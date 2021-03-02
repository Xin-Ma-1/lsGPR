ls.mcmc2pred.gpa<-function(a,a.te,U,U.te,Y,init,burn,nrun,prior){
  library(mvtnorm)
  library(MCMCpack)
  n=length(c(a))
  nn=length(c(a.te))
  a.full=c(a,a.te)
  U.full=append(U,U.te)
  Ea.full = Eu.full = matrix(rep(0,(n+nn)^2),(n+nn),(n+nn))
  for(i in 2:(n+nn)){
    for(j in 1:(i-1)){
      Ea.full[i,j]=Ea.full[j,i]=(a.full[i]-a.full[j])^2
      Eu.full[i,j]=Eu.full[j,i]=norm(U.full[[i]]-U.full[[j]],"F")^2
    }
  }
  Ea = Ea.full[1:n,1:n]
  sd.a = sd(Ea)
  Ea.full = Ea.full/sd.a
  Ea = Ea.full[1:n,1:n]
  Eu = Eu.full[1:n,1:n]
  sd.u = sd(Eu)
  Eu.full = Eu.full/sd.u
  Eu = Eu.full[1:n,1:n]

  if(missing(init)){init=c(1,1,1,1)}
  if(missing(burn)){burn=5000}
  if(missing(nrun)){nrun=8000}
  if(missing(prior)){prior=list(a.t=5,b.t=1,a.p1=2,b.p1=5,step=0.01)}
  tau=init[1]
  psi1=init[2]
  psia=init[3]
  psiu=init[4]
  phi=rep(0,n)
  PSI1=rep(NA,nrun); PSIA=rep(NA,nrun); PSIU=rep(NA,nrun); TAU=rep(NA,nrun)
  LOG=rep(NA,nrun)
  PHI = matrix(NA,n,nrun)
  PRED = matrix(NA,nn,nrun)
  for(g in 1:nrun){
    # update psi1
    EE = exp(-psia*Ea-psiu*Eu)
    a.psi1 = prior$a.p1 + 0.5*n
    b.psi1 = prior$b.p1 + 0.5*(matrix(phi,nrow=1)%*%solve(EE)%*%matrix(phi,ncol=1))
    psi1 = rinvgamma(1,shape=a.psi1,scale=b.psi1)
    PSI1[g] = psi1

    # update psia
    ev1 = diag(chol(EE))^2
    log.det = sum(log(ev1))
    quad = matrix(phi,nrow=1)%*%chol2inv(chol(EE))%*%matrix(phi,ncol=1)
    rand = exp(rnorm(1,log(psia),prior$step))
    E2 = exp(-rand*Ea-psiu*Eu)
    ev2 = diag(chol(E2))^2
    log.det2 = sum(log(ev2))
    quad2 = matrix(phi,nrow=1)%*%chol2inv(chol(E2))%*%matrix(phi,ncol=1)
    rho = min(-0.5*(log.det2 - log.det) - 0.5*(psi1^(-1))*(quad2 - quad), 0)
    acc=rbinom(1,1,exp(rho))
    if(acc==1){psia = rand}
    PSIA[g] = psia

    # update psiu
    EE = exp(-psia*Ea-psiu*Eu)
    ev1 = diag(chol(EE))^2
    log.det = sum(log(ev1))
    quad = matrix(phi,nrow=1)%*%chol2inv(chol(EE))%*%matrix(phi,ncol=1)
    rand = exp(rnorm(1,log(psiu),prior$step))
    E2 = exp(-psia*Ea-rand*Eu)
    ev2 = diag(chol(E2))^2
    log.det2 = sum(log(ev2))
    quad2 = matrix(phi,nrow=1)%*%chol2inv(chol(E2))%*%matrix(phi,ncol=1)
    rho = min(-0.5*(log.det2 - log.det) - 0.5*(psi1^(-1))*(quad2 - quad), 0)
    acc=rbinom(1,1,exp(rho))
    if(acc==1){psiu = rand}
    PSIU[g] = psiu

    # update tau
    a.tau = prior$a.t + 0.5*n
    b.tau = prior$b.t + 0.5*sum((c(Y)-c(phi))^2)
    tau = rgamma(1,shape=a.tau,rate=b.tau)
    TAU[g] = tau

    # update phi
    KK = psi1*exp(-psia*Ea-psiu*Eu)
    inv.KK = solve(KK)
    v.phi = solve(inv.KK + tau*diag(rep(1,n)))
    m.phi = v.phi%*%(tau*matrix(Y,ncol=1))
    phi = rmvnorm(1,mean=m.phi,sigma=v.phi)
    PHI[,g] = phi

    # log likelihood
    K = KK + tau^(-1)*diag(rep(1,n))
    inv.K = solve(K)
    LOG[g] = -0.5*log(det(K)) - 0.5*t(Y)%*%inv.K%*%Y

    # prediction
    cov.K = psi1*exp(-psia*Ea.full[1:n,(n+1):(n+nn)]-psiu*Eu.full[1:n,(n+1):(n+nn)])
    # pred = t(cov.K)%*%inv.K%*%matrix(Y,ncol=1)
    pred = t(cov.K)%*%inv.KK%*%matrix(phi,ncol=1)
    PRED[,g] = pred

    if(g%%1000==0){cat("iteration:",g,"tau:",tau,"psi1:",psi1,"psia:",psia,"psiu:",psiu,"loglik:",LOG[g],"\n")}
  }
  post.tau=mean(TAU[(burn+1):nrun])
  post.psi1=mean(PSI1[(burn+1):nrun])
  post.psia=mean(PSIA[(burn+1):nrun])
  post.psiu=mean(PSIU[(burn+1):nrun])
  post=c(post.tau,post.psi1,post.psia,post.psiu,sd.a,sd.u)
  return(list(post=post,post.loglik=LOG,post.mean=PHI,post.pred=PRED,TAU=TAU,PSI1=PSI1,PSIA=PSIA,PSIU=PSIU))
}

ls.mcmc3pred.gpa<-function(a,a.te,z,z.te,U,U.te,Y,init,burn,nrun,prior){
  library(mvtnorm)
  library(MCMCpack)
  n=length(c(a))
  nn=length(c(a.te))
  a.full=c(a,a.te)
  z.full=rbind(z,z.te)
  U.full=append(U,U.te)
  Ea.full = Ez.full = Eu.full = matrix(0,(n+nn),(n+nn))
  for(i in 2:(n+nn)){
    for(j in 1:(i-1)){
      Ea.full[i,j]=Ea.full[j,i]=(a.full[i]-a.full[j])^2
      Ez.full[i,j]=Ez.full[j,i]=sum((z.full[i,]-z.full[j,])^2)
      Eu.full[i,j]=Eu.full[j,i]=norm(U.full[[i]]-U.full[[j]],"F")^2
    }
  }
  Ea = Ea.full[1:n,1:n]
  sd.a = sd(Ea)
  Ea.full = Ea.full/sd.a
  Ea = Ea.full[1:n,1:n]

  Ez = Ez.full[1:n,1:n]
  sd.z = sd(Ez)
  Ez.full = Ez.full/sd.z
  Ez = Ez.full[1:n,1:n]

  Eu = Eu.full[1:n,1:n]
  sd.u = sd(Eu)
  Eu.full = Eu.full/sd.u
  Eu = Eu.full[1:n,1:n]

  if(missing(init)){init=c(1,1,0.1,0.1,0.01)}
  if(missing(burn)){burn=5000}
  if(missing(nrun)){nrun=8000}
  if(missing(prior)){prior=list(a.t=5,b.t=1,a.p1=2,b.p1=5,step=0.01)}
  tau=init[1]
  psi1=init[2]
  psia=init[3]
  psiz=init[4]
  psiu=init[5]
  phi=rep(0,n)
  PSI1=rep(NA,nrun); PSIA=rep(NA,nrun); PSIZ=rep(NA,nrun); PSIU=rep(NA,nrun); TAU=rep(NA,nrun)
  LOG=rep(NA,nrun)
  PHI = matrix(NA,n,nrun)
  PRED = matrix(NA,nn,nrun)
  for(g in 1:nrun){
    # update psi1
    EE = exp(-psia*Ea-psiz*Ez-psiu*Eu)
    a.psi1 = prior$a.p1 + 0.5*n
    b.psi1 = prior$b.p1 + 0.5*(matrix(phi,nrow=1)%*%solve(EE)%*%matrix(phi,ncol=1))
    psi1 = rinvgamma(1,shape=a.psi1,scale=b.psi1)
    PSI1[g] = psi1

    # update psia
    ev1 = diag(chol(EE))^2
    log.det = sum(log(ev1))
    quad = matrix(phi,nrow=1)%*%chol2inv(chol(EE))%*%matrix(phi,ncol=1)
    rand = exp(rnorm(1,log(psia),prior$step))
    E2 = exp(-rand*Ea-psiz*Ez-psiu*Eu)
    ev2 = diag(chol(E2))^2
    log.det2 = sum(log(ev2))
    quad2 = matrix(phi,nrow=1)%*%chol2inv(chol(E2))%*%matrix(phi,ncol=1)
    rho = min(-0.5*(log.det2 - log.det) - 0.5*(psi1^(-1))*(quad2 - quad), 0)
    acc=rbinom(1,1,exp(rho))
    if(acc==1){psia = rand}
    PSIA[g] = psia

    # update psiz
    EE = exp(-psia*Ea-psiz*Ez-psiu*Eu)
    ev1 = diag(chol(EE))^2
    log.det = sum(log(ev1))
    quad = matrix(phi,nrow=1)%*%chol2inv(chol(EE))%*%matrix(phi,ncol=1)
    rand = exp(rnorm(1,log(psiz),prior$step))
    E2 = exp(-psia*Ea-rand*Ez-psiu*Eu)
    ev2 = diag(chol(E2))^2
    log.det2 = sum(log(ev2))
    quad2 = matrix(phi,nrow=1)%*%chol2inv(chol(E2))%*%matrix(phi,ncol=1)
    rho = min(-0.5*(log.det2 - log.det) - 0.5*(psi1^(-1))*(quad2 - quad), 0)
    acc=rbinom(1,1,exp(rho))
    if(acc==1){psiz = rand}
    PSIZ[g] = psiz

    # update psiu
    EE = exp(-psia*Ea-psiz*Ez-psiu*Eu)
    ev1 = diag(chol(EE))^2
    log.det = sum(log(ev1))
    quad = matrix(phi,nrow=1)%*%chol2inv(chol(EE))%*%matrix(phi,ncol=1)
    rand = exp(rnorm(1,log(psiu),prior$step))
    E2 = exp(-psia*Ea-psiz*Ez-rand*Eu)
    ev2 = diag(chol(E2))^2
    log.det2 = sum(log(ev2))
    quad2 = matrix(phi,nrow=1)%*%chol2inv(chol(E2))%*%matrix(phi,ncol=1)
    rho = min(-0.5*(log.det2 - log.det) - 0.5*(psi1^(-1))*(quad2 - quad), 0)
    acc=rbinom(1,1,exp(rho))
    if(acc==1){psiu = rand}
    PSIU[g] = psiu

    # update tau
    a.tau = prior$a.t + 0.5*n
    b.tau = prior$b.t + 0.5*sum((c(Y)-c(phi))^2)
    tau = rgamma(1,shape=a.tau,rate=b.tau)
    TAU[g] = tau

    # update phi
    KK = psi1*exp(-psia*Ea-psiz*Ez-psiu*Eu)
    inv.KK = solve(KK)
    v.phi = solve(inv.KK + tau*diag(rep(1,n)))
    m.phi = v.phi%*%(tau*matrix(Y,ncol=1))
    phi = rmvnorm(1,mean=m.phi,sigma=v.phi)
    PHI[,g] = phi

    # log likelihood
    K = KK + tau^(-1)*diag(rep(1,n))
    inv.K = solve(K)
    LOG[g] = -0.5*log(det(K)) - 0.5*t(Y)%*%inv.K%*%Y

    # prediction
    cov.K = psi1*exp(-psia*Ea.full[1:n,(n+1):(n+nn)]-psiz*Ez.full[1:n,(n+1):(n+nn)]
                     -psiu*Eu.full[1:n,(n+1):(n+nn)])
    # pred = t(cov.K)%*%inv.K%*%matrix(Y,ncol=1)
    pred = t(cov.K)%*%inv.KK%*%matrix(phi,ncol=1)
    PRED[,g] = pred

    if(g%%1000==0){cat("iteration:",g,"tau:",tau,"psi1:",psi1,"psia:",psia,"psiz:",psiz,"psiu:",psiu,"loglik:",LOG[g],"\n")}
  }
  post.tau=mean(TAU[(burn+1):nrun])
  post.psi1=mean(PSI1[(burn+1):nrun])
  post.psia=mean(PSIA[(burn+1):nrun])
  post.psiz=mean(PSIZ[(burn+1):nrun])
  post.psiu=mean(PSIU[(burn+1):nrun])
  post=c(post.tau,post.psi1,post.psia,post.psiz,post.psiu,sd.a,sd.z,sd.u)
  return(list(post=post,post.loglik=LOG,post.mean=PHI,post.pred=PRED,TAU=TAU,PSI1=PSI1,PSIA=PSIA,PSIZ=PSIZ,PSIU=PSIU))
}

ls.mle2.gpa<-function(a,U,Y,init,maxiter,tol,param){
  n=length(c(a))
  Ea = Eu = matrix(rep(0,n^2),n,n)
  for(i in 2:n){
    for(j in 1:(i-1)){
      Ea[i,j]=Ea[j,i]=(a[i]-a[j])^2
      Eu[i,j]=Eu[j,i]=norm(U[[i]]-U[[j]],"F")^2
    }
  }
  sd.a = sd(Ea)
  sd.u = sd(Eu)
  Ea = Ea/sd.a
  Eu = Eu/sd.u
  if(missing(param)){param=list(a.p1=2, b.p1=5, a.a=2, b.a=5,a.u=2, b.u = 5, a.t = 5, b.t = 1)}
  if(missing(init)){init=c(1,1,1,1)}
  if(missing(maxiter)){maxiter=1000}
  if(missing(tol)){tol=10^(-8)}

  log.curr=log(init)
  init.loglik = -obj.log(vec=log.curr,Ea=Ea,Eu=Eu,Y=Y,param=param)
  iter=0
  lr=1
  for(g in 1:maxiter){
    iter = iter+1
    f0 = obj.log(vec=log.curr,Ea=Ea,Eu=Eu,Y=Y,param=param)
    df = grd.log(vec=log.curr,Ea=Ea,Eu=Eu,Y=Y,param=param)
    if(obj.log(vec=(log.curr-lr*df),Ea=Ea,Eu=Eu,Y=Y,param=param)<f0){
      log.curr = log.curr - lr*df
      if(abs(f0-obj.log(vec=log.curr,Ea=Ea,Eu=Eu,Y=Y,param=param))<tol){break}
      lr = 1.5*lr
      next
    }
    while(obj.log(vec=(log.curr-lr*df),Ea=Ea,Eu=Eu,Y=Y,param=param)>=f0){
      lr = 0.5*lr
    }
    log.curr = log.curr - lr*df
    if(abs(f0-obj.log(vec=log.curr,Ea=Ea,Eu=Eu,Y=Y,param=param))<tol){break}
  }
  update = exp(log.curr)
  loglik = -obj.log(vec=log.curr,Ea=Ea,Eu=Eu,Y=Y,param=param)
  cat("iteration step:",iter,"\n")
  return(list(vec=update,loglik=loglik,init.loglik=init.loglik,sd.a=sd.a,sd.u=sd.u))
}

ls.mle3.gpa<-function(a,z,U,Y,init,maxiter,tol,param){
  n=length(c(a))
  Ea = Ez = Eu = matrix(0,n,n)
  for(i in 2:n){
    for(j in 1:(i-1)){
      Ea[i,j]=Ea[j,i]=(a[i]-a[j])^2
      Ez[i,j]=Ez[j,i]=sum((z[i,]-z[j,])^2)
      Eu[i,j]=Eu[j,i]=norm(U[[i]]-U[[j]],"F")^2
    }
  }
  sd.a = sd(Ea)
  sd.z = sd(Ez)
  sd.u = sd(Eu)
  Ea = Ea/sd.a
  Ez = Ez/sd.z
  Eu = Eu/sd.u
  if(missing(param)){param=list(a.p1=2, b.p1=5, a.a=2, b.a=5,a.u=2, b.u = 5, a.t = 5, b.t = 1)}
  if(missing(init)){init=c(1,1,0.1,0.1,0.01)}
  if(missing(maxiter)){maxiter=1000}
  if(missing(tol)){tol=10^(-8)}

  log.curr=log(init)
  init.loglik = -obj.log2(vec=log.curr,Ea=Ea,Ez=Ez,Eu=Eu,Y=Y,param=param)
  iter=0
  lr=1
  for(g in 1:maxiter){
    iter = iter+1
    f0 = obj.log2(vec=log.curr,Ea=Ea,Ez=Ez,Eu=Eu,Y=Y,param=param)
    df = grd.log2(vec=log.curr,Ea=Ea,Ez=Ez,Eu=Eu,Y=Y,param=param)
    if(obj.log2(vec=(log.curr-lr*df),Ea=Ea,Ez=Ez,Eu=Eu,Y=Y,param=param)<f0){
      log.curr = log.curr - lr*df
      if(abs(f0-obj.log2(vec=log.curr,Ea=Ea,Ez=Ez,Eu=Eu,Y=Y,param=param))<tol){break}
      lr = 1.5*lr
      next
    }
    while(obj.log2(vec=(log.curr-lr*df),Ea=Ea,Ez=Ez,Eu=Eu,Y=Y,param=param)>=f0){
      lr = 0.5*lr
    }
    log.curr = log.curr - lr*df
    if(abs(f0-obj.log2(vec=log.curr,Ea=Ea,Ez=Ez,Eu=Eu,Y=Y,param=param))<tol){break}
  }
  update = exp(log.curr)
  loglik = -obj.log2(vec=log.curr,Ea=Ea,Ez=Ez,Eu=Eu,Y=Y,param=param)
  cat("iteration step:",iter,"\n")
  return(list(vec=update,loglik=loglik,init.loglik=init.loglik,sd.a=sd.a,sd.z=sd.z,sd.u=sd.u))
}



obj.log <- function(vec,Ea,Eu,Y,param){
  n=nrow(Ea)
  lt=vec[1]
  lp1=vec[2]
  lpa=vec[3]
  lpu=vec[4]
  E = lp1-exp(lpa)*Ea-exp(lpu)*Eu
  K = exp(E) + exp(-lt)*diag(rep(1,n))
  mLik = 0.5*log(det(K)) + 0.5*t(Y)%*%solve(K)%*%Y -param$a.t*lt + param$b.t*exp(lt) + param$a.p1*lp1 + param$b.p1*exp(-lp1)
}

grd.log <- function(vec,Ea,Eu,Y,param){
  n=nrow(Ea)
  lt=vec[1]
  lp1=vec[2]
  lpa=vec[3]
  lpu=vec[4]
  E = lp1-exp(lpa)*Ea-exp(lpu)*Eu
  K = exp(E) + exp(-lt)*diag(rep(1,n))
  inv.K = solve(K)
  alpha = inv.K%*%Y
  ka = inv.K-alpha%*%t(alpha)
  grd = rep(0,4)
  grd[1] = 0.5*sum(diag(-exp(-lt)*ka)) -param$a.t + param$b.t*exp(lt)
  grd[2] = 0.5*sum(diag(ka%*%exp(E))) +param$a.p1 - exp(-lp1)*param$b.p1
  grd[3] = 0.5*sum(diag(ka%*%(-exp(lpa)*exp(E)*Ea)))
  grd[4] = 0.5*sum(diag(ka%*%(-exp(lpu)*exp(E)*Eu)))
  grd
}

obj.log2 <- function(vec,Ea,Ez,Eu,Y,param){
  n=nrow(Ea)
  lt=vec[1]
  lp1=vec[2]
  lpa=vec[3]
  lpz=vec[4]
  lpu=vec[5]
  E = lp1-exp(lpa)*Ea-exp(lpz)*Ez-exp(lpu)*Eu
  K = exp(E) + exp(-lt)*diag(rep(1,n))
  mLik = 0.5*log(det(K)) + 0.5*t(Y)%*%solve(K)%*%Y
}

grd.log2 <- function(vec,Ea,Ez,Eu,Y,param){
  n=nrow(Ea)
  lt=vec[1]
  lp1=vec[2]
  lpa=vec[3]
  lpz=vec[4]
  lpu=vec[5]
  E = lp1-exp(lpa)*Ea-exp(lpz)*Ez-exp(lpu)*Eu
  K = exp(E) + exp(-lt)*diag(rep(1,n))
  inv.K = solve(K)
  alpha = inv.K%*%Y
  ka = inv.K-alpha%*%t(alpha)
  grd = rep(0,4)
  grd[1] = 0.5*sum(diag(-exp(-lt)*ka))
  grd[2] = 0.5*sum(diag(ka%*%exp(E)))
  grd[3] = 0.5*sum(diag(ka%*%(-exp(lpa)*exp(E)*Ea)))
  grd[4] = 0.5*sum(diag(ka%*%(-exp(lpz)*exp(E)*Ez)))
  grd[5] = 0.5*sum(diag(ka%*%(-exp(lpu)*exp(E)*Eu)))
  grd
}
