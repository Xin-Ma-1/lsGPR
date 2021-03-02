gpr.mcmc2pred<-function(X,X.te,Y,init,burn,nrun,prior){
  library(mvtnorm)
  library(MCMCpack)
  n=length(c(Y))
  nn=nrow(X.te)
  X.full=rbind(X,X.te)
  E.full = matrix(0,(n+nn),(n+nn))
  for(i in 2:(n+nn)){
    for(j in 1:(i-1)){
      E.full[i,j]=E.full[j,i]=sum((X.full[i,]-X.full[j,])^2)
    }
  }
  E = E.full[1:n,1:n]
  sd.x = sd(E)
  E.full = E.full/sd.x
  E = E.full[1:n,1:n]

  if(missing(prior)){prior=list(a.t=5,b.t=1,a.p1=2,b.p1=5,step=0.01)}
  tau=init[1]
  psi1=init[2]
  psi2=init[3]
  phi=rep(0,n)
  PSI1=rep(NA,nrun); PSI2=rep(NA,nrun); TAU=rep(NA,nrun)
  LOG=rep(NA,nrun)
  PHI = matrix(NA,n,nrun)
  PRED = matrix(NA,nn,nrun)
  for(g in 1:nrun){
    # update psi1
    EE = exp(-psi2*E) + diag(rep(1e-12,n))
    a.psi1 = prior$a.p1 + 0.5*n
    b.psi1 = prior$b.p1 + 0.5*(matrix(phi,nrow=1)%*%chol2inv(chol(EE))%*%matrix(phi,ncol=1))
    psi1 = rinvgamma(1,shape=a.psi1,scale=b.psi1)
    PSI1[g] = psi1

    # update psi2
    ev1 = diag(chol(psi1*EE + tau^(-1)*diag(rep(1,n))))^2
    log.det = sum(log(ev1))
    quad = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*EE + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
    rand = exp(rnorm(1,log(psi2),prior$step))
    E2 = exp(-rand*E)
    ev2 = diag(chol(psi1*E2 + tau^(-1)*diag(rep(1,n))))^2
    log.det2 = sum(log(ev2))
    quad2 = matrix(Y,nrow=1)%*%chol2inv(chol(psi1*E2 + tau^(-1)*diag(rep(1,n))))%*%matrix(Y,ncol=1)
    rho = min(-0.5*(log.det2 - log.det) - 0.5*(psi1^(-1))*(quad2 - quad), 0)
    acc=rbinom(1,1,exp(rho))
    if(acc==1){psi2 = rand}
    PSI2[g] = psi2

    # update tau
    a.tau = prior$a.t + 0.5*n
    b.tau = prior$b.t + 0.5*sum((Y-phi)^2)
    tau = rgamma(1,shape=a.tau,rate=b.tau)
    TAU[g] = tau

    # update phi
    KK = psi1*exp(-psi2*E) + diag(rep(1e-12,n))
    inv.KK = chol2inv(chol(KK))
    v.phi = chol2inv(chol(inv.KK + tau*diag(rep(1,n))))
    v.phi = (v.phi+t(v.phi))/2
    m.phi = v.phi%*%(tau*matrix(Y,ncol=1))
    phi = rmvnorm(1,mean=m.phi,sigma=v.phi)
    PHI[,g] = phi

    # in-sample prediction
    cov.K = psi1*exp(-psi2*E.full[1:n,(n+1):(n+nn)])
    pred = t(cov.K)%*%inv.KK%*%matrix(phi,ncol=1)
    PRED[,g] = pred

    if(g%%1000==0){cat("iteration:",g,"tau:",tau,"psi1:",psi1,"psi2:",psi2,"\n")}
  }
  post.tau=mean(TAU[(burn+1):nrun])
  post.psi1=mean(PSI1[(burn+1):nrun])
  post.psi2=mean(PSI2[(burn+1):nrun])
  post=c(post.tau,post.psi1,post.psi2,sd.x)
  return(list(post=post,post.mean=PHI,post.pred=PRED,TAU=TAU,PSI1=PSI1,PSI2=PSI2))
}


gpr.mle2<-function(X,Y,init,maxiter,tol,param){
  n=length(c(Y))
  E = matrix(rep(0,n^2),n,n)
  for(i in 2:n){
    for(j in 1:(i-1)){
      E[i,j]=E[j,i]=sum((X[i,]-X[j,])^2)
    }
  }
  sd.x = sd(E)
  E = E/sd.x
  if(missing(param)){param=list(a.p1=2, b.p1=5, a.p2=2, b.p2 = 5, a.t = 5, b.t = 1)}
  if(missing(init)){init=c(1,0.1,0.01)}
  if(missing(maxiter)){maxiter=100}
  if(missing(tol)){tol=10^(-8)}

  log.curr=log(init)
  init.loglik = -gpr.obj(vec=log.curr,E=E,Y=Y,param=param)
  iter=0
  lr=1
  for(g in 1:maxiter){
    iter = iter+1
    f0 = gpr.obj(vec=log.curr,E=E,Y=Y,param=param)
    df = gpr.grd(vec=log.curr,E=E,Y=Y,param=param)
    if(gpr.obj(vec=(log.curr-lr*df),E=E,Y=Y,param=param)<f0){
      log.curr = log.curr - lr*df
      if(abs(f0-gpr.obj(vec=log.curr,E=E,Y=Y,param=param))<tol){break}
      lr = 1.5*lr
      next
    }
    while(gpr.obj(vec=(log.curr-lr*df),E=E,Y=Y,param=param)>=f0){
      lr = 0.5*lr
    }
    log.curr = log.curr - lr*df
    if(abs(f0-gpr.obj(vec=log.curr,E=E,Y=Y,param=param))<tol){break}
  }
  update = exp(log.curr)
  loglik = -gpr.obj(vec=log.curr,E=E,Y=Y,param=param)
  cat("iteration step:",iter,"\n")
  return(list(vec=update,loglik=loglik,init.loglik=init.loglik,sd.x=sd.x))
}


gpr.obj <- function(vec,E,Y,param){
  n=nrow(E)
  lt=vec[1]
  lp1=vec[2]
  lp2=vec[3]
  EE = lp1-exp(lp2)*E
  K = exp(EE) + exp(-lt)*diag(rep(1,n))
  mLik = 0.5*log(det(K)) + 0.5*t(Y)%*%solve(K)%*%Y
  mLik
}

gpr.grd <- function(vec,E,Y,param){
  n=nrow(E)
  lt=vec[1]
  lp1=vec[2]
  lp2=vec[3]
  EE = lp1-exp(lp2)*E
  K = exp(EE) + exp(-lt)*diag(rep(1,n))
  inv.K = solve(K)
  alpha = inv.K%*%Y
  ka = inv.K-alpha%*%t(alpha)
  grd = rep(0,3)
  grd[1] = 0.5*sum(diag(-exp(-lt)*ka))
  grd[2] = 0.5*sum(diag(ka%*%exp(EE)))
  grd[3] = 0.5*sum(diag(ka%*%(-exp(lp2)*exp(EE)*E)))
  grd
}
