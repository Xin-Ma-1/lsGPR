func.mtv<-function(mat){
  lower=mat
  lower[upper.tri(lower,diag=T)]<-NA
  c.lower=c(lower)
  c.lower[!is.na(c.lower)]
}

func.vtm<-function(vec){
  v = length(c(vec))
  n = (1+sqrt(8*v+1))/2
  lower=matrix(rep(0,n^2),n)
  lower[lower.tri(lower, diag=FALSE)]=vec
  lower+t(lower)
}

EM_toU <- function(G, d=3, u1=0.5, sig.a=2, sig.u=2, tol=10^(-10), maxiter=200){
  library(igraph)
  G = c(G);v = length(G); p = (1+sqrt(8*v+1))/2
  sG = func.vtm(G)
  g = graph_from_adjacency_matrix(sG)
  D = distances(g, mode='all')
  D[D==Inf] = max(D[D!=Inf])+10
  # D[D==Inf] = p
  z0 = cmdscale(D,(d-1))

  u1=0.5
  U = cbind(u1,matrix(z0,ncol=(d-1)))
  # a = -log(1/mean(G)-1)
  a = 0
  step = 0; diff = 1
  while(step<=maxiter & diff>tol){
    # E step
    dotU = func.mtv(U%*%t(U))
    phi = a + dotU
    W = (exp(phi)-1)/(2*phi*(1+exp(phi)))
    mW = func.vtm(W)

    # M step
    # max a
    old.a = a
    a = sum(G-W*dotU-0.5)/(1/sig.a + sum(W))

    old.U = U

    # max U[k,-1]
    for(k in 1:p){
      tmp1 = t(U[-k,-1])%*%diag(mW[-k,k])%*%U[-k,-1] + (1/sig.u)*diag(rep(1,(d-1)))
      tmp2 = matrix(sG[-k,k]-mW[-k,k]*(a+u1^2)-0.5,nrow=1)%*%U[-k,-1]
      U[k,-1] = tmp2%*%solve(tmp1)
    }

    diff.U = max(abs(c(old.U-U)))
    diff.a = abs(old.a-a)
    diff = max(diff.U, diff.a)
    step = step + 1
  }
  list(U=U,a=a)
}

EM_toU_ver0 <- function(G, d=3, tol=10^(-10), maxiter=200){
  library(igraph)
  G = c(G);v = length(G); p = (1+sqrt(8*v+1))/2
  sG = func.vtm(G)
  g = graph_from_adjacency_matrix(sG)
  D = distances(g, mode='all')
  D[D==Inf] = max(D[D!=Inf])+1
  z0 = cmdscale(D,(d-1))

  u1=0.5
  U = cbind(u1,matrix(z0,ncol=(d-1)))
  a = -log(1/mean(G)-1)
  step = 0; diff = 1
  while(step<=maxiter & diff>tol){
    # E step
    dotU = func.mtv(U%*%t(U))
    phi = a + dotU
    W = (exp(phi)-1)/(2*phi*(1+exp(phi)))
    mW = func.vtm(W)

    # M step
    # max a
    old.a = a
    a = sum(G-W*dotU-0.5)/sum(W)

    old.U = U

    # max U[k,-1]
    for(k in 1:p){
      tmp1 = t(U[-k,-1])%*%diag(mW[-k,k])%*%U[-k,-1] + 0.5*diag(rep(1,(d-1)))
      tmp2 = matrix(sG[-k,k]-mW[-k,k]*(a+u1^2)-0.5,nrow=1)%*%U[-k,-1]
      U[k,-1] = tmp2%*%solve(tmp1)
    }

    diff.U = norm(old.U-U,"F")
    diff.a = abs(old.a-a)
    diff = max(diff.U, diff.a)
    step = step + 1
  }
  list(U=U,a=a)
}
