options(echo=TRUE) # if you want to see commands in output file

args = commandArgs(trailingOnly=T)
act = as.numeric(args[1])
act
ind = as.numeric(args[2])
ind

eff.d = 3
n=50; nn=50; p=100; v=p*(p-1)/2
library(MCMCpack)
library(glmnet)
library(doParallel)
library(foreach)
library(monomvn)
library(dimRed)
library(coda)
library(fastICA)
library(NMF)
library(igraph)
library(loe)
library(RSpectra)
library(RANN)
source("~/functions/lsEM.R")
source("~/functions/lsGPR_VS.R")
source("~/functions/GPRfuncs.R")
prior=list(a.t=0.1,b.t=80,a.p1=0.1,b.p1=80,step=0.01)
# generate network and scalar outcome
G = matrix(rep(0,v*(n+nn)),(n+nn),v)
for(i in 1:(n+nn)){
  nei = sample(4:12,1)
  g = watts.strogatz.game(1,100,nei,0.5)
  adj = as.matrix(as_adj(g))
  G[i,] = func.mtv(adj)
}
U0 = list()
a0 = rep(0,(n+nn))
for(i in 1:(n+nn)){
  EM = EM_toU(G[i,],d=eff.d,maxiter=500)
  U0[[i]] = EM$U
  a0[i] = EM$a
}
Ytr = rep(NA,n); Yte = rep(NA,nn)
b0 = 1/2
b = matrix(0,act,eff.d)
for(k in 1:(act/2)){b[k,2]=1}
for(k in (act/2+1):act){b[k,3]=-1}
for(i in 1:n){
  sum = 0
  for(k in 2:act){sum = sum + sum(b[k,]*U0[[i]][k,])}
  Ytr[i] = b0*a0[i] + exp(sum(b[1,]*U0[[i]][1,])) + sum + rnorm(1,sd=0.2)
}
for(i in 1:nn){
  sum = 0
  for(k in 2:act){sum = sum + sum(b[k,]*U0[[(i+n)]][k,])}
  Yte[i] = b0*a0[(i+n)] + exp(sum(b[1,]*U0[[(i+n)]][1,])) + sum + rnorm(1,sd=0.2)
}
m0 = mean(Ytr); sd0 = sd(Ytr)
Ytr = (Ytr-m0)/sd0; Yte = (Yte-m0)/sd0




###########################
### estimation        #####
###########################
mse = rep(NA,9)
coverage = rep(NA,5)
width = rep(NA,5)

## Lasso
cv.las.fit = cv.glmnet(G[1:n,],Ytr,alpha=1,family='gaussian',intercept=T)
las.fit=glmnet(G[1:n,],Ytr,alpha=1,family='gaussian',intercept=T,lambda=cv.las.fit$lambda.min)
las.pred = predict(las.fit, newx = G[(n+1):(n+nn),])
mse[1] = mean((Yte - las.pred)^2)/var(Yte)
cat("lasso:",mse[1],"\n")

## Ridge
cv.rid.fit = cv.glmnet(G[1:n,],Ytr,alpha=0,family='gaussian',intercept=T)
rid.fit=glmnet(G[1:n,],Ytr,alpha=0,family='gaussian',intercept=T,lambda=cv.rid.fit$lambda.min)
rid.pred = predict(rid.fit, newx = G[(n+1):(n+nn),])
mse[2] = mean((Yte - rid.pred)^2)/var(Yte)
cat("ridge:",mse[2],"\n")

## Elastic Net
a <- seq(0.05, 0.95, 0.05)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(G[1:n,],Ytr,family='gaussian', paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
}
cv.ela <- search[search$cvm == min(search$cvm), ]
ela.fit=glmnet(G[1:n,],Ytr,alpha=cv.ela$alpha,lambda=cv.ela$lambda.min,family='gaussian',intercept=T)
ela.pred=predict(ela.fit,newx=G[(n+1):(n+nn),])
mse[3]=mean((Yte-ela.pred)^2)/var(Yte)
cat("elastic.net:",mse[3],"\n")

## Bayesian Horseshoe "BHS"
hs.fit = bhs(G[1:n,],Ytr,T=5000,normalize=F,verb=0)
hs.mu = mean(hs.fit$mu[-seq(2500)])
hs.beta = apply(hs.fit$beta[-seq(2500),],2,mean)
hs.pred = apply(G[(n+1):(n+nn),],1,function(x) sum(x*hs.beta)+hs.mu)
mse[4] = mean((hs.pred-Yte)^2)/var(Yte)
cat("horseshoe:",mse[4],"\n","\n")
hs.edge.prop = colMeans(hs.fit$beta[-seq(2500),]!=0)

## GPR with full edge set
raw.fit=gpr.mle2(G[1:n,],Ytr,init=c(1,1,1),maxiter=1000)
raw.mcmc.fit = gpr.mcmc2pred(G[1:n,],G[(n+1):(n+nn),],Ytr,init=raw.fit$vec,5000,10000,prior=prior)
raw.mcmc.pred = rowMeans(raw.mcmc.fit$post.pred[,-seq(5000)])
mse[5] = mean((raw.mcmc.pred-Yte)^2)/var(Yte)
cat("rawGPR mse:",mse[5],"\n")
CI = HPDinterval(mcmc(t(raw.mcmc.fit$post.pred[,5001:10000])))
pred.in = (Yte>CI[,1])&(Yte<CI[,2])
coverage[1] = mean(pred.in)
cat("rawGPR coverage:",coverage[1],"\n")
width[1] = mean(CI[,2]-CI[,1])
cat("rawGPR width:",width[1],"\n","\n")

## PCA GPR
pca.G = prcomp(G)$x[,1:eff.d]
pca.fit=gpr.mle2(pca.G[1:n,],Ytr,init=c(1,1,1),maxiter=1000)
pca.mcmc.fit = gpr.mcmc2pred(pca.G[1:n,],pca.G[(n+1):(n+nn),],Ytr,init=pca.fit$vec,5000,10000,prior=prior)
pca.mcmc.pred = rowMeans(pca.mcmc.fit$post.pred[,-seq(5000)])
mse[6] = mean((pca.mcmc.pred-Yte)^2)/var(Yte)
cat("pcaGPR mse:",mse[6],"\n")
CI = HPDinterval(mcmc(t(pca.mcmc.fit$post.pred[,5001:10000])))
pred.in = (Yte>CI[,1])&(Yte<CI[,2])
coverage[2] = mean(pred.in)
cat("pcaGPR coverage:",coverage[2],"\n")
width[2] = mean(CI[,2]-CI[,1])
cat("pcaGPR width:",width[2],"\n","\n")

## manifold GPR
d.est = eff.d
leim=LaplacianEigenmaps(stdpars=list(ndim=d.est,sparse='knn',knn=20,eps=0.1,t=Inf,norm=T))
leim.data=dimRedData(G)
emb=leim@fun(leim.data, leim@stdpars)
red.G = emb@data@data
red.G = scale(red.G)
gp.fit=gpr.mle2(red.G[1:n,],Ytr,init=c(1,1,1),maxiter=1000)
gp.mcmc.fit = gpr.mcmc2pred(red.G[1:n,],red.G[(n+1):(n+nn),],Ytr,init=gp.fit$vec,5000,10000,prior=prior)
gp.mcmc.pred = rowMeans(gp.mcmc.fit$post.pred[,-seq(5000)])
mse[7] = mean((gp.mcmc.pred-Yte)^2)/var(Yte)
cat("mfGPR mse:",mse[7],"\n")
CI = HPDinterval(mcmc(t(gp.mcmc.fit$post.pred[,5001:10000])))
pred.in = (Yte>CI[,1])&(Yte<CI[,2])
coverage[3] = mean(pred.in)
cat("mfGPR coverage:",coverage[3],"\n")
width[3] = mean(CI[,2]-CI[,1])
cat("mfGPR width:",width[3],"\n","\n")

# our lsGPR
U = list()
a = rep(0,(n+nn))
for(i in 1:(n+nn)){
  EM = EM_toU(G[i,],d=eff.d,maxiter=500)
  U[[i]] = EM$U
  a[i] = EM$a
}

# lsGPR
ls.mle.fit = ls.mle2.gpa(a[1:n],U[1:n],Ytr,init=c(1,1,1,1),maxiter=1000)
ls.mcmc.fit = ls.mcmc2pred.gpa(a[1:n],a[(n+1):(n+nn)],U[1:n],U[(n+1):(n+nn)],Ytr,init=ls.mle.fit$vec,5000,10000,
                               prior=list(a.t=5,b.t=1,a.p1=2,b.p1=5,step=0.01))
ls.mcmc.pred = rowMeans(ls.mcmc.fit$post.pred[,-seq(5000)])
mse[8] = mean((ls.mcmc.pred-Yte)^2)/var(Yte)
cat("lsGPR mse:",mse[8],"\n")
CI = HPDinterval(mcmc(t(ls.mcmc.fit$post.pred[,-seq(5000)])),prob=0.95)
pred.in = (Yte>CI[,1])&(Yte<CI[,2])
coverage[4] = mean(pred.in)
cat("lsGPR coverage:",coverage[4],"\n")
width[4] = mean(CI[,2]-CI[,1])
cat("lsGPR width:",width[4],"\n","\n")

# sp-lsGPR
lsGPR.fit = lsGPR_VS(a=a[1:n],a.te=a[(n+1):(n+nn)],U=U[1:n],U.te=U[(n+1):(n+nn)],Y=Ytr,
                     prior=list(a.t=5,b.t=1,a.p1=2,b.p1=5,step=0.01,a.p=1,b.p=40))
lsGPR.pred = rowMeans(lsGPR.fit$post.pred[,-seq(5000)])
mse[9] = mean((lsGPR.pred-Yte)^2)/var(Yte)
cat("lsGPR mse:",mse[9],"\n")
CI = HPDinterval(mcmc(t(lsGPR.fit$post.pred[,-seq(5000)])),prob=0.95)
pred.in = (Yte>CI[,1])&(Yte<CI[,2])
coverage[5] = mean(pred.in)
cat("lsGPR coverage:",coverage[5],"\n")
width[5] = mean(CI[,2]-CI[,1])
cat("lsGPR width:",width[5],"\n","\n")
lsGPR.prop = rowMeans(lsGPR.fit$GAMMA[,-seq(5000)])

mse
coverage
width

save(mse,coverage,width,hs.edge.prop,lsGPR.prop,
     file=paste0("~/simulation/results/sce2_act",act,"_rep",ind,".Rdata"))
