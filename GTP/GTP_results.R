options(echo=TRUE) # if you want to see commands in output file

args = commandArgs(trailingOnly=T)
lam = as.numeric(args[1])
lam
ind = as.numeric(args[2])
ind

load(paste0("~/GTP/data/gradyfc_lam",lam,".Rdata"))

eff.d = 10
p=264; v=p*(p-1)/2
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
source("~/functions/lsEM_old.R")
source("~/functions/lsGPR_VS.R")
source("~/functions/GPRfuncs.R")
prior=list(a.t=0.1,b.t=80,a.p1=0.1,b.p1=80,step=0.01)
# load scalar outcome and other covariates
clinic = read.csv("~/GTP/data/grady_clinical_data.csv",header=T)
list = which(clinic$CDRISC10total!='')
Y = scale(clinic$CDRISC10total[list])
age = clinic$gtp_mri_participant_age
age = scale(age[list])
tei = clinic$TEI_TOTAL_TYPES_Experienced_som1 # missing subj33
tei[33] = mean(tei[-33])
tei = scale(tei[list])
ctq = clinic$CTQTOT
ctq = scale(log(ctq[list]))
z = cbind(age,tei,ctq)
G = G[list,]
n=round(0.5*length(list))
nn=length(list)-n

# sample train and test indices
load("~/GTP/data/grady_randomsplit_orders.Rdata")
order = orders[i,]
train = order[1:n]
test = order[-seq(n)]
Ytr = Y[train]; Yte = Y[test]
m0 = mean(Ytr); sd0 = sd(Ytr)
Ytr = (Ytr-m0)/sd0; Yte = (Yte-m0)/sd0
ztr = z[train,]; zte = z[test,]



###########################
### estimation        #####
###########################
mse = rep(NA,9)
coverage = rep(NA,5)
width = rep(NA,5)

## Lasso
cv.las.fit = cv.glmnet(cbind(ztr,cbind(ztr,G[train,])),Ytr,alpha=1,family='gaussian',intercept=T)
las.fit=glmnet(cbind(ztr,G[train,]),Ytr,alpha=1,family='gaussian',intercept=T,lambda=cv.las.fit$lambda.min)
las.pred = predict(las.fit, newx = cbind(zte,G[test,]))
mse[1] = mean((Yte - las.pred)^2)/var(Yte)
cat("lasso:",mse[1],"\n")

## Ridge
cv.rid.fit = cv.glmnet(cbind(ztr,G[train,]),Ytr,alpha=0,family='gaussian',intercept=T)
rid.fit=glmnet(cbind(ztr,G[train,]),Ytr,alpha=0,family='gaussian',intercept=T,lambda=cv.rid.fit$lambda.min)
rid.pred = predict(rid.fit, newx = cbind(zte,G[test,]))
mse[2] = mean((Yte - rid.pred)^2)/var(Yte)
cat("ridge:",mse[2],"\n")

## Elastic Net
a <- seq(0.05, 0.95, 0.05)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(cbind(ztr,G[train,]),Ytr,family='gaussian', paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
}
cv.ela <- search[search$cvm == min(search$cvm), ]
ela.fit=glmnet(cbind(ztr,G[train,]),Ytr,alpha=cv.ela$alpha,lambda=cv.ela$lambda.min,family='gaussian',intercept=T)
ela.pred=predict(ela.fit,newx=cbind(zte,G[test,]))
mse[3]=mean((Yte-ela.pred)^2)/var(Yte)
cat("elastic.net:",mse[3],"\n")

## Bayesian Horseshoe "BHS"
hs.fit = bhs(cbind(ztr,G[train,]),Ytr,T=5000,normalize=F,verb=0)
hs.mu = mean(hs.fit$mu[-seq(2500)])
hs.beta = apply(hs.fit$beta[-seq(2500),],2,mean)
hs.pred = apply(cbind(zte,G[test,]),1,function(x) sum(x*hs.beta)+hs.mu)
mse[4] = mean((hs.pred-Yte)^2)/var(Yte)
cat("horseshoe:",mse[4],"\n","\n")
hs.edge.prop = colMeans(hs.fit$beta[-seq(2500),]!=0)

## GPR with full edge set
raw.fit=gpr.mle2(cbind(ztr,G[train,]),Ytr,init=c(1,1,1),maxiter=1000)
raw.mcmc.fit = gpr.mcmc2pred(cbind(ztr,G[train,]),cbind(zte,G[test,]),Ytr,init=raw.fit$vec,5000,10000,prior=prior)
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
pca.fit=gpr.mle2(cbind(ztr,pca.G[train,]),Ytr,init=c(1,1,1),maxiter=1000)
pca.mcmc.fit = gpr.mcmc2pred(cbind(ztr,pca.G[train,]),cbind(zte,pca.G[test,]),Ytr,init=pca.fit$vec,5000,10000,prior=prior)
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
gp.fit=gpr.mle2(cbind(ztr,red.G[train,]),Ytr,init=c(1,1,1),maxiter=1000)
gp.mcmc.fit = gpr.mcmc2pred(cbind(ztr,red.G[train,]),cbind(zte,red.G[test,]),Ytr,init=gp.fit$vec,5000,10000,prior=prior)
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
ls.mle.fit = ls.mle3.gpa(atr,ztr,Utr,Ytr,init=c(1,1,1,1,1),maxiter=1000)
ls.mcmc.fit = ls.mcmc3pred.gpa(atr,ate,ztr,zte,Utr,Ute,Ytr,init=ls.mle.fit$vec,5000,10000,prior=prior)
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
lsGPR.fit = lsGPR_VS(a=a[train],a.te=a[test],U=U[train],U.te=U[test],
                     z=ztr,z.te=zte,Y=Ytr,
                     prior=list(a.t=0.1,b.t=80,a.p1=0.1,b.p1=80,step=0.01,a.p=1,b.p=1))
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
     file=paste0("~/GTP/results/GTP_lam",lam,"_rep",ind,".Rdata"))
