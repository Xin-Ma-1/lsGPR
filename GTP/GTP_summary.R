source("~/functions/lsEM.R")
setwd("~/GTP/results")
sample = 25
p = 264; v = p*(p-1)/2
mse.all = list()
cov.all = list()
width.all = list()

for(k in 1:3){
  lam = k*0.05
  MSE = matrix(NA,sample,9)
  COVERAGE = matrix(NA,sample,5)
  WIDTH = matrix(NA,sample,5)

  for(i in 1:sample){
    load(paste0("GTP_lam",lam,"_rep",i,".Rdata"))
    MSE[i,] = mse
    COVERAGE[i,] = coverage
    WIDTH[i,] = width
  }
  mse.all[[k]] = MSE
  cov.all[[k]] = COVERAGE
  width.all[[k]] = WIDTH
}

results = matrix(NA,9,9)
results[1,9:1] = colMeans(mse.all[[1]])
results[4,9:1] = colMeans(mse.all[[2]])
results[7,9:1] = colMeans(mse.all[[3]])

results[2,5:1] = colMeans(cov.all[[1]])
results[5,5:1] = colMeans(cov.all[[2]])
results[8,5:1] = colMeans(cov.all[[3]])

results[3,5:1] = colMeans(width.all[[1]])
results[6,5:1] = colMeans(width.all[[2]])
results[9,5:1] = colMeans(width.all[[3]])

colnames(results) = c("sp-lsGPR","lsGPR","mfGPR","pcaGPR","edgeGPR","Horseshoe","Elastic net","Ridge","Lasso")
row.names(results) = rep(c("MSE","Coverage","Width"),3)
options(knitr.kable.NA = '')
knitr::kable(results,digits=3,format="latex")


# test on MSE between mfGPR and lsGPR
t.test(mse.all[[1]][,7],mse.all[[1]][,8],alternative="greater")
t.test(mse.all[[2]][,7],mse.all[[2]][,8],alternative="greater")
t.test(mse.all[[3]][,7],mse.all[[3]][,8],alternative="greater")

# test on MSE between mfGPR and sp-lsGPR
t.test(mse.all[[1]][,7],mse.all[[1]][,9],alternative="greater")
t.test(mse.all[[2]][,7],mse.all[[2]][,9],alternative="greater")
t.test(mse.all[[3]][,7],mse.all[[3]][,9],alternative="greater")
