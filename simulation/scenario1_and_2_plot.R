source("~/functions/lsEM.R")
source("~/functions/auc_edge_to_node.R")
library(pROC)
library(ggsci)
library(ggplot2)
mypal = pal_npg()(9)
mypal[1] = "#8A4198FF"
source("~/functions/shared_legend.R")
library(gridExtra)

setwd("~/simulation/results")
sample = 25
p = 100; v = p*(p-1)/2


####################
## scenario1   #####
####################
mse.all = list()
cov.all = list()
width.all = list()
auc.all = list()
act.vec = c(10,30,50)

for(k in 1:3){
  act = act.vec[k]
  node.lab = rep(0,p)
  node.lab[1:act] = 1

  MSE = matrix(NA,sample,9)
  COVERAGE = matrix(NA,sample,5)
  WIDTH = matrix(NA,sample,5)
  AUC = matrix(NA,sample,2)

  for(i in 1:sample){
    load(paste0("sce1_act",act,"_rep",i,".Rdata"))
    MSE[i,] = mse
    COVERAGE[i,] = coverage
    WIDTH[i,] = width
    AUC[i,1] = auc_edge_to_node(node.lab, hs.edge.prop)
    AUC[i,2] = auc(node.lab, lsGPR.prop)
  }
  mse.all[[k]] = MSE
  cov.all[[k]] = COVERAGE
  width.all[[k]] = WIDTH
  auc.all[[k]] = AUC
}

name = c(rep("Lasso",sample),rep("Ridge",sample),rep("Elastic Net",sample),rep("Horseshoe",sample),
         rep("edgeGPR",sample),rep("pcaGPR",sample),rep("mfGPR",sample),rep("lsGPR",sample),rep("sp-lsGPR",sample))
active = c(rep(10,sample*9),rep(30,sample*9),rep(50,sample*9))
long = c(c(mse.all[[1]]),c(mse.all[[2]]),c(mse.all[[3]]))
data1 = data.frame(pred.MSE=long,active=as.factor(active),method=factor(rep(name,3),level=unique(name)))
plot1 = ggplot(data=data1, aes(x=active, y=pred.MSE)) + geom_boxplot(aes(fill=method),outlier.shape=NA) +
  scale_fill_manual(values=mypal) + ylim(0,2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title=element_text(size=18),axis.text=element_text(size=16),
        legend.text=element_text(size=16),legend.title=element_text(size=18),plot.title=element_text(size=20,face='bold'))


name = c(rep("edgeGPR",sample),rep("pcaGPR",sample),rep("mfGPR",sample),rep("lsGPR",sample),rep("sp-lsGPR",sample))
active = c(rep(10,sample*5),rep(30,sample*5),rep(50,sample*5))
long = c(c(cov.all[[1]]),c(cov.all[[2]]),c(cov.all[[3]]))
data2 = data.frame(coverage=long,active=as.factor(active),method=factor(rep(name,3),level=unique(name)))
plot2 = ggplot(data=data2, aes(x=active, y=coverage)) + geom_boxplot(aes(fill=method),outlier.shape=NA) +
  ylim(0,1) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title=element_text(size=18),axis.text=element_text(size=16),
        legend.text=element_text(size=16),legend.title=element_text(size=18),plot.title=element_text(size=20,face='bold')) +
  scale_fill_manual(values=mypal[c(5,6,7,8,9)])

name = c(rep("edgeGPR",sample),rep("pcaGPR",sample),rep("mfGPR",sample),rep("lsGPR",sample),rep("sp-lsGPR",sample))
active = c(rep(10,sample*5),rep(30,sample*5),rep(50,sample*5))
long = c(c(width.all[[1]]),c(width.all[[2]]),c(width.all[[3]]))
data3 = data.frame(interval.width=long,active=as.factor(active),method=factor(rep(name,3),level=unique(name)))
plot3 = ggplot(data=data3, aes(x=active, y=interval.width)) + geom_boxplot(aes(fill=method),outlier.shape=NA) +
  ylim(0,6) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title=element_text(size=18),axis.text=element_text(size=16),
        legend.text=element_text(size=16),legend.title=element_text(size=18),plot.title=element_text(size=20,face='bold')) +
  scale_fill_manual(values=mypal[c(5,6,7,8,9)])

name = c(rep("Horseshoe",sample),rep("sp-lsGPR",sample))
active = c(rep(10,sample*2),rep(30,sample*2),rep(50,sample*2))
long = c(c(auc.all[[1]]),c(auc.all[[2]]),c(auc.all[[3]]))
data4 = data.frame(AUC=long,active=as.factor(active),method=factor(rep(name,3),level=unique(name)))
plot4 = ggplot(data=data4, aes(x=active, y=AUC)) + geom_boxplot(aes(fill=method),outlier.shape=NA) +
  ylim(0,1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title=element_text(size=18),axis.text=element_text(size=16),
        legend.text=element_text(size=16),legend.title=element_text(size=18),plot.title=element_text(size=20,face='bold')) +
  scale_fill_manual(values=mypal[c(4,9)])



####################
## scenario2   #####
####################
mse.all = list()
cov.all = list()
width.all = list()
auc.all = list()
act.vec = c(10,30,50)

for(k in 1:3){
  act = act.vec[k]
  node.lab = rep(0,p)
  node.lab[1:act] = 1

  MSE = matrix(NA,sample,9)
  COVERAGE = matrix(NA,sample,5)
  WIDTH = matrix(NA,sample,5)
  AUC = matrix(NA,sample,2)

  for(i in 1:sample){
    load(paste0("sce2_act",act,"_rep",i,".Rdata"))
    MSE[i,] = mse
    COVERAGE[i,] = coverage
    WIDTH[i,] = width
    AUC[i,1] = auc_edge_to_node(node.lab, hs.edge.prop)
    AUC[i,2] = auc(node.lab, lsGPR.prop)
  }
  mse.all[[k]] = MSE
  cov.all[[k]] = COVERAGE
  width.all[[k]] = WIDTH
  auc.all[[k]] = AUC
}

name = c(rep("Lasso",sample),rep("Ridge",sample),rep("Elastic Net",sample),rep("Horseshoe",sample),
         rep("edgeGPR",sample),rep("pcaGPR",sample),rep("mfGPR",sample),rep("lsGPR",sample),rep("sp-lsGPR",sample))
active = c(rep(10,sample*9),rep(30,sample*9),rep(50,sample*9))
long = c(c(mse.all[[1]]),c(mse.all[[2]]),c(mse.all[[3]]))
data5 = data.frame(pred.MSE=long,active=as.factor(active),method=factor(rep(name,3),level=unique(name)))
plot5 = ggplot(data=data5, aes(x=active, y=pred.MSE)) + geom_boxplot(aes(fill=method),outlier.shape=NA) +
  scale_fill_manual(values=mypal) + ylim(0,2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title=element_text(size=18),axis.text=element_text(size=16),
        legend.text=element_text(size=16),legend.title=element_text(size=18),plot.title=element_text(size=20,face='bold'))


name = c(rep("edgeGPR",sample),rep("pcaGPR",sample),rep("mfGPR",sample),rep("lsGPR",sample),rep("sp-lsGPR",sample))
active = c(rep(10,sample*5),rep(30,sample*5),rep(50,sample*5))
long = c(c(cov.all[[1]]),c(cov.all[[2]]),c(cov.all[[3]]))
data6 = data.frame(coverage=long,active=as.factor(active),method=factor(rep(name,3),level=unique(name)))
plot6 = ggplot(data=data6, aes(x=active, y=coverage)) + geom_boxplot(aes(fill=method),outlier.shape=NA) +
  ylim(0,1) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title=element_text(size=18),axis.text=element_text(size=16),
        legend.text=element_text(size=16),legend.title=element_text(size=18),plot.title=element_text(size=20,face='bold')) +
  scale_fill_manual(values=mypal[c(5,6,7,8,9)])

name = c(rep("edgeGPR",sample),rep("pcaGPR",sample),rep("mfGPR",sample),rep("lsGPR",sample),rep("sp-lsGPR",sample))
active = c(rep(10,sample*5),rep(30,sample*5),rep(50,sample*5))
long = c(c(width.all[[1]]),c(width.all[[2]]),c(width.all[[3]]))
data7 = data.frame(interval.width=long,active=as.factor(active),method=factor(rep(name,3),level=unique(name)))
plot7 = ggplot(data=data7, aes(x=active, y=interval.width)) + geom_boxplot(aes(fill=method),outlier.shape=NA) +
  ylim(0,6) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title=element_text(size=18),axis.text=element_text(size=16),
        legend.text=element_text(size=16),legend.title=element_text(size=18),plot.title=element_text(size=20,face='bold')) +
  scale_fill_manual(values=mypal[c(5,6,7,8,9)])

name = c(rep("Horseshoe",sample),rep("sp-lsGPR",sample))
active = c(rep(10,sample*2),rep(30,sample*2),rep(50,sample*2))
long = c(c(auc.all[[1]]),c(auc.all[[2]]),c(auc.all[[3]]))
data8 = data.frame(AUC=long,active=as.factor(active),method=factor(rep(name,3),level=unique(name)))
plot8 = ggplot(data=data8, aes(x=active, y=AUC)) + geom_boxplot(aes(fill=method),outlier.shape=NA) +
  ylim(0,1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none",
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title=element_text(size=18),axis.text=element_text(size=16),
        legend.text=element_text(size=16),legend.title=element_text(size=18),plot.title=element_text(size=20,face='bold')) +
  scale_fill_manual(values=mypal[c(4,9)])

setEPS()
postscript("~/simulation/results/sce1sce2_boxplots.eps",width=16.67,height=10.42)
grid_arrange_shared_legend(plot1,plot2,plot3,plot4,
                           plot5,plot6,plot7,plot8,
                           nrow=2,ncol=4)
dev.off()
