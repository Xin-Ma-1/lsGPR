auc_edge_to_node<-function(node.lab,edge.prop){
  p=length(node.lab)
  v=length(edge.prop)
  if(v!=p*(p-1)/2){cat("Error: node and edge lengths need to match!","\n")}
  node1 = which(node.lab==1)
  node0 = which(node.lab==0)

  thres = seq(0,1,by=0.01)
  se = sp = rep(NA, length(thres))
  for(k in 1:length(thres)){
    edge.select = 1*(edge.prop>=thres[k])
    node.select = which(rowSums(func.vtm(edge.select))>0)
    se[k] = sum(node.select %in% node1)/length(node1)
    sp[k] = 1-sum(node.select %in% node0)/length(node0)
  }

  auc = sum(sp*diff(c(0, 1-se)))
  return(auc)
}

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
