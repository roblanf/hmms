library("foreach")
library("doParallel")
registerDoParallel(cores=5)
source("BW-models.R")
source("BW-mix_model.R")

root = "./AliSim/taxa4/"
tree_files=list.files(root,pattern = "*.sitelh")

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


vit.acc=c()
post.acc=c()
time.taken=c()
for (j in 1:4){
res<-foreach (i=1:length(tree_files), .combine='comb', .multicombine=TRUE,.init=list(list(),list(),list())) %dopar% {

  start.time <- Sys.time()
  pred=predict_tree(paste(root,tree_files[i],sep=''),j) #model and "post.prob.tree."
  end.time <- Sys.time()
  time.taken <- as.numeric (end.time - start.time, units = "secs")
  v.pred<- pred[[1]];p.pred<-pred[[2]]
  true=rep(1:3, each=10000)#### change for taxa4
  true=as.character(true)
  true=paste("T",true,sep='')
  vit.acc=length(which(true==v.pred))/length(true)
  post.acc=length(which(true==p.pred))/length(true)
  list(time.taken,vit.acc,post.acc)
}


time.taken[j]=mean(unlist(res[1]))
vit.acc[j]=mean(unlist(res[2]))
post.acc[j]=mean(unlist(res[3]))
}

