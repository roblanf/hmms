library("aphid") 
library("testit")

root = "./AliSim/taxa4/"
tree_files=list.files(root,pattern = "*.sitelh")

predict_tree <- function(file){
  
  data=read.table(file, header=TRUE, sep = ",")
  numTrees=(ncol(data)-5)/2
  print(numTrees)
  
  states <- c("Begin",paste("T",1:numTrees,sep=''))
  residues <- letters[1:(numTrees+3)]

  variable="log.like.tree."
  #variable="post.prob.tree."
  variables=paste(variable,as.character(1:numTrees),sep='')
  seq=head(residues,numTrees)[apply(data[variables], 1, which.max)]
  seq[data$sameParsimony==1]=letters[(numTrees+3)]
  seq[data$isInformative==0]=letters[(numTrees+2)]
  seq[data$isConstant==1]=letters[(numTrees+1)]
  print(table(seq))

  
  # Initial HMM
  #'   ### Define the transition probability matrix
  A <- matrix(c(rep(.01/(numTrees-1),(numTrees+1)**2)),nrow=numTrees+1,ncol=numTrees+1,byrow = TRUE)
  diag(A) <- .99
  A[1,]<-1/numTrees
  A[,1]<-0
  dimnames(A) <- list(from = states, to = states)
  #'   ### Define the emission probability matrix
  E=matrix(c(rep(.5/(numTrees+2),numTrees)),nrow=numTrees,ncol=numTrees+3,byrow = TRUE)
  diag(E) <- .5
  dimnames(E) <- list(states = states[-1], residues = residues)
  #'   ### Build the HMM object
  hmm <- structure(list(A = A, E = E), class = "HMM")

  
  # Baum-Welch
  warn=has_warning(bw <- train(hmm,seq,method = "BaumWelch",maxiter=10000,logspace = FALSE,quiet=TRUE))
  if (warn==TRUE){
    print("Not converged")
    #bw <- train(hmm,seq,method = "Viterbi",logspace = FALSE,quiet=TRUE)
  }
  
  # get viterbi path
  viterbi = Viterbi(bw,seq)
  viterbi.path=rownames(bw$E)[viterbi$path + 1]
  
  #get posterior highest state
  post.prob = posterior(bw,seq)
  post.path=tail(states,numTrees)[apply(post.prob, 2, which.max)]
  return(list(viterbi.path,post.path))
}

vit.acc=c()
post.acc=c()
time.taken=c()

#multiple files prediction and accuracy
for (i in 1:length(tree_files)){
  
  start.time <- Sys.time()
  pred=predict_tree(paste(root,tree_files[i],sep=''))
  end.time <- Sys.time()
  time.taken[i] <- end.time - start.time
  v.pred<- pred[[1]];p.pred<-pred[[2]]
  print(table(v.pred))
  print(table(p.pred))
  struct_file=paste(strsplit(tree_files[i],"\\.")[[1]][1],".regions.txt",sep="")
  struct=read.table(paste("./regions.txt/",struct_file,sep=""))
  true=c(rep(struct$V1[1],5000),rep(struct$V1[2],5000),rep(struct$V1[3],5000),rep(struct$V1[4],5000))
  true=as.character(true)
  true=paste("T",true,sep='')
  print(paste("Viterbi Accuracy:",tree_files[i],length(which(true==v.pred))/length(true)))
  print(paste("Posterior Accuracy:",tree_files[i],length(which(true==p.pred))/length(true)))
  vit.acc[i]=length(which(true==v.pred))/length(true)
  post.acc[i]=length(which(true==p.pred))/length(true)
}

for (i in 1:length(tree_files)){
  
  start.time <- Sys.time()
  pred=predict_tree(paste(root,tree_files[i],sep=''))
  end.time <- Sys.time()
  time.taken[i] <- as.numeric (end.time - start.time, units = "secs")
  v.pred<- pred[[1]];p.pred<-pred[[2]]
  print(table(v.pred))
  print(table(p.pred))
  true=rep(1:3, each=10000)
  true=as.character(true)
  true=paste("T",true,sep='')
  print(paste("Viterbi Accuracy:",tree_files[i],length(which(true==v.pred))/length(true)))
  print(paste("Posterior Accuracy:",tree_files[i],length(which(true==p.pred))/length(true)))
  vit.acc[i]=length(which(true==v.pred))/length(true)
  post.acc[i]=length(which(true==p.pred))/length(true)
}