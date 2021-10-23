library("aphid") 
library("testit")

tree_files=list.files("./tip4.sitelh",pattern = "*.sitelh")
#struct_files=list.files("./dataset",pattern = "*.txt")

predict_tree <- function(file){
  
  data=read.table(paste("./tip4.sitelh/",file,sep=''), header=TRUE, sep = ",")
  
  states <- c("Begin","T1", "T2")
  residues <- c("a","b")
  seq=head(residues,2)[apply(data[c("post.prob.tree.1","post.prob.tree.2")], 1, which.max)]
  #print(colSums(data != 0))
  print(table(seq))
  
  # Initial HMM
  #'   ### Define the transition probability matrix
  A <- matrix(c(0,0.5,0.5,0,.9,.1,0,.1,.9),nrow=3,ncol=3,byrow = TRUE)
  dimnames(A) <- list(from = states, to = states)
  #'   ### Define the emission probability matrix
  E=matrix(c(.9,.1,.1,.9),nrow=2,ncol=2,byrow = TRUE)
  dimnames(E) <- list(states = states[-1], residues = residues)
  #'   ### Build the HMM object
  hmm <- structure(list(A = A, E = E), class = "HMM")
  # hmm <- derive.HMM()
  
  # Baum-Welch
  warn=has_warning(bw <- train(hmm,seq,method = "BaumWelch",maxiter=100,logspace = FALSE,quiet=TRUE))
  if (warn==TRUE){print("Not converged")}
  
  # get viterbi path
  viterbi = Viterbi(bw,seq)
  viterbi.path=rownames(bw$E)[viterbi$path + 1]
  
  #get posterior highest state
  post.prob = posterior(bw,seq)
  post.path=tail(states,2)[apply(post.prob, 2, which.max)]
  return(list(viterbi.path,post.path))
}

vit.acc=c()
post.acc=c()

for (i in 1:length(tree_files)){
  
  pred=predict_tree(tree_files[i])
  v.pred<- pred[[1]];p.pred<-pred[[2]]
  print(table(v.pred))
  print(table(p.pred))
  struct_file=paste(strsplit(tree_files[i],"\\.")[[1]][1],".regions.txt",sep="")
  struct=read.table(paste("./regions.txt/",struct_file,sep=""))
  true=c(rep(struct$V1[1],5000),rep(struct$V1[2],5000),rep(struct$V1[3],5000),rep(struct$V1[4],5000))
  #true=c(rep(struct$V1[1],10000),rep(struct$V1[2],10000),rep(struct$V1[3],10000),rep(struct$V1[4],10000))
  true=as.character(true)
  true=paste("T",true,sep='')
  print(paste("Viterbi Accuracy:",tree_files[i],length(which(true==v.pred))/length(true)))
  print(paste("Posterior Accuracy:",tree_files[i],length(which(true==p.pred))/length(true)))
  vit.acc[i]=length(which(true==v.pred))/length(true)
  post.acc[i]=length(which(true==p.pred))/length(true)
}
