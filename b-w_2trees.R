library(HMM) 


tree_files=list.files("./dataset",pattern = "*.sitelh")
#struct_files=list.files("./dataset",pattern = "*.txt")

predict_tree <- function(file){
  
  data=read.table(paste("./dataset/",file,sep=''), header=TRUE, sep = ",")
  
  states=c('a','b')
  seq=states[apply(data[c("post.prob.tree.1","post.prob.tree.2")], 1, which.max)]

  
  # Initial HMM
  hmm = initHMM(c("T1","T2"),c("a","b"),
              transProbs=matrix(c(.9,.1,.1,.9),nrow=2,ncol=2,byrow = TRUE),
              emissionProbs=matrix(c(.9,.1,.1,.9),nrow=2,ncol=2,byrow = TRUE))
  
  # Baum-Welch
  bw = baumWelch(hmm,seq,30)

  # get viterbi path
  viterbi = viterbi(bw$hmm,seq)
  return(viterbi)
}

accuracy=c()

for (i in 1:length(tree_files)){
  
  predicted=predict_tree(tree_files[i])
  struct_file=paste(strsplit(tree_files[i],"\\.")[[1]][1],".regions.txt",sep="")
  struct=read.table(paste("./dataset/",struct_file,sep=""))
  true=c(rep(struct$V1[1],5000),rep(struct$V1[2],5000),rep(struct$V1[3],5000),rep(struct$V1[4],5000))
  true=as.character(true)
  true=paste("T",true,sep='')
  accuracy[i]=length(which(true==predicted))/length(true)
}

