library("aphid") 
library("testit")

predict_tree <- function(file,model,variable){
  if(missing(model)) {
    model=4
  }
  if(missing(variable)) {
    variable="log.like.tree."
    #variable="post.prob.tree."
  }
  
  data=read.table(file, header=TRUE, sep = ",")
  numTrees=(ncol(data)-5)/2
  #print(numTrees)
  
  states <- c("Begin",paste("T",1:numTrees,sep=''))

  
  variables=paste(variable,as.character(1:numTrees),sep='')
  
  seq=letters[1:numTrees][apply(data[variables], 1, which.max)]  
 
  
  # Initial HMM
  #'   ### Define the transition probability matrix
  A <- matrix(c(rep(.01/(numTrees-1),(numTrees+1)**2)),nrow=numTrees+1,ncol=numTrees+1,byrow = TRUE)
  diag(A) <- .99
  A[1,]<-1/numTrees
  A[,1]<-0
  dimnames(A) <- list(from = states, to = states)
  
  #'   ### Define the emission probability matrix
  if (model==1){
    residues <- letters[1:numTrees]
    E=matrix(c(rep(.5/(numTrees-1),numTrees)),nrow=numTrees,ncol=numTrees,byrow = TRUE)
  }
  
  else if(model==2){
    residues <- letters[1:(numTrees+1)]
    seq[data$isConstant==1]=letters[(numTrees+1)]
    E=matrix(c(rep(.5/(numTrees),numTrees)),nrow=numTrees,ncol=numTrees+1,byrow = TRUE)  
  }
  
  else if(model==3){
    residues <- letters[1:(numTrees+2)]
    seq[data$isInformative==0]=letters[(numTrees+2)]
    seq[data$isConstant==1]=letters[(numTrees+1)]
    E=matrix(c(rep(.5/(numTrees+1),numTrees)),nrow=numTrees,ncol=numTrees+2,byrow = TRUE)
  }
  else{
    residues <- letters[1:(numTrees+3)]
    seq[data$sameParsimony==1]=letters[(numTrees+3)]
    seq[data$isInformative==0]=letters[(numTrees+2)]
    seq[data$isConstant==1]=letters[(numTrees+1)]
    E=matrix(c(rep(.5/(numTrees+2),numTrees)),nrow=numTrees,ncol=numTrees+3,byrow = TRUE)
  }
  diag(E) <- .5
  dimnames(E) <- list(states = states[-1], residues = residues)
  print(table(seq))
  
  #'   ### Build the HMM object
  hmm <- structure(list(A = A, E = E), class = "HMM")
  
  # Baum-Welch
  warn=has_warning(bw <- train(hmm,seq,method = "BaumWelch",maxiter=10000,logspace = FALSE,cores="autodetect",quiet=TRUE))
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

