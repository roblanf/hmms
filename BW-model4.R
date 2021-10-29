library("aphid") 
library("testit")

tree_files=list.files("./tip7.spr1.sitelh.5k",pattern = "*.sitelh")
#struct_files=list.files("./dataset",pattern = "*.txt")

predict_tree <- function(file){
  
  data=read.table(paste("./tip7.spr1.sitelh.5k/",file,sep=''), header=TRUE, sep = ",")
  
  states <- c("Begin","T1", "T2")
  residues <- c("a","b","c","d","e")
  seq=head(residues,2)[apply(data[c("post.prob.tree.1","post.prob.tree.2")], 1, which.max)]
  seq[data$sameParsimony==1]='e'
  seq[data$isInformative==0]='d'
  seq[data$isConstant==1]='c'
  print(colSums(data != 0))
  print(table(seq))
  
  # Initial HMM
  #'   ### Define the transition probability matrix
  A <- matrix(c(0,0.5,0.5,0,.9,.1,0,.1,.9),nrow=3,ncol=3,byrow = TRUE)
  dimnames(A) <- list(from = states, to = states)
  #'   ### Define the emission probability matrix
  E=matrix(c(.7,.1,.05,.05,.1,.1,.7,.05,.05,.1),nrow=2,ncol=5,byrow = TRUE)
  dimnames(E) <- list(states = states[-1], residues = residues)
  #'   ### Build the HMM object
  hmm <- structure(list(A = A, E = E), class = "HMM")
  
  # Baum-Welch
  warn=has_warning(bw <- train(hmm,seq,method = "BaumWelch",maxiter=1000,logspace = FALSE,quiet=TRUE))
  if (warn==TRUE){
    print("Not converged, switching to Viterbi")
    #bw <- train(hmm,seq,method = "Viterbi",logspace = FALSE,quiet=TRUE)
  }
  
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

pred=predict_tree(paste("./tip7.spr1.sitelh.5k/",tree_files[1],sep=''))
v.pred<- pred[[1]];p.pred<-pred[[2]];df<-pred[[3]];seq<-pred[[4]];data=pred[[5]] 
true=c(rep("T2",5000),rep("T1",5000),rep("T2",5000),rep("T1",5000))

# mix graph
# df=data.frame(t(data.frame(c(df,data.frame(true)))))
# df <- tibble::rownames_to_column(df, "Array")
# df2<-df %>% 
#   melt(id.vars = "Array") %>%
#   mutate(variable = str_extract(variable, "[0-9]+")) %>%
#   mutate(value = case_when(
#     value == "a" ~ 1,
#     value == "b" ~ 2, 
#     TRUE ~ as.numeric(value)
#   )) %>%
#   mutate(variable = as.numeric(variable))
# 
# df2 %>% 
#   ggplot(aes(x = Array, y = variable, group = Array, fill = value)) +
#   geom_col() + coord_flip()


dl = data %>% pivot_longer(cols=c(post.prob.tree.1, post.prob.tree.2), names_to = "type", values_to = "measure")
sc=ggplot(dl, aes(x=site, y = measure, colour=type)) + geom_point() + geom_smooth(method='loess', span=0.03)


df=data.frame(true)
df$site <- seq.int(nrow(df))
df$path=as.numeric(factor(true))

r4=df %>% 
  ggplot(aes(x=site,y = path,col=path)) +
  geom_col() +coord_cartesian(ylim=c(0,1))

df=data.frame(v.pred)
df$site <- seq.int(nrow(df))
df$path=as.numeric(factor(v.pred))

r2=df %>% 
  ggplot(aes(x=site,y = path,col=path)) +
  geom_col() +coord_cartesian(ylim=c(0,1))


df=data.frame(seq)
df$site <- seq.int(nrow(df))
df$path=as.numeric(factor(seq))

r3=df %>% 
  ggplot(aes(x=site,y = path,col=path)) +
  geom_col() +coord_cartesian(ylim=c(0,1))

ggarrange(sc,r2, r3, r4,nrow=4)
