library(HMM) 
library(ggplot2)
library(tidyverse)
library(ggpubr)

data=read.table("data5.fas.sitelh", header=TRUE, sep = ",")

states=c('a','b','c')

seq=states[apply(data[c("log.like.tree.1","log.like.tree.2")], 1, which.max)]
#seq=states[apply(data[c("post.prob.tree.1","post.prob.tree.2")], 1, which.max)]
seq=states[apply(data[c("post.prob.tree.1","post.prob.tree.2","post.prob.tree.3")], 1, which.max)]
table(seq)
table(seq[1:20000])
table(seq[20001:40000])

# Initial HMM
#hmm = initHMM(c("T1","T2"),c("a","b"),
#              transProbs=matrix(c(.9,.1,.1,.9),nrow=2,ncol=2,byrow = TRUE),
 #             emissionProbs=matrix(c(.9,.1,.1,.9),nrow=2,ncol=2,byrow = TRUE))

hmm = initHMM(c("T1","T2","T3"),c("a","b","c"),
              transProbs=matrix(c(.8,.1,.1,.1,.8,.1,.1,.1,.8),nrow=3,ncol=3,byrow = TRUE),
              emissionProbs=matrix(c(.8,.1,.1,.1,.8,.1,.1,.1,.8),nrow=3,ncol=3,byrow = TRUE))
print(hmm)

# Baum-Welch
bw = baumWelch(hmm,seq,30)
print(bw$hmm)

# get viterbi path
viterbi = viterbi(bw$hmm,seq)
table(viterbi)
table(viterbi[1:5000])
table(viterbi[5001:15000])
table(viterbi[15001:50000])

t1=ggplot(data, aes(x=site,y = post.prob.tree.1)) + 
  geom_bar(stat = "identity") +
  scale_fill_gradient(low='white',high='black')

t2=ggplot(data, aes(x=site,y = post.prob.tree.2)) + 
  geom_bar(stat = "identity") +
  scale_fill_gradient(low='white',high='black')

t3=ggplot(data, aes(x=site,y = post.prob.tree.3)) + 
  geom_bar(stat = "identity") +
  scale_fill_gradient(low='white',high='black')


df=data.frame(viterbi)
df$site <- seq.int(nrow(df))
df$path=as.numeric(factor(viterbi))

r1=ggplot(df, aes(x=site,y = path,color = cut(path, breaks = c(0, 1, 2,3))))+geom_col()+coord_cartesian(ylim=c(0,1))+
  scale_colour_manual(values = c('red','green','blue'),limits=c('(0,1]', '(1,2]','(2,3]'))+
  guides(color = guide_legend(title = 'Tree'))+labs(title="Output")

ggarrange(t1, t2, t3, r1,nrow=4)

df=data.frame(seq)
df$site <- seq.int(nrow(df))
df$path=as.numeric(factor(seq))

r2=ggplot(df, aes(x=site,y = path,color = cut(path, breaks = c(0, 1, 2,3))))+geom_col()+coord_cartesian(ylim=c(0,1))+
  scale_colour_manual(values = c('red','green','blue'),limits=c('(0,1]', '(1,2]','(2,3]'))+
  guides(color = guide_legend(title = 'Tree'))+
  labs(title="Input")

 


# check with correct path
true=c(rep("T1",5000),rep("T2",10000),rep("T3",35000))
length(which(true==check))

df=data.frame(true)
df$site <- seq.int(nrow(df))
df$path=as.numeric(factor(true))

r3=ggplot(df, aes(x=site,y = path,color = cut(path, breaks = c(0, 1, 2,3))))+geom_col()+coord_cartesian(ylim=c(0,1))+
  scale_colour_manual(values = c('red','green','blue'),limits=c('(0,1]', '(1,2]','(2,3]'))+
  guides(color = guide_legend(title = 'Tree'))+labs(title="True")

ggarrange(r2, r1, r3,nrow=3)
