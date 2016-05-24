sim2introgress <-
function(x){
options(warn=-1)
parental1<-t(x[[1]])
parental2<-t(x[[2]])
n<-nrow(t(x[[1]]))
admix.gen <-matrix(nrow=n)
if (is.na(x[[3]]) == F){admix.gen<-cbind(admix.gen,t(x[[3]]))}
if (is.na(x[[4]]) == F){admix.gen<-cbind(admix.gen,t(x[[4]]))}
if (is.na(x[[5]]) == F){admix.gen<-cbind(admix.gen,t(x[[5]]))}
if (is.na(x[[6]]) == F){admix.gen<-cbind(admix.gen,t(x[[6]]))}
admix.gen<-admix.gen[,which(colMeans(is.na(admix.gen)) < 1)]

loci<-cbind(colnames(x[[1]]),rep("D",ncol(x[[1]])))
colnames(loci)<-c("locus","type")
introgress<-list(parental1,parental2,admix.gen,loci)
names(introgress)<-c("parental1","parental2","admix.gen","loci")
introgress

}
