plotclinical<-function(DATA, nam=names(DATA[[1]]), col=c("darkolivegreen3","darkorange1","lightsalmon3", "khaki1"))
{
N<-rep(NA,length(DATA)-1)
for (i in 1:(length(DATA)-1)) {N[i]<-(length(names(DATA[[i]]))==length(names(DATA[[i+1]])) && all(names(DATA[[i]])==names(DATA[[i+1]])))}

if (all(N))
{
 wid<-c(1,1,rep(2,length(nam)))
 hei<-c(rep(2,length(DATA)),2)
 #windows(height=length(DATA)+1, width=sum(wid))
 lay.mt<-matrix(c(1:((length(DATA)+1)*(length(nam)+2))), nrow=length(DATA)+1, ncol=length(nam)+2, byrow=T)
 lay<-layout(lay.mt,wid,hei)
 #layout.show(lay)
 par(mai=c(0.2,0,0.1,0.1))

 for (i in 1:(length(DATA)))
 {
  plot(1, type="n", axes=F, xlab="", ylab="")
  text(x=1,y=1,labels=names(DATA)[i], cex=2)
  plot(1, type="n", axes=F, xlab="", ylab="")
  text(x=1,y=1,labels=length(DATA[[i]][,1]), cex=2)

 for (j in nam)
  {
 if (all(is.na(DATA[[i]][,j])))  #chybajuci parameter ako vektor NA s rovnakym nazvom
      {plot(1, type="n", axes=F, xlab="", ylab="")} else {
  	if (is.factor(DATA[[i]][,j]) || is.character(DATA[[i]][,j]) )
      {temp<-as.matrix(table(as.factor(DATA[[i]][,j]),useNA="ifany"))
 	#rownames(temp)<-rownames(table(DATA[[i]][,j]))
 	barplot(temp/sum(temp)*100, beside=FALSE, col=col, horiz=TRUE, axes=F)
 	}
 	if (is.numeric(DATA[[i]][,j]) && i==length(DATA)) 
 	{boxplot(DATA[[i]][,j], frame.plot=FALSE, horizontal=TRUE, col="slateblue3", 
        ylim=c(min(DATA[[length(DATA)]][,j],na.rm=T),max(DATA[[length(DATA)]][,j],na.rm=T)) )}
	} 
      if (is.numeric(DATA[[i]][,j]) && i!=length(DATA)) {boxplot(DATA[[i]][,j], frame.plot=FALSE, axes=FALSE, horizontal=TRUE, col="slateblue3",  
        ylim=c(min(DATA[[length(DATA)]][,j],na.rm=T),max(DATA[[length(DATA)]][,j],na.rm=T)) )}
 }
 }
plot(1, type="n", axes=F, xlab="", ylab="")
plot(1, type="n", axes=F, xlab="", ylab="")
for (j in names(DATA[[i]]))
{
 plot(1, type="n", axes=F, xlab="", ylab="")
 if (is.factor(DATA[[i]][,j]) || is.character(DATA[[i]][,j]))
  {temp<-as.matrix(table(as.factor(DATA[[i]][,j]),useNA="ifany"))
   rownames(temp)[which(is.na(rownames(temp)))]<-"NA"
   legend(x="top", legend=rownames(temp),fill=col, bty="n", horiz=TRUE)}
   text(x=1, y=1, labels=j, cex=2)
}
  } else {stop("Names in list objects differ\n")}
}
