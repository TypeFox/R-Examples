plotevent <-
function(last.well,first.ill,data)
{
subsetdata <- data[!is.na(data[,first.ill]),c(last.well,first.ill)]
subsetdata$n <- seq(1,dim(subsetdata)[1])

plot(c(subsetdata[,1], subsetdata[,2]),rep(0,2*nrow(subsetdata)), bty="n", yaxt="n",
     type="n",
     xlim=c(round(min(subsetdata[,1],subsetdata[,2],na.rm=T)),round(max(subsetdata[,1],subsetdata[,2],na.rm=T))),
     ylim=c(-2, max(subsetdata[,3])),
     xlab="Time",ylab="Conversions", 
     main=paste("Times between ",last.well," and ",first.ill,"",sep=""))
mtext(seq(0,nrow(subsetdata),5)[-1],side=2,at=seq(0,nrow(subsetdata),5)[-1],las=1,line=0)
mtext("Eq Cl",side=2,at=-2,las=1, padj=0,col="blue",font=2)
segments(subsetdata[,1],subsetdata[,3],subsetdata[,2],subsetdata[,3],lwd=1)

left <- unique(subsetdata[,1])
right <- unique(subsetdata[,2])
names(left) <- rep("0",length(left))
names(right) <- rep("1",length(right))
MM <- sort(c(left,right))
type <- as.numeric(names(MM))
type2 <- c(type[-1],0)
diff <- type-type2
int <- MM[diff<0|diff>0]
mat <- matrix(int,length(int)/2,2,byrow=TRUE)
d <- mat[,1] - mat[,2]
segments(mat[,1][d!=0],-2,mat[,2][d!=0],-2,lwd=3,col="blue")
}
