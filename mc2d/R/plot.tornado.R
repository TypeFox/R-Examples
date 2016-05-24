#<<BEGIN>>
plot.tornado <- function(x,which=1,name=NULL,stat=c("median","mean"),xlab="method",ylab="",...)
#TITLE Draws a Tornado chart.
#DESCRIPTION
# Draws a Tornado chart as provided by \samp{tornado}.
#KEYWORDS hplot
#INPUTS
#{x}<<A \code{\link{tornado}} object or a \code{\link{tornadounc}} object.>>
#[INPUTS]
#{which}<<Which output to print -for multivariates output-.>>
#{name}<<Vector of name of input variables.
#If NULL, the name will be given from the name of the elements.>>
#{stat}<<The name of the statistics of the output to be considered. For a \samp{tornado} object: "median" or "mean". 
#For a \samp{tornadounc} object: the value should match one row name of the \samp{tornadounc} object.
#Alternatively, for a \samp{tornadounc} object, the number of the row may be used.>>
#{xlab}<<Label of the x axis. if "method", use the correlation method used in the \samp{tornado} object.>>
#{ylab}<<Label of the y axis.>>
#{\dots}<<Further arguments to be passed to the \samp{plot} function.>>
#VALUE
# NULL
#DETAILS
#A point is drawn at the estimate
#and the segment reflects the uncertainty around this estimate.
#SEE ALSO
#\code{\link{tornado}}
#EXAMPLE

#CREATED 07-08-01
#EXAMPLE
#data(ec)
#x <- evalmcmod(ec$modEC2, nsv=100, nsu=100, seed=666)
#tor <- tornado(x,7)
#plot(tor)
#REVISED 10-02-10
#--------------------------------------------
#
{
  val <- x$value[[which]]
  if(is.null(val)) stop("Invalid value for which")
  nc <- ncol(val)
  nr <- nrow(val)
  if(!is.null(name)) {colnames(val) <- (rep(name,nc))[1:nc]}

  if(xlab=="method") xlab <- c("Spearman's rho statistic","Kendall's tau statistic","Pearson correlation")[pmatch(x$method,c("spearman","kendall","pearson"))]

	plot(c(-1.5,1),c(1,nc),type="n",axes= FALSE, xlab=xlab,ylab=ylab,...)
  axis(1,at=c(-1,-0.5,0,0.5,1))

  stat <- match.arg(stat)
  stat <- ifelse(stat=="mean" && nr!=1, 2 ,1)
  val <- val[,order(abs(val[stat,])),drop=FALSE]
	if(nr==1){
		segments(0,1:nc,val,1:nc,lwd=2)
    points(val,1:nc, pch="|", lwd=2)
    }

  else {
		segments(0,1:nc,val[stat,],1:nc,lwd=2,col="grey")
    points(val[stat,],1:nc, pch="|", lwd=2)
    if(nr>3){
        val <- apply(val[3:nr,],2,range)
        segments(val[1,],1:nc,val[2,],1:nc,lwd=2)
        points(val[1,],1:nc, pch="|", lwd=2)
        points(val[2,],1:nc, pch="|", lwd=2)
            }
    }

		abline(v=0)
		text(-1.4,1:nc,labels=paste(colnames(val),sep=":"),cex=.8)
		}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<BEGIN>>
plot.tornadounc <- function(x,which=1, stat="median", name=NULL, xlab="method", ylab="",...)
#ISALIAS plot.tornado
#--------------------------------------------
#
{
  statposs <- rownames(x$value[[which]])
  
  if(is.character(stat)) stat <- pmatch(stat, rownames(x$value[[which]]))
  if(is.na(stat)) stop("stat should match with: ",paste(statposs,collapse=", ")) 

  x$value <- list(x$value[[which]][stat,,drop=FALSE])
	plot.tornado(x,which=1, stat="median", name=name, xlab=xlab,ylab=ylab,...)
 }

