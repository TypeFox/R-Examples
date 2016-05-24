#' Network plotting
#' 
#' Plot the graph of a sociomatrix 
#' 
#' @usage netplot(Y,X=NULL,xaxt="n",yaxt="n",xlab="",ylab="",
#'  lcol="gray",ncol="black",lwd=1,lty=1,pch=16,bty="n",plotnames=FALSE,
#'  seed=1,
#'  plot.iso=TRUE,directed=NULL,add=FALSE,...)
#' @param Y a sociomatrix 
#' @param X coordinates for plotting the nodes
#' @param xaxt x-axis type
#' @param yaxt y-axis type
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param lcol edge color
#' @param ncol node color (can be node-specific)
#' @param lwd  line width
#' @param lty  line type 
#' @param pch plotting character for nodes (can be node-specific) 
#' @param bty bounding box type
#' @param plotnames plot rownames of Y as node labels 
#' @param seed random seed
#' @param plot.iso include isolates in plot 
#' @param directed draw arrows
#' @param add add to an existing plot region
#' @param \ldots additional plotting parameters
#' @author Peter Hoff
#' @examples
#' data(addhealthc3)
#' Y<-addhealthc3$Y
#' X<-xnet(Y) 
#' netplot(Y,X) 
#' 
#' @export netplot
netplot<-function(Y,X=NULL,xaxt="n",yaxt="n",xlab="",ylab="",
                  lcol="gray",ncol="black",lwd=1,lty=1,pch=16,
                  bty="n",plotnames=FALSE,seed=1,plot.iso=TRUE,
                  directed=NULL,add=FALSE,...)
{
  if(nrow(Y)!=ncol(Y))
  {
    if(is.null(directed)){Y0<-el2sm(Y);directed<-(sum(Y0*t(Y0),na.rm=TRUE)!=0)}
    Y<-el2sm(Y,directed)
  }

  if(is.null(X)) { X<-xnet(Y,seed=seed) }
  if(is.null(directed)){ directed<-!all(Y==t(Y),na.rm=TRUE) }

  if(!plot.iso)
  {
    deg<-apply(Y,1,sum,na.rm=TRUE) + apply(Y,2,sum,na.rm=TRUE)
    Y<-Y[deg>0,deg>0]
    X<-X[deg>0,]
  }

  if(!add){plot(X,type="n",xaxt=xaxt,yaxt=yaxt,xlab=xlab,ylab=ylab,bty=bty,...)}

  if(!directed) {addlines(Y,X,lwd=lwd,lty=lty,col=lcol,...)}
  if(directed){addlines(Y,X,lwd=lwd,lty=lty,col=lcol,alength=.15,...) }

  if(!plotnames) {points(X[,1],X[,2],col=ncol,pch=pch,...) }
  if(plotnames)
  {
    if(is.null(rownames(Y))) { rownames(Y)<-1:nrow(Y)  }
    text(X[,1],X[,2],rownames(Y),col=ncol,...)
   }

}





