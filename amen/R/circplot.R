#' @title Circular network plot
#' 
#' @description Produce a circular network plot. 
#' 
#' @param Y (matrix) m by n relational matrix.
#' @param U (matrix) m by 2 matrix of row factors of Y.
#' @param V (matrix) n by 2 matrix of column factors of Y.
#' @param row.names (character vector) names of the row objects. 
#' @param col.names (character vector) names of the columns objects.
#' @param plotnames (logical) plot row and column names.
#' @param vscale (scalar) scaling factor for V coordinates.
#' @param pscale (scalar) scaling factor for plotting characters.
#' @param lcol (scalar or vector) line color(s) for the nonzero elements of Y.
#' @param rcol (scalar or vector) node color(s) for the rows.
#' @param ccol (scalar or vector) node color(s) for the columns.
#' @param pch (integer) plotting character.
#' @param lty (integer) line type.
#' @param jitter (scalar) a number to control jittering of nodes. 
#' @param bty (character) bounding box type.  
#' @param add (logical) add to existing plot 
#' 
#' @details 
#' This function creates a circle plot of a relational matrix or social network.
#' If not supplied via \code{U} and \code{V}, two-dimensional row factors and 
#' column factors are computed from the SVD of \code{Y}, scaled versions of 
#' which are used to plot positions on the outside edge (\code{U}) and inside
#' edge (\code{V}) of the circle plot. The magnitudes of the plotting characters
#' are determined by the magnitudes of the rows of \code{U} and \code{V}. 
#' Segments are drawn between each row object \code{i} and column object 
#' \code{j} for which \code{Y[i,j]!=0}. 
#'  
#' @return
#' NULL
#' 
#' @author Peter Hoff
#' 
#' @examples
#' data(IR90s) 
#' circplot(IR90s$dyadvars[,,1])
#' 
#' @export 
circplot<-function(Y,U=NULL,V=NULL,row.names=rownames(Y),col.names=colnames(Y),
                   plotnames=TRUE,vscale=.8,pscale=1.75,
                   lcol="gray",rcol="brown",ccol="blue",pch=16,lty=3,
                   jitter=.1*(nrow(Y)/(1+nrow(Y))) ,bty="n",add=FALSE )
{

  if(is.null(U))
  {
    a<-rowMeans(Y,na.rm=TRUE) ; b<-colMeans(Y,na.rm=TRUE)
    Y0<-Y ; Y0[is.na(Y)]<-(outer(a,b,"+"))[is.na(Y)] ; Y0<-Y0-mean(Y0)

    if(!all(Y==t(Y),na.rm=TRUE))
    {
      sY<-svd(Y0)
      u<-sY$u[,1:2] ; v<-sY$v[,1:2]
      mu<-sqrt( apply(u^2,1,sum) )
      mv<-sqrt( apply(v^2,1,sum) )
      u<-diag(1/mu)%*%u
      v<-diag(1/mv)%*%v*vscale
    }

    if( all(Y==t(Y),na.rm=TRUE) )
    {
      eY<-eigen(Y0)
      bv<-which(abs(eY$val)>= sort(abs(eY$val),decreasing=TRUE)[2])[1:2]
      u<-eY$vec[,bv]
      mu<-sqrt( apply(u^2,1,sum) )
      u<-diag(1/mu)%*%u
      mv<-mu ; v<-u
      ccol<-rcol
    }
  }

  if(!is.null(U))
  {
      if(is.null(V)){ V<-U ; ccol<-rcol ; vscale<-1 } 
      mu<-sqrt( apply(U^2,1,sum) )
      mv<-sqrt( apply(V^2,1,sum) )
      u<-diag(1/mu)%*%U
      v<-diag(1/mv)%*%V*vscale
  }


  ju<-1+jitter*(  rank(mu)/(nrow(Y)+1)  -.5 )
  u<-u*ju ; v<-v*ju


  rsum<-apply(abs(Y),1,sum,na.rm=TRUE)
  csum<-apply(abs(Y),2,sum,na.rm=TRUE)

  if(!add)
  {
    par(mfrow=c(1,1),mar=c(1,1,1,1) )
    plot(u*1.2,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty=bty) 
  }
  links<-which(Y!=0, arr.ind = TRUE)
  segments( u[links[,1],1],u[links[,1],2],
            v[links[,2],1],v[links[,2],2], col=lcol,lty=lty)

  if(plotnames)
  {
    if(is.null(row.names)){ row.names<-as.character(1:nrow(Y)) }
    if(is.null(col.names)){ col.names<-as.character(1:ncol(Y)) }
    text(u[rsum>0,] , row.names[rsum>0],cex=pscale*(mu[rsum>0])^.3,col=rcol)
    text(v[csum>0,] , col.names[csum>0],cex=pscale*(mv[csum>0])^.3,col=ccol)
  }

  if(!plotnames)
  {
    points(u[rsum>0,],cex=pscale*(mu[rsum>0])^.3,col=rcol,pch=pch)
    points(v[csum>0,],cex=pscale*(mv[csum>0])^.3,col=ccol,pch=pch)
  }
}





