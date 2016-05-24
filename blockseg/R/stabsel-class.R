##' Class "stab.blockSeg"
##'
##' Class of object returned by the \code{stab.blockSeg} function.
##'
##' @section Slots: \describe{
##' 
##' \item{\code{RowBreaks}: }{a vectors of length the number of rows. Each case contains the number of active variable identified
##' along the stability selection.}  
##'
##' \item{\code{ColBreaks}: }{a vectors of length the number of columns. Each case contains the number of active variable identified
##' along the stability selection.}
##' 
##' }
##'
##' @section Methods:
##' Specific plotting and predict methods are available and documented
##' (\code{\link{plot,stab.blockSeg-method}}, \code{\link{evolution,stab.blockSeg-method}}.
##'
##' @importFrom stats predict residuals deviance
##' @import Matrix
##'
##' @aliases print,stab.blockSeg-method show,stab.blockSeg-method, 
##'
##' @docType class
##'
##' @keywords class
##'
##' @seealso See also \code{\link{plot,stab.blockSeg-method}}, \code{\link{evolution,stab.blockSeg-method}}
##' \code{\link{print,blockSeg-method}} and \code{\link{stab.blockSeg}}.
##'
##' @name stab.blockSeg-class
##' @rdname stab.blockSeg-class
##'
##' @exportClass stab.blockSeg
##' @exportMethod print
##' @exportMethod show
##'
setClass(
  Class="stab.blockSeg",
  representation=representation(
    RowBreaks="numeric",
    ColBreaks="numeric")
)
##' Plot method for a stab.blockSeg object
##'
##' Produce a plot of two-dimensional segmentation of a \code{stab.blockSeg} fit.
##'
##' @param x an object of class \code{stab.blockSeg}.
##' @param y the observations data (or a transformation).
##' @param threshold the threshold used (percent the maximum value).
##' @param postprocessing the condition if plot used a post-processing (if $post=TRUE) or not. 
##' If there is a post-processing, post-processing$adjacent is the maximal distance between two points.
##' @param col colours of the graphics. By default, it is "GrayLevel" to black and white colours. 
##' If it is another "character", it is a level blue or red. Else, it is possible to propose a sequence with 
##' the colour (rgb format).
##' @param ... used for S4 compatibility.
##'
##' @seealso \code{\linkS4class{stab.blockSeg}}.
##'
##' 
##' @name plot,stab.blockSeg-method
##' @aliases plot,stab.blockSeg-method
##' @aliases plot.stab.blockSeg
##' @docType methods
##' @rdname plot.stab.blockSeg
##' @seealso \code{\linkS4class{stab.blockSeg}}.
##' 
##' @examples
##' \dontrun{
##' n <- 100
##' ## model parameters 
##' K <- 5
##' mu <- suppressWarnings(matrix(rep(c(1,0),ceiling(K**2/2)), K,K))
##' Y <- rblockdata(n,mu,sigma=.5)$Y
##' stab.out <- stab.blockSeg(Y, 100, 15)
##' plot(stab.out,Y)
##' }
##' @exportMethod plot
setMethod(
  f="plot",
  signature="stab.blockSeg",
  definition=function(x,y,threshold=40,postprocessing=list(post=TRUE,adjacent=2),col="GrayLevel",...){
    
    if (!is.numeric(threshold)){
        stop("threshold must be a percent strictly between 0 and 100")
    } else if ((threshold<=0)||(threshold>=100)||(length(threshold)!=1)){
        stop("threshold must be a percent strictly between 0 and 100")
    }
    if (!is.list(postprocessing)){
        stop("postprocessing must be a list")
    } else {
        if (!is.logical(postprocessing$post)){
            stop("postprocessing$post must be logical")
        } else if (postprocessing$post){
            if (!is.numeric(postprocessing$adjacent)){
                stop("postprocessing$adjacent must be a positive integer")
        } else if ((postprocessing$adjacent<=0)||(floor(postprocessing$adjacent)!=postprocessing$adjacent)){
            stop("postprocessing$adjacent must be a positive integer")
        }
        }
    }
    if (!(is.matrix(y)||(class(y)=="dgeMatrix"))){
        stop("y must be the observations data (or a transformation)")
    }
    if (!is.character(col)){
        stop("col must be a character")
    } else {
        if (length(col)==1){
            if (col=="GrayLevel"){
                couleur=gray(seq(0,1, length=256))
                couleur=couleur[length(couleur):1]
            }else{
                couleur=c(rgb((0:201)/201*200,(0:201)/201*200,255,maxColorValue = 255),"white",
                          rgb(255,(0:201)/201*200,(0:201)/201*200,maxColorValue = 255)[201:0])
            }
        }else{
            couleur=col
        }
    }   
    brek=seq(min(y),max(y),length=length(couleur)+1)
    
    ValSeuilline=threshold/100*max(x@RowBreaks)
    ValSeuilCol=threshold/100*max(x@ColBreaks)
    ColLine=x@RowBreaks>ValSeuilline
    ColCol=x@ColBreaks>ValSeuilCol
    RowBreak=which(ColLine)
    ColBreak=which(ColCol)
    
    n=nrow(y)
    d=ncol(y)
    
    if (postprocessing$post){
        nr=length(RowBreak)
        if (nr>1){
            dist=RowBreak[2:nr]-RowBreak[1:(nr-1)]
            Rowbreakbis=c()
            ind=1
            Rowcompl=c()
            nr=nr-1
            while(ind<nr){
                if (dist[ind]<postprocessing$adjacent){
                    inddeb=ind
                    while((dist[ind]<postprocessing$adjacent)&(ind<nr)){
                        ind=ind+1
                    }
                    indfin=ind
                    indref=inddeb-1+which.max(x@RowBreaks[RowBreak[inddeb:indfin]])
                    Rowbreakbis=c(Rowbreakbis,RowBreak[indref])
                    suite=inddeb:indfin
                    Rowcompl=c(Rowcompl,RowBreak[suite[suite!=indref]])
                }else{
                    Rowbreakbis=c(Rowbreakbis,RowBreak[ind])
                }
                ind=ind+1
            }
            if (dist[nr]>=postprocessing$adjacent){
                Rowbreakbis=c(Rowbreakbis,RowBreak[nr+1])
            }else{
                Rowcompl=c(Rowcompl,RowBreak[nr+1])
            }
            RowBreak=Rowbreakbis
            ColLine[Rowcompl]=ColLine[Rowcompl]+2
        }
        nc=length(ColBreak)
        if (nc>1){
            dist=ColBreak[2:nc]-ColBreak[1:(nc-1)]
            Colbreakbis=c()
            ind=1
            Colcompl=c()
            nc=nc-1
            while(ind<nc){
                if (dist[ind]<postprocessing$adjacent){
                    inddeb=ind
                    while((dist[ind]<postprocessing$adjacent)&(ind<nc)){
                        ind=ind+1
                    }
                    indfin=ind
                    indref=inddeb-1+which.max(x@ColBreaks[ColBreak[inddeb:indfin]])
                    Colbreakbis=c(Colbreakbis,ColBreak[indref])
                    suite=inddeb:indfin
                    Colcompl=c(Colcompl,ColBreak[suite[suite!=indref]])
                }else{
                    Colbreakbis=c(Colbreakbis,ColBreak[ind])
          }
                ind=ind+1
            }
            if (dist[nc]>=postprocessing$adjacent){
          Colbreakbis=c(Colbreakbis,ColBreak[nc+1])
            }else{
                Colcompl=c(Colcompl,ColBreak[nc+1])
            }
            ColBreak=Colbreakbis
            ColCol[Colcompl]=ColCol[Colcompl]+2
      }
    }
    
    #### Plot a proprement parle
    
    par(mfrow=c(2,2),oma=c(0,0,3,0))
    
    ### Affichage matrice originale
    
    image(1:d,1:n,t(y)[,n:1],xlab="",ylab="",xaxt="n",yaxt="n",main="Original data",col=couleur,breaks=brek)
    abline(v=ColBreak-0.5,col="purple")
    abline(h=n-RowBreak+1.5,col="purple")
    
    ## Affichage plot ligne
    
    plot(x@RowBreaks,n:1,col=ColLine+1,xlab="",ylab="n", ylim=c(1,n),yaxt="n")
    if (length(RowBreak)<=1){
        title(main=paste(as.character(length(RowBreak))," row break",sep=""))
    }else{
        title(main=paste(as.character(length(RowBreak))," row breaks",sep=""))
    }
    abline(v=ValSeuilline,col="purple")
    nax=floor(6*par("din")[2]/4)
    axis(2,at=floor(seq(0,n,length=nax)),labels=floor(seq(n,0,length=nax)))
    
    
    ## Affichage plot colonne
    
    plot(1:d,x@ColBreaks,col=ColCol+1,xlab="n", ylab="")
    if (length(ColBreak)<=1){
        title(main=paste(as.character(length(ColBreak))," column break",sep=""))
    }else{
        title(main=paste(as.character(length(ColBreak))," column breaks",sep=""))
    }
    abline(h=ValSeuilCol,col="purple")
    
    ## Affichage matrice resumee
    emplz=diff(unique(c(1,RowBreak,n+1)))
    z=bdiag(lapply(emplz,function(i) rep(1,i) ))
    nz=length(emplz)
    emplw=diff(unique(c(1,ColBreak,d+1)))
    w=bdiag(lapply(emplw,function(i) rep(1,i)))
    nw=length(emplw)
    resum=as.matrix(t(z)%*%y%*%w/(emplz%*%t(emplw)))
    if (is.matrix(resum)){
        if (nz!=1){
            image(x=unique(c(1,ColBreak,d)),y=unique(n-c(n,RowBreak[length(RowBreak):1],1)+1),
                  z=t(resum[nz:1,]),
                  xlab="",ylab="",xaxt="n",yaxt="n",main="Summarized data",col=couleur,breaks=brek)
        }else{
            image(x=unique(c(1,ColBreak,d)),y=unique(n-c(n,RowBreak[RowBreak:1],1)+1),
                  z=t(t(t(resum))),
                  xlab="",ylab="",xaxt="n",yaxt="n",main="Summarized data",col=couleur,breaks=brek)
        }
        abline(v=ColBreak,col="purple")
        abline(h=n-RowBreak+1,col="purple")
    }
    title(outer=TRUE,main=paste(as.character(threshold),"%",sep=""))
})

##' Plot method for a stab.blockSeg object
##'
##' Produce a plot of two-dimensional segmentation of a \code{stab.blockSeg} fit.
##'
##' @param x an object of class \code{stab.blockSeg}.
##' @param y the observations data (or a transformation).
##' @param thresholds the thresholds used (percent the maximum value). By default, thresholds = 10 * (8:1).
##' @param postprocessing the condition if plot used a post-processing (if $post=TRUE) or not. 
##' If there is a post-processing, post-processing$adjacent is the maximal distance between two points.
##' @param col colours of the graphics. By default, it is "GrayLevel" to black and white colours. 
##' If it is another "character", it is a level blue or red. Else, it is possible to propose a sequence with 
##' the colour (rgb format).
##' @param ask If \code{TRUE}, to hit will be necessary to see next plot.
##' @param ... used for S4 compatibility.
##' 
##' @name evolution
##' @aliases evolution
##' @docType methods
##' @rdname evolution
##' @seealso \code{\linkS4class{stab.blockSeg}}.
##' 
##' @examples
##' n <- 100
##' ## model parameters 
##' K <- 5
##' mu <- suppressWarnings(matrix(rep(c(1,0),ceiling(K**2/2)), K,K))
##' Y <- rblockdata(n,mu,sigma=.5)$Y
##' stab.out <- stab.blockSeg(Y, 100, 15)
##' evolution(stab.out,Y)
##' 
##' @exportMethod evolution
setGeneric("evolution",function(x,y,thresholds=10*(8:1),postprocessing=list(post=TRUE,adjacent=2),col="GrayLevel",ask=TRUE){standardGeneric("evolution")})

##' @rdname evolution
setMethod("evolution", "stab.blockSeg",
          definition=function(x,y,thresholds=10*(8:1),postprocessing=list(post=TRUE,adjacent=2),col="GrayLevel",ask=TRUE){
    
    if (!is.numeric(thresholds)){
      stop("thresholds must be a sequence percent strictly between 0 and 100")
    } else if (any(thresholds<=0)||any(thresholds>=100)){
      stop("thresholds must be a percent strictly between 0 and 100")
    } else if (length(thresholds)==1){
      warning("plot is used")
      plot(x=x,y=y,threshold=thresholds,postprocessing=postprocessing,col=col)
      opt=options(show.error.messages=FALSE)
      on.exit(opt) 
      stop() 
    }else{
      thresholds=unique(thresholds)
    }
    if (!is.list(postprocessing)){
      stop("postprocessing must be a list")
    }else{
      if (!is.logical(postprocessing$post)){
        stop("postprocessing$post must be logical")
      } else if (postprocessing$post){
        if (!is.numeric(postprocessing$adjacent)){
          stop("postprocessing$adjacent must be a positive integer")
        }else if ((postprocessing$adjacent<=0)||(floor(postprocessing$adjacent)!=postprocessing$adjacent)){
          stop("postprocessing$adjacent must be a positive integer")
        }
      }
    }
    if (!(is.matrix(y)||(class(y)=="dgeMatrix"))){
      stop("y must be the observations data (or a transformation)")
    }
    if (!is.character(col)){
      stop("col must be a character")
    }else{
      if (length(col)==1){
        if (col=="GrayLevel"){
          couleur=gray(seq(0,1, length=256))
          couleur=couleur[length(couleur):1]
        }else{
          couleur=c(rgb((0:201)/201*200,(0:201)/201*200,255,maxColorValue = 255),"white",
                    rgb(255,(0:201)/201*200,(0:201)/201*200,maxColorValue = 255)[201:0])
        }
      }else{
        couleur=col
      }
    }
    if (!is.logical(ask)){
      stop("ask must be logical")
    }
    
    brek=seq(min(y),max(y),length=length(couleur)+1)
    
    n=nrow(y)
    d=ncol(y)
    
    par(mfrow=c(2,2),oma=c(0,0,3,0))
    f=function(i){
      rep(1,i)
    }
    
    for (threshold in thresholds){
      
      ValSeuilline=threshold/100*max(x@RowBreaks)
      ValSeuilCol=threshold/100*max(x@ColBreaks)
      ColLine=x@RowBreaks>ValSeuilline
      ColCol=x@ColBreaks>ValSeuilCol
      RowBreak=which(ColLine)
      ColBreak=which(ColCol)
      
      if (postprocessing$post){
        nr=length(RowBreak)
        if (nr>1){
          dist=RowBreak[2:nr]-RowBreak[1:(nr-1)]
          Rowbreakbis=c()
          ind=1
          Rowcompl=c()
          nr=nr-1
          while(ind<=nr){
            if (dist[ind]<postprocessing$adjacent){
              inddeb=ind
              while((dist[ind]<postprocessing$adjacent)&(ind<nr)){
                ind=ind+1
              }
              indfin=ind
              indref=inddeb-1+which.max(x@RowBreaks[RowBreak[inddeb:indfin]])
              Rowbreakbis=c(Rowbreakbis,RowBreak[indref])
              suite=inddeb:indfin
              Rowcompl=c(Rowcompl,RowBreak[suite[suite!=indref]])
            }else{
              Rowbreakbis=c(Rowbreakbis,RowBreak[ind])
            }
            ind=ind+1
          }
          if (dist[nr]>=postprocessing$adjacent){
            Rowbreakbis=c(Rowbreakbis,RowBreak[nr+1])
          }else{
            Rowcompl=c(Rowcompl,RowBreak[nr+1])
          }
          RowBreak=Rowbreakbis
          ColLine[Rowcompl]=ColLine[Rowcompl]+2
        }
        nc=length(ColBreak)
        if (nc>1){
          dist=ColBreak[2:nc]-ColBreak[1:(nc-1)]
          Colbreakbis=c()
          ind=1
          Colcompl=c()
          nc=nc-1
          while(ind<=nc){
            if (dist[ind]<postprocessing$adjacent){
              inddeb=ind
              while((dist[ind]<postprocessing$adjacent)&(ind<nc)){
                ind=ind+1
              }
              indfin=ind
              indref=inddeb-1+which.max(x@ColBreaks[ColBreak[inddeb:indfin]])
              Colbreakbis=c(Colbreakbis,ColBreak[indref])
              suite=inddeb:indfin
              Colcompl=c(Colcompl,ColBreak[suite[suite!=indref]])
            }else{
              Colbreakbis=c(Colbreakbis,ColBreak[ind])
            }
            ind=ind+1
          }
          if (dist[nc]>=postprocessing$adjacent){
            Colbreakbis=c(Colbreakbis,ColBreak[nc+1])
          }else{
            Colcompl=c(Colcompl,ColBreak[nc+1])
          }
          ColBreak=Colbreakbis
          ColCol[Colcompl]=ColCol[Colcompl]+2
        }
      }
      
      #### Plot a proprement parle
      
      ### Affichage matrice originale
      
      image(1:d,1:n,t(y)[,n:1],xlab="",ylab="",xaxt="n",yaxt="n",main="Original data",col=couleur,breaks=brek)
      abline(v=ColBreak-0.5,col="purple")
      abline(h=n-RowBreak+1.5,col="purple")
      
      ### Affichage plot ligne
      
      plot(x@RowBreaks,n:1,col=ColLine+1,xlab="",ylab="n", ylim=c(1,n),yaxt="n")
      if (length(RowBreak)<=1){
        title(main=paste(as.character(length(RowBreak))," row break",sep=""))
      }else{
        title(main=paste(as.character(length(RowBreak))," row breaks",sep=""))
      }
      abline(v=ValSeuilline,col="purple")
      nax=floor(6*par("din")[2]/4)
      axis(2,at=floor(seq(0,n,length=nax)),labels=floor(seq(n,0,length=nax)))
      
      
      ### Affichage plot colonne
      
      plot(1:d,x@ColBreaks,col=ColCol+1,xlab="n", ylab="")
      if (length(ColBreak)<=1){
        title(main=paste(as.character(length(ColBreak))," column break",sep=""))
      }else{
        title(main=paste(as.character(length(ColBreak))," column breaks",sep=""))
      }
      abline(h=ValSeuilCol,col="purple")
      
      ### Affichage matrice resumee
      emplz=diff(unique(c(1,RowBreak,n+1)))
      z=bdiag(lapply(emplz,f))
      nz=length(emplz)
      emplw=diff(unique(c(1,ColBreak,d+1)))
      w=bdiag(lapply(emplw,f))
      nw=length(emplw)
      resum=as.matrix(t(z)%*%y%*%w/(emplz%*%t(emplw)))
      if (is.matrix(resum)){
        if (nz!=1){
          image(x=unique(c(1,ColBreak,d)),y=unique(n-c(n,RowBreak[length(RowBreak):1],1)+1),
                z=t(resum[nz:1,]),
                xlab="",ylab="",xaxt="n",yaxt="n",main="Summarized data",col=couleur,breaks=brek)
        }else{
          image(x=unique(c(1,ColBreak,d)),y=unique(n-c(n,RowBreak[RowBreak:1],1)+1),
                z=t(t(t(resum))),
                xlab="",ylab="",xaxt="n",yaxt="n",main="Summarized data",col=couleur,breaks=brek)
        }
        abline(v=ColBreak,col="purple")
        abline(h=n-RowBreak+1,col="purple")
      }
      title(outer=TRUE,main=paste(as.character(threshold),"%",sep=""))
      if (threshold==thresholds[1]){
        par(mfrow=c(2,2),oma=c(0,0,3,0),ask=ask)
      }
    }
    par(ask=FALSE)
  }
)
