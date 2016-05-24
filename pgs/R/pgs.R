### Source file for PGS package
### Version: 0.3-0
### Authors: Kien Kieu & Marianne Mora
### License: CeCILL (see COPYING file) 
### Copyright © 2010 INRA, Université Paris Ouest Nanterre la Défense

#######################################################################
#                        Generic functions                            #
#######################################################################

setGeneric("print")
setGeneric("plot")
setGeneric("scaling",def=function(x,s){standardGeneric("scaling")})
setGeneric("ltransf",def=function(x,m){standardGeneric("ltransf")})
setGeneric("content",def=function(x){standardGeneric("content")})
setGeneric("covariogram",def=function(x,f,sym=FALSE){standardGeneric("covariogram")})
setGeneric("dataCovariogram",def=function(x,coord){standardGeneric("dataCovariogram")})

#######################################################################
# Definition of the class "Figure" and derived classes (only bounded  #
# figures).                                                           #
#######################################################################

setClass("Figure",representation(dimspace="numeric",coord="matrix"))

## Scaling 
setMethod("scaling",signature(x="Figure",s="numeric"),
          function(x,s){
            x@coord<-x@coord*s
            x
          })

## Linear transform
setMethod("ltransf",signature(x="Figure",m="matrix"),
          function(x,m) {
            x@coord<-m%*%x@coord
            if(nrow(m)!=x@dimspace) x@dimspace <- nrow(m)
            x
          })

## Derived classes

setClass("Quadrat",representation("Figure"))
setClass("Segment",representation("Figure"))
setClass("PointPattern",representation("Figure"))

## Generators.
Quadrat<-function(hsize,vsize=hsize){ #Planar quadrat (rectangle)
  new("Quadrat",dimspace=2,coord=matrix(c(hsize,vsize),2,1))
}
Segment<-function(end){
  end <- as.matrix(end)
  new("Segment",dimspace=nrow(end),coord=as.matrix(end))
}
PointPattern<-function(coord){
  if(is.vector(coord)) coord<-as.matrix(coord)
  new("PointPattern",dimspace=nrow(coord),coord=coord)
}
PP2 <- function(n=4,h=1,v=h,h3=TRUE) {
  coord <- switch(n,
                  rep(0,2),
                  NaN,
                  NaN,
                  matrix(c(0,0, 1,0, 1,1, 0,1),2,4),
                  matrix(c(0,0, 1,0, 1,1, 0,1, .5,.5),2,5),
                  matrix(c(0,0, 0.5,0, 1,0,
                           0,1, 0.5,1, 1,1),2,6),
                  matrix(c(0.25,0, 0.75,0,
                           0,0.5, 0.5,0.5, 1,0.5,
                           0.25,1, 0.75,1),2,7),
                  matrix(c(0,0, 0.5,0, 1,0,
                           0.25,0.5, 0.75,0.5,
                           0,1, 0.5,1, 1,1),2,8),
                  matrix(c(0,0, 0.5,0, 1,0, 0,0.5,
                           0.5,0.5, 1,0.5,
                           0,1, 0.5,1, 1,1),2,9))
  if(any(is.null(coord)||is.nan(coord))) stop("Argument n must be in {1,4,5,6,7,8,9}")
  if(!h3) coord <- coord[2:1,]
  coord <- diag(c(h,v))%*%coord
  PointPattern(coord)
}

PP3 <- function(n=4,xp=1,yp=xp,h3=TRUE) {
  pp2 <- PP2(n,xp,yp,h3)
  coord <- rbind(pp2@coord,0)
  PointPattern(coord)
}

## Plot functions
setMethod("plot",signature(x="Quadrat",y="missing"),
          function(x,y,add=FALSE,origin=c(0,0),xlab="",ylab="",...) {
            if(add) {
              lines(origin[1]+c(0,x@coord[1],x@coord[1],0,0),origin[2]+c(0,0,x@coord[2],x@coord[2],0),...)
            } else {
              plot(origin[1]+c(0,x@coord[1],x@coord[1],0,0),origin[2]+c(0,0,x@coord[2],x@coord[2],0),type="l",
                 xlab=xlab,ylab=ylab,...)
            }
          })
setMethod("plot",signature(x="Segment",y="missing"),
          function(x,y,add=FALSE,origin=rep(0,x@dimspace),xlab="",ylab="",...) {
            if(x@dimspace==2) {
              if(add) {
                lines(origin[1]+c(0,x@coord[1]),origin[2]+c(0,x@coord[2]),...)
              } else {
                plot(origin[1]+c(0,x@coord[1]),origin[2]+c(0,x@coord[2]),type="l",xlab=xlab,ylab=ylab,...)
              }
            } else {
              stop(paste("Segment plot in ",x@dimspace,"D not implemented",sep=""))
            }
          }
         )
setMethod("plot",signature(x="Segment",y="matrix"),
          function(x,y,add=FALSE,origin=rep(0,x@dimspace),xlab="",ylab="",density=NULL,
                    angle=45,border=NULL,col=NA,lty=NULL,...) {
            if(x@dimspace==2) {
              if(!is.matrix(origin))
                origin <- as.matrix(origin)
              if(!add)
                plot(x=origin[1]+c(0,x@coord[1]),y=origin[2]+c(0,x@coord[2]),xlab=xlab,ylab=ylab,type="n",...)
              xlim <- par("usr")[1:2]
              ylim <- par("usr")[3:4]
              rect <- matrix(c(xlim[1],ylim[1],xlim[2],ylim[1],xlim[2],ylim[2],
                               xlim[1],ylim[2]),2,4)
              Ei <- svd(cbind(y,rep(0,2)),nv=0)$u
              bb <- t(apply(Ei%*%rect,1,range))
              ci <- Ei%*%(cbind(rep(0,2),x@coord)+kronecker(origin,matrix(1,1,2)))
              xi <- bb[1,]
              yi <- ci[2,]
              xyi <- rbind(rep(xi,length(yi)),rep(yi,rep(2,length(yi))))
              xy <- solve(Ei)%*%xyi
              xy[,3:4] <- xy[,4:3]
              polygon(xy[1,],xy[2,],density=density,angle=angle,border=border,col=col,lty=lty)
            } else {
            stop(paste("Segment plot in ",x@dimspace,"D not implemented",sep=""))
          }
        }
)
setMethod("plot",signature(x="PointPattern",y="missing"),
          function(x,y,add=FALSE,origin=rep(0,x@dimspace),xlab="",ylab="",...) {
            if(x@dimspace==2) {
              if(add) {
                points(x=origin[1]+x@coord[1,],y=origin[2]+x@coord[2,],...)
              } else {
                plot(x=origin[1]+x@coord[1,],y=origin[2]+x@coord[2,],xlab=xlab,ylab=ylab,...)
              }
            } else {
          stop(paste("Point Pattern plot in ",x@dimspace,"D not implemented",sep=""))
            }
          }
         )
setMethod("plot",signature(x="PointPattern",y="matrix"),
          function(x,y,add=FALSE,origin=rep(0,x@dimspace),xlab="",ylab="",...) {
            if(x@dimspace==2) {
              if(!is.matrix(origin))
                origin <- as.matrix(origin)
              if(!add)
                plot(x=origin[1]+x@coord[1,],y=origin[2]+x@coord[2,],xlab=xlab,ylab=ylab,type="n",...)
              xlim <- par("usr")[1:2]
              ylim <- par("usr")[3:4]
              rect <- matrix(c(xlim[1],ylim[1],xlim[2],ylim[1],xlim[2],ylim[2],
                               xlim[1],ylim[2]),2,4)
              Ei <- svd(cbind(y,rep(0,2)),nv=0)$u
              bb <- t(apply(Ei%*%rect,1,range))
              ci <- Ei%*%(x@coord+kronecker(origin,matrix(1,1,ncol(x@coord))))
              xi <- bb[1,]
              yi <- ci[2,]
              xyi <- rbind(rep(xi,length(yi)),rep(yi,rep(2,length(yi))))
              xy <- solve(Ei)%*%xyi
              for(i in seq(1,ncol(xy/2),by=2)) {
                lines(xy[1,c(i,i+1)],xy[2,c(i,i+1)],...)
              }
            } else {
            stop(paste("Point Pattern plot in ",x@dimspace,"D not implemented",sep=""))
          }
        }
)

## Print functions
setMethod("print","Quadrat",
          function(x,...) {
            cat("Quadrat\n")
            cat(paste("hsize=",x@coord[1],", vsize=",x@coord[2],"\n",sep=""))
          })
setMethod("print","Segment",
          function(x,...) {
            cat("Segment\n")
            cat("end:\n")
            print(x@coord)
            cat(paste("length=",format(content(x)),"\n",sep=""))
          })
setMethod("print","PointPattern",
          function(x,...) {
            cat("PointPattern\n")
            cat("coord:\n")
            print(x@coord)
          })

## Content functions
setMethod("content","Quadrat",
          function(x) {
            prod(x@coord)
          })
setMethod("content","Segment",
          function(x) {
            sqrt(sum(x@coord^2))
          })
setMethod("content","PointPattern",
          function(x) {
            ncol(x@coord)
          })

## Covariograms of the figures.
setMethod("covariogram",signature(x="Quadrat",f="function"),
          function(x,f,sym=FALSE){
            ucoord <- x@coord
            if(sym){
              integrand<-function(x){
                prod(ucoord-abs(x))*(f(x)+f(x*c(-1,1)))
              }
            } else {
              integrand<-function(x){
                prod(ucoord-abs(x))*(f(x)+f(x*c(-1,1))+f(x*c(1,-1))+f(-x))
              }
            }
            ifelse(sym,2,1)*cuhre(ndim=2,ncomp=1,integrand=integrand,
                                  flags=list(verbose=0),
                                  upper=ucoord,rel.tol=1e-2)$value
          })
setMethod("covariogram",signature(x="Segment",f="function"),
          function(x,f,sym=FALSE){
            l<-content(x)
            if(sym){
              integrand<-function(t){
                h<-outer(as.vector(x@coord)/l,t)
                (l-t)*f(h)
              }
            } else {
              integrand<-function(t){
                h<-outer(as.vector(x@coord)/l,t)
                (l-t)*(f(h)+f(-h))
              }
            }
            ifelse(sym,2,1)* integrate(integrand,0,l)$value
          })
setMethod("covariogram",signature(x="PointPattern",f="function"),
          function(x,f,sym=FALSE){
            n<-dim(x@coord)[2] # number of points
            g0 <- pointcov(x@coord)
            if(sym)
              sum(f(g0$ud)*g0$n)
            else
              f(g0$ud[,1,drop=FALSE])*g0$n[1]+sum(f(g0$ud[,-1,drop=FALSE])*g0$n[-1])/2+
                sum(f(-g0$ud[,-1,drop=FALSE])*g0$n[-1])/2
          })

#######################################################################
#                  Definition of class VecLat                         #
#######################################################################

setClass("VecLat",representation(dimspace="numeric",dimsupp="numeric",
                                 gmat="matrix",gmat0="matrix",det="numeric"))
## Generators
VecLat <- function(gmat) {
  gmat <- as.matrix(gmat)
  dimsupp <- ncol(gmat)
  detLat <- abs(prod(svd(gmat,0,0)$d))
  ## Check that the generating matrix is not singular. The criterion may be too
  ## strong.
  if(identical(all.equal(detLat^(1/nrow(gmat)),0,tol=.Machine$double.eps^0.75),TRUE))
    stop("The matrix gmat seems to be singular")
  gmat0 <- detLat^(-1/dimsupp)*gmat
  new("VecLat",dimspace=nrow(gmat),dimsupp=dimsupp,gmat=gmat,gmat0=gmat0,det=detLat)
}
RectLat2 <- function(h=1,v=h) {
  if(min(h,v)<=0) stop("Both h and v must be positive")
  VecLat(diag(c(h,v)))
}
HexLat2 <- function(det=1) {
  VecLat(sqrt(det)*matrix(sqrt(c(2/sqrt(3),0,1/(2*sqrt(3)),sqrt(3)/2)),2,2))
}
QcxLat2 <- function(d=1,dx=sqrt(2)*d) {
  if(d<dx/2) stop("d must be larger than dx/2")
  dy <- sqrt(4*d^2-dx^2)
  VecLat(matrix(c(dx,0,dx/2,dy/2),2,2))
}
RectLat3 <- function(dx=1,dy=dx,dz=dx) {
  if(min(dx,dy,dz)<=0) stop("dx, dy and dz must be positive")
  VecLat(diag(c(dx,dy,dz)))
}
BCRectLat3 <- function(dx=1,dy=dx,dz=dx) {
  if(min(dx,dy,dz)<=0) stop("dx, dy and dz must be positive")
  VecLat(cbind(c(dx,0,0),c(0,dy,0),c(dx/2,dy/2,dz)))
}
FCRectLat3 <- function(d=1,dx=sqrt(2)*d,dz=d) {
  if(d<dx/2) stop("d must be larger than dx/2")
  dy <- sqrt(4*d^2-dx^2)
  VecLat(matrix(c(dx,0,0,dx/2,dy/2,0,dx/2,0,dz),3,3))
}

## Print method
setMethod("print","VecLat",
          function(x,...) {
            cat("An object of class \"VecLat\"\n")
            cat("dimspace: ",format(x@dimspace),"\n")
            cat("dimsupp: ",format(x@dimsupp),"\n")
            cat("determinant: ",format(x@det),"\n")
            cat("generating matrix: \n")
            print(x@gmat)
            cat("determinant: ",format(x@det),"\n")
          })

## Scaling 
setMethod("scaling",signature(x="VecLat",s="numeric"),
          function(x,s){
            VecLat(s*x@gmat)
          })
          
## Linear transform 
setMethod("ltransf",signature(x="VecLat",m="matrix"),
          function(x,m){
            VecLat(m%*%x@gmat)
          })
          
#######################################################################
#                  Definition of class FigLat                         #
#######################################################################

setClass("FigLat",representation(vlat="VecLat",fig="Figure",lmat="matrix"))

## Generator
FigLat<-function(d,vlat,fig,lmat=matrix(0,nrow=d,ncol=1)){
  if(!is.matrix(lmat)) lmat<-as.matrix(lmat)
  ## Use a check method
  if(vlat@dimspace!=d | nrow(lmat)!=d) stop("The slot dimspace of vlat and the number of rows of gmat must be equal to d")
  diml<-qr(lmat)$rank
  if(any(lmat!=0) & diml!=ncol(lmat)) stop("The matrix lmat must be either a null vector or a full rank matrix")
  if (vlat@dimsupp+diml!=d) stop("The slot dimsupp of vlat  must be equal to d-dim(L)")
  if(!all.equal(t(vlat@gmat)%*%lmat,matrix(0,nrow=vlat@dimsupp,ncol=ncol(lmat))))
    stop("The column vectors in the slot gmat of vlat must be orthogonal to the column vectors in lmat")
  if(fig@dimspace!=d) {
    stop("The figure (fig) does not lie in a d-dimensional space")
  }
 new("FigLat",vlat=vlat,fig=fig,lmat=lmat)
}
PPRectLat2<-function(hl=1,vl=hl,n=1,hp=hl/5,vp=hp,h3=TRUE){
  FigLat(d=2,vlat=RectLat2(hl,vl),fig=PP2(n,hp,vp,h3))
}

PPHexLat2<-function(delta=sqrt(2/sqrt(3)),n=1,hp=delta/5,vp=hp,h3=TRUE){
  FigLat(d=2,vlat=HexLat2(delta^2*sqrt(3)/2),fig=PP2(n,hp,vp,h3))
}

PPQcxLat2 <- function(d=1,dx=sqrt(2)*d,n=1,hp=dx/5,vp=hp,h3=TRUE) {
  FigLat(d=2,vlat=QcxLat2(d,dx),fig=PP2(n,hp,vp,h3))
}

QRectLat2<-function(hl=1,vl=hl,hq=hl/5,vq=hq){
  FigLat(d=2,vlat=RectLat2(hl,vl),fig=Quadrat(hq,vq))
}

QHexLat2<-function(delta=sqrt(2/sqrt(3)),hq=delta/5,vq=hq){
  FigLat(d=2,vlat=HexLat2(delta^2*sqrt(3)/2),fig=Quadrat(hq,vq))
}

QQcxLat2 <- function(d=1,dx=sqrt(2)*d,hq=dx/5,vq=hq) {
  FigLat(d=2,vlat=QcxLat2(d,dx),fig=Quadrat(hq,vq))
}

SRectLat2<-function(hl=1,vl=hl,end=c(hl/5,0)){
  FigLat(d=2,vlat=RectLat2(hl,vl),fig=Segment(end))
}
  
SHexLat2<-function(delta=sqrt(2/sqrt(3)),end=c(delta/5,0)){
  FigLat(d=2,vlat=HexLat2(delta^2*sqrt(3)/2),fig=Segment(end))
}

SQcxLat2 <- function(d=1,dx=sqrt(2)*d,end=c(dx/5,0)) {
  FigLat(d=2,vlat=QcxLat2(d,dx),fig=Segment(end))
}

LLat2<-function(delta=1,theta=0){
  FigLat(d=2,vlat=VecLat(delta*c(-sin(theta),cos(theta))),PointPattern(rep(0,2)),lmat=c(cos(theta),sin(theta)))
}

PPRectLat3<-function(dx=1,dy=dx,dz=dx,n=1,xp=dx/5,yp=xp,h3=TRUE){
  FigLat(d=3,vlat=RectLat3(dx,dy),fig=PP3(n,xp,yp,h3))
}

PPBCRectLat3 <- function(dx=1,dy=dx,dz=dx,n=1,xp=dx/5,yp=xp,h3=TRUE) {
  FigLat(d=3,vlat=BCRectLat3(dx,dy,dz),fig=PP3(n,xp,yp,h3))
}

PPFCRectLat3 <- function(d=1,dx=sqrt(2)*d,dz=d,n=1,xp=dx/5,yp=xp,h3=TRUE) {
  FigLat(d=3,vlat=FCRectLat3(d,dx,dz),fig=PP3(n,xp,yp,h3))
}

## Plot method
setMethod("plot",signature(x="FigLat",y="missing"),
          function(x,y,add=FALSE,xlim,ylim,xlab="",ylab="",...) {
              if(x@vlat@dimspace!=2)
                stop("Method plot for FigLat objects only implemented for 2D lattices of figures")
              if(add) {
                xlim <- par("usr")[1:2]
                ylim <- par("usr")[3:4]
              }
              if(!add)
                plot(mean(xlim),mean(ylim),xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,type="n")
              if(x@vlat@dimsupp==2) {
                E <- x@vlat@gmat
                Ei <- solve(E) # Inverse of the generating matrix E
                ## rect: 2x4 matrix containing the Cartesian coordinates of the box corners 
                rect <- matrix(c(xlim[1],ylim[1],xlim[2],ylim[1],xlim[2],ylim[2],xlim[1],ylim[2]),2,4)
                ## bb: 2x2 matrix. R=[bb[1,1],bb[1,2]]x[bb[1,1],bb[1,2]] is the smallest
                ## rectangle such that ER is contained in the box defined by xlim and
                ## ylim.
                bb <- t(apply(Ei%*%rect,1,range))
                ## Round bb.
                bb[,1] <- floor(bb[,1])
                bb[,2] <- ceiling(bb[,2])
                ## Add safety margins.
                bb <- bb + matrix(c(-1,-1,1,1),2,2)
                x1 <- bb[1,1]:bb[1,2]
                x2 <- bb[2,1]:bb[2,2]
                ## Lattice vector coordinates inside the box.
                lat <- E%*%cubicgrid(list(x1,x2))
                ## Plot the figure lattice.
                for(j in 1:ncol(lat))
                  plot(x@fig,origin=lat[,j],add=TRUE,...)
              } else if(x@vlat@dimsupp==1) {
                ## complete the vector lattice basis
                E <- cbind(x@vlat@gmat,c(-1,1)*x@vlat@gmat[2:1])
                Ei <- solve(E) # Inverse of the generating matrix E
                ## rect: 2x4 matrix containing the Cartesian coordinates of the box corners 
                rect <- matrix(c(xlim[1],ylim[1],xlim[2],ylim[1],xlim[2],ylim[2],xlim[1],ylim[2]),2,4)
                ## bb: 2x2 matrix. R=[bb[1,1],bb[1,2]]x[bb[1,1],bb[1,2]] is the smallest
                ## rectangle such that ER is contained in the box defined by xlim and
                ## ylim.
                bb <- t(apply(Ei%*%rect,1,range))
                ## Round bb.
                bb[,1] <- floor(bb[,1])
                bb[,2] <- ceiling(bb[,2])
                ## Add safety margins.
                bb <- bb + matrix(c(-1,-1,1,1),2,2)
                ## Après, ça se complique. En particulier, dans le cas d'un
                ## lattice de droites, la figure est un point et on doit
                ## dessiner des droites ??? 
                x1 <- bb[1,1]:bb[1,2]
                x2 <- mean(bb[2,])
                ## Lattice vector coordinates inside the box.
                lat <- E%*%cubicgrid(list(x1,x2))
                ## Plot the figure lattice.
                for(j in 1:ncol(lat))
                  plot(x@fig,x@lmat,origin=lat[,j],add=TRUE,...)
              }
            }
          )

##Print method
setMethod("print","FigLat",
          function(x,...) {
            d <- nrow(x@lmat)
            if(all(x@lmat==0)) {
              type<-switch(class(x@fig),
                           "Quadrat"="quadrats",
                           "Segment"="segments",
                           "PointPattern"="point patterns")}
            else {
              type<-switch(class(x@fig),"PointPattern"="lines","Segment"="strips")
            }
            cat(paste(d,"D lattice of ",type,"\n",sep=""))
            cat("Vector Lattice:\n")
            print(x@vlat)
            cat("Figure: ")
            if(all(x@lmat==0))
              print(x@fig)
            else {
              switch(type,
                     "lines"=
                     cat("Line through the origin and the point (",format(x@lmat),")\n"),
                     "strips"=
                     cat("Strip parallel to the vector (",format(x@lmat),
                         ") and of width",format(content(x@fig)),"\n")
                     )
            }
          })

## Linear transform
setMethod("ltransf",signature(x="FigLat",m="matrix"),
          function(x,m){
            vlat <- ltransf(x@vlat,m)
            fig <- ltransf(x@fig,m)
            lmat <- m%*%x@lmat
            d <- vlat@dimspace
            FigLat(d,vlat,fig,lmat)
          })

## Scaling
setMethod("scaling",signature(x="FigLat",s="numeric"),
          function(x,s) {
            FigLat(x@vlat@dimspace,scaling(x@vlat,s),scaling(x@fig,s),x@lmat)
          })

#######################################################################
#                  Definition of class FigLatData                     #
#######################################################################

setClass("FigLatData",representation(flat="FigLat",data="array"))

## Generator
FigLatData <- function(figlat,array){
  ## check that array dimension is consistent with figlat
  if(length(dim(array))!=figlat@vlat@dimsupp+1) {
    if(length(dim(array))==figlat@vlat@dimsupp) {
      dim(array) <- c(dim(array),1)
    } else {
      stop("Dimension of argument array is not consistent with argument figlat")
    }
  }
  new("FigLatData",flat=figlat,data=array)
}
## Print
setMethod("print","FigLatData",
          function(x,...) {
            print(x@data)
          })
## Empirical Covariogram
setMethod("dataCovariogram",signature(x="FigLatData",coord="matrix"),
          function(x,coord){
            n <- x@flat@vlat@dimsupp
            if(missing(coord)) {
              coord <- cbind(rep(0,n),diag(1,n))
            }
            E <- x@flat@vlat@gmat
            ns <- dim(x@data)[-length(dim(x@data))]
            nrep <- dim(x@data)[length(dim(x@data))]
            res <- NULL
            for(i in 1:ncol(coord)) {
              shift <- coord[,i]
              l <- 1+pmax(0,shift)
              l <- c(l,1)
              u <- pmin(ns,ns+shift)
              u <- c(u,nrep)
              lu <- mapply(function(x,y) x:y,l,u,SIMPLIFY=FALSE)
              idx <- as.matrix(expand.grid(lu))
              idx.shift <- idx-rep(c(shift,0),rep(nrow(idx),n+1))
              ghat <- sum(x@data[idx]*x@data[idx.shift])/nrep
              ghat <- x@flat@vlat@det*ghat
              res <- rbind(res,c(E%*%shift,ghat))
            }
            res       
          })

#######################################################################
#                           Functions                                 #
#######################################################################

M <- function(vlat,fig,s,L) {
  if(!is(vlat,"VecLat")) stop("Argument vlat must be of class FigLat")
  if(!is(fig,"Figure")) stop("Argument fig must be of class Figure")
  ## If the support of vlat is strictly included in the space,
  ## "projection" of vlat.
  if(vlat@dimspace!=vlat@dimsupp) {
    om <- t(svd(vlat@gmat,nv=0)$u)
    vlat <- ltransf(vlat,om)
    fig <- ltransf(fig,om)
  }
  d <- vlat@dimspace # dimension of the space
  Es<-dualmat(vlat@gmat)
  pz <- Ezeta(s,VecLat(Es),L=L,prepare=TRUE)
  detE <- 1/pz$determ # pz$determ = determinant of the dual lattice
  fig0<-scaling(fig,detE^(-1/d))
  res<-covariogram(fig0,function(h){Ezeta(h=h,prepare=pz)},sym=TRUE)
  res/content(fig0)^2
}

area.mse <- function(x,B=1,L=3){
  ## x: a lattice of figures, object of class FigLat.
  ## B: the perimeter. Default 1.
  ## L: an integer, the criterion for stopping summation of the
  ## Epstein Zeta function. Default 3.

  detL <- x@vlat@det # determinant of Lambda
  B/(4*pi^3)*detL^(ifelse(all(x@lmat==0),3/2,3))*M(x@vlat,x@fig,3,L)
}

vol.mse <- function(x,S=1,L=3){
  ## x: a lattice of figures, object of class FigLat.
  ## S: the surface area. Default 1.
  ## L: an integer, the criterion for stopping summation of the
  ## Epstein Zeta function. Default 3.

  detL <- x@vlat@det # determinant of Lambda
  p <- x@vlat@dimsupp # dimension of L
  S/(8*pi^3)*detL^(ifelse(all(x@lmat==0),4/3,3/(3-p)))*M(x@vlat,x@fig,4,L)
}

dvol.mse <- function(x,S=1,L=3){
  ## x: a lattice of figures, object of class FigLat.
  ## S: the surface area. Default 1.
  ## L: an integer, the criterion for stopping summation of the
  ## Epstein Zeta function. Default 3.

  detL <- x@vlat@det # determinant of Lambda
  p <- x@vlat@dimsupp # dimension of L
  d <- x@vlat@dimspace # space dimension
  S/(2*pi^2*sphereSurface(d))*detL^(ifelse(all(x@lmat==0),(d+1)/d,d/(d-p)))*M(x@vlat,x@fig,d+1,L)
}

area.mse.est <- function(fldata,mse.only=TRUE,iso=FALSE,diff2use) {
  d <- fldata@flat@vlat@dimspace
  if(d!=2) {
    stop("Argument fldata does not involve a figure lattice in 2D")
  }
  p <- fldata@flat@vlat@dimsupp
  if(missing(diff2use)) {
    diff2use <- diag(p)
    if(!iso) {
      if(p==2) {
        diff2use <- cbind(diff2use,c(1,1),c(-1,1))
      } else { # p==1
        stop("Argument fldata must involve a two-dimensional vector lattice")
      }
    }
  }
  res <- dvol.mse.est(fldata,mse.only=FALSE,iso=iso,diff2use=diff2use)
  if(!iso) {
    longRes <- list(B.est=res$S.est,deformation=res$deformation,mse.est=res$mse.est)
  } else { # isotropic case
    longRes <- list(B.est=res$S.est,mse.est=res$mse.est)
  }
  if(mse.only) {
    res$mse.est
  } else {
    longRes
  }
}

vol.mse.est <- function(fldata,mse.only=TRUE,iso=FALSE,diff2use) {
  d <- fldata@flat@vlat@dimspace
  if(d!=3) {
    stop("Argument fldata does not involve a figure lattice in 3D")
  }
  p <- fldata@flat@vlat@dimsupp
  if(missing(diff2use)) {
    diff2use <- diag(p)
    if(!iso) {
      if(p==3) {
        diff2use <- cbind(diff2use,c(1,1,0),c(1,0,1),c(0,1,1))
      } else { # p==1 or 2
        stop("Argument fldata must involve a three-dimensional vector lattice")
      }
    }
  }
  dvol.mse.est(fldata,mse.only,iso,diff2use)
}

dvol.mse.est <- function(fldata,mse.only=TRUE,iso=FALSE,diff2use) {
  p <- fldata@flat@vlat@dimsupp
  d <- fldata@flat@vlat@dimspace
  gdiff <- dataCovariogram(fldata,cbind(rep(0,p),diff2use))
  y <- gdiff[1,p+1]-gdiff[-1,p+1]
  fig <- fldata@flat@fig
  x <- rep(NA,length=length(y))
  phi <- function(h) vec.norm(e-h)-vec.norm(h)
  for(i in 2:nrow(gdiff)) {
    e <- gdiff[i,1:p]
    x[i-1] <- covariogram(fldata@flat@fig,phi,sym=FALSE)
  }
  S.est <- (sum(x*y)/sum(x^2))*(sphereSurface(d)*gamma((d+1)/2)/pi^((d-1)/2))
  mse.est <- dvol.mse(fldata@flat,S=S.est)
  if(!iso) { # anisotropic boundary
    ypred <- function(ediff,L) {
      A <- t(L)
      AE <- A%*%ediff
      phi <- function(h) vec.norm(Ae-A%*%h)-vec.norm(A%*%h)
      x <- rep(NA,ncol(ediff))
      for(i in 1:ncol(ediff)) {
        Ae <- AE[,i]
        x[i] <- covariogram(fig,phi,sym=FALSE)
      }
      x
    }
    cov.fit <- function(param,y,E,diff2use,fig) {
      L <- matrix(0,d,d)
      coefs <- c(param[-length(param)],1)
      L[lower.tri(L,diag=TRUE)] <- coefs
      L[d,d] <- 1/prod(diag(L))
      S <- param[length(param)]
      constant <- (1/sphereSurface(d))*(pi^((d-1)/2))/gamma((d+1)/2)
      sum((y-constant*S*ypred(E%*%diff2use,L))^2)
    }
    A0 <- diag(d)
    L0 <- A0[lower.tri(A0,diag=TRUE)]
    L0 <- L0[-length(L0)]
    S0 <- S.est
    lower <- matrix(-Inf,d,d)
    diag(lower) <- 1e-4
    lower <- lower[lower.tri(lower,diag=TRUE)]
    lower[length(lower)] <- 0
    sol <- optim(c(L0,S0),method="L-BFGS-B",lower=lower,
                 fn=cov.fit,y=y,E=fldata@flat@vlat@gmat,diff2use=diff2use,fig=fldata@flat@fig)
    L <- matrix(0,d,d)
    coefs <- c(sol$par[-length(sol$par)],1)
    L[lower.tri(L,diag=TRUE)] <- coefs
    L[d,d] <- 1/prod(diag(L))
    Aeigen <- eigen(L%*%t(L))
    D <- diag(Aeigen$values)
    P <- t(Aeigen$vectors)
    A <- sqrt(D)%*%P
    S.iso.est <- sol$par[length(sol$par)]
    radius <- (S.iso.est/sphereSurface(d))^(1/(d-1))
    radii <- radius/diag(D)
    S.est <- ellipsoidSurface(radii)
    mse.est <- dvol.mse(ltransf(fldata@flat,A),S=S.iso.est)
    longRes <- list(S.est=S.est,deformation=A,mse.est=mse.est)
  } else { # isotropic case
    longRes <- list(S.est=S.est,mse.est=mse.est)
  }
  if(mse.only) {
    mse.est
  } else {
    longRes
  }
}


latscale <- function(x,A,shape,CE.n,upper,maxiter=100,tol=.Machine$double.eps^0.25,
                     lower=.Machine$double.eps ^ 0.5,L=3,only.root=TRUE) {
  ## Test if the vector lattice is unit.
  if(all.equal(x@vlat@det,1)!=TRUE)
    stop("The vector lattice x@vlat must be unit.")
  ## Difference between the MSE associated with the scaling parameter and the nominal MSE. 
  f <- function(u,x,A,shape,CE.n,L) {
     1/(4*pi^3)*shape*u^3*M(scaling(x@vlat,u),x@fig,3,L)-CE.n^2*A^(3/2)
  }
  ## Find the root of f.
  res <- uniroot(f,lower=lower,upper=upper,maxiter=maxiter,tol=tol,x=x,A=A,shape=shape,CE.n=CE.n,L=L)
  ## Return result.
  if(only.root) res$root
  else {
    list(scale=res$root,CE=sqrt((res$f.root+CE.n^2*A^(3/2))/A^2),iter=res$iter,prec=res$estim.prec)
  }
}

cubicgrid <- function(s,d=2) {
  ## Compute the coordinates of the points in a set of the form
  ## s1 x...x sd where si are ordered subsets of reals.
  ##
  ## s: a list of vectors or a vector. In the latter case, s
  ## is transformed into the  list list(s,...,s) where s is
  ## replicated d-times.
  ## d: an integer, space dimension. By default, 2.
  ## Value: a matrix with d lines. Each column contains the
  ## coordinates of a integral point in the cube.
  if (mode(s)!="list") {
    sf <- vector(mode="list",length=d)
    for(i in 1:d) {
      sf[[i]] <- s
    }} else sf<-s
  k <- expand.grid(sf) # Returns a data frame.
  k <- t(as.matrix(k)) # matrix with d lines
  dimnames(k) <- NULL
  k
}

Ezeta <- function(s,vlat,h=rep(0,vlat@dimspace),L=3,prepare=FALSE,norm=TRUE) {
  if(is.logical(prepare)) {
    d <- vlat@dimspace # space dimension
    ## Normalise the lattice so that we get a unit lattice
    determ <- vlat@det
    vlat <- scaling(vlat,determ^(-1/d))
    
    E <- vlat@gmat
    tE <- t(E)
    ev <- sqrt(eigen(t(E)%*%E,only.values=TRUE)$values)
    Es <- dualmat(vlat@gmat)
    ev <- range(ev) # extremas of the eigen values
    crit <- L^2
    
    ## Prepare the computation of the first sum on the half lattice
    kmax <- floor(L/ev[1])
    k <- cubicgrid(-kmax:kmax,d)
    k <- k[,1:floor(ncol(k)/2),drop=FALSE] # keep only the half-lattice, origin excluded
    ## x: lattice points in the "smallest" parallelepiped 
    x <- E%*%k
    nx2 <- apply(x^2,2,"sum")
    ## Use only the half-ellipsoid
    keep <- nx2<=crit
    x <- x[,keep,drop=FALSE]
    nx2 <- nx2[keep]
    
    ## Prepare the computation of the second part, dual lattice
    lmax <- floor(L*ev[2])
    l <- cubicgrid(-lmax:lmax,d)
    l <- l[,-ceiling(ncol(l)/2)] # remove the origin
    y <- Es%*%l # lattice points in the smallest cube
    ny2 <- apply(y^2,2,"sum")
    y <- y[,ny2<=crit,drop=FALSE] # points in the ellipsoid
    
    ## Compute parts independent of h
    ## Part of the summands on the lattice
    tos <- gaic(s/2,pi*nx2)/(pi*nx2)^(s/2)
    
    if(prepare) {
      return(list(d=d,s=s,determ=determ,tE=tE,Es=Es,tos=tos,x=x,y=y))
    }
  } else {
    d <- prepare$d
    s <- prepare$s
    determ <- prepare$determ
    tE <- prepare$tE
    Es <- prepare$Es
    tos <- prepare$tos
    x <- prepare$x
    y <- prepare$y
  }

  ## Part of the computation which depends on the phase

  if(!is.matrix(h)) h <- as.matrix(h)
  ## Normalise the phase
  if(norm) h <- normphase(h,Es,tE)
  
  ## ps: cross-scalar products
  ps <- kronecker(matrix(1,1,ncol(h)),x)*kronecker(h,matrix(1,1,ncol(x)))
  ps <- matrix(apply(ps,2,"sum"),ncol(x),ncol(h))
  tos <- tos*cos(2*pi*ps)
  ## Summation  on the lattice
  res <- -2/s + 2*apply(tos,2,"sum")
  
  ## Compute the summands on the dual lattice
  # cross-differences
  dyh <- kronecker(matrix(1,1,ncol(h)),y)+kronecker(h,matrix(1,1,ncol(y)))
  nyh2 <- matrix(apply(dyh^2,2,"sum"),ncol(y),ncol(h))
  # Summand for l=0
  ih0 <- apply(h==0,2,all) # indices of h columns where h is non zero
  if(any(ih0))
    res[ih0] <- res[ih0] + 2/(s-d)
  if(any(!ih0)) {
    h2 <- apply(h[,!ih0,drop=FALSE]^2,2,"sum")
    res[!ih0] <- res[!ih0] + gaic((d-s)/2,pi*h2)/(pi*h2)^((d-s)/2)
  }
  ## Summation outside the origin
  res <- res + apply(gaic((d-s)/2,pi*nyh2)/(pi*nyh2)^((d-s)/2),2,"sum")

  ## Final result
  res*pi^(s/2)/gamma(s/2)
}

normphase <- function(h,E,Einv=solve(E)) {
  E%*%((Einv%*%h)%%1)
}

gaic <- function(a,x) {
  ## vectorized computation of the incomplete gamma
  ## x may be a vector or an array. The result has
  ## same length and dim attributes as x.
  if(is.array(x)) {
    dimx <- dim(x)
  } else {
    dimx <- NULL
  }
  res <- sapply(x,gamma_inc,a=a)
  if(!is.null(dimx)) {
    dim(res) <- dim(x)
  }
  res
}

##Dual matrix
dualmat <- function(x) {
  ## Compute the dual of a square matrix.
  ## x: non-singular square matrix.
  ## Value: dual matrix.
  
  t(solve(x))
}

in.upper.halfspace <- function(x) {
  ## Test if the vector x is in the upper half-space.
  ## Upper half-plane: null vector or last non-zero
  ## coordinate is greater than 0
  ##
  ## x: a vector of coordinates.
  ##
  ## Value: TRUE if x is in the upper half-space, FALSE otherwise.
  
  ifelse(all(x==0),TRUE,sign(x)[max((1:length(x))[x!=0])]==1)
}

pointcov <- function(x,tol=.Machine$double.eps ^ 0.5) {
  x <- as.matrix(x)
  n <- ncol(x) # n: number of points
  res <- list(ud=NULL,n=NULL)
  ## Computation of all differences
  xdiff <- kronecker(matrix(1,1,ncol(x)),x)-kronecker(x,matrix(1,1,ncol(x)))
  ## Consider auto-differences
  res$ud <- matrix(0,nrow(x),1)
  res$n <- n
  ## Keep only half of the differences
  xdiff <- xdiff[,(1:(n^2))[upper.tri(diag(n),diag=FALSE)],drop=FALSE]
  ## Start counting
  while(ncol(xdiff)>0) {
    ## count: number of xdiff col equal to xdiff[,1]
    count <- 1
    ## index: vector of j's such that xdiff[,j]=xdiff[,1]
    index <- 1
    if(ncol(xdiff)>1) {
      for(j in 2:ncol(xdiff)) {
        if(identical(all.equal(xdiff[,1],xdiff[,j],tol),TRUE)) {
          count <- count + 1
          index <- c(index,j)
        }
      }
    }
    count <- 2*count
    if(all(xdiff[,1]==0)) {
      res$n[1] <- res$n[1] + count
    } else {
      res$ud <- cbind(res$ud,xdiff[,1])
      res$n <- c(res$n,count)
    }
    xdiff <- xdiff[,-index,drop=FALSE]
  }
  ## Gather close differences
  if(tol>0) {
    ud <- res$ud
    d2a <- kronecker(matrix(1,1,ncol(ud)),ud)-kronecker(ud,matrix(1,1,ncol(ud)))
    d2b <- kronecker(matrix(1,1,ncol(ud)),ud)+kronecker(ud,matrix(1,1,ncol(ud)))
    id <- apply((abs(d2a)<=tol),2,"all") | apply((abs(d2b)<=tol),2,"all")
    id <- matrix(id,ncol(ud),ncol(ud))&upper.tri(diag(ncol(ud)),diag=FALSE)
    if(any(id)) {
      gat <- 1:ncol(ud)
      for(i in 1:ncol(ud))
        gat[id[i,]] <- gat[i]
      ugat <- unique(gat)
      for(i in ugat) {
        res$ud[,gat==i] <- res$ud[,i]
        res$n[i] <- sum(res$n[gat==i])
      }
      res$ud <- res$ud[,ugat]
      res$n <- res$n[ugat]
    }
  }
  ## Change difference signs in lower half-space
  tochg <- !apply(res$ud,2,"in.upper.halfspace")
  if(any(tochg)) 
    res$ud[,tochg] <- -res$ud[,tochg]

  if(sum(res$n)!=ncol(x)^2)
    stop("Bug in pointcov: total count less than n^2!")
  res
}

sphereSurface <- function(d=2) {
  2*(pi)^(d/2)/gamma(d/2)
}

ellipsoidSurface <- function(a) {
  ### according to Garry Tee (2005) Surface area and capacity of
  ### ellipsoids in n dimensions. New Zealand Journal of Mathematics,
  ### 34(2), 165--198.
  n <- length(a)
  a <- sort(a,decreasing=TRUE)
  delta <- 1-a[n]^2/a^2
  delta <- delta[-n]
  a <- a[-n]
  integrand <- if(n%%2==0 && n>=4) {
    function(x)
      2*x^(n-2)*sqrt((2-x^2)^(n-3)/apply(1-outer(delta,(1-x^2)^2),2,prod))*apply((1-delta)/(1-outer(delta,(1-x^2)^2)),2,sum)
  } else {
    function(h)
      sqrt((1-h^2)^(n-3)/apply(1-outer(delta,h^2),2,prod))*apply((1-delta)/(1-outer(delta,h^2)),2,sum)
  }
  2*prod(a)*pi^((n-1)/2)/gamma((n+1)/2)*integrate(integrand,lower=0,upper=1)$value
}

vec.norm <- function(x) {
  if(is.vector(x)) {
    sqrt(sum(x^2))
  } else {
    sqrt(apply(x^2,2,sum)) # norms of matrix columns
  }
}
