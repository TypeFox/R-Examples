#-------------------------------------------------------
#
#  Created       : 29/10/02
#  Last Modified : Time-stamp: <2013-12-02 19:32:49 antoine>
#
#  Description   : Principal component analysis
#                  
#  Author        : Antoine Lucas
#                  lucas@toulouse.inra.fr
#
#  Licence       : GPL 
#
#-------------------------------------------------------


acp <- function(x,center=TRUE,reduce=TRUE,wI=rep(1,nrow(x)),wV=rep(1,ncol(x)))
{   
    x    <- as.matrix(x)
    if(center)
      x <- t(t(x) - as.vector(( wI %*% x)/sum(wI)))
##    x    <- scale(x ,center = center, scale = FALSE)
    if (reduce) 
      x    <- apply(x,2,function(u) { u/sd(u)}) 

    ##              Di.X'.Dv.X
    EIG  <- eigen( (t(x)* wI) %*% (x * wV) ,symmetric=FALSE) 
    V    <- EIG$vector    # ou bien: V=svd(x)$v

    EIG$values <- Re(EIG$values)
    V    <- V %*% diag(sign(EIG$values))
    val  <- sqrt(abs(EIG$values))

    scores <- x %*% V

    V      <- as.matrix(Re(V))
    scores <- as.matrix(Re(scores))

    dimnames(V)[[2]] <- paste("Comp",1:dim(x)[2])
    if(!is.null( dimnames(x)[[2]] ))
      dimnames(V)[[1]] <- dimnames(x)[[2]]
    if(!is.null(dimnames(x)[[1]]))
      dimnames(scores)[[1]] <- dimnames(x)[[1]]
    dimnames(scores)[[2]] <- paste("Comp",1:dim(x)[2])

    ##cmpr <- x %*% (sqrt(wV) * as.matrix(V))
    
    sdev   <- apply(scores,2,sd)    
    res  <- list(eig=val,sdev=sdev,scores=scores,loadings=V)
    class(res) <- "acp"
    res
}
pca <- acp


print.acp <- function(x, ...)
{
    #cat("Call:\n"); dput(x$call)
    cat("\nStandard deviations:\n")
    print(x$sdev, ...)
    cat("\nEigen values:\n")
    print(x$eig, ...)
    invisible(x)
}


# 
#   SECTION GRAPHIQUES
#

plot.acp <- function(x,i=1,j=2,text=TRUE,label='Composants',col='darkblue',main='Individuals PCA',variables=TRUE,individual.label=NULL,...)
{
    U    <- x$scores
    XLAB <- paste(label,i)
    YLAB <- paste(label,j)
    plot.new()
    plot.window(range(U[,i]),range(U[,j]))
    axis(1,labels=TRUE,tick=TRUE)
    axis(2,labels=TRUE,tick=TRUE)
    box()
    
    title(xlab=XLAB,ylab=YLAB,main=main)
    if(text){
      if(is.null(individual.label))
        {
          individual.label=dimnames(x$scores)[[1]]
        }
        text(labels=individual.label,U[,i],U[,j],col=col,...)   
    }
    else{
        points(U[,i],U[,j],col=col,...) 
    }
    if(variables)
      {
         par(new=TRUE)
         biplot.acp(x,circle=FALSE,label="",main="")
       }

}

biplot.acp <- function(x,i=1,j=2,label='Composants',col='darkblue',length=0.1,main='Variables PCA',circle=TRUE,...)
{
    U    <- x$loadings
    LIM  <- c(-1.3,1.3)
    XLAB <- paste(label,i)
    YLAB <- paste(label,j)

    # PLOT DES AXES
    plot.new()
    plot.window(LIM,LIM)
    axis(1,labels=TRUE,tick=TRUE)
    axis(2,labels=TRUE,tick=TRUE)
    box()
    title(xlab=XLAB,ylab=YLAB,main=main)


    # PLOT DU NOM DES FLECHES
    text(x=U[,i]*1.3,y=U[,j]*1.3,labels=dimnames(U)[[1]],col=col)   

    # PLOT DES FLECHES
    arrows(0,0,U[,i],U[,j],length = length,col=col)

    # CERCLE
    if(circle)
      {
        t2p <- 2 * pi * seq(0,1, length = 200)
        xc <- cos(t2p)
        yc <- sin(t2p)
        lines(xc,yc,col='darkblue')
      }
}

# Graphique: Eboulis des valeurs propres
plot2 <- function(x,pourcent=FALSE,eigen=TRUE,label='Comp.',col='lightgrey',main='Scree Graph',ylab='Eigen Values')
{
    if(eigen){ U <- x$eig }
    else { U <- x$sdev }

    if(pourcent){U <- U/sum(U) }
    n     <- length(U)
    names <- paste(label,1:n)
    barplot(U,main=main,ylab=ylab,col=col,names.arg=names)
}


plotAll <- function(x)
  {
    par(mfrow=c(2,2))
    plot2(x)
    ##    boxplot(as.list(as.data.frame(x$cmpr)))
    plot(x,variables=FALSE)
    biplot(x)
    plot(x,main="Both",variables=TRUE)
  }
