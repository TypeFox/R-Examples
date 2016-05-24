acfang.ltraj <- function(x, which=c("absolute","relative"),
                          nrep=999,lag=1, plot=TRUE,xlab="Lag",ylab="autocorrelation")
{
    if (!inherits(x,"ltraj"))
        stop("x should be of class ltraj")
    if (!is.regular(x))
        stop("x should be regular")
    if(length(lag)>1)
      stop("lag must be of length 1")
    whi.ang <- match.arg(which)

    foo <- function(df,whi.ang=whi.ang,nrepet=nrep,lag=1){
        if(whi.ang=="absolute"){
            y <- df[,"abs.angle"]
        }
        else if(whi.ang=="relative"){
            y <- df[,"rel.angle"]
        }
        ynoNA <- which(!is.na(y))
        res <- .C("acfangl", sim=as.double(rep(0,nrepet+1)),as.double(y),as.integer(length(y)),as.double(y),as.integer(length(ynoNA)),as.integer(ynoNA), as.integer(nrepet),as.integer(lag), NAOK=TRUE,PACKAGE = "adehabitatLT")
        return(c(obs=res$sim[1],quantile(res$sim,probs=c(0.025,0.5,0.975))))
    }
    res <- lapply(1:length(x),function(y) sapply(1:lag,function(z) foo(df=x[[y]],nrepet=nrep,whi.ang=whi.ang,lag=z)))
    res <- lapply(res,matrix,nrow=4,dimnames=list(c("obs","2.5%","50%","97.5%"),paste("lag",1:lag,sep=".")))
    if(plot){
      par(mfrow=n2mfrow(length(res)))
      for(i in 1:length(res)){
        plot(1:ncol(res[[i]]),res[[i]][1,],ylim=range(res),ty='n',xlab=xlab,ylab=ylab)
        polygon(c(1:ncol(res[[i]]),rev(1:ncol(res[[i]]))),c(res[[i]][2,],rev(res[[i]][4,])),col='lightgrey',border="grey")
        lines(1:ncol(res[[i]]),res[[i]][3,],col='darkgrey',ty="l",lty=2)
        points(1:ncol(res[[i]]),res[[i]][1,],ty='b',pch=ifelse(res[[i]][1,]< res[[i]][4,] & res[[i]][1,]> res[[i]][2,],1,15))
      }
      invisible(res)
    } else {
      return(res)
    }
}

acfdist.ltraj <- function(x,which=c("dist","dx","dy"),nrep=999,lag=1, plot=TRUE,xlab="Lag",ylab="autocorrelation")
{
    if (!inherits(x,"ltraj"))
      stop("x should be of class ltraj")
    if (!is.regular(x))
      stop("x should be regular")
    if(length(lag)>1)
      stop("lag must be of length 1")
    whi.lin <- match.arg(which)

    foo <- function(df,nrepet=nrep,lag,type=whi.lin){
      y <- df[,type]
      ynoNA <- which(!is.na(y))
      res <- .C("acfdist", sim=as.double(rep(0,nrepet+1)),as.double(y),as.integer(length(y)),as.double(y),as.integer(length(ynoNA)),as.integer(ynoNA), as.integer(nrepet),as.integer(lag), NAOK=TRUE,PACKAGE = "adehabitatLT")

      return(c(obs=res$sim[1],quantile(res$sim,probs=c(0.025,0.5,0.975))))
    }
    res <- lapply(1:length(x),function(burst) sapply(1:lag,function(z) foo(df=x[[burst]],nrepet=nrep,type=whi.lin,lag=z)))
    res <- lapply(res,matrix,nrow=4,dimnames=list(c("obs","2.5%","50%","97.5%"),paste("lag",1:lag,sep=".")))
    if(plot){
      par(mfrow=n2mfrow(length(res)))
      for(i in 1:length(res)){
        plot(1:ncol(res[[i]]),res[[i]][1,],ylim=range(res),ty='n',xlab=xlab,ylab=ylab)
        polygon(c(1:ncol(res[[i]]),rev(1:ncol(res[[i]]))),c(res[[i]][2,],rev(res[[i]][4,])),col='lightgrey',border="grey")
        lines(1:ncol(res[[i]]),res[[i]][3,],col='darkgrey',ty="l",lty=2)
        points(1:ncol(res[[i]]),res[[i]][1,],ty='b',pch=ifelse(res[[i]][1,]< res[[i]][4,] & res[[i]][1,]> res[[i]][2,],1,15))
      }
      invisible(res)
    } else {
      return(res)
    }
  }

