##' Tabulate grouped risks predicted by two different methods, models, algorithms
##'
##' All risks are multiplied by 100 before 
##' @title Risk reclassification table
##' @param list A list with two elements. Each element should either
##' be a vector with probabilities, or an object for which \code{predictStatusProb}
##' can extract predicted risk based on newdata.
##' @param newdata Passed on to \code{predictStatusProb}
##' @param cuts Risk quantiles to group risk
##' @param digits Number of digits to show for the predicted risks
##' @return reclassification table 
##' @seealso predictStatusProb
##' @examples
##'
#' set.seed(40)
#' N <- 40
#' X1 <- rnorm(N)
#' X2 <- rbinom(N,1,.4)
#' X3 <- rnorm(N)
#' expit <- function(x) exp(x)/(1+exp(x))
#' lp <- expit(X1 + X2 + X3)
#' Y <- factor(rbinom(N,1,lp))
#' dat <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3)
#' lm1 <- glm(Y~X1,data=dat,family="binomial")
#' lm2 <- glm(Y~X1+X2,data=dat,family="binomial")
#'
#' rc <- reclass(list("lrm.X1"=lm1,"lrm.X1.X2"=lm2),newdata=dat)
#' print(rc)
#' plot(rc)
#'
#' rc2 <- reclass(list("lrm.X1"=lm1,"lrm.X1.X2"=lm2),newdata=dat,cuts=c(0,5,10,50,100))
#' print(rc2)
#' plot(rc2)
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
reclass <- function(list,newdata,cuts=seq(0,100,25),digits=1){
    stopifnot(length(list)==2)
    predrisk <- lapply(list,function(x){
        if (class(x)[[1]]=="predictStatusProb")
            round(100*x,digits)
        else
            round(100*predictStatusProb(x,newdata),1)
    })
    retab <- table(cut(predrisk[[1]],
              cuts,
              include.lowest=TRUE,
              labels=paste(paste(cuts[-length(cuts)],cuts[-1],sep="-"),"%",sep="")),
          cut(predrisk[[2]],
              cuts,
              include.lowest=TRUE,
              labels=paste(paste(cuts[-length(cuts)],cuts[-1],sep="-"),"%",sep="")))
    nn <- names(list)
    if (!is.null(nn) & length(nn)==2)
        names(dimnames(retab)) <- nn
    out <- list(predRisk=predrisk,table=retab,cuts=cuts)
    class(out) <- "riskReclassification"
    out
}

##' @S3method print riskReclassification
print.riskReclassification <- function(x,...){
    print(x$table)
}

##' @S3method plot riskReclassification
plot.riskReclassification <- function(x,xlim=c(0,100),ylim=c(0,100),xlab,ylab,grid=TRUE,grid.col=gray(0.9),...){
    if (missing(xlab)) xlab <- paste("Risk (%):",names(dimnames(x$table))[[1]])
    if (missing(ylab)) ylab <- paste("Risk (%):",names(dimnames(x$table))[[2]])
    plot(x$predRisk[[1]],x$predRisk[[2]],axes=FALSE,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
    axis(1,at=x$cuts)
    axis(2,at=x$cuts)
    if (grid==TRUE)
        abline(h = x$cuts, v = x$cuts, col = gray(0.9))
}
