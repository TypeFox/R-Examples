#' Quantiles for Kaplan-Meier and Aalen-Johansen estimates.
#' 
#' Quantiles for Kaplan-Meier and Aalen-Johansen estimates.
#' 
#' 
#' @param x Object of class \code{"prodlim"}.
#' @param q Quantiles. Vector of values between 0 and 1.
#' @param cause For competing risks the cause of interest.
#' @param ... not used
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @keywords survival
#' @examples
#' library(lava)
#' set.seed(1)
#' d=SimSurv(30)
#' f=prodlim(Hist(time,status)~1,data=d)
#' f1=prodlim(Hist(time,status)~X1,data=d)
#' # default: median and IQR
#' quantile(f)
#' quantile(f1)
#' # median alone
#' quantile(f,.5)
#' quantile(f1,.5)
#'
#' # competing risks
#' set.seed(3)
#' dd = SimCompRisk(30)
#' ff=prodlim(Hist(time,event)~1,data=dd)
#' ff1=prodlim(Hist(time,event)~X1,data=dd)
#' ## default: median and IQR
#' quantile(ff)
#' quantile(ff1)
#' 
#' print(quantile(ff1),na.val="NA")
#' print(quantile(ff1),na.val="Not reached")
#' 
#' @export 
"quantile.prodlim" <- function(x,
                               q,
                               cause=1,
                               ...){
    ## require(stats)
    ## stopifnot(x$model=="survival")
    etype <- attr(x$model.response,"entry.type")
    if (!is.null(etype) && etype=="leftTruncated")
        stop("Don't know how to compute quantiles with delayed entry (left-truncation).")
    if(x$model=="survival"){
        if (missing(q)) q <- c(1,.75,0.5,.25,0)
        q <- 1-q ## since this is a survival function
        sumx <- summary(x,newdata=x$X,times=x$time,showTime=TRUE,verbose=FALSE)
        getQ <- function(sum){
            out <- do.call("cbind",lapply(c("surv","lower","upper"),function(w){
                                               notna=is.na(sum[,w])
                                               xxx=as.numeric(sum[,w,drop=TRUE][!notna])
                                               ttt=as.numeric(sum[,"time"][!notna])
                                               found <- 2+sindex(jump.times=xxx,eval.times=q,comp="greater",strict=FALSE)
                                               inner <- c(as.vector(c(0,ttt)[found]))
                                               inner
                                           }))
            out <- data.frame(out)
            out <- cbind(q,out)
            names(out) <- c("q","quantile","lower","upper")
            out}
        if (sumx$cotype==1) out <- list("quantiles.survival"=getQ(sumx$table))
        else out <- lapply(sumx$table,getQ)
        attr(out,"model") <- x$model
        class(out) <- "quantile.prodlim"
        out
    } else{
          ## cumulative incidence, competing risks
          if (missing(q)) q <- c(0,0.25,0.5,0.75,1)
          sumx <- summary(x,newdata=x$X,times=x$time,showTime=TRUE,verbose=FALSE,cause=cause)
          getQ <- function(sum){
              out <- do.call("cbind",lapply(c("cuminc","lower","upper"),function(w){
                                                 notna=is.na(sum[,w])
                                                 xxx=as.numeric(sum[,w][!notna])
                                                 ttt=as.numeric(sum[,"time"][!notna])
                                                 ## found <- 2+sindex(jump.times=xxx,eval.times=q,comp="greater",strict=FALSE)
                                                 found <- 2+sindex(jump.times=xxx,eval.times=q,comp="smaller",strict=FALSE)
                                                 inner <- c(as.vector(c(0,ttt)[found]))
                                                 inner
                                             }))
              out <- data.frame(out)
              out <- cbind(q,out)
              ## upper is lower and lower is upper
              names(out) <- c("q","quantile","upper","lower")
              out <- out[,c("q","quantile","lower","upper")]
              out}
          if (sumx$cotype==1)
              out <- list("quantiles.cuminc"=getQ(sumx$table[[1]]))
          else {
              out <- lapply(sumx$table[[1]],getQ)
          }
          attr(out,"model") <- x$model
          class(out) <- "quantile.prodlim"
          out
      }
}
      

