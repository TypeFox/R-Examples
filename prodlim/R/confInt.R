#' Add point-wise confidence limits to the graphs of Kaplan-Meier and
#' Aalen-Johansen estimates of survival and cumulative incidence.
#' 
#' This function is invoked and controlled by \code{plot.prodlim}.
#' 
#' This function should not be called directly. The arguments can be specified
#' as \code{Confint.arg} in the call to \code{plot.prodim}.
#' 
#' @param x an object of class `prodlim' as returned by the \code{prodlim}
#' function.
#' @param times where to compute point-wise confidence limits
#' @param newdata see \code{plot.prodim}
#' @param type Either \code{"cuminc"} or \code{"survival"} passed to
#' summary.prodlim as \code{surv=ifelse(type=="cuminc",FALSE,TRUE)}.
#' @param citype If \code{"shadow"} then confidence limits are drawn as colored
#' shadows.  Otherwise, dotted lines are used to show the upper and lower
#' confidence limits.
#' @param cause see \code{plot.prodim}
#' @param col the colour of the lines.
#' @param lty the line type of the lines.
#' @param lwd the line thickness of the lines.
#' @param density For \code{citype="shadow"}, the density of the shade. Default
#' is 55 percent.
#' @param \dots Further arguments that are passed to the function
#' \code{segments} if \code{type=="bars"} and to \code{lines} else.
#' @return Nil
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{plot.prodlim}}, \code{\link{atRisk}},
#' \code{\link{markTime}}
#' @keywords survival
#' @export
confInt <- function(x,
                    times,
                    newdata,
                    type,
                    citype,
                    cause,
                    col,
                    lty,
                    lwd,
                    density=55,
                    ...){
    ## if (citype=="shadow" && length(times)>100 && exact==FALSE)
    ## times <- seq(min(times),max(times),diff(range(times)/100))
    sumx <- summary(x,times=times,newdata=newdata,cause=cause,verbose=FALSE,surv=ifelse(type=="cuminc",FALSE,TRUE))$table
    if (x$model=="competing.risks" && x$covariate.type>1) sumx <- sumx[[1]]
    ## if (x$model=="survival" && x$covariate.type==1) sumx <- list(sumx)
    if (!is.list(sumx)) sumx <- list(sumx)
    nlines <- length(sumx)
    ci <- lapply(sumx,function(u){
        uu <- data.frame(u[,c("time","lower","upper"),drop=FALSE])
        uu=uu[!is.na(uu$lower),]
        # ----------remove confidence limits before the first event----------
        est <- u[!is.na(u[,"lower"]),type]
        cond <- est <1 & est>0
        uu=uu[((uu$upper-uu$lower)<1 | cond),]
        uu
    })
    nix <- lapply(1:nlines,function(i){
        if (NROW(ci[[i]])>0){
            switch(citype,
                   "bars"={
                       segments(x0=ci[[i]]$time,
                                x1=ci[[i]]$time,
                                y0=ci[[i]]$lower,
                                y1=ci[[i]]$upper,
                                lwd=lwd[i],
                                col=col[i],
                                lty=lty[i],
                                ...)
                   },
                   "shadow"={
                       cc <- dimColor(col[i],density=density)
                       ## ccrgb=as.list(col2rgb(col[i],alpha=TRUE))
                       ## names(ccrgb) <- c("red","green","blue","alpha")
                       ## ccrgb$alpha=density
                       ## cc=do.call("rgb",c(ccrgb,list(max=255)))
                       ttt <- ci[[i]]$time
                       nt <- length(ttt)
                       ttt <- c(ttt,ttt)
                       uuu <- c(0,ci[[i]]$upper[-nt],ci[[i]]$upper)
                       lll <- c(0,ci[[i]]$lower[-nt],ci[[i]]$lower)
                       neworder <- order(ttt)
                       uuu <- uuu[neworder]
                       lll <- lll[neworder]
                       ttt <- sort(ttt)
                       polygon(x=c(ttt,rev(ttt)),
                               y=c(lll,rev(uuu)),col=cc,border=NA)
                       ## xx=ci[[i]]$time
                       ## nix <- sapply(1:length(xx),function(b){
                       ## rect(xleft=xx[b],xright=xx[b+1],ybottom=ci[[i]]$lower[b],ytop=ci[[i]]$upper[b],col=cc,border=NA)
                       ## })
                   },{
                       lines(x=ci[[i]]$time,ci[[i]]$lower,type="s",lwd=lwd[i],col=col[i],lty=lty[i],...)
                       lines(x=ci[[i]]$time,ci[[i]]$upper,type="s",lwd=lwd[i],col=col[i],lty=lty[i],...)
                   })
        }
    })
}
