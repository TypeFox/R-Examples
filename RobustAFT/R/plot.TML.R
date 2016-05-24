plot.TML<-function(x, which = 1:3, caption = c("Residual QQ-plot",
  "Response vs. Fitted Values", "Standardized Residuals vs. Fitted Values"),
  panel = points, sub.caption = deparse(x$call$formula), main = "",
  ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...)
{
    if(!inherits(x, "TML"))
        stop("Use only with 'TML' objects")
    show <- rep(FALSE, 3)
    if(!is.numeric(which) || any(which < 1) || any(which > 3))
        stop("which must be in 1:3")
    show[which] <- TRUE
    r <- residuals(x)
    n <- length(r)
    sr <- r/x$v1
    yh <- fitted(x)
    one.fig <- prod(par("mfcol")) == 1
    if(ask){
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if(show[1]){
        if(x$errors == "Gaussian")
            qw <- qnorm(ppoints(sr))
        if(x$errors == "logWeibull")
            qw <- log(qweibull(ppoints(sr), shape = 1))
        plot(qw, sort(sr),ylim=c(-(max(abs(min(sr)),max(sr))+0.5),max(abs(min(sr)),max(sr))+0.5),
             xlab = "Theoretical Quantiles", ylab = "Standardized Residuals")
        abline(c(0, 1))
        abline(h = c(x$tl, x$tu), lty = 3)
        if (!is.null(caption[1]))
            mtext(paste(caption[1], " with ", sQuote(x$errors), " errors"), 3, 0.25)
        if(one.fig)
            title(sub = sub.caption, ...)
    }
    if(show[2]){
        y <- if(!is.null(x$model))
            model.response(x$model)
        else yh + r
        plot(yh, y, xlab = "Fitted Values", ylab = "Response",
            main = main, ...)
        mtext(caption[2], 3, 0.25)
        if(one.fig)
            title(sub = sub.caption, ...)
        abline(a = 0, b = 1)
    }
    if(show[3]){
        plot(yh, sr,ylim=c(-(max(abs(min(sr)),max(sr))+0.5),max(abs(min(sr)),max(sr))+0.5),
             xlab = "Fitted Values", ylab = "Standardized Residuals", main = main, ...)
        mtext(caption[3], 3, 0.25)
        if(one.fig)
            title(sub = sub.caption, ...)
        abline(h = c(x$tu, 0, x$tl), lty = 3)
        for(i in 1:length(sr)){
          if(weights(x)[i]==0){
            if(sr[i] < 0) text(yh[i],sr[i],label=i,pos=1,offset=0.5,cex=0.7)
            else text(yh[i],sr[i],label=i,pos=3,offset=0.5,cex=0.7)
          }
        }
    }
  invisible()
}
