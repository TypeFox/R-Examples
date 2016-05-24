plot.glinternet.cv = function(x, ...){
    #plot cv curve
    bestIndex = which.min(x$cvErr)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        cvErr=x$cvErr
        cvErrStd=x$cvErrStd
        lambdaIdx=1:length(x$lambda)
        plt = ggplot2::ggplot(data.frame(cvErr, cvErrStd, lambdaIdx), ggplot2::aes(x=lambdaIdx, y=cvErr)) + ggplot2::geom_line(size=1.5) + ggplot2::geom_vline(xintercept=bestIndex, linetype="dashed") + ggplot2::scale_x_discrete("Lambda index", breaks=seq(2,length(x$lambda),2)) + ggplot2::scale_y_continuous("CV error") + ggplot2::theme_bw() + ggplot2::geom_ribbon(ggplot2::aes(ymin=cvErr-cvErrStd, ymax=cvErr+cvErrStd), fill="black", alpha=0.5)
        print(plt)
    }
    else {
        barwidth = 0.25
        x = 1:length(x$lambda)
        y = x$cvErr
        delta = x$cvErrStd
        plot(x, x$cvErr, type="l", lwd=1.5, xlab="Lambda index", ylab="CV error", xaxt="n", ylim=c(min(y-delta), max(y+delta)))
        segments(x-barwidth, y+delta, x+barwidth, y+delta)
        segments(x-barwidth, y-delta, x+barwidth, y-delta)
        abline(v=bestIndex, lty=3)
        axis(1, at=seq(2, length(x), 2), las=1)
    }
}
