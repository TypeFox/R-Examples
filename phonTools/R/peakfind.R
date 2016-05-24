# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


peakfind = function (x, show = TRUE){
    rightbig = (diff(c(x[1]+1,x)) > 0)^2
    leftbig = rev((diff(c(rev(x)[1]+1,rev(x))) > 0)^2)

    peaks = leftbig + rightbig
    npeaks = sum(peaks == 2)
    if (npeaks==0) return (0)
    if (npeaks>0) peaks = which(peaks==2)
    if (show == TRUE & npeaks>0){
        plot(x, lwd = 2, type = "l", ylab = "Value", xaxs = "i")
        points(peaks, x[peaks], col = 2, pch = 17, cex = 1)
    }
    invisible(peaks)
}

