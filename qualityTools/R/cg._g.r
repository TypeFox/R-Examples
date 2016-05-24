.cg=function(x, target, tolerance=c(-1,1), ref.interval, facCg, facCgk)
{
    if (missing(x)) 
        stop("x must be given as a vector")
    if (missing(target)) 
        target = mean(x)
    if (missing(ref.interval)) 
        ref.interval = pnorm(3) - pnorm(-3)
    sd = sd(x)
    mean = mean(x)
    ref.ar = qnorm(ref.interval, mean, sd) - qnorm(1 - ref.interval, mean, sd)
    if (missing(facCg)) 
        facCg = 0.2
    if (missing(facCgk)) 
        facCgk = 0.1
    if (missing(tolerance)) 
        stop("tolerance is missing!")
    if (length(tolerance) != 2) 
        stop("tolerance has wrong length")                                         
 Cg = (facCg * (abs(diff(tolerance))))/ref.ar
 Cgk = (facCgk * (abs(diff(tolerance))) - abs(target - mean))/(ref.ar/2)
 return(list(Cg, Cgk))
}

cgRunChart=function(x, target, tolerance, ref.interval, facCg, facCgk, n = 0.2, type, col, pch, xlim, ylim,main, conf.level = 0.95, cex.val = 1, cgOut = TRUE)
{
    if (missing(x)) 
        stop("x must be given as a vector")
    if (missing(target)) {
        target = mean(x)
        targetmissing = FALSE
    }
    else targetmissing = TRUE
    if (missing(ref.interval)) 
        ref.interval = pnorm(3) - pnorm(-3)
    sd = sd(x)
    mean = mean(x)
    ref.ar = qnorm(ref.interval, mean, sd) - qnorm(1 - ref.interval, mean, sd)
    if (missing(facCg)) 
        facCg = 0.2
    if (missing(facCgk)) 
        facCgk = 0.1
    if (missing(tolerance)) 
        warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
    if (missing(tolerance)) {
        width = ref.ar/facCg
        tolerance = numeric(2)
        tolerance[1] = mean(x) - width/2
        tolerance[2] = mean(x) + width/2
    }
    quant1 = qnorm((1 - ref.interval)/2, mean, sd)
    quant2 = qnorm(ref.interval + (1 - ref.interval)/2, mean, sd)
    if (length(tolerance) != 2) 
        stop("tolerance has wrong length")
    if (missing(type)) 
        type = "b"
    if (missing(col)) 
        col = 1
    if (missing(pch)) 
        pch = 19
    if (missing(xlim)) 
        xlim = c(0, length(x))
    if (missing(ylim)) 
        ylim = c(min(x, target - n/2 * (abs(diff(tolerance))), quant1, quant2), max(x, target + n/2 * (abs(diff(tolerance))), quant1, quant2))
    if (missing(main))
        main="Run Chart"
    Cg = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[1]]
    Cgk = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[2]]
    par.old = par()$mar
    par(mar = c(5.1, 4.1, 4.1, 5.1))
     plot(x, type = type, col = col, main = main, pch = pch, xlim = c(xlim[1] - 0.05 * xlim[2], xlim[2]), ylim = ylim)
    abline(h = target)
    rect(xleft = min(xlim) - 0.09 * max(xlim), xright = -0.005 * max(xlim), ybottom = min(ylim), ytop = max(ylim), col = "white", border = "white")
    text(min(xlim) - 0.1 * max(xlim), target, pos = 4, srt = 0, "target", xpd = TRUE, bg = "white", cex = cex.val)
    if (targetmissing == TRUE) 
        abline(h = mean, lty = 2, col = "seagreen")
    abline(h = quant1, col = "seagreen", lty = 3)
    abline(h = quant2, col = "seagreen", lty = 3)
    abline(h = c(target + n/2 * (abs(diff(tolerance))), target - n/2 * (abs(diff(tolerance)))))
    lines(lowess(x), col = "red")
    temp = axis(4, at = c(target + n/2 * (abs(diff(tolerance))), target - n/2 * (abs(diff(tolerance))), mean(x)), labels = FALSE)
    text(max(xlim) + 0.075 * max(xlim), temp[3], srt = 0, adj = 0, substitute(x[tar] + a %.% b, list(a = round(n/2, 4), b = "T")), xpd = TRUE, cex = cex.val)
    text(max(xlim) + 0.075 * max(xlim), temp[2], srt = 0, adj = 0, expression(bar(x)), xpd = TRUE, cex = cex.val, col = "seagreen")
    text(max(xlim) + 0.075 * max(xlim), temp[1], srt = 0, adj = 0, substitute(x[tar] - a %.% b, list(a = round(n/2, 4), b = "T")), xpd = TRUE, cex = cex.val)
    lt1 = round(((1 - ref.interval)/2) * 100, 3)
    lt2 = round(((ref.interval + (1 - ref.interval)/2)) * 100, 3)
    text(max(xlim) + 0.075 * max(xlim), quant1, srt = 0, adj = 0, substitute(x[a * b], list(a = lt1, b = "%")), xpd = TRUE, cex = cex.val, col = "seagreen")
    text(max(xlim) + 0.075 * max(xlim), quant2, srt = 0, adj = 0, substitute(x[a * b], list(a = lt2, b = "%")), xpd = TRUE, cex = cex.val, col = "seagreen")
    if(cgOut == TRUE)
    {
     legend("topright",legend=(c(paste("Cg: ",Cg),paste("Cgk :",Cgk))),inset=c(0,0.06),bty="n")
    }
    par(mar = par.old)
invisible(list(Cg, Cgk))
}

cgHist=function(x, target, tolerance, ref.interval, facCg, facCgk, n = 0.2, col, xlim, ylim, main, conf.level = 0.95, cex.val = 1, cgOut = TRUE)
{
    if (missing(x)) 
        stop("x must be given as a vector")
    if (missing(target)) {
        target = mean(x)
        targetmissing = FALSE
    }
    else targetmissing = TRUE
    if (missing(ref.interval)) 
        ref.interval = pnorm(3) - pnorm(-3)
    sd = sd(x)
    mean = mean(x)
    ref.ar = qnorm(ref.interval, mean, sd) - qnorm(1 - ref.interval, mean, sd)
    if (missing(facCg)) 
        facCg = 0.2
    if (missing(facCgk)) 
        facCgk = 0.1
    if (missing(tolerance)) 
        warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
    if (missing(tolerance)) {
        width = ref.ar/facCg
        tolerance = numeric(2)
        tolerance[1] = mean(x) - width/2
        tolerance[2] = mean(x) + width/2
    }
    quant1 = qnorm((1 - ref.interval)/2, mean, sd)
    quant2 = qnorm(ref.interval + (1 - ref.interval)/2, mean, sd)
    if (length(tolerance) != 2) 
        stop("tolerance has wrong length")
    if (missing(col)) 
        col = "lightblue"
    if (missing(xlim)) 
        xlim = c(0, length(x))
    if (missing(ylim)) 
        ylim = c(min(x, target - n/2 * (abs(diff(tolerance))), quant1, quant2), max(x, target + n/2 * (abs(diff(tolerance))), quant1, quant2))
    if (missing(main))
        main=paste("Histogram of", deparse(substitute(x)), "- target")
    Cg = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[1]]
    Cgk = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[2]]
    
    x.c = x - target
    temp = hist(x.c, plot = FALSE)
 #   par.old = par()$mar
#    par(mar = c(3.1, 4.1, 0.55, 2.1))
    hist(x.c, freq = FALSE, ylim = c(0, max(density(x.c)$y, temp$density) + 0.4 * max(density(x.c)$y, temp$density)), col = col, main=main)
    test = t.test(x.c, mu = 0, conf.level = conf.level)
    lines(x = c(0, 0), y = c(0, max(density(x.c)$y, temp$density) + 0.3 * max(density(x.c)$y, temp$density)), col = "red")
    lines(x = c(test$conf.int[1], test$conf.int[1]), y = c(0, max(density(x.c)$y, temp$density) + 0.3 * max(density(x.c)$y, temp$density)), col = "blue", lty = 3)
    lines(x = c(test$conf.int[2], test$conf.int[2]), y = c(0, max(density(x.c)$y, temp$density) + 0.3 * max(density(x.c)$y, temp$density)), col = "blue", lty = 3)
    
    if(cgOut == TRUE)
    {
     legend("topright", legend = c("conf.int",paste("Cg: ",Cg),paste("Cgk :",Cgk)), col = c("blue",-1,-1), lty = c(3,-1,-1), inset = c(0.01, 0.06))
    }
    else
     legend("topright", legend = c("conf.int"), col = c("blue"), lty = c(3), inset = c(0.01, 0.06))
#    par(mar = par.old)
    legend("topleft", legend = c(expression(paste(H[0], " : Bias = 0")), paste("t-value: ", round(test$statistic, 3)), paste("p-value: ", round(test$p.value, 
        3))), inset = c(-0.01, 0.04), bty = "n")
    lines(density(x.c))
 #   par(mar = par.old)
 #   par(mar = par.old)
invisible(list(Cg, Cgk))
}

cgToleranceView=function(x, target, tolerance, ref.interval, facCg, facCgk, n = 0.2, type, col, pch, xlim, ylim, main, conf.level = 0.95, cex.val = 1, cgOut = TRUE)
{
    if (missing(x)) 
        stop("x must be given as a vector")
    if (missing(target)) {
        target = mean(x)
        targetmissing = FALSE
    }
    else targetmissing = TRUE
    if (missing(ref.interval)) 
        ref.interval = pnorm(3) - pnorm(-3)
    sd = sd(x)
    mean = mean(x)
    ref.ar = qnorm(ref.interval, mean, sd) - qnorm(1 - ref.interval, mean, sd)
    if (missing(facCg)) 
        facCg = 0.2
    if (missing(facCgk)) 
        facCgk = 0.1
    if (missing(tolerance)) 
        warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
    if (missing(tolerance)) {
        width = ref.ar/facCg
        tolerance = numeric(2)
        tolerance[1] = mean(x) - width/2
        tolerance[2] = mean(x) + width/2
    }
    quant1 = qnorm((1 - ref.interval)/2, mean, sd)
    quant2 = qnorm(ref.interval + (1 - ref.interval)/2, mean, sd)
    if (length(tolerance) != 2) 
        stop("tolerance has wrong length")
    if (missing(type)) 
        type = "b"
    if (missing(col)) 
        col = 1
    if (missing(pch)) 
        pch = 19
    if (missing(xlim)) 
        xlim = c(0, length(x))
    if (missing(ylim)) 
        ylim = c(min(x, target - n/2 * (abs(diff(tolerance))), quant1, quant2), max(x, target + n/2 * (abs(diff(tolerance))), quant1, quant2))
    if (missing(main))
        main="Tolerance View"   
    Cg = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[1]]
    Cgk = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[2]]    
    plot(x, type=type, col = col, main = main, pch = pch, xlim = xlim, ylim = c(min(x, tolerance[2], tolerance[1], target + n/2 * (tolerance[2] - 
        tolerance[1]), target - n/2 * (tolerance[2] - tolerance[1])), max(x, tolerance[2], tolerance[1], target + n/2 * (tolerance[2] - tolerance[1]), target - 
        n/2 * (tolerance[2] - tolerance[1]))))
    abline(h = c(tolerance[1], tolerance[2]), lty = 2, col = "red")
    abline(h = target)
    abline(h = c(target + n/2 * (tolerance[2] - tolerance[1]), target - n/2 * (tolerance[2] - tolerance[1])))
    if(cgOut == TRUE)
    {
     legend("topright",legend=(c(paste("Cg: ",Cg),paste("Cgk :",Cgk))),inset=c(0,0.06),bty="n")
    }
invisible(list(Cg, Cgk))
}

cg = function(x, target, tolerance, ref.interval, facCg, facCgk, n = 0.2, type, col, pch, xlim, ylim, conf.level = 0.95, cex.val = 1.5) {
    old.par <- par(no.readonly = TRUE)                                             #save old par settings
    if (missing(x)) 
        stop("x must be given as a vector")
    if (missing(target)) {
        target = mean(x)
        targetmissing = FALSE
    }
    else targetmissing = TRUE
    if (missing(ref.interval)) 
        ref.interval = pnorm(3) - pnorm(-3)
    sd = sd(x)
    mean = mean(x)
    ref.ar = qnorm(ref.interval, mean, sd) - qnorm(1 - ref.interval, mean, sd)
    if (missing(facCg)) 
        facCg = 0.2
    if (missing(facCgk)) 
        facCgk = 0.1
    if (missing(tolerance)) 
        warning("Missing tolerance! The specification limits are choosen to get Cg = 1")
    if (missing(tolerance)) {
        width = ref.ar/facCg
        tolerance = numeric(2)
        tolerance[1] = mean(x) - width/2
        tolerance[2] = mean(x) + width/2
    }
    quant1 = qnorm((1 - ref.interval)/2, mean, sd)
    quant2 = qnorm(ref.interval + (1 - ref.interval)/2, mean, sd)
    if (length(tolerance) != 2) 
        stop("tolerance has wrong length")
    if (missing(type)) 
        type = "b"
    if (missing(col)) 
        col = 1
    if (missing(pch)) 
        pch = 19
    if (missing(xlim)) 
        xlim = c(0, length(x))
    if (missing(ylim)) 
        ylim = c(min(x, target - n/2 * (abs(diff(tolerance))), quant1, quant2), max(x, target + n/2 * (abs(diff(tolerance))), quant1, quant2))
    Cg = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[1]]
    Cgk = .cg(x, target, tolerance, ref.interval, facCg, facCgk)[[2]]
    layout(matrix(data = c(1, 1, 1, 1, 1, 1, 2, 3, 4), ncol = 3, nrow = 3))
  #  par.old = par()$mar
  #  par(mar = c(5.1, 4.1, 4.1, 10.1))
    cgRunChart(x=x, target=target, tolerance=tolerance, ref.interval=ref.interval, 
               facCg=facCg, facCgk=facCgk, n = n, type=type, col=col, pch=pch, 
               xlim=xlim, ylim=ylim, main="Run Chart", conf.level = conf.level, 
               cex.val = cex.val, cgOut = FALSE)
 #   par(mar = par.old)
    plot(0:6, 0:6, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
    text(2.5, 5, expression(bar(x)), pos = 2, cex = cex.val)
    text(2.5, 5, paste("=", round(mean, 2)), pos = 4, cex = cex.val)
    text(2.5, 4, expression(s), pos = 2, cex = cex.val)
    text(2.5, 4, paste("=", round(sd, 2)), pos = 4, cex = cex.val)
    text(2.5, 3, expression(target), pos = 2, cex = cex.val)
    text(2.5, 3, paste("=", round(target, 5)), pos = 4, cex = cex.val)
    text(2.5, 2, expression(bold(C[g])), pos = 2, cex = cex.val)
    text(2.5, 2, paste("=", round(Cg, 2)), pos = 4, cex = cex.val)
    text(2.5, 1, expression(bold(C[gk])), pos = 2, cex = cex.val)
    text(2.5, 1, paste("=", round(Cgk, 2)), pos = 4, cex = cex.val)
    rect(0, 0, 6, 6)
    par(mar = c(3.1, 4.1, 0.55, 2.1))
    cgHist(x=x, target=target, tolerance=tolerance, ref.interval=ref.interval, 
           facCg=facCg, facCgk=facCgk, n = n, col="lightblue",xlim=xlim, 
           ylim=ylim, main=paste("Histogram of", deparse(substitute(x)), "- target"),
           conf.level = conf.level, cex.val = cex.val, cgOut = FALSE)
#    par(mar = par.old)
    par(mar = c(3.1, 4.1, 1.55, 2.1))
    cgToleranceView(x=x, target=target, tolerance=tolerance, ref.interval=ref.interval, 
                    facCg=facCg, facCgk=facCgk, n = n, type=type, col=col, pch=pch, 
                    xlim=xlim, ylim=ylim, main="Tolerance View", conf.level = conf.level, cex.val = cex.val, cgOut = FALSE)
    par(mar = old.par)
 invisible(list(Cg, Cgk))
} 


#set.seed(123)
#temp=rnorm(125,mean = 10.01 ,sd = 0.1)
#a=cg(temp, target = 10, tolerance = c(8,12),type="b")
#b=.cg(temp, target = 10, tolerance = c(8,12))
#c=cgRunChart(temp, target = 10, tolerance = c(8,12),n = 0.2, type="b", col=4, pch=3,conf.level = 0.95, cex.val = 1,main="Hello")
#d=cgHist(temp, target = 10, tolerance = c(8,12))
#e=cgToleranceView(temp, target = 10, tolerance = c(8,12))