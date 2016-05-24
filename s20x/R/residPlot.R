residPlot<-
function (lmfit, f = 0.5) 
## Fitted Values versus Residuals plot
{
    yhat <- fitted(lmfit)
    res <- resid(lmfit)
    y = yhat + res
    QQ = qr.Q(lmfit$qr)
    RR = qr.R(lmfit$qr)
    xx = QQ %*% RR
    newdata = data.frame(y, yhat^2, xx[, -1])

    xnam <- dimnames(newdata)[[2]][-1]
    form <- as.formula(paste("y ~", paste(xnam, collapse = "+")))
    foo <- lm(form, data = newdata)
    goo <- summary(foo)
    curvp <- goo$coefficients[dimnames(goo$coefficients)[[1]] == 
        "yhat.2", 4]
   ## sort the data so we can fit it okay
   
   yhat.sort=sort(yhat,index=TRUE)     
    newyhat=yhat.sort$x
    ix=yhat.sort$ix
    newres=res[ix]
    
 ## now get trend stuff from these sorted values
 
    trend <- lowess(newyhat, newres, f)
    e <- (newres - trend$y)^2
    scatter <- lowess(newyhat, e, f)
    
    uplim <- trend$y + sqrt(abs(scatter$y))
    lowlim <- trend$y - sqrt(abs(scatter$y))
    plot(yhat, res, pch = 1, xlab = "Fitted Values", ylab = "Residuals", 
        main = paste("Resids vs. Fitted ~ Test for Quadratic (p=", 
            signif(curvp, 3), ")"))
    lines(trend, col = "Blue")
    lines(scatter$x, uplim, lty = 2, col = "Red")
    lines(scatter$x, lowlim, lty = 2, col = "Red")

}

