"powerplot" <-
function(expr="x^2",xlab="x",ylab="y"){
   invtxt <- switch(expr, "x^2"="sqrt(y)","x^4"="y^0.25",
  "exp(x)"="log(x)","sqrt(x)"="y^2","x^0.25"="y^4",
  "log(x)"="exp(y)")
    x <- (1:60)/6
    y <- eval(parse(text=expr))
    form <- formula(paste("~",expr))
    dy <- deriv(form,"x")
    x0 <- min(x)+diff(range(x))*0.4
    y0 <- eval(parse(text=expr),list(x=x0))
    b <- eval(dy, list(x=x0))
    b <- attr(b,"gradient")
    plot(x, y, type = "n", xlab = "", ylab = "")
    lines(x,y,type="l",lwd=2,col=2)
    chh <- par()$cxy[2]
    theta <- atan(b*diff(par()$usr[1:2])/diff(par()$usr[3:4]))*180/pi
    mtext(side=1,line=2.5, xlab, cex=1)
    mtext(side=2,line=2.5, ylab, cex=1)
    funexpr <- parse(text=paste("y ==",expr))
    text(x0, y0+chh/2, funexpr, srt=theta,cex=1.5)
    invexpr <- parse(text=invtxt)[[1]]
    titletxt <- substitute(paste(tx, tilde(y) == invexpr),
    list(tx="Replace y by ", invexpr=invexpr))
    mtext(side=3,line=0.5,titletxt)
}
