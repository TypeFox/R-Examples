interval.plot <-
function(ll, ul)
   {y1 <- ll ; y2<-ul; n <- length(y1)
    plot(y1, type = "n", ylim=c(-.3,.3),xlab = " ",ylab = " ")
    condition <- (ll <= 0 & ul >= 0)
    segments((1:n)[y1<0&y2>0],y1[y1<0&y2>0],(1:n)[y1<0&y2>0],
    y2[y1<0&y2>0])
    segments((1:n)[y1>0],y1[y1>0],(1:n)[y1>0],y2[y1>0],col=17,
    lwd=8)
    segments((1:n)[y2<0],y1[y2<0],(1:n)[y2<0],y2[y2<0],col=17,
    lwd=8)
    SUM<-sum(condition)
    abline(h=0)
    list("Number of intervals that contain 0"=SUM)}

