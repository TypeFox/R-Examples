lines.attrib <- function(x, ...){
assign("y", x$obar.i)
assign("x", x$y.i)
do.call("lines", list(x, y, ...) )

}

####
lines.roc <- function(x, binormal = FALSE, ... ){

A <- roc.int(x$obs, x$pred, x$thres, binormal = binormal)
A<- as.data.frame(A)
if(binormal){
  dat  <- A
  names(dat) <- c("thres", "proby", "probn", "zH", "zF")
  dat <- dat[is.finite(dat$zH) & is.finite(dat$zF), ] ## reduce dat, get rid of nans and inf   
  new <-  as.data.frame( matrix(qnorm(seq(0.005, 0.995, 0.005 ) ), ncol = 1) )
  names(new) <- "zF"
  A <- lm(zH ~ zF, data = dat)$fitted.values
  B <- predict(lm(zH ~ zF, data = dat), newdata = new)

lines(pnorm(new$zF), pnorm(B), ...)
} else {lines(A$F, A$H, ...)}

}
