# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
 attribute.default<- function(x, obar.i,  prob.y=NULL, obar = NULL, class = "none", main = NULL,  CI = FALSE,  n.boot = 100, alpha = 0.05,  tck = 0.01, freq = TRUE, pred = NULL, obs = NULL, thres = thres, bins = FALSE, ...){
## attribute plot as displayed in Wilks, p 264.
## If the first object is a prob.bin class, information derived from that.

   old.par <- par(no.readonly = TRUE) # all par settings which
                                      # could be changed.
   # on.exit(par(old.par))

########################################
## if x is an verification object, bootstapping is possible.
   if(CI & class != "prob.bin" )  stop("x must be an 'prob.bin' object create by verify to create confidence intervals" )

########################################   

plot(x, obar.i,  col = 2, lwd = 2, type = "n",
     xlim = c(0,1), ylim = c(0,1),
     xlab =  expression( paste("Forecast probability, ", y[i] ) ),
     ylab = expression( paste("Observed relative frequency, ", bar(o)[1] ))
     )

###################  need to put down shading before anything else.

if(!is.null(obar)){
a   <- (1-obar)/2 + obar
b   <- obar / 2
x.p <- c(obar, obar, 1, 1, 0, 0)
y.p <- c(0, 1, 1, a, b, 0)

polygon(x.p, y.p, col = "gray")

text(0.6, obar + (a-b)*(0.6 - obar), "No skill", pos = 1,
     srt = atan( a - b )/(2*pi)*360 )

}

###########

ii <- is.finite(obar.i)
points(x[ii], obar.i[ii], type = "b", col = 2, lwd = 2)

####### bootstrap CI's
####### this causes a binding error since pred and obs is not introduced.
if(CI){   
n    <- length(pred)
OBAR <- matrix(NA, nrow = length(obar.i), ncol = n.boot)

for(i in 1:n.boot){
  ind      <- sample(1:n, replace = TRUE)
  YY       <- verify(obs[ind], pred[ind], show = FALSE, thresholds = thres, bins = bins)$obar.i    
  OBAR[,i] <- YY
  
} ## close 1:nboot

a<- apply(OBAR,1, quantile, alpha, na.rm = TRUE)
b<- apply(OBAR,1, quantile, 1-alpha, na.rm = TRUE)

for(i in 1:length(a) ){
 lines(rep(x[i], 2), c(a[i], b[i] ), lwd = 1)
 lines( c(x[i] - tck, x[i] + tck), rep(a[i],2),lwd = 1 )
 lines( c(x[i] - tck, x[i] + tck), rep(b[i],2), lwd = 1 )
}

rm(OBAR, a,b)

} ## close if CI

   
   
## plot relative frequency of each forecast
if(freq){
ind<- x< 0.5
text(x[ind], obar.i[ind], formatC(prob.y[ind], format = "f", digits = 3),
          pos = 3, offset = 2, srt = 90)
text(x[!ind], obar.i[!ind], formatC(prob.y[!ind], format = "f", digits = 3),
          pos = 1, offset = 2, srt = 90)
}
if(is.null(main)){title("Attribute Diagram")}else
{title(main)}

abline(0,1)

## resolution line
if(!is.null(obar)){
abline(h = obar, lty = 2)
abline(v = obar, lty = 2)
text( 0.6, obar, "No resolution", pos = 3)
}

invisible()
}
