# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
        # ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
conditional.quantile<- function(pred, obs, bins = NULL, thrs = c(10, 20), main = "Conditional Quantile Plot", ...){
## creates conditional probability plots as described in Murphy et al.
#set.seed(1)

#frcst <- round(runif(100, 20, 70))
#obs<- rnorm( 100, frcst, 10)
#bins <- seq(0,100,10)
#thrs<- c( 10, 20) # number of obs needed for a statistic to be printed #1,4 quartitle, 2,3 quartiles

  old.par <- par(no.readonly = TRUE) # all par settings which                                    # could be changed.
   on.exit(par(old.par))
  
# check bin sizes; issue warning
if(!is.null(bins)){
if( min(bins)> min (obs) | max(bins)< max(obs) ){warning("Observations outside of bin range. \n")}
if( min(bins)> min (pred) | max(bins)< max(pred) ){warning("Forecasts outside of bin range. \n")}
} else {
dat <- c(obs,pred); min.d <- min(dat); max.d <- max(dat)
  bins<- seq(floor(min.d), ceiling(max.d), length = 11)

}   ## close bin check  

## plot ranges
  lo<- min(bins); hi<- max(bins)
  
## if selected, the quasi-continuous data is subsetted into larger
## bins so that quantile statistics might be calculated.

  b<- bins[- length(bins)]
  labs<- b + 0.5*diff(bins)

obs.cut<- cut(obs, breaks = bins, include.lowest = TRUE, labels = labs)
obs.cut[is.na(obs.cut)]<- labs[1] # place anything below the limit into first bin. 
obs.cut<- as.numeric(as.character(obs.cut))

frcst.cut<- cut(pred, breaks = bins, include.lowest = TRUE, labels = labs)
frcst.cut[is.na(frcst.cut)]<- labs[1]
frcst.cut<- as.numeric(as.character(frcst.cut))

## calculate stats ext

n<- length(labs)

lng<- aggregate(obs, by = list(frcst.cut),length)
med<- aggregate(obs, by = list(frcst.cut),median)

q1 <- aggregate(obs, by = list(frcst.cut),quantile, 0.25)
q2 <- aggregate(obs, by = list(frcst.cut),quantile, 0.75)

q1$x[lng$x <= thrs[1]] <- NA
q2$x[lng$x <= thrs[1]] <- NA

q3 <- aggregate(obs, by = list(frcst.cut),quantile, 0.1)
q4 <- aggregate(obs, by = list(frcst.cut),quantile, 0.9)

q3$x[lng$x <= thrs[2]] <- NA
q4$x[lng$x <= thrs[2]] <- NA

par( mar = c(5,5,5,5) )
  
plot(frcst.cut, obs.cut, xlim = c(lo,hi), ylim = c(lo, hi), main = main,
     type = 'n', ylab = "Observed Value", xlab = "Forecast Value", ... )
mtext("Sample Size", side = 4, adj = -1)

#### legend
legend.txt<- c("Median", "25th/75th Quantiles", "10th/90th Quantiles")

legend(min(pred) + 0.55*diff(range(pred)),
       min(obs) + 0.25*diff(range(obs)), legend.txt, col = c(2,3,4),
       lty = c(1,2,3), lwd = 3, cex = 0.7 )

abline(0,1)
X <- as.numeric(as.character(med$Group.1))
  
lines(X, med$x, col = 2, lwd = 3)
lines(X, q1$x,
      col = 3, lty = 2, lwd = 3)
lines(X, q2$x,
      col = 3, lty = 2, lwd = 3)
lines(X, q3$x,
      col = 4, lty = 3, lwd = 3)
lines(X, q4$x,
      col = 4, lty = 3, lwd = 3)

pp<- par("plt")

par("plt" = c(pp[1], pp[2], pp[3], 0.2))

par(new = TRUE)

hist(frcst.cut, breaks = bins, col = "blue",
     main = "", axes = FALSE, xlim = c(lo, hi),
     xlab = " " , ylab = " ")
axis(4, line = 0)
# mtext("Sample Size", side = 4, line = 1)

}
