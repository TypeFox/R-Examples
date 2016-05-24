# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
"discrimination.plot" <- function(group.id, value, breaks = 11, main =
"Discrimination Plot", xlim = NULL, ylim = NULL, legend = FALSE,
leg.txt = paste("Model", sort(unique(group.id)) ),  marginal =
TRUE, cols = seq(2, length(unique(group.id)) + 1 ), xlab = "Forecast", ... ){
dat <- data.frame(group.id = group.id, value = value)
  groups <- sort(unique(group.id))  
n.group <- length(groups )

# test data  
# group.id <- dat$id2
# value     <- dat$allmaxsev
# breaks   <- 11
# main     <- "Comparison of Distributions"
# leg.txt  <- paste("Model", unique(dat$id2) )
# marginal <- TRUE
# cols     <-  seq(2, length(unique(dat$id2)) + 1 )

  old.par <- par(no.readonly = TRUE) # original parameters
  on.exit(par(old.par))

BRKS<- seq(min(value), max(value), length = breaks)
OUT <- matrix(NA, nrow = (breaks - 1), ncol = n.group)

for( i in 1:n.group){
  
XX     <- hist(value[group.id == groups[i] ], plot = FALSE, breaks = BRKS)
OUT[,i]<- XX$counts/sum(XX$counts)
}
## limits for plots

mx.1 <- max(value)
mn.1 <- min(value)
mx.2 <- max(OUT)
mn.2 <- min(OUT)

if(!is.null(xlim) ){
mx.1 <- xlim[2]
mn.1 <- xlim[1]
}

if(!is.null(ylim) ){
mx.2 <- ylim[2]
mn.2 <- ylim[1]
}

if(marginal){
par(oma = c(0,0,2,0))

layout(matrix(1:2, nrow = 2), heights = c(1,4) )
if(legend){par(mar = c(0,4,1,9) )} else
 {par(mar = c(0,4,1,1) ) }

boxplot(value~group.id, data = dat, horizontal = TRUE,  axes = FALSE,
        col = cols , ylim = c(mn.1, mx.1 ), ... )
axis(side = 2, at = 1:n.group, labels = leg.txt, las = 2 )

if(legend){par( mar = c(4,4,0,9))} else
{par(mar = c(4,4,0,1) )}
}else{if(legend){
  par(mar = c(4,4,4,8) )
} ## close if legend
                                        # par(mar = c(4,4,3,1) )
            }  ## close if marginal


plot(XX$mids, apply(OUT, 1, max) , type = "n", xlab = xlab, ylab
     = "Relative Frequency", ylim = c(0, mx.2), xlim = c(mn.1, mx.1 ), ...  )

for(i in 1:n.group){
points(XX$mids, OUT[,i], type = "b", col = cols[i], pch = 14+i )
}
if(marginal){mtext(main, outer = TRUE, ...)}else{title(main)}

abline(h = 0); abline(v=0)

if(legend){
  par(xpd = NA)
  xx <- mx.1 + 0.1*(mx.1 - mn.1)
  yy <- mean(c(mx.2, mn.2))
  legend(x= xx, y = yy, yjust = 0.5,  legend = leg.txt, col= cols, pch = seq(15,
         , 1, n.group), lty=1, merge=TRUE, cex = 0.8)
}# close if legend
}
