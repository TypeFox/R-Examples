"FileDrawer" <-
function(x, rc=.1)
 {
 k <- length (x$Rxy[!(is.na(x$Rxy))])
 rb <- rbar(x)
 rc <- rc 
 n <- k * (rb/rc - 1)
 mat <- matrix(n)
 colnames(mat) <- c("# of 'lost' studies needed")
 return(mat)
 }

