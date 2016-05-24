rysgran.hist <-
function (data, subset = NULL, which = NULL, ordered=TRUE) 
{
 phi <- data[1, ]
 if (sum(phi) > 550){
   phi<-(-log2(phi/1000))
   phi<-round(phi,digits=1)
 }

 phi <- as.vector(phi, mode = "numeric")
 phi <- factor(phi) 
 Weight <- as.matrix(data[-1,]) 
 colnames(Weight) <- phi
 Sum <- numeric(nrow(Weight))
 for(i in 1:length(Sum)) Sum[i] <- sum(Weight[i,]) 
 m <- matrix(nrow = nrow(Weight), ncol = ncol(Weight))
 for(i in 1:nrow(m)) 
 {
  for(j in 1:ncol(m)) m[i,j] <- Weight[i,j]*100 
 }
 Percent <- matrix(nrow = nrow(Weight), ncol = ncol(Weight))
 for(i in 1:nrow(Percent)) Percent[i,] <- m[i,]/Sum[i] 
 x <- rep(phi, times = nrow(Percent))
 Percent <- t(Percent)
 y <- as.vector(Percent, mode = "numeric")
 g <- rep(rownames(Weight), each = nlevels(x))
 i<-factor(g, levels = rownames(Weight))
 sub <- rep(subset, each = nlevels(x))

 if (is.null(subset)) 
 {
  if (ordered == FALSE)
  {
   bc <- barchart(y ~ x | g,
   horiz = FALSE, origin = 0,
   ylab = "%", xlab = expression(paste(, phi, )), col = "black",
   strip = strip.custom(bg = "grey90"),
   scales = list(x = list(rot = 90)),
   as.table=TRUE,
   panel = function(...)
   {
    panel.grid(h = -1, v = 0, lwd = 1.2)
    panel.barchart(..., border = "transparent")
   })
   return(bc)
  }
  else 
  {
   bc <- barchart(y ~ x | i,
   horiz = FALSE, origin = 0,
   ylab = "%", xlab = expression(paste(, phi, )), col = "black",
   strip = strip.custom(bg = "grey90"),
   scales = list(x = list(rot = 90)),
   as.table=TRUE,
   panel = function(...)
   {
    panel.grid(h = -1, v = 0, lwd = 1.2)
    panel.barchart(..., border = "transparent")
   })
   return(bc)
  }
 } 
 else 
 {
  if (ordered == FALSE)
  {
   bc <- barchart(y ~ x | g, subset = (sub == which),
   horiz = FALSE, origin = 0,
   ylab = "%", xlab = expression(paste(, phi, )), col = "black",
   strip = strip.custom(bg = "grey90"),
   scales = list(x = list(rot = 90)),
   as.table=TRUE,
   panel = function(...)
   {
    panel.grid(h = -1, v = 0, lwd = 1.2)
    panel.barchart(..., border = "transparent")
   })
   return(bc)
  }
  else 
  {
   bc <- barchart(y ~ x | i, subset = (sub == which),
   horiz = FALSE, origin = 0,
   ylab = "%", xlab = expression(paste(, phi, )), col = "black",
   strip = strip.custom(bg = "grey90"),
   scales = list(x = list(rot = 90)),
   as.table=TRUE,
   panel = function(...)
   {
    panel.grid(h = -1, v = 0, lwd = 1.2)
    panel.barchart(..., border = "transparent")
   })
   return(bc)
  } 
 } 
}

