`plot.aspect` <-
function(x, plot.type = "regplot", plot.var = c(1,2), xlab, ylab, main, type, ...)
{
# x         ... object of class "aspect"
# plot.type ... string with values "regplot", "transplot".
# plot.var  ... variables to be plotted (for regression plot only), either numeric or variable names (always vector of length 2) 


var1 <- plot.var[1]
if (is.numeric(var1)) var1 <- colnames(x$data)[var1]
var2 <- plot.var[2]
if (is.numeric(var2)) var2 <- colnames(x$data)[var2]

tab <- table(x$data[,var1], x$data[,var2])         #frequency table of the 2 variables
n <- dim(tab)[1]
m <- dim(tab)[2]


#------------------------------------ regplot -----------------------------------
#draws a before and after regression plot for a table and a set of scores

if (plot.type == "regplot") {

    #par(mfrow = c(1,2))
 
    if (missing(xlab)) xlab1 = var1 else xlab1 <- xlab      
    if (missing(ylab)) ylab1 = var2 else ylab1 <- ylab
    if (missing(main)) main1 <- "Unscaled Solution" else main1 <- main
    if (missing(main)) main2 <- "Scaled Solution" else main2 <- main
    if (missing(type)) type <- "b"
    
  #before lineals plot
    tau <- sum(tab)
    pr <- tab/tau                                     #relative frequencies
    r <- rowSums(pr)                                    #relative row margin
    c <- colSums(pr)                                    #relative column margins

    xave <- as.vector(as.matrix(pr)%*%1:m)/r            
    yave <- as.vector(1:n%*%as.matrix(pr))/c
    z <- c(1:n,1:m) 
    
    dev.new()
    plot(z, z, type = "n", xlab = paste(xlab1," categories"), ylab = paste(ylab1, " categories"), main = main1, 
    xaxt = "n", yaxt = "n", xlim = c(1,n), ylim = c(1,m),...)

    axis(1, at = 1:n, labels = rownames(tab))
    axis(2, at = 1:m, labels = colnames(tab))
    points(1:n, xave, type = type, col = "RED")
    points(yave, 1:m, type= type, col = "BLUE")
    abline(v=1:n, h=1:m, col = "lightgray", lty = 2 )
    for (i in 1:n) text(rep((1:n)[i],m),1:m,as.character(tab[i,]),cex=.8, col = "lightgray")

  #----------- scaled solution------------
    xa <- as.vector(x$catscores[[plot.var[1]]])   
    names(xa) <- rownames(x$catscores[[plot.var[1]]])       
    ya <- as.vector(x$catscores[[plot.var[2]]])
    names(ya) <- rownames(x$catscores[[plot.var[2]]])   

    xave <- as.vector(as.matrix(pr)%*%ya)/r
    yave <- as.vector(xa%*%as.matrix(pr))/c
    z <- c(xa,ya) 

    dev.new()
    plot(z, z, type = "n", xlab = paste(xlab1," scores"), ylab = paste(ylab1," scores"),main = main2,  
    xlim = range(xa), ylim = range(ya),...)
  
    points(xa[order(xa)], xave[order(xa)], type = type, col = "RED")
    points(yave[order(ya)], ya[order(ya)], type = type, col = "BLUE")
    abline(v = xa, h = ya, col = "lightgray", lty = 2)                                           #adds grid
    for (i in 1:n) text(rep(xa[i],m),ya,as.character(tab[i,]),cex=.8, col = "lightgray")
    axis(3, at = xa[order(xa)], labels = names(xa[order(xa)]), cex.axis = 0.6, col.axis = "lightgray", padj = 1)    #adds labels on x-axis (top)
    axis(4, at = ya[order(ya)], labels = names(ya[order(ya)]), cex.axis = 0.6, col.axis = "lightgray", padj = -1)   #adds labels on y-axis (right)

}

#-------------------------------- end regplot -------------------------------


#-------------------------------- transplot ------------------------------------

if (plot.type == "transplot") {
  if (missing(type)) type <- "b"
  if (missing(xlab)) xlab <- "categories"
  if (missing(ylab)) ylab <- "scores"

  for (i in 1:dim(x$data)[2])
  { 
    if (((i-1) %% 4) == 0) {
      dev.new()
      par(mfrow = c(2, 2))
    }
    if (missing(main)) main1 <- paste("Transformation Plot", colnames(x$data)[i]) else main1 <- main

    xa <- x$catscores[[i]]
    n <- length(xa)
    plot(1:n, xa, type = type, xlab = xlab, ylab = ylab, main = main1, xaxt = "n", pch = 1, ...)
    axis(1, at = 1:n, labels = rownames(xa))
    abline(v = 1:n, col = "lightgray", lty = 2)
  }
}

###--------------------------------- end transplot -------------------------------

}

