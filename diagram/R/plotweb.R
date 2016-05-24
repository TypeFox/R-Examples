
##==============================================================================
## plotweb      : plots a web
###==============================================================================

plotweb    <- function (flowmat, names=NULL, lab.size = 1.5, add  = FALSE,
     fig.size = 1.3, main = "", sub = "", sub2 = "",  log = FALSE,
     mar=c(2,2,2,2), nullflow = NULL, minflow = NULL, maxflow = NULL,
     legend=TRUE, leg.digit=5, leg.title=NULL, lcol= "black", arr.col= "black",
     val=FALSE, val.digit=5, val.size=0.6, val.col="red", val.title=NULL,
     val.ncol=1, budget = FALSE, bud.digit=5, bud.size=0.6,
     bud.title="budget", bud.ncol=1, maxarrow = 10, minarrow = 1,
     length=0.1, dcirc=1.2, bty = "o", ...)  {
   
  ##-----------------------------------------------------------------
  ## constructing the names,  flow matrix
  nm <- par("mar")
  if (ncol(flowmat) != nrow(flowmat))
    stop("flowmat has to be square")

  components <- names
  if (is.null(components))
    components <- colnames(flowmat )
  if (is.null(components))
    components <- rownames(flowmat )
  if (is.null(components))
    components <- as.character(1:ncol(flowmat))
  
  numcomp    <- length(components)  ## number of food web components
  if (ncol(flowmat) != numcomp)
    stop("flowmat and names not compatible")

  if (length(arr.col) ==1)
    arr.col  <- matrix(nrow=numcomp, ncol=numcomp, arr.col)

  flowmatrix <- flowmat
  if (!is.null(nullflow )) {
    flowmatrix[flowmatrix<nullflow[1] ] <- 0
    if (length(nullflow == 2))
      flowmatrix[flowmatrix>nullflow[2] ] <- 0
  }
  zero <- 0
  if (log) {
    flowmatrix <- log10(flowmatrix +1e-20)
    flowmatrix[flowmatrix==-20]<-0
  }

  if (is.null(maxflow))
    maxflow  <- max(flowmatrix)
  else if (log)
    maxflow <- log10(maxflow)
  if (is.null(minflow))
    minflow  <- min(flowmatrix[flowmatrix != zero])
  else if (log)
    minflow <- log10(minflow)

  ##-----------------------------------------------------------------
  ## new empty plot

  if (! add)  {

    figlim  <- c(-fig.size, fig.size)
    marg    <- par("mar")

    # if values written:shift left
    if (val)  mar <- mar + c(0,-2,0,2)
    mar <- pmax(mar,0)
#  was:  ifelse(val,  mar <- c(1, 0, 1, 2), mar<-c(1, 1, 1, 1))
        par(mar=mar)
        plot(c(0,  0),  type = "n",  ylab = "",  asp = 1,  xaxt = "n",
         yaxt = "n",  frame.plot = FALSE,  xlim = figlim,
         ylim = figlim, main=main, xlab="")
    mtext(side=3, line=-1, sub)
    mtext(side=1, adj=0.5, text=sub2)
  }
  
  ##-----------------------------------------------------------------  
  ## component labels positioned on a circle 
  
  alpha0 <- pi/2
  alpha  <- alpha0 - (1:numcomp) * 2 * pi/numcomp  # all angles
  xl     <- cos(alpha)                             # all x-positions
  yl     <- sin(alpha)                             # all y-positions

  # adjust depending on 
  for (i in 1:numcomp) {                           # write labels
		    if (    xl[i]  > 0     ) adjustx = 0       # text adjustments
		    if (    xl[i]  < 0     ) adjustx = 1
		    if (abs(xl[i]) < 0.0001) adjustx = 0.5
  	    if (    yl[i]  > 0     ) adjusty = 0
		    if (    yl[i]  < 0     ) adjusty = 1
		    if (abs(yl[i]) < 0.0001) adjusty = 0.5
        text(xl[i],  yl[i],  components[i],
              adj =c(adjustx, adjusty),  cex = par("cex") * lab.size)
                        }

  ##-----------------------------------------------------------------
  ## arrows representing the flows

	circle <- function (i, lwd, col) 	 {  # circular arrow
    cx   <- xl[i]*dcirc
    cy   <- yl[i]*dcirc
    r    <- 0.1  # radius
    x    <-c(seq(-pi, pi, by=0.01), pi)  # adding the last element ensures circle is closed
    lines(cx+r*sin(x), cy+r*cos(x), lwd=lwd, col=col)
  }

  par(lend=1)

  darrow   <- (maxarrow-minarrow)/(maxflow-minflow)
  dr       <- 0.02 

  xi       <- xl-dr* cos(alpha)   # arrow positions: shifted relative to labels
  yi       <- yl-dr* sin(alpha)
  iflow    <- 1
  offset   <- 1
  ltext    <- NULL
  for (i in 1:numcomp)  {

    x2 <- xi[i]
    y2 <- yi[i]
    for (j in 1:i)  {
      if (flowmatrix[i, j] >zero | flowmatrix[j, i] >zero) {
        Arr.col <- arr.col[i,j]
        x1 <- xi[j]
        y1 <- yi[j]
        dx <- x2-x1
        dy <- y2-y1

        ifelse (i == j, fsize<-flowmatrix[i, j],
                        fsize<-flowmatrix[i, j]-flowmatrix[j, i])

        if (fsize>0) {
           code <- 1
        } else {
           code<-2
           Arr.col <- arr.col[j,i]
        }
          
        size <- minarrow + darrow * (abs(fsize)-minflow) # arrow thickness
                              
 		    if (i != j)
           arrows (x1+dr*dx, y1+dr*dy, x2-dr*dx, y2-dr*dy,
                  length=length, code=code, lwd=size, col=Arr.col, ...)
			  if (i == j)
           circle (i, lwd=size, col=Arr.col)

        if (val)  {
          text(x=(x1+x2)*0.5, y=(y1+y2)*0.5, labels=iflow, offset=offset, col=val.col)
          ltext<- c(ltext, paste(iflow,  ":",  format.pval(abs(fsize), val.digit)))
        }
        iflow <- iflow+1
      }  # end flowmatrix>zero
    }  # end j
  }  # end i

  ##-----------------------------------------------------------------
  ## legends

  if (legend)  {
    sizeleg = par("cex") * lab.size

    ## size of largest and smallest arrows

    if (!log) {
      tmax <-    maxflow
      tmin <-    minflow
      title=leg.title
    } else  {
      tmax <- 10^maxflow
      tmin <- 10^minflow
      title=paste("logarithmic scale", leg.title)
    }

    legend("bottomright", legend=c(format.pval(tmax, leg.digit),
           format.pval(tmin, leg.digit)), cex=sizeleg, title=title,
           lwd=c(maxarrow, minarrow), bty=bty)
  }
  if (!val & !budget)
    return

  if (! add)   {
    par(mar=c(0, 0, 0, 0))
    par(new=TRUE)
    plot(c(0,  0),  type = "n",  ylab = "",  xaxt = "n",
         yaxt = "n",  frame.plot = FALSE, main="", xlab="")
  }
    
  if (val)
    legend("topright", legend=ltext, cex=val.size,
           title=val.title, ncol=val.ncol, bty=bty)

  if (budget) {
   rate <- NULL
    # sum all flows in - sum of flow out 
   for (i in 1:numcomp)
      rate <- c(rate, paste(components[i],  ":",
          format.pval(sum(flowmat[, i])-sum(flowmat[i, ]), bud.digit)))
   legend("topleft", legend=rate, cex=bud.size, title=bud.title,
          ncol=bud.ncol, bty=bty)
  }
      par("mar"=nm)
}
