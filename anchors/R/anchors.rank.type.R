#######################################################################
##
## Function: anchors.rank.type()
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created :  2002-08-14
##
## INPUT:
##   data  : object of class anchors.data
##
##
## OUTPUT:
##   B     : non-parametric values categorizing self-vignette
##   C     : non-parametric values categorizing self-vignette
##           This is an n x 2 matrix where:
##               FIRST column is the minimum value of C for the nth observation
##               SECOND column is the maximum value.
##               WHEN C is scalar-valued, the two columns will be identical
##
##   Cmax    =  2J+1 where J is the number of vignettes; indicates the maximum
##              value contained in C
##
##   n.interval  = Number of cases that have vector valued C (i.e., min and max C not same)
##
##   tfmat      = a logical matrix containing 0/1 for C category by pairwise
##                comparisons between each vignette and the self response.
##
## MODIFIED:
##    2007-09-01 : JW
##    - added B
## 
##    2008-04-20 : JW
##    - was anchors()
##    - now uses anchors.data object instead of processsing data here
##    
#######################################################################


anchors.rank.type <- function(data,type,debug=0) {

  ## extract data...
  z0 <- data$z0
  y0 <- data$y0
  
  ## count things...
  n.self<- NCOL(y0)
  n.vign<- NCOL(z0)
  n.obs <- NROW(z0)
  cmax  <- 2 * n.vign + 1
  bmax  <- n.vign + 1

  ## USAGE verification ....
  if (n.self != 1) {
    stop("You can only specify one self-response!")
  }

  
  ## logic matrix
#  cc <- matrix(FALSE, n.obs, cmax )
#  cc[, 1] <- y0 < z0[, 1]
#  cc[, 2] <- y0 == z0[, 1]
#  cc[, cmax ] <- z0[, n.vign] < y0
#  if (n.vign > 1) {
#    j <- 3
#    for (i in 1:(n.vign - 1)) {
#      v1 <- z0[, i]
#      v2 <- z0[, i + 1]
#      cc[, j]     <- v1 < y0 & y0 <  v2
#      cc[, j + 1] <-           y0 == v2
#      j <- j + 2
#    }
#  }

  if (debug>0) cat("anchors.rank.type: build cc\n")
  cc <- matrix(FALSE, nrow = n.obs, ncol = cmax,
                dimnames = list(rownames(z0), as.character(1:cmax)))
  for (i in 1:n.vign) cc[, 2*i   ] <-  y0 == z0[,i]
  ## MUST test for n.vign > 1
  if (n.vign > 1) {
    for (i in 2:n.vign) cc[,(2*i)-1] <- (y0 <  z0[,i]) & (y0 > z0[,(i-1)])
  }
  cc[,1] <- y0 < z0[,1]
  cc[,cmax] <- y0 > z0[,n.vign]
  
  jidx <- 1:cmax 
  get.range <- function( x ) {
    range(jidx[x])
  }
  
  if (debug>0) cat("anchors.rank.type: span cc\n")
  Cout <- t(apply(cc,1,get.range))
  rownames(cc) <- rownames(Cout) <- rownames(z0)
  colnames(Cout) <- c("Cs","Ce")
  
  ## how many cases have interval C,
  idx.interval <- apply(Cout,1, function(x) {x[1] != x[2] }) 
  n.interval <- sum( idx.interval )

  if (debug>0) cat("anchors.rank.type: build weight\n")
  Cwc <- matrix(as.numeric(cc), nrow=n.obs,ncol=cmax,
                dimnames = list(rownames(z0), as.character(1:cmax)))
  for (i in c(1:n.obs)[idx.interval]) {
    c1 <- Cout[i,1]
    c2 <- Cout[i,2] 
    Cwc[i, c1:c2 ] <- 1/(c2-c1+1)
  }
  

  
  ts <- Cout[,1]
  te <- Cout[,2]
  ## How to calculate B
  ss1 <- as.numeric((ts / 2) - floor(ts/2) > 0)
  Bs <- (floor(ts/2) + ss1)
  Be <- (floor(te/2) + 1)
  if (debug>0) cat("anchors.rank.type: build B\n")
  B <- as.data.frame(list(Bs=Bs,Be=Be))
  
  b.interval <- sum(Bs != Be)

  if (debug>0) cat("anchors.rank.type: build Bweight\n")
  Bwc <- matrix(0,
                nrow=n.obs,ncol=bmax,
                dimnames = list(rownames(z0), as.character(1:bmax)))
  for (i in c(1:n.obs)) {
    c1 <- Bs[i]
    c2 <- Be[i]
    Bwc[i, c1:c2 ] <- 1/(c2-c1+1)
  }


  if (debug>0) cat("anchors.rank.type: rownames\n")
  Cout <- as.data.frame(Cout)
  rownames(Cout) <- rownames(B) <- rownames(z0)
  
# C list currently excludes:  tfmat  =cc
  if (type == "C")
    out <- list(type   = "C",
                span   = Cout, 
                weight = Cwc,
                max    = cmax)
  if (type == "B")
    out <- list(type   = "B",
                span   = B,
                weight = Bwc,
                max    = bmax)

  class(out) <- "anchors.rank.type"
  return(out)
}
