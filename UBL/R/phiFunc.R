#==========================================================================#
# phi.R                                                                    #
#                                                                          #
#                                                                          #
# Rita P. Ribeiro  Sept 2011                                               #
#--------------------------------------------------------------------------#
## Objective: relevance definition                                         #
#==========================================================================#

# ======================================================================
# The phi function specifies the regions of interest in the target
# variable. It does so by performing a Monotone Cubic Spline
# Interpolation over a set of maximum and minimum relevance points. 
# The notion of RELEVANCE can be associated with rarity.
# Nonetheless, this notion may depend on the domain experts knowledge.
# ======================================================================

## keep up to date (in R and C)
phiMethods <- c("extremes","range")

# ======================================================================
# phi.setup
# This function does the control of loss parameters
# ======================================================================
phi.setup <- function(y, method = phiMethods,
                      extr.type = NULL, coef=1.5, control.pts = NULL) {

  method <- match.arg(method, phiMethods)

   
  control.pts <- do.call(paste("phi",method,sep="."),
                         c(list(y=y), extr.type = extr.type,
                           list(control.pts=control.pts),coef = coef))
  
  
  list(method = method, 
       npts = control.pts$npts, control.pts = control.pts$control.pts)
}

phi.control <- function(y, method="extremes", extr.type="both", coef=1.5, control.pts=NULL) {

  call <- match.call()
  phiP <- phi.setup(y, method, extr.type, coef, control.pts)
  
  phiP
  
}


## ======================================================================
## phi.extremes
## EXTREMES: 
## As we know that S1 and S2 have the same sign, the slope
## of the adjacent values will be the weighted harmonic mean between
## the median and the more extreme outlier.
## adjH = min{ y | y >= Q3 + coef * IQR}
## adjL = max{ y | y <= Q1 - coef * IQR}
## For adjL and adjH, coef = 1.5
## but for extreme outliers we can assign coef = 3
## stats =  lower whisker, the first quartile, the median,
##          the third quartile and the extreme of the upper whisker
## ======================================================================
phi.extremes <- function(y, extr.type = c("both","high","low"), control.pts,
                         coef=1.5) {
  extr.type <- match.arg(extr.type)
  
  control.pts <- NULL
  
  extr <- boxplot.stats(y,coef=coef)

  r <- range(y)
  
  if(extr.type %in% c("both","low") &&
     any(extr$out < extr$stats[1])) {

    ## adjL
    control.pts <- rbind(control.pts,c(extr$stats[1],1,0))
  
  } else {

    ## min
    control.pts <- rbind(control.pts,c(r[1],0,0))
  }

  ## median 
  if(extr$stats[3]!= r[1]){
    control.pts <- rbind(control.pts,c(extr$stats[3],0,0))
  }         
  if(extr.type %in% c("both","high") &&
     any(extr$out > extr$stats[5])) {
       
    ## adjH
    control.pts <- rbind(control.pts,c(extr$stats[5],1,0))
  
  } else {

    ## max
    if(extr$stats[3] != r[2]){
      control.pts <- rbind(control.pts,c(r[2],0,0))
    }
  }

  npts <- NROW(control.pts)

 
  list(npts = npts,
       control.pts = as.numeric(t(control.pts)))##,
           
}


## ======================================================================
## phi.range
## ======================================================================
phi.range <- function(y, extr.type, coef, control.pts, ...) {
  
  ## if it comes from pre-set env
  if(!is.null(names(control.pts))) 
    control.pts <- matrix(control.pts$control.pts,nrow=control.pts$npts,byrow=T)
  
  extr.type <- NULL
  coef <- NULL
  if(missing(control.pts) || !is.matrix(control.pts) ||
     (NCOL(control.pts) > 3 || NCOL(control.pts) < 2))
    stop('The control.pts must be given as a matrix in the form: \n',
         '< x, y, m > or, alternatively, < x, y >')

  npts <- NROW(control.pts)
  dx <- control.pts[-1L,1L] - control.pts[-npts,1L]
  
  if(any(is.na(dx)) || any(dx == 0))
    stop("'x' must be *strictly* increasing (non - NA)")
  
  if(any(control.pts[,2L] > 1 | control.pts[,2L] < 0))
    stop("phi relevance function maps values only in [0,1]")
  
  control.pts <- control.pts[order(control.pts[,1L]),]
  
  if(NCOL(control.pts) == 2) {
    
    ## based on "monoH.FC" method
    dx <- control.pts[-1L,1L] - control.pts[-npts,1L]
    dy <- control.pts[-1L,2L] - control.pts[-npts,2L]
    Sx <- dy / dx
    
    ## constant extrapolation
    m <- c(0, (Sx[-1L] + Sx[-(npts-1)]) / 2, 0)
    
    control.pts <- cbind(control.pts,m)
    
  }
  
 
  r <- range(y)
  npts <- NROW(matrix(control.pts,ncol=3))

 
  list(npts = npts,
       control.pts = as.numeric(t(control.pts)))##,

}
