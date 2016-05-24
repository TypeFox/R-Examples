#######################################################################
##
## Function: anchors.data.count
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created :  2008-04-20
##
## Extracted and refined from former chopit() and anchors()
#######################################################################
anchors.data.count <- function(data , options ) {

  no.self <- is.null(data$y0)
  count <- list()

  count$obs.cat <- unique(c(data$z0,data$y0))
  count$n.cat   <- max(count$obs.cat)
  
  count$nobs.y0  <- NROW(data$y0) # self q
  count$nobs.x0  <- NROW(data$x0) # self predictors
  count$nobs.z0  <- NROW(data$z0) # vign q
  count$nobs.v0s <- NROW(data$v0s) # tau predictors
  count$nobs.v0v <- NROW(data$v0v) # tau predictors

  ## better names...
  count$nobs.self <- count$nobs.y0 
  count$nobs.vign <- count$nobs.z0
  
  if (!no.self) {
    count$nobs.self.vec <- apply(data$y0,2,function(x) {sum(x>0) })
  } else {
    count$nobs.self.vec <- NULL
  }
  count$nobs.vign.vec <- apply(data$z0,2,function(x) {sum(x>0) })
  
  count$nvars.gamma   <- NCOL(data$v0v)
  count$nvars.gamma1  <- NCOL(data$v0v1)
  count$n.self        <- ifelse(is.null(data$y0), 0,  NCOL(data$y0) )

  count$n.vign.set    <- 1 # sum(substring(names(alphad),1,1) == "z")
  ## Q. How do vignettes and self-Q share taus?
  ## A. See NOTES in description of function
  ## number of vignette sets (currently fixed to 1!)
  vign.map      <- NULL
  if (is.null(vign.map))
    vign.map <- 1:count$n.vign.set
  ## Q. How many vignettes in each vignette set
  ## A. Previously loop over each set and counted, now just assume one set
  count$n.vign <- c( NCOL(data$z0) )


  ## 
  if (count$n.self > 0 & options$vign.cut == "hetero")
    stop("if n.self > 0 then only vign.cut == 'homo' allowed\n")

  if (count$n.self > 0) {
    count$n.tau.set <-  count$n.self
  } else if (options$vign.cut == "homo") {
    count$n.tau.set <- count$n.vign.set
  } else if (options$vign.cut == "hetero") {
    count$n.tau.set <- count$n.vign[1]
  } else {
    stop(paste("Invalid option given: vign.cut =",options$vign.cut))
  }

  tmp <- (count$n.cat-1)*count$n.tau.set
  count$tau.start.idx <- seq(1,tmp,by=count$n.cat-1) ## TT

  class(count) <- "anchors.data.count"
  return(count)
}
