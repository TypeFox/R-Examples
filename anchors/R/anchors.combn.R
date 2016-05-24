#######################################################################
##
## Function: anchors.combn()
## Author  : Jonathan Wand <wand@stanford.edu>
##           http://wand.stanford.edu
##           Daniel Hopkins and Olivia Lau linked in usage of cpolr
## Created :  2003-12-01
##
## MODIFIED:
## - was entropy()
## - no longer  require(combinat) -- now in library(utils)
## - no longer uses genoud
#######################################################################

anchors.combn <- function(adata, fdata, type, options) {
  ## single vignette case: no ties, Cmax=3
  mties <- NULL
  data2 <- adata
  n.vign<- NCOL(adata$z0)
  fdata <- trim.data( data=fdata, anchors=adata)

  if (is.null(adata$formula$cpolr))
    adata$formula$cpolr <- ~ 1
  if (type == "C") 
      fo <- as.formula( paste( "cbind(Cs,Ce) ~ ",  as.character( adata$formula$cpolr)[2] ))
  if (type == "B") 
      fo <- as.formula( paste( "cbind(Bs,Be) ~ ",  as.character( adata$formula$cpolr)[2] ))


  fn <- function(vidx) {
    data2$z0 <- as.matrix(adata$z0[,vidx])
    if (options$debug>0) cat("entropy: do anchors.rank\n")
    ra <- anchors.rank(adata=data2, fdata=fdata, type=type, options=options)
    if (options$debug>0) cat("entropy: fitted/minentropy\n")

    me <- fitted.anchors.rank(ra, ties="minentropy", average=TRUE)
    me <- entropy.calc( me )

    if (options$debug>0) cat("entropy: fitted/cpolr\n")
    cp <- fitted.anchors.rank(ra, ties="cpolr", average=TRUE, unconditional=FALSE)
    cp <- entropy.calc( cp )
    
    return( c( as.numeric(paste(vidx, sep="", collapse="")),
              cp, me,
              ra$summary$n.interval, ra$summary$avg.span, ra$summary$max) )
  }

  ## calculate the entropy for subsets of vignettes
  for (i in 1:n.vign) {
    if (options$verbose) cat("vign:",i,"\n")
    mties <- rbind(  mties, fn(i) )
  }
  for (i in (2:n.vign)) {
    r <- as.matrix(combn(1:n.vign,i))
    for (j in 1:ncol(r)) {
      if (options$verbose) cat("vign:",r[,j],"\n")
      mties <- rbind(mties, fn( r[,j] ) )
    }
  }

  mties <- as.data.frame(mties)
  
  oo <- rev(order(mties[,2]))
  mties <- mties[oo,]
  
  names(mties) <- c("Vignettes",
                    "Estimated entropy",
                    "Minimum entropy",
                    "Interval Cases",
                    "Span avg.",
                    "Max. rank")
  
  rownames(mties) <- 1:nrow(mties)
  class(mties) <- c(class(mties),"anchors.combn")
  return(mties)
}
