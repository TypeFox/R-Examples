#' Subset mpcross object
#'
#' Reduces an mpcross object down to a specified set of chromosomes, markers and/or lines
#' @S3method subset mpcross
#' @method subset mpcross
#' @param x Object of class \code{mpcross}
#' @param chr Selected chromosomes TO KEEP. Default is all
#' @param markers Selected markers TO KEEP. Default is all
#' @param lines Selected lines TO KEEP. Default is all
#' @param ... Additional arguments
#' @note Chromosomes can be input either as the character names of chromosomes or the index of the chromosomes in the map. Markers can be input as character names or the index in the matrix x$finals. Lines can be input as either character values (matching the rownames of x$finals) or indices of rows in that matrix. 
#' @return The original object with chromosomes/lines/markers removed which are not listed in the arguments.
#' @seealso \code{\link[mpMap]{mpcross.object}}
#' @examples
#' map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=map, pedigree=ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' sim.dat
#' red.dat <- subset(sim.dat, chr=1, lines=1:50)
#' red.dat

subset.mpcross <-
function(x, chr=NULL, markers=NULL, lines=NULL, ...)
{
  if (all(sapply(c(chr, markers, lines), length)==0)) return(x)
 
  output <- x

  if (!is.null(chr)) {
    if (is.numeric(chr)) chr <- names(x$map)[chr]
    output$map <- as.list(output$map[chr])

  for (ii in chr)
    markers <- unique(c(markers, names(x$map[[ii]])))
  } 

  if (!is.null(markers)) {
    if (is.character(markers)) 
      mrknum <- match(markers, colnames(x$finals)) else mrknum <- markers

    if (is.numeric(markers)) {
	    mrknum <- markers
	    markers <- colnames(x$finals)[mrknum]
    }

    output$founders <- as.matrix(output$founders[,mrknum])
    output$finals <- as.matrix(output$finals[,mrknum])
    colnames(output$founders) <- colnames(output$finals) <- markers

    if (!is.null(x$rf)) {
	m2 <- match(markers, colnames(output$rf$theta))
	output$rf$theta <- output$rf$theta[m2, m2]
	output$rf$lkhd <- output$rf$lkhd[m2, m2]
	output$rf$lod <- output$rf$lod[m2, m2]
    }
    if (!is.null(x$ld)) {
	m2 <- match(markers, colnames(output$ld$W))
	output$ld$W <- output$ld$W[m2, m2]
	output$ld$r2 <- output$ld$r2[m2, m2]
	output$ld$LewontinD <- output$ld$LewontinD[m2, m2]
	output$ld$delta2 <- output$ld$delta2[m2, m2]
    }
   ### should also be trimming down the map. 
   ### leave markers at estimated positions
   if (!is.null(x$map)) 
 	output$map <- lapply(x$map, function(y) return(y[which(names(y) %in% markers)]))
   output$map[which(lapply(output$map, length)==0)] <- NULL

    if (!is.null(x$lg)) {
	x$lg$groups[mrknum] <- NA
	x$lg$n.groups <- length(table(x$lg$groups))	
    }
  }

  if (!is.null(lines)) {
    if (is.character(lines)) linnum <- match(lines, rownames(x$finals)) else linnum <- lines
      output$finals <- output$finals[linnum,]
    if (!is.null(output$pheno))
      output$pheno <- output$pheno[linnum,]
  }
  output
}

