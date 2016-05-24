#' Add markers to a framework map using 3-point likelihoods
#' 
#' Based on a framework map and chromosome assignments, markers are inserted at the midpoints of intervals. Positions are chosen by maximizing the 3 point likelihood. 
#' @export
#' @useDynLib mpMap
#' @param mpcross Object of class \code{mpcross}, with framework markers
#' @param newmpcross Object of class \code{mpcross}, with additional markers
#' @param newchr Set of chromosome listings for additional markers
#' @param mapfx Map function to use - default is haldane
#' @return Combination of the two mpcross objects with new markers inserted at midpoints of intervals with maximum LOD values if LOD > 3. Note that if markers are inserted at the same midpoint, it will be necessary to locally reorder them and then re-estimate the length of the map by summing adjacent recombination fractions.
#' @details Note that the values in \code{newchr} need to correspond to the chromosomes already existing in \code{mpcross}.
#' @seealso \code{\link[mpMap]{mporder}}, \code{\link[mpMap]{mpadd}}
add3pt <- function(mpcross, newmpcross, newchr, mapfx=c("haldane", "kosambi"))
{
  if (missing(mapfx)) mapfx="haldane"
  if (mapfx=="kosambi") mf <- kosambiX2R  else mf <- haldaneX2R

  if (length(newchr)!=ncol(newmpcross$finals)) 
	stop("Must have a chromosome assignment for every marker in newmpcross")

  if (is.character(newchr)) newchr <- match(newchr, names(mpcross$map))

  if (nrow(newmpcross$finals) != nrow(mpcross$finals))
	stop("mpcross objects are not compatible")

  n.founders <- nrow(mpcross$founders)
  nobs <- nrow(mpcross$finals)
  fin <- newmpcross$finals 
  fio <- mpcross$finals 
  fon <- newmpcross$founders
  foo <- mpcross$founders

  lodmax <- list()
  newmap <- list()
  for (ii in 1:length(mpcross$map))
  {
    chrmap <- mpcross$map[[ii]]
    nmrk <- length(which(newchr==ii))
    if (nmrk>0) {
     lmap <- length(chrmap)
     indn <- which(newchr==ii)
     mrknam <- colnames(fin)[indn]
     indo <- match(names(chrmap), colnames(foo)) 
     fouchr <- cbind(fon[,indn], foo[,indo])
     finchr <- cbind(fin[,indn], fio[,indo])

     indices <- match(names(chrmap), colnames(finchr))
     dmap <- diff(chrmap)
     pos <- c(-5, chrmap[1:(lmap-1)]+dmap/2, chrmap[lmap]+5)
     
     left <- c(-1, indices[c(1:(lmap-1), lmap-1)])
     mid <- c(indices[1], rep(-1, lmap-1), indices[lmap])
     right <- c(indices[c(2, 2:lmap)], -1) 

     rl <- mf(c(5, dmap/2, dmap[lmap-1])/100)
     rr <- mf(c(dmap[1], dmap/2, 5)/100) 

     if (mapfx=="kosambi")
	ro <- (rl+rr)/(1+4*rl*rr) else ro <- rl+rr-2*rl*rr

     #### note, we also need to compute the base lkhd to get LOD values
     rl <- c(rl, rep(.5, length(left)-1), mf(dmap[lmap-1]/100))
     rr <- c(rr, mf(dmap[1]/100), rep(.5, length(left)-1))
     ro <- c(ro, .5, mf(dmap/100), .5)
     left <- c(left, left)
     mid <- c(mid, mid)
     right <- c(right, right)
   
     ### loop through each marker which is assigned to this chromosome
     chrmrk <- match(mrknam, colnames(finchr))

     finals <- cbind(as.numeric(as.character(mpcross$id)), finchr)
     finals[is.na(finals)] <- -1

     finvec <- as.vector(t(finals))
     fouvec <- as.vector(t(fouchr))

     npos <- length(left)

     lk3pt <- .C("add3pt",  as.integer(finvec), as.integer(fouvec), as.integer(mpcross$pedigree[,1]), as.integer(mpcross$pedigree[,2]), as.integer(mpcross$pedigree[,3]), as.integer(nobs), as.integer(n.founders), as.integer(nmrk), as.integer(nrow(mpcross$pedigree)), as.double(rl), as.double(rr), as.double(ro), as.integer(left), as.integer(mid), as.integer(right), as.integer(ncol(finchr)), as.integer(npos), out=integer(length=nmrk), lod=double(length=nmrk*npos/2), PACKAGE="mpMap")

     lod <- matrix(lk3pt$lod, nrow=nmrk, byrow=TRUE)
     lod[is.nan(lod)] <- 0
     if (nmrk>1) lodmax[[ii]] <- apply(lod, 1, max) else lodmax[[ii]] <-  max(lod)
     if (nmrk>1) wlod <- apply(lod, 1, which.max) else wlod <- which.max(lod)
     ### now take each cts marker and place it at the location indicated
     newmrk <- pos[wlod]
     names(newmrk) <- colnames(fon)[indn]
     newmrk <- newmrk[which(lodmax[[ii]]>3)]
     newmap[[ii]] <- sort(c(chrmap, newmrk))
     newmap[[ii]] <- newmap[[ii]] - min(newmap[[ii]])
     } else newmap[[ii]] <- chrmap
  } ## end of loop over chromosomes

  out <- mpadd(mpcross, newmpcross)
  out$oldmap <- mpcross$map
  names(newmap) <- names(mpcross$map) 

  class(newmap) <- "map"
  out$map <- newmap
  out <- maporder(out)
  out <- subset(out, markers=unlist(lapply(newmap, names)))
 
  return(out)
}
