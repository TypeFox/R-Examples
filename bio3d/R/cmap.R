cmap <- function(...)
  UseMethod("cmap")

cmap.default <- function(...)
  return(cmap.xyz(...))

cmap.xyz <-
function(xyz, grpby=NULL, dcut=4, scut=3, pcut=1, mask.lower = TRUE,
         ncore=1, nseg.scale=1, ...) {

  # Parallelized by parallel package (Mon Apr 22 16:32:19 EDT 2013)
  ncore <- setup.ncore(ncore)
  if(ncore > 1) {
     # Issue of serialization problem
     # Maximal number of cells of a double-precision matrix
     # that each core can serialize: (2^31-1-61)/8
     R_NCELL_LIMIT_CORE = 2.68435448e8
     R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE

     if(nseg.scale < 1) {
        warning("nseg.scale should be 1 or a larger integer\n")
        nseg.scale=1
     }
   }

  if (!(is.numeric(pcut) && pcut >= 0 && pcut <= 1)) {
    stop("Input 'pcut' should a number between 0 and 1")
  }

  xyz=as.xyz(xyz)
  
  if(nrow(xyz)>1) {
     if(is.null(grpby)) {
        nres <- ncol(xyz)/3
     } else {
        inds <- bounds(grpby, dup.inds = TRUE)
        nres <- nrow(inds)
     }
     if(ncore > 1) {
        ni = nrow(xyz)
        RLIMIT = floor(R_NCELL_LIMIT/(0.5*nres*(nres+1)))
        nDataSeg = floor((ni-1)/RLIMIT)+1
        nDataSeg = floor(nDataSeg * nseg.scale)
        lenSeg = floor(ni/nDataSeg)
        cmap.list <- NULL
        for(i in 1:nDataSeg) {
           istart = (i-1)*lenSeg + 1
           iend = if(i<nDataSeg) i*lenSeg else ni
           cmap.list <- c(cmap.list, mclapply(istart:iend, function(j) {
               dmat <- dm.xyz(xyz[j,], grpby, scut, mask.lower=TRUE)
               return(as.numeric(dmat[!lower.tri(dmat)] < dcut))
           }) )
        }
     } else {
        cmap.list <- lapply(1:nrow(xyz), function(j) {
            dmat <- dm.xyz(xyz[j,], grpby, scut, mask.lower=TRUE)
            return(as.numeric(dmat[!lower.tri(dmat)] < dcut))
        }) 
     }
     cmap.t <- rowMeans(do.call(cbind, cmap.list))
     cmap.t <- as.numeric(cmap.t >= pcut )
     cont.map <- matrix(NA, nrow=nres, ncol=nres)
     cont.map[!lower.tri(cont.map)] <- cmap.t
     if(!mask.lower) 
         cont.map[lower.tri(cont.map)] <- t(cont.map)[lower.tri(cont.map)]

  } else {

     ## Distance matrix (all-atom)
     dmat <- dm.xyz( xyz, grpby, scut, mask.lower = mask.lower, ncore=ncore)
     ## Contact map
     return(matrix(as.numeric(dmat < dcut),
                ncol = ncol(dmat),
                nrow = nrow(dmat)))

  }
  return (cont.map)
}

