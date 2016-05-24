filter.dccm <- function(x, cutoff.cij = 0.4, cmap = NULL, xyz = NULL, fac = NULL, 
                        cutoff.sims = NULL, collapse = TRUE, extra.filter = NULL, ...) {

   # check cij format
   cij <- x
   if("all.dccm" %in% names(cij)) {
      cij <- cij$all.dccm
   } else if(is.list(cij)) {
      cij <- array(unlist(cij), dim = c(dim(cij[[1]]), length(cij)))
   } else if(is.matrix(cij)) {
      cij <- array(cij, dim = c(dim(cij), 1))
   } else if(!is.array(cij)) {
      stop("Input x should be an array/list containing correlation matrices")
   }

   ## Check input is built of simmetric matrices
   if (dim(cij)[1] != dim(cij)[2]) {
     stop("Input 'x' should contain symmetric matrices.")
   }

   ## Check xyz and set cmap   
   if(is.null(cmap)) {
      if(is.null(xyz)) cmap = FALSE
      else cmap = TRUE
   }

   if(cmap) {

      # Inspect cij values with respect to cutoff.cij and contact map

      if(is.null(xyz)) stop("xyz coordinates or a 'pdbs' object must be provided")

      # check factor vector for multiple networks construction
      if(!is.null(fac)) {
         if(!is.factor(fac)) fac = as.factor(fac)
      } else {
         fac = factor(rep("a", dim(cij)[3L])) 
      } 

      # check xyz for contact map calculation
      if(inherits(xyz, "pdbs")) {
         gaps.pos <- gap.inspect(xyz$xyz)
         xyz <- xyz$xyz[, gaps.pos$f.inds]
      }
      if(nrow(xyz) != dim(cij)[3L] && nlevels(fac) > 1)
         stop("xyz matrix doesn't match x. Set fac=NULL for single network construction")
      
      # convert cij to upper.tri matrix for internal use
      pcij <- apply(cij, 3, function(x) x[upper.tri(x)])
    
      ncij <- tapply(1:dim(cij)[3L], fac, function(i) {
         
         # contact map
         if(nlevels(fac) > 1)
             cm <- cmap(xyz[i, ], ...) 
         else
             cm <- cmap(xyz, ...) 
   
         cij.min = apply(abs(pcij[, i]), 1, min)
         cij.max = apply(abs(pcij[, i]), 1, max)
   
         filter <- (cij.min >= cutoff.cij) | (cij.max >= cutoff.cij & cm[upper.tri(cm)]==1)
         
         if(!is.null(extra.filter))
            filter <- filter * extra.filter[upper.tri(extra.filter)] 
         
         ncij <- array(NA, dim=c(dim(cij[,,1]), length(i)))
         for(j in 1:dim(ncij)[3L]) {
            tcij <- cij[,,i[j]]
            tcij[upper.tri(tcij)] <- pcij[, i[j]] * filter
            tcij[lower.tri(tcij)] <- t(tcij)[lower.tri(tcij)]
            ncij[,,j] <- tcij
         }
#         if(length(i) == 1) ncij <- ncij[,,1]
         return(ncij)
      } )
     
      if(collapse) ncij <- lapply(ncij, rowMeans, dims = 2) 
      if(nlevels(fac)==1) ncij <- ncij[[1]]
      if(is.matrix(ncij)) class(ncij) = c("dccm", "matrix")

      return(ncij)

   } else {

      # Filter cijs with cutoff.sims and return mean dccm (dccm.mean())

      if(is.null(cutoff.sims)) cutoff.sims = dim(cij)[3L]
   
      if (cutoff.sims > dim(cij)[3L] || cutoff.sims < 0) {
        stop("The cutoff.sims should be a number between 0 and N, where N is the the number of simulations in the input matrix")
      }

      ## Filter by cutoff.cij and sum across simulations
      cut.cij.inds <- (abs(cij) < cutoff.cij)
      count <- array(NA, dim = dim(cij))
      count[!cut.cij.inds] = 1
      cij.sum <- apply(count, c(1:2), sum, na.rm = TRUE)
    
      ## Mask cij values below cutoff and average across simulations
      cij[cut.cij.inds] = NA
      cij.ave <- apply(cij, c(1:2), mean, na.rm = TRUE)
    
      ## Mask average values if below cutoff.sims
      cut.sims.inds <- (cij.sum < cutoff.sims)
      cij.ave[cut.sims.inds] = 0 ## Could use NA here
    
      class(cij.ave) = c("dccm", "matrix")
      return(cij.ave)
   }
}
