# radius of gyration
# xyz: length 3N, 1D array of atomic coordinates,
#      M*3N matrix, or a list object containing xyz
# mass: length N 1D array of atomic masses [amu],
#      or a PDB object having masses stored in the 
#      "B-factor" column
rgyr <- function(xyz, mass=NULL, ncore=1, nseg.scale=1)
{
   # Parallelized by parallel package
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
    
   # check xyz, vector, matrix, or list
   if(is.list(xyz)) 
      xyz <- xyz$xyz
   if(is.vector(xyz))
      xyz <- matrix(xyz, nrow=1)
   
   #check mass array and load masses
   if(is.list(mass)) 
     mass <- mass$atom[,"b"]
   if(is.null(mass))
     # assume carbon
     mass <- rep(12.01, ncol(xyz)/3)

   if(ncol(xyz)/3 != length(mass)) 
      stop("The length of masses doesn't match the number of atoms") 
   
   rg <- function (xyz, mass) {
      nAtom <- length(mass)
      mc <- matrix(xyz, 3, nAtom)
      v <- replicate(3, mass) * t(mc)

      com <- colSums(v)/sum(mass)

      recenteredpos <- mc - replicate(nAtom, com)
      rog_sq <- sum(colSums(recenteredpos**2) * mass)
      rog_sq <- rog_sq / sum(mass)
      return( sqrt(rog_sq) )
   }

   if(ncore > 1) {
      RLIMIT = R_NCELL_LIMIT
      nDataSeg = floor((nrow(xyz)-1)/RLIMIT)+1
      nDataSeg = floor(nDataSeg * nseg.scale)
      lenSeg = floor(nrow(xyz)/nDataSeg)
      rog <- NULL
      for(i in 1:nDataSeg) {
         istart = (i-1)*lenSeg + 1
         iend = if(i<nDataSeg) i*lenSeg else nrow(xyz)
         rog <- c(rog, mclapply(istart:iend, function(j)
                     rg(xyz[j,], mass)))
      }
      rog <- unlist(rog)
   } else {
      rog <- apply(xyz, 1, rg, mass)
   }
   return(rog)
}
