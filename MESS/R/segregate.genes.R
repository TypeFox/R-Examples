segregate.genes <- function(pedigree, maf) {

  if (! "pedigree" %in% class(pedigree))
    stop("requires a single pedigree as input")

  if (any(maf>1) | any(maf<0))
    stop("minor allele frequencies must be in the range [0,1]")

  depth <- kinship2::kindepth(pedigree)
  n <- length(depth)

  res <- 
    sapply(maf, function(x) {
      marker <- rep(0, n)
      
      # Make founders
      marker[depth==0] <- sample(2:0, size=sum(depth==0), replace=TRUE, prob=c(x**2, 2*x*(1-x), (1-x)**2))

      # Count the number of MZ twin pairs
      mztwins <- pedigree$relation[pedigree$relation$code=="MZ twin",]
      nmztwins <- nrow(mztwins)
      
      # Make the remaining generations
      if (max(depth)>0) {
        for (dep in 1:max(depth)) {
          # pick the individuals for this depth
          select <- (depth==dep)
          idx <- seq(1, n)[select]
          nselect <- sum(select)
          marker[select] <- rbinom(nselect, 1, marker[pedigree$findex[select]]/2) + rbinom(nselect, 1, marker[pedigree$mindex[select]]/2)
          # Fix the MZ twins
          # This is done by running through all MZ twins and copying the
          # alleles of the first found twin to his/her MZ relatives

          if (!is.null(nmztwins)) {
            # Check if they are in this generation
            inthisdepth <- pedigree$relation$indx1[pedigree$relation$indx1 %in% idx]
            # Copy the marker result from indx1 to indx2
            marker[pedigree$relation$indx2[pedigree$relation$indx1 %in% idx]] <- marker[inthisdepth]
          }          
        }
      }
      marker
    })

  if(max(depth)==0) {
    res <- as.data.frame(matrix(res, ncol=length(maf)))
  }
  colnames(res) <- paste("gene", seq(1,length(maf)), sep="")
  res
}

