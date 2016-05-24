#TO DO: - x is not needed when a pedigree is given, because unknown conributors must be specified in the pedigree
#       - When a pedigree is not specified, known_genotypes must not contain genotypes of non-contributors, only of the contributors!

simLR <- function(R, x, alleles, afreq, pDO, pDI, N, known_genotypes=NULL, ped=NULL, id.U=NULL, id.V=NULL) {
  # Input:
  # - R: evidence alleles
  # - x: number of unknown individuals
  # - alleles: allele names, must be vector of integers
  # - afreq: vector of allele frequencies with names(afreq)=alleles
  # - pDO: probability of drop-out, applied per allele
  # - pDI: probability of drop-in per locus
  # - N: number of simulations
  # - known_genotypes: list of genotypes of known contributors. If a pedigree is specified, the index of the contributor 
  #   must also be specified, e.g. individuals 3 and 4 with genotypes 1/2 and 2/2 are specified as list(c(3,1,2),c(4,2,2))
  # - ped: optional pedigree that specifies kinship between contributors
  # - id.U: index of untyped contributors in the pedigree. Only needed if a pedigree is specified
  # - id.V: index of typed non-contributors in the pedgree. Only needed if a pedigree is specified

  #Output:
  # - p.R: probability of the evidence
  
  fitEvid <- function(x,R){
    if(max(x) > 0) x=x[x>0]
    all(x%in%R) && all(R%in%x)
  }
  
  na <- length(alleles)
  
  #------- Genotypes for known and unknown contributors  --------
  #Without kinship
  if(is.null(ped)){
    #Genotypes of known contributors
    g <- numeric()
    if(length(known_genotypes)>0) {
      #Extract genotypes for known contributors
      gtKnown <- lapply(known_genotypes,function(x) { if(length(x)==2) x[1:2] else x[2:3] })
      g <- matrix(nrow=length(known_genotypes)*2,rep(unlist(gtKnown),N),byrow=FALSE)
    }
    #Genotypes of unknown contributors are sampled randomly from alleles
    if(x > 0) g <- rbind( g, matrix( sample(alleles, x*2*N, replace=TRUE, prob=afreq), nrow=x*2))
    
    #With kinshp
  } else {
      if (inherits(ped, c("linkdat", "singleton"))) 
        ped = list(ped)
      #Get index of individuals with known genotype
      all_typed = sapply(known_genotypes, "[", 1)
      #Index of known contributors 
      contrib_typed = setdiff(all_typed, id.V)
      #Index of unknown contributors
      contrib_untyped = id.U
      #Checks whether index of known individuals are in contrib_typed. Finds all unique alleles
      K = unique(unlist(lapply(known_genotypes, function(triple) if (triple[1] %in% contrib_typed) triple[2:3])))
      #Any alleles not explained by contributors?                                                                 
      R_not_masked = setdiff(R, K)
      if (length(alleles) == 1) 
        alleles = 1:alleles
      partialmarkers = lapply(ped, function(ped) {
        #Create empty marker
        m <- marker(ped, alleles = alleles, afreq = afreq)
        #Add known genotypes to marker
        for (tup in known_genotypes) if (tup[1] %in% ped$orig.ids) 
          m <- modifyMarker(ped, m, ids = tup[1], genotype = tup[2:3])
        m
      })
      #Genotypes of known contributors
      g <- numeric()
      if(length(contrib_typed)>0) {
        #Extract genotypes for known contributors
        gtKnown <- lapply(known_genotypes[which(all_typed%in%contrib_typed)],function(x) { if(length(x)==2) x[1:2] else x[2:3] })
        g <- matrix(nrow=length(contrib_typed)*2,rep(unlist(gtKnown),N),byrow=FALSE)
      } 
    #Genotypes for unknown contributors
    if(length(contrib_untyped)>0) {
        for( i in 1:length(ped)){ #APPLY HER???
          #Check if the unknown contributors are in the given pedigree
          available <- contrib_untyped[contrib_untyped %in% ped[[i]]$orig.ids] 
          if( length(available) > 0 ) {
            #set.seed(14)
            ysim <- markerSim(ped[[i]], N=N, available=available, partialmarker=partialmarkers[[i]],verbose=FALSE)
            #g <- rbind(g,matrix(nrow=2,unlist(lapply(ysim$markerdata, function(m) m[available,]))))
            g <- rbind(g,sapply(ysim$markerdata, function(m) t(m[available-min(ped[[i]]$orig.ids)+1,])))
          }
        }
      }
    }
    
  #------ Simulate drop-in and drop-out -------
  if(pDO==0 && pDI==0){
    p.R <- sum(apply(g,2, fitEvid, R=R))/N
    return(p.R)
  }
  #set.seed(14)
  dropInSample <- sample(0:na, 2*N, prob=c(1-pDI,afreq*pDI), replace=TRUE)
  dropInMat <- matrix(dropInSample, nrow=2, ncol=N)
  #Alleles of contributors plus drop-in alleles
  gDrop <- rbind(g,dropInMat)
  #set.seed(14)
  dropMat <- matrix(sample(c(0,1), prob=c(pDO,1-pDO), size=prod(dim(gDrop)), replace=TRUE), nrow=nrow(gDrop))
  data <- gDrop*dropMat
  #Only consider unique alleles
  sim <- apply(data,2,function(x) unique(x))
  #Count how many times the alleles fit the evidence
  if(is.matrix(sim)) p.R <- sum(apply(sim,2, fitEvid, R=R))/N 
  else p.R <- sum(sapply(sim, fitEvid, R=R))/N 
  
  return(p.R)
}
