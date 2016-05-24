"fit.xyz" <-
function(fixed,
         mobile,
         fixed.inds  = NULL,
         mobile.inds = NULL,
         verbose = FALSE,
         prefix = "",
         pdbext = "",
         outpath = "fitlsq",
         full.pdbs=FALSE,
         ncore=1,
         nseg.scale=1, # to resolve the memory problem in using multicore
         ...) {

  # Parallelized by parallel package (Tue Dec 11 17:41:08 EST 2012)
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
 
  ### Addation (Mon Jul 23 17:26:16 PDT 2007)
  if( is.null(fixed.inds) && is.null(mobile.inds) ) {
    if(is.list(mobile)) {
      fixed.inds <- intersect(which(!is.gap(fixed)),
                              gap.inspect(mobile$xyz)$f.inds )
    } else {
      fixed.inds <- intersect(which(!is.gap(fixed)),
                              gap.inspect(mobile)$f.inds )
    }
    mobile.inds <- fixed.inds
    warning(paste("No fitting indices provided, using the",
                  length(fixed.inds)/3,  "non NA positions\n"))
  }
  
  
  if (is.null(fixed.inds)) fixed.inds=which(!is.gap(fixed))
  if (is.null(mobile.inds))  mobile.inds=gap.inspect(mobile)$f.inds

  if (length(fixed.inds) != length(mobile.inds))
    stop("length of 'fixed.inds' != length of 'mobile.inds'")

#  if(!is.xyz(fixed) || !is.numeric(fixed))
  if(!is.numeric(fixed))
    stop("input 'fixed' should be a numeric 'xyz' vector or matrix")

  if(is.vector(mobile)) {   # INPUT is a single vector
    if(!is.numeric(mobile))
      stop("input 'mobile' should be numeric")

    if( any(is.na(fixed[fixed.inds])) ||
       any(is.na(mobile[mobile.inds])) ) {
      stop(" NA elements selected for fitting (check indices)")
    }
    
    fit <- rot.lsq(xx=mobile,
                   yy=fixed,
                   xfit=mobile.inds,
                   yfit=fixed.inds,
                   verbose=verbose)

    return(as.xyz(fit))
  } else {
    if(is.list(mobile)) {      # INPUT is a list object
      if(!is.numeric(mobile$xyz))
        stop("non numeric input 'mobile$xyz'")
      
      if( any(is.na(fixed[fixed.inds])) ||
         any(is.na(mobile$xyz[,mobile.inds])) ) {
        stop(" NA elements selected for fitting (check indices)")
      }

      if(ncore>1 && is.matrix(mobile$xyz) ) {        # Parallelized
         RLIMIT = floor(R_NCELL_LIMIT/ncol(mobile$xyz))
         nDataSeg = floor((nrow(mobile$xyz)-1)/RLIMIT)+1
         nDataSeg = floor(nDataSeg * nseg.scale)
         lenSeg = floor(nrow(mobile$xyz)/nDataSeg)
         fit = vector("list", nDataSeg)
         for(i in 1:nDataSeg) {
            istart = (i-1)*lenSeg + 1
            iend = if(i<nDataSeg) i*lenSeg else nrow(mobile$xyz)
            fit[[i]] <- mclapply(istart:iend, function(j) rot.lsq(xx=mobile$xyz[j,],
               yy=fixed, xfit=mobile.inds, yfit=fixed.inds, verbose=verbose), 
               mc.preschedule=TRUE)
         }
         fit <- matrix(unlist(fit), ncol=ncol(mobile$xyz), byrow=TRUE)
      } else {           # Single version
         fit <- t( apply(mobile$xyz, 1, rot.lsq,
                         yy = fixed,
                         xfit = mobile.inds,
                         yfit = fixed.inds,
                         verbose=verbose))
      }
       
      if(full.pdbs) {        # FULL PDB fitting and output
        core.inds.atom = mobile.inds[seq(3,length(mobile.inds),by=3)]/3
        dir.create(outpath, FALSE)
        
        full.files  <- paste(prefix, mobile$id, pdbext, sep="")
        
        mylapply <- lapply 
        if(ncore>1) mylapply <- mclapply

#        for(i in 1:length(mobile$id)) {
        mylapply(1:length(mobile$id), function(i) {
###          pdb  <- read.pdb( paste(pdb.path,"/",mobile$id[i],pdbext,sep=""), ... )
          pdb  <- read.pdb( full.files[i], ... )
          res.resno  <- mobile$resno[i,core.inds.atom]
          res.chains <- mobile$chain[i,core.inds.atom]
          chains <- unique(res.chains[!is.na(res.chains)])

          if(length(chains)==0) {
            ##string <- paste("///",
            ##          paste(mobile$resno[i,core.inds.atom],collapse = ","),
            ##                   "///CA/", sep="")
            
            inds <- atom.select(pdb, resno=res.resno, elety="CA",
                                verbose=verbose)$xyz

          } else {
            if(length(chains)==1) {
              #string <- paste("//",chains,"/",
              #                paste(res.resno, collapse = ","),
              #                "///CA/", sep="")
              
              inds <- atom.select(pdb, resno=res.resno, chain=chains, elety="CA", 
                                  verbose=verbose)$xyz
            } else {
              # indices for each chain
              inds <- NULL
              for(j in 1:length(chains)) {
                #string <- paste("//",chains[j],"/",
                #                paste(res.resno[ res.chains==chains[j] ],
                #                      collapse = ","),
                #                "///CA/", sep="")
                
                inds <- c(inds, atom.select(pdb, resno=res.resno[ res.chains==chains[j] ],
                                            chain=chains[j], elety="CA", verbose=verbose)$xyz)
              }
            }
          }
          pdb.xyz <- pdb$xyz
          #if (het) 
          #  pdb.xyz <- c(pdb.xyz,
          #               as.numeric(t(pdb$het[,c("x","y","z")])))

          if(length(inds) > length(fixed.inds)) {
            warning("Looks like we have a multi-chain pdb with no chain id: ignoring extra indices\n\t")
            inds <- inds[1:length(fixed.inds)]
          }
          
          xyz.fit <- rot.lsq(xx=pdb.xyz,
                             yy=fixed,
                             xfit=inds, # sort!!
                             yfit=fixed.inds)

          write.pdb(xyz = xyz.fit, pdb = pdb, ##het = het, 
                    file = file.path(outpath, paste(basename(mobile$id[i]),
                      "_flsq.pdb",sep = "")) )
          return (NULL)
        } )
      }
      return(as.xyz(fit))
    } else {
      if(full.pdbs)
        warning("Need 'mobile' list object for 'full.pdbs=TRUE'")
      if(is.matrix(mobile)) {       # INPUT is a matrix
        if(!is.numeric(mobile))
          stop("input 'mobile' should be numeric")

        if( any(is.na(fixed[fixed.inds])) ||
           any(is.na(mobile[,mobile.inds])) ) {
          stop("error: NA elements selected for fitting")
        }

        if(ncore > 1) {         # Parallelized
           RLIMIT = floor(R_NCELL_LIMIT/ncol(mobile))
           nDataSeg = floor((nrow(mobile)-1)/RLIMIT)+1
           nDataSeg = floor(nDataSeg * nseg.scale)
           lenSeg = floor(nrow(mobile)/nDataSeg)
           fit = vector("list", nDataSeg)
           for(i in 1:nDataSeg) {
              istart = (i-1)*lenSeg + 1
              iend = if(i<nDataSeg) i*lenSeg else nrow(mobile)
              fit[[i]] <- mclapply(istart:iend, function(j) rot.lsq(xx=mobile[j,],
                 yy=fixed, xfit=mobile.inds, yfit=fixed.inds, verbose=verbose), 
                 mc.preschedule=TRUE)
           }
           fit <- matrix(unlist(fit), ncol=ncol(mobile), byrow=TRUE)
        } else {         # Single version
           fit <- t( apply(mobile, 1, rot.lsq,
                           yy = fixed,
                           xfit = mobile.inds,
                           yfit = fixed.inds,
                           verbose=verbose))
        }
        return(as.xyz(fit))
      }
    }
  }
}
