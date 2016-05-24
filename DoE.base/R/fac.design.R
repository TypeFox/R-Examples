## full factorials of all kinds

## eventually in the wrapper package

fac.design <- function(nlevels=NULL, nfactors=NULL, factor.names = NULL, 
        replications=1, repeat.only = FALSE, randomize=TRUE, seed=NULL, 
        blocks=1, block.gen=NULL, block.name="Blocks", bbreps=replications, 
        wbreps=1, block.old.behavior = FALSE){
        ## nlevels either length 1 (if all equal) or numeric vector of length nfactors, 
        ## factor.names analogous to FrF2 (character vector or named list of levels)

        ## vector version of nlevels is sufficient
        ## list version of factor.names is sufficient
        ## scalar nlevels together with nfactors is sufficient
        
        ## if more than one of the entries are given:
        ## compatibility checks necessary
        
        ## factor levels are 1:entry of nlevels, except for 2-level factors only, where they become -1 and 1
        
        ## block generation:
        ## if only blocks is given and is a prime or a product of distinct primes, 
        ##      highest possible confounding is used, if block.gen=NULL
        ## if more than one combination of the same prime is needed, 
        ##      block.gen is required 
        ## block.gen must be a list or a matrix of numbers from 0 to p-1
        ##      or a vector (treated as one-row matrix)
        ##      appropriately matched to the columns of factorize.design

      ### check integer numbers
      creator <- sys.call()
      if (!is.null(nlevels)){ 
           if (!is.numeric(nlevels)) stop("nlevels must be numeric")
           if (!all(floor(nlevels)==nlevels)) 
              stop("nlevels must be an integer number or a vector of integer numbers.")
           if (any(nlevels < 2)) 
              stop("nlevels must not contain entries smaller than 2")
           }
      if (!is.null(nfactors)) if (!floor(nfactors)==nfactors) 
           stop("nfactors must be an integer number.")
      if (!is.null(seed)) if (!floor(seed)==seed) 
           stop("seed must be an integer number.")
      if (!floor(replications)==replications) 
           stop("replications must be an integer number.")
      if (identical(blocks,1) & !identical(wbreps,1)) 
           stop("wbreps must not differ from 1, if blocks = 1.")
      if (identical(blocks,1) & !is.null(block.gen)) 
           stop("block.gen must not be specified without specifiying the number of blocks.")
      if (bbreps > 1  & identical(blocks,1) & !replications > 1) 
        stop("Use replications, not bbreps, for specifying replications for unblocked designs.")
      ### check compatibilities of level number and factor number specifications
      ### and specify unspecified ones of these
      if (is.null(nlevels) & !is.list(factor.names)) 
             stop("If factor.names does not specify the factor levels, nlevels must be given!")
      if (is.null(nlevels) & is.list(factor.names)) if (!min(hilf <- sapply(factor.names,length))>1) 
             stop("If factor.names does not specify at least two levels for each factor, nlevels must be given!")
      if (!(is.null(nlevels) | is.null(nfactors))) if (length(nlevels)>1 & !nfactors==length(nlevels))
                          stop("nfactors does not match the length of nlevels.")
      if (is.null(nlevels)) {nlevels <- hilf
                      if (!is.null(nfactors)) if (!nfactors==length(nlevels))
                          stop("nfactors does not match the number of entries in factor.names.")}
      if (!(is.null(nlevels) | is.null(factor.names))) {
                      if (length(nlevels)>1 & !(length(factor.names)==length(nlevels)))
                          stop("length of factor.names and length of nlevels do not match.")
                      if (length(nlevels)==1) nlevels <- rep(nlevels,length(factor.names))
                      }
      if (is.null(nfactors)) nfactors <- length(nlevels)
      if (nfactors==1) stop("one factor only is not covered by fac.design")
      if (length(nlevels)==1) nlevels <- rep(nlevels, nfactors)
      if (is.list(factor.names)){ 
                             if (!(all(nlevels==sapply(factor.names,length) | sapply(factor.names,length)==1)))
                                 stop("Entries in nlevels do not match entries in factor.names.") 
            if (is.null(names(factor.names))){ if (nfactors<=50) names(factor.names) <- Letters[1:nfactors] 
                       else names(factor.names) <- paste("F",1:nfactors,sep="")
                           }}
      if (is.null(factor.names) | !is.list(factor.names)) {
                 ## null or character vector
                 hilf <- NULL
                 if (!is.null(factor.names)) hilf <- factor.names
                 factor.names <-  rep(list(numeric(0)),nfactors)
                 if (!is.null(hilf)) names(factor.names) <- hilf 
                 else if (nfactors<=50) names(factor.names) <- Letters[1:nfactors] 
                       else names(factor.names) <- paste("F",1:nfactors,sep="")
                 for (i in 1:nfactors) factor.names[i] <- list(1:nlevels[i])
             }
      if (is.list(factor.names)){ 
            if (is.null(names(factor.names))){ if (nfactors<=50) names(factor.names) <- Letters[1:nfactors] 
                       else names(factor.names) <- paste("F",1:nfactors,sep="")}
            if (any(sapply(factor.names,length)==1)) 
                 for (i in 1:nfactors) if (length(factor.names[[i]])==1) factor.names[[i]] <- 1:nlevels[i]
                 }
      ## make names valid under all circumstances
      names(factor.names) <- make.names(names(factor.names), unique=TRUE)
      
      ## check validity of blocking request
      if (!identical(blocks,1)){
        if (!is.numeric(blocks)) stop("blocks must be numeric")
        if (!round(blocks)==blocks) stop("blocks must be integer")
        if (!is.numeric(bbreps)) stop("bbreps must be an integer number.")
        if (!is.numeric(wbreps)) stop("wbreps must be an integer number.")
        ## pre-process needs for numbers of blocks and numbers of levels
        need.gen <- conf.design::factorize(blocks)
        ung <- unique(need.gen)   
        ## levels of pseudofactors,
        ## a list element for each factor
        hilfl <- conf.design::factorize(nlevels)
        names(hilfl) <- names(factor.names)
        ## numbers of pseudofactors for each factor
        lengths <- sapply(hilfl, length)
        ## levels of pseudofactors in natural order
        collevs <- unlist(hilfl)
        ## which element in collevs belongs to which original factor
        FUNC <- function(X, Y) rep(Y, X)
        pseudo.belongs <- unlist(mapply(FUNC, lengths, names(hilfl)))

        tab <- table(need.gen)   ## number of pseudo factors needed of each level

        if (is.null(block.gen)){
            ## unspecified block generators
            tab.greater3 <- tab[as.character(setdiff(unique(need.gen),c(2,3)))]
            if (length(tab.greater3)>0) 
            if (any(tab.greater3 > 1))
              stop("For this number of blocks, block.gen must be specified (see documentation)")
            ## use table of catalogued block generators for 2 and 3 levels
            ##    if up to 8 factors with 2 level contributors or 5 factors with 3 level contributors
            if (2 %in% ung){
                if (tab["2"] >= 8) 
                    stop("For this number of blocks, block.gen must be specified (see documentation)")
                else {
                    k.block2 <- tab["2"]
                    k2 <- sum(sapply(hilfl, function(obj) 2 %in% obj))
                    if (k2 > 8) stop("too many factors for 2-level related blocks with current method")
                    mult2 <- sapply(hilfl, function(obj) sum(2 == obj))
                    mult2 <- mult2[mult2>0]
                        if (k2 <= k.block2) 
                            stop("Too few factors with even number of factor levels for this number of blocks")
                }
            }
            if (3 %in% ung){
                if (tab["3"] >= 5 )
                    stop("For this number of blocks, block.gen must be specified (see documentation)")
                else {
                    k.block3 <- tab["3"]
                    k3 <- sum(sapply(hilfl, function(obj) 3 %in% obj))
                    if (k3 > 5) stop("too many factors for 3-level related blocks with current method")
                    mult3 <- sapply(hilfl, function(obj) sum(3 == obj))
                    mult3 <- mult3[mult3>0]
                        if (k3 <= k.block3) 
                           stop("Too few factors with number of factor levels a multiple of 3 for this number of blocks")
                }
            }
            ## now only one instance of each prime >3 needed
            ## and up to 8 instances for 2-level, 
            ## up to 5 instances for 3-level pseudo-factors
            
            ## too few pseudo factors with primes > 3
            for (i in setdiff(unique(need.gen),c(2,3))) 
                if (sum(sapply(hilfl, function(obj) i %in% obj)) %in% c(0,1))
                   stop("Number of blocks must not be a multiple of ", i, " for this full factorial.")
            
            ## block generators for pseudo factors with more than 3 levels
            ## NULL, if there are no such factors
            block.gen <- unlist(t(sapply(setdiff(need.gen,c(2,3)), function(obj) as.numeric(collevs==obj))))
           ### not needed, was duplicate check
           # for (i in 1:length(need.gen))
           #     if (length(unique(pseudo.belongs[which(block.gen[i,]>0)]))==1)
           #        stop("Blocks would be confounded with main effects")
            
            ## prepend block generators for three and two factors
            if (3 %in% ung){
               hilf <- Yates3[unlist(block.catlg3[which(block.catlg3$k==k3 & block.catlg3$k.block==k.block3),][paste("b",1:k.block3,sep="")])]
               ## k3 columns
               hilfmat <- do.call(rbind, hilf)
               if (k.block3==1) hilfmat <- matrix(hilfmat, nrow=1)
               hilfmat <- hilfmat[,1:k3, drop=FALSE]
               maxfacused <- colSums(hilfmat > 0)
               pseudo.belongs3 <- pseudo.belongs[which(collevs==3)]
               names(maxfacused) <- colnames(hilfmat) <- unique(pseudo.belongs3)
               ## mult3 was already defined and gives number of 3-level pseudo factors 
               ## for each factor
               
               ## length(collevs) columns
               hilf3 <- matrix(0,nrow=tab[["3"]],ncol=length(collevs))
               colnames(hilf3) <- pseudo.belongs
               for (i in 1:ncol(hilfmat)){
                   if (mult3[i]==1)
                      hilf3[,intersect(which(pseudo.belongs==(colnames(hilfmat)[i])), 
                          which(collevs==3))] <- hilfmat[,i]
                   else if (mult3[i]>1){
                   ## ith factor has more than one 3-level pseudo factor
                      
                      ## factors with multiple pseudo factors;
                      ## do not always use the same pseudo factor
                      ## no guarantee to be optimal in any way
                      hilfcols <- 
                          intersect(which(pseudo.belongs == colnames(hilfmat)[[i]]), which(collevs==3))
                      ### added with version 0.27
                      hilfrows <- which(hilfmat[,i]>0)  ## rows with non-zero entries of hilfmat
                      lr <- length(hilfrows)
                      lc <- length(hilfcols)
                      mal <- sapply(1:lr, function(obj) digitsBase(obj%%(3^lc), ndigits=lc, base=3))
                      hilf3[hilfrows, hilfcols] <- matrix(hilfmat[hilfrows,i],nrow=lr,ncol=lc)
# previous                hilf3[,hilfcols] <- matrix(hilfmat[,i],nrow=nrow(hilfmat),ncol=length(hilfcols))
                      if (!block.old.behavior){
                          ## cycle through different combinations of pseudo factors
                          ## with more than one 3-level factors
                          hilf3[hilfrows, hilfcols] <- hilf3[hilfrows,hilfcols]*t(mal)%%3
                          }
                      }
                   }
                      block.gen <- rbind(hilf3, block.gen)
                   }
            if (2 %in% ung){
               hilfmat <- matrix(0,nrow=k.block2, ncol=k2)
               hilf <- Yates[unlist(block.catlg[which(block.catlg$k==k2 & block.catlg$k.block==k.block2),][paste("b",1:k.block2,sep="")])]

               for (i in 1:nrow(hilfmat)) hilfmat[i, hilf[[i]]] <- 1
               maxfacused <- colSums(hilfmat > 0)
               pseudo.belongs2 <- pseudo.belongs[which(collevs==2)]

               names(maxfacused) <- colnames(hilfmat) <- unique(pseudo.belongs2)
               ## mult2 was already defined and gives number of 2-level pseudo factors 
               ## for each factor
               
               ## length(collevs) columns
               hilf2 <- matrix(0,nrow=tab[["2"]],ncol=length(collevs))
               colnames(hilf2) <- pseudo.belongs

               for (i in 1:ncol(hilfmat)){
                   if (mult2[i]==1)
                      hilf2[,intersect(which(pseudo.belongs==(colnames(hilfmat)[i])), which(collevs==2))] <- hilfmat[,i]
                   else if (mult2[i]>1){
                      ## ith factor has more than one 2-level pseudo factor
                       
                      ## changed with version 0.27:
                      ## factors with multiple pseudo factors;
                      ## do not always use the same pseudo factor

                          ## previously
                          ## factors with multiple pseudo factors always used the same pseudo factor
                          ## because of fears that otherwise there might be issues
                          ##     with aliasing from different cancelling out behavior
                      hilfcols <- 
                          intersect(which(pseudo.belongs == colnames(hilfmat)[[i]]), which(collevs==2))
                      hilfrows <- which(hilfmat[,i]>0)  ## rows with non-zero entries of hilfmat
                                                        ## in ith column
                      lr <- length(hilfrows)
                      lc <- length(hilfcols)
                      mal <- sapply(1:lr, function(obj) digitsBase(obj%%(2^lc), ndigits=lc))

                      hilf2[hilfrows, hilfcols] <- matrix(hilfmat[hilfrows,i], nrow=lr, ncol=lc)
                      # previous    hilf2[,hilfcols] <- matrix(hilfmat[,i], nrow=nrow(hilfmat), ncol=length(hilfcols))
                      
                      if (!block.old.behavior){
                          ## cycle through different single df contrasts for factors
                          ## with more than one 2-level factors
                          hilf2[hilfrows, hilfcols] <- hilf2[hilfrows,hilfcols]*t(mal)                          
                          }
                      }
                   }
                      block.gen <- rbind(hilf2, block.gen)
                   }
              colnames(block.gen) <- pseudo.belongs 
             }
        else{
            ## specified block generators
            if (! is.numeric(block.gen)) 
                stop("If given, block.gen must be a numeric matrix, or a numeric vector.")
            if (!is.matrix(block.gen)){
                nn <- names(block.gen) 
                block.gen <- matrix(block.gen, nrow=1)
                if (!is.null(nn)) colnames(block.gen) <- nn
                else colnames(block.gen) <- pseudo.belongs
                }
            if (!nrow(block.gen)==length(need.gen))
               stop(nrow(block.gen)," block generators specified, ", length(need.gen), " would be needed")
            if (!ncol(block.gen)==length(collevs))
               stop("coefficients for ", ncol(block.gen), " pseudo-factors were specified, ",
                  length(collevs), " would be needed")
               }
            
            ## continue of checking 

            ## identify relevant prime groups
            pg <- vector(mode="list", length=length(ung))
            for (i in 1:nrow(block.gen)){
                hilf <- block.gen[i,,drop=TRUE]
                chilf <- which(hilf>0)
                if (length(lev <- unique(collevs[chilf])) > 1) 
                    stop("each block generator must address pseudo factors with the same number of levels only")
                if (!lev %in% ung)
                    stop("the ", i, "th generator is not compatible with the requested number of blocks")
                pg[[which(ung==lev)]] <- rbind(pg[[which(ung==lev)]],hilf)
            }
            ## pg should now be list of separated generator matrices per prime
            ## check for correct number of elements 
            ## for correct element types
            ## and for too severe confounding
            for (i in 1:length(pg)){
                 if (!table(need.gen)[i] == nrow(pg[[i]]))
                     stop("something went wrong, wrong number of generators for prime ", ung[i])
                 if (any(!pg[[i]] %in% 0:(ung[i]-1)))
                     stop("wrong entries in block.gen for prime ", ung[i])
                 ### check confounding using conf.set function from conf.design
                 hilf <- conf.set(pg[[i]], ung[i])
                 if (!is.matrix(hilf)) hilf <- matrix(hilf, nrow=1)
                 if (!nrow(hilf)==(ung[i]^tab[as.character(ung[i])]-1)/(ung[i]-1))
                      stop("dependant block generators for prime ", ung[i])
                 for (j in 1:nrow(hilf)) {
                      hilf2 <- length(unique(pseudo.belongs[which(hilf[j,]>0)]))
                      if (hilf2==1) stop("confounding of blocks with main effects")
                      if (hilf2==2) warning("confounding of blocks with 2-factor interactions")
                 }
            }
        }
      
      nruns <- prod(sapply(factor.names,"length"))
     ### changed cat to message with version 0.27
      message("creating full factorial with ", nruns, " runs ...\n")

      design <- try(expand.grid(factor.names))
      if ("try-error" %in% class(design)) 
          stop("design with ", nruns, " runs could not be generated")
      row.names(design) <- 1:nruns 
      ## process blocking request (check feasibility for block.gen)
        if (!is.null(block.gen)){ 
              pseudo <- factorize.data.frame(design, long=TRUE)
              for (i in 1:length(pg)){
                  ## calculated block contributor for each relevant prime
                  pg[[i]] <- (pseudo%*%t(pg[[i]]))%%ung[i]
              }
      ## not elegant but works
      ## join does not necessarily work because pg elements can be matrices
       blockcol <- rep("",nruns)
               for (i in 1:length(pg))
                  for (j in 1:ncol(pg[[i]]))
                      blockcol <- paste(blockcol, pg[[i]][,j],sep="")
               blockcol <- as.factor(as.numeric(as.factor(blockcol)))
              ## attach this column to the design only later,
              ## because this prevents case distinctions
          }
      ## end of blocking
                desnum <- NULL
                quant <- sapply(factor.names, "is.numeric")
                for (i in 1:nfactors){
                    if (!is.factor(design[,i]))
                       design[,i] <- factor(design[,i],levels=factor.names[[i]]) 
                    if (nlevels[i]==2) contrasts(design[,i]) <- contr.FrF2(2)
                    else if (quant[i]) contrasts(design[,i]) <- contr.poly(nlevels[i],scores=factor.names[[i]])
                }
      
      ## prepend block column, if needed
      if (!identical(blocks, 1)){ 
         design <- cbind(blockcol, design)
         colnames(design) <- c(block.name, names(factor.names))
         design <- design[ord(data.frame(blockcol)),]
         }

      ## simple randomization situations
      if (identical(blocks, 1)){ 
          rand.ord <- rep(1:nrow(design),replications)
          if (replications > 1 & repeat.only) rand.ord <- rep(1:nrow(design),each=replications)
          if (randomize & !is.null(seed)) set.seed(seed)
          if (randomize & !repeat.only) for (i in 1:replications) 
                      rand.ord[((i-1)*nrow(design)+1):(i*nrow(design))] <- sample(nrow(design))
          if (randomize & repeat.only) rand.ord <- rep(sample(1:nrow(design)), each=replications)
          aus <- design[rand.ord,]
      }
      else{
         ## blocked randomization
          #### implement randomization and replication for blocks
          ## bbreps=replications, wbreps=1
          ## if ((!repeat.only) & !randomize)
              nblocks <- blocks
              blocksize <- nrow(design)%/%nblocks
              
              rand.ord <- rep(1:nruns, bbreps * wbreps)
              if ((!repeat.only) & !randomize)
                 for (i in 0:(nblocks-1))
                    for (j in 1:wbreps)
                    rand.ord[(i*blocksize*wbreps+(j-1)*blocksize+1):((i+1)*blocksize*wbreps+j*blocksize)] <- 
                         (i*blocksize+1):((i+1)*blocksize)
                    rand.ord <- rep(rand.ord[1:(nruns*wbreps)],bbreps)
              if (repeat.only & !randomize)
                    for (j in 1:wbreps)
                    rand.ord[(i*blocksize*wbreps + (j-1)*blocksize + 1) : 
                          (i*blocksize*wbreps + j*blocksize)] <- 
                               sample((i%%nblocks*blocksize+1):(i%%nblocks*blocksize+blocksize))
              if (wbreps > 1 & repeat.only) rand.ord <- rep(1:nruns,bbreps, each=wbreps)
              if ((!repeat.only) & randomize)
                 for (i in 0:(nblocks*bbreps-1))
                    for (j in 1:wbreps)
                    rand.ord[(i*blocksize*wbreps + (j-1)*blocksize + 1) : 
                          (i*blocksize*wbreps + j*blocksize)] <- 
                               sample((i%%nblocks*blocksize+1):(i%%nblocks*blocksize+blocksize))
                               
              if (repeat.only & randomize)
                 for (i in 0:(nblocks*bbreps-1))
                    rand.ord[(i*blocksize*wbreps + 1) : 
                          ((i+1)*blocksize*wbreps)] <- rep(sample((((i%%nblocks)*blocksize)+1):
                                        ((i%%nblocks+1)*blocksize)),each=wbreps)
              aus <- design[rand.ord,]
      }
      ## extract run number in standard order
      ## remove uniqueness appendix
      orig.no <- orig.no.rp <- sapply(strsplit(rownames(aus),".",fixed=TRUE),function(obj) obj[1])
      ## row added 27 01 2011 (for proper ordering of design)
      orig.no.levord <- sort(as.numeric(orig.no),index=TRUE)$ix
      if (blocks > 1) {
            orig.no.levord <- sort(100000*as.numeric(aus[[block.name]])+as.numeric(orig.no),index=TRUE)$ix
            orig.no <- paste(orig.no,aus[[block.name]], rep(1:blocksize,nblocks)[rand.ord], sep=".")
            }
      orig.no.rp <- orig.no
      if (bbreps * wbreps > 1){
           if (bbreps > 1) {
                ## !repeat.only covers all blocked cases and the repeat.only standard cases
                ## since bbreps stands in for replications
                if (repeat.only & identical(blocks,1))
                orig.no.rp <- paste(orig.no.rp, rep(1:bbreps,nruns),sep=".")
                else
                orig.no.rp <- paste(orig.no.rp, rep(1:bbreps,each=nruns*wbreps),sep=".")
           }
           if (wbreps > 1){
                ## blocked with within block replications
                if (repeat.only) 
                     orig.no.rp <- paste(orig.no.rp, rep(1:wbreps,nruns*bbreps),sep=".")
                else orig.no.rp <- paste(orig.no.rp, rep(1:wbreps,each=blocksize,nblocks*bbreps),sep=".")
           }
    }

      if (is.null(desnum)) desnum <- model.matrix(1:nrow(aus)~.,data=aus)[,-1,drop=FALSE] else
           desnum <- desnum[rand.ord,]
      rownames(aus) <- rownames(desnum) <- 1:nrow(aus)

      di <- list(type="full factorial", 
          nruns=nruns, nfactors=nfactors, nlevels=nlevels, factor.names=factor.names,
          replications=replications, repeat.only=repeat.only, 
          randomize=randomize, seed=seed, creator=creator)
      if (blocks>1){ 
            di$type <- "full factorial.blocked"
            di$block.name=block.name
            di$bbreps <- bbreps
            di$wbreps <- wbreps
            if (di$wbreps==1) di$repeat.only <- FALSE
            di$nblocks <- nblocks
            di$blocksize <- blocksize
            di$block.gen=block.gen
      }

      attr(aus,"desnum") <- desnum
      ## change 27 Jan 2011: leave orig.no as a factor, but with better-ordered levels
      orig.no <- factor(orig.no, levels=unique(orig.no[orig.no.levord]))
      attr(aus,"run.order") <- data.frame("run.no.in.std.order"=orig.no,"run.no"=1:nrow(aus),"run.no.std.rp"=orig.no.rp)
      attr(aus,"design.info") <- di
      class(aus) <- c("design", "data.frame")
      aus
}
