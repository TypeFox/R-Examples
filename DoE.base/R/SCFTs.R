SCFTs <- function(design, digits=3, all=TRUE, resk.only=TRUE, kmin=NULL, kmax=ncol(design), 
      regcheck=FALSE, arft=TRUE, cancors=FALSE, with.blocks=FALSE){
        ## returns list with SCFTs and ARFTs
        ## either for all full resolution subsets or for all subsets
        if ("design" %in% class(design)) {
            fn <- names(factor.names(design))
            if (with.blocks) fn <- c(fn, design.info(design)$block.name)
            design <- design[,fn]
            nfac <- length(fn)
            }
        else{
            nfac <- ncol(design)
            fn <- 1:nfac
        }
        
        nlev <- levels.no(design)
        dfs <- nlev-1
        ks <- which(round(GWLP(design, kmax=kmax)[-1],8)>0) 
        if (length(ks)==0) {
                     hilf <- list(list(SCFT=cbind(SC=0, frequency=sum(nlev) - kmax), 
                                       ARFT=cbind(aveR2=0, frequency=kmax), cancor1=0))
                     names(hilf) <- kmax
                     return(hilf)   
                     }
        k <- min(ks)   ## here, k is resolution
        if (k < 2) stop("resolution of design must be at least 2")
        kminset <- FALSE   ## kminset TRUE: user has chosen kmin 
        if (is.null(kmin)) kmin <- k
        else {
            kminset <- TRUE
            redu <- ks[ks >= kmin]   ## reduced set of levels
		message(paste("Check sets of sizes ", paste(redu, collapse=",")))
            if (length(redu)==0) return()
            kmin <- min(redu) 
            k <- kmin
            ks <- redu   ## added March 2016, after commenting out the interim 0 calculations
            }
        if (k >= kmin)
        {
        ## should always be reached
        if (!"design" %in% class(design)){
           ## make sure these are factors
           if (!is.data.frame(design)) design <- as.data.frame(design)
           faktoren <- sapply(design, is.factor)
           keinfaktor <- which(!faktoren)
           if (length(keinfaktor)>0) for (i in keinfaktor) design[[i]] <- as.factor(design[[i]])
        }
        
       ### calculations for k=kmin
        k <- kmin
        nGRindi <- choose(nfac-1,k-1)
        sel <- nchoosek(nfac-1,k-1)
        auswahl <- 1:nGRindi
        erg3 <- array(NA, c(nfac, nGRindi, max(nlev)-1))
        selproj <- nchoosek(nfac,k)
          ## added March 2016
          ## determine resolution k projections
          GWLPs <- round(apply(selproj, 2, function(obj) GWLP(design[,obj])[-1]),4) 
              ## column for each projection
        selproj <- apply(selproj, 2, function(obj) paste(obj, collapse=":"))
        ergproj <- rep(NA, length(selproj)) 
        dimnames(erg3) <- list(factor=fn, others=apply(sel, 2, function(obj) paste(obj, collapse=":")), 1:(max(nlev)-1))
        for (i in 1:nfac){
          ## added March 2016
          ## determine resolution k projections
          if (resk.only){
          reskproj <- apply(GWLPs, 2, function(obj) all(obj[-k]==0))
          if (all(!reskproj)){ 
                message("no projections with resolution ", k, " or higher")
                return()
             }
          }

          seli <- sel
          seli[seli>=i] <- seli[seli>=i]+1
          proji <- apply(seli, 2, function(obj) paste(sort(c(i,obj)),collapse=":")) 
          indexes <- sapply(proji, function(obj) which(selproj==obj))
          berechn <- t(sapply(auswahl, function(obj){
            spalten <- c(i,seli[,obj])
            hilf2 <- design[,spalten]
            mmX <- model.matrix(~., data=hilf2[,1,drop=FALSE])[,-1]  ## first factor
            
            if (k>2)
            mmOther <- model.matrix(formula(substitute(~.^km1, list(km1 = k-1))),
                     data=hilf2[,-1,drop=FALSE])[,-1]
            else
            mmOther <- model.matrix(~., data=hilf2[,-1,drop=FALSE])[,-1]
            hilfc <- cancor(mmOther, mmX)$cor
            if (length(hilfc) < dfs[i]) 
               hilfc <- c(hilfc, rep(0,dfs[i]-length(hilfc)))
            hilfc
         }))
          if (nrow(berechn)==1) berechn <- t(berechn)
          erg3[i,1:nGRindi,1:dfs[i]] <- berechn
            }
        rund <- round(erg3, digits)

          SCFT <- table(round(erg3^2, digits))
          SCFT <- cbind(SC=as.numeric(names(SCFT)), frequency=SCFT)
          rownames(SCFT) <- rep("",nrow(SCFT))

        aus <- list(SCFT=SCFT)
        if (arft) {
                aveR2s <- apply(erg3, c(1, 2), function(obj) mean(obj^2, 
                   na.rm = TRUE))
                ARFT <- table(round(aveR2s, digits))
                ARFT <- cbind(aveR2 = as.numeric(names(ARFT)), frequency = ARFT)
                rownames(ARFT) <- rep("", nrow(ARFT))
                aus <- c(aus, list(ARFT=ARFT))
            }
        if (cancors) aus <- c(aus, list(cancors=round(erg3^2, digits)))

        if (regcheck) if (length(setdiff(aus$SCFT[,1], c(0,1)))>0){
              message("The design is not regular")
              return(aus)
        }
        else message("The check for k=", k, " did not exclude regularity")
        }
        else aus <- list()
        
        ## calculations for larger k, if requested
        if (all){
        aus <- list(aus)
        names(aus) <- k
        ## March 2016: added condition
        ## is this really needed, even for not only resk.only ?
        ## it will provide interim 0 tables, which are unnecessary
        ## given one knows the GWLP
        ## therefore, commented out entirely
        ## if (!resk.only) ks <- kmin:kmax  ### entered in order to also obtain SCFTs for 
                                         ### numbers of factors with GWLP 0
        for (k in ks[ks>kmin & ks<=kmax]){
          nGRindi <- choose(nfac-1,k-1)  ## number of subsets excluding LHS
          sel <- nchoosek(nfac-1,k-1)    ## subsets excluding LHS
          auswahl <- 1:nGRindi           ## index vector of subsets exluding LHS
          erg3 <- array(NA, c(nfac, nGRindi, max(nlev)-1))
          
          selproj <- nchoosek(nfac,k)    ## k factor projections
          
          ## determine resolution k projections
          GWLPs <- round(apply(selproj, 2, function(obj) GWLP(design[,obj])[-1]),4) 
              ## column for each projection
          if (resk.only){
          reskproj <- apply(GWLPs, 2, function(obj) all(obj[-k]==0))
          if (all(!reskproj)){ 
                message("no projections with resolution ", k, " or higher")
                ks <- ks[ks<k]
                break
             }
          }
          
          ## make into character descriptor
          selproj <- apply(selproj, 2, function(obj) paste(obj, collapse=":"))
          ergproj <- rep(NA, length(selproj))
          ## prepare result 
          dimnames(erg3) <- list(factor=fn, 
                                 others=apply(sel, 2, function(obj) 
                                    paste(obj, collapse=":")), 1:(max(nlev)-1))
        for (i in 1:nfac){
          seli <- sel
          seli[seli>=i] <- seli[seli>=i]+1
          proji <- apply(seli, 2, function(obj) paste(sort(c(i,obj)),collapse=":")) 
          indexes <- sapply(proji, function(obj) which(selproj==obj))
          if (resk.only) needi <- sapply(proji, function(obj) reskproj[which(selproj==obj)])
          else needi <- rep(TRUE, nGRindi)
          if (any(needi)){
          berechn <- t(sapply(auswahl[needi], function(obj){
            spalten <- c(i,seli[,obj])
            hilf2 <- design[,spalten]
            mmX <- model.matrix(~., data=hilf2[,1,drop=FALSE])[,-1]  ## first factor
            if (k>2)
            mmOther <- model.matrix(formula(substitute(~.^km1, list(km1 = k-1))),
                     data=hilf2[,-1,drop=FALSE])[,-1]
            else 
            mmOther <- model.matrix(~., data=hilf2[,-1,drop=FALSE])[,-1]
            hilfc <- cancor(mmOther, mmX)$cor
            if (length(hilfc) < dfs[i]) 
               hilfc <- c(hilfc, rep(0,dfs[i] - length(hilfc)))
            hilfc
            }))
          if (nrow(berechn)==1) berechn <- t(berechn)
          erg3[i,(1:nGRindi)[needi],1:dfs[i]] <- berechn
            }
        }
        rund <- round(erg3, digits)
        SCFT <- table(round(erg3^2, digits))
        SCFT <- cbind(SC=as.numeric(names(SCFT)), frequency=SCFT)
        rownames(SCFT) <- rep("",nrow(SCFT))
        ausadd <- list(SCFT=SCFT)

        if (arft) {
                aveR2s <- apply(erg3, c(1, 2), function(obj) mean(obj^2, 
                   na.rm = TRUE))
                ARFT <- table(round(aveR2s, digits))
                ARFT <- cbind(aveR2 = as.numeric(names(ARFT)), frequency = ARFT)
                rownames(ARFT) <- rep("", nrow(ARFT))
                ausadd <- c(ausadd, list(ARFT=ARFT))
            }

        if (cancors) ausadd <- c(ausadd, list(cancors=round(erg3^2, digits)))
        aus <- c(aus, list(ausadd))

        if (regcheck) 
        if (length(setdiff(ausadd$SCFT[,1], c(0,1)))>0){
           message("The design is not regular")
           names(aus) <- ks[1:which(ks==k)]
           return(aus)
        }
           else message("The check for k=", k, " did not exclude regularity")
        }
	  names(aus) <- ks
        }
        aus
        
        }

