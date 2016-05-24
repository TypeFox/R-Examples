GRind <- function(design, digits=3, arft=TRUE, scft=TRUE, cancors=FALSE, with.blocks=FALSE){
        ## returns list of four or five 
        ## with first element GRs a vector of GR and GRind, 
        ## second element GR.i a matrix 2 x nfac elements
        	## with rows for GRtot.i and GRind.i
        ## third element the average R2 frequency table
        ## fourth element the squared canonical correlations table
        ## fifth element the canonical correlations
        	## fifth element is an nfac x choose(nfac-1, R-1) x max(nlev)-1 array 
            ## returned only if cancors is set to TRUE
            ## the canonical correlations are supplemented with 0es, 
            ## if there are fewer than the respective dfi
        ## version 0.27: make sure only factors are included in calculations
        ## with.blocks only used for class design
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
        ks <- which(round(GWLP(design)[-1],8)>0)
        if (length(ks)==0) return(list(GRtot=Inf, GRind=Inf, GRind.i=rep(Inf,nfac), cancor1=0))
        k <- min(ks) 
        ## 10/08/15: resolution requirement lowered to 2, avoided warning in case of Inf
        if (k<2) stop("resolution of design must be at least 2")
        
        ## now there is a finite resolution of at least 2   
        if (!"design" %in% class(design)){
           ## make sure these are factors
           if (!is.data.frame(design)) design <- as.data.frame(design)
           faktoren <- sapply(design, is.factor)
           keinfaktor <- which(!faktoren)
           if (length(keinfaktor)>0) for (i in keinfaktor) design[[i]] <- as.factor(design[[i]])
        }
        nGRindi <- choose(nfac-1,k-1)
        sel <- nchoosek(nfac-1,k-1)
        auswahl <- 1:nGRindi
        erg3 <- array(NA, c(nfac, nGRindi, max(nlev)-1))
        selproj <- nchoosek(nfac,k)
        selproj <- apply(selproj, 2, function(obj) paste(obj, collapse=":"))
        ergproj <- rep(NA, length(selproj)) 
        dimnames(erg3) <- list(factor=fn, others=apply(sel, 2, function(obj) paste(obj, collapse=":")), 1:(max(nlev)-1))
        for (i in 1:nfac){
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
            if (length(hilfc)<dfs[i]) 
               hilfc <- c(hilfc, rep(0,dfs[i]-length(hilfc)))
            hilfc
            }))
          if (nrow(berechn)==1) berechn <- t(berechn)
          erg3[i,1:nGRindi,1:dfs[i]] <- berechn
            }
        rund <- round(erg3, digits)
        aveR2s <- apply(erg3,c(1,2),function(obj) mean(obj^2, na.rm=TRUE))
        if (arft){
          ARFT <- table(round(aveR2s,digits))
          ARFT <- cbind(aveR2=as.numeric(names(ARFT)), frequency=ARFT)
          rownames(ARFT) <- rep("",nrow(ARFT))
          }
        if (scft){
          SCFT <- table(round(erg3^2, digits))
          SCFT <- cbind(SC=as.numeric(names(SCFT)), frequency=SCFT)
          rownames(SCFT) <- rep("",nrow(SCFT))
          }
        aus <- list(GRs=c(GR=k+1-round(sqrt(max(apply(erg3,c(1,2),
                                   function(obj) mean(obj^2, na.rm=TRUE)))),digits),
                          GRind=k+1-max(rund, na.rm=TRUE)), 
                    GR.i=rbind(GRtot.i=k+1-round(sqrt(apply(aveR2s,1, max)), digits),
                               GRind.i=k+1-apply(rund, 1, max, na.rm=TRUE)))
        if (arft) aus <- c(aus, list(ARFT=ARFT))
        if (scft) aus <- c(aus, list(SCFT=SCFT))
        if (cancors) aus <- c(aus, list(cancors=round(erg3^2, digits)))
        class(aus) <- c("GRind", class(aus))
        aus
        }

