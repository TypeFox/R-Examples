Dopt.design <- function(nruns, data=NULL, formula=~., 
    factor.names=NULL, nlevels=NULL, digits=NULL,
    constraint=NULL, center=FALSE, nRepeats=5, seed=NULL, 
    randomize=TRUE, blocks=1, block.name="Blocks", 
    wholeBlockData=NULL, ...){
    aufruf <- sys.call()

    ##wholeBlockData must be subsequently included into factor.names, nfactors, nlevels and quantitative, or not ???
    
    ## how to handle quantitative???
    
    if (!is.numeric(nruns)) stop("The number of runs (nruns) must be specified")
    if (!nruns%%1==0) stop("nruns must be an integer number")
    if (!is.numeric(blocks)) stop("blocks must be numeric")
    if (!all(blocks%%1==0)) stop("blocks must consist of integer numbers")
    if (length(blocks)==1 & !identical(blocks,1)){ 
         if (nruns%%blocks==0) blocks <- rep(nruns%/%blocks, blocks)
             else stop("the number of blocks is not a divisor of nruns")
    }
    if (!identical(blocks,1) & !sum(blocks)==nruns) stop("If specified, blocks must sum to nruns")
    
    ## transform character specification of formula
    if (is.character(formula)) formula <- try(as.formula(formula))
    if ("try-error" %in% class(formula)) stop("invalid character string for formula")
    
    ## provide nfactors
    if (!is.null(data)) nfactors <- ncol(data)
    else if (!is.null(factor.names)) nfactors <- length(factor.names)
       else if (!is.null(nlevels)) {
          nfactors <- length(nlevels)
          factor.names <- Letters[1:nfactors]
       }
    
    ## check digits
    if (!is.null(digits) & is.null(data)){ 
      if (!(length(digits)==nfactors | length(digits)==1))
         stop("if specified, digits must be a single number or must contain an entry for each factor")
         if (!is.numeric(digits)) stop("digits must be numeric")
         if (!all(floor(digits)==digits)) stop("digits must be integer")
         }
    
    if (is.null(data) & ((!is.list(factor.names)) & is.null(nlevels)))
       stop("Please provide at least data OR a valid formula and factor level information")
    if (!is.null(data) & (!is.null(factor.names) | !is.null(nlevels)))
       warning("factor.names and nlevels are ignored, since data has been specified")
    if (!(is.null(factor.names) | is.null(nlevels)) & !length(factor.names)==length(nlevels))
       stop("if both are specified, factor.names and nlevels must be of the same length")
       
    if (is.null(data) & is.list(factor.names) & !is.null(nlevels)){
       fnl <- sapply(factor.names,"length")
       if (!all(fnl==2 | fnl==nlevels))
          stop("factor.names and nlevels contradict each other")
       hilf <- which(fnl==2 & nlevels>2)
       for (i in hilf){
          if (!all(is.numeric(factor.names[[i]]))) stop("factor ", names(factor.names)[i], " is not numeric")
          factor.names[[i]] <- seq(factor.names[[i]][1],factor.names[[i]][2],
                length.out=nlevels[i])
          if (!is.null(digits)){
             if (length(digits)==length(nlevels)) factor.names[[i]] <- round(factor.names[[i]], digits=digits[i])
             else factor.names[[i]] <- round(factor.names[[i]], digits=digits)
          }
          }
    }
    ## candidate set provided as (within) data, if not already available
    if (is.null(data)){
        data <- fac.design(factor.names=factor.names, nlevels=nlevels, randomize=FALSE)
            ## character factor.names are taken care of in fac.design
        quantitative <- sapply(design.info(data)$factor.names, function(obj) all(is.numeric(obj)))
        data <- qua.design(data,quantitative)
        }
    else quantitative <- sapply(data.frame(data), is.numeric)
    
    ## apply constraint
    if (!is.null(constraint)){
        if (!is.character(constraint)) stop("constraint must be a character string")
        beding <- eval(parse(text=constraint), envir=data)
        if (!is.logical(beding))
            stop("evaluation of constraint not possible") 
        if (!length(beding)==nrow(data))
            stop("constraint does not produce a condition for each row of data")
        if (sum(beding) < nruns & identical(blocks,1)) 
            stop("The constraint reduces the candidate set to ", sum(beding), " rows.")
        data <- data[beding,]
        }
    
    ## populate factor.names, if not available
    if (is.null(factor.names)){
           factor.names <- as.list(colnames(data))
           names(factor.names) <- colnames(data)
            for (i in 1:nfactors) {
            if (quantitative[i]) factor.names[[i]] <- range(data[,i,drop=TRUE])
                else factor.names[[i]] <- names(table(data[,i,drop=TRUE]))}
         }
         
    if (!is.null(seed)) set.seed(seed)
    
    ## treat case without blocking
    if (identical(blocks,1)){
    plan <- optFederov(formula, data, nruns, augment=FALSE, center=center, nRepeats=nRepeats,
            approximate=FALSE,criterion="D",
           evaluateI=FALSE,space=NULL, ...)
        if (randomize) ord <- sample(nruns) else ord <- 1:nruns
        aus <- plan$design[ord,]
        class(aus) <- c("design", "data.frame")
        hilf <- model.matrix(formula(formula), data=aus)
        ## fixed August 2010
        if (ncol(hilf) < ncol(aus)){ 
             ## add missing columns
             if ("design" %in% class(data)) hilf2 <- desnum(data)[plan$rows,]
                 else hilf2 <- as.matrix(as.data.frame(lapply(data[plan$rows,], "as.numeric")))
             sp <- which(!colnames(hilf2) %in% colnames(hilf))
             hilf <- cbind(hilf, hilf2[,sp])
             }
       desnum(aus) <- hilf
       rownames(desnum(aus)) <- rownames(aus) <- 1:nrow(aus)
       ## changed 27.1.2011: make run.no.in.std.order a factor with proper ordering
       run.order(aus) <- data.frame(run.no.in.std.order=factor(plan$rows[ord],levels=1:nrow(data)), 
             run.no=1:nruns, run.no.std.rp=plan$rows[ord])
       }
    else{
      ## blocks with or without wholeBlockData
       if (is.null(wholeBlockData)) 
          plan <- optBlock(formula, data, blocks, center=center, criterion="D",nRepeats=nRepeats) 
       else
          plan <- optBlock(formula, data, blocks, wholeBlockData=wholeBlockData, center=center, 
            criterion="D",nRepeats=nRepeats)
        
        if (randomize) ord <- sapply(blocks,"sample",simplify=FALSE) else 
           ord <- sapply(blocks, function(obj) 1:obj, simplify=FALSE)
        for (i in 1:length(blocks)) 
           rownames(plan$Blocks[[i]]) <- 1:blocks[i] ## paste(i, 1:blocks[i], sep=".")   ## adjust row names
                                                     ## block number is handled by rbind (because of list names Bi)
        plan$Blocks <- mapply(function(matrix,select) as.data.frame(matrix[select,]), plan$Blocks, ord, SIMPLIFY=FALSE)

        ## stack all blocks beneath each other
        aus <- do.call(rbind, plan$Blocks)  ## deparse.level=1 prepends Block names to row names
        rownames(aus) <- substring(rownames(aus),2)  ## remove the "B"
             ### what would be better: call rbind with deparse.level=0; but I don't know how to
        
        ## prepend block factor
            aus <- cbind(block.name=as.factor(unlist(mapply(rep,1:length(blocks),each=blocks))),aus)
            colnames(aus)[1] <- block.name
        
        ## make ord into overall row numbers instead of block-specific ones
        ord <- unlist(ord) + unlist(mapply(rep,cumsum(c(0,blocks[-length(blocks)])),each=blocks))
        attr(aus, "desnum") <- cbind(model.matrix(~., data=aus[,1,drop=FALSE]),model.matrix(formula(formula), data=aus[,-1])[,-1])
        attr(aus, "run.order") <- data.frame(run.no.in.std.order = paste(plan$rows[ord], rownames(aus), sep="."), 
                                     run.no=1:nruns, 
                                     run.no.std.rp=paste(plan$rows[ord], rownames(aus), sep="."))
    }

    class(aus) <- c("design", "data.frame")
    design.info(aus) <- list(type="Dopt",nruns=nruns, nfactors=nfactors,factor.names=factor.names,nlevels=nlevels,
        replications=1,repeat.only=FALSE,randomize=randomize,seed=seed, creator=aufruf, quantitative=quantitative,
        digits=digits, constraint=constraint, formula=expand.formula(formula, names(factor.names)), 
        optimality.criteria=list(D=plan$D, Dea=plan$Dea, A=plan$A, G=plan$G), response.names=NULL)
    rownames(desnum(aus)) <- rownames(aus) <- 1:nrow(aus)
    if (!identical(blocks,1)){ 
        di <- design.info(aus)
        if (!is.null(wholeBlockData)){ 
            di$type <- "Dopt.splitplot"
            di$nfac.WP <- ncol(wholeBlockData)
            di$nfac.SP <- ncol(data)
            di$nWPs <- length(blocks)
            di$plotsize <- blocks
            di$plot.name <- block.name
            di$nfactors <- di$nfac.WP + di$nfac.SP
            ## prepend wp info to factor.names
            wpfn <- as.list(colnames(wholeBlockData))
            names(wpfn) <- colnames(wholeBlockData)
            di$factor.names <- c(wpfn, di$factor.names)
            
            for (i in 1:ncol(wholeBlockData)) 
                  di$factor.names[[i]] <- names(table(wholeBlockData[,i,drop=TRUE]))
            wpquan <- sapply(data.frame(wholeBlockData),is.numeric)
            names(wpquan) <- colnames(wholeBlockData)
            di$quantitative <- c(wpquan,quantitative)
        }
        else{
            di$type <- "Dopt.blocked"
            di$block.name <- block.name
            di$nblocks <- length(blocks)
            di$blocksize <- blocks
        }
      design.info(aus) <- di}
    aus
}
