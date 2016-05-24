length5 <- function(design, with.blocks=FALSE, J=FALSE, rela=FALSE){
    ## function to calculate generalized words of length 4
    ## according to Xu and Wu 2001 Annals
    ## it might be helpful to locate non-zeros
    ## this is so far not done
    if (rela & J) stop("rela and J must not be simultaneously TRUE")
    
    n <- nrow(design)
    if (!(is.data.frame(design) | is.matrix(design))) stop("design must be a data frame or a matrix")
    if (is.matrix(design)) design <- as.data.frame(design)
    
    if (rela) if (!(isTRUE(all.equal(length2(design), 0)) & isTRUE(all.equal(length3(design), 0)) & isTRUE(all.equal(length4(design),0)))) 
       stop("relative length5 is applicable for resolution V or higher designs only")
    
    if (!"design" %in% class(design)){
        for (i in 1:ncol(design)){
            design[,i] <- factor(design[,i])
            contrasts(design[,i]) <- "contr.XuWu"
        }
        nlevels <- sapply(as.list(design), function(obj) nlevels(obj))
        ##names(nlevels) <- colnames(design)
        fo <- formula("~.", data=design)
        }
    else{
        di <- design.info(design)
        nlevels <- di$nlevels
        if (is.null(nlevels)){
            if (length(grep("FrF2",di$type))>0 | length(grep("pb",di$type))>0 )
                nlevels <- rep(2,length(di$factor.name))
        }
        ## orthogonal contrasts
        design <- change.contr(design, "contr.XuWu")

        ## if blocked and requested, accomodate blocks
        if (with.blocks & !is.null(di$block.name)){
          if (!is.factor(design[[di$block.name]])) design[[di$block.name]] <- factor(design[[di$block.name]])
          contrasts(design[[di$block.name]]) <- "contr.XuWu"
          fo <- formula(paste("~",paste(c(di$block.name,names(di$factor.names)),collapse="+")), data=design)
          nlevels <- c(di$nblocks,nlevels)
        }
        else
          fo <- formula(paste("~",paste(names(di$factor.names),collapse="+")), data=design)
    }
    ## create model matrix

    mm <- model.matrix(fo,design)
    ## store column allocation to factors
    zuord <- attr(mm, "assign")
    ## no longer necessary, already done by contr.XuWu
    ## --> works for imbalanced designs as well
    ## normalize to sqrt(n)
    ## mm <- sqrt(n)*mm %*% diag(1/sqrt(colSums(mm^2)))
    ## remove intercept column
    mm <- mm[,-1]
    zuord <- zuord[-1]
    ##colnames(mm) <- zuord
    nfac <- max(zuord)
    if (nfac < 5) return(0)

    ## 5fi columns
    quintuples <- nchoosek(nfac, 5)
    anz <- 0
    for (i in 1:ncol(quintuples))
       anz <- anz+prod(nlevels[quintuples[,i]]-1)
    vec5 <- rep(NA, anz)
    zaehl <- 1

    ## same order as quintuples
    for (i in 1:(nfac-4)){
     icols <- which(zuord==i)
      for (j in (i+1):(nfac-3)){
        jcols <- which(zuord==j)
      for (k in (j+1):(nfac-2)){
        kcols <- which(zuord==k)
      for (l in (k+1):(nfac-1)){
        lcols <- which(zuord==l)
      for (m in (l+1):nfac){
        mcols <- which(zuord==m)
          for (a in icols){
           for (b in jcols){
             for (c in kcols){
             for (d in lcols){
             for (e in mcols){
              vec5[zaehl] <- sum(mm[,a]*mm[,b]*mm[,c]*mm[,d]*mm[,e])
              if (rela) vec5[zaehl] <- vec5[zaehl]/sqrt(min(nlevels[c(i,j,k,l,m)])-1)
              ## namvec[zaehl] <- paste(i,j,k,l,sep=":") ## update names
              zaehl <- zaehl+1
           }}}}}
       }
    }}}}
    if (J) {
        names(vec5) <- unlist(apply(quintuples,2,function(obj) rep(paste((1:nfac)[obj],collapse=":"), prod(nlevels[obj]-1))) )
        return(abs(vec5)/n)
        }
    else
    sum(vec5^2)/(n^2)
}