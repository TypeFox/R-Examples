length4 <- function(design, with.blocks=FALSE, separate=FALSE, J=FALSE, rela=FALSE){
    ## function to calculate generalized words of length 4
    ## according to Xu and Wu 2001 Annals
    if (separate & J) stop("separate and J must not both be TRUE")
    if (rela & J) stop("rela and J must not be simultaneously TRUE")
    ## this is a precaution; presumably, separate and rela would work together well
    if (separate & rela) stop("separate and rela must not both be TRUE")

    ## it might be helpful to locate non-zeros
    ## this is so far not done
    n <- nrow(design)
    if (!(is.data.frame(design) | is.matrix(design))) stop("design must be a data frame or a matrix")
    if (is.matrix(design)) design <- as.data.frame(design)
    
    if (rela) if (!(isTRUE(all.equal(length2(design), 0)) & isTRUE(all.equal(length3(design), 0)))) 
       stop("relative length4 is applicable for resolution IV or higher designs only")
    
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
    if (nfac < 4) return(0)

    ## 4fi columns
    quadruples <- nchoosek(nfac, 4)
    anz <- 0
    if (separate) {
       ## columns of mm4 to be assigned to quadruple 4fis
       anz.sep <- apply(quadruples, 2, function(obj) prod(nlevels[obj]-1))
       anz <- sum(anz.sep)  ## total column number of matrix mm4
       col.from <- c(1, (cumsum(anz.sep)+1)[-length(anz.sep)])
       col.to <- cumsum(anz.sep)
       rm(anz.sep)
       }
    else for (i in 1:ncol(quadruples))
       anz <- anz+prod(nlevels[quadruples[,i]]-1)
    vec4 <- rep(NA, anz)
  #  if (!separate) rm(quadruples)
    ##namvec <- rep("sp",ncol(mm4))  ## column names, only not attached as such
    ## takes up time, not needed because of quadruples

    zaehl <- 1

    ## same order as quadruples
    for (i in 1:(nfac-3)){
     icols <- which(zuord==i)
      for (j in (i+1):(nfac-2)){
        jcols <- which(zuord==j)
      for (k in (j+1):(nfac-1)){
        kcols <- which(zuord==k)
      for (l in (k+1):nfac){
        lcols <- which(zuord==l)
          for (a in icols){
           for (b in jcols){
             for (c in kcols){
             for (d in lcols){
              vec4[zaehl] <- sum(mm[,a]*mm[,b]*mm[,c]*mm[,d])
              if (rela) vec4[zaehl] <- vec4[zaehl]/sqrt(min(nlevels[c(i,j,k,l)])-1)
              ## namvec[zaehl] <- paste(i,j,k,l,sep=":") ## update names
              zaehl <- zaehl+1
           }}}}
       }
    }}}
    if (separate){
        fi2s <- nchoosek(nfac, 2)
        aus <- rep(NA, ncol(fi2s))
        for (i in 1:ncol(fi2s)){
              sel <- apply(quadruples, 2, function(obj) all(fi2s[,i] %in% obj))
                   ## logical selecting columns from quadruples = mm4
              ## try and expand to complete matrix ?
              cols <- unlist(mapply(":", col.from[sel], col.to[sel]))
              aus[i] <- sum(vec4[cols]^2)/(n^2)
              }
        names(aus) <- apply(fi2s, 2, "paste", collapse=":")
        return(aus)
    }
    else{
    if (J) {
        names(vec4) <- unlist(apply(quadruples,2,function(obj) rep(paste((1:nfac)[obj],collapse=":"), prod(nlevels[obj]-1))) )
        return(abs(vec4)/n)
        }
    else
    sum(vec4^2)/(n^2)}
}