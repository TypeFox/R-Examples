# conversion of object class 'gpData' to data.frame

gpData2data.frame <- function(gpData,trait=1,onlyPheno=FALSE,all.pheno=FALSE,all.geno=FALSE,repl=NULL,phenoCovars=TRUE,...){

      # check for class
      if(class(gpData)!="gpData") stop("object '",substitute(gpData),"' not of class 'gpData'")
      if(is.null(dimnames(gpData$pheno))) stop("There are not dimnames for the pheno slot of ", substitute(gpData))
      IDs <- dimnames(gpData$pheno)[[1]]
      reps <- dimnames(gpData$pheno)[[3]]
      pheno <- abind(gpData$pheno, matrix(1:dim(gpData$pheno)[1]+10**ceiling(log10(dim(gpData$pheno)[1])), ncol=dim(gpData$pheno)[3], nrow=dim(gpData$pheno)[1], byrow=FALSE), along=2)
      dimnames(pheno)[[2]][dim(pheno)[2]] <- "ID"
      # choose Traits
      if(all(is.numeric(trait))) trait <- (dimnames(pheno)[[2]])[trait]
      pheno <- array(pheno[, c("ID", trait), ], dim=c(dim(pheno)[1], length(trait)+1, dim(pheno)[3]))
      dimnames(pheno) <- list(IDs, c("ID", trait), reps)
      # look for covariables
      if(phenoCovars)
        if(!is.null(gpData$phenoCovars)){
          pheno <- abind(pheno, gpData$phenoCovars, along=2)
        }
       # append column for each replication
      if(dim(pheno)[3]>1){
        pheno <- abind(matrix(rep(reps, each=dim(gpData$pheno)[1]), ncol=dim(gpData$pheno)[3], nrow=dim(gpData$pheno)[1], byrow=FALSE), pheno, along=2)
        dimnames(pheno)[[2]][1] <- "repl"
      }
      if(!is.null(repl)){
        if(all(repl %in% 1:dim(pheno)[3])) repl <- dimnames(pheno)[[3]][repl]
        else if(!all(repl %in% dimnames(pheno)[[3]])) stop("wrong replication names used")
      } else repl <- dimnames(pheno)[[3]]
      pheno <- as.data.frame(apply(pheno[, ,repl], 2, cbind))
      if(phenoCovars)
        for(i in names(gpData$info$attrPhenoCovars)){
          if(gpData$info$attrPhenoCovars[i] == "numeric")
            pheno[, i] <- as(as.character(pheno[, i]), gpData$info$attrPhenoCovars[i])
        }
      if(!is.null(repl))
        for(i in colnames(pheno))
          if(class(pheno[, i]) == "factor")
            pheno[, i] <- as.factor(as.character(pheno[, i]))
      pheno$ID <- IDs[as.numeric(as.factor(as.character(pheno$ID)))]
      if(!onlyPheno){
        geno <- gpData$geno
        # merge genotypic and phenotypic data
        mergeData <- merge(pheno,geno,by.x="ID", by.y="row.names",all.x=all.pheno,all.y=all.geno)
        # omit row.names column
       } else mergeData <- pheno
       mergeData <- mergeData[, c("ID", trait, ("repl")["repl" %in% colnames(mergeData)], colnames(mergeData)[colnames(mergeData) %in% colnames(gpData$geno)], 
                                   colnames(mergeData)[colnames(mergeData) %in% dimnames(gpData$phenoCovars)[[2]]])]
     # sort by ID
     for(i in trait)
       mergeData[, i] <- as.numeric(as.character(mergeData[, i]))
     if(all(mergeData$ID %in% 1:nrow(mergeData))) mergeData$ID <- as.numeric(mergeData$ID)
     if(is.null(mergeData$repl)){
       mergeData <- orderBy(~ID,data=mergeData)
     } else {
       if(all(mergeData$repl %in% 1:dim(gpData$pheno)[3])){
         mergeData$repl <- as.numeric(as.character(mergeData$repl))
       }
       mergeData <- orderBy(~ID+repl,data=mergeData)
     }
     mergeData$ID <- as.character(mergeData$ID)
     rownames(mergeData) <- NULL
     return(mergeData)
}
