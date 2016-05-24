# adding markers to gpData object

add.markers <- function(gpData,geno,map=NULL){
  class(gpData$map) <- "data.frame"
  # check if markers are allready in data
  if (any(colnames(geno) %in% colnames(gpData$geno)) | any(rownames(map) %in% colnames(gpData$geno))){
    stop("some of the markers of ", substitute(geno)," are allready in ", substitute(gpData))
  }
  if(is.null(gpData$map) & !is.null(map)) stop("There is no map available for ", substitute(gpData), "!")
  if(!all(rownames(geno) %in% gpData$covar$id)) stop("You like to put new individuals into the data set!\nUse add.individuals() for that!")
  if(!all(rownames(map) %in% colnames(geno))) stop("There are markers in the map, which don't have information in ", substitute(geno), "!")
  # take names form map if available
  if(is.null(colnames(geno)) & !is.null(map))
    if(ncol(geno) == nrow(map)) colnames(geno) <- rownames(map) else
      stop("Check the colnames of", substitute(geno), " and the rownames of ", substitute(map), "!")

  # merge genotypic data
  nmiss <- sum(!rownames(gpData$geno) %in% rownames(geno))
  if(nmiss > 0){
    geno <- rbind(matrix(NA, nrow = nmiss, ncol = ncol(geno)), geno)
    rownames(geno)[1:nmiss] <- rownames(gpData$geno)[!rownames(gpData$geno) %in% rownames(geno)]
  }
  geno <- geno[match(rownames(gpData$geno), rownames(geno)), ]
  gpData$geno <- cbind(gpData$geno,geno)
  # first column as rownames and delete first column
  # merge map
  if(is.null(map)){
    #map <- gpData$map[1:ncol(geno),]
    map <- data.frame(chr=rep(NA,ncol(geno)),pos=rep(NA,ncol(geno)))
    rownames(map) <- colnames(geno)
  } else if(nrow(map) != ncol(geno)){
    map[colnames(geno)[!colnames(geno) %in% rownames(map)],] <- NA
  }
  map <- rbind(gpData$map,map)
  # same approach as in create.gpData
  map$sor <- substr(map$chr, nchar(as.character(map$chr)), nchar(as.character(map$chr)))
  if(!all(!unique(map$sor)[!is.na(unique(map$sor))] %in% 0:9)) map$sor <- 1
  map <- map[order(as.character(rownames(map))),]
  map <- orderBy(~sor+chr+pos,data=map)
  map$sor <- NULL
  gpData$map <- map
  gpData$geno <- gpData$geno[, match(rownames(gpData$map), colnames(gpData$geno))]
  gpData$info$codeGeno <- FALSE
  # sortcolumns in geno, too
  gpData$geno <- gpData$geno[,rownames(gpData$map)]
  # create new gpData object
  class(gpData$map) <- "GenMap"
  return(gpData)
}


# adding new individuals to gpData object

add.individuals <- function(gpData,pheno=NULL,geno=NULL,pedigree=NULL,covar=NULL,repl=NULL){
      # check if individuals are allready in data
      if (any(rownames(pheno) %in% gpData$covar$id) | any(rownames(geno) %in% gpData$covar$id) |
          any(pheno$ID %in% gpData$covar$id)){
        stop("some of the individuals of are already in ", substitute(gpData))
      }
      if(is.null(repl)) repl <- "repl"
      else if(!all(unique(pheno[, repl]) %in% dimnames(gpData$pheno)[[3]])) stop("Your values for replication is not in the dimnames of ", substitute(gpData$pheno))
      colnames(pheno)[colnames(pheno) == repl] <- "repl"
      # merge phenotypic data
      if(dim(gpData$pheno)[3] == 1) {
        if(!"ID" %in% colnames(pheno)) pheno$ID <- rownames(pheno)
        repl <- NULL
      } else {
        if(!"ID" %in% colnames(pheno))  stop("In", substitute(pheno), "the columns 'ID' and/or 'repl' are/is missing!")
      }
      if(!is.null(pheno))
        if(!all(!colnames(pheno)[colnames(pheno)!="ID"]%in%dimnames(gpData$pheno)[[2]] |
                !colnames(pheno)[colnames(pheno)!="ID"]%in%dimnames(gpData$phenoCovars)[[2]]))
          stop("different phenotypes (colnames) in '", substitute(gpData$pheno), "' or '",  substitute(gpData$phenoCovar),"' and '", substitute(pheno), "'")
      df.pheno <- gpData2data.frame(gpData, onlyPheno=TRUE, trait = dimnames(gpData$pheno)[[2]])
      if(!all(colnames(df.pheno) %in% colnames(pheno))) warning("Not all traits and covariates are available in the new data!")
      pheno[, colnames(df.pheno)[!colnames(df.pheno) %in% colnames(pheno)]] <- NA
      pheno <- pheno[, colnames(df.pheno)]
      pheno <- rbind(df.pheno, pheno)
      rm(df.pheno)
      if(dim(gpData$pheno)[3] == 1) {
        rownames(pheno) <- pheno$ID
        pheno$ID <- NULL
      }
      # merge genotypic data
      if(!is.null(geno)) if(any(colnames(geno)!=colnames(gpData$geno))) stop("different markers (colnames) in 'gpData$geno' and 'geno'")
      geno <- rbind(gpData$geno,geno)

      # merge pedigree
      pedigree <- rbind(gpData$pedigree,pedigree)
      # reorder if necessary
      if(!is.null(pedigree)) pedigree <- orderBy(~gener,pedigree)

      # merge covar  (not with first columns)

      if(any(colnames(covar) %in% c("genotyped","phenotyped","id"))) stop("not specify columns 'genotyped','phenotyped' and 'id' in 'covar' ")
      if(!is.null(covar)){
        cc <- data.frame(gpData$covar[,!colnames(gpData$covar) %in% c("genotyped","phenotyped","id")])
        rownames(cc) <- rownames(gpData$covar)
        colnames(cc) <- colnames(gpData$covar)[!colnames(gpData$covar) %in% c("genotyped","phenotyped","id")]
        covarUpdate <- rbind(cc,covar)
      }
      else covarUpdate <- gpData$covar

      # need id in rownames for create.gpData
      rownames(covarUpdate) <- c(as.character(gpData$covar$id),rownames(covar))

      # create new gpData object
      ret <- create.gpData(pheno=pheno,geno=geno,map=gpData$map,pedigree=pedigree,covar=covarUpdate,map.unit=gpData$info$map.unit,modCovar=dimnames(gpData$phenoCovars)[[2]],repeated=repl)
      return(ret)
}
