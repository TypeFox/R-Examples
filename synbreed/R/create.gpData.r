# read genomic prediction data
create.gpData <- function(pheno=NULL,geno=NULL,map=NULL,pedigree=NULL,family=NULL,covar=NULL,
                          reorderMap=TRUE,map.unit="cM",repeated=NULL,modCovar=NULL){

  # start with some checks on data
  # geno as matrix but not data.frame (storage)
  if(!map.unit %in% c("cM", "bp", "kb", "Mb")) warning("The measurement unit for the positions in the map should be either 'cM', 'bp', 'kb' or 'Mb'")
  if(!is.null(geno)){
    if(is.data.frame(geno)){
      geno <- as.matrix(geno)
      #if(any(duplicated(geno,MARGIN=1))) warning("individuals with duplicated genotypes")
    }
    if(!is.matrix(geno)) stop("geno must be a matrix or data.frame, not a ", class(geno))
    if(anyDuplicated(rownames(geno)))
      stop(paste("In", substitute(geno), " are duplicated names of genotypes in rownames!"))
    if(is.null(rownames(geno)) && is.null(pheno)) stop('rownames(geno) cannot be NULL unless a pheno of the same length is supplied')
    if(is.null(colnames(geno)) && is.null(map)) warning('colnames(geno) should not be NULL unless a map of the same length is supplied')
  }

  if(!is.null(map)){
    if(!is.data.frame(map)) stop('map must be a data.frame, not a ', class(map))
    # as a data.frame, map already has rownames
    if(!all(c("chr","pos") %in% colnames(map))) stop("colnames(map) must include 'chr' and 'pos'")
    # test if positions in map are numeric
    if(!is.numeric(map$pos)) stop("Position informations have to be numeric values!")

    # test if chr in map is numeric or character
    if(!(is.numeric(map$chr)|is.character(map$chr)))
      warning(paste("Chromosome information (in map$chr) should be numeric or character, not",class(map$chr)))
  } else map.unit <- NULL

  # match geno and map
  if(!is.null(geno) & !is.null(map)){
    if(!all(colnames(geno) %in% rownames(map))) {
      warning("not all markers in 'geno' mapped in 'map'. gaps filled with 'NA' \n")
    }
    if(ncol(geno) == nrow(map) && is.null(colnames(geno))){
      # assuming same markers in geno and map
      warning("missing colnames in 'geno': assuming to be identical as rownames in 'map' because of identical length \n")
      colnames(geno) <- rownames(map)
    }else{
      if(!identical(colnames(geno), rownames(map)))
        # markers must be reordered
        map <- data.frame(chr=map$chr[match(colnames(geno),rownames(map))],
                          pos=map$pos[match(colnames(geno),rownames(map))],
                          row.names=colnames(geno))
    }
  }

  phenoCovars <- NULL
  attrModCovars <- NULL
  if(!is.null(pheno)){
    if(2!=length(dim(pheno))) stop('pheno must be 2-dimensional, not ', length(dim(pheno)), "-dimensional")
    classList <- unlist(lapply(pheno, class))
    if(!all((classList[!names(classList) %in% repeated & !names(classList) %in% modCovar])[-1] %in% c("numeric", "integer"))) stop("Trait values have to be numeric!")
    # repeated measures? Use rownames of pheno as identifier for genotypes
    if(is.null(repeated)){
      if(is.null(rownames(pheno))) stop('rownames(pheno) cannot be null when not using repeated measures!')
      if(anyDuplicated(rownames(pheno))) warning("rownames in pheno should be unique when not using repeated measures!")
      add <- 10^ceiling(log10(nrow(pheno)))
      # if the rownames are numbers, then temporarily add a big number to allow alphabetical sorting.
      if(all(rownames(pheno) %in% 1:nrow(pheno))) rownames(pheno) <- add + as.numeric(rownames(pheno)) else add <- NULL
      if(dim(pheno)[2] ==1){# only a vector of traits
        phenoNames <- dimnames(pheno)
        arrPheno <- array(pheno[order(phenoNames[[1]]), ], dim = c(length(phenoNames[[1]]), 1, 1))
        dimnames(arrPheno) <- list(phenoNames[[1]][order(phenoNames[[1]])], phenoNames[[2]], "1")
      } else {# more than one trait, still unreplicated
        pheno <- pheno[order(rownames(pheno)), ]
        pheno.not.modCovar <- pheno[, !(colnames(pheno) %in% modCovar)]
        arrPheno <- array(as.matrix(pheno.not.modCovar), dim=c(dim(pheno.not.modCovar), 1))
        dimnames(arrPheno) <- list(rownames(pheno), colnames(pheno.not.modCovar), "1")
      }
      if(!is.null(add)) dimnames(arrPheno)[[1]] <- as.numeric(dimnames(arrPheno)[[1]]) - add
      if(!is.null(modCovar)){
        arrModCovars <- array(1,  dim=c(dim(pheno[, colnames(pheno) %in% modCovar]), 1))
        dimnames(arrModCovars) <- list(dimnames(arrPheno)[[1]], colnames(pheno)[colnames(pheno) %in% modCovar], "1")
        for(i in colnames(pheno)[colnames(pheno) %in% modCovar])
          arrModCovars[, i, 1] <- pheno[, i]
      }
    } else {# a vector with replication identifier is applied. The first column is the identifier for genotypes
      dim3 <- data.frame(unique(pheno[, repeated]))
      colnames(dim3) <- repeated
      dim3 <- orderBy(as.formula(paste("~", paste(repeated, collapse = " + "))), data = dim3)
      for(i in 1:ncol(dim3)) dim3[, i] <- as.character(dim3[, i])
      rownam <- sort(unique(pheno[, 1]))
      if(!is.null(modCovar)) repeated <- unique(c(repeated, modCovar))
      arrPheno <- array(NA, dim = c(length(rownam), ncol(pheno)-(1+length(repeated)), nrow(dim3)))
      dimnames(arrPheno) <- list(rownam, (colnames(pheno)[!colnames(pheno) %in% repeated])[-1], as.character(apply(dim3, 1, paste, collapse = "_")))
      for(i in 1:nrow(dim3)){
        vec.bool <- apply(as.matrix(pheno[, colnames(dim3)]) == as.matrix(dim3[rep(i, nrow(pheno)), ]), 1, all)
        arrPheno[as.character(pheno[vec.bool, 1]), , i] <- as.matrix(pheno[vec.bool, (colnames(pheno)[!colnames(pheno) %in% repeated])[-1]])
      }
      if(!is.null(modCovar)){ # strip out covariates
        arrModCovars <- arrPheno[,rep(1, length(modCovar)), ]
        dimnames(arrModCovars)[[2]] <- colnames(pheno)[colnames(pheno) %in% modCovar]
        for(i in 1:nrow(dim3)){
          vec.bool <- apply(matrix(pheno[, colnames(dim3)] == dim3[rep(i, nrow(pheno)), ], ncol=ncol(dim3)), 1, all)
          arrModCovars[as.character(pheno[vec.bool, 1]), , i] <- as.matrix(pheno[vec.bool, colnames(pheno)[colnames(pheno) %in% modCovar]])
        }
      }
    }
    if(!is.null(modCovar)){ # take the correct class of the covariates
      phenoCovars <- arrModCovars
      attrModCovars <- classList[dimnames(arrModCovars)[[2]]]
      for(i in names(attrModCovars)){
        if(attrModCovars[i] != "numeric")
          attrModCovars[ i] <- "factor"
      }
    }
    pheno <- arrPheno
    rm(arrPheno)
  }

  # match geno and pheno
  if(!is.null(geno) & !is.null(pheno)){
    if(is.null(dimnames(pheno)[[1]]) | is.null(rownames(geno))){
      if(dim(pheno)[1] == nrow(geno)){
        warning("assuming identical order of genotypes in 'pheno' and 'geno' because of identical length.\nControll the Output! There is no warranty of correctness!\n")
        if(is.null(dimnames(pheno)[[1]])) dimnames(pheno)[[1]] <- rownames(geno)
        else rownames(geno) <- dimnames(pheno)[[1]]
        if(is.null(dimnames(pheno)[[1]]) & is.null(rownames(geno))) dimnames(pheno)[[1]] <- rownames(geno) <- paste0("ID",10^ceiling(log10(nrow(geno)))+1:nrow(geno))
      }else stop("missing rownames (animal IDs) for either 'pheno' or 'geno' and lengths do not agree.")
      # now geno and pheno have rownames
    }
  }
  # sort geno by rownames (alphabetical order)
  if(!is.null(geno)){
    if(all(row.names(geno) %in% 1:nrow(geno)))
      geno <- geno[order(as.numeric(row.names(geno))), ]
    else
      geno <- geno[order(row.names(geno)),]
  }

  # sort markers by chromosome and position within chromosome
  if(!is.null(map)){
    if (reorderMap){
      map$sor <- substr(map$chr, nchar(as.character(map$chr)), nchar(as.character(map$chr)))
      if(any(unique(map$sor)[!is.na(unique(map$sor))] %in% 0:9)) map$sor <- 1
      # first order by rownames in alphabetical order (important for SNPs with the same position)
      map <- map[order(as.character(rownames(map))),]
      map <- orderBy(~sor+chr+pos,data=map)
      map$sor <- NULL
      # sortcolumns in geno, too
      geno <- geno[,rownames(map)]
    }
    map$pos[is.na(map$chr)] <- NA
    class(map) <- c("GenMap", "data.frame")
  }

  if(!is.null(pedigree)){
    if(!is(pedigree,"pedigree")) warning(paste("object pedigree should be of class 'pedigree', not class", paste(class(pedigree), collapse=' ')))
  }

  # return object
  obj <- list(covar=NULL,pheno=pheno,geno=geno,map=map,pedigree=pedigree,phenoCovars=phenoCovars)

  # add information to element covar
  # sort all available individuals
  ids <- sort(unique(c(dimnames(obj$pheno)[[1]],rownames(obj$geno),as.character(obj$pedigree$ID))))
  if(all(ids %in% 1:length(ids))) ids <- sort(as.numeric(ids))

  obj$covar <- data.frame(id=ids,
                          phenotyped=ids %in% dimnames(obj$pheno)[[1]],
                          genotyped=ids %in% rownames(obj$geno),
                          stringsAsFactors=FALSE)

  # family information for genotyped indviduals
  if(!is.null(family)){
    if(!is.data.frame(family) & !is.matrix(family)) stop('family must be either a data.frame or a matrix, not a ', class(family))
    colnames(family)[1] <- "family"
    family$id <- as.character(rownames(family))
    obj$covar <- merge(obj$covar,family,by="id",all=TRUE)
    obj$covar$genotyped[is.na(obj$covar$genotyped)] <- FALSE
    obj$covar$phenotyped[is.na(obj$covar$phenotyped)] <- FALSE
  } else obj$covar$family <- rep(NA, nrow(obj$covar))

  # add covar from arguments, if available
  if(!is.null(covar)){
    if(is.matrix(covar)){
      if(is.null(rownames(covar))) warning("the supplied covar's rownames will default to 1:nrow(covar), which is likely not correct.  Inspect the resulting $covar.")
      covar <- as.data.frame(covar)
    }
    if(!is.data.frame(covar)) stop('covar must be a data.frame, not a ',class(covar))
    # do not use any existing columns named 'genotyped', 'phenotyped' or 'id'
    covar <- covar[!colnames(covar) %in% c("genotyped","phenotyped","id","family")]
    # merge with existing data
    obj$covar <- merge(obj$covar,covar,by.x=1,by.y=0,all=TRUE)
  }

  if(!is.null(obj$pedigree))
    obj$covar <- obj$covar[match(obj$pedigree$ID, obj$covar$id), ]

  # further information
  obj$info$map.unit <- map.unit
  obj$info$codeGeno <- FALSE
  obj$info$attrPhenoCovars <- attrModCovars
  obj$info$version <- paste("gpData object was created by synbreed version", sessionInfo()$otherPkgs$synbreed$Version)

  # return object of class 'gpData'
  class(obj) <- "gpData"
  return(obj)
}
