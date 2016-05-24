setToBeDefinedPkgs <- function(pkgs, keywordstyles){
    if(!is.character(pkgs)) stop("Argument 'pkgs' must be a character")
    lP <- length(pkgs)
    if(missing(keywordstyles))
        keywordstyles <- getSweaveListingOption("Keywordstyle")
    keywordstyles <- rep(keywordstyles, length.out = lP)
    if(lP){
       tobeDefinedPkgs <- matrix(character(lP*2), nrow=lP, ncol =2)
       tobeDefinedPkgs[,1]  <- pkgs
       tobeDefinedPkgs[,2]  <- keywordstyles
       tbd <- if(is.null(.tobeDefinedPkgs)) tobeDefinedPkgs else
                 rbind(.tobeDefinedPkgs,tobeDefinedPkgs)
       .tobeDefinedPkgs <<- tbd
    }
    return(invisible())
}

lookUpKeywordStyles <- function(pkgs, defkws){
    lP <- length(pkgs)
    if(lP){
      kws <- character(lP)
      idx <- seq(along = pkgs)
      idxIn <- idx[pkgs %in% .tobeDefinedPkgs[,1]]
      idxNIn <- idx[!(pkgs %in% .tobeDefinedPkgs[,1])]
      idxRec <- NULL
      if(length(idxNIn)){
        if(missing(defkws))
           defkwsI <- setBaseOrRecommended(pkgs = pkgs[idxNIn])
        else{
           defkwsI <- if(length(idxNIn)<=length(defkws)) defkws[idxNIn] else defkws
           defkwsI <- defkwsI[!is.na(defkwsI)]
           defkwsI <- rep(defkwsI, length.out = length(idxNIn))
           defkwsI <- setBaseOrRecommended(pkgs = pkgs[idxNIn], kws = defkwsI)
        }
        kws[idxNIn] <- defkwsI
      }
      if(length(idxIn))
         kws[idxIn] <- .tobeDefinedPkgs[.tobeDefinedPkgs[,1] %in% pkgs[idxIn] ,2]
      return(kws)
    }
}

lstsetLanguage <- function(pkgs, posIdx, keywordstyles, overwrite = FALSE){
     genKWL <- function(pkg, kwd, kws, withcomma = TRUE){
        if (!overwrite)
             kwd <- kwd[!kwd%in%.keywordsR]
        num <- length(.alreadyDefinedPkgs)
        BaseOrRec <- ""
        if(!is.na(isBaseOrRecommended(pkg)))
            if(isBaseOrRecommended(pkg))
               BaseOrRec <- " [recommended or base] "
        cat("% --------------------------\n% Registration of package ",pkg,
             BaseOrRec,
            "\n% --------------------------\n",
            "morekeywords={[",num+2,"]", sep = "")
        ml <- length(kwd); m1 <- ml%/%5; m2 <- ml%%5
        kwd[ml] <- paste(kwd[ml],"%",sep="")
        if(length(kwd)){
           if(m2!=1){
              seps <-  c(rep(c(rep(",",4),",%\n"),ifelse(m2==0,m1-1,m1)),
                         c(rep(",",ifelse(m2==0, 4, m2-1))),"%")
           }else{
              seps <- c(rep(c(rep(",",4),",%\n"),m1),"%")
           }
           cat(kwd,sep = seps,  fill = FALSE); #cat("%\n")
           if(ml<=5) cat("\n")
           cat("},%\nkeywordstyle={[",num+2,"]",kws,"}",sep="")
           if(withcomma) cat(",")
           cat("%\n%\n", sep = "")
#           cat(paste("%\n%\n\n"))
           .alreadyDefinedPkgs <<- c(.alreadyDefinedPkgs,pkg)
        }
        return(invisible())
     }
     alreadyDefinedPkgs <- .alreadyDefinedPkgs
     seL <- .packages()
     if(!missing(pkgs)){
       posIdx0 <- sapply(pkgs, function(x){
                         w <- which(seL==x)
                         if( length(w)&& !x %in%alreadyDefinedPkgs)
                             return(w) else return(NA)}
                         )
       pkgs <- pkgs[!is.na(posIdx0)]
       posIdx0 <- posIdx0[!is.na(posIdx0)]
       if(missing(posIdx)) posIdx <- posIdx0
          else {posIdx <- posIdx[! .packages()[posIdx] %in% alreadyDefinedPkgs]
                posIdx <- unique(c(posIdx, posIdx0))}
     }else{
       if(missing(posIdx)) posIdx <- 1 : length(seL)
       posIdx <- posIdx[! .packages()[posIdx] %in% alreadyDefinedPkgs]
     }
     pkgs <- .packages()[posIdx]
     lP <- length(posIdx)
     if(missing(keywordstyles)){
#        keywordstyles <- getSweaveListingOption("Keywordstyle")
#        keywordstyles <- rep(keywordstyles, length.out = lP)
        keywordstyles <- setBaseOrRecommended(pkgs = pkgs)# ,
#                                              kws = keywordstyles)
     }else
        keywordstyles <- rep(keywordstyles, length.out = lP)

     if(lP) {
        if(is.null(.numberofRstyleDefs)) .numberofRstyleDefs <<- 1
        else .numberofRstyleDefs <<- .numberofRstyleDefs + 1
        withcomma <- TRUE
        pos1 <- 0
        for(i in 1: lP){
           pos0 <- which(pkgs[i] == gsub("package:","",search()))
           kwl <- ls(pos = pos0)
           kwl <- sort(kwl[grep("^[[:alpha:]]+\\w*",kwl,perl=TRUE)],
                       decreasing = TRUE)

           if(i==1){
              pos1 <- pos0 + length(.alreadyDefinedPkgs)
              cat("\n\\lstdefinestyle{RstyleO",.numberofRstyleDefs+1,"}{style=RstyleO",
                       .numberofRstyleDefs,",%\n",sep="")
           }

           if(i==lP) withcomma <- FALSE
           genKWL(pkg = pkgs[i], kwd = kwl,
                  kws = keywordstyles[i], withcomma = withcomma)
           if(i==lP) cat("}%\n\\lstdefinestyle{Rstyle}{style=RstyleO",
                          .numberofRstyleDefs+1,"}\n\n",
                          sep="")
           }
     }
     return(invisible())
}

changeKeywordstyles <- function(pkgs, keywordstyles){
     setkws <- function(num, kws){
        cat("%\n\\lstdefinestyle{RstyleO",.numberofRstyleDefs+1,"}{style=RstyleO",
                       .numberofRstyleDefs,",%\n", sep="")
        cat("{keywordstyle={[",num+1,"]",kws,"}\n}\n", sep = "")
        cat("\\lstdefinestyle{Rstyle}{style=RstyleO", .numberofRstyleDefs+1,"}\n",
            sep = "")
        }

     alreadyDefinedPkgs <- .alreadyDefinedPkgs

     if(!missing(pkgs)){
       numIdx0 <- sapply(pkgs, function(x){
                         w <- which(alreadyDefinedPkgs==x)
                         if( length(w))
                             return(w) else return(NA)}
                         )
       pkgs <- pkgs[!is.na(numIdx0)]
       numIdx <- numIdx0[!is.na(numIdx0)]
     }else{
       pkgs <- alreadyDefinedPkgs
       numIdx <- seq(along = length(pkgs))
     }

     lP <- length(pkgs)

     if(missing(keywordstyles))
        keywordstyles <- getSweaveListingOption("Keywordstyle")
     keywordstyles <- rep(keywordstyles, length.out = lP)

     if(lP){
        if(is.null(.numberofRstyleDefs)) .numberofRstyleDefs <<- 1
        else .numberofRstyleDefs <<- .numberofRstyleDefs + 1
        for(i in 1: lP)
            setkws(num = numIdx[i], kws = keywordstyles[i])
     }
     return(invisible())
}

require <- function(package, lib.loc = NULL, quietly = FALSE,
             warn.conflicts = TRUE,
             keep.source = getOption("keep.source.pkgs"),
             character.only = FALSE, save = TRUE, inSweave,
             keywordstyles, interm.keywordstyles, overwrite, intermediate){
             pkg <- as.character(substitute(package))
             mc <- as.list(match.call(expand.dots = FALSE))[-1]
             rerg <- .LibOrRequire(pkg, mc, base::require)
             return(invisible(rerg))
             }

library <- function(package, help, pos = 2, lib.loc = NULL,
    character.only = FALSE,
    logical.return = FALSE, warn.conflicts = TRUE,
    quietly = FALSE, keep.source = getOption("keep.source.pkgs"),
    verbose = getOption("verbose"), inSweave,
    keywordstyles, interm.keywordstyles, overwrite, intermediate){
             pkg <- as.character(substitute(package))
             mc <- as.list(match.call(expand.dots = FALSE))[-1]
             rerg <- .LibOrRequire(pkg, mc, base::library)
             if(class(rerg) == "libraryIQR" | class(rerg) == "packageInfo") return(rerg)
             if(logical.return) return(rerg)
             return(invisible(rerg))
    }

.LibOrRequire <- function(pkg, mc, fct){
             keywordstyles <- mc$keywordstyles
             interm.keywordstyles <- mc$"interm.keywordstyles"
             overwrite <- mc$overwrite
             intermediate <- mc$intermediate
             inSweave <- mc$inSweave
             if(is.null(overwrite))
                overwrite <-  getSweaveListingOption("overwrite")
             if(is.null(intermediate))
                intermediate <-  getSweaveListingOption("intermediate")
             if(is.null(inSweave))
                inSweave <-  getSweaveListingOption("inSweave")
             mc$keywordstyles <- NULL
             mc$"interm.keywordstyles" <- NULL
             mc$overwrite <- NULL
             mc$intermediate <- NULL
             mc$inSweave <- NULL
#
             if(inSweave){
                pold <- .packages()
                rerg <- suppressPackageStartupMessages(do.call(fct, args=mc))
                pnew <- .packages()
                pdiff <- pnew[! pnew %in% pold]
                pint <- pdiff[!pdiff %in% pkg]
#
                if(is.null(interm.keywordstyles) && length(pint))
                    interm.keywordstyles <- lookUpKeywordStyles(pkgs = pint,
                                     defkws = getSweaveListingOption("interm.Keywordstyle"))
                if(is.null(keywordstyles))
                    keywordstyles <- lookUpKeywordStyles(pkgs = pkg,
                                     defkws = getSweaveListingOption("Keywordstyle"))

                if(intermediate){
                   pkgs <- c(pkg, pint)
                   kws.int <- rep(interm.keywordstyles,length.out= length(pint))
                   kws.o <- rep(keywordstyles,length.out= length(pkg))
                   kws <- c(kws.o,kws.int)
                }else{
                   pkgs <- pkg
                   kws.o <- rep(keywordstyles,length.out= length(pkg))
                   kws <- kws.o
                }
                if(is.null(mc$help)&&!is.null(pkg))
                lstsetLanguage(pkgs = pkgs, keywordstyles = kws,
                               overwrite = overwrite)
             }else{
                rerg <- do.call(fct, args=mc)
             }
             return(rerg)
}

isBaseOrRecommended <- function(pkgs){
   sapply(pkgs, function(x){
      logres <- read.dcf(file = system.file("DESCRIPTION", package = x),
                         fields="Priority") %in% c("base","recommended")
      logres[is.na(logres)] <- FALSE
      any(logres)

   })
}


setBaseOrRecommended <- function(pkgs,kws =NULL){
      idxRec <- as.logical(isBaseOrRecommended(pkgs = pkgs))
      if(missing(kws)) {
         kws <- rep(getSweaveListingOption("Keywordstyle"),
                    length.out = length(pkgs))
         if(length(idxRec))
            kws[idxRec] <- rep(getSweaveListingOption("Recomd.Keywordstyle"),
                               length.out = length(idxRec))
      }else
         kws <- rep(kws, length.out = length(pkgs))

      if(any(idxRec)){
        kwsI <- kws[idxRec]

        grI <- grep("\\\\color\\{.+\\}", kwsI, invert = TRUE)
        kwsI <- gsub("\\\\color\\{.[^\\}]+\\}","\\\\color\\{RRecomdcolor\\}",
                           kwsI)
        kwsI[grI] <- gsub("\\}$","\\\\color\\{RRecomdcolor\\}\\}",kwsI[grI])
        kws[idxRec] <- kwsI
      }
      return(kws)
}
