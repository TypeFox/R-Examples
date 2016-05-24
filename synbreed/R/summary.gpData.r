# summary method fo class 'gpData'
summary.gpData <- function(object,...){
    obj <- object
    ans <- list()

    # summary of 'covar'
    ans$covar$n <- length(obj$covar$id)
    ans$covar$nphenotyped <- sum(obj$covar$phenotyped)
    ans$covar$ngenotyped <- sum(obj$covar$genotyped)


    # summary of 'pheno'
    ans$reps <- NULL
    if(is.null(obj$pheno)){
      ans$pheno <- NULL 
    } else {
      ans$pheno <- summary(data.frame(apply(obj$pheno, 2, as.numeric)))
      if(dim(obj$pheno)[[3]]>1) ans$reps <- dim(obj$pheno)[[3]]
    }
    ans$phenoCovars <- NULL
    if(!is.null(obj$phenoCovars)){
      df.phenoCovars <- as.data.frame(apply(obj$phenoCovars, 2, rbind))
      for(i in names(obj$info$attrPhenoCovars))
          df.phenoCovars[, i] <- as(df.phenoCovars[, i], obj$info$attrPhenoCovars[i])
       ans$phenoCovars <- summary(df.phenoCovars)
      ans$phenoCovarsMode <- unlist(lapply(df.phenoCovars, class))
      rm(df.phenoCovars)
    }
    
    #summary of 'geno'
    if(is.null(obj$geno)) ans$geno <- NULL
    else{
      geno <- obj$geno
      nobs <- ncol(geno)*nrow(geno)
      # table is very time consuming
      if(ncol(geno)*nrow(geno)<4e+06){
        ans$geno <- list(nMarkers=ncol(geno),genotypes=table(geno)/nobs,nNA=sum(is.na(geno))/nobs)
      }
      else{
         genotypes <- unique(as.vector(geno))
         genotypes <- genotypes[!is.na(genotypes)]
         frequencies <- rep(NA,length(genotypes))
         for (i in 1:length(genotypes)){
             frequencies[i] <- sum(geno==genotypes[i],na.rm=TRUE)/nobs
         }
         names(frequencies) <- genotypes
         ans$geno <- list(nMarkers=ncol(geno),genotypes=round(frequencies,3),nNA=sum(is.na(geno),na.rm=TRUE)/nobs)
      } 
      if(!is.null(obj$map)){
         ans$geno$markerChr <- table(obj$map$chr)#[unique(obj$map$chr)]     # keep same order as in map
         ans$geno$markerChr <- ans$geno$markerChr[!is.na(ans$geno$markerChr)]
         mapped <- !(is.na(obj$map$chr) | is.na(obj$map$pos))
         ans$geno$mappedMarkers <- sum(mapped)
      }


    }

    # summary of 'pedigree'
    if(is.null(obj$pedigree)) ans$pedigree <- NULL
    else ans$pedigree <- summary.pedigree(obj)
    
    class(ans) <- "summary.gpData"
    ans
}

# print method for summary method fo class 'gpData'
print.summary.gpData <- function(x,...){
    cat("object of class 'gpData' \n")
    cat("covar \n")
    cat("\t No. of individuals",x$covar$n,"\n" )
    cat("\t         phenotyped",x$covar$nphenotyped,"\n")
    cat("\t          genotyped",x$covar$ngenotyped,"\n")
    cat("pheno \n")
    cat("\t No. of traits:\t\t  ",dim(x$pheno)[[2]],"\n" )
    if(!is.null(x$reps)) cat("\t No. of repeated measures:", x$reps, "\n")
    cat("\n")
    print(x$pheno)
    cat("\n")
    if(!is.null(x$phenoCovars)){
      cat("phenoCovars \n")
      cat("\t No. of covariables: ",dim(x$phenoCovars)[[2]],"\n" )
      cat("\t Class of covariables: \n" )
      print(x$phenoCovarsMode)
      cat("\n")
      print(x$phenoCovars)
      cat("\n")
    }
    cat("geno \n")
    cat("\t No. of markers",x$geno$nMarkers,"\n")
    cat("\t genotypes",names(x$geno$genotypes),"\n")
    cat("\t frequencies",x$geno$genotypes,"\n")
    cat("\t NA's",format(x$geno$nNA*100,digits=3,nsmall=3),"%\n")
    cat("map \n")
    cat("\t No. of mapped markers ",x$geno$mappedMarkers,"\n")
    cat("\t No. of chromosomes    ",length(x$geno$markerChr),"\n\n")
    cat("\t markers per chromosome \n\t")
    print(x$geno$markerChr)
    cat("\n")
    cat("pedigree \n")
    print(x$pedigree)
}
