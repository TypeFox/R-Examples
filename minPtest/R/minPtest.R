minPtest <-
function(y,x,SNPtoGene,formula=NULL,cov=NULL,matchset=NULL,permutation=1000,
         seed=NULL,subset=NULL,parallel=FALSE,ccparallel=FALSE,trace=FALSE,
         aggregation.fun=min,adj.method=c("bonferroni","holm","hochberg","hommel","BH","BY","fdr","none"), ...){
  call <- match.call()
  adj.method = match.arg(adj.method)
  if (!is.null(subset)) {
    y <- y[subset]
    x <- x[subset, , drop=FALSE]
    if (!is.null(cov)) {cov <- cov[subset, , drop=FALSE]}
    if (!is.null(matchset)) {matchset <- matchset[subset]}
  }
  if(any(is.na(y))) {warning("response vector y includes NAs, subjects are excluded")}
  if(ncol(SNPtoGene)!=2){
    stop("mapping matrix SNPtoGene should have two columns")}
  if(any(is.na(SNPtoGene[,2]))) {stop("mapping matrix SNPtoGene includes NAs")}
  if(is.null(colnames(x))){
    warning("SNP matrix x has no colnames, not comparable to mapping matrix SNPtoGene")}
  if(any(!(SNPtoGene[,1]%in%colnames(x)))){
    stop("different names of SNPs in SNP matrix x and mapping matrix SNPtoGene")}
  if(any(!(colnames(x)%in%SNPtoGene[,1]))){
    stop("different names of SNPs in SNP matrix x and mapping matrix SNPtoGene")}
  if(nrow(x)!=length(y)){
    stop("length of response vector y should equal the number of rows of SNP matrix x")}
  snp.miss <- sum(is.na(x))
  if(!is.null(matchset)){
    if(length(y)!=length(matchset))
      {stop("length of response vector y should equal the length of the  matchset vector")}
    if(any(is.na(matchset))) {warning("matchset vector includes NA, subjects are excluded")}
    resp.match <- data.frame(y=y,matchset=matchset)
    matchnr.cases <- resp.match[resp.match[,1]==1,2]
    matchnr.controls <- resp.match[resp.match[,1]==0,2]
    cases.match <- lapply(seq_along(matchnr.cases), function(i){
      if(!is.na(matchnr.cases[i])){
        z <- which(matchnr.controls==matchnr.cases[i])
        z.miss <- length(z)>=1
        z.miss
      }
    })
    if(any(!unlist(cases.match))) {warning("not every case has controls")}
    control.match <- lapply(seq_along(matchnr.controls), function(i){
      if(!is.na(matchnr.controls[i])){
        z <- which(matchnr.cases==matchnr.controls[i])
        z.miss <- length(z)>=1
        z.miss
      }
    })
    if(any(!unlist(control.match))) {warning("not every control belongs to a case")}
  }
  nrsnp <- dim(x)[2]
  n <- length(y)
  nrgene <- length(unique(SNPtoGene[,2]))
  if(ccparallel){
    if(!require("snowfall")) {
      stop("package 'snowfall' must be installed and loaded, otherwise parallelization cannot be performed")
    }else{
      require("snowfall")
      if(!sfIsRunning()) sfInit(parallel=TRUE)
    }
  }
  if(parallel){
    if(!require("parallel")){
      stop("package 'parallel' must be installed and loaded, otherwise parallelization cannot be performed")
    }else{
      require("parallel")}}
  if(parallel & ccparallel){
    warning("parallelization is performed using 'snowfall'")
  }
##################################Cochran-Armitage trend test###################################################
  if(is.null(formula)){
    method <- "Cochran-Armitage Trend Test"
    if(!is.null(cov) | !is.null(matchset)){
      warning("missing formula, Cochran-Armitage Trend Test is performed")
    }
    if(snp.miss==0){
      if(length(unique(as.factor(x)))!=3){
        warning("SNP matrix x should contain three features")}
    }
    if(snp.miss!=0){
      if(length(unique(as.factor(x)))!=4){
        warning("SNP matrix x should contain three features")}
    }
    unique_snps <- unique(x)
    valmax <- max(unique_snps, na.rm=TRUE)
    valmin <- min(unique_snps, na.rm=TRUE)
    tsnps <- t(x)
    if(ccparallel){
      sfLibrary("scrime", character.only=TRUE)
      sfExport("valmax","valmin","tsnps")
    }
    cases <- rowTables(tsnps[, y==1, drop=FALSE],levels=valmin:valmax)
    controls <- rowTables(tsnps[, y==0, drop=FALSE],levels=valmin:valmax)
    p_value <- as.matrix(rowCATTs(cases,controls)$rawp)
    colnames(p_value) <- "p_value"
    perrfunc <- function(permut){
      if(trace) {cat("permutation",permut,"\n")}
      set.seed(seed[permut])
      permy <- sample(y)
      casesperm <- rowTables(tsnps[, permy==1, drop=FALSE],levels=valmin:valmax)
      controlsperm <- rowTables(tsnps[, permy==0, drop=FALSE],levels=valmin:valmax)
      pperr <- as.matrix(rowCATTs(casesperm,controlsperm)$rawp)
      return(pperr)
    }
##################################univariate logistische regression#############################################
  }else{
    snpvecnames <- colnames(x)
    if(!is.null(cov) & as.character(formula)[3]==1){
      warning("no specification of covariables in the formula, logistic regression is performed without covariables")
    }
    if(ccparallel){
      sfExport("formula","snpvecnames")
      if(!is.null(matchset)){
        sfExport("n","matchset")
        sfLibrary("Epi", character.only=TRUE)
      }
    }
##### logistic regression without covariables
    if(is.null(cov)){
      dat <- as.data.frame(cbind(y,x))
      colnames(dat)[1] <- as.character(formula)[2]
      if(as.character(formula)[3]!=1)
        {warning("no covariable matrix, logistic regression is performed without covariables")}
#unconditional logistic regression without covariables
      if(is.null(matchset)){
        method <- "unconditional logistic regression (glm)"
        actual.frame <- dat[,colnames(dat)[1], drop=FALSE]
        if(ccparallel){ sfExport("dat", "actual.frame") }
        p_valfunc <- function(i) {
          actual.frame$SNP <- dat[,snpvecnames[i]]
          new.form <- as.formula(paste(as.character(formula)[2],"SNP", sep = "~"))
          fit <- glm(new.form, binomial, actual.frame)
          fit.res <- summary(fit)$coefficients
          fit.res <- fit.res[dim(fit.res)[1],dim(fit.res)[2]]
          names(fit.res) <- snpvecnames[i]
          fit.res
        }
        perrfunc <- function(permut){
          if(trace) {cat("permutation", permut , "\n")}
          set.seed(seed[permut])
          permy <- sample(y)
          actual.frame$permy <- permy
          pperr <- apply(matrix(seq_along(snpvecnames)), MARGIN=1, function(i) {
            actual.frame$SNP <- dat[,snpvecnames[i]]
            new.form <- as.formula(paste("permy", "SNP", sep = "~"))
            fit <- glm(new.form, binomial, actual.frame)
            fit.res <- summary(fit)$coefficients
            fit.res <- fit.res[dim(fit.res)[1],dim(fit.res)[2]]
            fit.res
          })
          names(pperr) <- snpvecnames
          pperr <- as.matrix(pperr)
          return(pperr)
        }
      }else{
#conditional logistic regression without covariables
        method <- "conditional logistic regression"
        dat$matchset <- matchset
        actual.frame <- dat[,c(colnames(dat)[1],"matchset")]
        if(ccparallel){ sfExport("dat", "actual.frame") }
        p_valfunc <- function(i) {
          actual.frame$SNP <- dat[,snpvecnames[i]]
          new.form <- as.formula(paste(as.character(formula)[2],"SNP", sep = "~"))
          fit <- clogistic(new.form, strata=matchset, data=actual.frame)
	  coef <- fit$coefficients
          se <- sqrt(diag(fit$var))
          p <- 1 - pchisq((coef/se)^2, 1)
          fit.res <- p[length(p)]
          names(fit.res) <- snpvecnames[i]
          fit.res
        }
        perrfunc <- function(permut){
          if(trace) {cat("permutation", permut, "\n")}
          set.seed(seed[permut])
          perm <- sample(n)
          permy <- y[perm]
          permmatch <- matchset[perm]
          actual.frame$permy <- permy
          actual.frame$permmatch <- permmatch
          pperr <- apply(matrix(seq_along(snpvecnames)), MARGIN=1, function(i) {
            actual.frame$SNP <- dat[,snpvecnames[i]]
            new.form <- as.formula(paste("permy", "SNP", sep = "~"))
            fit <- clogistic(new.form, strata=permmatch, data=actual.frame)
	    coef <- fit$coefficients
	    se <- sqrt(diag(fit$var))
            p <- 1 - pchisq((coef/se)^2, 1)
            fit.res <- p[length(p)]
            fit.res
          })
          names(pperr) <- snpvecnames
          pperr <- as.matrix(pperr)
          return(pperr)
        }
      }
    }else{
##### logistic regression with covariable
      if(nrow(cov)!=length(y)){
        stop("length of response vector y should equal the number of rows of covariate matrix")
      }
      if(is.null(colnames(cov))){
        warning("no colnames in cov, not comparable to covariate names in formula")
      }
      w <- length(which(unlist(strsplit(as.character(formula)[3]," "))!="+"))
      if(ncol(cov)!=w){warning("number of covariates in formula differs from covariates in covariate matrix")}
      dat <- as.data.frame(cbind(y,cov,x))
      colnames(dat)[1] <- as.character(formula)[2]
      perr.formula <- paste("permy",as.character(formula)[3],sep="~")
      if(ccparallel){ sfExport("perr.formula") }
#unconditional logistic regression with covariables
      if(is.null(matchset)){
        method <- "unconditional logistic regression (glm)"
        actual.frame <- dat[,c(colnames(dat)[1], colnames(cov))]
        if(ccparallel) { sfExport("dat", "actual.frame") }
        p_valfunc <- function(i) {
          actual.frame$SNP <- dat[,snpvecnames[i]]
          new.form <- as.formula(paste(deparse(formula),"SNP", sep = "+"))
          fit <- glm(new.form, binomial, actual.frame)
          fit.res <- summary(fit)$coefficients
          fit.res <- fit.res[dim(fit.res)[1],dim(fit.res)[2]]
          names(fit.res) <- snpvecnames[i]
          fit.res
        }
        perrfunc <- function(permut){
          if(trace) {cat("permutation" , permut, "\n")}
          set.seed(seed[permut])
          permy <- sample(y)
          actual.frame$permy <- permy
          pperr <- apply(matrix(seq_along(snpvecnames)), MARGIN=1, function(i) {
            actual.frame$SNP <- dat[,snpvecnames[i]]
            new.form <- as.formula(paste(perr.formula, "SNP", sep = "+"))
            fit <- glm(new.form, binomial, actual.frame)
            fit.res <- summary(fit)$coefficients
            fit.res <- fit.res[dim(fit.res)[1],dim(fit.res)[2]]
            fit.res
          })
          names(pperr) <- snpvecnames
          pperr <- as.matrix(pperr)
          return(pperr)
        }
      }else{
# conditional logistic regression with covariables
        method <- "conditional logistic regression (clogistic)"
        dat$matchset <- matchset
        actual.frame <- dat[,c(colnames(dat)[1], colnames(cov), "matchset")]
        if(ccparallel) { sfExport("dat", "actual.frame") }
        p_valfunc <- function(i) {
          actual.frame$SNP <- dat[,snpvecnames[i]]
          new.form <- as.formula(paste(deparse(formula), "SNP", sep = "+"))
          fit <- clogistic(new.form, strata=matchset, data=actual.frame)
          coef <- fit$coefficients
          se <- sqrt(diag(fit$var))
          p <- 1 - pchisq((coef/se)^2, 1)
          fit.res <- p[length(p)]
          names(fit.res) <- snpvecnames[i]
          fit.res
        }
        perrfunc <- function(permut){
          if(trace) {cat("permutation", permut , "\n")}
          set.seed(seed[permut])
          perm <- sample(n)
          permy <- y[perm]
          permmatch <- matchset[perm]
          actual.frame$permy <- permy
          actual.frame$permmatch <- permmatch
          pperr <- apply(matrix(seq_along(snpvecnames)), MARGIN=1, function(i) {
            actual.frame$SNP <- dat[,snpvecnames[i]]
            new.form <- as.formula(paste(perr.formula, "SNP", sep = "+"))
            fit <- clogistic(new.form, strata=permmatch, data=actual.frame)
	    coef <- fit$coefficients
	    se <- sqrt(diag(fit$var))
            p <- 1 - pchisq((coef/se)^2, 1)
            fit.res <- p[length(p)]
            fit.res
          })
          names(pperr) <- snpvecnames
          pperr <- as.matrix(pperr)
          return(pperr)
        }
###################################################################################
      }
    }
    if(ccparallel){
      p_value <- sfClusterApplyLB(seq_along(snpvecnames),p_valfunc)
    }else{
      if(parallel){
        if(parallel>1){
          p_value <- mclapply(seq_along(snpvecnames),p_valfunc,mc.preschedule=FALSE,mc.cores=parallel)
        }else{
          p_value <- mclapply(seq_along(snpvecnames),p_valfunc,mc.preschedule=FALSE)
        }
      }else{
        p_value <- lapply(seq_along(snpvecnames),p_valfunc)
      }
    }
    p_value <- as.matrix(unlist(p_value))
    colnames(p_value) <- "p_value"
  }
  if(is.null(seed)){
    seed <- sample(1:1e7,size=permutation)
  }else{
    if(length(seed)!=permutation){
      stop("length of seed vector should equal the number of permutations")
    }
  }
  if(ccparallel){
    sfExport("y","seed","permutation","trace")
    p_valuesperr <- sfClusterApplyLB(1:permutation, perrfunc)
    sfStop()
  }else{
    if(parallel){
      if(parallel>1){
        p_valuesperr <- mclapply(1:permutation,perrfunc,mc.preschedule=FALSE,mc.cores=parallel)
      }else{
        p_valuesperr <- mclapply(1:permutation,perrfunc,mc.preschedule=FALSE)
      }
    }else{
      p_valuesperr <- lapply(1:permutation,perrfunc)
    }
  }
  psnpperr <- matrix(unlist(p_valuesperr),ncol=permutation,nrow=nrsnp)
  rownames(psnpperr) <- colnames(x)
  colnames(psnpperr) <- c(1:permutation)
  genes <- unique(SNPtoGene[,2])
  pgen <- lapply(seq_along(genes), function(i){
    snpnames <- SNPtoGene[which(SNPtoGene[,2]==genes[i]),1]
    minpgen <- aggregation.fun(p_value[snpnames,], ...)
    minpgen
  })
  names(pgen) <- genes
  pgen <- as.matrix(unlist(pgen))
  colnames(pgen) <- "test statistic"
  pgenperr <- lapply(p_valuesperr, function(x){
    p_list <- x
    minpper <- apply(matrix(seq_along(genes)), MARGIN=1, FUN=function(i){
      snpnames <- SNPtoGene[which(SNPtoGene[,2]==genes[i]), 1]
      minpgenperr <- aggregation.fun(p_list[snpnames,], ...)
      minpgenperr
    })
  })
  pgenperr <- matrix(unlist(pgenperr),ncol=permutation,nrow=length(genes))
  colnames(pgenperr) <- c(1:permutation)
  rownames(pgenperr) <- genes
  minp <- double(length(pgen))
  for(i in seq_along(pgen)){
    minp[i] <- mean(pgenperr[i,]<=pgen[i])    
  }
  names(minp) <- genes
################################################################################################################
## mutilple testing correction
  p.adj.minp <- matrix(p.adjust(p=minp, method=adj.method, n=nrgene),nrow=length(minp),ncol=1)
  rownames(p.adj.minp) <- genes
  colnames(p.adj.minp) <- c("p.adjust")
###############################################################################################################3
  minp <- as.matrix(minp)
  colnames(minp) <- "minP"
  p.adj.psnp <- matrix(p.adjust(p=as.vector(p_value), method=adj.method, n=nrsnp),nrow=length(p_value), ncol=1)	
  colnames(p.adj.psnp) <- c("p.adjust")
  rownames(p.adj.psnp) <- rownames(p_value)
  fit <- list(call=call,n=n,nrsnp=nrsnp,nrgene=nrgene,
              snp.miss=snp.miss,method=method,n.permute=permutation,
              SNPtoGene=SNPtoGene,psnp=p_value,psnpperm=psnpperr,
              zgen=pgen,zgenperm=pgenperr,minp=minp,
              p.adj.psnp=p.adj.psnp,p.adj.minp=p.adj.minp)
  class(fit) <- "minPtest"
  fit
}

print.minPtest <- function(x, ...){
  cat("\nUsed method:",x$method, "for", x$n, "subjects")
  cat("\nCall: ",deparse(x$call),"\n\n")
  cat("Number of genes:", x$nrgene, "\n")
  cat("Number of SNPs:", x$nrsnp, "\n")
  cat("Number of missings in the SNP matrix:", x$snp.miss, "\n")
  cat("Number of permutations:", x$n.permute, "\n")
  invisible()
}
