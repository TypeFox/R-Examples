###########################################################################################
scan.genome <- function(cross, pheno.col, pheno.parents, addcov, intcov,
                        threshold, method, ...)
{
  ## This currently scans one phenotype at a time.
  ## Want to extend to all phenotypes not in names(covM.dat).
  ## To do that, assume covariates are same for all i.
  
  n.pheno <- length(pheno.col)

  ## *** The following could be simplified.
  ## *** However, it involves a rethinking of myformula.
  ## *** Low priority for now.
  
  ## Design matrices for qtl/scanone. Must be numeric entries.
  ## Design matrix for parent phenotypes
  covM.dat <- pull.pheno.null(cross, pheno.parents)
  ## Design matrix for additive covariates.
  addcovM.dat <- covM.dat
  tmp <- unique(c(addcov, intcov))
  if(!is.null(tmp)) {
    tmp <- create.cov.matrix(cross, cov.names = tmp)
    if(length(tmp))
      addcovM.dat <- cbind(addcovM.dat, tmp)
  }
  ## Design matrix for interactive covariates.
  intcovM.dat <- create.cov.matrix(cross,cov.names=intcov)

  ## Phenotype matrices of actual data. Can be mix of numeric and factor.
  addcov.dat <- pull.pheno.null(cross, unique(c(addcov, intcov)))
  intcov.dat <- pull.pheno.null(cross, intcov)

  ## Call to scanone is the big time commitment.
  ss <- scanone.summary(cross, pheno.col, addcovM.dat, intcovM.dat, threshold, method)
  
  bic <- rep(NA, n.pheno)
  cross.type <- class(cross)[1]
  if(nrow(ss)) {
    ## Model may be different for each trait, depending on QTL.
    ## For loop is clumsy, but calculations are pretty fast.
    for(i in seq(n.pheno)) {
      ## Response for linear model.
      y <- cross$pheno[, pheno.col[i]]
      
      signif.lods <- (ss[ , 1 + 2 * i] >= threshold[i])
      le.markers <- sum(signif.lods)

      ## Build data.frame and formula for linear model fit.
      if(le.markers > 0){
        ## Need to extend this to multiple phenotypes.
        qtlo <- makeqtl(cross, chr = ss[signif.lods, 1],
                        pos = ss[signif.lods, 2 * i], what="prob")
        geno.dat <- hk.design.matrix(qtlo=qtlo, cross.type)[,-1, drop = FALSE]
        
        tmp <- paste("add", 1:le.markers, sep = "")
        if(cross.type == "f2")
          tmp <- as.vector(rbind(tmp, paste("dom", 1:le.markers, sep = "")))
        dimnames(geno.dat) <- list(NULL, tmp)
        
        form <- set.dat.form(y, covM.dat, addcov.dat, intcov.dat, geno.dat, le.markers,
                             cross.type)
        dat <- form$dat
        form <- form$form
      }
      else{
        form <- set.dat.form(y, covM.dat, addcov.dat, intcov.dat, cross.type)
        dat <- form$dat
        form <- form$form
      }
      
      ## Fit linear model.
      fm <- lm(form, dat)

      ## Record BIC.
      bic[i] <- AIC(fm, k = log(length(y)))[1]
    }
  }
  else {
    ## Same model for all traits (no QTL).
    for(i in seq(n.pheno)) {
      if(i == 1) {
        y <- cross$pheno[, pheno.col[1]]

        ## Need only do this once if no QTL.
        form <- set.dat.form(y, covM.dat, addcov.dat, intcov.dat, cross.type)
        dat <- form$dat
        form <- form$form
        
        ## Fit linear model.
        fm <- lm(form, dat)
      }
      else
        dat$y <- cross$pheno[, pheno.col[i]]
        
      ## Record BIC.
      bic[i] <- AIC(update(fm, data = dat), k = log(length(y)))[1]
    }
  }
  bic
}
###########################################################################################
scanone.summary <- function(cross, pheno.col, addcov, intcov, threshold,
                            method)
{
  ## This is the big time commitment.
  scan <- scanone(cross, pheno.col = pheno.col,
                  addcovar = addcov, intcovar = intcov, method = method)
  
  ## Mainly intersted in this summary to determine QTLs.
  summary(scan, format = ifelse(length(pheno.col) == 1, "onepheno", "allpeaks"),
                threshold = threshold)
}
###########################################################################################
set.dat.form <- function(y, covM.dat=NULL, addcov.dat=NULL, intcov.dat=NULL,
                         geno.dat=NULL, le.markers = 0, cross.type = "f2")
{
  ## Set up data.frame
  dat <- data.frame(y = y)
  if(!is.null(covM.dat))
    dat <- cbind.data.frame(dat, covM.dat)
  if(!is.null(addcov.dat))
    dat <- cbind.data.frame(dat, addcov.dat)
  if(!is.null(intcov.dat))
    dat <- cbind.data.frame(dat, intcov.dat)
  if(!is.null(geno.dat))
    dat <- cbind.data.frame(dat, geno.dat)

  ## Set up formula.
  form <- cov.formula(c(names(covM.dat),names(addcov.dat)),
                      names(intcov.dat),
                      le.markers, cross.type)

  list(dat = dat, form = form)
}
######################################################################
hk.design.matrix <- function(qtlo, cross.type="f2")
{
  nr <- nrow(qtlo$prob[[1]])
  ng <- length(qtlo$prob)
  if(cross.type == "f2"){
    tmp <- unlist(lapply(qtlo$prob,
                         function(x) {
                           if(ncol(x) == 3)
                             cbind(x[,1]-x[,3],x[,2])
                           else ## Must be X chr.
                             x[,1]-x[,2]
                         }))
    hkm <- matrix(tmp,nr,2*ng)
  }
  if(cross.type == "bc"){
    tmp <- unlist(lapply(qtlo$prob, function(x) x[,1]-x[,2]))
    hkm <- matrix(tmp,nr,ng)
  }
  cbind(rep(1,nr),hkm)
}
