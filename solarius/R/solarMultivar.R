# Run multivariate analysis and report results.
#
# @export
solarMultivar <- function(formula, data, traits, covlist = "1",
  na.rm = FALSE,
  ...,
  verbose = 0) 
{
  # Example: M <- solarMultivar(traits = c("trait1", "trait2"), covlist = list(c("SEX"), c("AGE")), data = dat40, na.rm = TRUE, verbose = 1)

  mc <- match.call()

  stopifnot(requireNamespace("Matrix"))
    
  ### Step 1: parse arguments
  if(!missing(formula)) {
    formula.str <- as.character(formula)

    traits <- formula.str[2]
    traits <- strsplit(traits, "\\+")[[1]]
    traits <- gsub(" ", "", traits)
    
    formula.polygenic <- formula
    
    traits.polygenic <- traits
    covlist.polygenic <- as.character(NA)    
  } else {
    stopifnot(!missing(traits))
    
    traits.polygenic <- traits
    
    if(length(covlist) == 1) {
      formula.polygenic <-  paste(paste(traits, collapse = "+"), "~", 
      paste(covlist, collapse = "+"))

      covlist.polygenic <- llply(1:length(traits), function(i) covlist)
      formula.polygenic <- as.formula(formula.polygenic)
    } else {
      stopifnot(length(covlist) == length(traits))
      
      covlist.polygenic <- covlist
      
      cov.polygenic <- llply(1:length(covlist), function(i) {
        cov <- covlist[[i]]
        cov <- cov[!(cov %in% c("1", ""))]
        
        paste0(cov, "(", traits[i], ")")
      })
      cov.polygenic <- unlist(cov.polygenic)

      if(length(cov.polygenic) == 0) {
        cov.polygenic <- "1"
      }

      formula.polygenic <-  paste(paste(traits, collapse = " + "), "~", 
      paste(cov.polygenic, collapse = " + "))

      formula.polygenic <- as.formula(formula.polygenic)
    }
  }
  
  ntraits <- length(traits.polygenic)
  stopifnot(ntraits == 2)
  
  # ID vars
  id.newnames <- matchIdNames(names(data), skip.sex = TRUE)
  # > id.newnames
  #     id   famid      mo      fa 
  #    "ID" "FAMID"    "MO"    "FA" 
  id.names <- names(id.newnames)
  
  id.var <- id.names[which(id.newnames == "ID")]
  tid.var <- paste0("t", id.var)

  data[, id.var] <- as.character(data[, id.var])
  
  ### Step 2: polygenic
  if(verbose) {
    cat(" * run polygenic with formula:\n")
    print(formula.polygenic)
  }  
  
  model <- solarPolygenic(formula.polygenic, data, ..., verbose = verbose)

  ### Step 3: loop over traits to construct `data2`
  traits <- traits.polygenic
  covlist <- unique(unlist(covlist.polygenic))

  covlist0 <- covlist
  covlist <- covlist[!(covlist %in% "1")]
  
  stopifnot(!any(covlist %in% id.names))

  # fix classes of some covariates: from `character` to `factor`
  ind.covlist <- laply(covlist, function(x) which(names(data) == x))
  classes.covlist <- laply(data[, ind.covlist], class)
  stopifnot(all(classes.covlist %in% c("numeric", "integer", "character", "factor")))

  ind.covlist.character <- ind.covlist[classes.covlist == "character"]
  for(i in ind.covlist.character) {
    data[, i] <- as.factor(data[, i])
  }

  classes.covlist <- laply(data[, covlist], class)
  stopifnot(all(classes.covlist %in% c("numeric", "integer", "factor")))

  out <- list()    

  for(i in 1:length(traits)) {
    t <- traits.polygenic[i]
    
    # response    
    tval <- data[, t]
    tval.scaled <- scale(tval, center = TRUE, scale = TRUE)
    tval.scaled <- tval.scaled + attributes(tval.scaled)[["scaled:center"]]
    
    # data frame `dt`
    dt <- data.frame(subset(data, select = id.names),
      tname = t, tnum = i - 1,
      trait = tval, strait = tval.scaled, 
      subset(data, select = covlist))
    
    covlist.t <- covlist.polygenic[[i]]
    for(f in covlist) {
      if(!(f %in% covlist.t)) {
        cl <- classes.covlist[f == covlist]
        if(cl %in% c("integer", "numeric")) {
          dt[, f] <- 0
        }
        else if(cl == "factor") {
          levels <- levels(dt[, f])
          n <- nrow(dt)
          dt[, f] <- factor(rep(1, n), levels = levels, labels = levels)
          
        } else {
          stop("erorr in classes.covlist")
        }
      }
    }      
    
    out[[t]] <- dt
  }

  data2 <- do.call(rbind, out)
  rownames(data2) <- NULL
  
  # process `na.rm`
  if(na.rm) {
    data2 <- na.omit(data2)
  }

  # intesect ids & fix `na.rm`
  tnum <- NULL # fix `no visible binding`
  ids.tnum <- llply(unique(data2$tnum), function(x) subset(data2, tnum == x, id.var, drop = TRUE))
  ids <- do.call(intersect, ids.tnum)
  
  ind <- (data2[, id.var] %in% ids)
  data2 <- data2[ind, ]
  
  # order of `ids`
  ids <- subset(data2, tnum == 0, id.var, drop = TRUE)
  
  # create column `tid.var`
  data2[, tid.var] <- paste0(data2[, id.var], "_", data2[, "tnum"])
  
  # `nobs`
  nobs.data <- nrow(data)  
  nobs <- nrow(data2)

  ### Step 4: `correlation` matrix
  if(verbose) cat(" * estimating kinship matrix...\n")

  kin2 <- solarKinship2(data)

  ids.kin2 <- rownames(kin2)
  stopifnot(all(ids %in% ids.kin2))

  ind <- laply(ids, function(x) which(ids.kin2 == x))
  kin2 <- kin2[ind, ind]

  # V = G + E
  if(verbose) cat(" * constructing `V = G + E` matrix...\n")

  rhog <- with(model$vcf, Estimate[varcomp == "rhog"])
  rhoe <- with(model$vcf, Estimate[varcomp == "rhoe"])
  
  h2r <- with(model$vcf, Estimate[grep("h2r\\(", varcomp)])
  e2 <- with(model$vcf, Estimate[grep("e2\\(", varcomp)])

  nids <- length(ids)

  G <- matrix(0.0, nrow = ntraits * nids, ncol = ntraits * nids)
  for(i in 1:ntraits) {
    ind <- (i - 1) * nids + 1:nids
    for(j in 1:ntraits) {
      jnd <- (j - 1) * nids + 1:nids
      
      if(i == j) {
        G[ind, jnd] <- h2r[i] * kin2
      } else {
        G[ind, jnd] <- rhog * sqrt(h2r[i]) * sqrt(h2r[j]) * kin2
      }
    }
  }

  E <- matrix(0.0, nrow = ntraits * nids, ncol = ntraits * nids)
  I <- diag(1.0, nids)
  for(i in 1:ntraits) {
    ind <- (i - 1) * nids + 1:nids
    for(j in 1:ntraits) {
      jnd <- (j - 1) * nids + 1:nids
      
      if(i == j) {
        E[ind, jnd] <- e2[i] * I
      } else {
        E[ind, jnd] <- rhoe * sqrt(e2[i]) * sqrt(e2[j]) * I
      }
    }
  }

  correlation <- G + E  
  rownames(correlation) <- data2[, tid.var]
  colnames(correlation) <- data2[, tid.var]  
  
  correlation <- Matrix::Matrix(correlation)
  
  ### Step 5: formula
  formula <- paste("trait ~ tnum",
   ifelse(length(covlist), paste("+", paste(covlist, collapse = " + ")), ""))
  formula <- as.formula(formula)
  
  ### output
  out <- list(call = mc, 
    formula.polygenic = formula.polygenic, 
    traits.polygenic = traits.polygenic, covlist.polygenic = covlist.polygenic,
    model = model, correlation = correlation,
    data = data2, id.var = id.var, tid.var = tid.var,
    nobs.data = nobs.data, nobs = nobs, ids = ids,
    traits = traits, covlist = covlist, covlist0 = covlist0,
    formula = formula
  )
  
  oldClass(out) <- c("solarMultivar")
  return(out)
}


#---------------------
# Class
#---------------------

#' S3 class solarMultivar.
#'
#' @name solarMultivarClass
#' @rdname solarMultivarClass
#'
#' @param x 
#'    An object of class \code{solarMultivar}.
#' @param ...
#'    Additional arguments.
#'
#' @exportClass solarMultivar

#' @rdname solarMultivarClass
#' @export
print.solarMultivar <- function(x, ...)
{
  cat("\nCall: ")
  print(x$call)
  
  cat("\nUnivariate polygenic models\n")
  cat(" * traits:", paste(x$traits.polygenic, collapse = " / "), "\n")
  cat(" * covlist:", paste(laply(x$covlist.polygenic, function(y) paste(y, collapse = ", ")), collapse = " / "), "\n")
  
  cat("\nMultivariate polygenic model\n")
  cat(" * traits:", paste(x$traits, collapse = ", "), "\n")
  cat(" * covlist:", paste(x$covlist, collapse = ", "), "\n")

  cat("\nData\n")
  cat(" * nobs.data ", x$nobs.data, ", nobs ", x$nobs, "\n", sep = "")
  cat(" * ids: ", paste(head(x$ids, 5), collapse = ", "), ", ... \n", sep = "")
  
  cat("\nCorrelation matrix (first 5x5 elements)\n")
  print(x$correlation[1:5, 1:5])
  
  cat("\nFormula\n")
  print(x$formula)  
}
