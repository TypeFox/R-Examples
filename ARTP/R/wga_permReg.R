# History Apr 05 2008  Initial coding
#         May 02 2008  Use cat instead of write
#                      Allow for Wald tests
#         May 07 2008  Add score tests
#         Jun 24 2008  Add out.file for observed results
#                      Rename interaction.vars to int.vars.
#                      Add main.vars
#         Jul 11 2008  Allow for different types of permutations
#         Aug 14 2008  Add new permutation method
#         Aug 29 2008  Remove sided option
#         Oct 21 2008  Make more efficient
#         Oct 22 2008  Use getModelData function
#         Jul 21 2009  Add option for output delimiter
#         Jul 24 2009  Use function writeVec
#                      Add option ndigits for formatted output
#         Oct 02 2009  Add option for single.marker.test function
#         Oct 02 2009  Make options consistent
#         Jun 02 2010  For test=4, add row to observed output for Fisher's test
#         Jun 14 2010  Update for no main effects
#         Jul 02 2010  Add observed beta, se for test = 4
#         Oct 08 2010  Rewrite for efficiency
#         Oct 13 2010  Add code for a large number of snps
#         Nov 02 2010  Add option to wait for snp temp data to be created (for running on a cluster)
#                      Do not gzip temp file.
#         Jan 11 2011  Add simplified function for the R package 
#                      Change the return flag from test 4, remove x from single.marker.test
#         Feb 01 2011  Add flag for continuous outcome in obs file of p-values
#         Mar 24 2011  Let min.count be applied to continuous outcomes
#         Feb 13 2012  Add option (mvn.k) to generate p-values from MVN distribution
#                      Add option to pass in the covariance matrix and save it

# Function to call glm for each SNP. 
permReg <- function(snp.list, pheno.list, op=NULL) {
  
  # pheno.list  
  #  response.var     Name of the response variable
  #                   No default
  #  main.vars        Character vector of variable names that
  #                   will be included as main effects.
  #                   The default is NULL
  #  int.vars         Character vector of variable names that
  #                   will interact with each SNP.
  #                   The default is NULL
  #              Not implemented yet
  ############################################################
  # op          A list with the folowing names:
  #  temp.list
  #  family           Family option for glm()
  #                   The default depends on response.var
  #  print            0 or 1
  #                   The default is 0
  #  observed         0 or 1 to compute the observed values
  #                   The default is 1
  #  test             1-4
  #                   1 = likelihood ratio
  #                   2 = Wald test
  #                   4 = Use the function single.marker.test. For this test,
  #                       interactions are not allowed (yet).
  #                   The default is 1.
  #  nperm            Number of permutations
  #                   The default is 0
  #  perm.method      1-3
  #                   1 is to permute the snps
  #                   2 is to generate new binary response vectors
  #                   3 is to generate new cont. response vector
  #                   The default is 1. 
  #  seed             If greater than 0, then the seed will be set
  #                   The default is -1
  #  out.delimiter    Delimiter used in output data sets
  #                   The default is a comma (",").
  #  ndigits          Number of digits in formatC() function.
  #                   The default is 6.
  #  test4.min.count      The default is 5
  #  test4.missing.rate   The default is 0.20
  #  obs.outfile      Output file for the observed values
  #                     The default is "obs.txt"
  #  perm.outfile     Output file for the permutation values
  #                     The default is "perm.txt
  #  out.file         Output file for the observed data only.
  #                   The file will contain 1 row for each snp.
  #                   The default is NULL.
  #  read.n           Integer for processing a large number of snps in chunks.
  #                   The default is 0
  #  read.n.file      Name of (temporary) file to hold the transformed genotype data.
  #                   The file will be deleted beforehand if waitForData = 0.
  #                   If the file name has .gz extention, then it will be gzipped.???
  #                   The default is NULL if read.n <= 0, else a temporary name will
  #                   be created from op$temp.list.
  #  read.n.miss      File for missing
  #  read.n.snps      File of snp names
  #  waitForData      0 or 1 to wait for the (temporary) genotype data set to be created.
  #                   The default is 0.
  #  mvn.k            Number of permutations to use to generate the correlation matrix
  #                   for generating p-values from a multi-variate normal distribution.
  #                   Set to 0 for generating p-values from the glm() function.
  #                   The default is 0
  #  mvn.method       "eigen", "svd", or "chol". The default is "chol"
  #  mvn.covFile      The .rda file containing the variance-covariance matrix to be used
  #                   when generating random vectors. This file must only contain 1 object.
  #                   The matrix must contain column names.
  #  mvn.covOut       Output file to save the covariance matrix using the save() function.
  #                   The default is NULL so that it will not be saved.

  # Local function
  fitModels <- function(vec, miss, response0) {

    vec <- as.numeric(vec)
    phenoData0[, snpcol] <- vec
    if (miss) {
      rows  <- !is.na(vec)
      nsubs <- sum(rows) 
      if (nsubs < 2) return(list(pvalue=NA, FisherFlag=NA, beta=NA, se=NA))
    } else {
      rows  <- allTRUE
      nsubs <- nsubjects 
    }

    # single.marker.test will take missing values into account
    if (test4) {
      if ((nsubjects-nsubs)/nsubjects >= miss.rate) return(list(pvalue=NA, FisherFlag=NA, beta=NA, se=NA))
      fit1 <- single.marker.test(response0[rows], phenoData0[rows,], 
               weights0[rows], offset0[rows], control, snpcol, min.count=min.count, 
               y.continuous=y.cont)
      if (print) print(fit1)
      return(list(pvalue=fit1[1], FisherFlag=fit1[2], beta=fit1[3], se=fit1[4]))
    }

    # Perform regression
    fit1 <- perm.glm(response0[rows], phenoData0[rows,], weights0[rows], offset0[rows],
              control, family) 

    # Check the return code
    if (fit1$converged) {
      if (print) print(summary(fit1))

      # Wald test 
      if (test2) {
        if (np1Flag) {
          tmp <- fit1$coefficients
          if (!is.na(tmp[snpcol])) {
            tmp <- summary(fit1)$coefficients 
            # Matrix can be reduced by other parameters with NA values
            ret <- list(pvalue=tmp[nrow(tmp), 4])
          } else {
            ret <- list(pvalue=NA)
          }
        } else {
          ret <- waldTest.main(fit1$coefficients, summary(fit1)$cov.scaled, snpcol)
        }
      } else {

        # Likelihood ratio test
        if (miss) {

          # Fit without the snp variables
          fit2 <- perm.glm(response0[rows], phenoData0[rows, -snpcol], weights0[rows],
                 offset0[rows], control, family) 
        } else {
          fit2 <- fit0
        }

        # Get the test and p-value
        ret <- likelihoodRatio(fit1, fit2)
      }
    } else {
      ret <- list(pvalue=NA)
    }

    ret

  } # END: fitModels

  # Function to load the correlation matrix
  loadCorr <- function(covfile) {

    ret <- loadFile(covfile, NULL)
    cnames <- colnames(ret)
    if (is.null(cnames)) stop("ERROR: covariance matrix in mvn.covFile must contain column names")
    rnames <- rownames(ret)
    if (is.null(rnames)) stop("ERROR: covariance matrix in mvn.covFile must contain row names")
    if (!all.equal(cnames, rnames)) stop("ERROR: row names and column names of mvn.covFile must be the same")

    temp <- snpNames %in% cnames
    if (sum(temp) != totalNsnps) {
      print("The following snps were not found in the covariance matrix:")
      print(snpNames[!temp])
      stop("ERROR: mvn.covFile does not contain all the required snps")
    }
    ret <- ret[snpNames, snpNames]

    ret

  } # END: loadCorr

  # Function to return the correlation matrix for the snps
  corrMatrix0 <- function(K) {

    # Matrix to hold p-values
    mat <- matrix(data=NA, nrow=K, ncol=totalNsnps)

    for (j in 1:K) {
      print(paste("perm ", j, sep=""))
      index <- 1
      if (chunkFlag) fid.list <- perm.openFiles(file.list)

      # Get the permutation
      temp      <- getPermutation(fit0, nsubjects, perm.method=perm.method)
      perm      <- temp$perm
      response0 <- temp$response

      while (1) {
        if (chunkFlag) {
          temp     <- perm.readData(fid.list, read.n)
          snpData  <- temp$snpData
          nsnps    <- length(snpData)
          rm(temp)
          gc()

          if (!nsnps) {
            perm.closeFiles(fid.list)
            break
          }
        }

        for (i in 1:nsnps) {
          # Permute the snp data
          temp <- getVecFromStr(snpData[i], delimiter=delimiter)
          temp <- temp[perm]
          temp <- fitModels(temp, missing[index], response0)

          out.pval[index] <- temp$pvalue
          index <- index + 1
        } 
        mat[j, ] <- out.pval 
        if (!chunkFlag) break

      } # END: while
    
    } # END: for (j in 1:K) 

    # Get the z-values
    mat <- qnorm(1-mat)

    # Get the correlation
    mat <- cor(mat)
    rownames(mat) <- snpNames
    colnames(mat) <- snpNames

    mat

  } # END: corrMatrix0

  # Function to return the correlation matrix for the snps
  corrMatrix <- function(K) {

    temp <- op[["mvn.covFile", exact=TRUE]]
    if (!is.null(temp)) {
      mat <- loadCorr(temp)
    } else {
      mat <- corrMatrix0(K)
    }

    temp <- op[["mvn.covOut", exact=TRUE]]
    if (!is.null(temp)) save(mat, file=temp)

    # Check for missing values
    temp <- !is.finite(mat)
    if (any(temp)) {
      nc   <- ncol(mat)
      temp <- rowSums(temp)
      temp <- temp == nc
      rows <- (1:nc)[temp]

      # Remove these rows/columns
      mat <- mat[-rows, -rows]
    } else {
      rows <- NULL
    }

    return(list(corr=mat, miss=rows)) 

  } # END: corrMatrix

  # Function to write out permutation p-values for mvn method
  mvn_pval <- function(fid, nperm) {

    nc   <- ncol(mvn.corr)
    mu   <- rep(0.0, nc)
    n    <- nperm 
    while (1) {
      # Try allocating memory for as many permutations as possible
      mat <- try(matrix(data=NA, nrow=n, ncol=totalNsnps), silent=TRUE)
      if (checkTryError(mat, conv=0)) {
        # Not enough memory available, free and set n smaller
        rm(mat)
        gc()
        n <- ceiling(n/2)
      } else {
        rm(mat)
        gc()
        break
      }
    }
    n  <- ceiling(n/6)  # For temp space
    ii <- 0
    nn <- n
    while (1) {
      mat <- myrmvnorm(nn, mean=mu, sigma=mvn.corr, method=mvn.method)
      mat <- 1 - pnorm(mat)
      mat <- formatC(mat, digits=ndigits, format="e")
      vec <- character(nn)
      for (j in 1:nn) vec[j] <- paste(mat[j,], collapse=out.sep, sep="")
      ii  <- ii + nn
      write(vec, file=fid.pvalue, ncolumns=1)
      flush(fid.pvalue)
      print(paste(ii, " permutations have been output", sep=""))
      if (ii >= nperm) break
      nn <- min(nn, nperm-ii)
      rm(vec)
      gc()
    }
    close(fid.pvalue)

  } # END: mvn_pval

  # Check the lists
  snp.list <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)

  if (is.null(op)) op <- list()
  op <- default.list(op, c("family", "y.cont"), list("ERROR", "ERROR"), error=c(1,1)) 
  op <- default.list(op, c("print", "observed", "nperm", "seed", 
          "test", "sided", "perm.method", "out.delimiter",
          "ndigits", "test4.min.count", "test4.missing.rate", "obs.outfile",
           "perm.outfile", "read.n", "waitForData", "mvn.k", "mvn.method"),
        list(0, 1, 0, -1, 1, 2, 1, ",", 6, 5, 0.20, "obs.txt", "perm.txt", 0, 0, 0, "chol"))
  temp.list <- op[["temp.list", exact=TRUE]]
  temp.list <- check.temp.list(temp.list)
  mvn.k     <- 0
  if (mvn.k) {
    mvn.covFile <- op[["mvn.covFile", exact=TRUE]]
    if (!is.null(mvn.covFile)) {
      if (!file.exists(mvn.covFile)) {
        temp <- paste("ERROR: The file ", mvn.covFile, " was not found.", sep="")
        stop(temp)
      }
    }
  }

  # Check for a response var
  pheno.list <- default.list(pheno.list, c("response.var"),
                list("ERROR"), error=c(1))

  # Keep only the variables we need
  temp <- c(pheno.list$response.var, pheno.list$id.var,
            pheno.list$main.vars, pheno.list$int.vars)
  pheno.list$keep.vars   <- unique(temp)
  pheno.list$remove.vars <- NULL
  pheno.list$make.dummy  <- 0
  pheno.list$remove.miss <- 1

  #  No streaming yet for permutation tests
  snp.list$stream <- 0

  # Determine if processing in chunks will be done
  read.n    <- op$read.n
  chunkFlag <- (read.n > 0)

  # Get the data vector of snps
  temp <- perm.getData(snp.list, pheno.list, temp.list, op=op) 

  if (chunkFlag) file.list  <- temp$file.list
  temp       <- temp$data
  snpData    <- temp$data 
  missing    <- temp$missing
  snpNames   <- temp$snpNames
  nsnps      <- length(snpNames)
  totalNsnps <- nsnps
  nsubjects  <- temp$nsubjects
  delimiter  <- snp.list$delimiter
  newVar     <- "SNP"

  # Get the phenotype data
  temp       <- temp$phenoData.list
  phenoData0 <- temp$data
  if (!is.null(temp$orig.id)) phenoData0[, temp$orig.id] <- NULL
  rm(temp, snp.list)

  # Remove the response variable from phenoData
  response0 <- phenoData0[, pheno.list$response.var]
  phenoData0[, pheno.list$response.var] <- NULL
  phenoData0[, pheno.list$id.var] <- NULL
  flag.main <- !is.null(pheno.list$main.vars)

  # Detemine the permutation method
  perm.method <- op[["perm.method", exact=TRUE]]

  # Check the response
  family <- op[["family", exact=TRUE]]
  y.cont <- op[["y.cont", exact=TRUE]]
  if (y.cont) {
    if (perm.method > 1) perm.method <- 3
  } else {
    if (perm.method > 1) perm.method <- 2
  }

  # Check for interactions
  temp <- pheno.list[["int.vars", exact=TRUE]]
  if (!is.null(temp)) {
    flag.inter   <- 1
    pheno0.inter <- removeOrKeepCols(phenoData0, temp, which=1)
    pheno0.inter <- createDummy(pheno0.inter)$data
    colnames(pheno0.inter) <- paste("v", 1:ncol(pheno0.inter), sep="")
    np1Flag <- 0
  } else {
    pheno0.inter <- NULL
    flag.inter   <- 0
    np1Flag      <- 1
  }

  rm(pheno.list)
  temp <- gc(verbose=FALSE)

  # Initialize
  print     <- op$print
  obs.pval  <- NULL
  test      <- op$test
  nperm     <- op$nperm
  test4     <- test == 4
  test2     <- test == 2

  # Make sure factor variables are converted to dummy variables
  if (flag.main) {
    phenoData0 <- createDummy(phenoData0)$data

    # Add intercept and column for SNP
    phenoData0 <- cbind(1, phenoData0, 0)
 
    # Convert to matrices
    phenoData0 <- as.matrix(phenoData0)

  } else {
    phenoData0 <- cbind(rep.int(1, nsubjects), rep.int(0, nsubjects))
  }
  snpcol <- ncol(phenoData0)
  colnames(phenoData0) <- paste("x", 1:snpcol, sep="")
  
  if (flag.inter) pheno0.inter <- as.matrix(pheno0.inter)


  # Define the output vectors
  out.pval <- rep(NA, times=totalNsnps)

  # Fit the model on all the phenotype data without any snp
  fit0 <- glm(response0 ~ phenoData0[, -snpcol] - 1, family=family, model=FALSE, y=TRUE, x=TRUE)
  if (!fit0$converged) stop("ERROR: model failed")
  print(summary(fit0))

  # See if out.file will be created
  out.file.flag <- (!is.null(op[["out.file", exact=TRUE]]) && op$observed)
  out.sep       <- op$out.delimiter

  ndigits <- op$ndigits
  if (test4) {
    min.count  <- op$test4.min.count
    miss.rate  <- op$test4.missing.rate
    out.fisher <- rep.int(0, times=totalNsnps) 
    out.beta   <- rep(NA, times=totalNsnps)
    out.se     <- rep(NA, times=totalNsnps)
  }

  weights0 <- rep.int(1, times=nsubjects)
  offset0  <- rep.int(0, times=nsubjects)
  control  <- glm.control()
  allTRUE  <- rep(TRUE, times=nsubjects)

  # Get the observed p-values
  if (op$observed) {
    index <- 1

    if (chunkFlag) fid.list <- perm.openFiles(file.list)
 
    while (1) {
      if (chunkFlag) {
        temp     <- perm.readData(fid.list, read.n)
        snpData  <- temp$snpData
        nsnps    <- length(snpData)
        rm(temp)
        gc()

        if (!nsnps) {
          perm.closeFiles(fid.list)
          break
        }
      }

      for (i in 1:nsnps) {
        temp <- getVecFromStr(snpData[i], delimiter=delimiter)
        temp <- fitModels(temp, missing[index], response0)

        out.pval[index] <- temp$pvalue
        if (test4)  {
          out.fisher[index] <- temp$FisherFlag
          out.beta[index]   <- temp$beta
          out.se[index]     <- temp$se
        } 
        index <- index + 1
      }
      if (!chunkFlag) break

    } # END: while 

    if (out.file.flag) {
      cnames <- "SNP"
      temp   <- snpNames
      cnames <- c(cnames, "Pvalue")
      temp   <- cbind(temp, out.pval) 
      
      if (test4) {
        cnames <- c(cnames, "Pvalue.flag", "Beta", "SE")
        temp   <- cbind(temp, out.fisher, out.beta, out.se)
      }

      colnames(temp) <- cnames
      write.table(temp, file=op$out.file, row.names=FALSE, sep=out.sep, quote=FALSE, col.names=TRUE)
    }

    # Check for all missing
    if (all(is.na(out.pval))) stop("ERROR: all observed p-values are missing")

    # Write out the results
    fid.pvalue <- writeVec(formatC(out.pval, digits=ndigits, format="e"),
                  op$obs.outfile, colnames=snpNames, sep=out.sep, isFID=0, close=0)
    if (test4) {
      writeVec(out.fisher, fid.pvalue, sep=out.sep, isFID=1, close=0)
      writeVec(out.beta, fid.pvalue, sep=out.sep, isFID=1, close=0)
      writeVec(out.se, fid.pvalue, sep=out.sep, isFID=1, close=0)
      write(y.cont, file=fid.pvalue, ncolumns=1)
    }
    close(fid.pvalue)

  } # END: if (op$observed)

  if (test4) {
    rm(out.fisher, out.beta, out.se)
    gc()
  }

  # Stop if no permutations are to be done
  if (!op$nperm) {
    #if ((chunkFlag) && (temp.list$delete)) perm.deleteFiles(file.list)
    return(0)
  }

  # Set the seed
  if (op$seed > 0) set.seed(op$seed)

  mvn.corr <- NULL
  if (mvn.k) {
    mvn.method   <- op$mvn.method
    temp         <- corrMatrix(mvn.k)
    mvn.corr     <- temp$corr
    mvn.miss     <- temp$miss
    mvn.missFlag <- !is.null(mvn.miss) 
    rm(temp)
    gc()
  }

  # Open the output files
  fid.pvalue <- writeVec(snpNames, op$perm.outfile,
                                type=3, close=0, sep=out.sep)
  rm(snpNames)
  gc()

  if (mvn.k) {
    rm(snpData, missing, phenoData0, response0, weights0, offset0, control, allTRUE, out.pval)
    gc()
    mvn_pval(fid.pvalue, nperm)
    return(0)
  }

  for (j in 1:nperm) {
    index <- 1
    if (chunkFlag) fid.list <- perm.openFiles(file.list)

    # Get the permutation
    temp      <- getPermutation(fit0, nsubjects, perm.method=perm.method)
    perm      <- temp$perm
    response0 <- temp$response

    while (1) {
      if (chunkFlag) {
        temp     <- perm.readData(fid.list, read.n)
        snpData  <- temp$snpData
        nsnps    <- length(snpData)
        rm(temp)
        gc()

        if (!nsnps) {
          perm.closeFiles(fid.list)
          break
        }
      }

      for (i in 1:nsnps) {
        # Permute the snp data
        temp <- getVecFromStr(snpData[i], delimiter=delimiter)
        temp <- temp[perm]
        temp <- fitModels(temp, missing[index], response0)

        out.pval[index] <- temp$pvalue
        index <- index + 1
      } 
      if (!chunkFlag) break

    } # END: while

    # Output the results
    writeVec(formatC(out.pval, digits=ndigits, format="e"), fid.pvalue, sep=out.sep, isFID=1, close=0)
    
  } # END: for (j in 1:op$nperm) 
  
  # Close file
  close(fid.pvalue)
   
  #if ((chunkFlag) && (temp.list$delete)) perm.deleteFiles(file.list)
 
  NULL

} # END: permReg

perm.glm <- function(y, x, wgt, off, control, fam) {

  fit <- glm.fit(x, y, weights=wgt, offset=off, family=fam, control=control) 

  class(fit) <- c("glm", "lm")

  fit

}

single.marker.test <- function(y, covariates, weights, offset, control, snpcol,
                     min.count=5, y.continuous=FALSE)
{
    # y is the outcome
    # Covariates is the design matrix which includes the intercept and snp
    # The snp must be the last column of covariates and must be coded as 0-1-2
    
    getReturnVec <- function(fit, flag) {
      if (!is.na(fit$coefficients[snpcol])) {
        scoef <- summary(fit)$coefficients 
        coef  <- scoef[nrow(scoef), ]
        ret   <- c(coef[4], flag, coef[1], coef[2])
      } else {
        ret <- c(NA, flag, NA, NA)
      }

      ret

    } # END: getReturnVec

    Flag <- 0  
 
    # y is continuous
    if (y.continuous)
    { 
      # Check that at least 2 cell counts are at least min.count
      x <- table(covariates[, snpcol])
      temp <- x >= min.count
      if (sum(temp) < 2) return(c(NA, Flag, NA, NA))
      temp <- lm(y~covariates-1)
      return(getReturnVec(temp, Flag))
    }
    else
    {
        x <- covariates[, snpcol]
        yeq1 <- y == 1
        yeq0 <- !yeq1
        xeqi <- x == 0
        r0   <- sum(xeqi & yeq1)
        s0   <- sum(xeqi & yeq0)
        xeqi <- x == 1
        r1   <- sum(xeqi & yeq1)
        s1   <- sum(xeqi & yeq0)
        xeqi <- x == 2
        r2 <- sum(xeqi & yeq1)
        s2 <- sum(xeqi & yeq0)
       
            if (min(cal.expect(matrix(c(r0,r1,r2,s0,s1,s2),nrow=2,byrow=TRUE)))>=min.count)
            {             
              temp <- perm.glm(y, covariates, weights, offset, control, binomial())
              return(getReturnVec(temp, Flag))
            }
            else
            {
                if (s0 >= s2) 
                { 
                    Flag <- 1
                    H <- matrix(c(r0,r1+r2,s0,s1+s2),ncol=2, byrow=TRUE) 
                    if (min(cal.expect(H))<min.count)
                    {
                        p.value <- fisher.test(H)$p.value
                        Flag <- -Flag
                        return(c(p.value, Flag, NA, NA))
                    }  
                    else
                    {   
                        x[x==2] <- 1
                        covariates[, snpcol] <- x
                        
                        temp <- perm.glm(y, covariates, weights, offset, control, binomial())
                        return(getReturnVec(temp, Flag))
                    } 
                }
                else 
                { 
                    Flag <- 2
                    H <- matrix(c(r0+r1,r2,s0+s1,s2),ncol=2,byrow=TRUE) 
                    if (min(cal.expect(H))<min.count)
                    {
                        p.value <- fisher.test(H)$p.value
                        Flag    <- -Flag
                        return(c(p.value, Flag, NA, NA))
                    }  
                    else
                    {   
                        x[x==1] <- 0
                        x[x==2] <- 1
                        covariates[, snpcol] <- x
                        
                        temp <- perm.glm(y, covariates, weights, offset, control, binomial())
                        return(getReturnVec(temp, Flag))
                    }
                }  
            }     
       
    }     
}

##
cal.expect <- function(H)
{
    a <- apply(H, 1, sum)
    b <- apply(H, 2, sum)
    n <- sum(H)
    
    kronecker(a,b)/n
}


# Function to transform and write out the snp data
perm.getData <- function(snp.list, pheno.list, temp.list, op=NULL) {

  read.n     <- op$read.n
  writeData  <- (read.n > 0)
  wait       <- op$waitForData
  tmpFlag    <- 0
  data.file0 <- op[["read.n.file", exact=TRUE]]
  miss.file  <- op[["read.n.miss", exact=TRUE]]
  snp.file   <- op[["read.n.snps", exact=TRUE]]
  if (is.null(data.file0)) {
    data.file <- getTempfile(temp.list$dir, paste("snp_", temp.list$id, "_", sep=""), 
                             ext=".ldat")
    wait      <- 0
    ret.f     <- data.file
  } else {
    # Store data in another temp file first
    data.file <- getTempfile(temp.list$dir, paste("tmp_", temp.list$id, "_", sep=""), 
                             ext=".ldat")
    tmpFlag   <- 1
    ret.f     <- data.file0
  }
  file.list <- list(file=ret.f, file.type=2, delimiter=snp.list$out.delimiter)

  tlist  <- list(include.row1=0, include.snps=0, return.type=1,
               missing=1, snpNames=1, orderByPheno=1, return.pheno=1)
  start  <- snp.list$start.vec
  stop0  <- snp.list$stop.vec
  if (stop0 < 1) stop0 <- Inf
  if ((read.n > 0) && (writeData)) {
    if (start + read.n >= stop0) read.n <- -1
  }

  if (wait) {
    # Just get the phenotype data
    snp.list$start.vec <- 1
    snp.list$stop.vec  <- 3
    writeData          <- 0
    read.n             <- -1

    # Wait for the data
    while (1) {
      if (file.exists(data.file0)) break

      # Check every minute
      Sys.sleep(60)
    }
  }

  if ((!writeData) || (read.n < 1)) {
    temp  <- try(getData.1(snp.list, pheno.list, temp.list, op=tlist),
               silent=TRUE)
    if (checkTryError(temp, conv=0)) {
      print(temp)
      stop("ERROR loading data")
    }

    if (wait) {
      # Get the missing vector and snp names
      temp$missing  <- scan(miss.file, what=integer(0), sep="\n")
      temp$snpNames <- scan(snp.file, what="character", sep="\n") 
    }

    if (!writeData) return(list(file.list=list(file.list), data=temp))

    write(temp$data, file=data.file, ncolumns=1)
    temp$data <- NULL
    gc()
    file.list <- list(file=data.file, file.type=2, delimiter=snp.list$out.delimiter)
    return(list(file.list=list(file.list), data=temp))
  } 

  stop      <- start + read.n
  missing   <- NULL
  snpNames  <- NULL
  
  fid <- file(data.file, "w")
  while (1) {
    snp.list$start.vec <- start
    snp.list$stop.vec  <- stop
    
    temp  <- try(getData.1(snp.list, pheno.list, temp.list, op=tlist),
               silent=TRUE)
    if (checkTryError(temp, conv=0)) {
      close(fid)
      print(temp)
      stop("ERROR loading data")
    }
    if (!length(temp$data)) {
      temp0 <- temp
      break
    }

    write(temp$data, file=fid, ncolumns=1)
    temp$data <- NULL
    gc()
    missing  <- c(missing, temp$missing)
    snpNames <- c(snpNames, temp$snpNames)

    start <- stop + 1
    stop  <- start + read.n
    if (start > stop0) {
      temp0 <- temp
      break
    }
    if (stop > stop0) stop <- stop0
  }
  close(fid)

  missing        <- as.integer(missing)
  temp0$missing  <- missing
  temp0$snpNames <- snpNames

  # Write out missing and snp names
  write(missing, file=miss.file, ncolumns=1)
  write(snpNames, file=snp.file, ncolumns=1)

  if (tmpFlag) {
    # Rename the file
    if (!file.rename(data.file, data.file0)) stop("ERROR: renaming file")
  }

  file.list <- list(file=ret.f, file.type=2, delimiter=snp.list$out.delimiter)
  list(file.list=list(file.list), data=temp0)

} # END: perm.getData

# Function to read data from files
perm.readData <- function(fid.list, read.n) {

  # Order should be data, names, missing
  snpData  <- scan(fid.list[[1]], what="character", sep="\n", nlines=read.n)
  #missing  <- scan(fid.list[[2]], what=integer(0), sep="\n", nlines=read.n)
  #snpNames <- scan(fid.list[[3]], what="character", sep="\n", nlines=read.n)

  list(snpData=snpData)

} # END: perm.readData

perm.openFiles <- function(file.list) {

  ret <- list()
  type <- c(2, 3, 3)
  for (i in 1:length(file.list)) {
    temp <- file.list[[i]]
    temp$file.type <- type[i]
    temp$open <- "r"
    ret[[i]] <- getFID(temp$file, temp)
  }
  ret

}
perm.closeFiles <- function(fid.list) {
 
  for (i in 1:length(fid.list)) close(fid.list[[i]])
  
}
perm.deleteFiles <- function(file.list) {
 
  for (i in 1:length(file.list)) file.remove(file.list[[i]]$file)
  
}

