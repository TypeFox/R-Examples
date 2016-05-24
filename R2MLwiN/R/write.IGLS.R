#' Writes MLwiN macros to fit models using the iterative generalized least
#' squares (IGLS) algorithm.
#' 
#' write.IGLS is an internal function which creates an MLwiN macro file to
#' fit a multilevel model using IGLS. 
#' 
#' @param indata A data.frame object containing the data to be modelled.
#' @param dtafile The name of the temporary file used to send the data to MLwiN, which
#' is in Stata format (i.e. with extension .dta).
#' @param oldsyntax Specified as \code{FALSE} if new syntax has been used in
#' \code{Formula} object, else specified as \code{TRUE} (to enable
#' back-compatibility).
#' @param resp A character string (vector) of the response variable name(s).
#' @param levID A character string (vector) of the specified level ID name(s). The
#' ID(s) should be sorted in the descending order of levels (e.g.
#' \code{levID = c('level2', 'level1')} where \code{'level2'} is the higher
#' level).
#' @param expl A character string (vector) of explanatory (predictor)
#' variable name(s).
#' @param rp A character string (vector) of random part of random variable name(s).
#' @param D A character string/vector specifying the type of distribution to be modelled, which
#' can include \code{'Normal'} (the default), \code{'Binomial'}, \code{'Poisson'},
#' \code{'Negbinom'}, \code{'Unordered Multinomial'}, \code{'Ordered Multinomial'},
#' \code{'Multivariate Normal'}, or \code{'Mixed'}. In the case of the latter, 
#' \code{'Mixed'} precedes the response types which also need to be listed in 
#' \code{D}, e.g. \code{c('Mixed', 'Normal', 'Binomial')}; these need to be
#' be listed in the same order to which they are referred to in the
#' \code{Formula} object (see \code{\link{runMLwiN}}, \code{\link{Formula.translate}},
#' \code{\link{Formula.translate.compat}}. \code{'Mixed'} combinations can consist of
#' \code{'Normal'} and \code{'Binomial'} or \code{'Normal'} and \code{'Poisson'}.
#' @param nonlinear A character vector specifying linearisation method for discrete
#' response models (see Chapter 9 of Rasbash et al 2012, and Goldstein 2011).
#' \code{N = 0} specifies marginal quasi-likelihood
#' linearization (MQL), whilst \code{N = 1} specifies penalised quasi-
#' likelihood linearization (PQL); \code{M = 1} specifies first order
#' approximation, whilst \code{M = 2} specifies second order approximation.
#' \code{nonlinear = c(N = 0, M = 1)} by default. First order marginal
#' quasi-likelihood (MQL1) only option for single-level discrete response
#' models.
#' @param categ Specifies categorical variable(s) as a matrix. Each column
#' corresponds to a categorical variable; the first row specifies the name(s)
#' of variable(s); the second row specifies the name(s) of reference group(s),
#' \code{NA}(s) if no reference group; the third row states the number of
#' categories for each variable. Supports back-compatibility with \pkg{R2MLwiN}
#' versions <0.8-0.
#' @param notation Specifies the model subscript notation to be used in the
#' MLwiN equations window. \code{'class'} means no multiple subscripts, whereas
#' \code{'level'} has multiple subscripts.
#' @param nonfp Removes the fixed part of random variable(s). \code{NA} if no
#' variable to be removed.
#' @param clre A matrix used to define which elements of the random effects matrix
#' to remove (i.e. hold constant at zero). Removes
#' from the random part at level <first row> the covariance matrix element(s)
#' defined by the pair(s) of rows <second row> <third row>. Each column
#' corresponds to a removed entry of the covariance matrix. See e.g. \code{demo(UserGuide07)}
#' for an example.
#' @param Meth Specifies which maximum likelihood estimation method to be used.
#' If \code{Meth = 0} estimation method is set to RIGLS. If \code{Meth = 1}
#' estimation method is set to IGLS (the default setting).
#' @param extra If \code{TRUE}, extra binomial, extra negative binomial,
#' extra Poisson or extra multinomial distributions assumed, else \code{FALSE}.
#' @param reset A vector of \code{length(levID)} specifying the action to be
#' taken, at each level, if a variance parameter is estimated at a particular
#' iteration to be negative during estimation. Values specified in
#' ascending order of level hierarchy: if \code{0} a negative variance
#' estimate is reset to zero and so are any associated covariances; if \code{1}
#' a negative variance estimate is reset to zero but not the associated
#' covariances; if \code{2} no resetting takes place. E.g. \code{reset = c(0, 1)}
#' to assign value \code{0} to level 1 and value \code{1} to level 2 of
#' two-level model.
#' @param rcon Matrix specifying constraints on the random parameters as
#' specified in \code{random.ui} and \code{random.ci} in the \code{constraints}
#' option within \code{estoptions} (see \code{\link{runMLwiN}}). \code{random.ci}
#' is appended to the bottom row of \code{random.ui}.
#' @param fcon Matrix specifying constraints on the fixed coefficients as
#' specified in \code{fixed.ui} and \code{fixed.ci} in the \code{constraints}
#' option within \code{estoptions} (see \code{\link{runMLwiN}}). \code{fixed.ci}
#' is appended to the bottom row of \code{fixed.ui}.
#' @param maxiter Numeric value specifying the maximum number of iterations, from
#' the start, before estimation halts.
#' @param convtol Numeric value specifying the convergence criterion, as
#' specified in the \code{tol} option within \code{estoptions} (see
#' \code{\link{runMLwiN}}). If value of \code{convtol} is m, estimation will be
#' deemed to have converged when the relative change in the estimate for all
#' parameters from one iteration to the next is less than 10(-m). Defaults to
#' value of \code{2} for m if not otherwise specified.
#' @param mem.init If calling write.IGLS directly, if wish to use defaults, value needs to be
#' specified as \code{'default'}, else specify a vector of length 5 corresponding
#' to the following order: number of levels; worksheet size in thousands of cells;
#' the number of columns; the number of explanatory variables; the number of group
#' labels. 
#' @param optimat This option instructs MLwiN to limit the maximum matrix size
#' that can be allocated by the (R)IGLS algorithm. Specify \code{optimat = TRUE}
#' if MLwiN gives the following error message 'Overflow allocating smatrix'.
#' This error message arises if one more higher-level units is extremely large
#' (contains more than 800 lower-level units). In this situation \pkg{R2MLwiN}'s
#' default behaviour is to instruct MLwiN to allocate a larger matrix size to
#' the (R)IGLS algorithm than is currently possible. Specifying
#' \code{optimat = TRUE} caps the maximum matrix size at 800 lower-level units,
#' circumventing the MLwiN error message, and allowing most MLwiN
#' functionality.
#' @param weighting A list of two items, one of which is a list called \code{weightvar}
#' the length of which corresponds
#' to the number of levels in the model, in descending order from highest level first.
#' The other is an option \code{standardised} which is \code{TRUE} or \code{FALSE}.
#' @param fpsandwich Specifies standard error type for fixed parameters. If
#' \code{fpsandwich = TRUE}, robust or `sandwich' standard errors based on raw
#' residuals are used, if \code{fpsandwich = FALSE} (default) then standard,
#' uncorrected, IGLS or RIGLS computation used.
#' @param rpsandwich  Specifies standard error type for random parameters. If
#' \code{rpsandwich = TRUE}, robust or `sandwich' standard errors based on raw
#' residuals are used, if \code{rpsandwich = FALSE} (default) then standard,
#' uncorrected, IGLS or RIGLS `plug in' estimates used.
#' @param macrofile A file name where the MLwiN macro file will be saved. The
#' default location is in the temporary folder.
#' @param IGLSfile A file name where the parameter estimates will be saved. The
#' default location is in the temporary folder.
#' @param resifile A file name where the residuals will be saved. The default
#' location is in the temporary folder.
#' @param resi.store A logical value to indicate if the residuals are to be
#' stored (\code{TRUE}) or not (\code{FALSE}).
#' @param resioptions A string vector to specify the various residual options.
#' The \code{'variances'} option calculates the posterior variances instead of
#' the posterior standard errors; the \code{'standardised'}, \code{'leverage'},
#' \code{'influence'} and \code{'deletion'} options calculate standardised,
#' leverage, influence and deletion residuals respectively; the
#' \code{'sampling'} option calculates the sampling variance covariance matrix
#' for the residuals; the \code{'norecode'} option prevents residuals with
#' values exceedingly close or equal to zero from being recoded to missing; the
#' reflate option returns unshrunken residuals. Note that the default option is
#' \code{resioptions = c('variance')}; \code{'variance'} cannot be used together
#' with the other options to calculate standardised, leverage, influence and
#' deletion residuals.
#' @param debugmode A logical value determining whether MLwiN is run in the
#' background or not. The default value is \code{FALSE}: i.e. MLwiN is run in
#' the background. If \code{TRUE} the MLwiN GUI is opened, and then pauses after the model
#' has been set-up, allowing user to check starting values; pressing 'Resume macro'
#' will then fit the model. Once fit, pressing 'Resume macro' once more will save
#' the outputs to the \code{workdir} ready to be read by \pkg{R2MLwiN}. Users can
#' instead opt to 'Abort macro' in which case the outputs are not saved to the
#' \code{workdir}. This option currently
#' works for 32 bit version of MLwiN only (automatically switches unless
#' \code{MLwiNPath} or \code{options(MLwiNPath)}
#' has been set directly to the executable).
#' @param startval A list of numeric vectors specifying the starting values.
#' \code{FP.b} corresponds to the estimates for the fixed
#' part; \code{FP.v} specifies the variance/covariance estimates for the fixed
#' part; \code{RP.b} specifies the variance estimates for the random part;
#' \code{RP.v} corresponds to the variance/covariance matrix of the variance
#' estimates for the random part.
#' @param namemap A mapping of column names to DTA friendly shorter names
#' @param saveworksheet A file name used to store the MLwiN worksheet after the
#' model has been estimated.
#' 
#' @return Outputs a modified version of namemap containing newly generated
#' short names.
#' 
#' @references
#' Goldstein, H. (2011) Multilevel Statistical Models. 4th Edition. London: John Wiley and Sons.
#' 
#' Rasbash, J., Steele, F., Browne, W.J. and Goldstein, H. (2012)
#' A User's Guide to MLwiN Version 2.26. Centre for Multilevel Modelling,
#' University of Bristol.
#' 
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#' 
#' @seealso \code{\link{write.MCMC}}
#' 
write.IGLS <- function(indata, dtafile, oldsyntax = FALSE, resp, levID, expl, rp, D = "Normal", nonlinear = c(0, 1), categ = NULL,
                         notation = NULL, nonfp = NA, clre = NULL, Meth = 1, extra = FALSE, reset, rcon = NULL, fcon = NULL, 
                         maxiter = 20, convtol = 2, mem.init = "default", optimat = FALSE, weighting = NULL, fpsandwich = FALSE, rpsandwich = FALSE, 
                         macrofile, IGLSfile, resifile, resi.store = FALSE, resioptions, debugmode = FALSE, startval = NULL,
                         namemap = sapply(colnames(indata), as.character), saveworksheet = NULL) {
  
  shortname <- function(...) {
    name <- paste0(...)
    if (!name %in% names(namemap)) {
      sname <- paste0("v", digest(name, algo = "xxhash64", serialize = FALSE))
      names(sname) <- name
      namemap <<- c(namemap, sname)
    }
    return(namemap[[name]])
  }
  
  nlev <- length(levID)
  
  nrp <- length(rp)
  if (nrp > 0) {
    rp.names <- names(rp)
    if (D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed") {
      for (i in 1:nrp) {
        temp <- rp.names[i]
        rp.names[i] <- paste("rp", as.integer(sub("rp", "", temp)) + 1, sep = "")
      }
    }
    names(rp) <- rp.names
  }
  num.expl.init <- function(p, nonfp, categ) {
    num_vars <- 0
    if (is.na(nonfp[1])) {
      if (is.null(categ) || sum(p == categ["var", ]) == 0) {
        num_vars <- num_vars + 1
      } else {
        if (is.na(categ["ref", which(p == categ["var", ])])) {
          num_vars <- num_vars + as.numeric(categ["ncateg", which(p == categ["var", ])])
        } else {
          num_vars <- num_vars + as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1
        }
      }
    } else {
      if (sum(p == nonfp) == 0) {
        if (is.null(categ) || sum(p == categ["var", ]) == 0) {
          num_vars <- num_vars + 1
        } else {
          if (is.na(categ["ref", which(p == categ["var", ])])) {
            num_vars <- num_vars + as.numeric(categ["ncateg", which(p == categ["var", ])])
          } else {
            num_vars <- num_vars + as.numeric(categ["ncateg", which(p == categ["var", ])]) - 1
          }
        }
      }
    }
    num_vars
  }
  
  if (is.list(expl)) {
    sep.coeff <- expl$sep.coeff
    if (is.list(nonfp)) {
      nonfp.sep <- nonfp$nonfp.sep
      nonfp.common <- nonfp$nonfp.common
    }
    if (!is.na(sep.coeff[1])) {
      num_vars <- sum(sapply(sep.coeff, function(x) num.expl.init(x, nonfp.sep, categ)))
    } else {
      num_vars <- 0
    }
    
    if (D[1] == "Multinomial") {
      nresp <- length(levels(indata[, resp])) - 1
      ## it is true if adding each expl variables separately
      num_vars <- num_vars * nresp
    }
    if (D[1] == "Multivariate Normal") {
      nresp <- length(resp)
      num_vars <- num_vars * nresp
    }
    
    common.coeff <- expl$common.coeff
    num_vars <- num_vars + sum(sapply(common.coeff, function(x) num.expl.init(x, nonfp$nonfp.common, categ)))
    common.coeff.id <- expl$common.coeff.id
    
  } else {
    num_vars <- sum(sapply(expl, function(x) num.expl.init(x, nonfp, categ)))
    if (D[1] == "Multinomial") {
      nresp <- length(levels(indata[, resp])) - 1
      num_vars <- num_vars * nresp
    }
    if (D[1] == "Multivariate Normal") {
      nresp <- length(resp)
      num_vars <- num_vars * nresp
    }
  }
  
  if (nrp > 0) {
    for (ii in 1:nrp) {
      rpx <- rp[[rp.names[ii]]]
      nrpx <- length(rpx)
      if (nrpx == 1) 
        num_vars <- num_vars + 1
      if (nrpx >= 2) 
        num_vars <- num_vars + nrpx * (nrpx - 1)/2 + nrpx
    }
  }
  
  ## Write into macro file
  wrt <- function(x) write(x, macrofile, append = TRUE)
  
  cat(file = macrofile)
  wrt("ECHO    0")
  wrt("NOTE    *****************************************************************")
  wrt("NOTE      MLwiN macro created by rtomlwin command")
  wrt("NOTE    *****************************************************************")
  wrt("")
  
  wrt("NOTE    Initialise MLwiN storage")
  if (mem.init[1] == "default") {
    wrt(paste("INIT    ", nlev + 1, " 6000 2500 ", num_vars + 10, " 30", sep = ""))
  } else {
    wrt(paste("INIT    ", mem.init[1], mem.init[2], mem.init[3], mem.init[4], mem.init[5]))
  }
  
  wrt("NOTE     Limit the maximum matrix size")
  if (optimat) {
    wrt("OPTS   1")
  } else {
    wrt("OPTS   0")
  }
  
  wrt("MONI    0")
  wrt("MARK    0")
  wrt("NOTE    Import the R data into MLwiN")
  wrt(paste("RSTA    '", dtafile, "'", sep = ""))
  
  wrt("NOTE    Correct column names")
  for (name in names(namemap)) {
    wrt(paste0("COLN '", namemap[[name]], "' b50"))
    wrt(paste0("DESC cb50 '", name, "'"))
    wrt(paste0("NAME cb50 '", name, "'"))
  }
  
  if ("class" %in% notation) {
    wrt("INDE 1")
  }
  if (!(D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed")) {
    wrt("NOTE   Specify the response variable")
    for (p in resp) wrt(paste("RESP    '", p, "'", sep = ""))
    wrt("")
  }
  
  if (D[1] == "Binomial") {
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:1
      if (!is.na(levID[ii])) 
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("")
    
    wrt("RDISt 1 0")
    wrt(paste("DOFFs 1 '", D[3], "'", sep = ""))
    if (D[2] == "logit") {
      wrt("LFUN 0")
      DD2 <- 0
    }
    if (D[2] == "probit") {
      wrt("LFUN 1")
      DD2 <- 1
    }
    if (D[2] == "cloglog") {
      wrt("LFUN 2")
      DD2 <- 2
    }
    
    interpos <- grep("\\:", expl)
    if (!isTRUE(oldsyntax) || length(interpos) == 0) {
      for (p in expl) {
        if (is.null(categ)) {
          wrt(paste("ADDT    '", p, "'", sep = ""))
        } else {
          if (sum(p == categ["var", ]) != 0) {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
            } else {
              wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                         ])]), sep = ""))
            }
          } else {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          }
        }
      }
    } else {
      exply <- expl[interpos]
      explx <- expl[-interpos]
      if (length(explx) > 0) {
        for (p in explx) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      }
      for (i in 1:length(exply)) {
        TT <- ""
        interx <- unlist(strsplit(exply[i], "\\:"))
        for (j in 1:length(interx)) {
          TT <- paste(TT, "'", interx[j], "' ", sep = "")
        }
        wrt(paste("ADDT    ", TT, sep = ""))
      }
      expl <- c(explx, exply)
    }
    wrt("")
  }
  
  if (D[1] == "Mixed") {
    nresp <- length(resp)
    for (ii in 1:nresp) wrt(paste("MVAR 1   '", resp[ii], "'", sep = ""))
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:2
      if (!is.na(levID[ii])) 
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("IDEN 1 'resp_indicator'")
    jj <- 1
    for (ii in 2:length(D)) {
      if (D[[ii]][1] == "Binomial") {
        wrt(paste("RDISt ", jj, " 0", sep = ""))
        wrt(paste("DOFFs ", jj, " '", D[[ii]][3], "'", sep = ""))
        if (D[[ii]][2] == "logit") {
          wrt("LFUN 0")
          DD2 <- 0
        }
        if (D[[ii]][2] == "probit") {
          wrt("LFUN 1")
          DD2 <- 1
        }
        if (D[[ii]][2] == "cloglog") {
          wrt("LFUN 2")
          DD2 <- 2
        }
      }
      if (D[[ii]][1] == "Poisson") {
        wrt(paste("RDISt ", jj, " 1", sep = ""))
        wrt("LFUN 3")
        DD2 <- 3
        DD2 <- 3
        if (!is.na(D[[ii]][3])) {
          wrt(paste("DOFFs ", jj, " '", D[[ii]][3], "'", sep = ""))
        }
      }
      jj <- jj + 1
    }
    wrt("")
    if (is.list(expl)) {
      if (!is.na(sep.coeff[1])) {
        interpos1 <- grep("\\:", sep.coeff)
        if (!isTRUE(oldsyntax) || length(interpos1) == 0) {
          for (x in 1:length(sep.coeff)) {
            p <- sep.coeff[x]
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        } else {
          exply <- sep.coeff[interpos1]
          explx <- sep.coeff[-interpos1]
          if (length(explx) > 0) {
            for (p in explx) {
              if (is.null(categ)) {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              } else {
                if (sum(p == categ["var", ]) != 0) {
                  if (is.na(categ["ref", which(p == categ["var", ])])) {
                    wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                  } else {
                    wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                               ])]), sep = ""))
                  }
                } else {
                  wrt(paste("ADDT    '", p, "'", sep = ""))
                }
              }
            }
          }
          for (i in 1:length(exply)) {
            TT <- ""
            interx <- unlist(strsplit(exply[i], "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          sep.coeff <- c(explx, exply)
        }
      }
      interpos2 <- grep("\\:", common.coeff)
      if (!isTRUE(oldsyntax) || length(interpos2) == 0) {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          
          p <- common.coeff[y]
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          p <- common.coeff[y]
          if (!(y %in% interpos2)) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          } else {
            TT <- ""
            interx <- unlist(strsplit(p, "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          wrt("RPAT")
          common.coeff[y] <- paste(p, ".", bb, sep = "")
        }
      }
    } else {
      interpos <- grep("\\:", expl)
      if (!isTRUE(oldsyntax) || length(interpos) == 0) {
        for (p in expl) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        exply <- expl[interpos]
        explx <- expl[-interpos]
        if (length(explx) > 0) {
          for (p in explx) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        }
        for (i in 1:length(exply)) {
          TT <- ""
          interx <- unlist(strsplit(exply[i], "\\:"))
          for (j in 1:length(interx)) {
            TT <- paste(TT, "'", interx[j], "' ", sep = "")
          }
          wrt(paste("ADDT    ", TT, sep = ""))
        }
        expl <- c(explx, exply)
      }
    }
    wrt("")
  }
  
  if (D[1] == "Multivariate Normal") {
    nresp <- length(resp)
    for (ii in 1:nresp) wrt(paste("MVAR 1   '", resp[ii], "'", sep = ""))
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:2
      if (!is.na(levID[ii])) 
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("IDEN 1 'resp_indicator'")
    wrt("LFUN 0")
    wrt("")
    if (is.list(expl)) {
      if (!is.na(sep.coeff[1])) {
        interpos1 <- grep("\\:", sep.coeff)
        if (!isTRUE(oldsyntax) || length(interpos1) == 0) {
          for (x in 1:length(sep.coeff)) {
            p <- sep.coeff[x]
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        } else {
          exply <- sep.coeff[interpos1]
          explx <- sep.coeff[-interpos1]
          if (length(explx) > 0) {
            for (p in explx) {
              if (is.null(categ)) {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              } else {
                if (sum(p == categ["var", ]) != 0) {
                  if (is.na(categ["ref", which(p == categ["var", ])])) {
                    wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                  } else {
                    wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                               ])]), sep = ""))
                  }
                } else {
                  wrt(paste("ADDT    '", p, "'", sep = ""))
                }
              }
            }
          }
          for (i in 1:length(exply)) {
            TT <- ""
            interx <- unlist(strsplit(exply[i], "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          sep.coeff <- c(explx, exply)
        }
      }
      interpos2 <- grep("\\:", common.coeff)
      if (!isTRUE(oldsyntax) || length(interpos2) == 0) {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          
          p <- common.coeff[y]
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          p <- common.coeff[y]
          if (!(y %in% interpos2)) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          } else {
            TT <- ""
            interx <- unlist(strsplit(p, "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          wrt("RPAT")
          common.coeff[y] <- paste(p, ".", bb, sep = "")
        }
      }
    } else {
      interpos <- grep("\\:", expl)
      if (!isTRUE(oldsyntax) || length(interpos) == 0) {
        for (p in expl) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        exply <- expl[interpos]
        explx <- expl[-interpos]
        if (length(explx) > 0) {
          for (p in explx) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        }
        for (i in 1:length(exply)) {
          TT <- ""
          interx <- unlist(strsplit(exply[i], "\\:"))
          for (j in 1:length(interx)) {
            TT <- paste(TT, "'", interx[j], "' ", sep = "")
          }
          wrt(paste("ADDT    ", TT, sep = ""))
        }
        expl <- c(explx, exply)
      }
    }
    wrt("")
  }
  
  if (D[1] == "Multinomial") {
    wrt("LINK 2 G21")
    wrt("NAME G21[1] 'resp' G21[2] 'resp_indicator'")
    wrt("LINK 0 G21")
    wrt(paste("MNOM ", as.numeric(D[4]), " '", resp, "' ", "'resp' 'resp_indicator' ", as.numeric(D[5]), sep = ""))
    wrt("RESP   'resp'")
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:c(nlev - 1)) {
      aa <- nlev:1
      if (!is.na(levID[ii])) 
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("IDEN 1 'resp_indicator'")
    wrt("")
    
    if (as.numeric(D[4]) == 0) 
      wrt("RDISt 1 4") else wrt("RDISt 1 5")
    wrt(paste("DOFFs 1 '", D[3], "'", sep = ""))
    if (D[2] == "logit") {
      wrt("LFUN 0")
      DD2 <- 0
    }
    if (D[2] == "probit") {
      wrt("LFUN 1")
      DD2 <- 1
    }
    if (D[2] == "cloglog") {
      wrt("LFUN 2")
      DD2 <- 2
    }
    
    wrt("")
    if (is.list(expl)) {
      if (!is.na(sep.coeff[1])) {
        interpos1 <- grep("\\:", sep.coeff)
        if (!isTRUE(oldsyntax) || length(interpos1) == 0) {
          for (x in 1:length(sep.coeff)) {
            p <- sep.coeff[x]
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        } else {
          exply <- sep.coeff[interpos1]
          explx <- sep.coeff[-interpos1]
          if (length(explx) > 0) {
            for (p in explx) {
              if (is.null(categ)) {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              } else {
                if (sum(p == categ["var", ]) != 0) {
                  if (is.na(categ["ref", which(p == categ["var", ])])) {
                    wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                  } else {
                    wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                               ])]), sep = ""))
                  }
                } else {
                  wrt(paste("ADDT    '", p, "'", sep = ""))
                }
              }
            }
          }
          for (i in 1:length(exply)) {
            TT <- ""
            interx <- unlist(strsplit(exply[i], "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          sep.coeff <- c(explx, exply)
        }
      }
      interpos2 <- grep("\\:", common.coeff)
      if (!isTRUE(oldsyntax) || length(interpos2) == 0) {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          
          p <- common.coeff[y]
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        for (y in 1:length(common.coeff)) {
          p <- common.coeff[y]
          len.common.id <- length(common.coeff.id[y, ])
          tt <- "RPAT    "
          aa <- 1:len.common.id
          partname <- aa[rep(1, len.common.id) == common.coeff.id[y, ]]
          bb <- ""
          for (ii in aa) bb <- paste(bb, partname[ii], sep = "")
          for (ii in 1:len.common.id) tt <- paste(tt, common.coeff.id[y, ii])
          wrt(tt)
          p <- common.coeff[y]
          if (!(y %in% interpos2)) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          } else {
            TT <- ""
            interx <- unlist(strsplit(p, "\\:"))
            for (j in 1:length(interx)) {
              TT <- paste(TT, "'", interx[j], "' ", sep = "")
            }
            wrt(paste("ADDT    ", TT, sep = ""))
          }
          wrt("RPAT")
          common.coeff[y] <- paste(p, ".", bb, sep = "")
        }
      }
    } else {
      interpos <- grep("\\:", expl)
      if (!isTRUE(oldsyntax) || length(interpos) == 0) {
        for (p in expl) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      } else {
        exply <- expl[interpos]
        explx <- expl[-interpos]
        if (length(explx) > 0) {
          for (p in explx) {
            if (is.null(categ)) {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            } else {
              if (sum(p == categ["var", ]) != 0) {
                if (is.na(categ["ref", which(p == categ["var", ])])) {
                  wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
                } else {
                  wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                             ])]), sep = ""))
                }
              } else {
                wrt(paste("ADDT    '", p, "'", sep = ""))
              }
            }
          }
        }
        for (i in 1:length(exply)) {
          TT <- ""
          interx <- unlist(strsplit(exply[i], "\\:"))
          for (j in 1:length(interx)) {
            TT <- paste(TT, "'", interx[j], "' ", sep = "")
          }
          wrt(paste("ADDT    ", TT, sep = ""))
        }
        expl <- c(explx, exply)
      }
    }
    wrt("")
  }
  
  if (D[1] == "Normal") {
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:1
      if (!is.na(levID[ii])) 
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("")
    
    wrt("NOTE   Specify covariate(s) used anywhere in the model")
    interpos <- grep("\\:", expl)
    if (!isTRUE(oldsyntax) || length(interpos) == 0) {
      for (p in expl) {
        if (is.null(categ)) {
          wrt(paste("ADDT    '", p, "'", sep = ""))
        } else {
          if (sum(p == categ["var", ]) != 0) {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
            } else {
              wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                         ])]), sep = ""))
            }
          } else {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          }
        }
      }
    } else {
      exply <- expl[interpos]
      explx <- expl[-interpos]
      if (length(explx) > 0) {
        for (p in explx) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      }
      for (i in 1:length(exply)) {
        TT <- ""
        interx <- unlist(strsplit(exply[i], "\\:"))
        for (j in 1:length(interx)) {
          TT <- paste(TT, "'", interx[j], "' ", sep = "")
        }
        wrt(paste("ADDT    ", TT, sep = ""))
      }
      expl <- c(explx, exply)
    }
    wrt("")
  }
  if (D[1] == "Poisson") {
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:1
      if (!is.na(levID[ii])) 
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("")
    
    wrt("RDISt 1 1")
    wrt("LFUN 3")
    DD2 <- 3
    if (!is.na(D[3])) {
      wrt(paste("DOFFs 1 '", D[3], "'", sep = ""))
    }
    interpos <- grep("\\:", expl)
    if (!isTRUE(oldsyntax) || length(interpos) == 0) {
      for (p in expl) {
        if (is.null(categ)) {
          wrt(paste("ADDT    '", p, "'", sep = ""))
        } else {
          if (sum(p == categ["var", ]) != 0) {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
            } else {
              wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                         ])]), sep = ""))
            }
          } else {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          }
        }
      }
    } else {
      exply <- expl[interpos]
      explx <- expl[-interpos]
      if (length(explx) > 0) {
        for (p in explx) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      }
      for (i in 1:length(exply)) {
        TT <- ""
        interx <- unlist(strsplit(exply[i], "\\:"))
        for (j in 1:length(interx)) {
          TT <- paste(TT, "'", interx[j], "' ", sep = "")
        }
        wrt(paste("ADDT    ", TT, sep = ""))
      }
      expl <- c(explx, exply)
    }
    wrt("")
  }
  
  
  if (D[1] == "Negbinom") {
    wrt("NOTE   Specify the level identifier(s)")
    for (ii in 1:nlev) {
      aa <- nlev:1
      if (!is.na(levID[ii])) 
        wrt(paste("IDEN ", aa[ii], "    '", levID[ii], "'", sep = ""))
    }
    wrt("")
    
    wrt("RDISt 1 2")
    wrt("LFUN 3")
    DD2 <- 3
    if (!is.na(D[3])) {
      wrt(paste("DOFFs 1 '", D[3], "'", sep = ""))
    }
    interpos <- grep("\\:", expl)
    if (!isTRUE(oldsyntax) || length(interpos) == 0) {
      for (p in expl) {
        if (is.null(categ)) {
          wrt(paste("ADDT    '", p, "'", sep = ""))
        } else {
          if (sum(p == categ["var", ]) != 0) {
            if (is.na(categ["ref", which(p == categ["var", ])])) {
              wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
            } else {
              wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                         ])]), sep = ""))
            }
          } else {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          }
        }
      }
    } else {
      exply <- expl[interpos]
      explx <- expl[-interpos]
      if (length(explx) > 0) {
        for (p in explx) {
          if (is.null(categ)) {
            wrt(paste("ADDT    '", p, "'", sep = ""))
          } else {
            if (sum(p == categ["var", ]) != 0) {
              if (is.na(categ["ref", which(p == categ["var", ])])) {
                wrt(paste("ADDT    '", p, "' ", -1e+07, sep = ""))
              } else {
                wrt(paste("ADDT    '", p, "' ", which(levels(indata[, p]) == categ["ref", which(p == categ["var", 
                                                                                                           ])]), sep = ""))
              }
            } else {
              wrt(paste("ADDT    '", p, "'", sep = ""))
            }
          }
        }
      }
      for (i in 1:length(exply)) {
        TT <- ""
        interx <- unlist(strsplit(exply[i], "\\:"))
        for (j in 1:length(interx)) {
          TT <- paste(TT, "'", interx[j], "' ", sep = "")
        }
        wrt(paste("ADDT    ", TT, sep = ""))
      }
      expl <- c(explx, exply)
    }
    wrt("")
  }
  
  if (is.list(nonfp)) {
    wrt("NOTE Turn off the fixed part of the explotary varible(s)")
    nonfp.sep <- nonfp$nonfp.sep
    nonfp.common <- nonfp$nonfp.common
    if (!is.na(nonfp.sep[1])) {
      interpos <- grep("\\:", nonfp.sep)
      if (!isTRUE(oldsyntax) || length(interpos) == 0) {
        for (p in nonfp.sep) wrt(paste("FPAR 0  '", p, "'", sep = ""))
      } else {
        for (i in 1:length(nonfp.sep)) {
          if (i %in% interpos && oldsyntax) {
            wrt(paste("FPAR 0  '", gsub("\\:", "\\.", nonfp.sep[i]), "'", sep = ""))
          } else {
            wrt(paste("FPAR 0  '", nonfp.sep[i], "'", sep = ""))
          }
        }
      }
    }
    if (!is.na(nonfp.common[1])) {
      interpos <- grep("\\:", nonfp.common)
      if (!isTRUE(oldsyntax) || length(interpos) == 0) {
        for (p in nonfp.common) wrt(paste("FPAR 0  '", p, "'", sep = ""))
      } else {
        for (i in 1:length(nonfp.common)) {
          if (i %in% interpos && oldsyntax) {
            wrt(paste("FPAR 0  '", gsub("\\:", "\\.", nonfp.common[i]), "'", sep = ""))
          } else {
            wrt(paste("FPAR 0  '", nonfp.common[i], "'", sep = ""))
          }
        }
      }
    }
  } else {
    if (!is.na(nonfp[1])) {
      wrt("NOTE Turn off the fixed part of the explotary varible(s)")
      for (p in nonfp) {
        if (oldsyntax){
          wrt(paste("FPAR 0  '", gsub("\\:", "\\.", p), "'", sep = ""))
        } else {
          wrt(paste("FPAR 0  '", p, "'", sep = ""))
        }
      }
    }
  }
  
  wrt("")
  wrt("NOTE   Specify random part covariate(s)")
  if (nrp > 0) {
    for (ii in 1:nrp) {
      for (p in rp[[ii]]) {
        if (oldsyntax) {
          p <- gsub("\\:", "\\.", p)
        }
        wrt(paste("SETV  ", as.numeric(sub("rp", "", rp.names[ii])), "   '", p, "'", sep = ""))
      }
    }
  }
  if (!is.null(clre)) {
    nclre <- ncol(clre)
    for (ii in 1:nclre) {
      if (oldsyntax) {
        wrt(paste("CLRE  ", as.numeric(clre[1, ii]), " '", gsub("\\:", "\\.", clre[2, ii]), "' '",
                  gsub("\\:", "\\.", clre[3, ii]), "'", sep = ""))
      } else {
        wrt(paste("CLRE  ", as.numeric(clre[1, ii]), " '",  clre[2, ii], "' '", clre[3, ii], "'", sep = ""))
      }
    }
  }
  
  nexpl <- length(expl)
  wrt("")
  
  if (!is.null(fcon)) {
    wrt("FCON b1001")
    wrt(paste("JOIN cb1001", paste(as.vector(fcon), collapse = " "), "cb1001"))
  }
  
  if (!is.null(rcon)) {
    wrt("RCON b1002")
    wrt(paste("JOIN cb1002", paste(as.vector(rcon), collapse = " "), "cb1002"))
  }
  
  wrt("NOTE   Set estimation method")
  if (Meth != 2) {
    wrt(paste("METH", Meth))
  }
  
  for (ii in 1:nlev) {
    wrt(paste("RESEt", ii, reset[ii]))
  }
  
  wrt(paste("LINE ", nonlinear[1], nonlinear[2]))
  if (isTRUE(extra)) {
    wrt("EXTRa 1")
  }
  wrt("")
  
  if (isTRUE(fpsandwich)) {
    wrt("NOTE   Turn on sandwich estimators for the fixed part parameter standard errors")
    wrt("FSDE 2")
  } else {
    wrt("FSDE 0")
  }
  wrt("")
  
  if (isTRUE(rpsandwich)) {
    wrt("NOTE   Turn on sandwich estimators for the random part parameter standard errors")
    wrt("RSDE 2")
  } else {
    wrt("RSDE 0")
  }
  wrt("")
  
  if (!is.null(weighting)) {
    for (i in 1:length(weighting$weightvar)) {
      wtlev <- (length(weighting$weightvar) - i) + 1
      if (!is.na(weighting$weightvar[[i]])) {
        wrt(paste("NOTE   Specify sampling weights at level", wtlev))
        wrt(paste("WEIG ", wtlev, " ", 1, " '", weighting$weightvar[[i]], "'", sep = ""))
      } else {
        wrt(paste("NOTE   Specify equal weights at level", wtlev))
        wrt(paste("WEIG ", wtlev, " ", 1, sep = ""))
      }
      wrt("")
      if (isTRUE(weighting$standardised)) {
        wrt("NOTE   Standardised weighting")
        wrt("LINK 1 G21")
        wrt(paste("WEIG ", wtlev, " ", 2, " G21[1]", sep = ""))
        wrt("FILL G21")
        wrt("LINK 0 G21")
        wrt("WEIG 2")
        wrt("")
      } else {
        wrt("NOTE   Raw weighting")
        wrt("LINK 1 G21")
        wrt(paste("WEIG ", wtlev, " ", 2, " G21[1]", sep = ""))
        wrt("FILL G21")
        wrt("LINK 0 G21")
        wrt("WEIG 1")
      }
    }
    wrt("")
    if (isTRUE(weighting$standardised)) {
      wrt("NOTE   Create the standardised weights")
      wrt("WEIG")
    } else {
      wrt("WEIG 0")
    }
    wrt("")
  }
  if (D[1] == "Normal") {
    wrt("PREF   0")
    wrt("POST   0")
  }

  if ("simple" %in% notation) {
    wrt("EXISt 'cons' b1000")
    wrt("SWITch b1000")
    wrt("CASE 0:")
    wrt("EXISt 'Intercept' b1000")
    wrt("SWITch b1000")
    wrt("CASE 1:")
    wrt("COLN 'Intercept' b1000")
    wrt("NAME cb1000 'cons'")
    wrt("ENDSWITch")
    wrt("ENDSWITch")
    wrt("NOTA 1")
  }
  
  wrt("NAME   c1098 '_FP_b'")
  wrt("NAME   c1099 '_FP_v'")
  wrt("NAME   c1096 '_RP_b'")
  wrt("NAME   c1097 '_RP_v'")
  
  wrt("NOTE   Fit the model")
  wrt("ECHO 1")
  wrt("BATC 1")
  wrt("MAXI 2")
  wrt("STAR")
  
  if (!is.null(startval)) {
    if (!is.null(startval$FP.b)) {
      wrt(paste("JOIN ", paste(startval$FP.b, collapse = " "), " '_FP_b'", sep = ""))
    }
    if (!is.null(startval$FP.v)) {
      wrt(paste("JOIN ", paste(startval$FP.v[!lower.tri(startval$FP.v)], collapse = " "), " '_FP_v'", sep = ""))
    }
    if (!is.null(startval$RP.b)) {
      wrt(paste("JOIN ", paste(startval$RP.b, collapse = " "), " '_RP_b'", sep = ""))
    }
    if (!is.null(startval$RP.v)) {
      wrt(paste("JOIN ", paste(startval$RP.v[!lower.tri(startval$RP.v)], collapse = " "), " '_RP_v'", sep = ""))
    }
  }

  if (!is.null(saveworksheet)) {
    wrt(paste0("STOR ", saveworksheet))
  }
  
  if (debugmode) {
    wrt("WSET 15 1")
    wrt("EXPA 3")
    wrt("ESTM 2")
    wrt("PAUS")
  }
  
  wrt(paste("TOLE", convtol))
  wrt(paste("MAXI", maxiter))
  wrt("NEXT")
  wrt("ECHO 0")
  wrt("MONI 1")
  wrt("ITNU 0 b21")
  wrt("CONV b22")
  wrt("")
  
  
  
  if (debugmode) {
    wrt("NOTE   Open the equations window")
    wrt("WSET 15 1")
    wrt("EXPA 3")
    wrt("ESTM 2")
    wrt("PAUS")
  }
  wrt("NOTE    *****************************************************************")
  wrt("")
  
  
  wrt("NOTE    *****************************************************************")
  wrt("NOTE       Export the model results to R")
  wrt("NOTE    *****************************************************************")
  wrt("LINK 1 G21")
  wrt("NAME   G21[1] '_Stats'")
  if (D[1] == "Multinomial" || D[1] == "Multivariate Normal" || D[1] == "Mixed") {
    wrt("NOBS 2 b31 b32")
  } else {
    wrt("NOBS 1 b31 b32")
  }
  wrt("EDIT 1 '_Stats' b31")
  wrt("EDIT 2 '_Stats' b32")
  if (D[1] == "Normal" || D[1] == "Multivariate Normal") {
    wrt("LIKE   b100")
  }
  wrt("EDIT 3 '_Stats' b100")
  wrt("EDIT 7 '_Stats' b21")
  wrt("EDIT 8 '_Stats' b22")
  wrt("NAME   c1098 '_FP_b'")
  wrt("NAME   c1099 '_FP_v'")
  wrt("NAME   c1096 '_RP_b'")
  wrt("NAME   c1097 '_RP_v'")
  wrt("NAME   c1094 '_esample'")
  wrt("SUM '_esample' b1")
  wrt("EDIT 9 '_Stats' b1")
  wrt(paste("PSTA '", IGLSfile, "' ", "'_FP_b' ", "'_FP_v' ", "'_RP_b' ", "'_RP_v' ", "'_Stats'", sep = ""))
  
  calcresiduals <- function(level, displevel, rpx, resioptions, clre = clre) {
    wrt("")
    if (!("norecode" %in% resioptions)) {
      wrt("MISR 0")
    }
    
    len.rpx <- length(rpx)
    if (level == 2 & D[[1]][1] == "Mixed") {
      for (i in 2:length(D)) {
        if (D[[i]][1] == "Binomial" | D[[i]][1] == "Poisson") {
          len.rpx <- len.rpx + 1
        }
      }
    }
    
    wrt(paste("LINK", len.rpx, "G21"))
    for (k in 1:len.rpx) {
      wrt(paste0("NAME G21[", k, "] '", shortname("lev_", displevel, "_resi_est_", rpx[k]), "'"))
      wrt(paste0("DESC G21[", k, "] ", "'residual estimates'"))
    }
    
    wrt(paste("LINK", len.rpx, "G22"))
    if ("variance" %in% resioptions) {
      for (k in 1:len.rpx) {
        wrt(paste0("NAME G22[", k, "] '", shortname("lev_", displevel, "_resi_var_", rpx[k]), "'"))
        wrt(paste0("DESC G22[", k, "] ", "'residual variance'"))
      }
    } else {
      for (k in 1:len.rpx) {
        wrt(paste0("NAME G22[", k, "] '", shortname("lev_", displevel, "_resi_se_", rpx[k]), "'"))
        wrt(paste0("DESC G22[", k, "] ", "'residual standard error'"))
      }
    }
    wrt("RFUN")
    wrt("ROUT G21 G22")
    wrt("")
    
    wrt(paste("RLEV   ", level, sep = ""))
    wrt("RCOV   1")
    
    outgroups <- c("G21", "G22")
    
    if ("standardised" %in% resioptions || "deletion" %in% resioptions || "leverage" %in% resioptions) {
      wrt(paste("LINK", len.rpx, "G23"))
      for (k in 1:len.rpx) {
        wrt(paste0("NAME G23[", k, "] '", shortname("lev_", displevel, "_std_resi_est_", rpx[k]), "'"))
        wrt(paste0("DESC G23[", k, "] ", "'std standardised residual'"))
      }
      
      wrt("RTYP   0")
      wrt("RESI")
      for (k in 1:len.rpx) {
        wrt(paste0("CALC G23[", k, "]=G21[", k, "]/sqrt(G22[", k, "])", sep = ""))
      }
      if ("standardised" %in% resioptions) {
        outgroups <- c(outgroups, "G23")
      }
    }
    
    # NOTE leverage depends on residuals with rtype 0
    if ("leverage" %in% resioptions || "influence" %in% resioptions) {
      # influence requires leverage to be calculated
      wrt(paste("LINK", len.rpx, "G24"))
      for (k in 1:len.rpx) {
        wrt(paste0("NAME G24[", k, "] '", shortname("lev_", displevel, "_resi_leverage_", rpx[k]), "'"))
        wrt(paste0("DESC G24[", k, "] ", "'leverage residual'"))
      }
      outgroups <- c(outgroups, "G24")
      
      for (k in 1:len.rpx) {
        wrt("LINK 1 G25")
        wrt(paste0("OMEGa ", level, " '", rpx[k], "' G25[1]"))  # retrieve variance for corresponding random parameter
        wrt("PICK 1 G25[1] b50")
        wrt("ERASe G25[1]")
        wrt("LINK 0 G25")
        wrt(paste0("CALC G24[", k, "]=1-sqrt(G22[", k, "])/sqrt(b50)"))
      }
      
      if (!("standardised" %in% resioptions) && !("deletion" %in% resioptions)) {
        wrt("ERAS G23")
        wrt("LINK 0 G23")
        outgroups <- outgroups[outgroups != "G23"]
      }
    }
    
    wrt("RTYP   1")  # Compute comparative variances
    wrt("RESI")
    
    if (!("variance" %in% resioptions)) {
      for (k in 1:len.rpx) {
        wrt(paste0("CALC G22[", k, "]=sqrt(G22[", k, "])"))  # Convert the variances to standard errors
      }
    }
    
    if ("deletion" %in% resioptions || "influence" %in% resioptions) {
      wrt(paste("LINK", len.rpx, "G25"))
      for (k in 1:len.rpx) {
        wrt(paste0("NAME G25[", k, "] '", shortname("lev_", displevel, "_resi_deletion_", rpx[k]), "'"))
        wrt(paste0("DESC G25[", k, "] ", "'deletion residual'"))
      }
      outgroups <- c(outgroups, "G25")
      
      wrt(paste("NOBS ", level, "b31 b32"))
      for (k in 1:len.rpx) {
        wrt(paste0("CALC G25[", k, "]=G23[", k, "]/ sqrt((b31 - 1 - G23[", k, "]^2)/(b31 - 2))"))
      }
      
      if (!("standardised" %in% resioptions) && !("leverage" %in% resioptions)) {
        wrt("ERASe G23")
        wrt("LINK 0 G23")
        outgroups <- outgroups[outgroups != "G23"]
      }
    }
    
    if ("influence" %in% resioptions) {
      wrt(paste("LINK", len.rpx, "G26"))
      for (k in 1:len.rpx) {
        wrt(paste0("NAME G26[", k, "] '", shortname("lev_", displevel, "_resi_influence_", rpx[k]), "'"))
        wrt(paste0("DESC G26[", k, "] ", "'influence residual'"))
      }
      outgroups <- c(outgroups, "G26")
      
      for (k in 1:len.rpx) {
        wrt(paste0("SUM    G24[", k, "]", " b50"))
        wrt(paste0("CALC   G26[", k, "]=G24[", k, "]/b50"))
        wrt(paste0("CALC   G26[", k, "]=sqrt(G26[", k, "]/(1-G26[", k, "]))*abso(G25[", k, "])"))
      }
      
      if (!("deletion" %in% resioptions)) {
        wrt("ERASE G25")
        wrt("LINK 0 G25")
        outgroups <- outgroups[outgroups != "G25"]
      }
      if (!("leverage" %in% resioptions)) {
        wrt("ERASE G26")
        wrt("LINK 0 G26")
        outgroups <- outgroups[outgroups != "G26"]
      }
    }
    
    if ("sampling" %in% resioptions) {
      numele <- (len.rpx * (len.rpx + 1))/2
      wrt(paste0("LINK ", numele, " G26"))
      k <- 1
      for (i in 1:len.rpx) {
        for (j in 1:i) {
          if (i == j) {
            wrt(paste0("NAME G26[", k, "] '", shortname("lev_", displevel, "_resi_var_", rpx[i]), "'"))
            wrt(paste0("DESC G26[", k, "] ", "'sampling variance'"))
          } else {
            wrt(paste0("NAME G26[", k, "] '", shortname("lev_", displevel, "_resi_cov_", rpx[i], "_", rpx[j]), "'"))
            wrt(paste0("DESC G26[", k, "] ", "'sampling covariance'"))
          }
          k <- k + 1
        }
      }
      outgroups <- c(outgroups, "G26")
      
      
      
      wrt(paste0("LINK ", len.rpx, " G27"))
      wrt(paste0("LINK ", 1, " G28"))
      wrt(paste0("LINK ", 1, " G29"))
      
      wrt(paste0("NAME G29[1] ", paste0("'lev_", displevel, "_resi_cov'")))
      wrt(paste0("DESC G29[1] ", paste0("'sampling var(cov)'")))
      outgroups <- c(outgroups, "G29")
      
      
      wrt("RFUN")
      wrt("ROUT G27 G29")
      wrt("RCOV 2")
      wrt("RESI")
      wrt("")
      
      # NOTE: This is square rooted, as the residual covariances are sometimes negative
      wrt(paste0("NOBS ", level, " b31 b32"))
      wrt(paste0("CODE ", numele, " 1 b31 G28[1]"))
      wrt("SPLIt G29 G28 G26")
      wrt("ERAS G28 G29")
      wrt("LINK 0 G28")
      wrt("LINK 0 G29")
      wrt("ERAS G27")
      wrt("LINK 0 G27")
      outgroups <- outgroups[outgroups != "G29"]
    }
    
    wrt("")
    wrt(paste0("NOBS ", level, " b30 b31"))
    wrt("LINK 1 G29")
    wrt("GENE 1 b30 1 G29[1]")
    wrt(paste0("NAME G29[1] 'lev_", level, "_residualid'"))
    outgroups <- c(outgroups, "G29")
    if (!("norecode" %in% resioptions)) {
      wrt("MISR 1")
    }
    wrt("")
    wrt("LINK 0 G30")
    for (i in 1:length(outgroups)) {
      wrt(paste("GSET 2 G30", outgroups[i], "G30"))
      wrt(paste("LINK 0", outgroups[i]))
    }
  }
  
  
  if (resi.store && nrp > 0) {
    for (j in nrp:1) {
      rpx <- rp[[j]]
      len.rpx <- length(rp[[j]])
      wrt(paste("NOTE Calculate level ", as.numeric(sub("rp", "", rp.names[j])), " residuals", sep = ""))
      levtt <- as.numeric(sub("rp", "", rp.names[j]))
      displevel <- levtt
      if (D[1] == "Multivariate Normal" || D[1] == "Mixed" || D[1] == "Multinomial") {
        displevel <- displevel - 1
      }
      calcresiduals(levtt, displevel, rpx, resioptions, clre = clre)
      wrt(paste("PSTA '", resifile[j], "' G30", sep = ""))
    }
    wrt("ERAS G30")
    wrt("LINK 0 G30")
    
  }
  
  if (debugmode) {
    wrt("WPMT 'Do you want to close MLwiN?' b50")
    wrt("SWITch b50")
    wrt("  CASE 1:")
    wrt("    EXIT ")
  } else {
    wrt("EXIT")
  }
  return(namemap)
} 
