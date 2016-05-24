#' An internal function, allowing back-compatibility, which translates a model
#' formula from a formula object or character string into an R list object.
#'
#' Supports Formula syntax as used in earlier (<0.8-0) versions of \pkg{R2MLwiN}. A model
#' formula, as a formula object (or a character string) is translated into
#' an R list object. Called by \code{\link{runMLwiN}} if \code{oldsyntax = TRUE}
#' (when user specifies \code{levID} not \code{NULL} in \code{\link{runMLwiN}} function
#' call). For corresponding function supporting new syntax, see
#' \code{\link{Formula.translate}}.
#'
#' @param Formula A formula object (or a character string) specifying a
#' multilevel model. See \code{Value} for details.
#' @param levID A character (vector) specifying the level ID(s).
#' @param D A character string/vector specifying the distribution to be
#' modelled, which can include \code{'Normal'} (the default), \code{'Binomial'},
#' \code{'Poisson'}, \code{'Negbinom'}, \code{'Unordered Multinomial'},
#' \code{'Ordered Multinomial'}, \code{'Multivariate Normal'}, or \code{'Mixed'}.
#' @param indata A data.frame object containing the data to be modelled.
#'
#' @details If \code{Formula} is a character string, then the following
#' syntax applies:
#' \itemize{
#' \item \code{~} A tilde is used to separate response variable(s) and
#' explanatory variable(s).
#' \item \code{()} Round brackets are used to specify each random
#' variable in the model together with its fixed/random part information.
#' \item \code{|} Separates explanatory variable(s) (placed to the right of
#' \code{|}) from the fixed/random part information (placed to the left of
#' \code{|}) when placed within \code{()}.
#' \item \code{[]} When placed
#' immediately after an explanatory variable, indicates that the variable is
#' categorical. The string in the \code{[]} represents the reference category;
#' if empty, no reference category is used; See note.
#' \item \code{:} Indicates an interaction term: i.e. the variables adjacent to \code{:}, and
#' separated by it, are interacted with each other.
#' \item \code{0} When placed to the left of \code{|} within \code{()} indicates that the variables
#' to the right of \code{|} within the same \code{()} are to be added to the
#' fixed part of the model.
#' \item \code{1} When placed to the left of
#' \code{|} within \code{()} indicates that the coefficients of the variables
#' placed to the right of \code{|} within the same \code{()} are to be allowed
#' to randomly vary at level 1 (and so on for \code{2} for level 2, \code{3} for level 3, etc.)
#' \item \code{0s/0c} When placed to the left of \code{|} within \code{()}
#' indicates that separate (hence \code{s}) / common (hence \code{c})
#' coefficients for the variables to the right of \code{|} within the same
#' \code{()} are to be added to the fixed part (hence \code{0}) of multivariate
#' normal, multinomial and mixed responses models.
#' \item \code{2s/2c} When
#' placed to the left of \code{|} within \code{()} indicates that separate
#' (hence \code{s}) / common (hence \code{c}) coefficients for the variables to
#' the right of \code{|} within the same \code{()} are to be added to the
#' random part of the model, and allowed to vary at level 2; applies to
#' multivariate normal, multinomial and mixed responses models only.
#' \item \code{\{\}} gives a vector of binary indicators specifying a
#' common coefficient. 1 is to include the component at the corresponding
#' positions; zero otherwise. These digits are separated by commas; applies to
#' multivariate normal, multinomial and mixed responses models only.
#' \item \code{.} Used for adding a separate coefficient for a particular
#' component at a specific level; applies to multivariate normal, multinomial
#' and mixed responses models only
#' }
#'
#' If \code{Formula} is a formula object, \code{0s/0c}, \code{2s/2c}, .... and
#' \code{\{\}} have to be replaced by \code{`0s`/`0c`}, \code{`2s`/`2c`}, ....
#' and \code{()} respectively. Other syntax remains the same.
#'
#' @return Outputs an R list object, which is then used as the input for
#' \code{\link{write.IGLS}} and/or \code{\link{write.MCMC}}.
#'
#' @note Note that some characters listed above have special meanings in the
#' formula, so avoid using them when you name the random variable. Alphanumeric
#' characters (i.e. \code{[:alnum:]}) are recommended for naming the random
#' variable. They are also recommended for naming a reference category, inside
#' \code{[]}.  Note: use \code{[]} notation only in the fixed part when
#' there is no categorical variable in the random effects. If there is one in
#' the random part, the categorical variable has to be converted into a set of
#' binary variables (e.g., using \code{\link{Untoggle}}).
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso
#' \code{\link{runMLwiN}}, \code{\link{write.IGLS}}, \code{\link{write.MCMC}}, \code{\link{Formula.translate}}
#'

Formula.translate.compat <- function(Formula, levID, D = "Normal", indata) {
  nlev <- length(levID)
  cc <- c(0:nlev)
  if (is.character(Formula)) {
    Formula <- gsub("\\{", "\\(", Formula)
    Formula <- gsub("\\}", "\\)", Formula)
    Formula <- gsub("[[:space:]]", "", Formula)
    if (sum(grepl("\\(+[[:digit:]]+[[:alpha:]]+\\|", Formula)) > 0) {
      for (i in cc) {
        Formula <- sub(paste(i, "s\\|", sep = ""), paste("\\`", i, "s`\\|", sep = ""), Formula)
        Formula <- sub(paste(i, "c\\|", sep = ""), paste("\\`", i, "c`\\|", sep = ""), Formula)
      }
    }
    Formula <- as.formula(Formula)
  }
  Terms <- terms.formula(Formula, keep.order = TRUE)
  resp <- rownames(attr(Terms, "factors"))[attr(Terms, "response")]
  resp <- gsub("[[:space:]]", "", resp)
  left <- attr(Terms, "term.labels")
  left <- gsub("\\(", "\\{", left)
  left <- gsub("\\)", "\\}", left)
  left <- gsub("[[:space:]]", "", left)
  if (sum(grepl("\\`+[[:digit:]]+[[:alpha:]]+\\`+\\|", left)) > 0) {
    for (i in cc) {
      left <- sub(paste("\\`", i, "s`\\|", sep = ""), paste(i, "s\\|", sep = ""), left)
      left <- sub(paste("\\`", i, "c`\\|", sep = ""), paste(i, "c\\|", sep = ""), left)
    }
  }
  
  if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial" || D[1] == "Multivariate Normal" || D[1] ==
        "Mixed") {
    nlev <- length(levID)
    cc <- c(0:nlev)
    cflag <- 0
    if (sum(grepl("\\(+[[:digit:]]+[[:alpha:]]+\\|", left)) > 0)
      cflag <- 1
    if (cflag == 0) {
      for (i in cc) {
        left <- sub(paste(i, "\\|", sep = ""), paste(i, "s\\|", sep = ""), left)
      }
    }
    
    if (D[1] == "Multivariate Normal") {
      resp <- sub("c\\(", "", resp)
      resp <- sub("\\)", "", resp)
      resp <- strsplit(resp, ",")[[1]]
    }
    if (D[1] == "Unordered Multinomial") {
      D <- rep(NA, 5)
      resp <- regmatches(resp, regexec("([[:alnum:]]+)\\(([[:alnum:][:space:],_]+)\\)", resp))[[1]]
      link <- resp[2]
      resp <- strsplit(resp[3], ",")[[1]]
      D[1] <- "Unordered Multinomial"
      names(D)[1] <- "distr"
      D[2] <- link
      names(D)[2] <- "link"
      D[3] <- resp[2]
      names(D)[3] <- "denominator"
      D[4] <- 0
      names(D)[4] <- "mode"
      if (resp[3] %in% indata[[resp[1]]]) {
        if (is.factor(indata[[resp[1]]])) {
          D[5] <- which(resp[3] == levels(indata[[resp[1]]]))
        } else {
          D[5] <- which(resp[3] == sort(unique(indata[[resp[1]]])))
        }
      }
      names(D)[5] <- "ref.cat"
      resp <- resp[1]
    }
    if (D[1] == "Ordered Multinomial") {
      D <- rep(NA, 5)
      resp <- regmatches(resp, regexec("([[:alnum:]]+)\\(([[:alnum:][:space:],_]+)\\)", resp))[[1]]
      link <- resp[2]
      resp <- strsplit(resp[3], ",")[[1]]
      D[1] <- "Ordered Multinomial"
      names(D)[1] <- "distr"
      D[2] <- link
      names(D)[2] <- "link"
      D[3] <- resp[2]
      names(D)[3] <- "denominator"
      D[4] <- 1
      names(D)[4] <- "mode"
      if (resp[3] %in% indata[[resp[1]]]) {
        if (is.factor(indata[[resp[1]]])) {
          D[5] <- which(resp[3] == levels(indata[[resp[1]]]))
        } else {
          D[5] <- which(resp[3] == sort(unique(indata[[resp[1]]])))
        }
      }
      names(D)[5] <- "ref.cat"
      resp <- resp[1]
    }
    if (D[1] == "Mixed")
      D <- as.list(D)
    if (D[[1]] == "Mixed") {
      resp <- sub("^c\\(", "", resp)
      resp <- sub("\\)$", "", resp)
      resp <- strsplit(resp, ",")[[1]]
      lenD <- length(D) - 1
      resp2 <- rep(NA, lenD)
      j <- 1
      for (i in 1:length(resp)) {
        if (!(grepl("\\(", resp[i])) && !(grepl("\\)", resp[i]))) {
          resp2[j] <- resp[i]
          j <- j + 1
        } else {
          if (grepl("\\(", resp[i])) {
            ts <- paste(resp[i], ",", sep = "")
          }
          if (grepl("\\)", resp[i])) {
            ts <- paste(ts, resp[i], sep = "")
            resp2[j] <- ts
            j <- j + 1
          }
        }
      }
      resp <- resp2
      for (i in 1:length(resp)) {
        respx <- resp[i]
        if (D[[i + 1]] == "Normal") {
          resp[i] <- respx
        } else if (D[[i + 1]] == "Binomial") {
          D[[i + 1]] <- rep(NA, 3)
          respx <- regmatches(respx, regexec("([[:alnum:]]+)\\(([[:alnum:][:space:],_]+)\\)", respx))[[1]]
          link <- respx[2]
          D[[i + 1]][1] <- "Binomial"
          D[[i + 1]][2] <- link
          respx <- strsplit(respx[3], ",")[[1]]
          D[[i + 1]][3] <- respx[2]
          resp[i] <- respx[1]
        } else if (D[[i + 1]] == "Poisson") {
          respx <- regmatches(respx, regexec("([[:alnum:]]+)\\(([[:alnum:][:space:],_]+)\\)", respx))[[1]]
          link <- respx[2]
          respx <- strsplit(respx[3], ",")[[1]]
          D[[i + 1]] <- rep(NA, 3)
          D[[i + 1]][1] <- "Poisson"
          D[[i + 1]][2] <- link
          if (length(respx) == 2) {
            D[[i + 1]][3] <- respx[2]
          }
          resp[i] <- respx[1]
        }
      }
    }
    
    nleft <- length(left)
    
    categ <- NULL
    leftsplit <- strsplit(left, "(\\+)|(\\|)")
    categstr <- unique(unlist(sapply(leftsplit, function(x) {
      unlist(regmatches(x, gregexpr("([[:alnum:]]*(\\_|\\-|\\^|\\&)*[[:alnum:]])+(\\[+\\]|\\[+[[:print:]]+\\])",
                                    x)))
    })))
    
    left <- sapply(regmatches(left, gregexpr("\\[[^]]*\\]", left), invert = TRUE), function(x) paste(x, collapse = ""))
    
    ncategstr <- length(categstr)
    if (ncategstr > 0) {
      categ <- matrix(, nrow = 3, ncol = ncategstr)
      rownames(categ) <- c("var", "ref", "ncateg")
      for (i in 1:ncategstr) {
        cvx <- unlist(strsplit(categstr[i], "\\["))
        categ[1, i] <- cvx[1]
        cvy <- sub("\\]", "", cvx[2])
        if (cvy == "")
          cvy <- NA
        categ[2, i] <- cvy
        categ[3, i] <- length(levels(indata[[categ[1, i]]]))
      }
    }
    
    fixs.no <- grep("0s+\\|", left)
    fixs <- left[fixs.no]
    fixc.no <- grep("0c+\\|", left)
    fixc <- left[fixc.no]
    if (length(fixc) != 0) {
      fixc <- unlist(strsplit(fixc, "\\|"))
      fixc <- unlist(strsplit(fixc[2], "\\+"))
      cidmat <- matrix(, nrow = length(fixc), ncol = 2)
      for (i in 1:length(fixc)) {
        if (length(grep("\\{", fixc[i])) == 0) {
          cidmat[i, 1] <- fixc[i]
        } else {
          tt <- unlist(strsplit(fixc[i], "\\{"))
          cidmat[i, 1] <- tt[1]
          cidmat[i, 2] <- sub("\\}", "", tt[2])
        }
      }
      fixc <- cidmat[, 1]
    } else {
      cidmat <- NULL
    }
    if (length(fixs) != 0) {
      fixs <- unlist(strsplit(fixs, "\\|"))
      fixs <- unlist(strsplit(fixs[2], "\\+"))
    }
    rands.no <- rep(NA, nlev)
    randc.no <- rep(NA, nlev)
    effect.lev <- nlev:1
    for (i in 1:(nlev)) {
      t1 <- grep(paste(effect.lev[i], "s+\\|", sep = ""), left)
      if (length(t1) != 0)
        rands.no[i] <- t1
      t2 <- grep(paste(effect.lev[i], "c+\\|", sep = ""), left)
      if (length(t2) != 0)
        randc.no[i] <- t2
    }
    randS <- left[rands.no]
    randC <- left[randc.no]
    rands <- randc <- randcpos <- list()
    
    if (nlev > 0) {
      for (i in 1:(nlev)) {
        if (length(randS[i]) != 0) {
          rands[[i]] <- unlist(strsplit(randS[[i]], "\\|"))
          rands[[i]] <- unlist(strsplit(rands[[i]][2], "\\+"))
        }
        if (length(randC[i]) != 0) {
          randc[[i]] <- unlist(strsplit(randC[[i]], "\\|"))
          randc[[i]] <- unlist(strsplit(randc[[i]][2], "\\+"))
          randcpos[[i]] <- rep(NA, length(randc[[i]]))
          if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
            nresp <- length(unique(indata[[resp]]))
          } else {
            nresp <- length(resp)
          }
          for (j in 1:length(randc[[i]])) {
            if (!is.na(randc[[i]][1])) {
              if (length(grep("\\{", randc[[i]][j])) == 0) {
                cidmat <- rbind(cidmat, c(randc[[i]][j], NA))
              } else {
                randcc <- unlist(strsplit(randc[[i]][j], "\\{"))
                randc[[i]][j] <- randcc[1]
                tempid <- sub("\\}", "", randcc[2])
                randcpos[[i]][j] <- tempid
                cidmat <- rbind(cidmat, c(randc[[i]][j], tempid))
              }
            }
          }
        }
      }
      randS <- unique(na.omit(unlist(rands)))
      for (i in 1:length(randS)) {
        if (sum(grepl("\\.", randS[i])) > 0) {
          ttemp <- unlist(strsplit(randS[i], "\\."))
          ttemp <- paste(ttemp[-length(ttemp)], collapse = "")
          if (ttemp %in% randS) {
            randS <- randS[-i]
          }
        }
      }
      if (length(fixs) == 0) {
        temps <- randS
        nonfps <- NULL
        if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
          if (is.factor(indata[[resp]])) {
            names.resp <- levels(indata[[resp]])
          } else {
            names.resp <- as.character(sort(unique(indata[[resp]])))
          }
          names.resp <- names.resp[-as.numeric(D["ref.cat"])]
          refcatint <- as.numeric(D["ref.cat"])
          if (!is.na(refcatint)) {
            if (refcatint == 1)
              usign <- ">" else usign <- "<"
          } else {
            usign <- ""
          }
        } else {
          names.resp <- resp
        }
        for (i in 1:length(temps)) {
          if (D[1] == "Ordered Multinomial") {
            if (sum(grepl("\\.", temps[i])) == 0) {
              nonfps <- c(nonfps, sapply(names.resp, function(x) paste(temps[i], ".(", usign, "=", x, ")",
                                                                       sep = "")))
            } else {
              ttemp <- unlist(strsplit(temps[i], "\\."))
              ttemp <- ttemp[length(ttemp)]
              if (ttemp %in% sapply(names.resp, function(x) paste(".(", usign, "=", x, ")", sep = ""))) {
                nonfps <- c(nonfps, temps[i])
              } else {
                nonfps <- c(nonfps, sapply(names.resp, function(x) paste(temps[i], ".(", usign, "=", x, ")",
                                                                         sep = "")))
              }
            }
          } else {
            if (sum(grepl("\\.", temps[i])) == 0) {
              nonfps <- c(nonfps, sapply(names.resp, function(x) paste(temps[i], x, sep = ".")))
            } else {
              ttemp <- unlist(strsplit(temps[i], "\\."))
              ttemp <- ttemp[length(ttemp)]
              if (ttemp %in% names.resp) {
                nonfps <- c(nonfps, temps[i])
              } else {
                nonfps <- c(nonfps, sapply(names.resp, function(x) paste(temps[i], x, sep = ".")))
              }
            }
          }
        }
        fixs <- temps
      } else {
        temps <- randS
        nonfps <- NULL
        for (i in 1:length(temps)) {
          if (sum(grepl("\\.", temps[i])) > 0) {
            ttemp0 <- unlist(strsplit(temps[i], "\\."))
            ttemp0 <- paste(ttemp0[-length(ttemp0)], collapse = "")
            temps[i] <- ttemp0
          }
        }
        temps <- randS[!(temps %in% fixs)]
        if (length(temps) != 0) {
          nonfps <- NULL
          if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
            if (is.factor(indata[[resp]])) {
              names.resp <- levels(indata[[resp]])
            } else {
              names.resp <- as.character(sort(unique(indata[[resp]])))
            }
            names.resp <- names.resp[-as.numeric(D["ref.cat"])]
            refcatint <- as.numeric(D["ref.cat"])
            if (!is.na(refcatint)) {
              if (refcatint == 1)
                usign <- ">" else usign <- "<"
            } else {
              usign <- ""
            }
          } else {
            names.resp <- resp
          }
          for (i in 1:length(temps)) {
            if (D[1] == "Ordered Multinomial") {
              if (sum(grepl("\\.", temps[i])) == 0) {
                nonfps <- c(nonfps, sapply(names.resp, function(x) paste(temps[i], ".(", usign, "=", x, ")",
                                                                         sep = "")))
              } else {
                ttemp <- unlist(strsplit(temps[i], "\\."))
                ttemp <- ttemp[length(ttemp)]
                if (ttemp %in% sapply(names.resp, function(x) paste(".(", usign, "=", x, ")", sep = ""))) {
                  nonfps <- c(nonfps, temps[i])
                } else {
                  nonfps <- c(nonfps, sapply(names.resp, function(x) paste(temps[i], ".(", usign, "=", x,
                                                                           ")", sep = "")))
                }
              }
            } else {
              if (sum(grepl("\\.", temps[i])) == 0) {
                nonfps <- c(nonfps, sapply(names.resp, function(x) paste(temps[i], x, sep = ".")))
              } else {
                ttemp <- unlist(strsplit(temps[i], "\\."))
                ttemp <- ttemp[length(ttemp)]
                if (ttemp %in% names.resp) {
                  nonfps <- c(nonfps, temps[i])
                } else {
                  nonfps <- c(nonfps, sapply(names.resp, function(x) paste(temps[i], x, sep = ".")))
                }
              }
            }
          }
          fixs <- c(fixs, temps)
        } else {
          nonfps <- character(0)
        }
      }
      
      randC <- na.omit(unlist(randc))
      randCC <- rep(NA, length(randC))
      
      if (length(randC) != 0) {
        na.pos <- which(is.na(cidmat[, 2]))
        non.na.pos <- which(!is.na(cidmat[, 2]))
        for (i in na.pos) {
          if (cidmat[i, 1] %in% cidmat[non.na.pos, 1]) {
            cidmat[i, ] <- c(NA, NA)
          } else {
            tempid <- 1:nresp
            if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
              tempid <- tempid[-as.numeric(D["ref.cat"])]
            }
            cidmat[i, 2] <- paste(tempid, collapse = ",")
          }
        }
        if (sum(is.na(cidmat[, 1])) > 0)
          cidmat <- cidmat[-which(is.na(cidmat[, 1])), ]
        common.coeff <- unique(paste(cidmat[, 1], cidmat[, 2], sep = "@"))
        lencom <- length(common.coeff)
        tt.id <- unlist(strsplit(common.coeff, "\\@"))[(1:lencom) * 2]
        tt.names <- unlist(strsplit(common.coeff, "\\@"))[(1:lencom) * 2 - 1]
        common.coeff <- sub("\\@", "\\.", common.coeff)
        if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial")
          lencol <- nresp - 1 else lencol <- nresp
        ccid.mat <- matrix(0, nrow = length(tt.names), ncol = lencol)
        rownames(ccid.mat) <- tt.names
        
        for (i in 1:length(tt.id)) {
          nonrefcatpos <- as.numeric(unlist(strsplit(tt.id[i], ",")))
          if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
            refcatint <- as.numeric(D["ref.cat"])
            if (!is.na(refcatint)) {
              nonrefcatpos[which(nonrefcatpos > refcatint)] <- nonrefcatpos[which(nonrefcatpos > refcatint)] -
                1
            }
          }
          ccid.mat[i, nonrefcatpos] <- 1
        }
        
        randC <- unique(na.omit(unlist(randc)))
        randCC <- rep(NA, length(randC))
        for (i in 1:length(randC)) {
          randCC[i] <- grep(randC[i], common.coeff)
        }
        randCC <- common.coeff[randCC]
        randCC <- gsub(",", "", randCC)
        if (length(fixc) == 0) {
          nonfpc <- randCC
        } else {
          nonfpc <- randCC[!(randC %in% fixc)]
        }
        fixc <- tt.names
      } else {
        ccid.mat <- NULL
        nonfpc <- NULL
      }
      lenfixc <- length(fixc)
      
      if (lenfixc != 0 && length(randC) == 0) {
        if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial")
          lencol <- nresp - 1 else lencol <- nresp
        ccid.mat <- matrix(0, nrow = lenfixc, ncol = lencol)
        rownames(ccid.mat) <- fixc
        
        for (i in 1:lenfixc) {
          if (is.na(cidmat[i, 2])) {
            tempid <- 1:nresp
            if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
              tempid <- tempid[-as.numeric(D["ref.cat"])]
            }
            cidmat[i, 2] <- paste(tempid, collapse = ",")
          }
          
          nonrefcatpos <- as.numeric(unlist(strsplit(cidmat[i, 2], ",")))
          if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
            refcatint <- as.numeric(D["ref.cat"])
            if (!is.na(refcatint)) {
              nonrefcatpos[which(nonrefcatpos > refcatint)] <- nonrefcatpos[which(nonrefcatpos > refcatint)] -
                1
            }
          }
          ccid.mat[i, nonrefcatpos] <- 1
        }
      }
      
      rp <- list()
      rp.names <- NULL
      
      if (D[1] == "Mixed") {
        for (j in 2:length(D)) {
          if (D[[j]][1] == "Binomial" | D[[j]][1] == "Poisson") {
            rp[["rp1"]] <- c(rp[["rp1"]], paste0("bcons.", j - 1))
          }
        }
      }
      
      for (i in 1:length(rands)) {
        if (!is.na(rands[[i]][1])) {
          rptemp <- NULL
          for (j in 1:length(rands[[i]])) {
            rp.name <- paste("rp", effect.lev[i], sep = "")
            rp.names <- c(rp.names, rp.name)
            if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
              if (is.factor(indata[[resp]])) {
                names.resp <- levels(indata[[resp]])
              } else {
                names.resp <- as.character(sort(unique(indata[[resp]])))
              }
              names.resp <- names.resp[-as.numeric(D["ref.cat"])]
              if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial")
                refcatint <- as.numeric(D["ref.cat"])
            } else {
              names.resp <- resp
            }
            if (D[1] == "Ordered Multinomial") {
              if (!is.na(refcatint)) {
                if (refcatint == 1)
                  usign <- ">" else usign <- "<"
              } else {
                usign <- ""
              }
              if (sum(grepl("\\.", rands[[i]][j])) == 0) {
                rptemp <- c(rptemp, sapply(names.resp, function(x) paste(rands[[i]][j], ".(", usign, "=",
                                                                         x, ")", sep = "")))
              } else {
                ttemp <- unlist(strsplit(rands[[i]][j], "\\."))
                ttemp <- ttemp[length(ttemp)]
                if (ttemp %in% sapply(names.resp, function(x) paste(".(", usign, "=", x, ")", sep = ""))) {
                  rptemp <- c(rptemp, rands[[i]][j])
                } else {
                  rptemp <- c(rptemp, sapply(names.resp, function(x) paste(rands[[i]][j], ".(", usign, "=",
                                                                           x, ")", sep = "")))
                }
              }
            } else {
              if (sum(grepl("\\.", rands[[i]][j])) == 0) {
                rptemp <- c(rptemp, sapply(names.resp, function(x) paste(rands[[i]][j], x, sep = ".")))
              } else {
                ttemp <- unlist(strsplit(rands[[i]][j], "\\."))
                ttemp <- ttemp[length(ttemp)]
                if (ttemp %in% names.resp) {
                  rptemp <- c(rptemp, rands[[i]][j])
                } else {
                  rptemp <- c(rptemp, sapply(names.resp, function(x) paste(rands[[i]][j], x, sep = ".")))
                }
              }
            }
          }
          rp[[rp.name]] <- c(rp[[rp.name]], rptemp)
        }
      }
      
      for (i in 1:length(randc)) {
        if (!is.na(randc[[i]][1])) {
          rp.name <- paste("rp", effect.lev[i], sep = "")
          if (!(rp.name %in% rp.names)) {
            rp.names <- c(rp.names, rp.name)
          }
          rptemp <- NULL
          for (j in 1:length(randc[[i]])) {
            randcid <- randcpos[[i]][j]
            rptemp <- c(rptemp, paste(randc[[i]][j], gsub(",", "", randcid), sep = "."))
          }
          if (is.null(rp[[rp.name]])) {
            rp[[rp.name]] <- rptemp
          } else {
            rptt <- rp[[rp.name]]
            rptemp <- c(rptt, rptemp)
            rp[[rp.name]] <- rptemp
          }
        }
      }
    } else {
      lenfixc <- length(fixc)
      if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
        nresp <- length(unique(indata[[resp]]))
      } else {
        nresp <- length(resp)
      }
      if (lenfixc != 0) {
        if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial")
          lencol <- nresp - 1 else lencol <- nresp
        ccid.mat <- matrix(0, nrow = lenfixc, ncol = lencol)
        rownames(ccid.mat) <- fixc
        
        for (i in 1:lenfixc) {
          if (is.na(cidmat[i, 2])) {
            tempid <- 1:nresp
            if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
              tempid <- tempid[-as.numeric(D["ref.cat"])]
            }
            cidmat[i, 2] <- paste(tempid, collapse = ",")
          }
          
          nonrefcatpos <- as.numeric(unlist(strsplit(cidmat[i, 2], ",")))
          if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
            refcatint <- as.numeric(D["ref.cat"])
            if (!is.na(refcatint)) {
              nonrefcatpos[which(nonrefcatpos > refcatint)] <- nonrefcatpos[which(nonrefcatpos > refcatint)] -
                1
            }
          }
          ccid.mat[i, nonrefcatpos] <- 1
        }
      } else {
        ccid.mat <- NULL
      }
      rp <- nonfps <- nonfps <- nonfpc <- NULL
    }
    
    if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial")
      D[1] <- "Multinomial"
    invars <- new.env()
    if (length(rp) != 0)
      rp <- rp[order(names(rp), decreasing = TRUE)]
    if (length(randC) == 0 && length(fixc) == 0) {
      invars$resp <- resp
      invars$expl <- fixs
      if (length(rp) != 0)
        invars$rp <- rp
      if (length(nonfps) != 0)
        invars$nonfp <- nonfps
      if (!is.null(categ))
        invars$categ <- categ
      invars$D <- D
    } else {
      invars$resp <- resp
      if (length(fixs) != 0)
        invars$expl$sep.coeff <- fixs else invars$expl$sep.coeff <- NA
      if (length(fixc) != 0)
        invars$expl$common.coeff <- fixc else invars$expl$common.coeff <- NA
      if (length(rp) != 0)
        invars$rp <- rp
      if (length(ccid.mat) != 0)
        invars$expl$common.coeff.id <- ccid.mat
      if (length(nonfps) != 0)
        invars$nonfp$nonfp.sep <- nonfps else invars$nonfp$nonfp.sep <- NA
      if (length(nonfpc) != 0)
        invars$nonfp$nonfp.common <- nonfpc else invars$nonfp$nonfp.common <- NA
      if (!is.null(categ))
        invars$categ <- categ
      invars$D <- D
    }
    invars <- as.list(invars)
    
  }
  if (D[1] == "Binomial") {
    D <- rep(NA, 3)
    names(D) <- c("Distr", "link", "denominator")
    D[1] <- "Binomial"
    nlev <- length(levID)
    
    resp <- regmatches(resp, regexec("([[:alnum:]]+)\\(([[:alnum:][:space:],_]+)\\)", resp))[[1]]
    link <- resp[2]
    D[2] <- link
    resp <- strsplit(resp[3], ",")[[1]]
    D[3] <- resp[2]
    resp <- resp[-2]
    
    nleft <- length(left)
    
    categ <- NULL
    leftsplit <- strsplit(left, "(\\+)|(\\|)")
    categstr <- unique(unlist(sapply(leftsplit, function(x) {
      unlist(regmatches(x, gregexpr("([[:alnum:]]*(\\_|\\-|\\^|\\&)*[[:alnum:]])+(\\[+\\]|\\[+[[:print:]]+\\])",
                                    x)))
    })))
    
    left <- sapply(regmatches(left, gregexpr("\\[[^]]*\\]", left), invert = TRUE), function(x) paste(x, collapse = ""))
    
    ncategstr <- length(categstr)
    if (ncategstr > 0) {
      categ <- matrix(, nrow = 3, ncol = ncategstr)
      rownames(categ) <- c("var", "ref", "ncateg")
      for (i in 1:ncategstr) {
        cvx <- unlist(strsplit(categstr[i], "\\["))
        categ[1, i] <- cvx[1]
        cvy <- sub("\\]", "", cvx[2])
        if (cvy == "")
          cvy <- NA
        categ[2, i] <- cvy
        categ[3, i] <- length(levels(indata[[categ[1, i]]]))
      }
    }
    
    fixs.no <- grep("0+\\|", left)
    fixs <- left[fixs.no]
    if (length(fixs) != 0) {
      fixs <- unlist(strsplit(fixs, "\\|"))
      fixs <- unlist(strsplit(fixs[2], "\\+"))
    }
    
    rp <- list()
    rp[["rp1"]] <- "bcons.1"
    
    effect.lev <- nlev:1
    rands.no <- rep(NA, nlev)
    for (i in 1:nlev) {
      t1 <- grep(paste(effect.lev[i], "+\\|", sep = ""), left)
      if (length(t1) != 0)
        rands.no[i] <- t1
    }
    randS <- left[rands.no]
    rands <- list()
    
    for (i in 1:nlev) {
      if (length(randS[i]) != 0) {
        rands[[i]] <- unlist(strsplit(randS[[i]], "\\|"))
        rands[[i]] <- unlist(strsplit(rands[[i]][2], "\\+"))
      }
    }
    randS <- unique(na.omit(unlist(rands)))
    if (length(fixs) == 0) {
      nonfps <- randS
      fixs <- randS
    } else {
      temps <- randS[!(randS %in% fixs)]
      if (length(temps) != 0) {
        nonfps <- temps
        fixs <- c(fixs, temps)
      } else {
        nonfps <- character(0)
      }
    }
    
    rp.names <- NULL
    for (i in 1:length(rands)) {
      if (!is.na(rands[[i]][1])) {
        rp.name <- paste("rp", effect.lev[i], sep = "")
        rp.names <- c(rp.names, rp.name)
        rptemp <- rp[[rp.name]]
        
        for (j in 1:length(rands[[i]])) {
          rptemp <- c(rptemp, rands[[i]][j])
        }
        rp[[rp.name]] <- rptemp
      }
    }
    
    invars <- new.env()
    invars$resp <- resp
    invars$expl <- fixs
    if (length(rp) != 0)
      rp <- rp[order(names(rp), decreasing = TRUE)]
    if (length(rp) != 0)
      invars$rp <- rp
    if (length(nonfps) != 0)
      invars$nonfp <- nonfps
    invars$D <- D
    if (!is.null(categ))
      invars$categ <- categ
    invars <- as.list(invars)
  }
  
  if (D[1] == "Poisson" || D[1] == "Negbinom") {
    nlev <- length(levID)
    resp <- regmatches(resp, regexec("([[:alnum:]]+)\\(([[:alnum:][:space:],_]+)\\)", resp))[[1]]
    link <- resp[2]
    resp <- strsplit(resp[3], ",")[[1]]
    DD <- D[1]
    D <- rep(NA, 3)
    D[1] <- DD
    D[2] <- link
    if (length(resp) == 2) {
      D[3] <- resp[2]
      resp <- resp[-2]
    }
    
    nleft <- length(left)
    
    categ <- NULL
    leftsplit <- strsplit(left, "(\\+)|(\\|)")
    categstr <- unique(unlist(sapply(leftsplit, function(x) {
      unlist(regmatches(x, gregexpr("([[:alnum:]]*(\\_|\\-|\\^|\\&)*[[:alnum:]])+(\\[+\\]|\\[+[[:print:]]+\\])",
                                    x)))
    })))
    left <- sapply(regmatches(left, gregexpr("\\[[^]]*\\]", left), invert = TRUE), function(x) paste(x, collapse = ""))
    
    ncategstr <- length(categstr)
    if (ncategstr > 0) {
      categ <- matrix(, nrow = 3, ncol = ncategstr)
      rownames(categ) <- c("var", "ref", "ncateg")
      for (i in 1:ncategstr) {
        cvx <- unlist(strsplit(categstr[i], "\\["))
        categ[1, i] <- cvx[1]
        cvy <- sub("\\]", "", cvx[2])
        if (cvy == "")
          cvy <- NA
        categ[2, i] <- cvy
        categ[3, i] <- length(levels(indata[[categ[1, i]]]))
      }
    }
    
    fixs.no <- grep("0+\\|", left)
    fixs <- left[fixs.no]
    if (length(fixs) != 0) {
      fixs <- unlist(strsplit(fixs, "\\|"))
      fixs <- unlist(strsplit(fixs[2], "\\+"))
    }
    rp <- list()
    if (D[1] == "Poisson") {
      rp[["rp1"]] <- "bcons.1"
    }
    if (D[1] == "Negbinom") {
      rp[["rp1"]] <- c("bcons.1", "bcons2.1")
    }
    
    effect.lev <- nlev:1
    rands.no <- rep(NA, nlev)
    for (i in 1:nlev) {
      t1 <- grep(paste(effect.lev[i], "+\\|", sep = ""), left)
      if (length(t1) != 0)
        rands.no[i] <- t1
    }
    randS <- left[rands.no]
    rands <- list()
    
    for (i in 1:nlev) {
      if (length(randS[i]) != 0) {
        rands[[i]] <- unlist(strsplit(randS[[i]], "\\|"))
        rands[[i]] <- unlist(strsplit(rands[[i]][2], "\\+"))
      }
    }
    randS <- unique(na.omit(unlist(rands)))
    if (length(fixs) == 0) {
      nonfps <- randS
      fixs <- randS
    } else {
      temps <- randS[!(randS %in% fixs)]
      if (length(temps) != 0) {
        nonfps <- temps
        fixs <- c(fixs, temps)
      } else {
        nonfps <- character(0)
      }
    }
    
    rp.names <- NULL
    for (i in 1:length(rands)) {
      if (!is.na(rands[[i]][1])) {
        rp.name <- paste("rp", effect.lev[i], sep = "")
        rp.names <- c(rp.names, rp.name)
        rptemp <- rp[[rp.name]]
        
        for (j in 1:length(rands[[i]])) {
          rptemp <- c(rptemp, rands[[i]][j])
        }
        rp[[rp.name]] <- rptemp
      }
    }
    
    invars <- new.env()
    invars$resp <- resp
    invars$expl <- fixs
    if (length(rp) != 0)
      rp <- rp[order(names(rp), decreasing = TRUE)]
    if (length(rp) != 0)
      invars$rp <- rp
    if (length(nonfps) != 0)
      invars$nonfp <- nonfps
    invars$D <- D
    if (!is.null(categ))
      invars$categ <- categ
    invars <- as.list(invars)
  }
  
  if (D[1] == "Normal") {
    nlev <- length(levID)
    nleft <- length(left)
    
    categ <- NULL
    leftsplit <- strsplit(left, "(\\+)|(\\|)")
    categstr <- unique(unlist(sapply(leftsplit, function(x) {
      unlist(regmatches(x, gregexpr("([[:alnum:]]*(\\_|\\-|\\^|\\&)*[[:alnum:]])+(\\[+\\]|\\[+[[:print:]]+\\])",
                                    x)))
    })))
    left <- sapply(regmatches(left, gregexpr("\\[[^]]*\\]", left), invert = TRUE), function(x) paste(x, collapse = ""))
    ncategstr <- length(categstr)
    if (ncategstr > 0) {
      categ <- matrix(, nrow = 3, ncol = ncategstr)
      rownames(categ) <- c("var", "ref", "ncateg")
      for (i in 1:ncategstr) {
        cvx <- unlist(strsplit(categstr[i], "\\["))
        categ[1, i] <- cvx[1]
        cvy <- sub("\\]", "", cvx[2])
        if (cvy == "")
          cvy <- NA
        categ[2, i] <- cvy
        categ[3, i] <- length(levels(indata[[categ[1, i]]]))
      }
    }
    
    fixs.no <- grep("0+\\|", left)
    fixs <- left[fixs.no]
    if (length(fixs) != 0) {
      fixs <- unlist(strsplit(fixs, "\\|"))
      fixs <- unlist(strsplit(fixs[2], "\\+"))
    }
    effect.lev <- nlev:1
    rands.no <- rep(NA, nlev)
    for (i in 1:nlev) {
      t1 <- grep(paste(effect.lev[i], "+\\|", sep = ""), left)
      if (length(t1) != 0)
        rands.no[i] <- t1
      
    }
    randS <- left[rands.no]
    rands <- list()
    
    for (i in 1:nlev) {
      if (length(randS[i]) != 0) {
        rands[[i]] <- unlist(strsplit(randS[[i]], "\\|"))
        rands[[i]] <- unlist(strsplit(rands[[i]][2], "\\+"))
      }
    }
    randS <- unique(na.omit(unlist(rands)))
    if (length(fixs) == 0) {
      nonfps <- randS
      fixs <- randS
    } else {
      temps <- randS[!(randS %in% fixs)]
      if (length(temps) != 0) {
        nonfps <- temps
        fixs <- c(fixs, temps)
      } else {
        nonfps <- character(0)
      }
    }
    
    rp <- list()
    rp.names <- NULL
    for (i in 1:length(rands)) {
      if (!is.na(rands[[i]][1])) {
        rp.name <- paste("rp", effect.lev[i], sep = "")
        rp.names <- c(rp.names, rp.name)
        rptemp <- NULL
        
        for (j in 1:length(rands[[i]])) {
          rptemp <- c(rptemp, rands[[i]][j])
        }
        rp[[rp.name]] <- rptemp
      }
    }
    
    invars <- new.env()
    invars$resp <- resp
    invars$expl <- fixs
    if (length(rp) != 0)
      rp <- rp[order(names(rp), decreasing = TRUE)]
    if (length(rp) != 0)
      invars$rp <- rp
    if (length(nonfps) != 0)
      invars$nonfp <- nonfps
    invars$D <- D
    if (!is.null(categ))
      invars$categ <- categ
    invars <- as.list(invars)
  }
  return(invars)
}
