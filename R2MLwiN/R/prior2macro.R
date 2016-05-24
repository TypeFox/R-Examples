#' Translates informative prior information into a concise MLwiN macro.
#'
#' An internal function which takes an R list object containing informative
#' prior information for a multilevel model and translates it into a concise
#' vector object to be used in an MLwiN macro.
#'
#' @param prior An R list object containing prior information for a multilevel
#' model. See `Details' below.
#' @param formula A formula object; see
#' \code{\link{Formula.translate}} or \code{\link{Formula.translate.compat}}.
#' @param levID A string vector specifying the level ID(s).
#' @param D A character string/vector specifying the distribution of the
#' current model.
#' @param indata A data.frame object containing the variables to be modelled.
#'
#' @details
#' The \code{prior} list can contain the following:
#' \itemize{
#' \item \code{fixe}: For the fixed
#' parameters, if proper normal priors are used for some parameters, a list of
#' vectors of length two is provided, each of which specifies the mean and the
#' standard deviation. If not given, default ('flat' or 'diffuse') priors are
#' used for the parameters.
#' \item \code{fixe.common}: For multivariate normal,
#' multinomial and mixed response models, if common coefficients are added, use
#' \code{fixe.common} rather than \code{fixe}.
#' \item \code{fixe.sep}: If the common
#' coefficients are added, use \code{fixe.sep} for the separate coefficients.
#' \item \code{rp<level number>}: A list object specifying the Wishart or gamma prior for the
#' covariance matrix or scalar variance at the levels specified, e.g. \code{rp1} for
#' level 1, \code{rp2} for level 2, etc. Consists of: (1)
#' \code{estimate} -- an estimate for the true value of the inverse of the
#' covariance matrix; (2) \code{size} -- the number of rows in the covariance
#' matrix. Note that this is a weakly-informative prior and the default prior
#' is used if missing.
#' }
#'
#' @return A long vector is returned in the format of MLwiN macro language. This
#' includes all the specified prior parameters.
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso \code{\link{runMLwiN}}
#'
prior2macro <- function(prior, formula, levID, D, indata) {
  ## translation from prior information to MLwiN macro
  
  # fixed effect
  get.terms <- function(fstr) {
    if (is.na(fstr)) {
      return(fstr)
    }
    anycurly <- any(grepl("\\{|\\}", fstr))
    if (anycurly) {
      fstr <- gsub("\\{", "\\[", fstr)
      fstr <- gsub("\\}", "\\]", fstr)
    }
    form <- as.formula(paste0("~", fstr))
    Terms <- terms(form, keep.order = TRUE)
    xterm <- attr(Terms, "term.labels")
    fstr <- unlist(strsplit(fstr, "\\+"))
    fstr <- gsub("[[:space:]]", "", fstr)
    if (any(fstr == "1")) {
      xterm <- c("1", xterm)
    }
    if (anycurly) {
      xterm <- gsub("\\[", "\\{", xterm)
      xterm <- gsub("\\]", "\\}", xterm)
    }
    gsub("[[:space:]]", "", xterm)
  }
  
  get.Idata <- function(left, indata) {
    anycurly <- any(grepl("\\{|\\}", left))
    if (anycurly) {
      left <- gsub("\\{", "\\[", left)
      left <- gsub("\\}", "\\]", left)
    }
    svec <- sapply(left, function(x) unlist(strsplit(x, "\\|"))[2])
    nn <- length(svec)
    Iterms <- character(0)
    for (ii in 1:nn) {
      xform <- as.formula(paste0("~", svec[ii]))
      tterms <- terms(xform, keep.order = TRUE)
      ttermsLabs <- unlist(sapply(attr(tterms, "term.labels"), function(x) unlist(strsplit(x, "\\:"))))
      is_Ifunc <- grepl("^[[:alpha:]]{1}[[:alnum:]]*\\({1}[[:print:]]+\\)+", ttermsLabs)
      if (any(is_Ifunc)) {
        Iterms <- c(Iterms, ttermsLabs[is_Ifunc])
      }
    }
    if (length(Iterms) > 0) {
      Iterms <- sapply(regmatches(Iterms, gregexpr("\\[{1}([[:digit:]]|\\,|[[:space:]])*\\]{1}", Iterms), invert = TRUE),
                       function(x) paste(x, collapse = ""))
      tform <- as.formula(paste0("~0+", paste(Iterms, collapse = "+")))
      dataplus <- model.frame(formula = tform, data = indata, na.action = NULL)
      dataplus.names <- names(dataplus)
      if (anycurly) {
        dataplus.names <- gsub("\\[", "\\{", dataplus.names)
        dataplus.names <- gsub("\\]", "\\}", dataplus.names)
      }
      names(dataplus) <- gsub("[[:space:]]", "", dataplus.names)
      indata <- cbind(indata, dataplus)
    }
    indata
  }
  
  get.polydata <- function(left, indata) {
    is_polyfunc <- grepl("(poly|polym)\\([[:print:]]+\\)", left)
    if (!any(is_polyfunc)) {
      return(list(newleft = character(0), indata = indata))
    }
    anycurly <- any(grepl("\\{|\\}", left))
    if (anycurly) {
      left <- gsub("\\{", "\\[", left)
      left <- gsub("\\}", "\\]", left)
    }
    lvec <- sapply(left, function(x) unlist(strsplit(x, "\\|"))[1])
    svec <- sapply(left, function(x) unlist(strsplit(x, "\\|"))[2])
    nn <- length(svec)
    newleft <- rep(NA, nn)
    newsvec <- rep(NA, nn)
    polyterms <- character(0)
    for (ii in 1:nn) {
      ttermsLabs <- get.terms(svec[ii])
      is_polyfunc <- grepl("(poly|polym)\\([[:print:]]+\\)", ttermsLabs)
      if (any(is_polyfunc)) {
        for (jj in 1:length(ttermsLabs)) {
          xlabs <- unlist(strsplit(ttermsLabs[jj], "\\:"))
          pos_polyfunc <- grepl("(poly|polym)\\([[:print:]]+\\)", xlabs)
          if (any(pos_polyfunc)) {
            xployterms <- xlabs[pos_polyfunc]
            xployterms <- sapply(regmatches(xployterms, gregexpr("\\[{1}([[:digit:]]|\\,|[[:space:]])*\\]{1}",
                                                                 xployterms), invert = TRUE), function(x) paste(x, collapse = ""))
            labs.common <- sapply(regmatches(xployterms, gregexpr("\\[{1}([[:digit:]]|\\,|[[:space:]])*\\]{1}",
                                                                  xployterms)), function(x) paste(x, collapse = ""))
            is_labs.common <- !(labs.common == "")
            otherterms <- xlabs[!pos_polyfunc]
            for (kk in 1:length(xployterms)) {
              tform <- as.formula(paste0("~0+", paste(xployterms[kk], collapse = "+")))
              dataplus <- as.data.frame(model.frame(formula = tform, data = indata, na.action = NULL)[[1]])
              dataplus.names <- paste(make.names(names(dataplus)), xployterms[kk], sep = "_")
              dataplus.names <- gsub("[[:space:]]", "", dataplus.names)
              ndp <- sapply(dataplus.names, nchar)
              if (any(ndp > 32)) {
                dataplus.names <- gsub("x\\=", "", dataplus.names)
                dataplus.names <- gsub("degree\\=", "", dataplus.names)
                dataplus.names <- gsub("coefs\\=", "", dataplus.names)
                dataplus.names <- gsub("raw\\=", "", dataplus.names)
                ndp <- sapply(dataplus.names, nchar)
                if (any(ndp > 32)) {
                  dataplus.names <- gsub("[[:punct:]]", "_", dataplus.names)
                  dataplus.names <- abbreviate(dataplus.names, minlength = 31)
                }
              }
              if (anycurly) {
                dataplus.names <- gsub("\\[", "\\{", dataplus.names)
                dataplus.names <- gsub("\\]", "\\}", dataplus.names)
              }
              names(dataplus) <- dataplus.names
              if (is_labs.common[kk]) {
                tmpnames0 <- paste0(names(dataplus), labs.common[kk])
                if (kk == 1) {
                  dataplus.names0 <- tmpnames0
                } else {
                  dataplus.names0 <- paste(dataplus.names0, tmpnames0, sep = ":")
                }
              } else {
                if (kk == 1) {
                  dataplus.names0 <- names(dataplus)
                } else {
                  dataplus.names0 <- paste(dataplus.names0, names(dataplus), sep = ":")
                }
              }
              indata <- cbind(indata, dataplus)
            }
            if (length(otherterms) > 0) {
              otherTerms <- paste(otherterms, collapse = ":")
              ttermsLabs[jj] <- paste(otherTerms, dataplus.names0, sep = ":")
            } else {
              ttermsLabs[jj] <- paste(dataplus.names0, collapse = "+")
            }
          }
        }
      }
      newsvec[ii] <- paste(ttermsLabs, collapse = "+")
      newleft[ii] <- paste(lvec[ii], newsvec[ii], sep = "|")
    }
    if (anycurly) {
      newleft <- gsub("\\[", "\\{", newleft)
      newleft <- gsub("\\]", "\\}", newleft)
    }
    list(newleft = newleft, indata = indata)
  }
  
  get_categstr <- function(left, indata) {
    categstr0 <- NULL
    categstr1 <- NULL
    leftsplit <- strsplit(left, "\\|")
    for (ii in 1:length(leftsplit)) {
      leftsplit[[ii]] <- get.terms(leftsplit[[ii]][2])
    }
    tmpcategstr <- unique(unlist(leftsplit))
    tmpcategstr <- gsub("\\{{1}([[:digit:]]|\\,)*\\}{1}", "", tmpcategstr)
    lfcol <- sapply(indata, is.factor)
    for (ii in 1:length(tmpcategstr)) {
      ttcategstr <- unlist(strsplit(tmpcategstr[ii], "\\:"))
      lttcateg <- ttcategstr %in% names(indata)[lfcol]
      if (any(lttcateg)) {
        categstr0 <- c(categstr0, tmpcategstr[ii])
        categstr1 <- c(categstr1, ttcategstr[lttcateg])
      }
    }
    categstr0 <- unique(categstr0)
    categstr1 <- unique(categstr1)
    ncategstr0 <- length(categstr0)
    ncategstr1 <- length(categstr1)
    
    categstr2 <- vector("list", ncategstr1)
    names(categstr2) <- categstr1
    categstr3 <- vector("list", ncategstr0)
    names(categstr3) <- categstr0
    
    if (ncategstr0 > 0) {
      # extend data
      for (ii in 1:ncategstr1) {
        f.ext <- as.formula(eval(paste("~0+", categstr1[ii])))
        contrMat <- attr(indata[[categstr1[ii]]], "contrasts")
        if (is.null(contrMat)) {
          data.ext <- model.matrix(f.ext, indata)[, -1, drop = FALSE]
          categstr2[[categstr1[ii]]] <- colnames(data.ext)
        } else {
          keeppos <- rowSums(contrMat) > 0
          data.ext <- model.matrix(f.ext, indata)[, keeppos, drop = FALSE]
          categstr2[[categstr1[ii]]] <- colnames(data.ext)
        }
        indata <- cbind(indata, as.data.frame(data.ext))
      }
      for (ii in 1:ncategstr0) {
        ttcategstr <- unlist(strsplit(categstr0[ii], "\\:"))
        if (length(ttcategstr) == 1) {
          categstr3[[categstr0[ii]]] <- categstr2[[ttcategstr]]
        } else {
          ttlist <- vector("list", length(ttcategstr))
          names(ttlist) <- ttcategstr
          for (jj in 1:length(ttcategstr)) {
            if (ttcategstr[jj] %in% categstr1) {
              ttlist[[ttcategstr[jj]]] <- categstr2[[ttcategstr[jj]]]
            } else {
              ttlist[[ttcategstr[jj]]] <- ttcategstr[jj]
            }
          }
          ttcombs <- apply(expand.grid(ttlist), 1, function(x) paste(x, collapse = ":"))
          categstr3[[categstr0[ii]]] <- ttcombs
        }
      }
    }
    list(categstr = categstr3, indata = indata)
  }
  
  left2leftsc <- function(left, nlev, nresp, D) {
    is_cvar <- grepl("(\\+|\\-|\\*|\\/|\\:|\\|){1}\\(*[[:alnum:]]{1}[[:graph:]]*\\[{1}(c\\(|\\)|\\,|\\:|[[:digit:]])*\\]{1}\\)*(\\+|\\-|\\*|\\/|\\:|$)",
                     left)
    if (!any(is_cvar)) {
      return(left)
    }
    cc <- c(0:nlev)
    respid <- 1:nresp
    if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
      refcatint <- as.integer(D["ref.cat"])
      respid <- respid[respid != refcatint]
    }
    left <- sort(left)
    nleft <- length(left)
    newleft <- character(0)
    for (ii in cc) {
      for (jj in 1:nleft) {
        lev_found <- grepl(paste0("^", ii, "\\|"), left[jj])
        cvar_found <- grepl("(\\+|\\-|\\*|\\/|\\:|\\|){1}\\(*[[:alnum:]]{1}[[:graph:]]*\\[{1}(c\\(|\\)|\\,|\\:|[[:digit:]])*\\]{1}\\)*(\\+|\\-|\\*|\\/|\\:|$)",
                            left[jj])
        if (lev_found && cvar_found) {
          leftjj <- sub(paste0("^", ii, "\\|"), "", left[jj])
          leftjj <- get.terms(leftjj)
          is_cvar_last <- grepl("^[[:alnum:]]{1}[[:graph:]]*\\[{1}(c\\(|\\)|\\,|\\:|[[:digit:]])*\\]{1}$", leftjj)
          is_cvar <- grepl("^[[:alnum:]]{1}[[:graph:]]*\\[{1}(c\\(|\\)|\\,|\\:|[[:digit:]])*\\]{1}", leftjj)
          for (uu in 1:length(leftjj)){
            if (!is_cvar_last[uu] && is_cvar[uu]){
              #move brackets to the end of each term
              tempjju1 <- sapply(regmatches(leftjj[uu], gregexpr("\\[{1}(c\\(|\\)|\\,|\\:|[[:digit:]])*\\]{1}",
                                                               leftjj[uu]), invert = TRUE), function(x) paste(x, collapse = ""))
              tempjju2 <- unlist(regmatches(leftjj[uu], gregexpr("\\[{1}(c\\(|\\)|\\,|\\:|[[:digit:]])*\\]{1}",
                                                               leftjj[uu])))[1]
              leftjj[uu] <- paste0(tempjju1, tempjju2, collapse="")
            }
          }
          svar <- leftjj[!is_cvar]
          if (length(svar) > 0) {
            newsvarjj <- paste0(ii, "s|", paste(svar, collapse = "+"))
            newleft <- c(newleft, newsvarjj)
          }
          leftjj2 <- leftjj[is_cvar]
          leftjj2 <- gsub("\\[{1}(c\\(|\\)|\\,|\\:|[[:digit:]])*\\]{1}\\:{1}", ":", leftjj2)
          for (kk in 1:length(leftjj2)) {
            tmpstr <- unlist(strsplit(leftjj2[kk], "\\["))
            tmpstr1 <- tmpstr[1]
            tmpstr2 <- sub("\\]", "", tmpstr[2])
            if (tmpstr2 != "") {
              tmpvec <- eval(parse(text = tmpstr2))
              tmpstr2 <- paste(tmpvec[tmpvec %in% respid], collapse = ",")
            }
            leftjj2[kk] <- paste0(tmpstr1, "[", tmpstr2, "]")
          }
          newcvarjj <- paste0(ii, "c|", paste(leftjj2, collapse = "+"))
          tmpstr2 <- paste(respid, collapse = ",")
          newcvarjj <- gsub("\\[\\]", paste0("[", tmpstr2, "]"), newcvarjj)
          newleft <- c(newleft, newcvarjj)
        } else {
          if (lev_found) {
            leftjj <- sub(paste0("^", ii, "\\|"), "", left[jj])
            newsvarjj <- paste0(ii, "s|", leftjj)
            newleft <- c(newleft, newsvarjj)
          }
        }
      }
    }
    return(newleft)
  }
  
  lefts2leftc <- function(left, nlev, nresp, D) {
    is_cvar <- grepl("^[[:digit:]]{1,2}c\\|", left)
    if (!any(is_cvar)) {
      return(left)
    }
    cc <- c(0:nlev)
    respid <- 1:nresp
    if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
      refcatint <- as.integer(D["ref.cat"])
      respid <- respid[respid != refcatint]
    }
    left <- sort(left)
    nleft <- length(left)
    
    pos_svar <- grep("^[[:digit:]]{1,2}s\\|", left)
    if (length(pos_svar) > 0) {
      slist <- vector("list", length(pos_svar))
      for (ii in 1:length(pos_svar)) {
        leftii <- unlist(strsplit(left[pos_svar[ii]], "\\|"))
        names(slist)[ii] <- leftii[1]
        slist[[ii]] <- get.terms(leftii[2])
      }
    } else {
      return(left)
    }
    slistnames <- names(slist)
    targetnames <- gsub("s", "c", slistnames)
    pos_target <- sapply(targetnames, function(x) grep(paste0("^", x), left))
    is_target_found <- as.logical(sapply(pos_target, length))
    
    newleft <- vector("list", sum(!is_target_found))
    names(newleft) <- targetnames[!is_target_found]
    
    for (ii in cc) {
      cvar_found <- grepl(paste0("^", ii, "c\\|"), left)
      if (any(cvar_found)) {
        leftjj <- sub(paste0("^", ii, "c\\|"), "", left[cvar_found])
        leftjj <- get.terms(leftjj)
        if (length(pos_svar) > 0) {
          for (kk in 1:length(leftjj)) {
            for (jj in 1:length(slist)) {
              ele <- leftjj[kk]
              valid_ele <- grepl("\\[[[:digit:]]{1}\\]$", ele)
              ele <- unlist(strsplit(ele, "\\["))[1]
              is_match <- ele == slist[[jj]]
              if (any(is_match) && valid_ele) {
                slist[[jj]] <- slist[[jj]][!is_match]
                newele <- paste0(ele, "[", respid, "]")
                if (length(pos_target[[jj]]) > 0) {
                  left[pos_target[[jj]]] <- paste(c(left[pos_target[[jj]]], newele), collapse = "+")
                } else {
                  newleft[[targetnames[jj]]] <- paste(c(newleft[[targetnames[jj]]], newele), collapse = "+")
                }
              }
            }
          }
        }
      }
    }
    for (ii in 1:length(slist)) {
      newele <- slist[[ii]]
      if (length(newele) > 0) {
        left[pos_svar[ii]] <- paste0(slistnames[ii], "|", paste(newele, collapse = "+"))
      } else {
        left[pos_svar[ii]] <- paste0(targetnames[ii], "|", newleft[[targetnames[ii]]])
      }
    }
    if (length(newleft) > 0) {
      for (ii in 1:length(newleft)) {
        xcvar_found <- any(grepl(paste0("^", names(newleft)[ii], "\\|"), left))
        xele <- newleft[[names(newleft)[ii]]]
        if (!xcvar_found && length(xele) > 0) {
          xele.add <- paste0(names(newleft)[ii], "|", xele)
          left <- c(left, xele.add)
        }
      }
    }
    is_empty <- grepl("^[[:digit:]]+(c|s){1}\\|$", left)
    left[!is_empty]
  }
  
  nlev <- length(levID)
  cc <- c(0:nlev)
  Formula <- formula
  is_str_form <- is.character(Formula)
  if (is_str_form) {
    Formula <- gsub("\\{", "\\(", Formula)
    Formula <- gsub("\\}", "\\)", Formula)
    Formula <- gsub("[[:space:]]", "", Formula)
    if (sum(grepl("\\({1}[[:digit:]]+\\|{2}", Formula)) > 0) {
      for (i in cc) {
        Formula <- sub(paste(i, "\\|{2}", sep = ""), paste("\\`", i, "c`\\|", sep = ""), Formula)
        Formula <- sub(paste(i, "\\|", sep = ""), paste("\\`", i, "s`\\|", sep = ""), Formula)
      }
    }
    if (sum(grepl("\\({1}[[:digit:]]+[[:alpha:]]{1}\\|", Formula)) > 0) {
      for (i in cc) {
        Formula <- sub(paste(i, "s\\|", sep = ""), paste("\\`", i, "s`\\|", sep = ""), Formula)
        Formula <- sub(paste(i, "c\\|", sep = ""), paste("\\`", i, "c`\\|", sep = ""), Formula)
      }
    }
    Formula <- as.formula(Formula)
  }
  
  tempfstr <- as.character(Formula)[3]
  tempfstr <- unlist(strsplit(tempfstr, "\\+"))
  tempfstr <- gsub("[[:space:]]", "", tempfstr)
  
  if (!any(D %in% c("Normal", "Multivariate Normal"))) {
    Formula <- update(Formula, ~. + (0 | l1id))
  }
  
  Terms <- terms.formula(Formula, keep.order = TRUE)
  resp <- rownames(attr(Terms, "factors"))[attr(Terms, "response")]
  resp <- gsub("[[:space:]]", "", resp)
  left <- attr(Terms, "term.labels")
  if (is_str_form) {
    left <- gsub("\\(", "\\{", left)
    left <- gsub("\\)", "\\}", left)
  }
  left <- gsub("[[:space:]]", "", left)
  
  if (any(tempfstr == "1")) {
    # if(!all(grepl('\\|', left)) && as.logical(attr(Terms,'intercept'))){
    left <- c("1", left)
  }
  
  if (is.null(levID)) {
    charposlevID <- grepl("\\|{1,2}[[:alpha:]]{1}[[:graph:]]*$", left)
    vlpos <- grepl("\\|", left)
    nonzeropos <- !grepl("\\|{1,2}0{1}(s|c)*$", left)
    vlpos <- vlpos & nonzeropos
    if (any(charposlevID) && (sum(charposlevID) == sum(vlpos))) {
      levID <- sub("[[:graph:]]+\\|{1,2}", "", left[which(vlpos)])
      nlev <- length(levID)
      cc <- c(0:nlev)
      for (ii in 1:nlev) {
        pos_first <- grep(paste0("\\|{1,2}", levID[ii], "$"), left)
        if (length(pos_first) > 0) {
          kk <- 0
          for (jj in 1:length(pos_first)){
            if (grepl("\\|{2}",left[pos_first[jj]])){
              left[pos_first[jj]] <- sub(paste0("\\|{2}", levID[ii], "$"), "", left[pos_first[jj]])
              left[pos_first[jj]] <- paste0(paste0(c(nlev:1)[ii], "||"), left[pos_first[jj]])
            }else{
              if (kk == 0){
                left[pos_first[jj]] <- sub(paste0("\\|{1}", levID[ii], "$"), "", left[pos_first[jj]])
                left[pos_first[jj]] <- paste0(paste0(c(nlev:1)[ii], "|"), left[pos_first[jj]])
                kk <- kk + 1
              }
            }
          }
        }
      }
      onevlzero <- left == "0|1"
      onevdlzero <- left == "0||1"
      delpos <- onevlzero | onevdlzero
      left <- left[!(delpos)]
    } else {
      stop("levID cannot be determined based on the formula")
    }
  } else {
    charposlevID <- grepl("\\|{1,2}[[:alpha:]]{1}[[:graph:]]*$", left)
    if (any(charposlevID)) {
      for (ii in 1:nlev) {
        pos_first <- grep(paste0("\\|{1,2}", levID[ii], "$"), left)
        if (length(pos_first) > 0) {
          kk <- 0
          for (jj in 1:length(pos_first)){
            if (grepl("\\|{2}",left[pos_first[jj]])){
              left[pos_first[jj]] <- sub(paste0("\\|{2}", levID[ii], "$"), "", left[pos_first[jj]])
              left[pos_first[jj]] <- paste0(paste0(c(nlev:1)[ii], "||"), left[pos_first[jj]])
            }else{
              if (kk == 0){
                left[pos_first[jj]] <- sub(paste0("\\|{1}", levID[ii], "$"), "", left[pos_first[jj]])
                left[pos_first[jj]] <- paste0(paste0(c(nlev:1)[ii], "|"), left[pos_first[jj]])
                kk <- kk + 1
              }
            }
          }
        }
      }
      onevlzero <- left == "0|1"
      onevdlzero <- left == "0||1"
      delpos <- onevlzero | onevdlzero
      left <- left[!(delpos)]
    }
  }
  if (sum(grepl("^[[:digit:]]+\\|{2}", left)) > 0) {
    for (i in cc) {
      left <- sub(paste(i, "\\|{2}", sep = ""), paste("\\`", i, "c`\\|", sep = ""), left)
      left <- sub(paste(i, "\\|", sep = ""), paste("\\`", i, "s`\\|", sep = ""), left)
    }
  }
  if (sum(grepl("^\\`{1}[[:digit:]]+[[:alpha:]]{1}\\`{1}\\|", left)) > 0) {
    for (i in cc) {
      left <- sub(paste("\\`", i, "s`\\|", sep = ""), paste(i, "s\\|", sep = ""), left)
      left <- sub(paste("\\`", i, "c`\\|", sep = ""), paste(i, "c\\|", sep = ""), left)
    }
  }
  non0pos <- !grepl("\\|", left)
  if (sum(non0pos) > 0) {
    pos0s <- grepl("^0s\\||0\\|", left)
    if (sum(pos0s) == 1) {
      left[pos0s] <- paste(c(left[pos0s], left[non0pos]), collapse = "+")
      left <- left[!(non0pos)]
    }
    if (sum(pos0s) == 0) {
      mergeterm <- paste0("0|", paste(left[non0pos], collapse = "+"))
      left <- left[!(non0pos)]
      left <- c(left, mergeterm)
    }
    if (sum(pos0s) > 1) {
      stop("allow a 0s/0 term in the formula only")
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
    if (is_str_form) {
      if (resp[3] %in% indata[[resp[1]]]) {
        if (is.factor(indata[[resp[1]]])) {
          D[5] <- which(resp[3] == levels(indata[[resp[1]]]))
        } else {
          D[5] <- which(resp[3] == sort(unique(indata[[resp[1]]])))
        }
      }
    } else {
      if (length(resp) == 3) {
        D[5] <- as.integer(resp[3])
      } else {
        if (length(resp) == 2) {
          # take the first category as base
          D[5] <- 1
        }
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
    if (is_str_form) {
      if (resp[3] %in% indata[[resp[1]]]) {
        if (is.factor(indata[[resp[1]]])) {
          D[5] <- which(resp[3] == levels(indata[[resp[1]]]))
        } else {
          D[5] <- which(resp[3] == sort(unique(indata[[resp[1]]])))
        }
      }
    } else {
      if (length(resp) == 3) {
        D[5] <- as.integer(resp[3])
      } else {
        if (length(resp) == 2) {
          # take the first category as base
          D[5] <- 1
        }
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
  
  if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") {
    names.resp <- levels(indata[[resp]])
    names.resp <- names.resp[-as.numeric(D["ref.cat"])]
    nresp <- length(unique(indata[[resp]]))
  } else {
    names.resp <- resp
    nresp <- length(resp)
  }
  
  if (!is_str_form) {
    left <- left2leftsc(left, nlev, nresp, D)
    left <- lefts2leftc(left, nlev, nresp, D)
    left <- gsub("\\[", "\\{", left)
    left <- gsub("\\]", "\\}", left)
  }
  cflag <- 0
  if (sum(grepl("\\(+[[:digit:]]+[[:alpha:]]+\\|", left)) > 0)
    cflag <- 1
  if (cflag == 0) {
    for (i in cc) {
      left <- sub(paste(i, "\\|", sep = ""), paste(i, "s\\|", sep = ""), left)
    }
  }
  nleft <- length(left)
  
  categ <- NULL
  if (is_str_form) {
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
  } else {
    indata <- get.Idata(left, indata)
    xpoly <- get.polydata(left, indata)
    if (length(xpoly$newleft) > 0) {
      left <- xpoly$newleft
      indata <- xpoly$indata
    }
    rm(xpoly)
    categobj <- get_categstr(left, indata)
    categstr3 <- categobj$categstr
    categstr0 <- names(categstr3)
  }
  TT <- NULL
  
  if (cflag == 1) {
    fixs.no <- grep("0s\\|", left)
    fixs <- left[fixs.no]
    if (length(fixs) != 0) {
      fixs <- unlist(strsplit(fixs, "\\|"))
      fixs <- get.terms(fixs[2])
    }
    if (length(categstr3) != 0) {
      if (length(fixs) != 0) {
        delpos <- numeric(0)
        addterms <- numeric(0)
        for (ii in 1:length(fixs)) {
          replacepos <- fixs[ii] %in% categstr0
          if (replacepos) {
            taddterms <- categstr3[[fixs[ii]]]
            delpos <- c(delpos, ii)
            addterms <- c(addterms, taddterms)
          }
        }
        if (length(delpos) > 0) {
          fixs <- fixs[-delpos]
          fixs <- c(fixs, addterms)
        }
      }
    }
    for (i in 1:length(fixs)) {
      if (fixs[i] %in% categ[1, ]) {
        pos <- which(fixs[i] == categ[1, ])
        fixsa <- fixs[1:i]
        fixsa <- fixsa[-i]
        fixsb <- fixs[i:length(fixs)]
        fixsb <- fixsb[-1]
        if (is.na(categ[2, pos])) {
          fixs <- c(fixsa, levels(indata[[categ[1, pos]]]), fixsb)
        } else {
          refx <- categ[2, pos]
          categx <- levels(indata[[categ[1, pos]]])
          categx <- categx[!(refx == categx)]
          fixs <- c(fixsa, categx, fixsb)
        }
      }
    }
    fixS <- NULL
    for (i in 1:length(names.resp)) {
      fixS <- c(fixS, paste(fixs, names.resp[i], sep = "."))
    }
    fixs <- fixS
    fixc.no <- grep("0c\\|", left)
    fixc <- left[fixc.no]
    if (length(fixc) != 0) {
      fixc <- unlist(strsplit(fixc, "\\|"))
      fixc <- get.terms(fixc[2])
      fixcc <- rep(NA, nrow = length(fixc))
      for (i in 1:length(fixc)) {
        if (length(grep("\\{", fixc[i])) == 0) {
          fixcc[i] <- fixc[i]
        } else {
          tt <- unlist(strsplit(fixc[i], "\\{"))
          fixcc <- tt[1]
        }
      }
      fixc <- fixcc
    }
    if (length(categstr3) != 0) {
      if (length(fixc) != 0) {
        delpos <- numeric(0)
        addterms <- numeric(0)
        for (ii in 1:length(fixc)) {
          replacepos <- fixc[ii] %in% categstr0
          if (replacepos) {
            taddterms <- categstr3[[fixc[ii]]]
            delpos <- c(delpos, ii)
            addterms <- c(addterms, taddterms)
          }
        }
        if (length(delpos) > 0) {
          fixc <- fixc[-delpos]
          fixc <- c(fixc, addterms)
        }
      }
    }
    if (length(fixs) > 0) {
      fpps <- prior$fixe.sep
      for (i in 1:length(fixs)) {
        if (fixs[i] %in% names(fpps)) {
          tname <- fixs[i]
          TT <- c(TT, 1, fpps[[tname]][1], fpps[[tname]][2])
        } else {
          TT <- c(TT, 0)
        }
      }
    }
    if (length(fixc) > 0) {
      fppc <- prior$fixe.common
      for (i in 1:length(fixc)) {
        if (fixc[i] %in% names(fppc)) {
          tname <- fixc[i]
          TT <- c(TT, 1, fpps[[tname]][1], fpps[[tname]][2])
        } else {
          TT <- c(TT, 0)
        }
      }
    }
  } else {
    fixs.no <- grep("0s\\|", left)
    fixs <- left[fixs.no]
    if (length(fixs) != 0) {
      fixs <- unlist(strsplit(fixs, "\\|"))
      fixs <- get.terms(fixs[2])
    }
    if (length(categstr3) != 0) {
      if (length(fixs) != 0) {
        delpos <- numeric(0)
        addterms <- numeric(0)
        for (ii in 1:length(fixs)) {
          replacepos <- fixs[ii] %in% categstr0
          if (replacepos) {
            taddterms <- categstr3[[fixs[ii]]]
            delpos <- c(delpos, ii)
            addterms <- c(addterms, taddterms)
          }
        }
        if (length(delpos) > 0) {
          fixs <- fixs[-delpos]
          fixs <- c(fixs, addterms)
        }
      }
    }
    for (i in 1:length(fixs)) {
      if (fixs[i] %in% categ[1, ]) {
        pos <- which(fixs[i] == categ[1, ])
        fixsa <- fixs[1:i]
        fixsa <- fixsa[-i]
        fixsb <- fixs[i:length(fixs)]
        fixsb <- fixsb[-1]
        if (is.na(categ[2, pos])) {
          fixs <- c(fixsa, levels(indata[[categ[1, pos]]]), fixsb)
        } else {
          refx <- categ[2, pos]
          categx <- levels(indata[[categ[1, pos]]])
          categx <- categx[!(refx == categx)]
          fixs <- c(fixsa, categx, fixsb)
        }
      }
    }
    if (length(names.resp) > 1) {
      fixS <- NULL
      for (i in 1:length(names.resp)) {
        fixS <- c(fixS, paste(fixs, names.resp[i], sep = "."))
      }
      fixs <- fixS
    }
    if (length(fixs) > 0) {
      fpps <- prior$fixe
      for (i in 1:length(fixs)) {
        if (fixs[i] %in% names(fpps)) {
          tname <- fixs[i]
          TT <- c(TT, 1, fpps[[tname]][1], fpps[[tname]][2])
        } else {
          TT <- c(TT, 0)
        }
      }
    }
  }
  
  if (D[1] == "Normal" || D[1] == "Multivariate Normal" || D[1] == "Mixed") {
    efflev <- nlev:1
  } else {
    if (nlev > 1)
      efflev <- nlev:2 else efflev <- NULL
  }
  
  if (!is.null(efflev)) {
    rp.names <- paste("rp", efflev, sep = "")
    for (i in 1:length(rp.names)) {
      if (rp.names[i] %in% names(prior)) {
        tname <- rp.names[i]
        mat <- prior[[tname]]$estimate
        mat[upper.tri(mat)] <- mat[lower.tri(mat)]
        tt <- c(as.vector(mat[upper.tri(mat, diag = TRUE)]), prior[[tname]]$size)
        TT <- c(TT, 1, tt)
      } else {
        TT <- c(TT, 0)
      }
    }
  }
  if (!(D[1] == "Normal" || D[1] == "Multivariate Normal" || D[1] == "Mixed"))
    TT <- c(TT, 0)
  TT
}
