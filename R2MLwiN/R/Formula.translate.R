#' An internal function to translate an R formula into an R list object.
#'
#' A model formula, as a formula object written in R-type syntax, is translated
#' into an R list object.
#'
#' @param Formula A \code{\link[stats]{formula}} object specifying a multilevel model.
#' See \code{\link[stats]{formula}} for notes on general usage, and \code{\link{runMLwiN}} for
#' further details.
#' @param D A character string/vector specifying the type of distribution to be modelled, which
#' can include \code{'Normal'} (the default), \code{'Binomial'}, \code{'Poisson'},
#' \code{'Negbinom'}, \code{'Unordered Multinomial'}, \code{'Ordered Multinomial'},
#' \code{'Multivariate Normal'}, or \code{'Mixed'}. In the case of the latter,
#' \code{'Mixed'} precedes the response types which also need to be listed in
#' \code{D}, e.g. \code{c('Mixed', 'Normal', 'Binomial')}; these need to be
#' be listed in the same order to which they are referred to in the
#' \code{Formula} object.
#' @param indata A data.frame object containing the data to be modelled.
#' Optional (can \code{\link[base]{attach}} as an alternative) but recommended.
#'
#' @return Outputs an R list object, which is then used as the input for
#' \code{\link{write.IGLS}} or \code{\link{write.MCMC}}.
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#' @seealso
#' \code{\link{runMLwiN}}, \code{\link{write.IGLS}}, \code{\link{write.MCMC}}; for
#' function allowing back-compatibility with Formula syntax used in older
#' versions of \pkg{R2MLwiN} (<0.8.0) see \code{\link{Formula.translate.compat}}.
#'
#' @examples
#' \dontrun{
#' # NB: See demo(packge = 'R2MLwiN') for a wider range of examples.
#' library(R2MLwiN)
#' # NOTE: if MLwiN not saved in location R2MLwiN defaults to, specify path via:
#' # options(MLwiN_path = 'path/to/MLwiN vX.XX/')
#' # If using R2MLwiN via WINE, the path may look like this:
#' # options(MLwiN_path = '/home/USERNAME/.wine/drive_c/Program Files (x86)/MLwiN vX.XX/')
#'
#' # Two-level random intercept model with student (level 1) nested within
#' # school (level 2) and standlrt added to the fixed part.
#' # Importantly, the ordering of school and student reflects their hierarchy,
#' # with the highest level (school) specified first.
#' # E.g. see demo(UserGuide04)
#' data(tutorial, package = 'R2MLwiN')
#' (mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student),
#'                      data = tutorial))
#'
#' # Adding a random slope
#' (mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school)
#'                      + (1 | student), data = tutorial))
#'
#' # Exploring complex level 1 variation
#' # E.g. see demo(UserGuide07)
#' (mymodel3 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school)
#'                       + (1 + standlrt | student), data = tutorial))
#'
#' # Logit link with cons specified as denominator
#' # Note level 1 ID not explicitly specified
#' # E.g. see demo(UserGuide09)
#' data(bang, package = 'R2MLwiN')
#' (mymodel4 <- runMLwiN(logit(use, cons) ~ 1 + lc + age + (1 | district),
#'                       D = 'Binomial', data = bang))
#'
#' # Mixed response model
#' # Note using MCMC estimation (EstM = 1)
#' # Normal (english) and Bernoulli (behaviour) distributed responses
#' # probit link modelling behaviour with cons as denominator
#' # E.g. see demo(MCMCGuide19)
#' data(jspmix1, package = 'R2MLwiN')
#' (mymodel <- runMLwiN(c(english, probit(behaviour, cons)) ~
#'                      1 + sex + ravens + fluent[1] + (1 | school) + (1[1] | id),
#'                      D = c('Mixed', 'Normal', 'Binomial'),
#'                      estoptions = list(EstM = 1,
#'                      mcmcMeth = list(fixM = 1, residM = 1, Lev1VarM = 1)),
#'                      data = jspmix1))
#' }
#'
Formula.translate <- function(Formula, D = "Normal", indata) {
  
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
  
  get.offset <- function(Formula, indata) {
    tterms <- terms.formula(Formula, keep.order = TRUE)
    pos_offset <- attr(tterms, "offset")
    if (length(pos_offset) > 0) {
      tlabs <- as.character(attr(tterms, "variables"))[-1]
      left0 <- tlabs[attr(tterms, "offset")]
    } else {
      return(list(offset.label = character(0), indata = indata))
    }
    is_num <- sapply(indata, is.numeric)
    # use the first numeric column as the temporary response
    trespname <- names(indata)[is_num]
    trespname <- trespname[1]
    
    tform <- as.formula(paste0(trespname, "~", paste(left0, collapse = "+")))
    tframe <- model.frame(formula = tform, data = indata, na.action = NULL)
    myoffset <- model.offset(tframe)
    if (is.null(myoffset)) {
      return(list(offset.label = character(0), indata = indata))
    } else {
      indata[["_OFFSET"]] <- myoffset
      return(list(offset.label = "_OFFSET", indata = indata))
    }
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
  
  get.Interdata <- function(left, indata) {
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
      ttermsLabs <- attr(tterms, "term.labels")
      is_inter <- grepl("\\:", ttermsLabs)
      if (any(is_inter)) {
        Iterms <- c(Iterms, ttermsLabs[is_inter])
      }
    }
    if (length(Iterms) > 0) {
      Iterms <- sapply(regmatches(Iterms, gregexpr("\\[{1}([[:digit:]]|\\,|[[:space:]])*\\]{1}", Iterms), invert = TRUE), 
                       function(x) paste(x, collapse = ""))
      tform <- as.formula(paste0("~0+", paste(Iterms, collapse = "+")))
      na_act <- options("na.action")[[1]]
      options(na.action = "na.pass")
      dataplus <- model.matrix(object = tform, data = indata)
      options(na.action = na_act)
      dataplus <- as.data.frame(dataplus)
      dataplus.names <- names(dataplus)
      if (anycurly) {
        dataplus.names <- gsub("\\[", "\\{", dataplus.names)
        dataplus.names <- gsub("\\]", "\\}", dataplus.names)
      }
      names(dataplus) <- gsub("\\.", "\\_", dataplus.names)
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
        na_act <- options("na.action")[[1]]
        options(na.action = "na.pass")
        if (is.null(contrMat)) {
          data.ext <- model.matrix(f.ext, indata)[, -1, drop = FALSE]
        } else {
          keeppos <- rowSums(contrMat) > 0
          data.ext <- model.matrix(f.ext, indata)[, keeppos, drop = FALSE]
        }
        options(na.action = na_act)
        colnames(data.ext) <- gsub("\\.", "\\_", colnames(data.ext))
        categstr2[[categstr1[ii]]] <- colnames(data.ext)
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
  
  nlev <- 0
  cc <- c(0:nlev)
  
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
  left <- gsub("[[:space:]]", "", left)
  
  if (any(tempfstr == "1")) {
    # if(!all(grepl('\\|', left)) && as.logical(attr(Terms,'intercept'))){
    left <- c("1", left)
  }
  
  
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
              #match the first position for this levID
              #levID can have the same name but represent different levels
              left[pos_first[jj]] <- sub(paste0("\\|{1}", levID[ii], "$"), "", left[pos_first[jj]])
              left[pos_first[jj]] <- paste0(paste0(c(nlev:1)[ii], "|"), left[pos_first[jj]])   
              kk <- kk +1     
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
  
  if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial" || D[1] == "Multivariate Normal" || D[1] == 
        "Mixed") {
    nlev <- length(levID)
    cc <- c(0:nlev)
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
      if (length(resp) > 1) {
        D[3] <- resp[2]
      } else {
        # Create vector of ones for the denominator
        D[3] <- "_denom"
        indata[["_denom"]] <- rep(1, nrow(indata))
      }
      names(D)[3] <- "denominator"
      D[4] <- 0
      names(D)[4] <- "mode"
      if (length(resp) > 2) {
        D[5] <- as.integer(resp[3])
      } else {
        # take the first category as base
        D[5] <- 1
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
      if (length(resp) > 1) {
        D[3] <- resp[2]
      } else {
        # Create vector of ones for the denominator
        D[3] <- "_denom"
        indata[["_denom"]] <- rep(1, nrow(indata))
      }
      names(D)[3] <- "denominator"
      D[4] <- 1
      names(D)[4] <- "mode"
      if (length(resp) > 2) {
        D[5] <- as.integer(resp[3])
      } else {
        # take the first category as base
        D[5] <- 1
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
          if (length(respx) > 1) {
            D[[i + 1]][3] <- respx[2]
          } else {
            # Create vector of ones for the denominator
            D[[i + 1]][3] <- "_denom"
            indata[["_denom"]] <- rep(1, nrow(indata))
          }         
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
      nresp <- length(unique(indata[[resp]]))
    } else {
      nresp <- length(resp)
    }
    if (D[[1]] == "Mixed") {
      for (i in 1:nresp) {
        if (D[[i + 1]] == "Poisson" && is.na(D[[i + 1]][3])) {
          myoffset <- get.offset(Formula, indata)
          if (length(myoffset$offset.label) > 0) {
            D[[i + 1]][3] <- myoffset$offset.label
            indata <- myoffset$indata
          }
        }
      }
    }
    
    left <- left2leftsc(left, nlev, nresp, D)
    left <- lefts2leftc(left, nlev, nresp, D)
    left <- gsub("\\[", "\\{", left)
    left <- gsub("\\]", "\\}", left)
    
    cflag <- 0
    if (sum(grepl("\\(+[[:digit:]]+[[:alpha:]]+\\|", left)) > 0) 
      cflag <- 1
    if (cflag == 0) {
      for (i in cc) {
        left <- sub(paste(i, "\\|", sep = ""), paste(i, "s\\|", sep = ""), left)
      }
    }
    nleft <- length(left)
    
    indata <- get.Idata(left, indata)
    indata <- get.Interdata(left, indata)
    xpoly <- get.polydata(left, indata)
    if (length(xpoly$newleft) > 0) {
      left <- xpoly$newleft
      indata <- xpoly$indata
    }
    rm(xpoly)
    categobj <- get_categstr(left, indata)
    categstr3 <- categobj$categstr
    categstr0 <- names(categstr3)
    
    fixs.no <- grep("0s+\\|", left)
    fixs <- left[fixs.no]
    fixc.no <- grep("0c+\\|", left)
    fixc <- left[fixc.no]
    if (length(fixc) != 0) {
      fixc <- unlist(strsplit(fixc, "\\|"))
      fixc <- get.terms(fixc[2])
      cidmat <- matrix(, nrow = length(fixc), ncol = 2)
      for (i in 1:length(fixc)) {
        if (length(grep("\\{", fixc[i])) == 0) {
          cidmat[i, 1] <- fixc[i]
        } else {
          tt <- unlist(strsplit(fixc[i], "\\{"))
          cidmat[i, 1] <- tt[1]
          if (length(grep( "\\}$",tt[2]))==1) {
            cidmat[i, 2] <- sub("\\}", "", tt[2])
          } else{
            ttt <- unlist(strsplit(tt[2], "\\}"))
            cidmat[i, 2] <- ttt[1]
            cidmat[i, 1] <- paste0(cidmat[i, 1], ttt[2])
          }
        }
      }
      fixc <- cidmat[, 1]
      fixcid <- cidmat[, 2]
    } else {
      cidmat <- NULL
    }
    if (length(fixs) != 0) {
      fixs <- unlist(strsplit(fixs, "\\|"))
      fixs <- get.terms(fixs[2])
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
          rands[[i]] <- get.terms(rands[[i]][2])
        }
        if (length(randC[i]) != 0) {
          randc[[i]] <- unlist(strsplit(randC[[i]], "\\|"))
          randc[[i]] <- get.terms(randc[[i]][2])
          randcpos[[i]] <- rep(NA, length(randc[[i]]))
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
        if (length(temps) > 0) {
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
          cidmat <- cidmat[!(is.na(cidmat[, 1])), , drop=FALSE]
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
        if (length(fixc) != 0) {
          fixC <- paste(fixc, fixcid, sep = ".")
          nonfpc <- gsub("\\,", "", common.coeff[!(common.coeff %in% fixC)])
        } else {
          nonfpc <- gsub("\\,", "", common.coeff)
        }
        fixc <- tt.names
        # randCC=integer(0) for (i in 1:length(randC)){ tmpitem <- which(randC[i]==tt.names) randCC <- c(randCC, tmpitem)
        # } randCC=common.coeff[randCC] randCC=gsub(',','',randCC) if (length(fixc)==0){ nonfpc= randCC }else{
        # nonfpc=randCC[!(randC%in%fixc)] }
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
    
    # substitute categorical variables
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
      if (length(nonfps) != 0) {
        delpos <- numeric(0)
        addterms <- numeric(0)
        for (ii in 1:length(nonfps)) {
          replacepos <- nonfps[ii] %in% categstr0
          if (replacepos) {
            taddterms <- categstr3[[nonfps[ii]]]
            delpos <- c(delpos, ii)
            addterms <- c(addterms, taddterms)
          }
        }
        if (length(delpos) > 0) {
          nonfps <- nonfps[-delpos]
          nonfps <- c(nonfps, addterms)
        }
      }
      
      if (length(fixc) != 0) {
        delpos <- numeric(0)
        delpos2 <- numeric(0)
        addterms <- numeric(0)
        addccid <- numeric(0)
        for (ii in 1:length(fixc)) {
          replacepos <- fixc[ii] %in% categstr0
          repeated <- sum(fixc[ii] == fixc) > 1
          if (replacepos) {
            taddterms <- categstr3[[fixc[ii]]]
            if (!repeated) {
              len.taddterms <- length(taddterms)
              taddccid <- matrix(rep(ccid.mat[fixc[ii], ], len.taddterms), nrow = len.taddterms, byrow = TRUE)
              rownames(taddccid) <- taddterms
              addccid <- rbind(addccid, taddccid)
              delpos2 <- c(delpos2, ii)
            }
            delpos <- c(delpos, ii)
            addterms <- c(addterms, taddterms)
          }
        }
        if (length(delpos) > 0) {
          fixc <- fixc[-delpos]
          fixc <- c(fixc, addterms)
        }
        if (length(delpos2) > 0) {
          ccid.mat <- ccid.mat[-delpos2, ]
          ccid.mat <- rbind(ccid.mat, addccid)
        }
      }
      if (length(nonfpc) != 0) {
        delpos <- numeric(0)
        addterms <- numeric(0)
        for (ii in 1:length(nonfpc)) {
          if (grepl("\\.{1}([[:digit:]]+|[[:alpha:]]{1}[[:graph:]]*)$", nonfpc[ii])) {
            tmprp <- sub("\\.{1}([[:digit:]]+|[[:alpha:]]{1}[[:graph:]]*)$", "", nonfpc[ii])
            tmprp.common <- sapply(regmatches(nonfpc[ii], gregexpr("\\.{1}([[:digit:]]+|[[:alpha:]]{1}[[:graph:]]*)$", 
                                                                   nonfpc[ii])), function(x) paste(x, collapse = ""))
            replacepos <- tmprp %in% categstr0
            if (replacepos) {
              taddterms <- paste0(categstr3[[tmprp]], tmprp.common)
              delpos <- c(delpos, ii)
              addterms <- c(addterms, taddterms)
            }
          } else {
            replacepos <- nonfpc[ii] %in% categstr0
            if (replacepos) {
              taddterms <- categstr3[[nonfpc[ii]]]
              delpos <- c(delpos, ii)
              addterms <- c(addterms, taddterms)
            }
          }
        }
        if (length(delpos) > 0) {
          nonfpc <- nonfpc[-delpos]
          nonfpc <- c(nonfpc, addterms)
        }
      }
      
      if (length(rp) != 0) {
        for (ii in 1:length(rp)) {
          delpos <- numeric(0)
          addterms <- numeric(0)
          for (jj in 1:length(rp[[ii]])) {
            if (grepl("\\.{1}([[:digit:]]+|[[:alpha:]]{1}[[:graph:]]*)$", rp[[ii]][jj])) {
              tmprp <- sub("\\.{1}([[:digit:]]+|[[:alpha:]]{1}[[:graph:]]*)$", "", rp[[ii]][jj])
              tmprp.common <- sapply(regmatches(rp[[ii]][jj], gregexpr("\\.{1}([[:digit:]]+|[[:alpha:]]{1}[[:graph:]]*)$", 
                                                                       rp[[ii]][jj])), function(x) paste(x, collapse = ""))
              replacepos <- tmprp %in% categstr0
              if (replacepos) {
                delpos <- c(delpos, jj)
                taddterms <- paste0(categstr3[[tmprp]], tmprp.common)
                addterms <- c(addterms, taddterms)
              }
            } else {
              replacepos <- rp[[ii]][jj] %in% categstr0
              if (replacepos) {
                delpos <- c(delpos, jj)
                taddterms <- categstr3[[rp[[ii]][jj]]]
                addterms <- c(addterms, taddterms)
              }
            }
          }
          if (length(delpos) > 0) {
            rp[[ii]] <- rp[[ii]][-delpos]
            rp[[ii]] <- c(rp[[ii]], addterms)
          }
        }
      }
    }
    
    if (D[1] == "Ordered Multinomial" || D[1] == "Unordered Multinomial") 
      D[1] <- "Multinomial"
    invars <- new.env()
    invars$levID <- levID
    if (length(rp) != 0) 
      rp <- rp[order(names(rp), decreasing = TRUE)]
    if (length(randC) == 0 && length(fixc) == 0) {
      invars$resp <- resp
      invars$expl <- fixs
      if (length(rp) != 0) 
        invars$rp <- rp
      if (length(nonfps) != 0) 
        invars$nonfp <- nonfps
      invars$indata <- categobj$indata
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
      invars$indata <- categobj$indata
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
    if (length(resp) > 1) {
      D[3] <- resp[2]
    } else {
      # Create vector of ones for the denominator
      D[3] <- "_denom"
      indata[["_denom"]] <- rep(1, nrow(indata))
    }
    resp <- resp[-2]
    
    nleft <- length(left)
    
    indata <- get.Idata(left, indata)
    indata <- get.Interdata(left, indata)
    xpoly <- get.polydata(left, indata)
    if (length(xpoly$newleft) > 0) {
      left <- xpoly$newleft
      indata <- xpoly$indata
    }
    rm(xpoly)
    categobj <- get_categstr(left, indata)
    categstr3 <- categobj$categstr
    categstr0 <- names(categstr3)
    
    fixs.no <- grep("0+\\|", left)
    fixs <- left[fixs.no]
    if (length(fixs) != 0) {
      fixs <- unlist(strsplit(fixs, "\\|"))
      fixs <- get.terms(fixs[2])
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
        rands[[i]] <- get.terms(rands[[i]][2])
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
    
    # substitute categorical variables
    
    if (length(categstr3) != 0) {
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
      if (length(nonfps) != 0) {
        delpos <- numeric(0)
        addterms <- numeric(0)
        for (ii in 1:length(nonfps)) {
          replacepos <- nonfps[ii] %in% categstr0
          if (replacepos) {
            taddterms <- categstr3[[nonfps[ii]]]
            delpos <- c(delpos, ii)
            addterms <- c(addterms, taddterms)
          }
        }
        if (length(delpos) > 0) {
          nonfps <- nonfps[-delpos]
          nonfps <- c(nonfps, addterms)
        }
      }
      if (length(rp) != 0) {
        for (ii in 1:length(rp)) {
          delpos <- numeric(0)
          addterms <- numeric(0)
          for (jj in 1:length(rp[[ii]])) {
            replacepos <- rp[[ii]][jj] %in% categstr0
            if (replacepos) {
              delpos <- c(delpos, jj)
              taddterms <- categstr3[[rp[[ii]][jj]]]
              addterms <- c(addterms, taddterms)
            }
          }
          if (length(delpos) > 0) {
            rp[[ii]] <- rp[[ii]][-delpos]
            rp[[ii]] <- c(rp[[ii]], addterms)
          }
        }
      }
    }
    
    invars <- new.env()
    invars$levID <- levID
    invars$resp <- resp
    invars$expl <- fixs
    if (length(rp) != 0) 
      rp <- rp[order(names(rp), decreasing = TRUE)]
    if (length(rp) != 0) 
      invars$rp <- rp
    if (length(nonfps) != 0) 
      invars$nonfp <- nonfps
    invars$D <- D
    invars$indata <- categobj$indata
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
    
    if (is.na(D[3])) {
      myoffset <- get.offset(Formula, indata)
      if (length(myoffset$offset.label) > 0) {
        D[3] <- myoffset$offset.label
        indata <- myoffset$indata
      }
    }
    indata <- get.Idata(left, indata)
    indata <- get.Interdata(left, indata)
    xpoly <- get.polydata(left, indata)
    if (length(xpoly$newleft) > 0) {
      left <- xpoly$newleft
      indata <- xpoly$indata
    }
    rm(xpoly)
    categobj <- get_categstr(left, indata)
    categstr3 <- categobj$categstr
    categstr0 <- names(categstr3)
    
    fixs.no <- grep("0+\\|", left)
    fixs <- left[fixs.no]
    if (length(fixs) != 0) {
      fixs <- unlist(strsplit(fixs, "\\|"))
      fixs <- get.terms(fixs[2])
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
        rands[[i]] <- get.terms(rands[[i]][2])
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
    
    # substitute categorical variables
    
    if (length(categstr3) != 0) {
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
      if (length(nonfps) != 0) {
        delpos <- numeric(0)
        addterms <- numeric(0)
        for (ii in 1:length(nonfps)) {
          replacepos <- nonfps[ii] %in% categstr0
          if (replacepos) {
            taddterms <- categstr3[[nonfps[ii]]]
            delpos <- c(delpos, ii)
            addterms <- c(addterms, taddterms)
          }
        }
        if (length(delpos) > 0) {
          nonfps <- nonfps[-delpos]
          nonfps <- c(nonfps, addterms)
        }
      }
      if (length(rp) != 0) {
        for (ii in 1:length(rp)) {
          delpos <- numeric(0)
          addterms <- numeric(0)
          for (jj in 1:length(rp[[ii]])) {
            replacepos <- rp[[ii]][jj] %in% categstr0
            if (replacepos) {
              delpos <- c(delpos, jj)
              taddterms <- categstr3[[rp[[ii]][jj]]]
              addterms <- c(addterms, taddterms)
            }
          }
          if (length(delpos) > 0) {
            rp[[ii]] <- rp[[ii]][-delpos]
            rp[[ii]] <- c(rp[[ii]], addterms)
          }
        }
      }
    }
    
    invars <- new.env()
    invars$levID <- levID
    invars$resp <- resp
    invars$expl <- fixs
    if (length(rp) != 0) 
      rp <- rp[order(names(rp), decreasing = TRUE)]
    if (length(rp) != 0) 
      invars$rp <- rp
    if (length(nonfps) != 0) 
      invars$nonfp <- nonfps
    invars$D <- D
    invars$indata <- categobj$indata
    invars <- as.list(invars)
  }
  
  if (D[1] == "Normal") {
    nlev <- length(levID)
    nleft <- length(left)
    
    indata <- get.Idata(left, indata)
    indata <- get.Interdata(left, indata)
    xpoly <- get.polydata(left, indata)
    if (length(xpoly$newleft) > 0) {
      left <- xpoly$newleft
      indata <- xpoly$indata
    }
    rm(xpoly)
    categobj <- get_categstr(left, indata)
    categstr3 <- categobj$categstr
    categstr0 <- names(categstr3)
    
    fixs.no <- grep("0+\\|", left)
    fixs <- left[fixs.no]
    if (length(fixs) != 0) {
      fixs <- unlist(strsplit(fixs, "\\|"))
      fixs <- get.terms(fixs[2])
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
        rands[[i]] <- get.terms(rands[[i]][2])
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
    
    # substitute categorical variables
    
    if (length(categstr3) != 0) {
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
      if (length(nonfps) != 0) {
        delpos <- numeric(0)
        addterms <- numeric(0)
        for (ii in 1:length(nonfps)) {
          replacepos <- nonfps[ii] %in% categstr0
          if (replacepos) {
            taddterms <- categstr3[[nonfps[ii]]]
            delpos <- c(delpos, ii)
            addterms <- c(addterms, taddterms)
          }
        }
        if (length(delpos) > 0) {
          nonfps <- nonfps[-delpos]
          nonfps <- c(nonfps, addterms)
        }
      }
      if (length(rp) != 0) {
        for (ii in 1:length(rp)) {
          delpos <- numeric(0)
          addterms <- numeric(0)
          for (jj in 1:length(rp[[ii]])) {
            replacepos <- rp[[ii]][jj] %in% categstr0
            if (replacepos) {
              delpos <- c(delpos, jj)
              taddterms <- categstr3[[rp[[ii]][jj]]]
              addterms <- c(addterms, taddterms)
            }
          }
          if (length(delpos) > 0) {
            rp[[ii]] <- rp[[ii]][-delpos]
            rp[[ii]] <- c(rp[[ii]], addterms)
          }
        }
      }
    }
    
    invars <- new.env()
    invars$levID <- levID
    invars$resp <- resp
    invars$expl <- fixs
    if (length(rp) != 0) 
      rp <- rp[order(names(rp), decreasing = TRUE)]
    if (length(rp) != 0) 
      invars$rp <- rp
    if (length(nonfps) != 0) 
      invars$nonfp <- nonfps
    invars$D <- D
    invars$indata <- categobj$indata
    invars <- as.list(invars)
  }
  return(invars)
} 
