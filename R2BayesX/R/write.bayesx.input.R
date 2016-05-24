write.bayesx.input <- function(object)
{
  if(is.null(object) || missing(object))
    stop("nothing to do!")
  if(class(object) != "bayesx.input")
    stop("object must be of class bayesx.input")
  data.file <- NULL
  if(is.null(object$data))
    object$data <- environment(object$formula)
  else {
    if(is.character(object$data))
      data.file <- object$data
  }
  if(is.null(object$outfile)) {
    if(.Platform$OS.type == "windows") {
      tmp <- splitme(tempdir())
      for(i in 1L:length(tmp))
        if(tmp[i] == "\\")
          tmp[i] <- "/"
      object$outfile <- paste(resplit(tmp), "/bayesx", sep = "")
    } else object$outfile <- paste(tempdir(), "/bayesx", sep = "")
  }
  if(!is.character(object$outfile))
    stop("argument outfile must be a character!")
  if(is.null(object$hlevel))
    is.h <- FALSE
  else {
    if(object$hlevel > 1L)
      is.h <- TRUE
    else
      is.h <- FALSE
  }
  if(!file.exists(object$outfile) || object$replace) {
    if(!is.null(object$outfile) && length(fdir <- list.files(object$outfile))) {
      file.remove(paste(object$outfile, "/", fdir, sep = "")) 
    }
    dir.create(object$outfile, showWarnings = FALSE, recursive = TRUE)
  } else {
    if(!is.h && !object$replace) {
      wok <- TRUE
      count <- 0L
      while(wok) {
        count <- count + 1L
        if(!file.exists(paste(object$outfile, count, sep = ""))) {
          object$outfile <- paste(object$outfile, count, sep = "")
          wok <- FALSE
        }
      }
      dir.create(object$outfile, showWarnings = FALSE)
      cat("Note: created new output directory \'", object$outfile, "\'!\n", sep = "")
    } else {
      if(is.h)
        object$outfile <- object$houtfile
    }
  }
  wd <- as.character(getwd())
  if(is.null(object$hlevel) || object$hlevel < 2L)
    setwd(object$outfile)
  prg.file <- paste(object$model.name, ".input.prg", sep = "")
  thismodel <- NULL
  if(file.exists(prg.file) && !is.h) {
    wok <- TRUE
    thismodel <- 0L
    while(wok) {
      thismodel <- thismodel + 1L
      prg.file <- paste(object$model.name, ".input", thismodel, ".prg", sep = "")
      if(!file.exists(prg.file))
        wok <- FALSE
    }
  }
  if(!is.h) {
    cat(paste("% usefile ", object$outfile, "/", prg.file, "\n\n", sep = ""), 
      file = prg.file, append = FALSE)
    cat(paste("logopen using ", object$outfile, "/", prg.file, ".log", "\n\n", sep = ""),
      file = prg.file, append = TRUE)
    if(object$method == "MCMC") {
      if(is.null(object$mcmcreg) && object$first && !object$hmcmc) {
        cat("bayesreg b\n\n", file = prg.file, append = TRUE)
      }
      if((!is.null(object$mcmcreg) && object$first) || object$hmcmc) {
        cat("mcmcreg b\n\n", file = prg.file, append = TRUE)
      }
    }
    if(object$method == "REML")
      cat("remlreg b\n\n", file = prg.file, append = TRUE)
    if(object$method == "STEP")
      cat("stepwisereg b\n\n", file = prg.file, append = TRUE)   
  }

  ## hierarchical models set here 
  if(!is.null(object$h.random))
    for(k in 1:length(object$h.random)) {
      object$h.random[[k]]$houtfile <- object$outfile
      write.bayesx.input(object$h.random[[k]])
    }
  dropvars <- NULL
  if(is.null(data.file)) {
    tf <- as.character(object$formula)
    cnd <- colnames(object$data)
#    for(k in 1L:ncol(object$data)) {
#      if(is.factor(object$data[,k])){
#        if(is.null(object$contrasts[[cnd[k]]]))
#          contrasts(object$data[,k]) = contr.sum(nlevels(object$data[,k]))
#      }
#    }
    if(grepl(tf[2L], tf[3L], fixed = TRUE)) {
      dat <- model.matrix(as.formula(paste(tf[1L], tf[3L])), 
        data = object$data, contrasts.arg = object$contrasts)
    } else {
      dat <- model.matrix(object$formula, data = object$data,
        contrasts.arg = object$contrasts)
    }
    nc <- ncol(dat)
    if(!is.null(model.offset(object$data))) {
      nc <- ncol(dat)
      dat <- cbind(dat, model.offset(object$data))
      colnames(dat)[nc + 1L] <- "ModelOffset"
    }
    if(!is.null(model.weights(object$data))) {
      nc <- ncol(dat)
      dat <- cbind(dat, model.weights(object$data))
      colnames(dat)[nc + 1L] <- "ModelWeights"
    }
    nd <- rmf(names(object$data))
    names(object$data) <- nd
    colnames(dat) <- rmf(colnames(dat))
    intcptcheck <- intcpt <- attr(object$terms, "intercept") == 1L
    for(sf in nd) {
      if(is.factor(object$data[[sf]])) {
        if(intcpt) {
          ff <- as.data.frame(eval(parse(text = paste("model.matrix(~ -1 +", sf,")", 
            sep = "")), envir = object$data))
        } else {
          ff <- as.data.frame(eval(parse(text = paste("model.matrix(~ 1 +", sf,")", 
            sep = "")), envir = object$data))
        }
        lf <- colnames(ff)
        vars <- colnames(dat)
        if(!all(lf %in% vars)) {
          mf <- lf[!lf %in% vars]
          dat <- cbind(dat, as.matrix(ff[mf]))
        }
      }
    }
    object$Yn <- rmf(object$Yn)
    dat <- cbind(as.vector(object$data[[object$Yn]]), dat)
    colnames(dat)[1L] <- object$Yn
    vars <- colnames(dat)
    nc <- ncol(dat)
    intcpt <- FALSE
    if("(Intercept)" %in% vars && nc > 1L) {
      dat <- matrix(dat[, !vars %in% "(Intercept)"], ncol = (nc - 1L))
      colnames(dat) <- vars[!vars %in% "(Intercept)"]
      intcpt <- TRUE
    }
    vars <- rmf(colnames(dat))
    for(i in 1L:length(vars)) {
      vars[i] <- gsub("(offset)", "ModelOffset", vars[i], fixed = TRUE)
      vars[i] <- gsub("(weights)", "ModelWeights", vars[i], fixed = TRUE)
    }
    colnames(dat) <- vars
    vars <- vars[vars != object$Yn]
    if(intcpt)
      vars <- c("(Intercept)", vars)
    if(!is.null(object$hlevel)) {
      if(object$hlevel > 1L) {
        object$model.name <- paste(object$model.name, "_hlevel", 
          object$hlevel, "_RANDOM_EFFECTS", sep = "")
      } else {
        object$model.name <- paste(object$model.name, "_hlevel", 
          object$hlevel, "_MAIN_REGRESSION", sep = "")
      }
    }
    data.file <- paste(object$outfile, "/", object$model.name, ".data.raw", sep = "")
    if(!is.null(object$begin.vec)) {
      dn <- c(colnames(dat), object$begin)
      dat <- cbind(dat, object$begin.vec)
      colnames(dat) <- dn
    }
    if(!file.exists(data.file))
      write.table(dat, data.file, col.names = TRUE, row.names = FALSE, quote = FALSE)
    else {
      wok <- TRUE
      i <- 1L
      while(wok) {
        data.file <- paste(object$outfile, "/", object$model.name, ".data", i, ".raw", sep = "")
        i <- i + 1L
        if(!file.exists(data.file)) {
          write.table(dat, data.file, col.names = TRUE, row.names = FALSE, quote = FALSE)
          wok <- FALSE
        }
      }
    }
    terms <- attr(object$terms, "term.labels")
    infofile <- rdafile <- paste(object$outfile, "/", object$model.name, thismodel, sep = "")
    infofile <- paste(infofile, ".terms.info", sep = "")
    rdafile <- paste(rdafile, ".formula.rda", sep = "")
    write.term.info(infofile, terms, object$data, object, 
      contrasts.arg = object$contrasts, rdafile = rdafile)
  } else {
    vdf <- sub(" +$", "", readLines(data.file, n = 1L))
    vars <- strsplit(vdf, " ")[[1L]]
    ok <- FALSE
    if(length(vars) > 1L)
      if(all(attr(terms(object$formula), "term.labels") %in% vars))
        ok <- TRUE
    if(!ok)
      warning("variable names in specified data file do not match with formula variable names!")
  }
  st <- split.terms(attr(object$terms, "term.labels"), vars, object$data, dropvars, intcptcheck)
  bayesx.prg <- write.prg.setup(rmf(object$Yn), object, prg.file, data.file,
    thismodel, st)
  if(!is.h)
    cat("logclose \n", file = prg.file, append = TRUE)
  if(is.null(object$hlevel) || object$hlevel < 2L)
    setwd(wd)

  return(invisible(list(prg = bayesx.prg, prg.name = prg.file,
    model.name = object$model.name, file.dir = object$outfile)))
}

