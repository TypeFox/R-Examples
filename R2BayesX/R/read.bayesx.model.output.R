read.bayesx.model.output <- function(dir, model.name)
{
  if(is.null(dir))
    stop("no directory specified!")
  if(is.null(model.name))
    stop("no model name specified!")
  files <- dir.files <- list.files(dir)
  if(!any(grep(model.name,files)))
    stop(paste("no model results existing for ", model.name, "!", sep = ""))
  else {
    fileext <- function(x) {
      pos <- regexpr("\\.([[:alnum:]]+)$", x)
      ifelse(pos > -1L, substring(x, pos + 1L), "")
    }
    rval <- list()
    files <- grep(model.name, files, value = TRUE, fixed = TRUE)
    files <- files[fileext(files) != "ps"]
    filep <- grep(paste(model.name, ".", sep = ""), files, value = TRUE, fixed = TRUE)  
    files <- c(grep(paste(model.name, "_", sep = ""), files, value = TRUE, fixed = TRUE), filep)  
    info <- paste(model.name, ".terms.info", sep = "")
    mformula <- paste(model.name, ".formula.rda", sep = "")
    minfo <- NULL
    if(info %in% dir.files) {
      minfo <- readLines(file.path(dir, info))
      minfo <- eval(parse(text = minfo[length(minfo)]))
    }
    ## search for data
    data <- isnadata <- NULL
    if(length(i <- grep(paste(model.name, ".data.raw", sep = ""), files))) {
      data <- as.matrix(read.table(file.path(dir, files[i]), header = TRUE))
      if(any(is.na(data))) {
        isnadata <- rowSums(is.na(data) * 1) > 0
        data <- data[!isnadata, ]
      }
    }
    for(char in c("_predict.raw", "_predictmean.raw", "_predict.res"))
      if(length(i <- grep(char, files)))
        data <- cbind(data, as.matrix(read.table(file.path(dir, files[i]), header = TRUE)))

    ## set response and predictor
    N <- NA
    response <- eta <- residuals <- NULL
    if(!is.null(data)) {
      data[data == "."] <- "NA"
      warn <- getOption("warn")
      options("warn" = -1)
      mode(data) <- "numeric"
      options("warn" = warn)
      data <- as.data.frame(data)
      data$intnr <- NULL
      dn <- unique(names(data))
      data <- data[dn]
      response <- data[[1L]]
      N <- length(response)
      rval$fitted.values <- eta <- get.eta(data)
      if(!is.null(minfo)) {
        if(!is.null(minfo$order)) {
          ooo <- 1:NROW(data) ## order(eval(parse(text = minfo$order))) FIXME!
          data <- data[ooo, ]
          rownames(data) <- 1:NROW(data)
          if(!is.null(response))
            response <- response[ooo]
          if(!is.null(eta)) {
            if(!is.null(isnadata) && length(eta) == length(isnadata)) {
              if(is.matrix(eta))
                eta <- eta[!isnadata,]
              else
                eta <- eta[!isnadata]
            }
            if(is.matrix(eta) || is.data.frame(eta)) {
              if(nrow(eta) == length(ooo))
                eta <- eta[ooo,]
            } else {
              if(length(eta) == length(ooo))
                eta <- eta[ooo]
            }
          }
          if(is.matrix(eta)) {
            if(nrow(eta) != nrow(data))
              eta <- NULL
            else
              rownames(eta) <- 1:NROW(eta)
          } else {
            if(nrow(data) != length(eta))
              eta <- NULL
          }
        }
        if(!is.null(minfo$YLevels)) {
          response <- as.factor(response)
          YLevels <- eval(parse(text = minfo$YLevels))
          levels(response) <- YLevels
          if(is.matrix(eta)) {
            nYLevels <- eval(parse(text = minfo$nYLevels))
            cne <- cne <- colnames(eta)
            for(k in 1:length(cne)) {
              tmp1 <- strsplit(cne[k], "eta")[[1L]]
              tmp2 <- strsplit(cne[k], "mu")[[1L]]
              if(length(tmp1) > 1) {
                if(tmp1[2L] %in% nYLevels) {
                  cne[k] <- paste("eta:", YLevels[nYLevels == tmp1[2L]], sep = "")
                }
              }
              if(length(tmp2) > 1) {
                if(tmp2[2L] %in% nYLevels)
                  cne[k] <- paste("mu:", YLevels[nYLevels == tmp2[2L]], sep = "")
              }
            }
            colnames(eta) <- cne
            rval$fitted.values <- eta
          }
        }
      } 
      if(!is.factor(response))
        rval$residuals <- response - eta
      rval$response <- response
    }

    ## get smooth and random effects
    rval <- c(rval, find.smooth.random(dir, files, data, response, eta, 
      model.name, minfo, file.path(dir, info)))

    ## get fixed effects
    rval <- find.fixed.effects(dir, files, data, response, eta, model.name, 
      rval, minfo, file.path(dir, info))

    ## get scale estimate
    rval$variance <- get.scale(files, dir)

    ## search for other results
    model.results <- mf <- mf2 <- NULL
    method <- ""
    if(any(grep("deviance.raw", files))) {
      mf <- grep("deviance.raw", files, value = TRUE)
      mf <- read.table(file.path(dir, mf), header = TRUE)
      pd <- mf$unstandardized_deviance[length(mf$unstandardized_deviance) - 1L]
      DIC <- mf$unstandardized_deviance[length(mf$unstandardized_deviance)]
      mf <- list(DIC = DIC, pd = pd)
    }
    if(any(grep(".tex", files))) {
      sm <- readLines(file.path(dir, grep(".tex", files, value = TRUE)))
      model.results <- search.bayesx.tex(sm)
      method <- model.results$method
    }
    if(any(grep("modelfit.raw", files))) {
      mf <- grep("modelfit.raw", files, value = TRUE)
      mf <- chacol(read.table(file.path(dir, mf), header = TRUE))
      mf <- mf2 <- as.list(mf)
    }
    if(!is.null(model.results)) {
      if(!is.null(mf2)) {
        n1 <- names(mf2)
        n2 <- names(model.results)
        for(i in 1L:length(mf2))
          if(!n1[i] %in% n2)
            eval(parse(text = paste("model.results$", n1[i], "<- mf2[[i]]", sep = "")))
      }
      mf <- model.results
      if(!is.null(mf$smooth.hyp.step) & is.null(rval$smooth.hyp)) {
        rval$smooth.hyp <- mf$smooth.hyp.step
        if(length(log <- grep(".log", files, fixed = TRUE, value = TRUE))) {
          log <- readLines(file.path(dir, log))
          if(length(i <- grep("Final Model:", log))) {
            i <- i[length(i)]
            log <- log[i:length(log)][3]
            if(length(splitme(log))) {
              log <- strsplit(log, " + ", fixed = TRUE)[[1L]]
              if(length(log <- grep("(random", log, fixed = TRUE, value = TRUE))) {
                for(k in 1L:length(log)) {
                  tmp <- strsplit(log[k], ",", fixed = TRUE)[[1L]]
                  if(length(tmp <- grep("df=", tmp, fixed = TRUE, value = TRUE))) {
                    tmp <- strsplit(tmp, "df=", fixed = TRUE)[[1L]][2L]
                    if(length(splitme(tmp))) {
                      tmp <- as.numeric(tmp)
                      cn <- colnames(rval$smooth.hyp)
                      if(any(cn == "df")) {
                        tmpv <- rval$smooth.hyp[cn == "df"]
                        rn <- rownames(rval$smooth.hyp)
                        rn[tmpv == tmp] <- paste(gsub("s(", "sx(", rn[tmpv == tmp], fixed = TRUE),
                          "re", sep = ":")
                        rownames(rval$smooth.hyp) <- rn
                      }
                    }
                  }
                }
              }
            }
          }
        }
        if(!is.null(rval$smooth.hyp))
          rownames(rval$smooth.hyp) <- gsub("s(", "f(", rrmfs(rownames(rval$smooth.hyp)), fixed = TRUE)
        mf$smooth.hyp.step <- NULL
      }
    } else mf <- c(list(method = method, N = N), mf)
    if(!is.null(mf$smooth.hyp.step))
      mf$smooth.hyp.step <- NULL
    if(mf$method == " ")
      mf$method <- "NA"
    if(!is.null(mf$logLik))
      mf$logLik <- mf$logLik/(-2)
    rval$model.fit <- mf

    ## reordering and naming
    if(info %in% dir.files)
      rval$effects <- term.reorder(rval$effects, file.path(dir, info))
    rval$effects <- delete.NULLs(rval$effects)
    if(any(duplicated(names.eff <- names(rval$effects)))) {
      for(k in names.eff) {
        if(sum(which <- names.eff == k) > 1L) {
          for(j in 1:length(which))
            if(which[j]) {
              cw <- gsub(".bayesx", "", class(rval$effects[[j]])[1L], fixed = TRUE)
              cw <- gsub("random", "re", cw)
              cw <- paste(names.eff[j], cw, sep = ":")
              names.eff[j] <- attr(rval$effects[[j]], "specs")$label <- cw
            } 
        }
      }
      names(rval$effects) <- names.eff
    }

    ## was there a .prg/.log file?
    if(length(prg <- grep(".prg", files, value = TRUE))) {
      for(j in prg) {
        if(grepl(".log", j))
          rval$bayesx.run <- list("log" = readLines(file.path(dir, j)))
        else
          rval$bayesx.prg <- list("prg" = readLines(file.path(dir, j)))
      }
    }

    ## search for additional info
    rval$model.fit <- smi(file.path(dir, info), rval$model.fit)

    ## new with HMCMC
    if(any(grep("_DIC.res", files))) {
      dic <- grep("_DIC.res", files, value = TRUE)
      dic <- chacol(read.table(file.path(dir, dic), header = TRUE))
      rval$model.fit <- c(rval$model.fit, as.list(dic))
    }

    ## reformate output
    rval <- bayesx.reformate(rval)

    ## get long formulas
    if(mformula %in% files) {
      nenv <- new.env()
      load(file.path(dir, mformula), envir = nenv)
      rval$model.fit$formula <- get("f", envir = nenv)
      unlink(file.path(dir, mformula))
    }

    ## get log file
    if(length(log <- grep(".log", files, fixed = TRUE, value = TRUE)))
      rval$logfile <- readLines(file.path(dir, log[1]))

    return(rval)
  }
}
