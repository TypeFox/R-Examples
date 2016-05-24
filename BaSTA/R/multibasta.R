# Function to run and compare multiple BaSTA models:
multibasta <- function(object, studyStart, studyEnd, models, 
                               shapes, ...) {
  argList <- list(...)
  if (is.null(argList$nsim)) nsim <- 2 else nsim <- argList$nsim
  mods <- rbind(data.frame(model = "EX", shape = "simple"),
                expand.grid(model = c("GO", "WE", "LO"), 
                            shape = c("simple", "Makeham", "bathtub")))
  
  modsAbr <- rbind(data.frame(model = "Ex", shape = "Si"),
                   expand.grid(model = c("Go", "We", "Lo"), 
                               shape = c("Si", "Ma", "Bt")))
  if (missing(models)) models <- c("EX", "GO", "WE", "LO")
  if (missing(shapes)) shapes <- c("simple", "Makeham", "bathtub")
  idIncl <- which(mods$model %in% models & mods$shape %in% shapes)
  mods <- mods[idIncl, ]
  modsAbr <- modsAbr[idIncl, ]
  nmod <- nrow(mods)
  modNames <- paste(modsAbr$model, modsAbr$shape, sep = ".")
  if (nsim > 1) DIC <- TRUE else DIC <- FALSE
  if (DIC) {
    dics <- matrix(NA, nmod, 5, dimnames = 
                     list(modNames, c("D.ave", "D.mode", "pD", "k", "DIC")))
  }
  runList <- list()
  for (mod in 1:nmod) {
    updt <- sprintf("Run number %s, model: %s", mod, modNames[mod])
    cat(paste("\n", paste(rep("-", nchar(updt)), collapse = ""), sep = ""))
    cat(sprintf("\n%s\n", updt))
    cat(paste(paste(rep("-", nchar(updt)), collapse = ""), "\n", sep = ""))
    if (is.null(argList$nsim)) {
      out <- basta(object, studyStart, studyEnd, model = mods$model[mod],
                   shape = mods$shape[mod], nsim = nsim, ...)
    } else {
      out <- basta(object, studyStart, studyEnd, model = mods$model[mod],
                   shape = mods$shape[mod],  ...)
    }
    runList[[modNames[mod]]] <- out
    if (DIC & out$DIC[1] != "Not calculated") {
      dics[mod, ] <- out$DIC
    }
  }
  if (DIC) {
    idNa <- which(is.na(dics[, "DIC"]))
    idnNa <- which(!is.na(dics[, "DIC"]))
    idModFit <- c(idnNa[sort.int(dics[idnNa, 5], index.return = TRUE)$ix], idNa)
    if (is.na(idModFit[1])) idModFit <- 0
    modFit <- data.frame(mods, dics)[idModFit, ]
    DICdiff <- modFit$DIC - modFit$DIC[1]
    Rank <- 1:nmod
    modFit <- data.frame(modFit, DICdiff, Rank)
  } else {
    modFit <- data.frame(DICs = "Not Calculated")
    idModFit <- NA
  }
  results <- list(runs = runList[idModFit], DICs = modFit, models = mods)
  class(results) <- "multibasta"
  return(results)
}


print.multibasta <- function(x, ...) {
  argList <- list(...)
  if (is.null(argList$digits)) digits <- 3 else digits <- argList$digits
  cat("\nDICs:\n")
  print(x$DICs, digits = digits)
}

summary.multibasta <- function(object, ...) {
  argList <- list(...)
  if (is.null(argList$digits)) digits <- 3 else digits <- argList$digits
  if ("version" %in% names(object$runs[[1]])) {
    cat(sprintf("\nBaSTA version %s\n", object$runs[[1]]$version))
  }
  
  cat("\nCall:\n")
  cat(paste("Covars. structure \t\t: ", object$runs[[1]]$modelSpecs[3], "\n", 
            sep = ""))
  cat(paste("Minimum age       \t\t: ", object$runs[[1]]$modelSpecs[4], "\n", 
            sep = ""))
  cat(paste("Cat. covars.      \t\t: ", object$runs[[1]]$modelSpecs[5], "\n", 
            sep = ""))
  cat(paste("Cont. covars.     \t\t: ", object$runs[[1]]$modelSpecs[6], "\n", 
            collapse = ""))
  
  cat("\nModel settings:\n")
  print(object$runs[[1]]$set)

  cat("\nDICs:\n")
  print(object$DICs, digits = digits)
}

coef.multibasta <- function(object, showAll = FALSE, ...) {
  argList <- list(...)
  if (is.null(argList$digits)) digits <- 3 else digits <- argList$digits
  ans <- list()
  if (class(showAll) == "logical") {
    if (showAll) nmod <- 1:length(object$runs) else nmod <- 1
  } else if (class(showAll) == 'numeric') {
    if (length(showAll) == 1) {
      nmod <- 1:showAll
    } else {
      nmod <- showAll
    }
    nmod <- nmod[nmod <= length(object$runs)]
  } else {
    stop("\nArgument 'showAll' should be logical or an integer.\n", 
         call. = FALSE)
  }
  modNames <- names(object$runs)
  for (mod in nmod) {
    updt <- sprintf("DIC rank %s, model: %s", as.character(object$DICs$Rank)[mod], 
                    modNames[mod])
    cat(paste("\n", paste(rep("-", nchar(updt)), collapse = ""), sep = ""))
    cat(sprintf("\n%s\n", updt))
    cat(paste(paste(rep("-", nchar(updt)), collapse = ""), "\n", sep = ""))
    print.default(object$runs[[mod]]$coefficients, digits = digits)
    ans$coefficients[[modNames[mod]]] <- object$runs[[mod]]$coefficients
  }
  ans$DICs <- object$DICs
  return(invisible(ans))
}


