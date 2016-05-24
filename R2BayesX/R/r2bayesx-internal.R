.print_bayesx <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  if(!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  } else {
    if(!is.null(x$model.fit$formula)) {
      cat("Formula:\n")
      if(is.character(x$model.fit$formula))
        cat(x$model.fit$formula, "\n")
      else
        print(x$model.fit$formula)
    }
  }
  if(!is.null(x$model.fit)) {
    cat("Summary:\n")
    mfn <- names(x$model.fit)
    mfn <- mfn[mfn != "formula" & mfn != "order" & 
      mfn != "YLevels" & mfn != "nYLevels" & 
      mfn != "model.name"]
    step <- 5L
    for(i in 1L:length(mfn)) {
      txt <- x$model.fit[[mfn[i]]]
      if(is.numeric(txt))
        txt <- round(txt, digits)
      txt <- deparse(txt)
      if(i < step) {
        if(!is.null(txt) && txt != "") {
          if(mfn[i] != "step.final.model")
            cat(mfn[i], "=", gsub('\"', "", txt, fixed = TRUE), " ")
#          else {
#            cat("\n\n")
#            cat("Stepwise final model:\n")
#            cat(gsub('\"', "", txt, fixed = TRUE))
#            cat("\n\n")
#          }
        }
      }
      if(i == step) {
        if(i != length(mfn))
          cat("\n")
        step <- step + step
      }
    }
  cat("\n")
  }

  return(invisible(NULL))
}

.print_summary_bayesx <- function(x, digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...)
{
  if(!is.null(x$model.fit))
    if(!is.null(x$model.fit$model.name))
      if(length(grep("_hlevel", x$model.fit$model.name))) {
        hlevel <- splitme(strsplit(x$model.fit$model.name, "_hlevel")[[1]][2])
        go <- TRUE
        hl <- NULL
        for(i in 1:length(hlevel)) {
          if(hlevel[i] == "_")
            go <- FALSE
          if(go)
            hl <- c(hl, hlevel[i])
        }
        hlevel <- as.integer(resplit(hl))
        if(hlevel > 1)
          cat("Hierarchical random effects model results: stage", hlevel, "\n")
        else {
          cat("Main effects model results: stage", hlevel, "\n")
          cat("\n")
        }
      }
  if(!is.null(x$call)) {
    cat("Call:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  } else {
    if(!is.null(x$model.fit$formula)) {
      cat("Formula:\n")
      if(is.character(x$model.fit$formula))
        cat(x$model.fit$formula, "\n")
      else
        print(x$model.fit$formula, showEnv = FALSE)
    }
  }
  liner <- ""
  fc <- FALSE
  if(!is.null(x$fixed.effects)) {
    fc <- TRUE
    if(nrow(x$fixed.effects) < 2L) {
      if(!any(is.na(x$fixed.effects)) && all(x$fixed.effects[1L,] == 0))
        fc <- FALSE
    } else {
      if(!all(as.character(x$fixed.effects) == "NaN") && all(x$fixed.effects[1L,] == 0)) {
        m <- ncol(x$fixed.effects)
        nc <- colnames(x$fixed.effects)
        nr <- rownames(x$fixed.effects)[2L:nrow(x$fixed.effects)]
        x$fixed.effects <- matrix(x$fixed.effects[2L:nrow(x$fixed.effects),], ncol = m)
        colnames(x$fixed.effects) <- nc
        rownames(x$fixed.effects) <- nr
      }
    }
    x$fixed.effects <- round(x$fixed.effects, digits)
  }
  if(fc || (!is.null(x$smooth.hyp))) {
    cat(liner, "\n")
    cat("Fixed effects estimation results:\n")
    cat("\n")
  }
  if(fc) {
    cat("Parametric coefficients:\n")
    printCoefmat(x$fixed.effects)
  }
  if(!is.null(x$smooth.hyp)) {
    if(fc)
      cat("\n")
    if(x$model.fit$method == "MCMC" || x$model.fit$method == "HMCMC") 
      cat("Smooth terms variances:\n")
    else
      cat("Smooth terms:\n")
    ls <- ncol(x$smooth.hyp)
    terms <- colnames(x$smooth.hyp)	
    rn <- rownames(x$smooth.hyp)
    x$smooth.hyp <- round(x$smooth.hyp, digits)
    printCoefmat(x$smooth.hyp)
  }
  cat(liner, "\n")
  if(!is.null(x$random.hyp)) {
    cat("Random effects variances:\n")
    x$random.hyp <- round(x$random.hyp, digits)
    printCoefmat(x$random.hyp)		
    cat(liner, "\n")
  }		
  if(!is.null(x$model.fit)) {
    if(x$model.fit$method == "MCMC") {
      if(!is.null(x$variance)) {
        cat("Scale estimate:\n")
        x$variance <- round(x$variance, digits)
        printCoefmat(x$variance)
        cat(liner, "\n")
      }
    } else {
      if(!is.null(x$variance)) {
        cat("Scale estimate:", round(as.numeric(x$variance)[1], digits), "\n")
        cat(liner, "\n")
      }
    }
    x$model.fit <- delete.NULLs(x$model.fit)
    mfn <- names(x$model.fit)
    step <- 5L
    mfn <- mfn[!is.null(x$model.fit)]
    mfn <- mfn[mfn != "model.name"]
    mfn <- mfn[mfn != "formula"]
    mfn <- mfn[mfn != "step.final.model"]
    mfn <- mfn[mfn != "YLevels"]
    mfn <- mfn[mfn != "nYLevels"]
    mfn <- mfn[mfn != "order"]
    if(is.na(x$model.fit$N))
      x$model.fit$N <- "NA"
    if(x$model.fit$method == "")
      x$model.fit$method <- "NA"
    for(i in 1L:length(mfn)) {
        if(!is.null(x$model.fit[[mfn[i]]]) && !is.na(x$model.fit[[mfn[i]]] != "")) {
          if(length(splitme(as.character(x$model.fit[[mfn[i]]])))) {
            if(i < step)
              cat(mfn[i], "=", x$model.fit[[mfn[i]]], " ")
            if(i == step) {
              if(i != length(mfn))
                cat("\n")
              cat(mfn[i], "=", x$model.fit[[mfn[i]]], " ")
              step <- step + step
            }
          }
      }
    }
    cat("\n")
  }
#  if(!is.null(x$model.fit$step.final.model)) {
#    cat(liner,"\n")
#    cat("Stepwise final predictor choosen:\n")
#    cat("\n")
#    cat(x$model.fit$step.final.model, "\n")
#  }

  return(invisible(NULL))
}

