################################################################################
######################## summary method and its print method ###################
################################################################################

################# summary method for multiPIM objects ##########################

summary.multiPIM <- function(object,
                             type = c("both", "statistical", "time"),
                             use.plug.in.se = is.null(object$boot.param.array),
                             alternative.se.matrix = NULL,
                             two.sided.p.vals = TRUE,
                             bf.multiplier = object$num.exp * object$num.out,
                             by.exposure = TRUE,
                             digits = 4,
                             ...) {

  ## add some error checking to this function eventually

  type = match.arg(type)

  ## do statistical summary stuff
  
  if(type != "time") {

    sum.attributes <- c("param.estimate", "stand.error", "test.stat",
                        "p.val", "p.val.bon.adj")

    ## instantiate an array to hold summary info

    summary.array <- array(0, dim = c(dim(object$param.estimates),
                                      length(sum.attributes)))

    ## param estimates

    summary.array[,,1] <- object$param.estimates

    ## stand.errs

    if(!is.null(alternative.se.matrix)) {
      summary.array[,,2] <- alternative.se.matrix
      stand.err.type <- "alternative"
    } else if(use.plug.in.se) {

      if(identical(object$plug.in.stand.errs, NA)
         && object$estimator == "G-COMP")
        stop("When using the G-COMP estimator, no plug-in standard errors are",
             "\navailable and thus no summary table can be provided.",
             "\nStandard errors can be estimated by bootstrapping\n",
             "e.g. by using the multiPIMboot function instead of multiPIM.\n",
             "Or, to see a breakdown of running time by super learner candidate, use\n",
             'summary(multiPIM.result, type = "time")')

      summary.array[,,2] <- object$plug.in.stand.errs
      stand.err.type <- "plug.in"
    } else if(!use.plug.in.se) {
      summary.array[,,2] <- apply(object$boot.param.array, c(2,3), sd)
      stand.err.type <- "bootstrap"
    } else stop("internal error #1")

    ## test stats = param estimates / stand errs

    summary.array[,,3] <- abs(summary.array[,,1] / summary.array[,,2])

    ## calculate p values

    multiplier.1 <- ifelse(two.sided.p.vals, 2, 1)
    summary.array[,,4] <- multiplier.1 * pnorm(summary.array[,,3],
                                               lower.tail = FALSE)

    ## get bonferroni p vals

    summary.array[,,5] <- bf.multiplier * summary.array[,,4]
    summary.array[,,5][summary.array[,,5] > 1] <- 1

    ## add names

    dimnames(summary.array) <- c(dimnames(object$param.estimates),
                                 list(sum.attributes))

    stat.list <- list(summary.array = summary.array,
                       two.sided.p.vals = two.sided.p.vals,
                       stand.err.type = stand.err.type,
                       bf.multiplier = bf.multiplier,
                       by.exposure = by.exposure)

  } else stat.list <- vector("list", length = 0)

  ## do timing summary stuff
  
  if(type != "statistical") {
    
    ## get summary of time taken broken down by g vs. Q modeling

    g.Q.time.frame <- data.frame(c(object$g.method, object$Q.method),
                                 c(object$g.time, object$Q.time),
                                 c(object$g.time/object$main.time * 100,
                                   object$Q.time/object$main.time * 100))

    names(g.Q.time.frame) <- c("method", "seconds", "%.of.total")
    rownames(g.Q.time.frame) <- c("g modeling", "Q modeling")

    if(object$estimator != "G-COMP" && object$g.method == "sl") {
  
      g.sl.xval.time.mat <- cbind(object$g.sl.cand.times,
                                  object$g.sl.cand.times * 100
                                  / object$g.sl.time)

      dimnames(g.sl.xval.time.mat) <- list(object$g.sl.cands,
                                           c("seconds", "%.of.g.x-val.time"))

    } else g.sl.xval.time.mat <- NA

    if(object$estimator != "IPCW" && object$Q.method == "sl") {
    
      Q.sl.xval.time.mat <- cbind(object$Q.sl.cand.times,
                                  object$Q.sl.cand.times * 100
                                  / object$Q.sl.time)

      dimnames(Q.sl.xval.time.mat) <- list(object$Q.sl.cands,
                                           c("seconds", "%.of.Q.x-val.time"))

    } else Q.sl.xval.time.mat <- NA

    time.list <- list(main.time = object$main.time,
                      g.time = object$g.time,
                      Q.time = object$Q.time,
                      g.Q.time.frame = g.Q.time.frame,
                      g.sl.time = object$g.sl.time,
                      Q.sl.time = object$Q.sl.time,
                      g.sl.xval.time.mat = g.sl.xval.time.mat,
                      Q.sl.xval.time.mat = Q.sl.xval.time.mat)

  } else time.list <- vector("list", length = 0)

  summary.list <- c(list(type = type, digits = digits, call = object$call),
                    stat.list,
                    time.list)

  class(summary.list) <- "summary.multiPIM"
  return(summary.list)
}

###################### print method for summary objects ########################

print.summary.multiPIM <- function(x, by.exposure, digits, ...) {

  if(missing(digits)) digits <- x$digits

  ## print the call

  cat("\nThe call was:\n\n")
  print(x$call)
  
  ## print the statistical summary part

  if(x$type != "time") {

    if(missing(by.exposure)) by.exposure <- x$by.exposure

    cat("\n")

    if(all(dim(x$summary.array)[c(1,2)] == 1)) {
      ## i.e. if num.exposures==1 and num.outcomes==1 

      cat("Results for the exposure \"", dimnames(x$summary.array)[[1]],
          "\" vs the outcome \"", dimnames(x$summary.array)[[2]], "\"\n\n",
          sep = "")
      print(x$summary.array[1, 1, ], digits = digits, ...)
      cat("\n")

    } else if(dim(x$summary.array)[1] == 1) {
      ## i.e. if num.exposures==1 

      cat("Results for the exposure \"", dimnames(x$summary.array)[[1]],
          "\" vs the outcomes listed on the left:\n", sep = "")
      print(x$summary.array[1,,], digits = digits, ...)
      cat("\n")

    } else if(dim(x$summary.array)[2] == 1) {
      ## i.e. if num.outcomes==1
    
      cat("Results for the exposures listed on the left vs the outcome \"",
          dimnames(x$summary.array)[[2]], "\":\n", sep = "")
      print(x$summary.array[,1,], digits = digits, ...)
      cat("\n")

    } else if(by.exposure) {

      for(i in 1:dim(x$summary.array)[1]) {
        cat("Results for the exposure \"", dimnames(x$summary.array)[[1]][i],
            "\" vs the outcomes listed on the left:\n", sep = "")
        print(x$summary.array[i,,], digits = digits, ...)
        cat("\n")
      }

    } else {

      for(i in 1:dim(x$summary.array)[2]) {
        cat("Results for the exposures listed on the left vs the outcome \"",
            dimnames(x$summary.array)[[2]][i], "\":\n", sep = "")
        print(x$summary.array[,i,], digits = digits, ...)
        cat("\n")
      }
    }
  } # end statistical part
  
  ## Print the time summaries

  if(x$type != "statistical") {

    cat("\nTotal time for main multiPIM run:\n\n")
    cat(x$main.time, "seconds\n\n")
    cat("Breakdown by g vs. Q modeling:\n\n")
    print(x$g.Q.time.frame, digits = digits)

    if(is.matrix(x$g.sl.xval.time.mat)) {

      cat("\nTotal time for g model super learner",
          "cross validation (x-val):\n\n")
      cat(x$g.sl.time, "seconds\n\n")
      cat("Breakdown by candidate:\n\n")
      print(x$g.sl.xval.time.mat, digits = digits)
    }

    if(is.matrix(x$Q.sl.xval.time.mat)) {

      cat("\nTotal time for Q model super learner",
          "cross validation (x-val):\n\n")
      cat(x$Q.sl.time, "seconds\n\n")
      cat("Breakdown by candidate:\n\n")
      print(x$Q.sl.xval.time.mat, digits = digits)
    }
  } # end timing part

  invisible(x)
}
