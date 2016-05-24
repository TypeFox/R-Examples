geekin <- function(formula,
                   family=gaussian, 
                   data,
                   weights,
                   subset,
                   id,
                   na.action,
                   control=geese.control(...),
                   varlist,
#                   type=c("pearson", "OR"),
                   ...
                   ) {

  # TODO
  # Pass everything correctly to geeglm

  # Check that the input is correct
  if (is.matrix(varlist))
    varlist <- list(varlist)
  if (!is.list(varlist))
    stop("varlist must be a list")

#  type <- match.arg(type)

  # Setup the initial glm call for a model without the correlation structure taken into account
  # This is used to check for missing data and is not necessary otherwise
  Call <- match.call()
  glmcall <- Call
  glmcall$id <- NULL
  glmcall$id <- glmcall$control <- glmcall$corstr <- glmcall$zcor <- glmcall$type <- NULL
  glmcall$varlist <- NULL
  glmcall[[1]] <- as.name("glm")
  glmFit <- eval(glmcall, parent.frame())
  mf <- model.frame(glmFit)

  
  # Extract clusters from formula
  # Can only accept ONE simple random effect
#  modterm <- attr(terms(formula), "term.labels")
#  clusters <- modterm[grep("\\|", attr(terms(formula), "term.labels"))]
#  if (length(clusters) != 1) {
#    stop("must provide exactly one random effect that determines the clusters")
#  }
#  if (length(grep("^1 \\| ", clusters))!=1) {
#    stop("only accepts a simple random effect to identify the clusters")
#  }

#  id <- clusters
  id <- id

  # Eval det rigtige sted

  # Check subset på variansmatricerne
 
#  print(get(id))
#  print(id)
  
  # Check that data are in the correct order
  if (!ordered.clusters(id)) {
    stop("the clusters must appear as contiguous blocks")
  }
  
  # The updated cluster index
  id.fixed <- id

  # The vector of observations to be removed
  missing.index <- as.numeric(glmFit$na.action)
  if (length(missing.index>0))
    id.fixed <- id.fixed[-missing.index]

  # for each element k in the varlist
  # extract the lower triangular matrix for use with the zcorr specification in geeglm
  varelements <- lapply(varlist, function(x) {
    if (length(missing.index>0)) {
        x <- x[-missing.index,-missing.index]
    }
    lower.tri.vector(x, cluster=id.fixed, diag=FALSE)
  })

  #  print(id.fixed)

  # Setup the proper call to geeglm
  geecall <- Call
  geecall$varlist <- NULL
  geecall$zcor <- do.call("cbind", varelements)
  geecall$corstr <- "userdefined"
  geecall[[1]] <- as.name("geeglm")
  geecall$control <- control
  geecall$id <- id.fixed
  geecall$data <- mf
  
  geeFit <- eval(geecall, parent.frame())

  # Cheat the output so the varlist and data don't look too horrible
  geeFit$call <- Call
  geeFit$na.action <- glmFit$na.action
  class(geeFit) <- c("geekin", class(geeFit))
  geeFit
}

         
  
print.geekin <- function (x, digits = NULL, quote = FALSE, prefix = "", ...) 
{
    xg <- x$geese
    if (is.null(digits)) {
        digits <- options()$digits
    } else { options(digits = digits) }
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(unclass(x$coefficients), digits = digits)
    cat("\nDegrees of Freedom:", length(x$y), "Total (i.e. Null); ", 
        x$df.residual, "Residual\n")
    if (!xg$model$scale.fix) {
        cat("\nScale Link:                  ", xg$model$sca.link)
        cat("\nEstimated Scale Parameters:  ")
        print(as.numeric(unclass(xg$gamma)), digits = digits)
    } else cat("\nScale is fixed.\n")
    cat("\nCorrelation:  Structure =", xg$model$corstr, " ")
    if (pmatch(xg$model$corstr, "independence", 0) == 0) {
        cat("  Link =", xg$model$cor.link, "\n")
        cat("Estimated Correlation Parameters:\n")
        print(unclass(xg$alpha), digits = digits)
    }
    cat("\nNumber of clusters:  ", length(xg$clusz), "  Maximum cluster size:", 
        max(xg$clusz), "\n")
    if (nzchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")

    cat("\n")
    invisible(x)
}
