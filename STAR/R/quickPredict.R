##############################################################################
##############################################################################
##############################################################################
quickPredict <- function(object,
                         include=object$terms$labels[2],
                         se.fit=TRUE,
                         length.out,
                         otherTermsFct=median
                         )
{
######################################################################
### Function quickPredict: an "easy" interface to predict.ssanova and
### predict.ssanova0
### Designed to quickly compute the effect of a _single_ model term.
### This term can correspond to a single variable effect or to the
### interaction of two variables
### -----------------------------------------------------------------
### Arguments:
###  object: an object inheriting from ssanova and ssanova0 (gssanova
###          and gssanova0 objects are therefore suitable).
###  include: a character, the model term for which one wants the
###           the prediction.
###  se.fit: see predict.ssanova and predict.ssanova0
###  length.out: a positive integer, the number of points on the
###              variable / term definition domain at which the
###              "prediction" will be made.
###  otherTermsFct: a function used to set the value of the other
###                 model terms in the data frame fed to
###                 predict.ssanova / predict.ssanova0
### -------------------------------------------------------------------
### Value:
###  A "quickPredict" object. This object is a list with the following
###  components:
###  xx: a numeric vector with the values of the selected term at
###      which prediction was made. When an interaction term was
###      selected the values of the first variable are stored here.
###  yy: values of the second variable (for interaction terms) at
###      which prediction was made. NULL for none interaction terms.
###  include: the value of the argument with this name.
###  call: the matched call.
###  est.mean: a numeric vector, the estimated mean value of the term,
###            or a matrix for interaction terms.
###  est.sd: a numeric vector or NULL, the estimated standard
###          deviation of the term. It is NULL if argument "se.fit"
###          is set to FALSE. It is a matrix for interaction terms.
##########################################################################
  
  ## Make sure that object inherits from ssanova or ssanova0
  if (!inherits(object,c("ssanova","ssanova0")))
    stop("object should be a ssanova or a ssanova0 obbject.")

  allTerms <- object$terms$labels[-1]
  if (length(include) != 1) stop("include should be a character of length 1.")
  if (!(include %in% allTerms)) stop("include should be one of the model terms.")

  ## Find out if the term is an interaction term
  isInteraction <- !(include %in% names(object$mf))

  if (missing(length.out)) {
    if (isInteraction) length.out <- 101
    else length.out <- 501
  } ## End of conditional on missing(length.out)
  
  if (!isInteraction) {
    xx <- seq(from=object$terms[[include]]$rk$env$env$min,
              to=object$terms[[include]]$rk$env$env$max,
              length.out=length.out)
    yy <- NULL
    newdata <- data.frame(xx)
    names(newdata) <- include
    otherTerms <- allTerms[!(allTerms %in% include) &
                           (allTerms %in% names(object$mf))
                           ]
  } else {
    vNames <- strsplit(include,":")[[1]]
    xx <- seq(from=object$terms[[include]]$rk$env$rk[[1]]$env$min,
              to=object$terms[[include]]$rk$env$rk[[1]]$env$max,
              length.out=length.out)
    yy <- seq(from=object$terms[[include]]$rk$env$rk[[2]]$env$min,
              to=object$terms[[include]]$rk$env$rk[[2]]$env$max,
              length.out=length.out)
    newdata <- data.frame(rep(xx,length.out),
                          rep(yy,rep(length.out,length.out))
                          )
    names(newdata) <- vNames
    otherTerms <- allTerms[!(allTerms %in% vNames) &
                           (allTerms %in% names(object$mf))
                           ]
  } ## End of conditional on !isInteraction

  nbOtherTerms <- length(otherTerms)
  if (nbOtherTerms > 0) {
    otherVal <- sapply(otherTerms, function(n) otherTermsFct(object$mf[,n]))
    for (idx in 1:nbOtherTerms) {
      newdata[[otherTerms[idx]]] <- rep(otherVal[idx],dim(newdata)[1])
    } ## End of for loop on idx
  } ## End of conditional on nbOtherTerms > 0 

  est <- predict(object,newdata,se.fit=se.fit,include=include)

  result <- list(xx=xx,
                 yy=yy,
                 call=match.call(),
                 include=include)
  if (se.fit) {
    if (isInteraction) {
      result$est.mean <- matrix(est$fit,length.out,length.out)
      result$est.sd <- matrix(est$se.fit,length.out,length.out)
    } else {
      result$est.mean <- est$fit
      result$est.sd <- est$se.fit
    }                  
  } else {
    if (isInteraction) {
      result$est.mean <- matrix(est,length.out,length.out)
    } else {
      result$est.mean <- est
    }
    result$est.sd <- NULL
  } ## End of conditional on se.fit

  class(result) <- "quickPredict"
  result
}
##############################################################################
##############################################################################
##############################################################################
"%qp%" <- function(object,
                   include
                   )
{
#######################################################################
### Function %qp%
### Binary version of quickPredict
### ----------------------------------------------------------
### Arguments:
###  object: an object inheriting from ssanova and ssanova0 (gssanova
###          and gssanova0 objects are therefore suitable).
###  include: a character, the model term for which one wants the
###           the prediction.
### The other arguments of quickPredict are set to their default values.
### -------------------------------------------------------------------
### Value: see quickPredict.
#######################################################################
  x <- deparse(substitute(object))
  y <- deparse(substitute(include))
  cmd <- paste("quickPredict(",
               x,
               ",",
               y,
               ")",sep="")
  cmd <- parse(text=cmd)
  eval(cmd)

}
