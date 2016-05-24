getCallHL <- function(object) {
  if (is.null(call <- attr(object,"corrHLfitcall"))) {
    if (is.null(call <- attr(object,"HLCorcall"))) call <- getCall(object)
  }
  call
}

update.HLfit <- function (object, formula., ..., evaluate = TRUE) 
{
  if (is.null(call <- getCallHL(object))) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) {
    predictor <- formula(object) ## formula.default gets formula from $call, not from $predictor
    if (inherits(predictor,"predictor")) {
      form <- update.formula(attr(predictor,"oriFormula"),formula.) ## LOSES ALL ATTRIBUTES 
    } else  form <- update.formula(predictor,formula.) 
    ## !!!! FR->FR does not handle etaFix$beta !!!!
    if (! is.null(findOffset(formula.))) {off <- NULL} else { off <- attr(predictor,"offsetObj")$total }
    predArgs <- list(formula=form,
                     LMatrix=attr(predictor,"LMatrix"),
                     AMatrix=attr(predictor,"AMatrix"),
                     ZALMatrix=attr(predictor,"ZALMatrix"), ## all again from $call, not from $predictor
                     offset=off)
    ## attributes BinDenForm and oriFormula will be reconstructed:
    call$formula <- do.call("Predictor",predArgs) ## reconstructs oriFormula... otherwise we have a predictor without it...
  }
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call))) ## which to replace and which to add to the call
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]] ## replace
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing]) ## add
      call <- as.call(call)
    }
  }
  if (evaluate) 
    eval(call, parent.frame())
  else call
}


#`update.HLfit` <- function(object,formula.,...) {update.HL(object=object,formula.=formula.,...)}
#`update.HLCor` <- function(object,formula.,...) {update.HL(object=object,formula.=formula.,...)}

