## initialiser offset a O ?? pas actuellement pcq bcp de tests qu'il est nul
Predictor <- function (formula, offset=NULL, LMatrix = NULL,  AMatrix = NULL, ZALMatrix = NULL) {
  if (inherits(formula,"predictor")) {
    pastefrom("Do not call 'Predictor' on a predictor object.",prefix="(!) From ")
  }
  oriFormula <- formula
  if (substr((as.character(formula[2])),1,5)=="cbind") { 
    ## FR->FR e.g. strsplit("cbind(npos,ntot-npos)","[(,)-]") gives a list which element [[1]] is "cbind" "npos"  "ntot"  "npos"
    positives <- strsplit(strsplit(as.character(formula[2]),"\\(")[[1]][2],",")[[1]][1] ## names of variables
    negatives <- strsplit(strsplit(as.character(formula[2]),",")[[1]][2],"\\)")[[1]][1] ## names of variables
    if (length(positives)>1 && length(negatives)>1) {
      stop("For binomial data, please use cbind(<pos>,<neg>) where at least one of <pos> and <neg> is a variable from the data frame")
    } ## because problem bootstrap otherwise...
    ## cbind syntax for binomial model => is converted to alternative syntax; cbind would have failed in HLframes <10/04/2014 
    ## HLframes was modified on 10/04/2014 but the present code as not been revised following that change whichmay not be sufficient
    formula <- as.formula(paste(positives,"~",as.character(formula[3])))
    BinDenForm <- paste(positives,"+",negatives)       
  } else {BinDenForm <-NULL}
  if ( ! ( is.null(AMatrix) || is.list(AMatrix)) ) {
    ## assumes any of a number of matrix classes. Not clear how to test compactly for the soup of Matrix classes
    AMatrix <- list(AMatrix) ## further code expects a list (should ultimately be the same for all matrices...) 
  }
  if ( ! is.null(offset) ) {
    offterm <- findOffset(formula)
    if ( ! is.null(offterm) ) stop("in 'Predictor', offset should be given EITHER as $formula term OR as $offset element")
  } 
  if ( ! is.null(LMatrix)) {
    if ( ! is.list(LMatrix)) LMatrix <- list(dummyid=LMatrix)
    LMatrix <- lapply(LMatrix, function(lmatrix) {
      ranefs <- attr(lmatrix,"ranefs")
      if (is.null(ranefs)) {
        ranefs <- parseBars(formula) ## FR->FR or oriformula ???
      } else {
        ranefs <- unlist(lapply(ranefs,function(term) {parseBars(as.formula(paste("bla~",term)))}))      
      }
      attr(lmatrix,"ranefs") <- ranefs
      attr(lmatrix,"userLfixed") <- TRUE ## else remains NULL...
      lmatrix
    })
  }
  res <- formula
  attr(res,"oriFormula") <- oriFormula
  attr(res,"ZALMatrix") <- ZALMatrix
  attr(res,"AMatrix") <- AMatrix
  attr(res,"BinDenForm") <- BinDenForm
  attr(res,"LMatrix") <- LMatrix
  attr(res,"offsetObj") <- list(offsetArg=offset,nonZeroInfo= !is.null(offset))
  class(res) <- c("predictor",class(res))
  return(res)
}

## modified from getAnywhere(print.formula)
print.predictor <- function (x, showAttr = FALSE, ...) 
{
  .x <- x ## .x is original, returned invisibly
  if (! showAttr) attributes(x) <- NULL
  print.default(unclass(x), ...) ## from print.formula...
  invisible(.x)
}



