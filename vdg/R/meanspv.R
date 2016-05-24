#' Compute Mean Spherical SPV
#' 
#' Computes the matrix of spherical region moments for a given model formula and a vector of radii, and uses this to
#' calculate the mean spherical SPV for each of the radii. The function \code{expmat} calculates
#' the matrix containing the exponents of each model factor within each model term as columns. 
#' Only simple formulae are allowed for. Only products of terms should be included in
#' calls to \code{\link{I}}. The power operator \code{\link{^}} should be used instead
#' of \code{\link{sqrt}}. Models should contain only monomial terms. 
#' 
#' @param formula model formula
#' @param radii numeric vector or radii at which to calculate the matrix of spherical region moments
#' @param FtF.inv inverse of F'F, where F is the design matrix
#' @param n integer giving the number of design runs
# ' @return a matrix of variables by terms containing the exponents of the varaibles in every term
#' @examples
#' f1 <- formula(~ x1*x2)
#' expmat(f1)
#' f2 <- update(f1, ~ . + I(x1^2) + I(x2^2))
#' expmat(f2)
#' f3 <- update(f2, ~ . + I(x2^0.4))
#' expmat(f3)
#' f4 <- update(f3, ~ . + I(x1^2):I(x2^2))
#' expmat(f4)
#' f5 <- update(f4, ~ . + I(x1^3*x2^0.5))
#' expmat(f5)
meanspv <- function(formula, radii, FtF.inv, n){
  expmat <- expmat(formula = formula)
  nterms <- ncol(expmat)
  
  sigfun <- function(delta, r){
    if(any(delta %% 2 != 0)) return(rep(0, length(r)))
#     if(all(delta == 0)) return(rep(1, length(r)))
    m <- length(delta)
    sdelta <- sum(delta)
    out <- r^sdelta * gamma(m/2) * prod(gamma((delta + 1)/2)) / (pi^(m/2) * 
           gamma((sdelta + m)/2))
    return(out)
  }
  
  lowers <- t(apply(combn(nterms, 2), 2, function(x) sigfun(expmat[,x[1]] + 
                                                     expmat[,x[2]], r = radii)))
  diags <- t(apply(2*expmat, 2, sigfun, r = radii))
  smom <- mapply(as.data.frame(lowers), as.data.frame(diags), 
            FUN = function(x, y){tmp <- matrix(0, length(y), length(y)); 
                              tmp[lower.tri(tmp)] <- x;
                              tmp <- tmp + t(tmp);
                              diag(tmp) <- y; 
                              dimnames(tmp) <- rep(list(colnames(expmat)), 2);
                              return(tmp)
            }, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  mspv <- mapply(FUN = function(x, y) n * sum(diag(x %*% y)), smom, 
                 list(FtF.inv))
  out <- structure(list(Radius = radii, SPV = mspv), 
                   class = c("meanspv", "list"))
  return(out)
}
#' @rdname meanspv
#' @export
expmat <- function(formula){
  vars <- all.vars(formula)
  terms <- terms(formula)
  factors <- attr(terms, "factors")
  test <- tnames <- colnames(factors)
  
  for(i in vars) test <- gsub(i, replacement = "", x = test)
  test <- gsub("[[:digit:]^()I.:* -]", "", test)
  if(!all(test == "")) 
    stop("Characters other that variable names and the regexp class '[[:digit:]^()I.:* -]' found in formula")
  
  varrows <- na.omit(match(vars, rownames(factors)))
  Irows <- seq_len(nrow(factors))[-varrows]
  Inames <- rownames(factors)[Irows]
  expmat <- factors[varrows, , drop = FALSE]
  expmat[expmat == 2] <- 1
  
  getExp <- function(term, var){
    if(!grepl(var, term)) return(0L)
    sterm <- gsub("^I\\(", "", term)
    sterm <- gsub(")$", "", sterm) 
    sterm <- unlist(strsplit(sterm, "*", fixed = TRUE))
    sterm <- sterm[grepl(var, sterm)]
    if(length(sterm) != 1) stop("A variable should only occur once in each term")
    tmatch <- match(var, sterm) 
    if(!is.na(tmatch)) if(tmatch == 1) return(1L)
    sterm <- gsub("[()]", "", sterm)
    sterm <- unlist(strsplit(sterm, "^", fixed = TRUE))
    if(length(sterm) != 2) stop("Unkown use of '^' operator")
    return(as.numeric(sterm[2]))
  }
  
  whichI <- match(Inames, tnames)
  for(i in whichI) for(j in seq_along(vars[varrows])) 
    expmat[j, i] <- getExp(tnames[i], vars[varrows][j])
  whichInt <- setdiff(seq_along(tnames), match(rownames(factors), tnames))
  whichIntI <- whichInt[grep("I\\(", tnames[whichInt])]
  nrInt <- length(whichIntI)
  if(nrInt){
    sint <- strsplit(tnames[whichIntI], ":", fixed = TRUE)
    intPos <- lapply(sint, match, Inames)
    for(i in seq_along(intPos)){
      index <- intPos[[i]][!is.na(intPos[[i]])]
      expmat[, whichIntI[i]] <- expmat[, whichIntI[i]] + 
        rowSums(expmat[, Inames[index], drop = FALSE])
    } 
  }
  if(attr(terms, "intercept")) expmat <- cbind("(Intercept)" = 0, expmat)
  return(expmat)
}