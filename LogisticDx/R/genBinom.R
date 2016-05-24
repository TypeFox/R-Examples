#' @name genBinom
#' @rdname genBinom
#'
#' @title Generate data for binomial regression
#' @description
#' Generates a \code{data.frame} or \code{data.table}
#' with a binary outcome, and a logistic model to
#' describe it.
#'
#' @param b The number of \bold{b}inomial variables
#' (the number of predictors
#' which are binary).
#' \cr
#' These are limited to \eqn{0} or \eqn{1}.
#' @param f The number of \bold{f}actor predictors.
#' \cr
#' The number of predictors
#' which are \code{factor}s.
#' @param c The number of \bold{c}ontinuous predictors.
#' \cr
#' the number of predictors which are
#' continuous.
#' @param n The \bold{n}umber of observations (rows) in the
#' \code{data.frame} or \code{data.table}.
#' @param nlf The \bold{n}umber of \bold{l}evels
#' in a \bold{f}actor.
#' @param pb The \bold{p}robability for \bold{b}inomnial
#' predictors:
#' \cr
#' the probability of binomial predictors being \eqn{=1}.
#' \cr
#' E.g. if \code{pb=0.3}, \eqn{30\%} will be \eqn{1}s,
#'  \eqn{70\%} will be \eqn{0}s
#' @param rc The \bold{r}atio for \bold{c}ontinuous variables.
#' \cr
#' The ratio of levels of
#' continuous variables to the total number of
#' observations \code{n}.
#' \cr
#' E.g. if \code{rc=0.8} and \code{n=100},
#' it will be in the range \eqn{1} to \eqn{80}.
#' @param py The \bold{r}atio for \bold{y},
#' the ratio of \eqn{1}s to the total number of observations
#' for the
#' binomial predictors.
#' \cr
#' E.g. if \code{ry=0.5},
#' 50\% will be \eqn{1}s, \eqn{50\%} will be \eqn{0}s.
#' @param asFactor If \code{asFactor=TRUE} (the default),
#' predictors given as \code{factor}s
#' will be converted to \code{factor}s
#' in the data frame before the model
#' is fit.
#' @param model If \code{model=TRUE},
#' will also return a model fitted with
#' \code{stats::glm} or \code{speedglm::speedglm}
#' @param timelim function will timeout after \code{timelim} secs.
#' This is present to prevent duplication of rows.
#' @param speedglm If \code{speedglm=TRUE},
#' return a model fitted with \code{speedglm}
#' instead of \code{glm}. See:
#' ?speedglm::speedglm
#' @return If \code{model=TRUE}: a list with the following values:
#'  \item{df or dt}{A \code{data.frame} (for \code{genBinomDf})
#' or \code{data.table}
#' (for \code{genBinomDt}).
#' \cr
#' Predictors are labelled \eqn{x1, x2, ..., xn}.
#' \cr
#' The response is \eqn{y}.
#' \cr
#' Rows represent to \eqn{n} observations}
#'  \item{model}{A model fit with \code{stats::glm}
#' or \code{speedglm::speedglm}}
#' If \code{model=FALSE} a \code{data.frame}
#' or \code{data.table} as above.
#'
#' @note  \code{genBinomDt} is faster
#' and more efficient for large datasets.
#' \cr \cr
#' Using \code{asFactor=TRUE} with \code{factor}s
#' which have a large number of \code{levels}
#' (e.g. \code{nlf > 30})
#' on large datasets (e.g. \eqn{n > 1000})
#' can cause fitting to be excessively slow.
#' @keywords datagen
NULL
#' @rdname genBinom
#' @export genBinomDf
#'
#' @examples
#' set.seed(1)
#' genBinomDf(speedglm=TRUE)
#'
genBinomDf <- function(b=2L, f=2L, c=1L, n=20L,
                       nlf=3L,
                       pb=0.5, rc=0.8, py=0.5,
                       asFactor=TRUE,
                       model=FALSE,
                       timelim=5,
                       speedglm=FALSE) {
    if ( (pb|rc|py) > 1 ) stop("Ratios should be <1")
    if (nlf <= 2) stop("Factors should have at least 3 levels")
    sbcf1 <- sum(c(b, f, c))
    if (sbcf1 <= 0) stop("Need at least one predictor")
    if (sbcf1 >= n) stop("Need n to be larger for this no. of predictors")
    ## prevent taking more than (timelimit)
    setTimeLimit(elapsed=timelim, transient=TRUE)
    ## define frame to hold values
    df1 <- as.data.frame(matrix(0L, ncol=(sbcf1 + 1L), nrow=n))
    xnam1 <- paste("x", 1:(ncol(df1)-1), sep="")
    colnames(df1) <- c(xnam1, "y")
    ## repeat until no dupicated columns
    repeat{
        if(!b==0){
            df1[, 1L:b] <- sample(x=c(0L, 1L),
                                  size=b*n,
                                  replace=TRUE,
                                  prob=c(pb, 1-pb))}
        if(!f==0){
            if(asFactor){
                df1[, (b+1L):(b+f)] <- as.factor(
                    sample(x=seq.int(1L, nlf),
                           size=f*n,
                           replace=TRUE))
            } else {
                 df1[, (b+1L):(b+f)] <- sample(x=seq.int(1L, nlf),
                                               size=f*n,
                                               replace=TRUE)
             }
        }
        if(!c==0){
            df1[, (b+f+1L):(b+f+c)] <- sample(x=seq.int(1, (rc*n)),
                                              size=c*n,
                                              replace=TRUE)}
        df1[, ncol(df1)] <- sample(x=c(0L, 1L),
                                   size=n,
                                   replace=TRUE,
                                   prob=c(py, 1-py))
### check no duplicate columns (prevent overfitting)
        if (!all(duplicated(df1))){break}
    }
###
    if(!model) return(df1)
### generate formula
    fmla <- stats::as.formula(paste("y ~ ", paste(xnam1, collapse= "+")))
### use speeedglm?
    if(speedglm){
        f1 <- speedglm::speedglm(formula=fmla,
                                 family=stats::binomial,
                                 data=df1)
    } else {
        f1 <- stats::glm(formula=fmla,
                         family=stats::binomial,
                         data=df1)
    }
    res <- list(df=df1,
                model=f1)
    return(res)
}
###
###----------------------------------------
###
#' @rdname genBinom
#' @export genBinomDt
#'
#' @examples
#' genBinomDt(b=0, c=2, n=100L, rc=0.7, model=FALSE)
#'
genBinomDt <- function(b=2L, f=2L, c=1L, n=20L,
                       nlf=3L,
                       pb=0.5, rc=0.8, py=0.5,
                       asFactor=TRUE,
                       model=FALSE,
                       timelim=5,
                       speedglm=FALSE) {
    if ( (pb|rc|py) > 1 ) stop("Ratios should be <1")
    if (nlf <= 2) stop("Factors should have at least 3 levels")
    sbcf1 <- sum(c(b, f, c))
    if (sbcf1 <= 0) stop("Need at least one predictor")
    if (sbcf1 >= n) stop("Need n to be larger for this no. predictors")
### prevent taking more than (timelimit)
    setTimeLimit(elapsed=timelim, transient=TRUE)
### define frame to hold values
    dt1 <- data.table::data.table(matrix(0, ncol=(sbcf1+1L), nrow=n))
    xnam1 <- paste("x", 1:(ncol(dt1)-1), sep="")
    data.table::setnames(dt1, c(xnam1, "y"))
    repeat{
        if(!b==0){
            for (i in 1L:b){
                data.table::set(dt1, j=i, value=sample(x=c(0L, 1L),
                                          size=n,
                                          replace=TRUE,
                                          prob=c(pb, 1-pb)
                                          )
                                )
            }
        }
        if(!f==0){
            for (i in (b+1L) : (b+f)){
                if (asFactor){
                    data.table::set(dt1, j=i, value=as.factor(sample(x=seq.int(1L, nlf),
                                              size=n,
                                              replace=TRUE))
                                    )
                } else {
                    data.table::set(dt1, j=i, value=sample(x=seq.int(1L, nlf),
                                              size=n,
                                              replace=TRUE)
                                    )
                }
            }
        }
        if(!c==0){
            for (i in (b+f+1L) : (b+f+c)){
                data.table::set(dt1, j=i, value=sample(x=seq.int(1, (rc*n)),
                                          size=n,
                                          replace=TRUE)
                                )
            }
        }
        data.table::set(dt1, j=as.integer(sbcf1+1L), value= sample(x=c(0L, 1L),
                                         size=n,
                                         replace=TRUE,
                                         prob=c(py, 1-py)
                                         ))
### check no duplicate columns (prevent overfitting)
        if (!all(duplicated(dt1))){break}
    }
###
    if(!model) return(dt1)
### generate model
### generate formula
    fmla <- stats::as.formula(paste("y ~ ", paste(xnam1, collapse= "+")  ))
### use speeedglm?
    if(speedglm){
        f1 <- speedglm::speedglm(formula=fmla, family=stats::binomial, data=dt1)
    } else {
        f1 <- stats::glm(formula=fmla, family=stats::binomial, data=dt1)
    }
    res <- list(dt=dt1,
                model=f1)
    return(res)
}
