#'Data Imputation
#'
#'The fill.missing function uses the \code{\link{transcan}} function from the
#'\pkg{Hmisc} package to impute values for the given data.frame.
#'
#'The fill.missing function will fill the missing values within a data.frame
#'with the values imputed with the \code{\link{transcan}} function.  An idcol may be
#'specified to prevent including the use of IDs in the imputation.  In addition
#'for every column that contains missing data, a new column will be attached to
#'the data.frame containing an indicator of missingness.  A "1" indicates that
#'the value was missing and has been imputed.
#'
#'@aliases fill.missing fill.missing,data.frame-method
#'@param x A data.frame object.  It should have missing values.
#'@param seed Seed provided for random-number generation.  Default value of
#'101.
#'@param simplify logical: whether to remove duplicate missingness columns.
#'@param idcol An integer value or character string.  Indicates the column
#'containing IDs, specified as column index or column name.  Defaults to "id",
#'or NA, when not found.
#'@param \dots Additional arguments, potentially passed to \code{\link{transcan}}.
#'@return data.frame with imputed values
#'@exportMethod fill.missing
#'@author Cole Beck
#'@seealso \code{\link{transcan}}
#'@examples
#'
#'set.seed(1)
#'df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
#'df[sample(seq_len(nrow(df)), ceiling(nrow(df)*0.1)), 2] <- NA
#'df <- fill.missing(df)
#'

setGeneric("fill.missing", function(x, seed=101, simplify=TRUE, idcol="id", ...) standardGeneric("fill.missing"))
setMethod("fill.missing", "data.frame", function(x, seed=101, simplify=TRUE, idcol="id", ...) {
    if(exists(".Random.seed", envir = .GlobalEnv)) {
        save.seed <- get(".Random.seed", envir= .GlobalEnv)
        on.exit(assign(".Random.seed", save.seed, envir = .GlobalEnv))
    } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    miss.cols <- which(apply(x, MARGIN=2, FUN=function(y) { xm<-is.na(y);(any(xm) & !all(xm))}))
    if(length(miss.cols) > 0) {
        missingness <- is.na(x[,miss.cols, drop=FALSE])+0
        dup <- duplicated(missingness, MARGIN=2)
        # simplify will remove duplicate "missing" columns
        if(simplify == TRUE && any(dup)) {
            missingness <- missingness[,!dup, drop=FALSE]
        }
        colnames(missingness) <- sprintf("%s.missing", colnames(missingness))
    } else {
      # no missing data
      return(x)
    }
    if(!is.numeric(seed)) seed <- 101
    set.seed(seed)
    cnames <- names(x)
    if(!is.null(idcol) && length(idcol) == 1) {
        if(is.character(idcol)) {
            idcol <- match(idcol, cnames)
        }
        if(!is.na(idcol) && is.numeric(idcol) && idcol > 0 && idcol <= ncol(x)) {
            cnames <- cnames[-idcol]
        }
    }
    fmla <- as.formula(paste("~", paste(cnames, collapse= "+")))
    # run until no errors
    fails <- 0
    while(fails < 10) {
        if(!is.null(tryCatch(imputer <- transcan(fmla, data=x, transformed=FALSE, pl=FALSE, pr=FALSE, imputed=TRUE, ...),  error = function(e){}))) break
        fails <- fails+1
    }
    if(fails >= 10) stop("cannot impute values")
    invisible(lapply(names(imputer$imputed), FUN=function(cname) {
        if(length(imputer$imputed[[cname]]) > 0L) {
            imputed.vals <- imputer$imputed[[cname]]
            x[names(imputed.vals), cname] <<- imputed.vals
        }
    }))
    cbind(x, missingness)
})
