#'Generate a Distance Matrix
#'
#'The gendistance function creates an \eqn{(N+K)}x\eqn{(N+K)} distance matrix
#'from an \eqn{N}x\eqn{P} covariates matrix, where \eqn{N} is the number
#'of subjects, \eqn{P} the number of covariates, and \eqn{K} the number of
#'phantom subjects requested (see \code{ndiscard} option). Provided the
#'covariates' covariance matrix is invertible, the distances computed are
#'Mahalanobis distances, or if covariate weights are provided, Reweighted
#'Mahalanobis distances (see \code{weights} option and Greevy, et al.,
#'Pharmacoepidemiology and Drug Safety 2012).
#'
#'Given a data.frame of covariates, generate a distance matrix.  Missing values
#'are imputed with \code{\link{fill.missing}}.  For each column with missing
#'data, a missingness indicator column will be added.  Phantoms are fake
#'elements that perfectly match all elements.  They can be used to discard a
#'certain number of elements.
#'
#'@aliases gendistance gendistance,data.frame-method
#'@param covariate A data.frame object, containing the covariates of the data
#'set.
#'@param idcol An integer or column name, providing the index of the column
#'containing row ID's.
#'@param weights A numeric vector, the length should match the number of
#'columns.  This value determines how much weight is given to each column when
#'generating the distance matrix.
#'@param prevent A vector of integers or column names, providing the index of
#'columns that should be used to prevent matches.  When generating the distance
#'matrix, elements that match on these columns are given a maximum distance.
#'@param force An integer or column name, providing the index of the column
#'containing information used to force pairs to match.
#'@param rankcols A vector of integers or column names, providing the index of
#'columns that should have the rank function applied to them before generating
#'the distance matrix.
#'@param missing.weight A numeric value, or vector, used to generate the weight
#'of missingness indicator columns.  Missingness indicator columns are created
#'if there is missing data within the data set.  Defaults to 0.1.  If a single
#'value is supplied, weights are generating by multiplying this by the original
#'columns' weight.  If a vector is supplied, it's length should match the
#'number of columns with missing data, and the weight is used as is.
#'@param ndiscard An integer, providing the number of elements that should be
#'allowed to match phantom values.  The default value is 0.
#'@param singular.method A character string, indicating the function to use
#'when encountering a singular matrix.  By default, \code{\link{solve}} is called.
#'The alternative is to call \code{\link{ginv}} from the \pkg{MASS} package.
#'@param talisman An integer or column name, providing location of talisman column.
#'The talisman column should only contains values of 0 and 1.  Records with zero
#'will match phantoms perfectly, while other records will match phantoms at max distance.
#'@param prevent.res.match An integer or column name, providing location of the column
#'containing assigned treatment groups.  This is useful in some settings, such as
#'trickle-in randomized trials.  When set, non-NA values from this column are
#'replaced with the value 1.  This prevents records with previously assigned
#'treatments (the \sQuote{reservior}) from matching each other.
#'@param \dots Additional arguments, not used at this time.
#'@return a list object with several elements
#'
#'  \item{dist}{generated distance matrix}
#'
#'  \item{cov}{covariate matrix used to generate distances}
#'
#'  \item{ignored}{ignored columns from original covariate matrix}
#'
#'  \item{weights}{weights applied to each column in covariate matrix}
#'
#'  \item{prevent}{columns used to prevent matches}
#'
#'  \item{mates}{index of rows that should be forced to match}
#'
#'  \item{rankcols}{index of columns that should use rank}
#'
#'  \item{missing.weight}{weight to apply to missingness indicator columns}
#'
#'  \item{ndiscard}{number of elements that will match phantoms}
#'@exportMethod gendistance
#'@author Cole Beck
#'@seealso \code{\link{distancematrix}}
#'@examples
#'
#'set.seed(1)
#'df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
#'# add some missing data
#'df[sample(seq_len(nrow(df)), ceiling(nrow(df)*0.1)), 2] <- NA
#'df.dist <- gendistance(df, idcol=1, ndiscard=2)
#'# up-weight the second column
#'df.weighted <- gendistance(df, idcol=1, weights=c(1,2,1), ndiscard=2, missing.weight=0.25)
#'df[,3] <- df[,2]*2
#'df.sing.solve <- gendistance(df, idcol=1, ndiscard=2)
#'df.sing.ginv <- gendistance(df, idcol=1, ndiscard=2, singular.method="ginv")
#'

setGeneric("gendistance", function(covariate, idcol=NULL, weights=NULL,
           prevent=NULL, force=NULL, rankcols=NULL, missing.weight=0.1,
           ndiscard=0, singular.method='solve', talisman=NULL,
           prevent.res.match=NULL, ...) standardGeneric("gendistance"))
setMethod("gendistance", "data.frame", function(covariate, idcol=NULL,
          weights=NULL, prevent=NULL, force=NULL, rankcols=NULL,
          missing.weight=0.1, ndiscard=0, singular.method='solve',
          talisman=NULL, prevent.res.match=NULL, ...) {
    nr <- nrow(covariate)
    nc <- ncol(covariate)
    stopifnot(nr > 0)
    myrownames <- seq_len(nr)
    mycolnames <- names(covariate)
    mateIDs <- integer(0)
    bad.data <- NULL

    # columns that aren't numeric should be marked bad
    badcol <- which(sapply(1:nc, FUN=function(x) suppressWarnings(!is.numeric(covariate[,x]))))
    # if all values in a column are the same, mark as bad column
    badcol <- union(badcol, which(sapply(1:nc, FUN=function(x) length(unique(covariate[,x])) == 1L)))

    if(!is.null(idcol) && length(idcol) == 1L) {
        if(is.character(idcol)) idcol <- match(idcol, mycolnames)
        idcol <- as.integer(idcol)
        if(!is.na(idcol) && idcol > 0 && idcol <= nc) {
            myrownames <- as.character(covariate[,idcol])
            row.names(covariate) <- myrownames
            badcol <- union(badcol, idcol)
        }
    }

    # warning for factor variables
    factorcols <- setdiff(which(sapply(covariate, is.factor)), idcol)
    if(length(factorcols)) {
      warning(sprintf("consider converting factor variables [%s] before calling gendistance", paste(mycolnames[factorcols], collapse=', ')))
    }

    if(is.null(weights)) {
        weights <- rep(1, nc)
    } else {
        weights <- as.numeric(weights)
        weights <- ifelse(is.na(weights) | weights < 0, 0, weights)
        weights <- c(weights, numeric(nc-length(weights)))
    }
    if(all(weights == 0)) stop("At least one column weight must be non-zero")

    if(!is.null(prevent)) {
        if(is.character(prevent)) prevent <- match(prevent, mycolnames)
        prevent <- as.integer(prevent)
        # in the past, prevent was ignored if found in badcol
        # but it could make sense to prevent unequal character strings
        prevent <- na.omit(ifelse(prevent < 1 | prevent > nc, NA, prevent))
        if(length(prevent) >= 1L) {
            badcol <- union(badcol, prevent)
            # convert column index to column name
            prevent <- mycolnames[prevent]
        }
    }

    tal.fail <- numeric()
    if(!is.null(talisman) && length(talisman) == 1L) {
        if(is.character(talisman)) talisman <- match(talisman, mycolnames)
        talisman <- as.integer(talisman)
        if(!is.na(talisman) && talisman > 0 && talisman <= nc) {
            badcol <- union(badcol, talisman)
            # non-zero values fail talisman
            tal.fail <- which(covariate[,talisman] != 0)
        }
    }

    if(!is.null(prevent.res.match) && length(prevent.res.match) == 1L) {
        if(is.character(prevent.res.match)) prevent.res.match <- match(prevent.res.match, mycolnames)
        prevent.res.match <- as.integer(prevent.res.match)
        if(!is.na(prevent.res.match) && prevent.res.match > 0 && prevent.res.match <= nc) {
            badcol <- union(badcol, prevent.res.match)
            prevent <- union(prevent, mycolnames[prevent.res.match])
            # set non-NA to value of 1
            nna.ix <- which(!is.na(covariate[,prevent.res.match]))
            if(length(nna.ix) == nr) stop(sprintf("invalid prevent.res.match argument, all values in column %s are non-missing", mycolnames[prevent.res.match]))
            covariate[nna.ix, prevent.res.match] <- 1
        }
    }

    if(!is.null(force)) {
        if(is.character(force)) force <- match(force, mycolnames)
        force <- as.integer(force)
        force <- setdiff(na.omit(ifelse(force < 1 | force > nc, NA, force)), badcol)
        # while this could adapt to a vector of columns, only accept one column
        if(length(force) == 1L) {
            badcol <- union(badcol, force)
            mateIDs <- as.integer(covariate[,force])
            # ensure ids are valid and not duplicated
            mateIDs <- ifelse(mateIDs < 1 | mateIDs > length(mateIDs) | duplicated(mateIDs), NA, mateIDs)
            # ensure ids are reflexive (1->2, 2->1)
            mateIDs <- ifelse(mateIDs[mateIDs] == seq_along(mateIDs), mateIDs, NA)
        }
    }
    badcol <- badcol[order(badcol)]

    # remove bad columns
    if(length(badcol) >= 1L) {
        weights <- weights[-badcol]
        # selecting one column from a data.frame creates a vector; using drop keeps a data.frame
        bad.data <- covariate[,badcol, drop=FALSE]
        covariate <- covariate[,-badcol, drop=FALSE]
    }
    # validate rankcols
    if(!is.null(rankcols)) {
        if(is.character(rankcols)) rankcols <- match(rankcols, mycolnames)
        rankcols <- as.integer(rankcols)
        rankcols <- setdiff(na.omit(ifelse(rankcols < 1 | rankcols > nc, NA, rankcols)), badcol)
        # bad columns have been removed, so if there are any rank columns, their column index has changed
        if(length(rankcols) >= 1L) {
            rankcols <- sapply(rankcols, FUN=function(y) { y - sum((y > badcol)*1) })
        }
    }
    # validate missing.weight and ndiscard
    if(!is.numeric(missing.weight) || (length(missing.weight) == 1 && missing.weight < 0)) missing.weight <- 0.1
    if(!is.numeric(ndiscard) || ndiscard < 0 || (nr - ndiscard) < 2) ndiscard <- 0

    # impute any missing values
    if(any(is.na(covariate))) {
        orig.colnames <- names(covariate)
        covariate <- fill.missing(covariate)
        new.colnames <- names(covariate)
        # calculate new weights
        weight.lookups <- sapply(sub(".missing", "", new.colnames[(length(orig.colnames)+1):length(new.colnames)]), FUN=function(y) { which(orig.colnames == y)  })
        if(length(missing.weight) > 1L) {
            if(length(missing.weight) != length(weight.lookups)) {
                missing.weight <- rep(missing.weight, length.out=length(weight.lookups))
                warning('the number of elements in missing.weight does not equal the number of variables with missingness - weights will be recycled as necessary')
            }
            # weight lookups are ignored when vector of missing.weights is provided
            weights <- append(weights, missing.weight)
        } else {
            weights <- append(weights, weights[weight.lookups]*missing.weight)
        }
        if(length(rankcols) >= 1L) {
            rankcols <- append(rankcols, (length(orig.colnames)+1):length(new.colnames))
        } else {
            rankcols <- (length(orig.colnames)+1):length(new.colnames)
        }
    }

    # Define your matrix of covariates covariate
    X <- as.matrix(covariate)
    if(nrow(X) < 2) stop('covariates data.frame must have at least two rows')
    if(ncol(X) < 1) stop('covariates data.frame must have at least one column')
    if(length(rankcols) >= 1L) {
        for(i in rankcols) {
            X[,i] <- rank(covariate[,i])
        }
    }

    # Create the covariance matrix and invert it
    X.cov <- cov.wt(X)
    # use pseudo-inverse if matrix is singular
    Sinv <- tryCatch(solve(X.cov$cov), error=function(e) { warning(e[[1]]); NULL })
    if(is.null(Sinv)) {
        # options for singular.method: [solve, ginv]
        if(is.null(singular.method) || !is.character(singular.method)) singular.method <- 'solve'
        if(singular.method=='ginv') {
            Sinv <- tryCatch(ginv(X.cov$cov), error=function(e) { warning(e[[1]]); NULL})
            if(!is.null(Sinv)) {
                warning("The covariance matrix is not invertible. The Moore-Penrose pseudoinverse (generalized inverse) was used to compute distances.")
            }
        }
        if(is.null(Sinv)) {
            Sinv <- solve(diag(diag(X.cov$cov)))
            warning("The covariance matrix is not invertible and/or the option for Euclidean distances was selected. Standardized Euclidean distances were used to compute distances.")
        }
    }
    # prevent negative distances
    for(i in seq_len(nrow(Sinv))) {
        Sinv[i,] <- Sinv[i,]*weights[i]
        Sinv[,i] <- Sinv[,i]*weights[i]
    }

    # Create the distance matrix mdists
    mdists <- sapply(seq_len(nr), FUN=function(x) { y <- X - X[rep(x, nr),]; rowSums(y %*% Sinv * y) })
    # Define a function to create the distance matrix using Sinv and X, previously done this way
    # mdistmaker <- function(row1, row2) { t(X[row1,]-X[row2,]) %*% Sinv %*% (X[row1,]-X[row2,]) }
    # mdists <- sapply(seq_len(nr), FUN=function(x) mapply(mdistmaker, x, seq_len(nr)))
    # take the square root
    mdists <- sqrt(mdists)

    maxval <- Inf
    # add back row names
    dimnames(mdists) <- list(myrownames, myrownames)

    # penalize "prevent" columns with matching values by setting distance to maxval
    # prevent is a vector of column names found in bad.data
    if(length(prevent) >= 1L) {
        indeces <- match(prevent, names(bad.data))
        for(colnum in indeces) {
            for(rownum in seq_len(nr)) {
                # get the row number of all rows that have the same value in the given column
                eqrows <- setdiff(which(bad.data[rownum,colnum] == bad.data[,colnum]), colnum)
                if(length(eqrows) > 0L) {
                    mdists[rownum,eqrows] <- maxval
                    mdists[eqrows,rownum] <- maxval
                }
            }
        }
    }

    # need to add phantom rows/columns of zero distance
    GROUPS <- 2
    nphantoms <- ndiscard + ((GROUPS - (nr - ndiscard) %% GROUPS) %% GROUPS)
    if(nphantoms > 0L) {
        mdists <- make.phantoms(mdists, nphantoms, maxval=maxval)
        # failed talisman receive max distance
        if(length(tal.fail)) {
            pcols <- seq(nr + 1, nr + nphantoms)
            mdists[tal.fail, pcols] <- mdists[pcols, tal.fail] <- maxval
        }
    }

    # forced matches/mates should receive a distance of zero -- non-matches receive maxval
    if(length(mateIDs)) {
        for(i in seq_along(mateIDs)) {
            j <- mateIDs[i]
            if(!is.na(j)) {
                mdists[i,] <- mdists[,i] <- maxval
                mdists[i,j] <- mdists[j,i] <- 0
            }
        }
    }

    diag(mdists) <- maxval
    # convert matrix to data frame
    mdists <- as.data.frame(mdists)

    list(dist=mdists, cov=covariate, ignored=bad.data, weights=weights, prevent=prevent, mates=mateIDs, rankcols=rankcols, missing.weight=missing.weight, ndiscard=nphantoms)
})
