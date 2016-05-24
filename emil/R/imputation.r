#' Support function for identifying missing values
#' 
#' @param data Fitting and testing data sets, as returned by
#'   \code{\link{pre_split}}.
#' @return Data frame containing row and column indices of missing values or
#'   \code{NULL} if the data doesn't contain any.
#' @examples
#' x <- as.matrix(iris[-5])
#' y <- iris$Species
#' x[sample(length(x), 10)] <- NA
#' cv <- resample("crossvalidation", y)
#' sets <- pre_split(x, y, cv[[1]])
#' sets <- pre_remove(sets, 3L)
#' na_index(sets)
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
na_index <- function(data){
    # `unname` is needed to avoid problems with duplicate names
    fit.na  <- which(unname(is.na(data$fit$x)), arr.ind=TRUE)
    if(any(fit.na)){
        fit.na <- data.frame(set = "fit", fit.na)
        fit.na$original_row <- index_fit(data$fold)[fit.na$row]
        fit.na$original_col <- data$feature_selection[fit.na$col]
    } else {
        fit.na <- NULL
    }
    test.na <- which(unname(is.na(data$test$x)), arr.ind=TRUE)
    if(any(test.na)){
        test.na <- data.frame(set = "test", test.na)
        test.na$original_row <- index_test(data$fold)[test.na$row]
        test.na$original_col <- data$feature_selection[test.na$col]
    } else {
        test.na <- NULL
    }
    rbind(fit.na, test.na)
}

#' Basic imputation
#'
#' This solution is optimized for the scenario that the dataset is very large
#' but only contains missing values in a small number of columns.
#'
#' @param data Fitting and test datasets, as returned by \code{\link{pre_split}}
#'   or any other standard pre-processing function.
#' @param fun Function for calculating imputation values. Should take a vector
#'  and return a scalar.
#' @param ... Sent to \code{fun}.
#' @return A pair of fitting and testing datasets.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
pre_impute <- function(data, fun, ...){
    na.ind <- na_index(data)
    if(is.null(na.ind)) return(data)

    na.feats <- unique(na.ind$col)
    m <- rep(NA, ncol(data$fit$x))
    m[na.feats] <- apply(data$fit$x[,na.feats,drop=FALSE], 2, fun, ...)

    na.ind %<>%
        split(na.ind$set) %>%
        lapply(function(x) as.matrix(x[,c("row", "col")]))
    if(!is.null(na.ind$fit))  data$fit$x[na.ind$fit]   <- m[na.ind$fit[,"col"]]
    if(!is.null(na.ind$test)) data$test$x[na.ind$test] <- m[na.ind$test[,"col"]]

    impute_failed <- intersect(which(is.na(m)), na.feats)
    if(any(impute_failed)){
        data <- pre_remove(data, impute_failed)
        warning(sprintf("Could not impute %i features.", length(impute_failed)))
    }
    data
}
#' @rdname pre_impute
#' @export
pre_impute_median <- function(data){
    pre_impute(data, fun=median, na.rm=TRUE)
}
#' @rdname pre_impute
#' @export
pre_impute_mean <- function(data){
    pre_impute(data, fun=mean, na.rm=TRUE)
}

#' Nearest neighbors imputation
#' 
#' Nearest neighbor methods needs to have a distance matrix of the dataset it works on.
#' When doing repeated model fittings on subsets of the entire dataset it is
#' unnecessary to recalculate it every time, therefore this function requires
#' the user to manually calculate it prior to resampling and supply it in a
#' wrapper function.
#'
#' Features with fewer than \code{k} non-missing values will be removed
#' automatically.
#' 
#' @param data Fitting and testing data sets, as returned by
#'   \code{\link{pre_split}}.
#' @param k Number of nearest neighbors to calculate mean from. Set to < 1 to
#'   specify a fraction.
#' @param distance_matrix A matrix, \code{\link{dist}} object or
#'   \code{"auto"}. Notice that \code{"auto"} will recalculate the distance
#'   matrix in each fold, which is only meaningful in case the features of
#'   \code{x} vary between folds. Otherwise you are just wasting time.
#'   
#' @examples
#' x <- iris[-5]
#' x[sample(nrow(x), 30), 3] <- NA
#' my.dist <- dist(x)
#' evaluate(modeling_procedure("lda"), x = x, y = iris$Species,
#'     pre_process = function(...){
#'         pre_split(...) %>% pre_impute_knn(k = 4, distance_matrix = my.dist)
#'     }
#' )
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
pre_impute_knn <- function(data, k=.05, distance_matrix){
    na_ind <- na_index(data)
    if(is.null(na_ind)) return(data)

    if(k < 1) k <- max(1, round(.05*length(index_fit(data$fold))))
    if(k > sum(data$fold > 0, na.rm=TRUE))
        stop("k is larger than number of fitting observations.")

    # Exclude features with fewer non-NAs than k
    na.count <- na_ind %>%
        filter_("set == 'fit'") %>%
        group_by_("col") %>%
        summarize_(n = "length(row)")
    impute_failed <- which(nrow(data$fit$x) - na.count$n < k)
 
    if(length(impute_failed) > 0){
        warning(sprintf("Could not impute %s.", 
            example_string("feature", names(impute_failed))))
        data <- pre_remove(data, impute_failed)
        na_ind %<>% filter_(~!col %in% impute_failed)
        na_ind$col <- na_ind$col - findInterval(na_ind$col, impute_failed)
    }

    # Check that the distance matrix is in order
    if(missing(distance_matrix))
        stop("You must supply a distance matrix, see `?pre_impute_knn` for details.")
    if(is.character(distance_matrix) && distance_matrix == "auto"){
        distance_matrix <- matrix(nrow = length(data$fold), ncol = length(data$fold))
        ind <- c(index_fit(data$fold, allow_oversample=FALSE),
                 index_test(data$fold))
        distance_matrix[ind, ind] <- as.matrix(dist(rbind(data$fit$x, data$test$x)))
    } else if(!is.matrix(distance_matrix)){
        distance_matrix <- as.matrix(distance_matrix)
    }
    if(any(length(data$fold) != dim(distance_matrix)))
        stop("Distance matrix does not match dataset.")

    # Perform the imputation
    diag(distance_matrix) <- NA
    NN <- as.data.frame(apply(distance_matrix[na_ind$original_row, index_fit(data$fold), drop=FALSE], 1, order))
    na_ind$fill <- mapply(function(i, col){
        x <- data$fit$x[i, col]
        mean(x[!is.na(x)][1:k])
    }, NN, na_ind$col)

    if(any(na_ind$set == "fit"))
        data$fit$x[as.matrix(na_ind[na_ind$set == "fit", c("row", "col")])] <-
            na_ind$fill[na_ind$set == "fit"]
    if(any(na_ind$set == "test"))
        data$test$x[as.matrix(na_ind[na_ind$set == "test", c("row", "col")])] <-
            na_ind$fill[na_ind$set == "test"]
    data
}

#' Regular imputation
#'
#' If you want to impute, build model and predict you should use
#' \code{\link{pre_impute_median}} or \code{\link{pre_impute_knn}}.
#' This function imputes using all observations
#' without caring about cross-validation folds.
#'
#' For additional information on the parameters see \code{\link{pre_impute_knn}}
#' and \code{\link{pre_impute}}.
#' 
#' @param x Dataset.
#' @param k Number of nearest neighbors to use.
#' @param distance_matrix Distance matrix.
#' @return An imputed matrix.
#' @examples
#' x <- matrix(rnorm(36), 6, 6)
#' x[sample(length(x), 5)] <- NA
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @name impute
#' @seealso \code{\link{emil}}, \code{\link{pre_process}},
#'   \code{\link{pre_impute_knn}}, \code{\link{pre_impute_median}}
{}
#' @examples
#' impute_knn(x)
#' @rdname impute
#' @export
impute_knn <- function(x, k=.05, distance_matrix="auto"){
    if(identical(distance_matrix, "auto"))
        distance_matrix <- dist(x)
    pre_split(x, y=NULL) %>%
        pre_impute_knn(k=k, distance_matrix=distance_matrix) %>%
        (function(data) data$fit$x)
}
#' @examples
#' impute_median(x)
#' @rdname impute
#' @export
impute_median <- function(x){
    pre_split(x, y=NULL) %>%
        pre_impute_median %>%
        (function(data) data$fit$x)
}

#' Impute a data frame
#' 
#' This function imputes each column of data frames univariately with different
#' functions depending on their class.
#' 
#' @param data Pre-processed data set with features in a data frame.
#' @param class_fun List of functions to use for imputating specific feature
#'   classes.
#' @param default_fun Function to use for imputation features of classes not
#'   listed in \code{class_fun}.
#' @param na.rm Whether to remove missing values.
#' @examples
#' x <- iris
#' x[sample(150, 3), 1] <- NA                                                                          
#' x[sample(150, 1), 3] <- NA
#' x[sample(150, 5), 5] <- NA
#' y <- gl(2, 75)
#' fold <- resample("holdout", y, nfold=1)[[1]]
#' data <- pre_split(x, y, fold) %>%
#'     pre_impute_df
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
pre_impute_df <- function(data,
        class_fun = list(numeric = function(x) median(x, na.rm=TRUE),
                         integer = function(x) median(x, na.rm=TRUE)),
        default_fun = function(x) mode(x, na.rm=TRUE, allow_multiple=FALSE)){
    stopifnot(inherits(data$fit$x, "data.frame"))
    na_ind <- na_index(data)
    m <- Map(function(x, any_na){
        if(any_na){
            fun_i <- match(class(x), names(class_fun))
            fun_i <- fun_i[!is.na(fun_i)]
            if(any(fun_i)) class_fun[[fun_i[1]]](x)
            else default_fun(x)
        } else NULL
    }, data$fit$x, seq_len(ncol(data$fit$x)) %in% unique(na_ind$col))

    impute_failed <- which(vapply(m, identical, logical(1), NA))
    if(length(impute_failed) > 0){
        warning(sprintf("Could not impute %s.", 
            example_string("feature", names(impute_failed))))
        data <- pre_remove(data, impute_failed)
        na_ind %<>% filter_(~!col %in% impute_failed)
        na_ind$col <- na_ind$col - findInterval(na_ind$col, impute_failed)
        m <- m[-impute_failed]
    }

    data$fit$x[TRUE] <- Map(function(x, value){
        if(is.null(value)) x
        else x[is.na(x)] <- value
    }, data$fit$x, m)
    data$test$x[TRUE] <- Map(function(x, value){
        if(is.null(value)) x
        else x[is.na(x)] <- value
    }, data$test$x, m)
    data
}

