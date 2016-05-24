##' Distributed gradient boosting based on the \pkg{mboost} package.
##' Gaussian, Binomial and Poisson families are currently supported.
##'
##' Generally gradient boosting offers more flexibility and better
##' predictive performance than random forests, but is usually not
##' used for large data sets because of its iterative
##' nature. \code{parboost} is designed to scale up component-wise
##' functional gradient boosting in a distributed memory environment
##' by splitting the observations into disjoint subsets, or
##' alternatively by bootstrapping the original data. Each cluster
##' node then fits a boosting model to its subset of the data. These
##' boosting models are combined in an ensemble, either with equal
##' weights, or by fitting a (penalized) regression model on the
##' predictions of the individual models on the complete data. The
##' motivation behind \code{parboost} is to offer a boosting
##' framework that is as easy to parallelize and thus scalable as
##' random forests.
##'
##' If you want to modify the boosting parameters, please take a look
##' at the \code{\link{mboost}} package documentation and pass the
##' appropriate parameters to \code{tree_control} and
##' \code{boost_control}.
##'
##' @title Distributed gradient boosting based on the \pkg{mboost} package.
##' @param cluster_object Cluster object from the \pkg{parallel} package to
##' carry out distributed computations.
##' @param mc.cores If not \code{NULL}, \code{parboost} uses mclapply for shared memory parallelism.
##' @param data A data frame containing the variables in the
##' model. It is recommended to use path_to_data instead for
##' IO efficiency. Defaults to NULL
##' @param path_to_data A string pointing to the location of the data. \code{parboost} assumes that the data is located at
##' the same location on every cluster node. This parameter is
##' ignored if you pass a data frame to the data argument.
##' @param data_import_function Function used to import data. Defaults
##' to \code{read.csv}. This parameter is ignored if you pass a data
##' frame to the data argument.
##' @param split_data String determening the way the data should be
##' split. \code{disjoint} splits the data into disjoint
##' subsets. \code{bootstrap} draws a bootstrap sample instead.
##' @param nsplits Integer determining the number of disjoint sets the
##' data should be split into. If \code{split_data} is set
##' to \code{bootstrap}, \code{nsplits} determines the number of
##' bootstrap samples.
##' @param preprocessing Optional preprocessing function to apply to
##' the data. This is useful if you cannot modify the data on the
##' cluster nodes.
##' @param seed Integer determining the random seed value for
##' reproducible results.
##' @param formula Formula to be passed to mboost.
##' @param baselearner Character string to determine the type of
##' baselearner to be used for boosting.
##' See \code{\link[mboost]{mboost}} for details.
##' @param family A string determining the family. Currently
##' gaussian, binomial and poisson are implemented.
##' @param control An object of type \code{boost_control} for
##' controlling \code{\link[mboost]{mboost}}.  See
##' \code{\link[mboost]{boost_control}} in the \code{\link{mboost}}
##' for details.
##' @param tree_controls Optional object of type \code{TreeControl}.
##' See \code{ctree_control} in the \code{party}
##' documentation for detailos. Used to set hyperparameters for tree
##' base learners.
##' @param cv Logical to activate crossvalidation to determine
##' \eqn{m_{stop}}. Defaults to \code{TRUE}.
##' @param cores_cv Integer determining the number of CPU cores used
##' for cross-validation on each node (or locally). Defaults to
##' maximum available using \code{\link[parallel]{detectCores}}.
##' @param folds Integer determening the number of folds used during
##' cross-validation on each cluster node. Defaults to 8.
##' It is computationally more efficient to set the
##' value of of folds to a multiple of the number of cores on each
##' cluster node.
##' @param stepsize_mstop Integer denoting the stepsize used during
##' cross-validation for tuning the value of \eqn{m_{stop}}.
##' @param postprocessing String to set the type of
##' postprocessing. Defaults to \code{"none"} for a simple average of
##' the ensemble components.
##' @return An object of type \code{parboost} with \code{print},
##' \code{summary}, \code{predict}, \code{coef} and \code{selected} methods.
##' @author Ronert Obst
##' @examples ## Run parboost on a cluster (example not run)
##' # data(friedman2)
##' # library(parallel)
##' # cl <- makeCluster(2)
##' # parboost_model <- parboost(cluster_object = cl, data = friedman2,
##' #                            nsplits = 2, formula = y ~ .,
##' #                            baselearner="bbs", postprocessing = "glm",
##' #                            control = boost_control(mstop=10))
##' # stopCluster(cl)
##' # print(parboost_model)
##' # summary(parboost_model)
##' # head(predict(parboost_model))
##' #
##' # ## Run parboost serially for testing/debugging purposes
##' # parboost_model <- parboost(data = friedman2, nsplits = 2, formula
##' # = y ~ ., baselearner="bbs", postprocessing = "glm", control =
##' # boost_control(mstop=10))
##' @export
parboost <- function(cluster_object=NULL, mc.cores=NULL, data=NULL, path_to_data="",
                     data_import_function = NULL,
                     split_data = c("disjoint", "bagging"), nsplits,
                     preprocessing = NULL, seed = NULL,
                     formula, baselearner = c("bbs", "bols", "btree", "bss", "bns"),
                     family = c("gaussian", "binomial", "poisson"),
                     control = boost_control(),
                     tree_controls = NULL,
                     cv = TRUE, cores_cv = detectCores(), folds = 8, stepsize_mstop = 1,
                     postprocessing = c("none", "glm", "lasso", "ridge", "elasticnet")) {

    if (is.null(data_import_function)) {
        ### determine classes on subsample to spead up data import
        data_import_function <- function(path) {
            sample_data <- read.csv(file = path, header = TRUE, nrows = 100)
            classes <- sapply(sample_data, class)
            read.csv(file = path, header = TRUE, colClasses = classes)
        }
    }

    if (is.null(data)) {
        export <- FALSE
    } else {
        export <- TRUE
    }

    postprocessing <- match.arg(postprocessing)
    split_data <- match.arg(split_data)
    family <- match.arg(family)

    if (is.character(baselearner)) {
        baselearner <- match.arg(baselearner)
    }

    if(is.null(data)) {
        data <- data_import_function(path_to_data)
    }

    if(!is.null(preprocessing)) {
        data <- preprocessing(data)
    }

    response <- data[, as.character(formula[2])]

    if(!is.null(seed)) set.seed(seed)
    if(split_data == "disjoint") {
        list_of_subsample_indices <- createFolds(y = response,
                                                 k = nsplits,
                                                 list = TRUE)
    } else {
        list_of_subsample_indices <- createResample(y = response,
                                                    times = nsplits,
                                                    list = TRUE)
    }

    if(is.null(cluster_object) && is.null(mc.cores)) {
        ### Estimate models serially
        list_of_models <- lapply(X = list_of_subsample_indices, parboost_fit, data,
                                 path_to_data, data_import_function, preprocessing, seed,
                                 formula, baselearner, family, control, tree_controls, cv, cores_cv = cores_cv, folds, stepsize_mstop)
    } else if (is.null(cluster_object) && !is.null(mc.cores) &&
               Sys.info()[1] != "Windows") {
        list_of_models <- mclapply(X = list_of_subsample_indices, parboost_fit, data,
                                   path_to_data, data_import_function, preprocessing, seed,
                                   formula, baselearner, family, control, tree_controls, cv, cores_cv = cores_cv, folds, stepsize_mstop, mc.cores = mc.cores)
    }
    else if (export) {
        ### Estimade models on the cluster and return the results
        clusterEvalQ(cl = cluster_object, library(parboost))
        list_of_models <- parLapply(cl = cluster_object,
                                    X = list_of_subsample_indices, parboost_fit,
                                    data, path_to_data, data_import_function,
                                    preprocessing, seed,
                                    formula, baselearner, family, control, tree_controls,
                                    cv, cores_cv = cores_cv, folds, stepsize_mstop)

    } else {
        ### Estimade models on the cluster and return the results
        clusterEvalQ(cl = cluster_object, library(parboost))
        list_of_models <- parLapply(cl = cluster_object,
                                    X = list_of_subsample_indices, parboost_fit, data = NULL,
                                    path_to_data, data_import_function,
                                    preprocessing, seed,
                                    formula, baselearner, family, control, tree_controls,
                                    cv, cores_cv = cores_cv, folds, stepsize_mstop)
    }

    ### Carry out postprocessing
    final_model<- postprocess(response, family, list_of_models,
                              postprocessing, cores_cv = cores_cv)


    ### Set up output
    ret <- list()
    ret$models <- final_model$models
    ret$predict <- final_model$postprocessed_model$predict
    ret$fitted <- final_model$postprocessed_model$fitted
    ret$call <- match.call()
    ret$family <- family
    ret$stepsize_mstop <- stepsize_mstop
    ret$split_data <- split_data
    ret$folds <- folds
    ret$postprocessing <- postprocessing
    ret$ensemble_weights <- final_model$postprocessed_model$ensemble_weights
    class(ret) <- "parboost"
    return(ret)
}

##' Internal function to fit mboost model on a subset of the data
##'
##' Fits a mboost model on each subset of the data
##' @title Fit individual parboost component using mboost
##' @param subsample_indices A numeric vector containing the indices
##' of the subsample
##' @param data A data frame containing the variables in the
##' model. It is recommended to use path_to_data instead for
##' IO efficiency. Defaults to NULL
##' @param path_to_data A string with the path to the data.
##' @param data_import_function What function should be used to import
##' the data?
##' @param preprocessing Optional preprocessing function to apply to
##' the data passed from parboost
##' @param seed Set a seed for reproducible results.
##' @param formula Formula for mboost.
##' @param baselearner Character string determining the type of base learner.
##' @param family mboost family
##' @param control mboost control
##' @param tree_controls party control
##' @param cv Cross-validate?
##' @param cores_cv Number of cores to use during cv.
##' @param folds Number of folds to use for cv.
##' @param stepsize_mstop Stepsize used for optimizing mstop.
##' @return The fitted submodel and its predictions
##' @author Ronert Obst
##' @keywords internal
parboost_fit <- function(subsample_indices, data = NULL, path_to_data, data_import_function,
                         preprocessing, seed, formula, baselearner, family,
                         control, tree_controls, cv, cores_cv = detectCores(), folds, stepsize_mstop) {


    if(is.null(data)) {
        data <- data_import_function(path_to_data)
    }

    if(!is.null(preprocessing)) {
        data <- preprocessing(data)
    }
    data_subset <- data[subsample_indices, ]

    ### Create family object for boosting
    if (family == "gaussian") boost_family <- Gaussian()
    if (family == "binomial") boost_family <- Binomial()
    if (family == "poisson") boost_family <- Poisson()

    if (!is.null(tree_controls)) {
        fitted_submodel <- blackboost(formula = formula, data = data_subset,
                                  family = boost_family, control = control, tree_controls = tree_controls)
    } else {
        if (is.character(baselearner)) {
            fitted_submodel <- mboost(formula = formula, data = data_subset,
                                      baselearner = baselearner,
                                      family = boost_family, control = control)
        } else {
            fitted_submodel <- mboost(formula = formula, data = data_subset,
                                      family = boost_family, control = control)
        }
    }
    if(!is.null(seed)) set.seed(seed)
    subsample_folds <- createFolds(y = data_subset[, as.character(formula[2])], k = folds, list = TRUE)
    subsample_folds <- t(laply(subsample_folds, function(fold, vec_length) ifelse(1:vec_length %in% fold, 0, 1), length(data_subset[, as.character(formula[2])])))

    if (cv) {
        if (stepsize_mstop == 1) {
            m_stop <- cv_subsample(fitted_submodel, folds = subsample_folds, cores_cv = cores_cv)
        } else {
            grid <- seq(0, control$mstop, stepsize_mstop)
            grid[1] <- 1
            m_stop <- cv_subsample(fitted_submodel, folds = subsample_folds, grid = grid, cores_cv = cores_cv)
        }
        fitted_submodel[m_stop]
    }
    response_index <- match(as.character(formula[[2]]), names(data))
    predictions <- predict(fitted_submodel, newdata = data[, -response_index], type = "response")

    return(list(model = fitted_submodel, predictions = predictions))
}
