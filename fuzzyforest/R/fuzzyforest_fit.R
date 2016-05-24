#' Fits fuzzy forest algorithm.
#'
#' Fits fuzzy forest algorithm.  Returns
#' fuzzy forest object.
#' @export
#' @param X                 A data.frame.
#'                          Each column corresponds to a feature vectors.
#' @param y                 Response vector.  For classification, y should be a
#'                          factor.  For regression, y should be
#'                          numeric.
#' @param Z                 A data.frame. Additional features that are not to be
#'                          screened out at the screening step.
#' @param module_membership A character vector giving the module membership of
#'                          each feature.
#' @param screen_params     Parameters for screening step of fuzzy forests.
#'                          See \code{\link[fuzzyforest]{screen_control}} for
#'                          details. \code{screen_params} is an object of type
#'                          \code{screen_control}.
#' @param select_params     Parameters for selection step of fuzzy forests.
#'                          See \code{\link[fuzzyforest]{select_control}} for details.
#'                          \code{select_params} is an object of type
#'                          \code{select_control}.
#' @param final_ntree       Number of trees grown in the final random forest.
#'                          This random forest contains all selected features.
#' @param num_processors    Number of processors used to fit random forests.
#' @param nodesize          Minimum terminal nodesize. 1 if classification.
#'                          5 if regression.  If the sample size is very large,
#'                          the trees will be grown extremely deep.
#'                          This may lead to issues with memory usage and may
#'                          lead to significant increases in the time it takes
#'                          the algorithm to run.  In this case,
#'                          it may be useful to increase \code{nodesize}.
#' @param test_features     A data.frame containing features from a test set.
#'                          The data.frame should contain the features in both
#'                          X and Z.
#' @param test_y            The responses for the test set.
#' @return An object of type \code{\link[fuzzyforest]{fuzzy_forest}}.  This
#' object is a list containing useful output of fuzzy forests.
#' In particular it contains a data.frame with list of selected features.
#' It also includes the random forest fit using the selected features.
#' @references
#' Leo Breiman (2001). Random Forests. Machine Learning, 45(1), 5-32.
#'
#' Daniel Conn, Tuck Ngun, Christina M. Ramirez (2015). Fuzzy Forests: a New
#' WGCNA Based Random Forest Algorithm for Correlated, High-Dimensional Data,
#' Journal of Statistical Software, Manuscript in progress.
#'
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17
#' @examples
#' #ff requires that the partition of the covariates be previously determined.
#' #ff is handy if the user wants to test out multiple settings of WGCNA
#' #prior to running fuzzy forests.
#' library(WGCNA)
#' library(randomForest)
#' library(fuzzyforest)
#' data(ctg)
#' y <- ctg$NSP
#' X <- ctg[, 2:22]
#'
#' #set tuning parameters for WGCNA
#' net = blockwiseModules(X, power = 6, minModuleSize = 1, nThreads = 1)
#'
#'
#' #extract module membership for each covariate
#' module_membership <- net$colors
#'
#' #set tuning parameters
#' mtry_factor <- 1; min_ntree <- 500;  drop_fraction <- .5; ntree_factor <- 1
#' screen_params <- screen_control(drop_fraction = drop_fraction,
#'                                 keep_fraction = .25, min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' select_params <- select_control(drop_fraction = drop_fraction,
#'                                 number_selected = 5,
#'                                 min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#'
#' #fit fuzzy forests
#' \donttest{
#' ff_fit <- ff(X, y, module_membership = module_membership,
#'                 screen_params = screen_params,
#'                 select_params = select_params,
#'                 final_ntree = 500)
#'
#' #extract variable importance rankings
#' vims <- ff_fit$feature_list
#'
#' #plot results
#' modplot(ff_fit)
#' }
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
ff <- function(X, y, Z=NULL, module_membership,
                        screen_params = screen_control(min_ntree=5000),
                        select_params = select_control(min_ntree=5000),
                        final_ntree = 5000,
                        num_processors=1, nodesize, test_features=NULL,
                        test_y=NULL) {
  CLASSIFICATION <- is.factor(y)
  if ( !((mode(y)=="numeric") || is.factor(y)) ) {
    stop("y must be a numeric vector or factor")
  }
  if( (!CLASSIFICATION) && (length(unique(y)) < 5) ) {
    warning("y has 5 or fewer unique values?  In this case, we recommend
            classification instead of regression.  For classification,
            y must be a factor.")
  }
  if(!is.data.frame(X)) {
    stop("X must be a data.frame.")
  }
  if(!is.null(Z)) {
    if (!is.data.frame(Z)) {
      stop("Z must be a data.frame.")
    }
  }
  if(CLASSIFICATION == TRUE) {
    if(missing(nodesize)){
      nodesize <- 1
    }
  }
  if(CLASSIFICATION == FALSE) {
    if(missing(nodesize)){
      nodesize <- 5
    }
  }
  screen_control <- screen_params
  select_control <-  select_params
  module_list <- unique(module_membership)
  if(num_processors > 1) {
    #set up parallel backend
    cl = parallel::makeCluster(num_processors)
    parallel::clusterCall(cl, library, package = "randomForest", character.only = TRUE)
    doParallel::registerDoParallel(cl)
    #close parallel backend on exit
    on.exit(try(parallel::stopCluster(cl), silent=TRUE))
  }
  survivors <- vector('list', length(module_list))
  drop_fraction <- screen_control$drop_fraction
  mtry_factor <- screen_control$mtry_factor
  ntree_factor <- screen_control$ntree_factor
  min_ntree <- screen_control$min_ntree
  keep_fraction <- screen_control$keep_fraction
  if(ncol(X)*keep_fraction < select_control$number_selected){
    warning(c("ncol(X)*keep_fraction < number_selected", "\n",
              "number_selected will be set to floor(ncol(X)*keep_fraction)"))
              select_control$number_selected <- max(floor(ncol(X)*keep_fraction), 1)
  }

  for (i in 1:length(module_list)) {
    module <- X[, which(module_membership == module_list[i]), drop=FALSE]
    num_features <- ncol(module)
    #TUNING PARAMETER mtry_factor
    if(CLASSIFICATION == TRUE) {
      mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
      if(missing(nodesize)){
        nodesize <- 1
      }
    }
    if(CLASSIFICATION == FALSE) {
      mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
      if(missing(nodesize)){
        nodesize <- 5
      }
    }
    #TUNING PARAMETER ntree_factor
    ntree <- max(num_features*ntree_factor, min_ntree)
    #TUNING PARAMETER keep_fraction
    target = ceiling(num_features * keep_fraction)
    while (num_features >= target){
      if(num_processors > 1) {
        rf = foreach(ntree = rep(ntree/num_processors, num_processors),
                   .combine = combine, .packages = 'randomForest') %dorng% {
                   randomForest(module, y, ntree = ntree, mtry = mtry,
                   importance = TRUE, scale = FALSE, nodesize=nodesize) }
      }
      if(num_processors == 1) {
        rf <- randomForest(module, y, ntree = ntree, mtry = mtry,
                           importance = TRUE, scale = FALSE,
                           nodesize = nodesize)
      }
      var_importance <- importance(rf, type=1, scale=FALSE)
      var_importance <- var_importance[order(var_importance[, 1],
                                             decreasing=TRUE), ,drop=FALSE]
      reduction <- ceiling(num_features*drop_fraction)
      if(num_features - reduction > target) {
          trimmed_varlist <- var_importance[1:(num_features - reduction), ,drop=FALSE]
          features <- row.names(trimmed_varlist)
          module <- module[, which(names(module) %in% features)]
          num_features <- length(features)
          if(CLASSIFICATION == TRUE) {
            mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
          }
          if(CLASSIFICATION == FALSE) {
            mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
          }
          ntree <- max(num_features*ntree_factor, min_ntree)
        }
      else {
          num_features <- target - 1
          mod_varlist <- var_importance[, 1][1:target]
          features <- row.names(var_importance)[1:target]
          survivors[[i]] <- cbind(features, mod_varlist)
          row.names(survivors[[i]]) <- NULL
          survivors[[i]] <- as.data.frame(survivors[[i]])
          survivors[[i]][, 1] <- as.character(survivors[[i]][, 1])
          survivors[[i]][, 2] <- as.numeric(as.character(survivors[[i]][, 2]))
        }
    }
  }
  survivor_list <- survivors
  names(survivor_list) <- module_list
  survivors <- do.call('rbind', survivors)
  survivors <- as.data.frame(survivors, stringsAsFactors = FALSE)
  survivors[, 2] <- as.numeric(survivors[, 2])
  names(survivors) <- c("featureID", "Permutation VIM")
  X_surv <- X[, names(X) %in% survivors[,1]]
  if(!is.null(Z)) {
    X_surv <- cbind(X_surv, Z, stringsAsFactors=FALSE)
  }
  select_args <- list(X_surv, y, num_processors, nodesize)
  select_args <- c(select_args, select_control)
  names(select_args)[1:4] <- c("X", "y", "num_processors", "nodesize")
  select_results <- do.call("select_RF", select_args)
  final_list <- select_results[[1]]
  selection_list <- select_results[[2]]
  final_list[, 2] <- round(as.numeric(final_list[, 2]), 4)
  row.names(final_list) <- NULL
  colnames(final_list) <- c("feature_name", "variable_importance")
  final_list <- as.data.frame(final_list, stringsAsFactors=FALSE)
  final_list[, 2] <- as.numeric(final_list[, 2])
  final_list <- cbind(final_list, rep(".", dim(final_list)[1]),
                     stringsAsFactors=FALSE)
  names(final_list)[3] <- c("module_membership")
  select_X <- names(X)[which(names(X) %in% final_list[, 1])]
  select_mods <- module_membership[which(names(X) %in% final_list[,1])]
  select_order <- final_list[, 1][which(final_list[,1] %in% names(X))]
  select_mods <- select_mods[match(select_order, select_X)]
  final_list[, 3][final_list[, 1] %in% names(X)] <- select_mods
  final_X <- X[, names(X) %in% final_list[, 1], drop=FALSE]
  if(!is.null(Z)) {
    final_X <- cbind(final_X, Z[, names(Z) %in% final_list[, 1], drop=FALSE],
                     stringsAsFactors=FALSE)
  }
  current_p <- dim(final_X)[2]
  if(CLASSIFICATION == TRUE) {
    final_mtry <- min(ceiling(select_control$mtry_factor*sqrt(current_p)),
                      current_p)
  }
  if(CLASSIFICATION == FALSE) {
    final_mtry <- min(ceiling(select_control$mtry_factor*current_p/3),
                      current_p)
  }
  if(!is.null(test_features)) {
    test_features <- test_features[, which(names(test_features) %in%
                                      names(final_X))]
  }
  final_rf <- randomForest(x=final_X, y=y, mtry=final_mtry, ntree=final_ntree,
                           importance=TRUE, nodesize=nodesize,
                           xtest=test_features, ytest=test_y)
  module_membership <- as.data.frame(cbind(names(X), module_membership),
                                     stringsAsFactors=FALSE)
  names(module_membership) <- c("feature_name", "module")
  out <- fuzzy_forest(final_list, final_rf, module_membership,
                      survivor_list=survivor_list,
                      selection_list=selection_list)
  return(out)
}


#' Fits WGCNA based fuzzy forest algorithm.
#'
#' Fits fuzzy forest algorithm using WGCNA.  Returns
#' fuzzy forest object.
#' @export
#' @param X                 A data.frame. Each column corresponds to a feature
#'                          vector.  WGCNA will be used to cluster the
#'                          features in X.  As a result, the features should be
#'                          all be numeric.  Non-numeric features may be input
#'                          via Z.
#' @param y                 Response vector.  For classification, y should be a
#'                          factor.  For regression, y should be
#'                          numeric.
#' @param Z                 Additional features that are not to be screened out
#'                          at the screening step.  WGCNA is not carried out on
#'                          features in Z.
#' @param WGCNA_params      Parameters for WGCNA.
#'                          See \code{\link[WGCNA]{blockwiseModules}} and
#'                          \code{\link[fuzzyforest]{WGCNA_control}} for details.
#'                          \code{WGCNA_params} is an object of type
#'                          \code{WGCNA_control}.
#' @param screen_params     Parameters for screening step of fuzzy forests.
#'                          See \code{\link[fuzzyforest]{screen_control}} for details.
#'                          \code{screen_params} is an object of type
#'                          \code{screen_control}.
#' @param select_params     Parameters for selection step of fuzzy forests.
#'                          See \code{\link[fuzzyforest]{select_control}} for details.
#'                          \code{select_params} is an object of type
#'                          \code{select_control}.
#' @param final_ntree       Number of trees grown in the final random forest.
#'                          This random forest contains all selected features.
#' @param num_processors    Number of processors used to fit random forests.
#' @param nodesize          Minimum terminal nodesize. 1 if classification.
#'                          5 if regression.  If the sample size is very large,
#'                          the trees will be grown extremely deep.
#'                          This may lead to issues with memory usage and may
#'                          lead to significant increases in the time it takes
#'                          the algorithm to run.  In this case,
#'                          it may be useful to increase \code{nodesize}.
#' @param test_features     A data.frame containing features from a test set.
#'                          The data.frame should contain the features in both
#'                          X and Z.
#' @param test_y            The responses for the test set.
#' @return An object of type \code{\link[fuzzyforest]{fuzzy_forest}}.  This
#' object is a list containing useful output of fuzzy forests.
#' In particular it contains a data.frame with list of selected features.
#' It also includes the random forest fit using the selected features.
#' @references
#' Leo Breiman (2001). Random Forests. Machine Learning, 45(1), 5-32.
#'
#' Daniel Conn, Tuck Ngun, Christina M. Ramirez (2015). Fuzzy Forests: a New
#' WGCNA Based Random Forest Algorithm for Correlated, High-Dimensional Data,
#' Journal of Statistical Software, Manuscript in progress.
#'
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17
#' @examples
#' library(WGCNA)
#' library(randomForest)
#' library(fuzzyforest)
#' data(ctg)
#' y <- ctg$NSP
#' X <- ctg[, 2:22]
#' WGCNA_params <- WGCNA_control(p = 6, minModuleSize = 1, nThreads = 1)
#' mtry_factor <- 1; min_ntree <- 500;  drop_fraction <- .5; ntree_factor <- 1
#' screen_params <- screen_control(drop_fraction = drop_fraction,
#'                                 keep_fraction = .25, min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' select_params <- select_control(drop_fraction = drop_fraction,
#'                                 number_selected = 5,
#'                                 min_ntree = min_ntree,
#'                                 ntree_factor = ntree_factor,
#'                                 mtry_factor = mtry_factor)
#' \donttest{
#' wff_fit <- wff(X, y, WGCNA_params = WGCNA_params,
#'                 screen_params = screen_params,
#'                 select_params = select_params,
#'                 final_ntree = 500)
#'
#' #extract variable importance rankings
#' vims <- wff_fit$feature_list
#'
#' #plot results
#' modplot(wff_fit)
#' }
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
wff <- function(X, y, Z=NULL, WGCNA_params=WGCNA_control(power=6),
                        screen_params=screen_control(min_ntree=5000),
                        select_params=select_control(min_ntree=5000),
                        final_ntree=500, num_processors=1, nodesize,
                        test_features=NULL, test_y=NULL) {
  if ( !("package:WGCNA" %in% search()) ) {
    stop("WGCNA must be loaded and attached. Type library(WGCNA) to do so.")
  }
  numeric_test <- sapply(X, is.numeric)
  if (sum(numeric_test) != dim(X)[2]) {
    stop("To carry out WGCNA, all columns of X must be numeric.")
  }
  CLASSIFICATION <- is.factor(y)
  if(CLASSIFICATION == TRUE) {
    if(missing(nodesize)){
      nodesize <- 1
    }
  }
  if(CLASSIFICATION == FALSE) {
    if(missing(nodesize)){
      nodesize <- 5
    }
  }
  WGCNA_control <- WGCNA_params
  screen_control <- screen_params
  select_control <-  select_params
  WGCNA_args <- list(X,WGCNA_control$power)
  WGCNA_args <- c(WGCNA_args, WGCNA_control$extra_args)
  names(WGCNA_args) <- c("datExpr", "power", names(WGCNA_control$extra_args))
  bwise <- do.call("blockwiseModules", WGCNA_args)
  module_membership <- bwise$colors
  screen_drop_fraction <- screen_control$drop_fraction
  screen_keep_fraction <- screen_control$keep_fraction
  screen_mtry_factor <- screen_control$mtry_factor
  screen_ntree_factor <- screen_control$ntree_factor
  screen_min_ntree <- screen_control$min_ntree
  out <- ff(X, y, Z, module_membership,
                    screen_control, select_control, final_ntree,
                    num_processors, nodesize=nodesize,
                    test_features=test_features, test_y=test_y)
  out$WGCNA_object <- bwise
  return(out)
}






