#Copyright (c) 2016 Jussi Korpela (Finnish Institute of Occupational Healt, jussi.korpela@ttl.fi) and Andreas Henelius (Finnish Institute of Occupational Healt, andreas.henelius@iki.fi)
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

## -----------------------------------------------------------------------------
## The CoCoReg algorithm
## -----------------------------------------------------------------------------

#' Create a path
#' 
#' @description 
#' Given a set of numbers and a starting point, create a path from
#' the starting point through each point. If cyclic = FALSE, the path
#' stops at p[length(p)] and if cyclic = TRUE the patch goes back to ind.
#'
#' @param ind [1,1] int, The starting point.
#' @param p [1,m] int, A set of numbers through which the path must go,
#'        excluding the starting point.
#' @param cyclic boolean, If TRUE, the path leads back to the starting point.
#'        If FALSE, the path stops at last element of p.
#'                
#' @return A list containing the path definition.
#'
#' @keywords internal
pathify <- function(ind, p, cyclic = FALSE) {
    if (cyclic)
        p  <- c(ind, rep(p, each = 2), ind)
    else
        p  <- c(p[1], rep(p[-1], each = 2), ind)

    p2 <- split(p, ceiling(seq_along(p) / 2))
    path_lst <- lapply(p2, function(i) paste(as.character(i), collapse = "-"))
    attr(path_lst, "startDS") <- p[1]
    path_lst
}


#' Generate all/some paths between points
#'
#' @param ind [1,2] int, The starting and ending point c(start, end).
#' @param n [1,1] int, Number of points in the whole set.
#' @param cyclic boolean, Should the path be cyclic (1-2-1) or noncyclic (2-1).
#' @param sample_paths boolean, If FALSE, all possible paths are generated.
#'        If true one path per ending point is selected.
#'
#' @return A list of lists containing the paths.
#'
#' @export
generate_paths <- function(ind, n, cyclic = FALSE, sample_paths = FALSE) {
    if (cyclic){
        paths <- generate_paths_cyclic(ind, n)
    } else {
        paths <- generate_paths_noncyclic(ind, n, sample_paths = sample_paths)
    }
}


#' Generate cyclic paths
#'
#'@description
#' From a set of n numbers, generate all possible paths starting from
#' and ending on a given number.
#'
#' @param ind [1,1] ind, The starting dataset (equals to ending dataset because of cycle).
#' @param n [1,1] ind, The number of datasets.
#'
#' @return A list of lists containing the paths.
#'
#' @importFrom combinat permn
#'
#' @export
generate_paths_cyclic <- function(ind, n) {
    if (n == 2) {
        if (ind == 1)
            res <- list(list("1-2", "2-1"))
            attr(res[[1]], "startDS") <- 1
        if (ind == 2)
            res <- list(list("2-1", "1-2"))
            attr(res[[1]], "startDS") <- 2
    }

    if (n > 2) {
        plist <- combinat::permn(setdiff(seq.int(n), ind))
        res   <- lapply(plist, function(i) pathify(ind, i, cyclic = TRUE))
    }

    res
}


#' Generate non-cyclic paths
#'
#' @description 
#' From a set of n numbers, generate all possible paths starting from
#' and ending on a given number.
#'
#' @param ind The starting dataset
#' @param n The number of datasets.
#' @param sample_paths boolean, If FALSE, all possible paths are generated.
#'        If TRUE one path per ending point is selected.
#'
#' @return A list of lists containing the paths.
#'
#' @importFrom combinat permn
#'
#' @export
generate_paths_noncyclic <- function(ind, n, sample_paths = FALSE) {
    if (n == 2) {
        if (ind == 1){
            res <- list(list("2-1"))
            attr(res[[1]], "startDS") <- 2
        }
        if (ind == 2){
            res <- list(list("1-2"))
            attr(res[[1]], "startDS") <- 1
        }
    }

    if (n > 2) {
      if (sample_paths & (n>3)){
        ## pick one path for each possible starting index
        plist <- lapply(setdiff(seq.int(n),ind), function(i) c(i, sample(setdiff(seq.int(n), c(i, ind)))))
        ## Note: sample(3) does not sample from vector c(3) but from 1:3!
        ## This is why sample_paths = TRUE has effect only for data collections
        ## with K > 3.
      } else {
        ## all paths for all possible starting indices
        plist <- combinat::permn(setdiff(seq.int(n), ind))  
      }
      res <- lapply(plist, function(i) pathify(ind, i, cyclic = FALSE))
    }
    res
}


#' Generate a mapping function between two datasets
#' 
#' @description 
#' Generate a mapping function between two datasets,
#' using some method, such as linear regression (lm),
#' or some classifier such as a random forest (randomForest).
#' 
#' Wraps the function as well as data into a single object.
#'
#' @param method A function to be used in the mapping. A function object.
#'
#' @return Returns a function object that does the mapping between two datasets i.e.
#' from dataset 1 to dataset 2.
#'
#' @export
generate_mapping_function <- function(method = lm) {
    function(data1, data2) {
        f_list <- lapply(seq.int(ncol(data2)),
                         function(i) method(data2[,i] ~ 0 + . , data = data1))
        function(x) {
            tmp           <- lapply(f_list, function(i) predict(i, newdata = x))
            tmp           <- data.frame(do.call(cbind, tmp))
            colnames(tmp) <- colnames(data2)
            tmp
        }
    }
}


## Create mapping functions:

#' Mapping stats::lm
#' 
#' @param data1 data.frame, Dataset 1, the independent variables
#' @param data2 data.frame, Dataset 2, the dependent variables
#' 
#' @export
mapping_lm  <- generate_mapping_function(method = stats::lm)

#' Mapping MASS::rlm
#' 
#' @param data1 data.frame, Dataset 1, the independent variables
#' @param data2 data.frame, Dataset 2, the dependent variables
#' 
#' @importFrom MASS rlm
#' @export
mapping_rlm  <- generate_mapping_function(method = MASS::rlm)

#' Mapping randomForest
#' 
#' @param data1 data.frame, Dataset 1, the independent variables
#' @param data2 data.frame, Dataset 2, the dependent variables
#' 
#' @export
mapping_rf  <- generate_mapping_function(method = randomForest::randomForest)

#' Mapping svm
#' 
#' @param data1 data.frame, Dataset 1, the independent variables
#' @param data2 data.frame, Dataset 2, the dependent variables
#' 
#' @importFrom e1071 svm
#' @export
mapping_svm <- generate_mapping_function(method = e1071::svm)

#' SVM using sigmoid kernel
#' 
#' @param ... Further arguments passed on to e1071::svm()
#' 
#' @importFrom e1071 svm
#' @export
svm_sigmoid <- function(...) {
  e1071::svm(..., kernel = "sigmoid", coef0 = 0)
}


#' Mapping svm using sigmoid
#' 
#' @param data1 data.frame, Dataset 1, the independent variables
#' @param data2 data.frame, Dataset 2, the dependent variables
#' 
#' @export
mapping_svm_sigmoid <- generate_mapping_function(method = svm_sigmoid)

#' Define a mapping function using MASS::lm.ridge
#' 
#' @param data1 data.frame, Dataset 1, the independent variables
#' @param data2 data.frame, Dataset 2, the dependent variables
#' 
#' @return Returns a function object that does the mapping between two datasets.
#' 
#' @importFrom MASS lm.ridge
#' @export
mapping_lmridge <- function(data1, data2) {
  f_list <- lapply(seq.int(ncol(data2)), function(i) MASS::lm.ridge(data2[,i] ~ . , data = data1))
  
  function(x) {
    tmp           <- lapply(f_list, function(model){
      matrix( scale(as.matrix(x), center = F, scale = model$scales) %*% matrix(model$coef) + model$ym
        , nrow=nrow(data2))} )
    tmp           <- data.frame(do.call(cbind, tmp))
    colnames(tmp) <- colnames(data2)
    tmp
  }
}



#' Define a mapping function using pls::pcr
#' 
#' @param data1 data.frame, Dataset 1, the independent variables
#' @param data2 data.frame, Dataset 2, the dependent variables
#' 
#' @return Returns a function object that does the mapping between two datasets.
#' 
#' @importFrom pls pcr
#' @export
mapping_pcr <- function(data1, data2) {
  f_list <- lapply(seq.int(ncol(data2)), function(i) pls::pcr(data2[,i] ~ . , data = data1))
  
  function(x) {
    tmp           <- lapply(f_list, function(model){
                matrix(predict(model, x)[,1,2], nrow=nrow(data2))} )
    tmp           <- data.frame(do.call(cbind, tmp))
    colnames(tmp) <- colnames(data2)
    tmp
  }
}


#' Define a mapping function using glmnet::glmnet
#'
#' @param data1 data.frame, Dataset 1, the independent variables
#' @param data2 data.frame, Dataset 2, the dependent variables
#' 
#' @return Returns a function object that does the mapping between two datasets.
#' 
#' @importFrom glmnet glmnet predict.glmnet
#' @export
mapping_glmnet <- function(data1, data2) {
  f_list <- lapply(seq.int(ncol(data2)), function(i) glmnet::glmnet(as.matrix(data1), data2[,i]) )
  
  function(x) {
    tmp           <- lapply(f_list, function(i) predict.glmnet(i, newx = as.matrix(x), s=49))
    tmp           <- data.frame(do.call(cbind, tmp))
    colnames(tmp) <- colnames(data2)
    tmp
  }
}


#' Generate all possible pairwise mappings between the given multivariate
#' datasets.
#' 
#' @description 
#' The following naming convention is used in the output:
#' '1-2' means '1' mapped to '2', i.e.,, '2' explained by '1'.
#'
#' @param data A list of data.frames (the datasets)
#' @param mapping_function (optional) Default is \code{mapping_lm}.
#'
#' @return A named list containing the pairwise mapping functions.
#'
#' @export
create_mappings <- function(data, mapping_function = mapping_lm) {
    n   <- length(data)
    out <- vector(mode = "list", length = ((n * n) - n))
    ind <- 1

    for (i in seq.int(n)) {
        for (j in seq.int(n)) {
            if (i != j) {
                out[[ind]]      <-  mapping_function(data[[i]], data[[j]])
                names(out)[ind] <- paste(i, j, sep = "-")
                ind             <- ind + 1
            }
        }
    }

    out
}


#' Calculate the composition formed by applying all functions
#' in the given path to a dataset.
#'
#' @param x A data frame or vector
#' @param path A list describing the path.
#' @param mappings A list containing the mapping functions described in the path.
#'                  Usually a list of all M*M-M available mappings between the M
#'                  data sets.
#'
#' @return A vector containing the result of the composition.
#'
#' @export
compose <- function(x, path, mappings) {
    for (i in path){
      #str(x)
      x <- mappings[[i]](x)
    }
        
    x
}


#' Helper function to get the starting dataset based on
#' a path.
#'
#' @param p A path.
#'
#' @return A number indicating the starting dataset
get_starting_dataset <- function(p) {
    as.numeric(unlist(strsplit(p[[1]], "-"))[1])
}


#' Calculate the average of the composition formed by applying all functions
#' in all possible paths to a dataset.
#'
#' @param x A list of data frames.
#' @param paths A list of list with paths.
#' @param mappings A list containing the mapping functions described in the path.
#'
#' @return A vector containing the average composition from all the paths.
#'
#' @importFrom abind abind
#'
#' @export
compose_all <- function(x, paths, mappings) {
    tmp <- vector(mode = "list", length = length(paths))
    for (i in seq.int(length(paths))){
        #cat(sprintf("\n\ncompose_all, round:%d",i))
        tmp[[i]] <- compose(x[[get_starting_dataset(paths[[i]])]], paths[[i]], mappings)
    }

    tmp <- lapply(tmp, as.matrix)
    tmp <- do.call(abind::abind, c(tmp, along = 3)) #catenated along the third dimension
    as.data.frame(apply(tmp, c(1, 2), mean)) # average out 3rd dimension
}


#' The Common Components by Regression (CoCoReg) algorith.
#' 
#' @description 
#' An algorithm that extracts common variation between datasets using regression.
#'
#' @param data   [1,K] list of data.frames.
#' @param cyclic boolean, Operation mode: cyclic or non-cyclic
#' @param mapping_function function, The function to use in mappings.
#'                          See mapping_lm() for an example.
#' @param sample_paths boolean, If FALSE all paths are computed. If TRUE a
#'        subset of paths is taken: one (random) path for each starting point.
#'        Currently implemented only for cyclic=F.
#' @param center_data boolean, Should the data be centered?
#' @param scale_data boolean, Should the data be scaled?
#'
#'
#' @return A list with elements:
#' \item{$data:}{[1,K] list of data.frames containing the joint information,
#'                organised identically to the input data.}
#' \item{$mappings:}{[1,K*K-K] list of functions, mappings between datasets}
#' \item{$paths:}{[(K-1)(K-2)!, K] list of lists, paths for each data set}
#' \item{$cyclic:}{input cyclic as is}
#' \item{$sample_paths:}{boolean, TRUE if paths have been sampled, FALSE otherwise.}
#' \item{$dataid:}{string, Dataset identifier string}
#' \item{$method:}{string, Analysis method identifier string}
#' \item{$wall_time_taken:}{[1,1] double, Time taken to run the analysis in seconds}
#'
#' @examples
#' dc <- create_syn_data_toy()
#' ccr <- cocoreg(dc$data)
#' 
#' ggplot_dflst(dc$data, ncol=1)
#' ggplot_dflst(ccr$data, ncol=1)
#' 
#' \dontrun{
#' ggplot_dclst(list(orig = dc$data, ccr = ccr$data)) 
#' ggplot_dclst(list(orig = dc$data, shared = ccr$data), legendMode = 'none')
#' ggplot_dclst(list(orig = dc$data, ccr = ccr$data), legendMode = 'all')
#' }
#'
#' @export
cocoreg <- function(data, cyclic = FALSE, mapping_function = mapping_lm,
                    sample_paths = FALSE, center_data = T, scale_data = F) {
  start_time <- Sys.time()
  
  ## Validate data and convert to internal format
  if ("id" %in% names(attributes(data)) ){
    dataid <- attr(data, "id")
  } else {
    dataid <- "NA"
  }
  data  <- validate_data(data)
  dmeta <- get_dc_meta(data)
  data  <- rename_variables(data)
  if (center_data | scale_data){
    data  <- lapply(data, function(el){df_scale(el, center=center_data, scale=scale_data)})  
  }
  
  ## Compute cocoreg
  mappings <- create_mappings(data, mapping_function = mapping_function)
  ccr_data <- lapply(seq.int(length(data)), function(i) compose_all(data, generate_paths(i, n = length(data), cyclic = cyclic, sample_paths=sample_paths), mappings))
  paths <- sapply(seq.int(length(data)), function(i)  generate_paths(i, n = length(data), cyclic = cyclic, sample_paths=sample_paths))
  #TODO: test if the use of vapply() would not leave out the dim attribute?
  if (length(data)==2){
    dim(paths) <- c(1, 2) #prevent from becoming dimensionless  
  }
  
  ## Convert back to original data format
  ccr_data <- apply_dc_meta(ccr_data, dmeta) 
  
  ## Return
  res_lst <- list(
    "data"         = ccr_data,
    "mappings"     = mappings,
    "paths"        = paths,
    "cyclic"       = cyclic,
    "sample_paths" = sample_paths,
    "dataid"       = dataid,
    "method"       = 'ccr')
  
  res_lst        <- c(res_lst, list("wall_time_taken" = as.numeric(difftime(Sys.time(), start_time, units = "sec"))))
  class(res_lst) <- "cocoreg" ## Needed for effective recursive traversing of results in a nested list
  res_lst
}




## -----------------------------------------------------------------------------
## Helper functions
## -----------------------------------------------------------------------------


#' Computes the R^2 (variance explained) between two lists of data.frames
#' 
#' @param df_orig_lst List of original data.frames
#' @param df_est_lst List of estimated data.frames
#' 
#' @return Returns a data.frame with R2 values, one value for each data set and variable. Molten/long format.
#' 
#' @importFrom reshape melt melt.data.frame
#' 
#' @export

## 
average_R2_dflst <- function(df_orig_lst, df_est_lst){
  df_var_explained <- function(df1, df2){
    as.data.frame(mapply(function(x,y){var_explained(x,y)$R2},
                         df1, df2, SIMPLIFY = F))
  }
  tmp <- mapply(df_var_explained, df_orig_lst, df_est_lst, SIMPLIFY = F)
  tmp <- reshape::melt(tmp)
  tmp$dsvar <- tmp$variable
  tmp$variable <- "R2"
  tmp$ds <- as.factor(tmp$L1)
  tmp$L1 <- NULL
  tmp
}


#' Validate a data collection for use with cocoreg
#' 
#' @description 
#' Check that the data collection has all the required properties.
#' 
#' @param df_list list of data.frames, The data collection to validate
#' 
#' @return A list of data.frames that conform to the requirements
#' 
#' @export
validate_data <- function(df_list){  
  
  ## Add dataset names if missing
  if (is.null(names(df_list))){
    warning('validate_data: datasets do not have names, creating names...')
    names(df_list) <- sprintf("DS%d", seq.int(length(df_list)))
  }
  
  ## Make sure the datasets are data.frames
  df_match <- vapply(df_list, is.data.frame, logical(1))
  if (!all(df_match)){
    warning('validate_data: some datasets are not data.frames, converting...')
    df_list <- lapply(df_list, as.data.frame)
  }
  
  #TODO: Maybe check that variable names follow the COCOREG convention and rename if necessary
  ## Rename variables
  #if (rename_vars){
  #  df_list <- rename_variables(df_list)
  #}
  
  df_list
}


#' Rename variables of a data collection
#' 
#' @description 
#' Rename variable in all datasets such that the data.frame list conforms to the 
#' requirements of CoCoReg.
#'
#' @param df_list  list of data.frames, The datasets to process
#'
#' @return A list of data.frames with changed variable names. Original dimension names are stored
#'          as an attribute.
#'
#' @export
rename_variables <- function(df_list){
  for (i in 1:length(df_list)){
    attr(df_list[[i]],'orig_dimnames') <- dimnames(df_list[[i]])
    ## CoCoReg dimension names
    tmp_colnames <- paste(sprintf("x%d",i), as.character(seq.int(ncol(df_list[[i]]))), sep="_")
    dimnames(df_list[[i]]) <- list((1:nrow(df_list[[i]])), tmp_colnames)
  }
  df_list
}


#' Extract important properties of data collection
#'
#' @param df_list  list of data.frames, The data collection to process
#' @param type string, If 'current' then data collection metadata is collected from the
#'        data collection itself. If 'original' then metadata is collected from special attributes.
#'
#' @return A list with elements:
#' \item{$dcnames:}{Dataset names}
#' \item{$dcdimnames:}{A list of dataset dimension names for each dataset}
#'
#' @export
get_dc_meta <- function(df_list, type = 'current'){
  
  if (type == 'current'){
    meta <- list( dcnames = names(df_list),
                  dcdimnames = lapply(df_list, dimnames) )
  } else if (type == 'original'){
    meta <- list( dcnames = names(df_list),
                  dcdimnames = lapply(df_list, function(el){attr(el, 'orig_dimnames')}) )
  } else {stop()}
  meta
}


#' Apply extracted properties of a data collection to a data collection (restore)
#' 
#' @param df_list  list of data.frames, The data collection to process
#' @param meta, list, Output of get_dc_meta()
#' 
#' @return A list of data.frames, the data collection with updated metadata
#' 
#' @export
apply_dc_meta <- function(df_list, meta){
  names(df_list) <- meta$dcnames
  for (ind in seq.int(length(df_list))){
    dimnames(df_list[[ind]]) <- meta$dcdimnames[[ind]]
  }
  df_list
}


#' Extract R2 values from a list of mappings using summary()
#'
#' @param mappings [M*M-M] list, Exhaustive list of mappings between the M datasets
#' @param n_datasets integer, Number of datasets i.e. M
#' @param aggfun function, A function to apply when aggregating R2 values over variables in
#'                a multivariate dataset
#'
#' @return A [M,M] matrix of R2 values stored such that the R2 value for mapping
#'          a->b is read from row a and column b.
#'
#' @examples
#' \dontrun{
#' ccr <- cocoreg(data_collection)
#' R2mat <- mappings_R2_matrix(ccr$mappings, length(ccr$data))
#' }
#'
#' @export
mappings_R2_matrix <- function(mappings, n_datasets, aggfun = mean){

    ## Compute values
    mappings_summaries <- lapply(mappings, function(el){summary(el)})

    ## Extract and average R2 values
    if ("r.squared" %in% names(summary(mappings[[1]]))) {
        ## univariate datasets (?)
        mappings_R2 <- lapply(mappings_summaries, function(el){el$r.squared})
    } else {
        ## multivariate datasets
        sbf_r2_fun <- function(summ) {
            aggfun(sapply(summ, function(el){el$r.squared}, simplify=T))
        }
        mappings_R2 <- lapply(mappings_summaries, sbf_r2_fun)
    }

    ## Create respective R2 matrix
    R2mat <- matrix(NA, nrow=n_datasets, ncol=n_datasets)
    for (i in 1:length(mappings_R2)) {
        i_inds <- as.integer(strsplit(names(mappings_R2)[[i]],"-")[[1]])
        R2mat[i_inds[1], i_inds[2]] <- mappings_R2[[i]]
    }

    R2mat
}


#' Compute D_joint for dataset i separately for all paths
#' Can be used to study path variability
#'
#' @param data_orig list of data.frames, Original data collection
#' @param ccr list, output of cocoreg(data_orig)
#' @param ds_ind integer, index of the dataset to process
#'
#' @return A list of data.frames, D_joint corresponding to each path that ends
#'         at 'ds_ind' i.e. paths defined by ccr$paths[,ds_ind].
#'         Dimensions of the matrices are the same as for the data.frames
#'         in data_orig.
#'
#' @examples
#' \dontrun{
#' ccr <- cocoreg(data_list)
#' jibp <- cocoreg_by_path(ccr, 1)
#' }
#' @export
cocoreg_by_path <- function(data_orig, ccr, ds_ind){
    data_orig <- validate_data(data_orig)
    data_orig <- rename_variables(data_orig)
    
    ##paths <- generate_paths(ds_ind, n=length(ccr$data), cyclic=ccr$cyclic)
    paths <- ccr$paths[,ds_ind]

    ## Compute path joint informations
    tmp <- vector(mode = "list", length = length(paths))
    for (i in seq.int(length(paths))){
        tmp[[i]] <- compose(data_orig[[ attr(paths[[i]],'startDS') ]],
                            paths[[i]], ccr$mappings)
    }
    names(tmp) <- vapply(paths, function(x){attr(x,'startDS')}, numeric(1))
    tmp
}


#' Compute "variance" of the vectors var()
#' 
#' @param vec_lst Data to process as a list of numeric vectors
#' @param mean_vec (optional) Desired mean vector as a numeric vector
#' 
#' @return Variance of data values around mean as a numeric matrix
#' 
#' @importFrom abind abind
#' @export
vector_variability <- function(vec_lst,
                               mean_vec = apply(abind::abind(vec_lst, along = 2), 1, mean) ){
    ## Variability around mean
    diff_mat <- abind::abind(lapply(vec_lst, function(x){x-mean_vec}), along=2)
    apply(diff_mat, 2, var)
}


#' Compute "variance" of the matrices using Frobenius norm.
#' Variance is by default computed with respect to the mean of the matrices.
#'
#' @param mat_lst A [1,M] list of [I,J] matrices from which variability should
#'        be computed
#' @param mean_mat A [I,J] matrix describing the mean observation for mat_lst.
#'
#' @return A list with elements
#' \item{fbnorm}{Frobenius norm values for each of the matrices}
#' \item{diff}{[I,J,M] matrix of differences mat_lst-mean_mat}
#'
#' @importFrom abind abind
#' @export
matrix_variability <- function(mat_lst,
                               mean_mat = apply(abind::abind(mat_lst, along = 3), c(1,2), mean) ){
    ## Variability around mean
    diff_mat <- abind::abind(lapply(mat_lst, function(x){x-mean_mat}), along=3)
    list(fbnorm=apply(diff_mat, 3, norm, type="F"),
         diff=diff_mat)
}


#' Determine the variability of matrices under row suffling
#'
#' @param mat_lst A list of matrices
#' @param B integer, Number of times to sample (shuffle)
#'
#' @return [B,K] matrix, Frobenius norm vectors corresponding to row suffling
#'
#' @export
row_suffle_variability <- function(mat_lst, B=50){
    nd_var = matrix(NA, nrow=B, ncol=length(mat_lst))
    for (i in 1:B){
        ## suffle
        mat_lst_suff <- lapply(mat_lst, function(x){apply(x,2,sample,size=nrow(x))})
        ## compute variability
        nd_var[i,] <- matrix_variability(mat_lst_suff)$fbnorm
    }
    nd_var
}


#' Compute ds_variability for all datasets in a data collection
#' 
#' @param data Unprocessed dataset as a list of data.frames
#' @param ccr Processed dataset as a list of data.frames, output of cocoreg()$data
#' 
#' @return Path variability as a data.frame
#' 
#' @export
dc_variability <- function(data, ccr){
  data <- validate_data(data)
  data <- rename_variables(data)
  x <- lapply(seq.int(length(data)),
              function(i){ds_variability(data, ccr, i)})
  names(x) <- names(data)
  x_df <- dflst2df(x, id_var_name="ds")
}


#' Compute variability_components for a dataset
#'
#' Note: might not work for cyclic ccr <TODO>
#' 
#' @param ds_ind Starting dataset of the set of paths to analyze as [1,1] integer
#' @inheritParams dc_variability
#' 
#' @return Path variability as a data.frame
#' 
#' @importFrom reshape melt melt.data.frame
ds_variability <- function(data, ccr, ds_ind){
    
    paths <- ccr$paths[,ds_ind]

    ## Compute path joint informations
    tmp <- vector(mode = "list", length = length(paths))
    for (i in seq.int(length(paths))){
        tmp[[i]] <- compose(data[[get_starting_dataset(paths[[i]])]],
                            paths[[i]], ccr$mappings)
    }
    names(tmp) <- vapply(paths, function(x){attr(x,'startDS')}, numeric(1))


    path_arr <- dflst2array(tmp, dim=3)
    myfun <- function(vec){
        variability_components(vec, grp=dimnames(path_arr)[[3]], fun=ss)
    }
    ss_arr <- apply(path_arr, c(1,2), myfun)
    ss_arr <- apply(ss_arr, c(1,3), mean) #average out time
    names(dimnames(ss_arr)) <- c("varcomp","variable")
    ss_df <- reshape::melt(ss_arr)
    ss_df
}

