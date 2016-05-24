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
## This file contains common tools such as data creation and plotting
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
## Manipulation of lists of data.frames
## -----------------------------------------------------------------------------
# Manipulations that are independet for data.frames in the list can be done
# using lapply(). The more complicated ones are stored here.

#' Remove rows with NA values from a list of data.frames
#' 
#' @description 
#' Same rows removed from all data frames in the list.
#'
#' @param df_lst [1,m] list of data.frames, A list of data.frames to process
#'
#' @return [1,m] list of data.frames, Data without rows that contain NA,
#' same rows removed from all data.frames in the input list.
#'
#' @export
dl_remove_NA <- function(df_lst){
  df <- abind(df_lst, along=2) ## collect all dfs in one giant one
  is.complete <- complete.cases(df)
  df_lst <- lapply(df_lst, function(el){el[is.complete,]})
  df_lst
}


#' Catenate a list of data.frames vertically to a single data.frame
#' 
#' @description  
#' Assumes equal variables for all datasets!
#' Ouput has columns: <variables>, "dataset"
#' Preserves list element names in column "dataset".
#' For a more generic approach see dflst2dfmelt (uses reshape::melt(df_lst))
#'
#' @param df_lst [1,m] list of data.frames, A list of data.frames to process
#' @param id_var_name string, Column name for the dataset id variable
#'
#' @return A data.frame with elements of df_lst catenated vertically, An extra column with dataset id
#' is added.
#'
#' @export
dflst2df <- function(df_lst, id_var_name="dataset"){

  if (is.null(names(df_lst))){
    warning('dflst2df: df_lst entries not named. Making them up...')
    names(df_lst) <- sprintf('dataset%d', 1:length(df_lst))
  }
  
  df_lst <- mapply(function(x,y){x[,id_var_name] <- as.factor(y)
                                    #x[,id_var_name] <- as.factor(x[,id_var_name])
                                    x}, df_lst, names(df_lst),
                    SIMPLIFY=F, USE.NAMES=F)
  ggdf <- do.call(rbind, df_lst)
}


#' Combine a list of data.frames to a single molten data.frame
#' 
#' @description  
#' Output maximally "molten" with columns "dataset", "obs", "variable", "value"
#' Preserves list element names in column "dataset".
#'
#' @param df_lst [1,m] list of data.frames, A list of data.frames to process
#'
#' @return A data.frame with elements of df_lst combined using reshape::melt().
#' 
#' Extra columns:
#' 
#' \describe{
#' \item{dataset}{dataset name}
#' \item{obs}{running observation index (time)}
#' }
#'
#' @importFrom reshape melt.list melt.data.frame
#'
#' @export
dflst2dfmelt <- function(df_lst){
  
  #df <- lapply(df_lst, reshape::melt.data.frame, id.vars = NULL)
  #df <- reshape::melt.list(df, id.vars = c('variable','value'))
  df <- reshape::melt(df_lst) #won't work with importFrom
  
  names(df)[names(df) %in% "L1"] <- "dataset"
  #names(df)[names(df) %in% "L2"] <- "datacollection"
  df$obs <- rep(1:nrow(df_lst[[1]]), length(df_lst))
  df
}


#' Catenate a list of data.frames to a matrix along dim
#'
#' @param df_lst [1,m] list of data.frames, A list of data.frames to process
#' @param dim [1,1] int, Dimension to apply over
#'
#' @return A matrix with elements of df_lst converted to matrix and catenated along dim.
#'
#' @importFrom abind abind
#' @export
dflst2array <- function(df_lst, dim=2){
  df_lst <- dflst_dsnames2varnames(df_lst) 
  tmp <- lapply(df_lst, as.matrix)
  tmp <- do.call(abind::abind, c(tmp, along = dim)) #catenated along dimension dim
}


#' Append dataset names to variable names of the respective dataset
#'
#' @param dflst [1,m] list of data.frames, A list of data.frames to process
#' @param sep string, Separator to use
#'
#' @return A [1,m] list of data.frames with modified variable names.
#'
#' @export
dflst_dsnames2varnames <- function(dflst, sep = "_"){
  for (i in seq.int(length(dflst))){
    names(dflst[[i]]) <- paste(names(dflst)[i], names(dflst[[i]]), sep=sep)
  }
  dflst
}


#' Catenate a set of data collections (lists of data.frames) into a single molted data.frame.
#' 
#' @description 
#' Can be used e.g. to prepare data for plotting with ggplot().
#'
#' @param ... Several lists of data.frames to catenate
#' @param id.vars [1,m] string, ID variables for reshape::melt
#'
#' @return A data.frame with elements of ... melted and catenated vertically into a single data.frame.
#' 
#' Extra columns created:
#' \item{ds:}{dataset id within data collection}
#' \item{dc:}{data collection id}
#' 
#' @examples 
#' df_lst <- list(df1 = iris[,2:3], df2 = iris[2:3])
#' data_collections2ggdf(dc1 = df_lst, dc2 = df_lst)
#'
#' @importFrom reshape melt melt.data.frame
#' 
#' @export
data_collections2ggdf <- function(..., id.vars=NULL){
  input_lst <- list(...)
  ## melt all data.frames
  if (!is.null(id.vars)){
    input_lst <- lapply(input_lst, function(d){
                      lapply(d,function(d2){reshape::melt(d2, id.vars=id.vars)}) })
  } else {
    input_lst <- lapply(input_lst, function(d){ lapply(d,function(d2){reshape::melt(d2)}) })
  }
  ## collapse data collections to data frames
  input_lst <- lapply(input_lst, dflst2df, id_var_name="ds")
  ## collapse into data.frame
  dflst2df(input_lst, id_var_name="dc")
}


#' Apply fun to the bottom level of a nested list structure
#' 
#' @description 
#' Used to batch process computation results that are stored into a nested list structure.
#' Analysis results are stored as lists but with class attribute changed. This signals that the
#' recursion into the list structure should end and fun should be applied instead. Can be used e.g. 
#' to pick out results from a complex list structure.
#'
#' @param lst nested list, A nested list structure to process
#' @param fun function object, The function to apply at the bottom level
#' @param exclude_names string array, Names of list elements to skip at any level
#' @param ... Furhter parameters passed on to fun
#'     
#' @return A list outputs generated when applying fun to the bottom level of input lst.
#' Bottom level is considered reached when something other than class == 'list' is 
#' encountered.      
#'     
#' @export
traverse_nested_list <- function(lst, fun, exclude_names = NULL, ...){
  if (class(lst)=="list"){
    include_match <- !(names(lst) %in% exclude_names)
    lapply(lst[include_match],
           function(el){traverse_nested_list(el, fun=fun, exclude_names=exclude_names, ...)})
  } else {
    # The hack comes in here: results of analyses are lists but their class attribute has been
    # changed to something else. Hence they end up here.
    fun(lst, ...)
  }
}


#' Apply PCA to the data after catenating data.frames horizontally
#'
#' @param df_lst [1,m] list of data.frames, A list of data.frames to process
#' @param center boolean, TRUE -> center data, FALSE -> do nothing
#' @param scale boolean, TRUE -> scale data, FALSE -> do nothing
#' 
#' @return A list with elements:
#' \item{pcdf:}{data.frame, PCA components (prcomp()$x)}
#' \item{model:}{list, Output of prcomp()}
#'
#' @import stats
#' @export
dflst_pca <- function(df_lst, center = F, scale = F){
  d <- dflst2array(df_lst, dim = 2)
  fit <- stats::prcomp(d, center = center, scale = scale)
  pcdf <- as.data.frame(fit$x)
  list(pcdf = pcdf, model = fit)
}


#' Add a data.frame (dataset) to a list of data.frames (datasets)
#'
#' @param dflst [1,m] list of data.frames, A list of data.frames
#' @param df data.frame, Data frame to add
#' @param dsname string, Dataset name for the data.frame to add
#'
#' @return A list of data.frames, A new list of data.frames with one new dataset in the end
#'
#' @export
dflst_add_ds <- function(dflst, df, dsname){
  df <- list(df)
  names(df) <- dsname
  dflst <- c(dflst, df)
  dflst
}


## -----------------------------------------------------------------------------
## Manipulation of data.frames
## -----------------------------------------------------------------------------

#' Apply scale on a numeric data.frame
#'
#' @param df data.frame, A numeric data.frame to process
#' @param ... arguments to scale()
#' 
#' @return data.frame, A scaled data.frame with attributes preserved
#'
#' @export
df_scale <- function(df,...){
  att <- attributes(df)
  df <- as.data.frame(scale(df,...))
  attributes(df) <- att
  df
}

#' Scales variables in data.frame dfx using ordinary least squares such 
#' 
#' @description
#' Scales variables in data.frame dfx using ordinary least squares such that the scaled
#' result explains as much of the variance in dfy as possible.
#' Scaling is done separately for each variable (i.e. no linear mixing of variables).
#' Assumes data.frames dfx and dfy to be of identical structure.
#' Intended use: to scale up cocoreg projections to account for the lost variance.
#'
#' @param dfx data.frame, Data frame to use as independent variable
#' @param dfy data.frame, Data frames to use as dependent variable 
#' 
#' @return data.frame, A rescaled version of dfx with dimnames from dfy.
#' 
#' @examples
#' \dontrun{
#' dc <- create_syn_data_toy()
#' ccr <- cocoreg(dc$data)
#' dfLst <- mapply(df_scale_ols, ccr$data, dc$data , SIMPLIFY=F)
#' }
#' 
#' @export
df_scale_ols <- function(dfx, dfy){
  dfy_est <- as.data.frame(mapply(
    function(x1,x2){
      mapping_lm(as.data.frame(x1),as.data.frame(x2))(as.data.frame(x1))
    }, dfx, dfy))
  dimnames(dfy_est) <- dimnames(dfy)
  dfy_est
}


## -----------------------------------------------------------------------------
## Statistics
## -----------------------------------------------------------------------------

#' Compute Euclidean norm of vector
#' 
#' @description
#' Convenience function for use with e.g. lapply
#' 
#' @param x [1,m] double, A vector of data
#' 
#' @return [1,1] double, Euclidiean norm of x
#' 
#' @export
vecnorm <- function(x){
  base::norm(as.matrix(x),'F')
}

#' Make vector of unit norm
#' 
#' @param x [1,m] double, A vector of data
#' 
#' @return [1,m] double, Same vector normalized to unit Euclidean norm
#' 
#' @export
to_unit_vec <- function(x){
  x / vecnorm(x)
}

#' Compute RMSE between vectors v1 and v2
#'
#' @param v1 [1,m] numeric, First data vector
#' @param v2 [1,m] numeric, Second data vector
#' @param relative boolean, If TRUE, relate the rmse value to the rmse of v1.
#'        If FALSE, just compute RMSE between v1 and v2
#'
#' @return [1,1] double, RMSE value
#'
#' @export
rmse <- function(v1, v2, relative=F){
  if (relative){
    rmse(v1,v2)/rmse(v1, rep(0,length(v1)) )
  } else {
    sqrt( mean((v1-v2)^2, na.rm = TRUE) )  
  }
}


#' Compute RMSE between data.matrices dm1 and dm2
#' 
#' @description
#' 
#' A data.matrix has observations as rows and variables as columns
#'
#' @param dm1 [N,M] numeric, First data.matrix
#' @param dm2 [N,M] numeric, Second data.matrix
#'
#' @return [1,M] numeric, A vector of RMSE values, one per variable.
#'
#' @examples 
#' \dontrun{
#'  dm1 <- matrix(rep(1,6),nrow=2)
#'  dm2 <- matrix(rep(3,6),nrow=2)
#'  data_matrix_rmse(dm1, dm2)
#'  
#'  first = list(dm1, dm1)
#'  second = list(dm2, dm2)
#'  (tmp = mapply(data_matrix_rmse, first, second, SIMPLIFY=FALSE))
#'  }
#'
#' @export
data_matrix_rmse <- function(dm1, dm2){
  ## average the rows of data.matrix, columns (variables) remain
  sqrt( apply((dm1-dm2)^2, 2, mean , na.rm = TRUE) )
}


#' Standard error of mean
#'
#' @param x [1,M] numeric, data vector
#' @param na.rm procedure for NA's, passed on to sd(), default: na.rm = T
#'
#' @return [1,1] numeric, standard error of mean
#'
#' @export
se <- function(x, na.rm = T){
  sd(x, na.rm=na.rm)/sqrt(length(x))
}


#' Get [min, max] of a list of numeric objects
#'
#' @param dataList [1,m] list of numeric objects
#' 
#' @return [1,2] double, [min, max] of the input 
#'
#' @export
get_range_datalist <- function(dataList){
  ## TODO: rename to dl_range()
  dl_max <- max(as.numeric( lapply(dataList, function(x){max(x)}) ))
  dl_min <- min(as.numeric( lapply(dataList, function(x){min(x)}) ))
  c(dl_min, dl_max)
}


#' A rotation matrix
#'
#' @param angle_deg [1,1] numeric, Angle in degrees
#' 
#' @return [2,2] matrix, Rotation matrix for making angle_deg 2D rotation 
#' 
#' @export
rotation_matrix <- function(angle_deg){
  angle_rad = (angle_deg/180)*pi
  matrix(c(cos(angle_rad), -sin(angle_rad),
           sin(angle_rad), cos(angle_rad)), ncol=2 ,byrow=T)
}


#' Sum-of-squares values showing what portion of variance in dvec is explained
#' by dvec_est
#' 
#' @description
#'
#' Computation as in:
#' http://en.wikipedia.org/wiki/Fraction_of_variance_unexplained
#'
#' ss_est becomes zero if dvec_est equals dvec_0=rep(mean(dvec),length(dvec)).
#' If dvec_est is better estimate than dvec_0, R2 is positive.
#' If dvec_est is worse than dvec_0, R2 is negative.
#'
#' @param dvec [1,m] numeric, data vector
#' @param dvec_est [1,m] numeric, data vector, an estimate of dvec
#' 
#' @return A list with elements:
#' \item{ss_tot:}{Sum of squares in dvec}
#' \item{ss_est:}{Sum of squares in dvec_est}
#' \item{ss_err:}{Sum of squares of dvec - dvec_est}
#' \item{R2:}{Percentage of variance explained i.e. 1 - ss_err/ss_tot}
#'
#' @export
var_explained <- function(dvec, dvec_est){
  
  ss_tot <- ss(dvec - mean(dvec))
  ss_est <- ss(dvec_est - mean(dvec))
  ss_err <- ss(dvec - dvec_est)
  R2 = 1 -  ss_err/ss_tot

  list(ss_tot=ss_tot, ss_est=ss_est, ss_err=ss_err, R2=R2)
}


#' Compute total, within group and between group variability using fun
#' 
#' @description
#' The function used the definition: gvar = tvar - wgvar
#'
#' @param vec [1,M] numeric, Data vector
#' @param grp [1,M] integer/character vector, Some grouping of vec
#' @param fun function, Function to use when quantifying the variability
#'
#' @return A list with elements:
#' \item{tvar:}{Total variability}
#' \item{bgvar:}{Between groups variability, tvar - sum(wgvar_*)}
#' \item{wgvar_<groupname>:}{Within group variability for each group}
#' \item{wg_rel:}{sum(wgvar)/tvar}
#' \item{bg_rel:}{bgvar/tvar}
#'
#' @examples
#' vec <- rnorm(10)
#' grp <- rep(c("a","b","c"), c(3,3,4))
#' variability_components(vec, grp, ss)
#' 
#' @export
variability_components <- function(vec, grp, fun){
  grpmean_vec <- tapply(vec, grp, mean)
  tvar <- fun(vec - mean(vec))
  wgvar <- vapply(seq.int(length(grpmean_vec)),
                  function(i){ fun(vec[grp %in% names(grpmean_vec)[i]] - grpmean_vec[i]) },
                  numeric(1))
  names(wgvar) <- names(grpmean_vec)
  bgvar <- tvar - sum(wgvar) ## this term cannot be computed directly as
  ## there might be "interactions" present
  
  res <-        c(tvar  , bgvar  , wgvar                              , sum(wgvar)/tvar, bgvar/tvar)
  names(res) <- c("tvar", "bgvar", paste0("wgvar_",names(grpmean_vec)), "wg_rel"       , "bg_rel")
  res
}



#' Sum of squares
#' 
#' @param x [1,m] numeric, A data vector
#'
#' @return Sum of squares of x
#' 
#' @export
ss <- function(x){sum(x^2)}


## -----------------------------------------------------------------------------
## Data generation
## -----------------------------------------------------------------------------

#' Make 2D gauss data (maybe obsolete)
#'
#' @param n [1,1] int, Number of observations
#' @param var [1,2] numeric, Variances
#' @param angle_deg [1,1] numeric, Rotation angle
#' @param scale boolean, Scale data? T -> scale, F -> do not scale
#' @param seed [1,1] int, Random seed
#' 
#' @return Matrix of 2D gaussian data
#' 
#' @export
make_data_gauss_2d <- function(n, var, angle_deg, scale=T, seed=42){
  R <- rotation_matrix(angle_deg) #rotation matrix
  L = diag(var) #variances
  set.seed(seed)
  S = matrix(rnorm(n*2), ncol=n)
  d = t(R %*% L %*% S) 
  if (scale){
    d = scale(d)  
  } else { d }
}


#' Add notch-like gaussian snippets to an existing signal x
#'
#' @param x [1,N] numeric, Original data
#' @param pos [1,m] integer, Positions to add notches to
#' @param sd [1,1] numeric, (optional) Desired width of the Gaussian notch
#' @param amplitude [1,1] numeric, (optional) Desired amplitude for the notches
#' 
#' @return [1,N] numeric, Modified signal with notches
#'
#' @export
add_notches <- function(x, pos, sd=0.01*length(x), amplitude=1){
  
  N <- length(x)
  mid <- floor(N/2)
  
  nv <- dnorm(-mid:(mid-1), sd=sd)
  if (length(nv) != N) {
    # can only be one shorter: for uneven N
    nv <- c(nv,0)
  }
  nv <- (amplitude/max(nv))*nv #scale peak height
  
  # Add notches
  for (p in pos){
    if (p > mid){
      by = p - mid -1
    } else {
      by = -(mid-p+1)
    }
    x <- x + cshift(nv, by=by)    
  }
  x
}


## -----------------------------------------------------------------------------
## General utilities
## -----------------------------------------------------------------------------

#' Circularly shift vector elements
#'
#' @param x [1,N] numeric, A vector
#' @param by [1,1] integer, How many positions to shift.
#'        by > 0 -> shift to right
#'        by = 0 -> no shift
#'        by < 0 -> shift to left
#' 
#' @return [1,N] numeric, Circularly shifted signal
#'
#' @export
cshift <- function(x, by){
  N <- length(x)
  if (by > 0){
    # to right
    c(x[(N-by+1):N], x[1:(N-by)])  
  } else if (by==0){ x }
  else {
    # to left
    c(x[abs(by)+1:(N-abs(by))], x[1:(abs(by))])
  }
}


#' Replicate matrix to create a larger one
#' 
#' @description
#' From:
#' http://haky-functions.blogspot.fi/2006/11/repmat-function-matlab.html (accessed 27.3.2015)
#'
#' @param X A [I,J] matrix or J element vector, Matrix used as such, vector coerced to a row
#'        matrix with dim(X)=[1,J].
#' @param m [1,1] integer, Replication count vertically
#' @param n [1,1] integer, Replication count horizontally
#' 
#' @return [m*I,n*J] matrix, Replicated data
#'
#' @export
repmat = function(X, m, n){
  ##R equivalent of repmat (matlab)
  if (class(X) != "matrix"){
    X = matrix(X, ncol=length(X))
  }
  
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}


## -----------------------------------------------------------------------------
## Wrappers
## -----------------------------------------------------------------------------

#' Run BGFA by Klami et al using data format conventions of this repo
#'
#' @param df_list [1,m] list of data.frames, Input data to GFA in COCOREG format
#' @param K [1,1] int, The (maximum) number of components; should be high enough to capture all of the components. This can be recognized by at least a few of the components being shut down
#' @param Nrep [1,1] int, Number of random initializatin used for learning the model
#' 
#' @return A list, The output of CCAGFA::GFA()
#'
#' @import CCAGFA      
#' @export
wrapper_BGFA <- function(df_list, K=8, Nrep=2){
  dmat_list <- lapply(df_list, as.matrix)
  
  opts <- CCAGFA::getDefaultOpts()
  ## These options are from the examples in CCAGFA
  opts$iter.crit <- 1e-6     # Need to make this more strict if
  opts$lbfgs.factr <- 1e5    # A bit less strict convergence criterion,
  opts$iter.max <- 1e3
  #   having a large sample size
  #   should be enough for our simple dataset
  opts$verbose <- 1          # Looks a bit nicer
  opts$low.mem <- FALSE
  
  #model <- CCAGFA::GFAexperiment(dmat_list, K, opts, Nrep=Nrep)
  model <- CCAGFA::GFA(dmat_list, K, opts)
}


#' Project BGFA components common to all datasets back to the original space
#' 
#' @param model Output of CCAGFA::GFA()
#' @param threshold [1,1] double, GFA model trimming threshold
#' 
#' @return A list of data.frames, Original data reconstructed using only 
#' those latent components that are active in all datasets
#' 
#' @import CCAGFA     
#' @export
BGFA_joint_info <- function(model, threshold=0.001){
  model.trim <- CCAGFA::GFAtrim(model, threshold=threshold) #0.001 GFAtrim default threshold
  all_active_match <- colSums(model.trim$active)==nrow(model.trim$active)
  lapply(model.trim$W, function(w){ as.data.frame(model.trim$Z[,all_active_match] %*%
                                                      t(w[,all_active_match])) })
}

#' Apply GFA using the same interface as cocoreg()
#' 
#' @description 
#' Note: if K is too high GFA() might not converge in a meaningful time
#' or the computation may mysteriously crash.
#' 
#' @param df_list [1,m] list of data.frames, Input data to GFA in COCOREG format
#' @param K [1,1] int, (Maximum) number of GFA components 
#' @param Nrep [1,1] int, Number of random initializatin used for learning the model
#' @param threshold [1,1] double, GFA model trimming threshold
#' 
#' @return A list with elements:
#' \item{$data:}{[1,m] list of data.frames, Original data reconstructed using only 
#'         those latent components that are active in all datasets}
#' \item{$model:}{a list, Non-trimmed output of CCAGFA::GFA()}      
#' \item{$dataid:}{string, Dataset identifier string}
#' \item{$method:}{string, Analysis method identifier string}
#' \item{$wall_time_taken:}{[1,1] double, Time taken to run the analysis in seconds}
#' 
#' @export
BGFA_cocoreg_interface <- function(df_list, K=8, Nrep=2, threshold=0.001){
  start_time <- Sys.time()
  
  if ('id' %in% names(attributes(df_list)) ){
    dataid <- attr(df_list,'id') #get before its gone  
  } else {
    dataid <- 'NA'
  }
  df_list <- validate_data(df_list)
  dmeta <- get_dc_meta(df_list)

  df_list  <- lapply(df_list, function(el){df_scale(el, center=T, scale=F)})
  
  ## A quick-n-dirty fix for one special dataset that triggers a bug in GFA
  if ( length(grep("syn", dataid)) == 1 ){
    #synthetic data -> reduce K to be able to compute
    model <- wrapper_BGFA(df_list, K=4, Nrep=Nrep)  
  } else {
    model <- wrapper_BGFA(df_list, K=K, Nrep=Nrep)  
  }
  data <- BGFA_joint_info(model, threshold=threshold)

  ## Convert back to original data format
  data <- apply_dc_meta(data, dmeta)   

  res_lst <- list(data=data,
                  model=model,
                  dataid=dataid,
                  method='gfa')
  res_lst <- c(res_lst, list(wall_time_taken = as.numeric(difftime(Sys.time(), start_time, units = "sec"))))
  class(res_lst) <- "gfa" ## Needed for effective recursive traversing of results in a nested list
  res_lst
}


#' PCA projection using cocoreg interface
#' 
#' @param df_list [1,m] list of data.frames, Input data to GFA in COCOREG format
#' @param prc_th [1,1] double, Threshold in precentage of cumulative variance explained
#'        PCA components are included until cumulative explained variance reaches prc_th.
#'       
#' @return A list with elements:         
#' \item{$data}{[1,m] list of data.frames, Original data reconstructed using only 
#'         those latent components that are active in all datasets}
#' \item{$dataid}{string, Dataset identifier string}
#' \item{$method}{string, Analysis method identifier string}
#' \item{$wall_time_taken}{[1,1] double, Time taken to run the analysis in seconds}
#' 
#' @export
PCA_cocoreg_interface <- function(df_list, prc_th=0.9){
  start_time <- Sys.time()
  dataid <- attr(df_list,'id') #get before its gone
  df_list <- validate_data(df_list)
  dmeta <- get_dc_meta(df_list)
  
  df_list  <- lapply(df_list, function(el){df_scale(el, center=T, scale=F)})
  
  df <- do.call(cbind, df_list)
  varnames <- lapply(seq.int(length(df_list)),
                     function(i){
                       paste0(names(df_list[i]),".",colnames(df_list[[i]]))
                       })
  
  fit <- stats::prcomp(df, center=T, scale=T) #SVD
  ##summary(fit) # print variance accounted for 
  pc_inds <- which(summary(fit)$importance[3,] < prc_th)
  if (length(pc_inds)==0){
    pc_inds = 1
  } else {
    pc_inds <- as.numeric(c(pc_inds, max(pc_inds)+1))  
  }
  pcd <- as.matrix(df) %*% fit$rotation[,pc_inds]
  pd <- pcd %*% t(fit$rotation[,pc_inds])
  pd <- lapply(varnames, function(el){as.data.frame(pd[,el])})
  
  ## Convert back to original data format
  pd <- apply_dc_meta(pd, dmeta) 
  
  res_lst <- list(data=pd,
                  dataid=dataid,
                  method='pca',
                  wall_time_taken = as.numeric(difftime(Sys.time(), start_time, units = "sec")) )
  class(res_lst) <- "pca" ## Needed for effective recursive traversing of results in a nested list
  res_lst
}


#' SCA projection using cocoreg interface
#' 
#' @param df_list, [1,m] list of data.frames, Input data to GFA in COCOREG format
#' @param nfac, [1,1] int, see multiway::sca() for details
#' @param type, string, Type of analysis, see multiway::sca() for details
#' @param rotation, string, see multiway::sca() for details
#' @param nstart, [1,1] int, see multiway::sca() for details
#' 
#' @return A list with elements:
#' \item{$data}{[1,m] list of data.frames, Original data reconstructed using only 
#'         those latent components that are active in all datasets}
#' \item{$model}{list, The output of multiway::sca()}
#' \item{$dataid}{string, Dataset identifier string}
#' \item{$method}{string, Analysis method identifier string}
#' \item{$wall_time_taken}{[1,1] double, Time taken to run the analysis in seconds}
#' 
#' @importFrom multiway sca
#' @export

SCA_cocoreg_interface <- function(df_list,
                                  nfac=1, type="sca-p",
                                  rotation="none", nstart=10){
  start_time <- Sys.time()
  
  ## Extract metadata
  dataid <- attr(df_list,'id') #get before its gone
  df_list <- validate_data(df_list)
  dmeta <- get_dc_meta(df_list)
  
  df_list  <- lapply(df_list, function(el){df_scale(el, center=T, scale=F)})
  
  varnames <- lapply(seq.int(length(df_list)),
                     function(i){
                       paste0(names(df_list[i]),".",colnames(df_list[[i]]))
                     })
  
  ## SCA
  dm_list <- lapply(df_list, as.matrix)
  scamod <- multiway::sca(dm_list, nfac=nfac, type=type, rotation=rotation, nstart=nstart)
  
  # Projection onto original variables
  sca_comp <- vector(length=length(df_list), mode="list")
  for (i in seq.int(length(df_list))){
    sca_comp[[i]] <- scamod$D[[i]] %*% t(scamod$B)
  }
  
  ## Convert back to original data format
  sca_comp <- lapply(sca_comp, as.data.frame)
  sca_comp <- apply_dc_meta(sca_comp, dmeta) 
  
  res_lst <- list(data=sca_comp,
                  model=scamod,
                  dataid=dataid,
                  method='sca',
                  wall_time_taken = as.numeric(difftime(Sys.time(), start_time, units = "sec")) )
  class(res_lst) <- "sca" ## Needed for effective recursive traversing of results in a nested list
  res_lst
}


#' COCOREG style analysis using RGCCA projection
#' 
#' @description
#' COCOREG interface used for both input and output.
#'
#' @param dflst, [1,m] list of data.frames, Input data to GFA in COCOREG format
#' @param tauArr, [1,m] double, See RGCCA::rgcca() for details
#' 
#' @return A list with elements:
#' \item{$data}{[1,m] list of data.frames, Original data reconstructed using only 
#'         those latent components that are active in all datasets}
#' \item{$model}{list, The output RGCCA::rgcca()}
#' \item{$dataid}{string, Dataset identifier string}
#' \item{$method}{string, Analysis method identifier string}
#' \item{$wall_time_taken}{[1,1] double, Time taken to run the analysis in seconds}
#' 
#' @import RGCCA
#' @export
RGCCA_cocoreg_interface <- function(dflst, tauArr = rep(0.5, length(dflst))) {
  start_time <- Sys.time()
  
  ## Extract metadata
  dataid <- attr(dflst,'id') #get before its gone
  df_list <- validate_data(dflst)
  dmeta <- get_dc_meta(dflst)
  
  # Get the first RGCCA component
  res.rgcca = RGCCA::rgcca(  A = lapply(dflst, as.matrix),
                      ncomp = rep(1,length(dflst)),
                      tau = tauArr,
                      verbose = F)
  rgccaProjs <- lapply(res.rgcca$Y, as.data.frame)
  rgccaEstimate <- mapply(df_scale_ols, rgccaProjs, dflst , SIMPLIFY=F)
  
  res_lst <- list(data = rgccaEstimate,
                  model = res.rgcca,
                  dataid = dataid,
                  method = 'rgcca',
                  wall_time_taken = as.numeric(difftime(Sys.time(), start_time, units = "sec")) )
  class(res_lst) <- "rgcca" ## Needed for effective recursive traversing of results in a nested list
  res_lst
}

