######################################################################
# Uplift random forest
######################################################################

# x = a data frame of predictors
# y = a binary response (numeric) vector 
# ct = a binary (numeric) vector representing treatment (1) and control group (0)
# mtry = the number of variables to be tested in each node;  default is floor(sqrt(ncol(x))),
# ntree = the number of trees to generate in the forest; default is ntree = 100
# split_method = the split criteria used at each node of each tree; Possible values are: "ED", "Chisq" and "KL"
# minsplit =  the minimum number of observations that must exist in a node in order for a split to be attempted.
# minbucket_ct0 = the minimum number of control observations in any terminal <leaf> node
# minbucket_ct1 = the minimum number of treatment observations in any terminal <leaf> node 
# bag.fraction = the fraction of the training set observations randomly selected for the purpose of fitting each
                 #tree in the forest
# keep.inbag = Should an nrow(x) by ntree matrix be returned that keeps track of which samples are
               #“in-bag” in which trees
# verbose = print status messages?

  upliftRF <- function(x, ...) UseMethod("upliftRF")  

  upliftRF.default <- function(x,  
                               y,  
                               ct, 
                               mtry = floor(sqrt(ncol(x))),
                               ntree = 100, 
                               split_method = c("ED", "Chisq", "KL", "L1", "Int"),
                               interaction.depth = NULL,
                               bag.fraction = 0.5,
                               minsplit = 20, 
                               minbucket_ct0 = round(minsplit/4), 
                               minbucket_ct1 = round(minsplit/4), 
                               keep.inbag = FALSE,
                               verbose = FALSE,
                               ...) {
                     
  ### check classes of training data 
  if(!is.data.frame(x))
    stop("uplift: x must be data frame. Aborting...")
  
  if(!is.numeric(y))
    stop("uplift: y must be a numeric vector. Aborting...")
  
  if(!is.numeric(ct))
    stop("uplift: ct must be a numeric vector. Aborting...")
  
  ### make sure we have valid data
  if (any(is.na(as.vector(x))) |
      any(is.infinite(as.matrix(x))) |
      any(is.na(as.vector(y))) |
      any(is.infinite(as.vector(y))) |
      any(is.na(as.vector(ct))) |
      any(is.infinite(as.vector(ct))))
    stop("uplift: training data contains NaNs or Inf values. Please correct and try again. Aborting...")
  
  ### check split method
  out.method <- charmatch(tolower(split_method), 
                          c("ed", "chisq", "kl", "l1", "int"))
  if (is.na(out.method))
    stop("uplift: split_method must be one of 'ED', 'Chisq', 'KL', 'L1' or 'Int'. Aborting...")
  
  ### check bag fraction
  if (bag.fraction <= 0 || bag.fraction > 1)
    stop("uplift: bag.fraction must be greater than 0 and equal or less than 1. Aborting...")
  
  ### check interaction depth
  if (!is.null(interaction.depth) && interaction.depth < 1)
    stop("uplift: interaction.depth must be greater than 0. Aborting...")
  
  ### check number of variables to split on at each node
  if (mtry < 1 || mtry > ncol(x)) 
    stop("uplift: invalid mtry: reset to within valid range. Aborting...")
  
  ### check number of classes
  if (length(unique(y)) != 2)
    stop("uplift: upliftRF supports only binary response variables. Aborting...")
 
  ### check response coding
  if (!all(unique(y) %in% c(0,1)))
    stop("uplift: y must be coded as 0/1. Aborting...")
  
  ### check treatment coding
  if (!all(unique(ct) %in% c(0,1)))
    stop("uplift: ct must be coded as 0/1. Aborting...")
  
  ### check number of treatments
  if (length(unique(ct)) != 2)
    stop("uplift: upliftRF supports only 2 treatments at the moment. Aborting...")
  
  ### check length of training data 
  if (length(y) != nrow(x) || nrow(x) != length(ct))
   stop("uplift: length of x, y, and ct must be similar. Aborting...")
  
  ### check number of factor levels is no more than 32
  xlevels <- lapply(x, mylevels)
  ncat <- sapply(xlevels, length)
  maxcat <- max(ncat)
  if (maxcat > 32)
    stop("uplift: can not handle categorical predictors with more than 32 categories. Aborting...")
  
  ### Check minbucket
  if (minbucket_ct0 < 1 || minbucket_ct1 < 1)
    stop("uplift: minbucket_ct0 and minbucket_ct1 must be greater than 0. Aborting...")
  
  ### check verbose setting 
  if (verbose)
    message("uplift: status messages enabled; set \"verbose\" to false to disable")
      
  ### Store all data in a data frame
  dframe <- cbind(x, ct, y, obs.index = 1:nrow(x))
    
  ### Misc. arguments passed to buildTree (similar to all trees)
  nr_samples <- nrow(dframe)
  nr_vars <- ncol(x)
  dframe_sp <- split(dframe, list(dframe$y, dframe$ct))
  split_len <- length(dframe_sp)
  if (split_len != 4) 
    stop("uplift: each level of treatment ct must have positive (y=1) and negative (y=0) responses. Aborting...")
  nr_in_samples_sp <- lapply(dframe_sp, function(x) floor(nrow(x) * bag.fraction))
  nr_in_samples <- sum(unlist(nr_in_samples_sp))
  nr_nodes <- 2 * nr_in_samples + 1 #total number of possible nodes
  trees <- vector("list", ntree)
  b.ind.m <- matrix(nrow = nr_in_samples, ncol = ntree, 
                    dimnames = list(NULL, paste('tree', 1:ntree,sep='')))
  
  
  if (verbose)
    cat("upliftRF: starting.",date(),"\n")
  
  for (i in 1:ntree) {  
    
    if (verbose)
      if ((i %% 10) == 0 && i < ntree) message( "", i, " out of ", ntree, " trees so far...")    
      
    b.ind <- lapply(1:split_len, function(k) sample(1:nrow(dframe_sp[[k]]),
                                                    nr_in_samples_sp[[k]], replace = FALSE))
    dframe_sp_s <- lapply(1:split_len, function(k) dframe_sp[[k]][b.ind[[k]], ])
    b.dframe <- do.call("rbind", dframe_sp_s)
    
    trees[[i]] <- buildTree(b.dframe, 
                            mtry, 
                            split_method,
                            interaction.depth,
                            minsplit, 
                            minbucket_ct0, 
                            minbucket_ct1, 
                            nr_vars,
                            nr_in_samples,
                            nr_nodes);
    
    if (keep.inbag) 
      b.ind.m[, i] <- b.dframe$obs.index
  }     
  
  
  cl <- match.call()
  var.names <- colnames(x)
  var.class <- sapply(x, class)
  
  res.trees <- list(call = cl,
                    trees = trees, 
                    split_method = split_method,
                    ntree = ntree,
                    mtry = mtry, 
                    var.names = var.names, 
                    var.class = var.class, 
                    inbag = b.ind.m)
                                     
  class(res.trees) <- "upliftRF"
  return(res.trees) 
}

### END FUN
