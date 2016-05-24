######################################################################
# Causal conditional inference forests
######################################################################

# x = a data frame of predictors
# y = a binary response (numeric) vector 
# ct = a binary (numeric) vector representing treatment (1) and control group (0)
# mtry = the number of variables to be tested in each node;  default is floor(sqrt(ncol(x))),
# ntree = the number of trees to generate in the forest; default is ntree = 100
# split_method = the split criteria used at each node of each tree; Possible values are: "ED", "Chisq" and "KL"
# interaction.depth = The maximum depth of the interactions
# pvalue = the maximum acceptable pvalue required in order to make a split
# bonferroni = apply bonferroni adjustment to pvalue
# minsplit =  the minimum number of observations that must exist in a node in order for a split to be attempted.
# minbucket_ct0 = the minimum number of control observations in any terminal <leaf> node
# minbucket_ct1 = the minimum number of treatment observations in any terminal <leaf> node 
# keep.inbag = Should an nrow(x) by ntree matrix be returned that keeps track of which samples are
               #“in-bag” in which trees
# verbose = print status messages?
# (...) arguments passed to independence_test{coin}

  ccif <- function(x, ...) UseMethod("ccif")  

  ccif.default <- function(x,  
                           y,  
                           ct, 
                           mtry = floor(sqrt(ncol(x))),
                           ntree = 100, 
                           split_method = c("ED", "Chisq", "KL", "L1", "Int"),
                           interaction.depth = NULL,
                           pvalue = 0.05,
                           bonferroni = FALSE,
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
  
  ### check interaction depth
  if (!is.null(interaction.depth) && interaction.depth < 1)
    stop("uplift: interaction.depth must be greater than 0. Aborting...")
  
  ### check pvalue
  if (!(pvalue >= 0 & pvalue <= 1))
    stop("uplift: pvalue must be between 0 and 1. Aborting...")
  
  ### check bonferroni
  if (!is.logical(bonferroni)) 
    stop("uplift: bonferroni must either be TRUE or FALSE. Aborting...")
  
  ### check number of variables to split on at each node
  if (mtry < 1 || mtry > ncol(x)) 
    stop("uplift: invalid mtry: reset to within valid range. Aborting...")
  
  ### check number of classes
  if (length(unique(y)) != 2)
    stop("uplift: ccif supports only binary response variables. Aborting...")
 
  ### check response coding
  if (!all(unique(y) %in% c(0,1)))
    stop("uplift: y must be coded as 0/1. Aborting...")
  
  ### check treatment coding
  if (!all(unique(ct) %in% c(0,1)))
    stop("uplift: ct must be coded as 0/1. Aborting...")
  
  ### check number of treatments
  if (length(unique(ct)) != 2)
    stop("uplift: ccif supports only 2 treatments at the moment. Aborting...")
  
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
  
  ### Define transformed response variable    
  ### Store data in a data frame
  z <- factor(ifelse((ct == 1 & y == 1) | (ct == 0 & y == 0), 1, 0))
  dframe <- cbind(x, ct, y, z, obs.index = 1:nrow(x))
    
  ### Misc. arguments passed to buildTree (similar to all trees)
  nr_samples <- nrow(dframe)
  nr_vars <- ncol(x)
  var_names <- colnames(x)
  dframe_sp <- split(dframe, list(dframe$y, dframe$ct))
  split_len <- length(dframe_sp)
  if (split_len != 4) 
    stop("uplift: each level of treatment ct must have positive (y=1) and negative (y=0) responses. Aborting...")
  nr_in_samples_sp <- lapply(dframe_sp, nrow)
 
  ### Prob(ct=1) = Prob(ct=0) = 1/2
  p.ct <- prop.table(table(dframe$ct))
  if (p.ct[2] > p.ct[1]) {
    
    nr_in_samples_sp[[3]] <- round(nr_in_samples_sp[[3]] * p.ct[1] / p.ct[2])
    nr_in_samples_sp[[4]] <- round(nr_in_samples_sp[[4]] * p.ct[1] / p.ct[2])
    
  } else 
    
     if (p.ct[2] < p.ct[1]) {
    
    nr_in_samples_sp[[1]] <- round(nr_in_samples_sp[[1]] * p.ct[2] / p.ct[1])
    nr_in_samples_sp[[2]] <- round(nr_in_samples_sp[[2]] * p.ct[2] / p.ct[1])
  }
  
  nr_in_samples <- sum(unlist(nr_in_samples_sp))
  nr_nodes <- 2 * nr_in_samples + 1 #total number of possible nodes
  trees <- vector("list", ntree)
  b.ind.m <- matrix(nrow = nr_in_samples, ncol = ntree, 
                    dimnames = list(NULL, paste('tree', 1:ntree,sep='')))
   
  if (verbose)
    cat("ccif: starting.",date(),"\n")
  
  for (i in 1:ntree) {  
    
    if (verbose)
      if ((i %% 10) == 0 && i < ntree) message( "", i, " out of ", ntree, " trees so far...")    
      
    b.ind <- lapply(1:split_len, function(k) sample(1:nrow(dframe_sp[[k]]),
                                                    nr_in_samples_sp[[k]], replace = TRUE))
    dframe_sp_s <- lapply(1:split_len, function(k) dframe_sp[[k]][b.ind[[k]], ])
    b.dframe <- do.call("rbind", dframe_sp_s)
    
    trees[[i]] <- build_ccit(b.dframe, 
                             mtry, 
                             split_method,
                             interaction.depth,
                             pvalue,
                             bonferroni,
                             minsplit, 
                             minbucket_ct0, 
                             minbucket_ct1, 
                             nr_vars,
                             var_names,
                             nr_in_samples,
                             nr_nodes,
                             ...);
    
    if (keep.inbag) 
      b.ind.m[, i] <- b.dframe$obs.index
  }     
  
  
  cl <- match.call()
  var_class <- sapply(x, class)
  
  res.trees <- list(call = cl,
                    trees = trees, 
                    split_method = split_method,
                    ntree = ntree,
                    mtry = mtry, 
                    var.names = var_names, 
                    var.class = var_class, 
                    inbag = b.ind.m)
                                     
  class(res.trees) <- "ccif"
  return(res.trees) 
}

### END FUN
