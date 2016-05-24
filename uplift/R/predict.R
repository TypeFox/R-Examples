######################################################################
# Predict upliftRF
######################################################################

predict.upliftRF <- function(object, newdata, n.trees = object$ntree, predict.all = FALSE, ...) {
  
  if (!inherits(object, "upliftRF"))
    stop("uplift: object not of class upliftRF")
  if (is.null(object$ntree))
    stop("uplift: no trees in the object")
  if (nrow(newdata) == 0)
    stop("uplift: newdata has 0 rows")
  if (any(!object$var.names %in% colnames(newdata)))
    stop("uplift: variables in the training data missing in newdata")
  if (n.trees > object$ntree)
    stop("uplift: n.trees cannot be greater than the number of fitted trees")
  if (n.trees < 1)
    stop("uplift: n.trees must be greater than 0")
  
  
  # rewrite newdata to match training data variable indexes
  newdata <- model.frame(terms(reformulate(object$var.names)),
                         newdata, na.action = na.pass) 
  
  if (any(is.na(newdata)))
    stop("uplift: missing values in newdata")
  
  var.class.new <- sapply(newdata, class)
  if (!all(object$var.class == var.class.new))
    stop("uplift: type of predictors in new data do not match that of the training data")
  
  nr_samples_t <- nrow(newdata)
  lx <- n.trees 
  pred.trees <- vector("list", lx)
  
  for (i in 1:lx) {
    obs_node_t <- rep(1, nr_samples_t) # initilize which obs belong to which node
    ### Assign observations to terminal nodes
    for (curr_node_t in 1:object$trees[[i]]$total_nr_nodes) { 
      obs_curr_node.ind_t <- which(obs_node_t == object$trees[[i]]$s_curr_node[curr_node_t]) 
      if (object$trees[[i]]$s_node_status[curr_node_t] == 1) {
        if (is.numeric(object$trees[[i]]$s_bs.x.value[[curr_node_t]])) {
          obs_node_t[obs_curr_node.ind_t] <- ifelse(newdata[obs_curr_node.ind_t, object$trees[[i]]$s_bs.var[curr_node_t]] <= 
                                                      object$trees[[i]]$s_bs.x.value[[curr_node_t]],
                                                    2 * object$trees[[i]]$s_curr_node[curr_node_t],
                                                    2 * object$trees[[i]]$s_curr_node[curr_node_t] + 1)
        } else { # if logical
          obs_node_t[obs_curr_node.ind_t] <- ifelse(newdata[obs_curr_node.ind_t, object$trees[[i]]$s_bs.var[curr_node_t]] %in%
                                                      names(object$trees[[i]]$s_bs.x.value[[curr_node_t]][object$trees[[i]]$s_bs.x.value[[curr_node_t]] == TRUE]),
                                                    2 * object$trees[[i]]$s_curr_node[curr_node_t],
                                                    2 * object$trees[[i]]$s_curr_node[curr_node_t] + 1)     
        }
      }
    }
    
    ### Matrix of response probabilities
    obs_node_tt <- sapply(obs_node_t, function(x) which(x == object$trees[[i]]$s_curr_node))
    pred <- cbind(pr.y1_ct1 = object$trees[[i]]$s_pr.y1_ct1[obs_node_tt], pr.y1_ct0 = object$trees[[i]]$s_pr.y1_ct0[obs_node_tt])
    pred.trees[[i]] <- list(obs_node_t = obs_node_t, 
                            pred = pred)
  }
  ### compute average prediction over all trees  
  pred.trees.t <- lapply(pred.trees, function(x) x$pred) # extract predictions from all trees
  pred.sum <- matrix(rep(0, nr_samples_t) * 2, nrow = nr_samples_t, ncol = 2, 
                     dimnames = list(NULL, c("pr.y1_ct1", "pr.y1_ct0")))
  for (i in 1:length(pred.trees.t)) {
    pred.sum.temp <- pred.trees.t[[i]]
    pred.sum <- pred.sum + pred.sum.temp 
  }
  pred.avg <- pred.sum / length(pred.trees.t)
  if (predict.all) {
    all.res <- list(individual = pred.trees.t, pred.avg = pred.avg) } else {
      all.res <- pred.avg   
    }
  return(all.res)
}

### END FUN


######################################################################
# Predict ccif
######################################################################

predict.ccif <- function(object, newdata, n.trees = object$ntree, predict.all = FALSE, ...) {
  
  if (!inherits(object, "ccif"))
    stop("uplift: object not of class ccif")
  if (is.null(object$ntree))
    stop("uplift: no trees in the object")
  if (nrow(newdata) == 0)
    stop("uplift: newdata has 0 rows")
  if (any(!object$var.names %in% colnames(newdata)))
    stop("uplift: variables in the training data missing in newdata")
  if (n.trees > object$ntree)
    stop("uplift: n.trees cannot be greater than the number of fitted trees")
  if (n.trees < 1)
    stop("uplift: n.trees must be greater than 0")
  
  
  # rewrite newdata to match training data variable indexes
  newdata <- model.frame(terms(reformulate(object$var.names)),
                         newdata, na.action = na.pass) 
  
  if (any(is.na(newdata)))
    stop("uplift: missing values in newdata")
  
  var.class.new <- sapply(newdata, class)
  if (!all(object$var.class == var.class.new))
    stop("uplift: type of predictors in new data do not match that of the training data")
  
  nr_samples_t <- nrow(newdata)
  lx <- n.trees 
  pred.trees <- vector("list", lx)
  
  for (i in 1:lx) {
    obs_node_t <- rep(1, nr_samples_t) # initilize which obs belong to which node
    ### Assign observations to terminal nodes
    for (curr_node_t in 1:object$trees[[i]]$total_nr_nodes) { 
      obs_curr_node.ind_t <- which(obs_node_t == object$trees[[i]]$s_curr_node[curr_node_t]) 
      if (object$trees[[i]]$s_node_status[curr_node_t] == 1) {
        if (is.numeric(object$trees[[i]]$s_bs.x.value[[curr_node_t]])) {
          obs_node_t[obs_curr_node.ind_t] <- ifelse(newdata[obs_curr_node.ind_t, object$trees[[i]]$s_bs.var[curr_node_t]] <= 
                                                      object$trees[[i]]$s_bs.x.value[[curr_node_t]],
                                                    2 * object$trees[[i]]$s_curr_node[curr_node_t],
                                                    2 * object$trees[[i]]$s_curr_node[curr_node_t] + 1)
        } else { # if logical
          obs_node_t[obs_curr_node.ind_t] <- ifelse(newdata[obs_curr_node.ind_t, object$trees[[i]]$s_bs.var[curr_node_t]] %in%
                                                      names(object$trees[[i]]$s_bs.x.value[[curr_node_t]][object$trees[[i]]$s_bs.x.value[[curr_node_t]] == TRUE]),
                                                    2 * object$trees[[i]]$s_curr_node[curr_node_t],
                                                    2 * object$trees[[i]]$s_curr_node[curr_node_t] + 1)     
        }
      }
    }
    
    ### Matrix of response probabilities
    obs_node_tt <- sapply(obs_node_t, function(x) which(x == object$trees[[i]]$s_curr_node))
    pred <- cbind(pr.y1_ct1 = object$trees[[i]]$s_pr.y1_ct1[obs_node_tt], pr.y1_ct0 = object$trees[[i]]$s_pr.y1_ct0[obs_node_tt])
    pred.trees[[i]] <- list(obs_node_t = obs_node_t, 
                            pred = pred)
  }
  ### compute average prediction over all trees  
  pred.trees.t <- lapply(pred.trees, function(x) x$pred) # extract predictions from all trees
  pred.sum <- matrix(rep(0, nr_samples_t) * 2, nrow = nr_samples_t, ncol = 2, 
                     dimnames = list(NULL, c("pr.y1_ct1", "pr.y1_ct0")))
  for (i in 1:length(pred.trees.t)) {
    pred.sum.temp <- pred.trees.t[[i]]
    pred.sum <- pred.sum + pred.sum.temp 
  }
  pred.avg <- pred.sum / length(pred.trees.t)
  if (predict.all) {
    all.res <- list(individual = pred.trees.t, pred.avg = pred.avg) } else {
      all.res <- pred.avg   
    }
  return(all.res)
}

### END FUN
