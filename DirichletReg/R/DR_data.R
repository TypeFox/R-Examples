DR_data <- function(Y,                                   # response (compositional variable)
                    trafo = sqrt(.Machine$double.eps),   # transform (compress) the data?
                    base  = 1,                           # base variable for the reparametrized (alternative) model
                    norm_tol = sqrt(.Machine$double.eps) # tolerance for normalization [0, ?]
                   ){


# initialization
  force.norm <- FALSE   # was normalization forced?
  state.tran <- FALSE   # was Y transformed?
  force.tran <- FALSE   # was transformation forced?

# set all rows containing NAs to NA
  if(any(is.na(Y))){ Y[which(rowSums(is.na(Y)) > 0),] <- NA }

# check for negative values in Y
  if(any(na.delete(Y) < 0)) stop('"Y" contains values < 0.')

# save the original data for reference
  Y.original <- Y



### CONVENIENTLY HANDLE BETA-DISTRIBUTED VARIABLES
  if((!is.matrix(Y) && !is.data.frame(Y)) || ifelse(is.null(ncol(Y)), FALSE, ncol(Y) == 1L)){
    if(any((na.delete(Y) < 0) || (na.delete(Y) > 1))){
      stop('only one variable supplied with values outside [0, 1].\nbeta distribution cannot safely be assumed.\nprepare your data first.')
    }
    Y <- cbind(1.0-Y, Y)

    .name <- deparse(match.call()$Y)
    .name <- gsub(".*\\$", "", .name)   # eliminate references to the object the variable comes from
    .name <- gsub("\\[.*", "", .name)   # eliminate indices
    if(length(.name) == 0L) .name <- "Y"
    colnames(Y) <- c(paste("1 -", .name), .name)

    message('only one variable in [0, 1] supplied - beta-distribution assumed.\ncheck this assumption.')
  }



### CHECKS
  if(!is.matrix(Y) && !is.data.frame(Y)) stop('"Y" must be either a matrix or a data.frame.')
  if(ncol(Y) <= 1) stop('"Y" must at least have two columns.')
  if((base < 1) || (base > ncol(Y))) stop('"base" must be in the range of variables.')
  if((base %% 1) != 0) stop('"base" must be an integer value.')
  if(is.na(trafo) || (is.numeric(trafo) && ((length(trafo) != 1L) || (trafo < 0) || (trafo >= .5))) ){
    stop('"trafo" must either be a logical or a (small) numeric value > 0. See ?DR_data')
  }
  if(is.null(colnames(Y))) colnames(Y) <- paste0("v", seq_len(ncol(Y)))



### NORMALIZATION - forced only if rowSums != 1 w/tolerance = norm_tol
  row.sums <- rowSums(Y)

  if( !isTRUE(all.equal( na.delete(row.sums), rep(1.0, length(na.delete(row.sums))), tolerance = norm_tol, check.attributes = FALSE)) ){
    Y <- Y/row.sums
    force.norm <- TRUE
  }



### TRANSFORMATION
  if( (is.logical(trafo) && trafo) ||
      (is.numeric(trafo) && any(na.delete(Y) < trafo, na.rm=TRUE) || any(na.delete(Y) > (1-trafo), na.rm=TRUE)) ){
    n.obs <- length(na.delete(row.sums))
    Y     <- (Y * (n.obs - 1.0) + 1.0/ncol(Y)) / n.obs
    state.tran <- TRUE
    if(is.logical(trafo) && trafo){
      force.tran <- FALSE
    } else {
      force.tran <- TRUE
    }
  }

  if(any(na.delete(Y) <= 0.0, na.rm=TRUE) || any(na.delete(Y) >= 1.0, na.rm=TRUE)){
    stop('"trafo" was suppressed, yet values on the boundary of the support are present (0 and 1).\nConsider setting "trafo" = TRUE or to a threshold.\nSee ?DR_data')
  }



### OBJECT DEFINITION
  res <- structure(
    as.matrix(Y),
    "Y.original"  = as.data.frame(Y.original),
    "dims"        = ncol(Y),
    "dim.names"   = colnames(Y),
    "obs"         = nrow(Y),
    "valid_obs"   = length(na.delete(row.sums)),
    "normalized"  = force.norm,
    "transformed" = state.tran,
    "base"        = base,
    "class" = "DirichletRegData"
  )



### warnings
  if(force.norm && force.tran){
    warning("not all rows sum up to 1 => normalization forced\n  some entries are 0 or 1 => transformation forced", immediate.=TRUE)
  } else if(force.norm){
    warning("not all rows sum up to 1 => normalization forced", immediate.=TRUE)
  } else if(force.tran){
    warning("some entries are 0 or 1 => transformation forced", immediate.=TRUE)
  }



  return(res)

}
