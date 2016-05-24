
validate.setup <- function(setup){
  
  if(is.character(setup)){ # path of setup file
    ld <- try(load(setup), silent = TRUE)
    if(error.try(ld)){
      msg <- paste('Cannot load', setup)
      stop(msg)
    }
  }
  
  if(!is.list(setup)){
    msg <- 'setup should be a list'
    stop(msg)
  }
  
  if(!all(c('norm.stat', 'options') %in% names(setup))){
    msg <- 'norm.stat or options is not found in setup'
    stop(msg)
  }
  
  
  for(i in 1:length(setup$norm.stat$V)){
    rs <- sort(names(setup$norm.stat$score0[[i]]))
    if(length(rs) == 1){
      # the chromosome contains only one SNP.
      # in that case, setup$norm.stat$V[[i]] might not be a matrix if version between v0.8.10 and v0.8.14 was used
      # we need to reformat it as a 1x1 matrix, and also add row/column names to it
      setup$norm.stat$V[[i]] <- matrix(setup$norm.stat$V[[i]], dimnames = list(rs, rs))
    }else{
      setup$norm.stat$score0[[i]] <- setup$norm.stat$score0[[i]][rs]
      setup$norm.stat$V[[i]] <- setup$norm.stat$V[[i]][rs, rs, drop = FALSE]
    }
    
  }
  
  setup$options$only.setup <- NULL
  setup$options$save.setup <- NULL
  setup$options$path.setup <- NULL
  
  setup
  
}
