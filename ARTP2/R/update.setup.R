
update.setup <- function(setup, nperm, lambda, nthread){
  
  if(!is.null(nperm)){
    setup$options$nperm <- nperm
  }
  
  if(!is.null(nthread)){
    setup$options$nthread <- nthread
  }
  
  if(is.null(lambda)){
    lambda <- 1.0
  }else{
    setup$options$lambda <- setup$options$lambda * lambda
    for(i in 1:length(setup$norm.stat$V)){
      setup$norm.stat$V[[i]] <- setup$norm.stat$V[[i]] / lambda
      setup$norm.stat$score0[[i]] <- setup$norm.stat$score0[[i]] / lambda
    }
  }
  
  tmp <- .C("check_nthread", nthread = as.integer(setup$options$nthread), PACKAGE = "ARTP2")
  setup$options$nthread <- tmp$nthread
  
  setup
  
  
}
