odds_ratios <- function(x){
  
  this_call <- match.call()
  
  if(is.table(x)){
    if(!all(dim(x) == c(2L, 2L))) stop('if "x" is supplied as a table, it must be a 2 x 2 table (i.e., 2 variables with 2 categories each).')
    tt <- list(x)
  } else if(!is.data.frame(x) || (ncol(x) < 2L)){
    stop('"x" must be either a 2x2 table or a data.frame with at least 2 variables.')
  } else if(any(unlist(lapply(tt <- combn(seq_len(ncol(x)), 2L, function(i){
                                table(x[, i])
                              }, simplify=FALSE), function(tbl){
                                !all(dim(tbl) == c(2L, 2L))
                              })))){
     stop('all variables in the data frame "x" must have two categories.')
  }
  
  res <- list("call"=this_call, "x"=x, "tables"=tt, "comps"=combn(colnames(x), 2, c, simplify=FALSE))
  
  no_zero_warning <- TRUE
  
  res$ORs <- lapply(tt, function(tbl){
    if(any(tbl == 0) && no_zero_warning){
      warning('one or more frequencies equal to 0 encountered. there will be no results for at least one table.', call.=FALSE, immediate.=TRUE)
      no_zero_warning <- FALSE
      return(list("or"=NA, "lor"=NA, "se"=NA, "z.value"=NA, "p.value"=NA))
    }
    or <- (tbl[1L, 1L] * tbl[2L, 2L]) / (tbl[1L, 2L] * tbl[2L, 1L])
    lor <- log(or)
    se <- sqrt(sum(1/(tbl)))
    z <- lor/se
    p <- 2*pnorm(-abs(z))
    return(list("or"=or, "lor"=lor, "se"=se, "z.value"=z, "p.value"=p))
  })
  
  class(res) <- "REdaS_ORs"
  
  return(res)
  
}
