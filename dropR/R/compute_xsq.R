
#' Compute Chisq test for a list 
#' 
#' Compute chisq given a list of data.frames
#' 
#' @param li a list
#' @param pos integer position within the respective data.frane
#' @param participants known integer
#' @param sel integer position in th result table
#' @export
#'  
compute_xsq <- function(li,pos,participants,sel){
  drop_outs <- unlist(lapply(li,'[',pos))
  remain <- participants - drop_outs
  
  M <- cbind(remain,drop_outs)
  dimnames(M) <- list(condition = names(li),
                      status = c('remaining','dropouts') )
  
  chisq.test(as.table(M[names(li)[sel],]),correct = T)
  
}




