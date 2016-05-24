#' random.range - supporting function for other rangemodel fucntions; no
#'                practical use
#' @description random.range is used within other rangemodel functions to
#'              radomly place given number of species occurences
#' @param uid a vector of unique ids for selection
#' @param nb a neighbour object similar to generated from 'shp2nb'
#' @param range.size a vector of number of sites occupied by each species
#' @param var an optional vector of variables for constraining the randomization
#' @param first If true, var is used while choosing the first occurence as well.
#'        if var is null, first is always set FALSE
#' @details this function is not intended for any direct use but is called
#'           within other functions of 'rangemodelR'.
#' @return a numeric vector specifying selected possitions in 'uid'
#' @export

random.range <- function(uid,nb,range.size,var,first){
  suppressWarnings(if(is.na(nb)){
    sel.vec <- sample(uid,range.size,prob = var)
  }else{
    sel.vec <- NULL
    sel.nb <- NULL # objects to store seleted cells and neighbours

    if(first == T){
      sel.vec <- c(sel.vec,sample(uid,1,prob = var)) # select first cell

    }else{sel.vec <- c(sel.vec,sample(uid,1))} # select first cell)

    if(range.size == 1){
      sel.vec
      }else{
        sel.nb <- nb[which(uid%in%sel.vec[1])] # query the nb object
        # with the possitionat which selected cell appears in the uid vector
        sel.nb.vec <- unlist(sel.nb)
        if(is.null(var)){
          for(i in 2:range.size){
            sel.vec[i] <- sample(uid[unique(sel.nb.vec)],1)
            sel.nb <- nb[which(uid%in%sel.vec)]
            sel.nb.vec <- unlist(sel.nb)
            sel.nb.vec <- sel.nb.vec[!sel.nb.vec%in%which(uid%in%sel.vec)]
          }  #extend the selction to desired length
        }else{
          for(i in 2:range.size){
            sel.vec[i] <- sample(uid[unique(sel.nb.vec)],1,
                                 prob = var[unique(sel.nb.vec)] )
            sel.nb <- nb[which(uid%in%sel.vec)]
            sel.nb.vec <- unlist(sel.nb)
            sel.nb.vec <- sel.nb.vec[!sel.nb.vec%in%which(uid%in%sel.vec)]
          }  #extend the selction to desired length
        }
      }
      }
  )
  return(sel.vec)
}
