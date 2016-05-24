#' Function to do pairwise deletion of missing rows. 
#' 
#' The aim in pair-wise deletions is to retain the largest 
#' number of available data pairs with all non-missing data.
#' 
#' @param x vector of x data
#' @param y vector of y data
#' @return 
#' \item{newx}{new vector x after removing pairwise missing data} 
#' \item{newy}{new vector y after removing pairwise missing data} 
## @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @examples
#' 
#' \dontrun{
#' x=sample(1:10);y=sample(1:10);x[2]=NA; y[3]=NA
#' napair(x,y)}
#' 
#' @export
 
napair <-
function(x,y){
#author:  H D Vinod, Fordham University, 2013
ava.x=which(!is.na(x))#ava means available
ava.y=which(!is.na(y))#ava means non-missing
ava.both=intersect(ava.x,ava.y)
list(newx=x[ava.both],#delete NAs from x
newy=y[ava.both])#delete NAs from y
}
