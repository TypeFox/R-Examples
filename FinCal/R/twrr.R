#' Computing TWRR, the time-weighted rate of return
#' 
#' @param ev ordered ending value list 
#' @param bv ordered beginning value list
#' @param cfr ordered cash flow received list
#' @seealso \code{\link{hpr}}
#' @export
#' @examples
#' twrr(ev=c(120,260),bv=c(100,240),cfr=c(2,4))
twrr <- function(ev,bv,cfr){
  r=length(ev)
  s=length(bv)
  t=length(cfr)
  wr=1
  if(r != s || r != t || s != t){
    stop("Different number of values!")
  }else{
    for(i in 1:r){
      wr = wr * (hpr(ev[i],bv[i],cfr[i]) + 1)
    }
    return(wr^(1/r) - 1)
  }
}

