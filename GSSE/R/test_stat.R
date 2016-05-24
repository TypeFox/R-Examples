
#' @export

## function to compute test statistic 

test_stat = function(F1, F2){
  Fest = rbind(F1, F2);
  difference <- matrix(apply(Fest,2,diff), ncol=1);
  if(dim(Fest)[2]==1){
    ## single point
    out <- apply(abs(difference),1,max)
  } else {
    out <- apply(abs(difference),2,max)
  }
  return(out)
}

