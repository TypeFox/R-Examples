
#-------------------------------------------------------------------------------#
#------------------------- SOME USEFUL DATA CHECK FUNCTIONS --------------------#

check.rep.num <-function(rep.num)
{
  res = as.integer(10) 
  if( is.vector(rep.num) ) 
  { 
    rep.num = as.integer(rep.num[1])
    if( !is.na(rep.num) && !is.nan(rep.num) )
      if( 0 < rep.num && rep.num < 100 )  res = rep.num
  }
  return(res)
}

check.subset.ratio <- function(data, subset.ratio)
{
  if( !is.numeric(subset.ratio) ) subset.ratio=0.75
  if( subset.ratio > 1 || subset.ratio <= 0 ) subset.ratio=0.75

  obj.num = dim(data)[1]
  if( obj.num*subset.ratio < 1 ) subset.ratio = 1/obj.num
  return(subset.ratio)
}

cut.cl.num <-function(data, cl.num, ratio)
{
  max.cls.num = dim(data)[1]*ratio
  if( max(cl.num) > max.cls.num ) cl.num = cl.num[cl.num <= max.cls.num]
  
  if( length(cl.num) < 1 )
  {
    cat("The tripple <data, cls.num, subset.ratio> is incorrect.\n")
    cat("There are three solutions:\n")
    cat("1.increase number of objects in data,\n")
    cat("2.increase subset.ratio,\n")
    cat("3.decrease values in cl.num vector.\n")
    stop("There is not enough data to be clustered.")
  } 
  return(cl.num)
}

cut.cl.num.pred <-function(data, cl.num, ratio)
{
  max.cls.num = dim(data)[1] * min(ratio, 1 - ratio) 
  if( max(cl.num) > max.cls.num ) cl.num = cl.num[cl.num <= max.cls.num]
  
  if( length(cl.num) < 1 )
  {
    cat("The tripple <data, cls.num, subset.ratio> is incorrect.\n")
    cat("There are three solutions:\n")
    cat("1.increase number of objects in data,\n")
    cat("2.change subset.ratio\n")
    cat("3.decrease values in cl.num vector.\n")
    stop("Subsets produced during execution are too small to be clustered.")
  }

  return(cl.num)
}
