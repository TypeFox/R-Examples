predict.yai <- 
function(object,newdata,...)
{
  if (missing(newdata)) impute.yai(object,...) else 
  { 
    al <- list(...)
    impute.yai(newtargets(object,newdata,al$k,al$ann),...)
  }
}

