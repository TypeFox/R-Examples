### crm function for finding next dose level and the updated parameter 'a'
crm <- function(target,prior,ptdata,model=1,a0=1,b=3){
  if(model != 1 && model != 2){
    stop("Error: model must be 1 or 2\n")
  }
  if(target<0 || target > 1){
    stop("Error: target must be greater than 0 and less than 1\n")
  }
  if(any(prior<0) || any(prior>1)){
    stop("Error: All elements in prior must be greater than 0 and less than 1\n")
  }
  ptdt = matrix(ptdata,ncol=2)
  callCRM = .C("CRM",as.integer(model),as.double(target),as.double(prior),as.integer(length(prior)),
    as.double(a0),as.double(b),as.integer(ptdt),as.integer(nrow(ptdt)),mtd=as.integer(0),
    aMean=as.double(0),PACKAGE="CRM")
  return(list(MTD=callCRM$mtd,a=callCRM$aMean))
}

