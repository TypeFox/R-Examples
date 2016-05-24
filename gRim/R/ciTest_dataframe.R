
## ########################################################
##
## Testing for condional independence in a dataframe
## <x>  : dataframe
## <set>: NULL, a vector or a formula
##
## ########################################################

ciTest_df <- function(x, set=NULL,...){
  if (is.null(set)){
    set <- names(x)
  } else {
    if (inherits(set,c("formula","character"))){
      set <- unlist(rhsFormula2list(set))
      set <- names(x)[pmatch(set, names(x))]
    }
  }
  
  wdata       <- x[,set]
  varTypes    <- uniquePrim(unlist(lapply(wdata, class)))

  has.factor  <- "factor" %in% varTypes
  has.numeric <- any(c("integer","numeric") %in% varTypes)

  if (has.factor & has.numeric){
    .ciTest_df_internal(wdata, set,...)
  } else {
    if (has.factor){
      ciTest_table(xtabs(~., data=wdata),set=set,...)
    } else {
      if (has.numeric){
        ciTest_mvn(cov.wt(wdata,method="ML"), set=set,...)
      } else {
        stop("Strange error...\n")
      }
    } 
  }
}

## ########################################################
##
## Testing for condional independence in mixed data
## <x>  : dataframe
## <set>: NULL, a vector or a formula
##
## ########################################################

.ciTest_df_internal <- function(x,set=NULL,...){
  ##cat("CHK: ciTestmixed\n")

  obj <- mmod(list(set), data=x)

  ans <- testdelete(obj, set[1:2])
  ans <- ans[c(1,3,2)] ## FIXME: This is fragile
  ans$method   <- "CHISQ"
  ans$statname <- "DEV"
  ans$varNames <- set
  class(ans)   <- "citest"
  ans
}



