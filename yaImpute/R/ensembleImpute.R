ensembleImpute <- function (imputes,method="mean",...)
{
  cl = match.call()
  posM = c("mean","median")
  if (!(method %in% posM)) stop ('method="',method,'" must be one of ',
    paste0('"',posM,'"',collapse=", "))

  for (i in 1:length(imputes)) if (("yai" %in% class(imputes[[i]])))
    imputes[[i]] = impute.yai(imputes[[i]],...)

  colns = unique(unlist(lapply(imputes,function(x) colnames(x))))
  ctf = if (method!="mean") median else mean
  rowns = sort(unique(unlist(lapply(imputes,function (x) rownames(x)))))
  ave = list()
  sd = list()
  N = list()
  methods = list()
  for (cl in colns)
  {
    one = matrix(unlist(lapply(imputes,function (x,rowns,cl) 
          x[rowns,cl],rowns,cl)), nrow=length(rowns))
    n = apply(one,1,function (x) sum(!is.na(x)))
    if (any(n != length(imputes))) N[[cl]] = n
    if (mode(one) == "character") 
    {
      ave[[cl]] = apply(one,1,function(x) 
        {
          x = na.omit(x)
          if (length(x) == 0) return(NA)
          x = table(x)
          x = x+(runif(length(x))*.01)
          names(x)[which.max(x)]
        })
      ave[[cl]] = as.factor(ave[[cl]])
      methods[[cl]] = "mode"
    } else {
      ave[[cl]] = apply(one,1,ctf,na.rm=TRUE)
      ave[[cl]][is.nan(ave[[cl]])] = NA
      sd [[cl]] = apply(one,1,function (x) 
        {
          x = na.omit(x)
          if (length(x) > 1) sd(x) else 0
        })
      methods[[cl]] = method
    }
  }
  ans = as.data.frame(ave)
  rownames(ans) = rowns
  class(ans) = c("impute.yai","data.frame")
  if (length(sd)>0) 
  {
    sumsgtz = unlist(lapply(sd,sum)) > 0
    if (any(sumsgtz))
    {
      sd = as.data.frame(sd[sumsgtz])
      rownames(sd) = rowns
      attr(ans,"sd") = sd
    }
  } 
  attr(ans,"N") = if (length(N)>0) as.data.frame(N) else length(imputes)
  attr(ans,"methods") = unlist(methods)
  ans  
}  

