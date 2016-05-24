interv_covariate <- function(n, tau, delta){
  stopifnot( #Are the arguments valid?
    n%%1==0 & n>=0,
    length(tau)==length(delta),
    tau==Inf | ( tau%%1==0 & tau>=0 & tau<=n),
    0<=delta & delta<=1
  )
  single_intervention <- function(k) ifelse((1:n)>=tau[k], delta[k]^((1:n)-tau[k]), 0) #generates covariate vector for intervention k
  result <- vapply(seq(along=tau), single_intervention, FUN.VALUE=numeric(n)) #matrix with the k intervention covariates in the columns
  if(length(tau)>0) colnames(result) <- paste("interv_", seq(along=tau), sep="")
  return(result)
}
