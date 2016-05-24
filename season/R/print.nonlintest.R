## print.nonlintest.R
## Prints basic results from nonlintest
print.nonlintest<-function(x, ...){
  ## Check
  if (class(x) != "nonlintest"){
    stop("Object must be of class 'nonlintest'")
  } 

# 1) Stats on third order moment
diff_l=x$region*0
diff_u=x$region*0
index=which(x$region>0,arr.ind=T)
if(length(index)>0){diff_u[index]=x$region[index]}
index=which(x$region<0,arr.ind=T)
if(length(index)>0){diff_l[index]=x$region[index]}

  cat('Largest and smallest co-ordinates of the third-order moment outside the test limits\n');
# Upper limit
if(max(diff_u)>0){
   index=which(diff_u==max(diff_u),arr.ind=T);
   cat('Largest positive difference at lags:\n');
   cat(index-1,'\n');
   cat('Size of largest largest difference:\n');
   cat(max(diff_u),'\n');
}
if(max(diff_u)==0){
   cat('Largest difference is zero\n');
}
# Lower limit
if(min(diff_l)<0){
   index=which(diff_l==min(diff_l),arr.ind=T);
   cat('Largest negative difference at lags:\n');
   cat(index-1,'\n');
   cat('Size of largest negative difference:\n');
   cat(min(diff_l),'\n');
}
if(min(diff_l)==0){
   cat('Smallest difference is zero\n');
}


# 2) Bootstrap stats
  cat('\nBootstrap test of non-linearity using the third-order moment\n');
  cat('Statistics for areas outside test limits:\n');
  cat('observed     obs/SD     median-null    95%-null    p-value\n');
  cat(x$stats$outside,x$stats$stan,x$stats$median,x$stats$upper,x$stats$pvalue,'\n');
} # end of function
