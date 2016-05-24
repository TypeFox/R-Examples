functional <-
function(type_kernel, y1,x1,bandwidth,pvalue)
   {
  a<-kde(type_kernel=type_kernel, vec_data=x1, y=y1, bw=bandwidth)$Estimated_values
  b<-pvalue
  return(a-b)
  }
