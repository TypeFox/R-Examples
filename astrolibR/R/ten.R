ten=function(dd,mm=0,ss=0) {
  np = nargs()

  testneg = NULL
  sign = 1
  if((np==1 && is.character(dd)) ){
      temp = gsub(':',' ',dd)
      testneg = grep('-',temp)
      if(length(testneg)==1) temp = sub('-','',temp)
      if(length(testneg)==1) sign = -1
      vector = as.double(strsplit(temp,' ')[[1]])
    }                                                     #  colon-separated input

  else{
      vector= as.double(c(dd,mm,ss))  # comma-separated input
}

  fac=c(1,60,3600)
  vector = abs(vector)
  
  return(sign*sum(vector/fac))
}
