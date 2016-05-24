mdiffIRT=function(y,ai,vi,ter,sd2A,sd2V,A,W,model,eps){
  res=.C("mdiffQQ",as.double(y), as.integer(length(y)), as.double(ai),as.double(vi), as.double(sd2A),
  as.double(sd2V),as.double(ter),as.integer(model), as.double(A),as.double(W),as.integer(length(A)), as.double(eps),as.double(matrix(-999,length(y))))
  return(res[[13]])
}

