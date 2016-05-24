mat_split<- function(mat,nchunk){
  mat.story=list()
  n=nrow(mat)
  #case when input is number of chunks 
  if (length(nchunk) == 1) {
    nc=floor(n/nchunk)
    seq.chu=seq(from=1,to=n,by=nc)
    seq.chu[nchunk+1] = n
  }
  #case when input is specific block sizes
  if (length(nchunk) > 1){
  #  nc = length(nchunk)
    seq.chu=0
    seq.chu=c(seq.chu,cumsum(nchunk))+1
    #fix last value
    seq.chu[length(seq.chu)] = n
  }
  
  kk=length(seq.chu)
  
  if (kk > 2)
  {
    mat.story[[1]]=mat[1:seq.chu[2],]
  
    for(k in 2:(kk-1)){
      a=(seq.chu[k]+1)
      b=(seq.chu[k+1])
      mat.story[[k]]=mat[a : b, ]
      }
    
  } else { #case of a single chunk
    mat.story[[1]] = mat  
  }
  out=list()
  out$splitMat=mat.story
  
  out
}


#     mat.story[[2]]=mat2[92 : 181, ]
#print(mat.story)
#