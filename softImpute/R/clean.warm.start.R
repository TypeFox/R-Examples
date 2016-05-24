clean.warm.start=function(a){
  if(is.null(a))return(NULL)
  d=a$d
  if(is.null(d))return(NULL)
  if(any(d>0)){
    if(length(d)==1){
      a$u=matrix(a$u,ncol=1)
      a$v=matrix(a$v,ncol=1)
    }
    a
  }
  else NULL
}

            
