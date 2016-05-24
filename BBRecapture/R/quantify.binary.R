quant.binary.integer=function(x){

  if(is.character(x)){
    binarystring=x
    bh=as.numeric(unlist(strsplit(binarystring,split="")))
  }else{
    bh=x
  }
  if(length(setdiff(unique(bh),c(0,1)) )>0){
    stop('Binary history contains non binary 0/1 content')
  }

  n=length(bh)
  
  if(n==0){
    out=0
  }else{
  out=sum(bh*(2^(0:(n-1))))
}
  
  return(out)

}	

quant.binary=function(x){

  if(is.character(x)){
    binarystring=x
    bh=as.numeric(unlist(strsplit(binarystring,split="")))
  }else{
    bh=x
  }
  if(length(setdiff(unique(bh),c(0,1)) )>0){
    stop('Binary history contains non binary 0/1 content')
  }

  n=length(bh)
  if(n==0){
    out=0
  }else{
    
  num=sum(bh*(2^(0:(n-1))))
  den=sum(rep(1,n)*(2^(0:(n-1))))

  out=num/den
}
  
  return(out)

}	


quant.binary.markov=function(x,markov.ord){

  if(is.character(x)){
    binarystring=x
    bh=as.numeric(unlist(strsplit(binarystring,split="")))
  }else{
    bh=x
  }
  if(length(setdiff(unique(bh),c(0,1)) )>0){
    stop('Binary history contains non binary 0/1 content')
  }

  bh=c(rep(0,markov.ord-1),bh)
  n=length(bh)

  if(n==0){
    out=0
  }else{
    
  num=sum(bh*(2^(0:(n-1))))
  den=sum(rep(1,n)*(2^(0:(n-1))))

  out=num/den
}
  
  return(out)

}	


quant.binary.counts=function(x){

  if(is.character(x)){
    binarystring=x
    bh=as.numeric(unlist(strsplit(binarystring,split="")))
  }else{
    bh=x
  }
  if(length(setdiff(unique(bh),c(0,1)) )>0){
    stop('Binary history contains non binary 0/1 content')
  }

  n=length(bh)
  if(n==0){
    out=0
  }else{
    
  num=sum(bh)
  den=length(bh)

  out=num/den
}
  return(out)

}	


quant.binary.counts.integer=function(x){

  if(is.character(x)){
    binarystring=x
    bh=as.numeric(unlist(strsplit(binarystring,split="")))
  }else{
    bh=x
  }
  if(length(setdiff(unique(bh),c(0,1)) )>0){
    stop('Binary history contains non binary 0/1 content')
  }

  n=length(bh)
  if(n==0){
    out=0
  }else{
    
  out=sum(bh)

}
  return(out)

}	
