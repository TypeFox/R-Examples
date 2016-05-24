quick.norm <-
function(A,mod=2){
  if (mod==1){
    t(apply(A,1,function(i){
      if (sum(i)==0) {return(i)}
      else {i/sqrt(sum(i^2))}}))
  } else if (mod==2){
    t(apply(A,1,function(i){
      if (sum(i)==0) {return(i)}
      else {i/sum(i)}}))
  }
}

