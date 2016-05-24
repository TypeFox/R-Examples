logit.one <-
function(u){

 if(is.finite(u)){
 if(0<u & u<1){
  out <- log(u)-log(1-u)
 }
 else if( u ==0){
  u   <- u +.0001
  out <- log(u)-log(1-u)
  }
  else if(u==1){
  u   <- u -.0001
  out <- log(u)-log(1-u)
  }else {
  warning("error..negative value passed to logit function")
  out <- NA
  }
  }else{

     warning("error..NaN value passed to logit function")
  out <- NA

  }

 out

}
