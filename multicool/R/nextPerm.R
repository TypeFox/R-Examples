nextPerm = function(mcObj){
  r = mcObj$mc$nextPerm()
  
  if(r$b){
    return(mcObj$elements[r$set])
  }else{
    return(FALSE)
  }
}
