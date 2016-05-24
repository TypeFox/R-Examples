getCov <-
function(obj){
  if(length(obj)!=3) print("Please input the correct formula in the form of y~ x1+x2+...")
  term=obj[[3]]
  res=all.vars(term)
  return(res)
}
