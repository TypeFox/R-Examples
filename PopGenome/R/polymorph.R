get.polymorph <- function(matr){

#count1 <<- 0
poly  <- apply(matr,2,function(x){

  # cat(count1,"\n")
  # count1 <<- count1 + 1 
  return(length(unique(x))!=1)

})

return(which(poly))

}
