make.col.vec <-
function(x, colors){
#input a vector (x) of 0, 1, 2 outputs a vector of "col1", "col2" , "col3"
# as defined by the vector of three colors (col)
result<-NULL
for( i in 1:length(x)){
  result = c(result, colors[x[i]+1])
  
}
return(result)
}

