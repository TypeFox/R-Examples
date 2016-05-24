spanrow <-
function(col1){
  b3=c(which(!is.na(col1) & col1!=""),length(col1)+1)
  b4=b3[-1]-b3[-length(b3)]-1
  rs=rep(0,length(col1))
  rs[b3[-length(b3)]]=b4+1
  return(rs)
}
