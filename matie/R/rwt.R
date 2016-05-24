rwt <-
function(d){
  
  tRank <- function(x){
    # n<-length(x)
    rt <- rank(x, ties.method= "random")
    return(rt)
  }
  
  # rd <- apply(t(d),1,tRank)
  rd <- apply(d,2,tRank)
  return(rd)
}
