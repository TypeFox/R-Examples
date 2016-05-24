ancestors <-
function(x, ID, sources) {
  transchain <- x
  cur <- x
  while(1) {
    cur <- sources[which(ID==cur)]
    transchain <- c(transchain, cur)
    if (cur == 0) {
      break
    }
  }
  return(as.numeric(transchain))
}
