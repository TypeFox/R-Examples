knots.order <- function(penden.env) {
  help.seq <- seq(1,length(get("knots",penden.env)))
  val <- ceiling(median(help.seq))
  if(get("q",penden.env)==2) x.order <- c(0,help.seq[val]/(2^get("d",penden.env)+1),1) else x.order <- c(0,1)
  ind <- 1
  #ungerade <- mod(help.seq,2)
  ungerade <- help.seq[which((help.seq-floor(help.seq/2)*2)!=0)]
  if(get("q",penden.env)==2) ungerade[ungerade>val] <- ungerade[ungerade>val]+1
  for(i in 0:(get("d",penden.env)-1)) {
    help.seq <- seq(1,((2)^(i))) #2^(i+ind)
    if(get("q",penden.env)==1) x.order <- c(x.order,ungerade[help.seq]/((2)^(i+1))) #(2+ind)^
    if(get("q",penden.env)==2) x.order <- c(x.order,ungerade[help.seq]/((2)^(i+1)+1)) #(2+ind)^
  }
  return(x.order*2^get("d",penden.env)+1)#(2+ind)^
}
