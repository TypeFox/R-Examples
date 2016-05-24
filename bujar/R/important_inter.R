#Given a pxp matrix with some interaction scores, to rank these scores
# create an index of the values in descending order
import.inter <- function(pred.data,cross.tab){
  n.preds <- ncol(pred.data)
  pred.names <- colnames(pred.data)
  search.index <- ((n.preds^2) + 1) - rank(cross.tab, ties.method = "first")

#  n.important <- max(2,round(0.1 * ((n.preds^2)/2),0))
  n.important <- n.preds
  var1.names <- rep(" ",n.important)
  var1.index <- rep(0,n.important)
  var2.names <- rep(" ",n.important)
  var2.index <- rep(0,n.important)
  int.size <- rep(0,n.important)

  for (i in 1:n.important) {

    index.match <- match(i,search.index)

    j <- min(n.preds,trunc(index.match/n.preds) + 1)
    var1.index[i] <- j
    var1.names[i] <- pred.names[j]

    k <- index.match%%n.preds
    if (k > 0) {   #only do this if k > 0 - otherwise we have all zeros from here on 
      var2.index[i] <- k
      var2.names[i] <- pred.names[k]

      int.size[i] <- cross.tab[k,j]
    }

  }

  rank.list <- data.frame(var1.index,var1.names,var2.index,var2.names,int.size)

  return(list(rank.list = rank.list, interactions = cross.tab))
}

