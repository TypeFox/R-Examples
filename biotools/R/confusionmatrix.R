confusionmatrix <-
function(obs, predict)
{
   obs <- as.vector(obs)
   predict <- as.vector(predict)
   if (length(obs) != length(predict))
      stop("incompatible dimensions!")
   aux <- data.frame(gr = obs, newgr = predict)
   gr <- NULL
   ng <- length(unique(aux[["gr"]]))
   lev <- levels(as.factor(aux[["gr"]]))
   newlab <- paste("new", lev)
   m <- matrix(0, ng, ng,
      dimnames = list(lev, newlab))
   equal <- function(x, y) sum(x == y)
   for(i in 1:ng) {
      for(j in 1:ng) {
         m[i, j] <-
         equal(subset(aux, gr == lev[i])[["newgr"]], lev[j])
      }
   }
   return(m)
}
