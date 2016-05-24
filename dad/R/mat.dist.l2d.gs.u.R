mat.dist.l2d.gs.u <-
function(x)  
{
  lot = x[, ncol(x)]
  x = x[, -ncol(x)]
  distances = diag(0, nrow = nlevels(lot))
  dimnames(distances) = list(levels(lot), levels(lot))
  for (i in 2:nlevels(lot))  
     {for (j in 1:(i-1))  
       {i.lot = which(lot == levels(lot)[i])
       j.lot = which(lot == levels(lot)[j])
       distances[i, j] = distances[j, i] = dist.l2d.gs.u(x[i.lot], x[j.lot])
       }
     }
  as.dist(distances)
}
