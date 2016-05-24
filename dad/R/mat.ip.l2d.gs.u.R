mat.ip.l2d.gs.u <-
function(x)  
{
  lot = as.factor(x[, ncol(x)])
  x = x[, -ncol(x)]
  meanL<-by(x, lot, mean);
  varL<-by(x, lot, var);
  W = diag(0, nrow = length(meanL))
  dimnames(W) = list(names(meanL), names(meanL))
  W[1, 1] = l2d.gp.u(meanL[[1]],varL[[1]],meanL[[1]],varL[[1]])
  for (i in 2:length(meanL))
    {W[i, i] = l2d.gp.u(meanL[[i]],varL[[i]],meanL[[i]],varL[[i]])
    for (j in 1:(i-1))
      {W[i,j] = W[j, i] = l2d.gp.u(meanL[[i]],varL[[i]],meanL[[j]],varL[[j]])
      }
    }
  W
}
