mat.ip.l2d.gp <-
function(meanL, varL)
{
  W = diag(0, nrow = length(meanL))
  dimnames(W) = list(names(meanL), names(meanL))
  W[1, 1] = l2d.gp(meanL[[1]],varL[[1]],meanL[[1]],varL[[1]])
  for (i in 2:length(meanL))  
    {W[i, i] = l2d.gp(meanL[[i]],varL[[i]],meanL[[i]],varL[[i]])
     for (j in 1:(i-1))  
      {W[i, j] = W[j, i] = l2d.gp(meanL[[i]],varL[[i]],meanL[[j]],varL[[j]])
      }
    }
  W
}
