tchisq <-
function (obs) 
{esp = matrix(0, 2, dim(obs)[2])
 for (i in 1:2) {
     for (j in 1:dim(obs)[2]) {
         esp[i, j] = (sum(obs[i, ]) * sum(obs[, j]))/sum(obs)
         if (esp[i, j] == 0) {
             esp[i, j] = 0.1
         }
     }
 }
 test = sum(((obs - esp)^2)/esp)    
}
