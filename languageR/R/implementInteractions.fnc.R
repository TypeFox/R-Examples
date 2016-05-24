`implementInteractions.fnc` <-
function(m) {
   nams =  strsplit(colnames(m), ":")
   for (i in 1:length(nams)) {
     if (length(nams[[i]]) > 1) {
       m[,i] = m[,nams[[i]][1]]
       for (j in 2:length(nams[[i]])) {
         m[,i] = m[,i]*m[,nams[[i]][j]]
       }
     }
   }
   return(m)
}

