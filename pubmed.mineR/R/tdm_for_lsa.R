tdm_for_lsa =
 function (object, y) 
 {
     TDM_matrix_new = NULL
     testa2 = NULL
    for (i in 1:length(y)) {
         tempO = gregexpr(y[i], object@Abstract, fixed=T)
         tempa = unlist(lapply(tempO, function(x){ if (x[1] != -1) return(length(x)) else return(0)}   ))
        testa2 = matrix(tempa, nrow = 1, ncol = length(object@Abstract))
        rownames(testa2) <- y[i]
         TDM_matrix_new = rbind(TDM_matrix_new, testa2)
     }
   return(TDM_matrix_new)
 }


