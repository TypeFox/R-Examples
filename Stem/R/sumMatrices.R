`sumMatrices` <-
function(matrices){
         if (length(matrices) > 2) matrices[[1]] + Recall(matrices[-1])
         else matrices[[1]] + matrices[[2]]
}

