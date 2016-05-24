get_original_term=function(m,n){

lapply(n, function(y){resAA = y$Words; resBB = lapply(resAA, function(x){ resCC = searchabsL(m, include = x); resDD = gregexpr(x, resCC@Abstract[1], ignore.case = TRUE); resEE = attr(resDD[[1]], "match.length"); resFF = substr(resCC@Abstract[1], start =  resDD[[1]][1] , stop =  resDD[[1]][1] + (resEE[1] - 1)) ; return(resFF)});return(unlist(resBB))})}


