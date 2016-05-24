`kin.cohort` <-
function(...,method=c("marginal","mml","chatterjee", "moments","km","watcholder")){
if (tolower(method[1]) %in% c("marginal","mml","chatterjee") )
    o<- kc.marginal(...)
else if (tolower(method[1]) %in% c("moments","km","watcholder") ) 
    o<- kc.moments(...)
else  
  stop("please provide a valid method")
o
}

