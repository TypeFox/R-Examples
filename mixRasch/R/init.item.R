`init.item` <-
function(it.stats){

 dimset <- dim(it.stats$S.ih)
 delta.i <- array(0, dim=dimset[2],dimnames=list(colnames(it.stats$S.ih)))  
 tau     <- array(0, dim=dimset,dimnames=list(NULL,colnames(it.stats$S.ih)))
 tau[is.na(it.stats$S.ih)] <- NA
 delta   <- t(apply(tau,1, function(XXX) XXX + delta.i))
 list(delta.i = delta.i, tau = tau, delta = delta)
} # end init.item

