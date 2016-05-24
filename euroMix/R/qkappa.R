qkappa=function(kappa=c(0,1,0),q=NULL)
kappa[1]*q[[1]]+kappa[2]*q[[2]]+kappa[3]*diag(q[[3]])