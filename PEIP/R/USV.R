USV <-
function(G)
{

uu = svd(G)

##  assign("U", uu$u, envir = .GlobalEnv)
## assign("S", uu$d, envir = .GlobalEnv)
## assign("V", uu$v, envir = .GlobalEnv)

return(list(U=uu$u, S=uu$d, V= uu$v))

}
