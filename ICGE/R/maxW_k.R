maxW_k <- function(x_d,v_d, pert){
######################### Used in INCAnumclu #################################
#                  When number of clusters k>=2
# Input: 
# x_d <- distance matrix of the partition considered as populations
# v_d <- distance from each individual of the testing cluster to the rest.
# pert <- partition of individual considered as populations.
# Output:
# contar: number of individuals of the testing cluster considered as atipical or #         well classified.
##########################################################################

if (is.null(dim(v_d))){dim(v_d) <- c(1, length(v_d))}
nt <- dim(v_d)[1]
k <- max(pert)
if (is.null(dim(x_d))){dim(x_d) <- c(1, length(x_d))}
nn <- dim(x_d)[1]
atipico <- matrix(0,nt,1)




vg <- vgeo(x_d, pert);


delta <- deltas_simple(x_d,vg,pert)


# Calculo de W para los individuos de xx/d_x

TW <- matrix(0,nn,1)+999;
for (l in 1:nn){   
   
    phi <- proxi_simple(x_d[l,],vg,pert)
    
    TW[l] <- estW_simple(vg,delta, phi)$Wvalue
    
}


# Calculate  maxumum of  W for ind. in  xx
M <- max(TW)

############################################
# Calculate W for ind. in  v_d and evalutate whether is atypical or not

for (ind in 1:nt){
    phi <- proxi_simple(v_d[ind,],vg,pert)
    W0 <- estW_simple(vg,delta,phi)$Wvalue
    
    # Clasificacion en atipico o no atipico
    
    if (W0<=M){
        atipico[ind] <-0
    }
    if (W0 > M){
        atipico[ind] <-1
    }
    
} #for ind=1:nt


contar <- 0 # To count atypical units in v_d
contar <- sum(atipico)
  

#out <- list(Wx=TW, Ati=atipico, suma=contar)
out<-contar

return(out)

} # end function







