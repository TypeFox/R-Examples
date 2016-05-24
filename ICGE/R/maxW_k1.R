maxW_k1 <- function(x_d,v_d){
######################### Used in INCAnumclu #################################
#                  When number of clusters k=1
# Input:
# x_d <- distance matrix of the partition considered as populations
# v_d <- distance from each individual of the testing cluster to the rest.
# pert <- partition of individuals considered as populations.
# Output:
# contar: number of individuals of the testing cluster considered as atipical or #         well classified.
##########################################################################
if (is.null(dim(v_d))){dim(v_d) <- c(1, length(v_d))}
nt <- dim(v_d)[1]
k <- 1
if (is.null(dim(x_d))){dim(x_d) <- c(1, length(x_d))}
nn <- dim(x_d)[1]
atipico <- matrix(0,nt,1)


# As there is an unique cluster in population
pert <- matrix(1,nn,1)
vg <- vgeo(x_d, pert);


# Calculate W for ind. in xx/d_x

TW <- matrix(0,nn,1)+999;
for (l in 1:nn){   
   
    phib <- proxi_simple(x_d[l,],vg,pert);
    
    TW[l] <- phib              #this case W is the proximity function    
}


# Calculate the maximum of W for ind. in xx
M <- max(TW)

############################################
# Calculate W for ind. in v_d and evaluter whether is atypical or not
 
for (ind in 1:nt){
    W0 <- proxi_simple(v_d[ind,],vg,pert)
    # Classification    
    if (W0<=M){
        atipico[ind] <-0
    }
    if (W0 > M){
        atipico[ind] <-1
    } 
} #for ind=1:nt



contar <- 0 # To count atypicals in  v_d
contar <- sum(atipico)
  
out<-contar

return(out)

} # end function
