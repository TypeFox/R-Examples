distrW <- function(d, pert, np, frec, frecacum){
################# Bootstrap distribution of W under H0 #####################
# Input:
# d: Distance matrix
# pert: integer vector indicating the group each individual belongs to. 
# np: sample size for the bootstrap distribution
# Output:
# TW: null distribution of W
############################################################################

n<-dim(d)[1] # there are n elements
k<-max(pert) # there are k populations


TW <- matrix(0, np,1)+999

indor<-order(pert,c(1:n)) # <- Names of INDIVIDUALS ordered by pert
  
    for (l in 1:np){
        #choose one ind. known to be from one initial group, and calculate its W 
        #indatestar=round(runif(1,0.5,n+0.5))
        indatestar <- sample(1:n, 1)
        
        indelegidos <- matrix(0, 1,n)
         # Obtain sample with size n_pob (frec[pob]) from the original groups 
         # They are in the vector "indelegidos"
        for ( pob in 1:k){
          indelegidos[1,(frecacum[pob]+1):frecacum[pob+1]] <- sample((frecacum[pob]+1):frecacum[pob+1], frec[pob], replace=TRUE)
        }


        
        # Construct the new distance matrix
        dboot <- d[indelegidos, indelegidos]
       

        # Calculate geom. var. and deltas 
        # pert keeps the cluster each generated ind. belongs to.
        vg <- vgeo(dboot, pert)
        delta<-deltas_simple(dboot,vg, pert)

        #calculate distances from indatestar to clusters
        dx0b <- d[indatestar, indelegidos]
    
        phib <- proxi_simple(dx0b,vg, pert)
        
        TW[l] <- estW_simple(vg,delta,phib)$Wvalue
        
    } # for l in 1:np
    
  
return(TW)
}
