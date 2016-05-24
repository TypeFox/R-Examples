#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified (unused)


ehh <- function(hapldata,n1,n2){

if(missing(n1) & missing(n2)){
   n1 <- floor(dim(hapldata)[2]/2)
   n2 <- n1
}

core   <- hapldata[,n1:n2,drop=FALSE] 
values <- my_unique(core)

# return(values)

n        <- values$numHap # number of haplotypes
nsites   <- dim(hapldata)[2] # number of sites

data_left  <- hapldata[, 1:(n1-1),drop=FALSE]
data_right <- hapldata[, (n2+1):nsites,drop=FALSE]

m2 <- nsites - (n2-n1)  # size of core (sites)

h              <- matrix(1,nrow=n,ncol=m2)
rownames(h)    <- rownames(values$uniquematrix)
hvar           <- matrix(0,nrow=n,ncol=m2)
rownames(hvar) <- rownames(values$uniquematrix)
core_p         <- rep(0,n)
names(core_p)  <- rownames(values$uniquematrix)

for(k in 1:n){

  idx        <- values$idx==k # Welche sind Haplotyp k ?
  core_p[k]  <- mean(idx)     # Frequenz des Haplotypen k ---> (Anzahl k)/(Alle Sequenzen)
  
  dl <- data_left [idx,,drop=FALSE]      # 
  dr <- data_right[idx,,drop=FALSE]      #
  
  # Check the left side of the core
  for(i in 1:(n1-1)){
      dlfenster <- dl[,i:dim(dl)[2],drop=FALSE]      
     # print(dlfenster)
      val <- my_unique(dlfenster) # von links nach rechts --> verkleiner die Fenster
      z   <- sum(val$sizHap)
      if(z>1){
         p <- val$sizHap/z # relative Frequency of Haplotypes
         a <- sum(p^2)
         h[k,i]    <- (a-1/z)/(1-1/z) #
         hvar[k,i] <- (2*(z-1)/z^3)*(2*(z-2)*(sum(p^3)-a^2) + a - a^2 )
      }   
  }  
  # Check the ride side of the core 
  for(j in (n2+1):m2){
      val <- my_unique(dr[,1:(j-n1),drop=FALSE]) # 
      z <- sum(val$sizHap)
      if(z>1){
        p <- val$sizHap/z
        a <- sum(p^2)
        h[k,j]    <- (a-1/z)/(1-1/z)
        hvar[k,j] <-(2*(z-1)/z^3)*(2*(z-2)*(sum(p^3)-a^2) + a - a^2 )
      }
  }

} # End of for 1:n


return(list(h=h,hvar=hvar,core_p=core_p))

}






