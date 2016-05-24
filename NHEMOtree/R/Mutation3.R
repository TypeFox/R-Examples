Mutation3 <-
function(Chromosom, m){ 	   
  j   <- 0
  L   <- length(Chromosom)
  FLIP<- c()

  for (i in 1:L){
     # Ziehen einer Zufallszahl, die kleiner als Mutationswsk sein muss, damit geflipt wird!
       temp<- sort(runif(n=1, min=0, max=1))

       if (temp <= m){
           if (Chromosom[i]==0) Chromosom[i]<- 1 else Chromosom[i]<- 0 # Vertauschen!
           j      <- j+1
           FLIP[j]<- i
  }}
  # Ausgabe:
    return(list(Chromosom_mutiert=Chromosom, Genes_flipped=FLIP)) 
}
