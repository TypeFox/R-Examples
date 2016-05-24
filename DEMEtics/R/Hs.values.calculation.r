# Function used within H.out.r
Hs.values.calculation <- function(Hj.values2.onelocus){
Hs.one.locus<-Hs(as.numeric(as.vector(Hj.values2.onelocus$Hj.value)))
        
# The Hs value is calculated for the actual locus from the Hj values
# of all the populations for the actual locus.
          
Hs.values.onelocus<-cbind(Hs.one.locus,as.character((Hj.values2.onelocus$Locus)[1]))
invisible(Hs.values.onelocus)
          
}
