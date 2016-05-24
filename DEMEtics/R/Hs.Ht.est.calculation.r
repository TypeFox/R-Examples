# Function called within H.out.r

Hs.Ht.est.calculation <- function(Hs.values2.onelocus,Ht.values2.onelocus,sample.sizes2.onelocus){
# The harmonic mean is calculated from the sample sizes...
n <- length(split(sample.sizes2.onelocus,sample.sizes2.onelocus$population))
sample.size <- as.numeric(as.vector(sample.sizes2.onelocus$sample.size))
N<-harmonic(sample.size,n)
        
# ... and ascribed to the variable N.

                      
Hs.actual<-as.numeric(as.vector(Hs.values2.onelocus$Hs.value))
Ht.actual<-as.numeric(as.vector(Ht.values2.onelocus$Ht.value))
                                          
Hs.est.actual <- Hs.est(Hs.actual,N)
Ht.est.actual <- Ht.est(Ht.actual,Hs.est.actual,N,n)
                                          
Hs.est.values<-cbind(Hs.est.actual,as.character((sample.sizes2.onelocus$locus)[1]))
Ht.est.values<-cbind(Ht.est.actual,as.character((sample.sizes2.onelocus$locus)[1]))

                                        # The Ht.est and Hs.est values
                                                    # for each locus
                                                    # are combined and
                                                    # the locus names
                                                    # are added to the
                                                    # values they
                                                    # belong to
H.est.values <- list(Hs.est.values,Ht.est.values)
invisible(H.est.values)
}
