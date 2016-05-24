# ' compare signs of the coefficients in two vectors
# ' @param trueA first vector
# ' @param Aalgo second vector
compare_sign<-function(trueA=trueA,Aalgo=Aalgo){
      quivrai0=which(trueA==0)
      nbbon0=length(which(Aalgo[quivrai0]==0))
      nbbon1=length(which(Aalgo[-quivrai0]!=0))
      nbfaux0=length(which(Aalgo[-quivrai0]==0))
      nb0mank=length(quivrai0)-nbbon0
      
      #comptage des signes +
      quivraiplus=which(trueA>0)
      quivraimoins=which(trueA<0)
      #attention , on n'a pas une partition (a cause des 0)
      nbbonplus=length(which(Aalgo[quivraiplus]>0))
      nbbonmoins=length(which(Aalgo[quivraimoins]<0))
      nbfauxplus=length(which(Aalgo[quivraimoins]>0))
      nbfauxmoins=length(which(Aalgo[quivraiplus]<0))
      
      return(list(truepositive=nbbonplus,truenegative=nbbonmoins,falsepositive=nbfauxplus,falsenegative=nbfauxmoins))    
      
   
}