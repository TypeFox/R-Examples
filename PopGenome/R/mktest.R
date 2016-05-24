#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified

########################################################
###################### MKTEST ##########################
########################################################
# unused 

mktest <- function(seq,taxset1,taxset2){

# In der Funktion DsPsDnPn wird gesplittet !

O <- DsPsDnPn(seq,taxset1,taxset2)
Ds <- O$Ds
Ps <- O$Ps
Dn <- O$Dn
Pn <- O$Pn



#p <- fisherextest(round(Ds),round(Ps),round(Dn),round(Pn))



#cat("Ds - Syn. Divergences; Ps - Syn. Polymorphisms","\n")
#cat("Dn - Nonsyn. Divergences; Pn - Nonsyn. Polymorphisms","\n")
#cat("\n")
#cat("Fisher''s exact test","\n","P-value (two tailed):",p,"-->","significant ? :",sigtag(p))
#cat("\n")

X <- gtest(Ds,Ps,Dn,Pn)
P <- X$P
G <- X$G

#cat("G TEST, G value:",G,"\n")
#cat("P-value:",P)
 #if(P!="NaN"){
  # cat("--> significant ? :",sigtag(P))
 #}
#cat("\n")

 if(Ps!=0 && Dn!=0 && Ds!=0){

   NI <- (Pn/Dn)/(Ps/Ds)
   alphax <- 1-(Pn*Ds)/(Ps*Dn)
 }else{
   NI <- 0
   alphax <- "NaN"
 }

#cat("Neutrality Index (NI), [(Pn/Dn)/(Ps/Ds)]:",NI,"\n")
#cat("Alpha [1-(Pn*Ds)/(Ps*Dn)]:",alphax,"\n")    

return(list(Ps=Ps,Pn=Pn,Ds=Ds,Dn=Dn,NI=NI,alphax=alphax,G=G,P=P))

}
