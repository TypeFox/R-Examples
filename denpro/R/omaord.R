omaord<-function(values,recs,frekv=NULL){
#Jarjestaa paloittain vakion funktion palat funktion arvojen
#mukaan suuruusjarjestykseen
#
#palvak on lnum*(1+2*xlkm)-matriisi, missa lnum on laatikkojen lkm,
#matriisin ensimmainen sarake sisaltaa estimaatin values laatikoittain,
#naitten mukaan matriisin rivit jarjestetaan.
#Muut sarakkeet sis laatikoitten maaritykset, ts jokaista
#muuttujaa kohden vaihteluvali 
#ep on toleranssiparametri yhtasuuruuden testauksessa
#
#kutsuu: omaind
#
lnum<-length(values)        #length(recs[,1])     #laatikoitten lkm
ordrecs<-recs               #tahan oikea jarjestys
ordvalues<-values
if (is.null(frekv)){
 ordfrekv<-NULL
 i<-1
 while (i<=lnum){
   pienin<-omaind(values)
   ordrecs[i,]<-recs[pienin,]
   ordvalues[i]<-values[pienin]
   values[pienin]<-NA       #NA on plus aareton
   i<-i+1
 }
}
else{
ordfrekv<-frekv
i<-1
 while (i<=lnum){
   pienin<-omaind(values)
   ordrecs[i,]<-recs[pienin,]
   ordvalues[i]<-values[pienin]
   ordfrekv[i]<-frekv[pienin]
   values[pienin]<-NA       #NA on plus aareton
   i<-i+1
 }
}
return(list(values=ordvalues,recs=ordrecs,frekv=ordfrekv))
}





