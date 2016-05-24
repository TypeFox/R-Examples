declevnew<-function(rindeksit,linkit,n){
#Splits level set to disconnected subsets
#
#rindeksit is m-vector, links to observations, levset is union of m atoms
#linkit is n*n-matrix, n is the total number of atoms,
#  describes which atoms touch each other
#
#Returns sublevnum*n-matrix, describes disconnected parts of levset  
#
m<-length(rindeksit)
tulos<-matrix(0,m,n)      #in the worst case we split to m parts
#                         #blokit !!!!!!!!!!!!!!!!
merkatut<-matrix(0,m,1) #laitetaan 1 jos on jo sijoitettu johonkin tasojouk.
pino<-matrix(0,m,1) #pinoon laitetaan aina jos koskettaa, max kosketuksia m
pinind<-0           #pinossa viitataan rindeksin elementteihin
curleima<-1
i<-1   #i ja j viittavat rindeksit-vektoriin, jonka alkiot viittavat atomeihin
while (i<=m){
  if (merkatut[i]==0){  #jos ei viela merkattu niin 
    pinind<-pinind+1    #pannaan pinoon
    pino[pinind]<-i    
    while (pinind>0){
      curviite<-pino[pinind]  #otetaan pinosta viite rindeksit-vektoriin
                              #jossa puolestaan viitteet itse palloihin
      curpallo<-rindeksit[curviite]  #haetaan viite palloon
      pinind<-pinind-1  
      tulos[curleima,curpallo]<-1    #laitetaan pallo ko tasojoukkoon
#      merkatut[curviite]<-1          #merkataan kaytetyksi
      j<-1
      while (j<=m){          #pannnaan linkeista pinoon
        ehdokas<-rindeksit[j]         #kaydaan ko tasojoukon atomit lapi
        touch<-(linkit[curpallo,ehdokas]==1)
        if ((touch) && (merkatut[j]==0)){
          pinind<-pinind+1      
          pino[pinind]<-j  
          merkatut[j]<-1
        }
        j<-j+1
      }
    }
    curleima<-curleima+1   #uusi leima
  }
  i<-i+1
} 
tullos<-matrix(0,(curleima-1),n)
tullos[1:(curleima-1),]<-tulos[1:(curleima-1),] #poistetaan ylimaaraiset
return(tullos)
}


















