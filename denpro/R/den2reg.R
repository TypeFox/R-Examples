den2reg<-function(dendat,binlkm,kantaja)
{
#muuntaa tiheysdatan regressiodataksi

#dendat on tiheysdatan sisaltava n*xlkm matriisi
#binlkm on niitten lokeroitten maara joihin yksi muuttuja ositetaan
#otosavaruus ositetaan binlkm^p lokeroon
#kantaja on 2*xlkm vaakavektori, sisaltaa kantajan kullekin 
#muuttujalle, olet etta kaikki dendat:n hav todella sisaltyvat kantajaan.

#palauttaa: listan vast
#vast.dep on saatu*1 vektori, sisaltaa frekvenssit,
#(saatu on erillisten diskretoitujen havaintojen lkm) 
#vast.ind on saatu*xlkm matriisi, sisaltaa niitten lokeroitten 
#keskipisteet joita (diskretoidussa) datassa esiintyy vahintaan yksi, 
#vast.hila on xlkm*3 matriisi, 
#vast.hila:n i:s rivi sis. i:nnelle muuttujalle 
#pisteiden maaran
#ensimmaisen regressiodatan havaintopisteen,
#viimeisen regressiodatan havaintopisteen.

xlkm<-length(dendat[1,]) #dendat:in sarakkeitten lkm on muuttujien lkm
n<-length(dendat[,1])    #dendat:in rivien lkm on havaintojen maara
hila<-matrix(1,xlkm,3)   #hila on xlkm*3 matriisi
binhila<-hila            #binhila on xlkm*3 matriisi
valpit<-matrix(1,xlkm,1)      #hilan valien pituudet
hila[,1]<-binlkm*hila[,1]
binhila[,1]<-binlkm*binhila[,1]
i<-1
while (i<=xlkm){
     binhila[i,2]<-kantaja[i,1]   #min(dendat[,i])-epsi 
     binhila[i,3]<-kantaja[i,2]     #max(dendat[,i])+epsi 
     valpit[i]<-(binhila[i,3]-binhila[i,2])/binlkm
       #binhila:n i:s rivi sis. i:nnelle muuttujalle 
       #lokeroinnin alkupisteen, arvoalueen loppupisteen   
       #valpit: arvoalueen pituus jaettuna lokeroitten lkm:lla
       #eli yhden lokeron leveys. 
     i<-i+1
}
hila[,2]<-binhila[,2]+valpit/2
hila[,3]<-binhila[,3]-valpit/2
       #hila:n i:s rivi sis. i:nnelle muuttujalle 
       #ensimmaisen regressiodatan havaintopisteen,
       #viimeisen regressiodatan havaintopisteen
#if (valpit<=0) stop("in some variable there is no variation")
hiladat<-matrix(1,n,xlkm) #muunnetaan dendat hiladat:iksi
                          #ts. pyoristetaan havainnot 
                          #hiladat sis. diskretoidut havainnot
i<-1
while (i<=n){       #kaydaan lapi aineisto
    #pyoristetaan i:s havainto hilapisteeseen
    j<-1
    while (j<=xlkm){
       alavali<-floor((dendat[i,j]-binhila[j,2])/valpit[j])
           #alavali ilmaisee monennessako lokerossa hav. sijaitsee
       hiladat[i,j]<-binhila[j,2]+alavali*valpit[j]+valpit[j]/2
       j<-j+1
    }
i<-i+1
}
xtulos<-matrix(0,n,xlkm) #periaatteessa mahdollista etta kaikki n
                         #havaintoa ovat eri lokeroissa, siksi
                         #laitetaan xtulos matriisin rivien maaraksi n
ytulos<-matrix(0,n,1)
xtulos[1,]<-hiladat[1,]  #hiladat:in ensimmainen rivi esiintyy ainakin kerran
ytulos[1]<-1             #sen frekvenssi ainakin yksi
saatu<-1                 #toistaiseksi yksi erillinen havainto
i<-1
while (i<n){  #kaydaan lapi aineisto
   i<-i+1     #aloitetaan kakkosesta
   lippu<-0   #apriori kyseessa uusi lajityyppi
   j<-1
   while ((j<=saatu) && (lippu==0)){  #kaydaan lapi keratyt havinnot
       if (all(hiladat[i,]==xtulos[j,])){  #jos on jo saatu
            lippu<-1     #liputetaan etta havaittiin toisto
            jind<-j      #merkataan indeksi frekvenssin paivitysta varten
       }
       j<-j+1       
   }
   if (lippu==1) ytulos[jind]<-ytulos[jind]+1 
            #jos saatiin toisto, paivitetaan frekvenssi 
      else{ 
         saatu<-saatu+1      #jos saatiin uusi, lisataan saatu:un yksi ja
         xtulos[saatu,]<-hiladat[i,]  #merkitaan uusi lajityyppi muistiin
         ytulos[saatu]<-1    #uuden lajityypin frekvenssi on aluksi yksi
      }
}
xtulos<-xtulos[1:saatu,]
ytulos<-ytulos[1:saatu]
ytulos<-t(t(ytulos))
if (xlkm==1) xtulos<-t(t(xtulos))
return(list(dep=ytulos,ind=xtulos,hila=hila))
}

