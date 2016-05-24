findsplitter3<-function(x,rec,n,method,recdown)   
{
# x is m*d matrix of data
# rec is 2*d-vector

l<-dim(x)[1]
d<-dim(x)[2]

xx<-matrix(0,l-1,d)  # x contains split points 
for (i in 1:d){
    ordi<-order(x[,i])
    ala<-x[ordi[1:(l-1)],i]
    yla<-x[ordi[2:l],i]
    xx[,i]<-(ala+yla)/2
}

ssrvec1<-matrix(1,d,1)  #ssrvec1:ssa pienin ssr-arvo (yli ko. muuttujan mahdollisiin
                        #jakopisteisiin liittyvista ssr arvoista).
valvec<-matrix(1,d,1)   #kullekin muuttujalle jakopiste
valvecint<-matrix(1,d,1)   

i<-1
while (i<=d){           #kaydaan muuttujat lapi
  jpislkm<-l-1
  if (jpislkm>=1){   #jos voidaan jakaa
    ssrvec2<-matrix(1,jpislkm,1) #ssrvec2:een talletet. kuhunkin mahd.
                  #jakopisteeseen liittyva ssr arvo i:nnelle muuttujalle 
    for (j in 1:jpislkm){ #kayd i:nnelle muuttujalle mahd. jakopist. lapi
        jakopiste<-xx[j,i]
        #jakpisteint<-recdown[i]+j
        # hila jaetaan vasempaan ja oikeaan hilaan
        leftrec<-rec 
        leftrec[2*i-1]<-rec[2*i-1]    
        leftrec[2*i]<-jakopiste
        rightrec<-rec
        rightrec[2*i-1]<-jakopiste                  
        rightrec[2*i]<-rec[2*i]  
        #
        leftobslkm<-j
        rightobslkm<-jpislkm-j+1
        #
        volumeleft<-massone(leftrec)
        volumeright<-massone(rightrec)
        #
        ssrvec2[j]<-denssr(volumeleft,leftobslkm,n,method)+denssr(volumeright,rightobslkm,n,method)
    }

    minvali<-omaind(-ssrvec2) #indeksi, jossa ssr:n suurin arvo i. muuttuj.
    valvec[i]<-xx[minvali,i]
    valvecint[i]<-recdown[i]+minvali
    ssrvec1[i]<-ssrvec2[minvali]        # min(ssrvec2)
  }

  else ssrvec1[i]<-NA

  i<-i+1  
}
vec<-omaind(-ssrvec1)     # sen muuttujan numero joka halkaistu
val<-valvec[vec]          # halkaisupiste
valio<-valvecint[vec]

#leftrec<-rec
#leftrec[2*vec]<-val   #vas rec:n loppupiste
#rightrec<-rec
#rightrec[2*vec-1]<-val  #oik rec:n alkupiste

resu<-list(val=val,vec=vec,valio=valio)
return(resu)
}
    






