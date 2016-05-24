findsplitter2<-function(x,rec,n,method)
{
# x is m*d matrix of split points (and boundaries)
# rec is 2*d-vector of pointers to x

d<-length(x[1,])        
ssrvec1<-matrix(1,d,1)  #ssrvec1:een talletetaan kullekin muuttujalle
                        #pienin ssr-arvo (yli ko. muuttujan mahdollisiin
                        #jakopisteisiin liittyvista ssr arvoista).
valvec<-matrix(1,d,1)   #kullekin muuttujalle jakopiste

i<-1
while (i<=d){                         #kaydaan muuttujat lapi
  jpislkm<-rec[2*i]-rec[2*i-1]-1
  if (jpislkm>=1){   #jos voidaan jakaa
    ssrvec2<-matrix(1,jpislkm,1) #ssrvec2:een talletet. kuhunkin mahd.
                  #jakopisteeseen liittyva ssr arvo i:nnelle muuttujalle 
    for (j in 1:jpislkm){ #kayd i:nnelle muuttujalle mahd. jakopist. lapi
        jakopiste<-rec[2*i-1]+j
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
        trueleftrec<-leftrec
        for (dd in 1:d){ 
            trueleftrec[2*dd-1]<-x[leftrec[2*dd-1],dd]
            trueleftrec[2*dd]<-x[leftrec[2*dd],dd]
        }
        truerightrec<-leftrec
        for (dd in 1:d){ 
            truerightrec[2*dd-1]<-x[rightrec[2*dd-1],dd]
            truerightrec[2*dd]<-x[rightrec[2*dd],dd]
        }
        volumeleft<-massone(trueleftrec)
        volumeright<-massone(truerightrec)
        #
        ssrvec2[j]<-denssr(volumeleft,leftobslkm,n,method)
                   +denssr(volumeright,rightobslkm,n,method)
    }

    minvali<-omaind(-ssrvec2) #indeksi, jossa ssr:n suurin arvo i. muuttuj.
    valvec[i]<-rec[2*i-1]+minvali
    ssrvec1[i]<-ssrvec2[minvali]        # min(ssrvec2)
  }

  else ssrvec1[i]<-NA

  i<-i+1  
}
vec<-omaind(-ssrvec1)     # sen muuttujan numero joka halkaistu
val<-valvec[vec]          # halkaisupiste

#leftrec<-rec
#leftrec[2*vec]<-val   #vas rec:n loppupiste
#rightrec<-rec
#rightrec[2*vec-1]<-val  #oik rec:n alkupiste

resu<-list(val=val,vec=vec)
return(resu)
}
    






