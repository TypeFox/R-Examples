findsplitter<-function(grid,x,rec,n,method,minobs,recdown,rechigh)   
{
# x is m*d matrix of data
# rec is 2*d-vector

l<-dim(x)[1]
d<-dim(x)[2]

ssrvec1<-matrix(1,d,1)  #ssrvec1:ssa pienin ssr-arvo (yli ko. muuttujan mahdollisiin
                        #jakopisteisiin liittyvista ssr arvoista).
valvec<-matrix(1,d,1)   #kullekin muuttujalle jakopiste
valvecint<-matrix(1,d,1)   

i<-1
while (i<=d){           #kaydaan muuttujat lapi
  jpislkm<-rechigh[i]-recdown[i]-1
  if (jpislkm>=1){   #jos voidaan jakaa
    ssrvec2<-matrix(1,jpislkm,1) #ssrvec2:een kuhunkin jpistessr-arvo i:lle muutt.
    for (j in 1:jpislkm){ 
        jakopisteint<-recdown[i]+j
        jakopiste<-grid[jakopisteint,i]
        # hila jaetaan vasempaan ja oikeaan hilaan
        leftrec<-rec 
        leftrec[2*i-1]<-rec[2*i-1]    
        leftrec[2*i]<-jakopiste
        rightrec<-rec
        rightrec[2*i-1]<-jakopiste                  
        rightrec[2*i]<-rec[2*i]  
        #
        leftobslkm<-sum(x[,i]<jakopiste)
        rightobslkm<-sum(x[,i]>jakopiste)
        #
        volumeleft<-massone(leftrec)
        volumeright<-massone(rightrec)
        #
        if ((leftobslkm==0)||(rightobslkm==0))
        ssrvec2[j]<-NA
        else
        ssrvec2[j]<-denssr(volumeleft,leftobslkm,n,method)+denssr(volumeright,rightobslkm,n,method)
    }

    minvali<-omaind(-ssrvec2) #indeksi, jossa ssr:n suurin arvo i. muuttuj.
    valvecint[i]<-recdown[i]+minvali
    valvec[i]<-grid[valvecint[i],i]
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
    






