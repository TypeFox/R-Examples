findsplitter.dyadic<-function(grid,x,rec,n,method,minobs,recdown,rechigh)   
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
     j<-ceiling(jpislkm/2)

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
       ssr2<-NA
     else
       ssr2<-denssr(volumeleft,leftobslkm,n,method)+denssr(volumeright,rightobslkm,n,method)

     minvali<-j
     valvecint[i]<-recdown[i]+minvali
     valvec[i]<-grid[valvecint[i],i]
     ssrvec1[i]<-ssr2
  }

  else ssrvec1[i]<-NA

  i<-i+1  
}
vec<-omaind(-ssrvec1)     # sen muuttujan numero joka halkaistu
val<-valvec[vec]          # halkaisupiste
valio<-valvecint[vec]

resu<-list(val=val,vec=vec,valio=valio)
return(resu)
}
    






