findsplitG<-function(dendat,currec,curbeg,curend,obspoint,suppo,n,method){
#Finds a splitting point.

#dendat is n*d data-matrix
#currec is 2*d-vector
#curbeg, curend in 1:n, curbeg<curend, pointers to pointers
#obspoint is n-vector, points to rows of dendat

#Returns
#list(val,vec,leftrec,rightrec,leftbeg,leftend,rightbeg,rightend,obspoint)

#mahd. puolituspisteet ovat a+(b-a)*j/(n+1), j=1,...,n
#currec koostuu valeista muotoa [a0,b0]
#a0=a+(b-a)*j0/(n+1), b0=a+(b-a)*j1/(n+1), 1<= j0 < j1 <= n
#siis mahd jakopisteiden lkm on [(n+1)*(b0-a0)/(b-a)]-1 = j1-j0-1 >= 0 

d<-length(dendat[1,])        #x-muuttujien lkm
n<-length(dendat[,1])

step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(suppo[2*i]-suppo[2*i-1])/n
}

obs<-obspoint[curbeg:curend]
obslkm<-curend-curbeg+1

ssrvec1<-matrix(1,d,1)  #ssrvec1:een talletetaan kullekin muuttujalle
                        #pienin ssr-arvo (yli ko. muuttujan mahdollisiin
                        #jakopisteisiin liittyvista ssr arvoista).
dordobs<-matrix(0,d,obslkm)   #kullekin muuttujalle ordered obs
dleftend<-matrix(0,d,1)       #kullekin muuttujalle end of left 
valvec<-matrix(1,d,1)         #kullekin muuttujalle jakopiste
suppvolume<-massone(suppo)
i<-1
while (i<=d){     #kaydaan muuttujat lapi
  supplen<-suppo[2*i]-suppo[2*i-1]    #alkup. jaettavan valin pituus
  valipit<-(currec[2*i]-currec[2*i-1])*step[i]

  jpislkm<-currec[2*i]-currec[2*i-1]-1   #floor((n+1)*(valipit/supplen))-1

  if (jpislkm>=1){   #jos voidaan jakaa
    ssrvec2<-matrix(1,jpislkm,1) #ssrvec2:een talletet. kuhunkin mahd.
                     #jakopisteeseen liittyva ssr arvo i:nnelle muuttujalle 

    ordobs<-makeorder(obs,dendat,i) #order pointers acc. to the i:th coord.
    dordobs[i,]<-ordobs

    lefends<-matrix(1,jpislkm,1)
    j<-1 
    while (j<=jpislkm){  #kayd i:nnelle muuttujalle mahd. jakopist. lapi

      jakopiste<-j    #supplen*j/(n+1)

      #hila jaetaan vasempaan ja oikeaan hilaan
      leftrec<-currec 
      leftrec[2*i-1]<-currec[2*i-1]    
      leftrec[2*i]<-currec[2*i-1]+jakopiste 
      
      rightrec<-currec
      rightrec[2*i-1]<-currec[2*i-1]+jakopiste                  
      rightrec[2*i]<-currec[2*i]  
  
      leftrecend<-suppo[2*i-1]+leftrec[2*i]*step[i]
      leftend<-findobsR(dendat,ordobs,leftrecend,i)
      leftobslkm<-leftend
      lefends[j]<-leftend
      rightobslkm<-obslkm-leftobslkm
      
      volumeleft<-1
      for (ji in 1:d)
         volumeleft<-volumeleft*(leftrec[2*ji]-leftrec[2*ji-1])*step[ji]
      volumeright<-1
      for (ji in 1:d)
         volumeright<-volumeright*(rightrec[2*ji]-rightrec[2*ji-1])*step[ji]
      
      #meanleft<-denmean(volumeleft,leftobslkm,n) #vas hilan estim arvo
      #meanright<-denmean(volumeright,rightobslkm,n) #oik hilan estim arvo

      ssrvec2[j]<-denssr(volumeleft,leftobslkm,n,method)+denssr(volumeright,rightobslkm,n,method)
      j<-j+1
    }

    minvali<-omaind(-ssrvec2) #indeksi, jossa ssr:n suurin arvo i. muuttuj.

    valvec[i]<-minvali                  #currec[2*i-1]+supplen*minvali/(n+1)
    ssrvec1[i]<-ssrvec2[minvali]        #min(ssrvec2)
    dleftend[i]<-lefends[minvali]
  }

  else ssrvec1[i]<-NA

  i<-i+1  
}
vec<-omaind(-ssrvec1)     #sen muuttujan numero joka halkaistu
val<-currec[2*vec-1]+valvec[vec]          #halkaisupiste
leftrec<-currec
leftrec[2*vec]<-val   #vas rec:n loppupiste
rightrec<-currec
rightrec[2*vec-1]<-val  #oik rec:n alkupiste

if (dleftend[vec]==0){
  leftbeg<-0
  leftend<-0
  rightbeg<-curbeg
  rightend<-curend
  }
else if (dleftend[vec]==obslkm){
  leftbeg<-curbeg
  leftend<-curend
  rightbeg<-0
  rightend<-0
  }
else{
  leftbeg<-curbeg
  leftend<-curbeg+dleftend[vec]-1
  rightbeg<-leftend+1
  rightend<-curend  
}
obspoint[curbeg:curend]<-dordobs[vec,]

    #123 127
    #kapu<-length(obspoint)
    #juo<-1
    #while (juo<=kapu){
    #   if (obspoint[juo]==0){
    #         stop(format(supplen))
    #   }
    #   juo<-juo+1
    #}

resu<-list(val=val,vec=vec,leftrec=leftrec,rightrec=rightrec,leftbeg=leftbeg,
leftend=leftend,rightbeg=rightbeg,rightend=rightend,obspoint=obspoint)
return(resu)
}
    










