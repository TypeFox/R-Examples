findsplitpenaR<-function(dendat,currec,curbeg,curend,obspoint,
mcdendat,mccurbeg,mccurend,mcobspoint,mix,suppo,step)
{

d<-length(dendat[1,])        #x-muuttujien lkm
n<-length(dendat[,1])
mcn<-length(mcdendat[,1])

obs<-obspoint[curbeg:curend]
obslkm<-curend-curbeg+1

mcobs<-mcobspoint[mccurbeg:mccurend]
mcobslkm<-mccurend-mccurbeg+1

ssrvec1<-matrix(1,d,1)  #ssrvec1:een talletetaan kullekin muuttujalle
                        #pienin ssr-arvo (yli ko. muuttujan mahdollisiin
                        #jakopisteisiin liittyvista ssr arvoista).
dordobs<-matrix(0,d,obslkm)   #kullekin muuttujalle ordered obs
mcdordobs<-matrix(0,d,mcobslkm)   #kullekin muuttujalle ordered obs

dleftend<-matrix(0,d,1)       #kullekin muuttujalle end of left 
mcdleftend<-matrix(0,d,1)   

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

    mcordobs<-makeorder(mcobs,mcdendat,i)
    mcdordobs[i,]<-mcordobs

    lefends<-matrix(1,jpislkm,1)
    mclefends<-matrix(1,jpislkm,1)
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
 
      mcleftend<-findobsR(mcdendat,mcordobs,leftrecend,i)
      mcleftobslkm<-mcleftend
      mclefends[j]<-mcleftend
      mcrightobslkm<-mcobslkm-mcleftobslkm
     
      volumeleft<-1
      for (ji in 1:d)
         volumeleft<-volumeleft*(leftrec[2*ji]-leftrec[2*ji-1])*step[ji]
      volumeright<-1
      for (ji in 1:d)
         volumeright<-volumeright*(rightrec[2*ji]-rightrec[2*ji-1])*step[ji]
      
      meanleft<-denmean(volumeleft,leftobslkm,n) #vas hilan estim arvo
      meanright<-denmean(volumeright,rightobslkm,n) #oik hilan estim arvo

      ssrvec2[j]<-denssr(volumeleft,leftobslkm,n,method="projec",mix=mix)+denssr(volumeright,rightobslkm,n,method="projec",mix=mix)+2*(mix-1)*mix*mcleftobslkm*meanleft/mcn+2*(mix-1)*mix*mcrightobslkm*meanright/mcn

      j<-j+1
    }

    minvali<-omaind(-ssrvec2) #indeksi, jossa ssr:n suurin arvo i. muuttuj.

    valvec[i]<-minvali                  #currec[2*i-1]+supplen*minvali/(n+1)
    ssrvec1[i]<-ssrvec2[minvali]        #min(ssrvec2)
    dleftend[i]<-lefends[minvali]
    mcdleftend[i]<-mclefends[minvali]

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

if (mcdleftend[vec]==0){
  mcleftbeg<-0
  mcleftend<-0
  mcrightbeg<-mccurbeg
  mcrightend<-mccurend
  }
else if (mcdleftend[vec]==mcobslkm){
  mcleftbeg<-mccurbeg
  mcleftend<-mccurend
  mcrightbeg<-0
  mcrightend<-0
  }
else{
  mcleftbeg<-mccurbeg
  mcleftend<-mccurbeg+mcdleftend[vec]-1
  mcrightbeg<-mcleftend+1
  mcrightend<-mccurend  
}

obspoint[curbeg:curend]<-dordobs[vec,]
mcobspoint[mccurbeg:mccurend]<-mcdordobs[vec,]

resu<-list(val=val,vec=vec,leftrec=leftrec,rightrec=rightrec,
leftbeg=leftbeg,leftend=leftend,rightbeg=rightbeg,rightend=rightend,
obspoint=obspoint,
mcleftbeg=mcleftbeg,mcleftend=mcleftend,mcrightbeg=mcrightbeg,mcrightend=mcrightend,
mcobspoint=mcobspoint
)
return(resu)
}
    










