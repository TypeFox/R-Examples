findsplit<-function(x,rec,beg,end,obspoint,suppo,n,method)
{
#Finds a splitting point.

#x is n*d data-matrix
#rec is 2*d-vector
#beg, end in 1:n, beg<end, pointers to pointers
#obspoint is n-vector, points to rows of x

#Returns
#list(val,vec,leftrec,rightrec,leftbeg,leftend,rightbeg,rightend,obspoint)

#mahd. puolituspisteet ovat a+(b-a)*j/(n+1), j=1,...,n
#rec koostuu valeista muotoa [a0,b0]
#a0=a+(b-a)*j0/(n+1), b0=a+(b-a)*j1/(n+1), 1<= j0 < j1 <= n
#siis mahd jakopisteiden lkm on [(n+1)*(b0-a0)/(b-a)]-1 = j1-j0-1 >= 0 

d<-length(x[1,])        #x-muuttujien lkm
n<-length(x[,1])
obs<-obspoint[beg:end]

obslkm<-end-beg+1

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
  valipit<-rec[2*i]-rec[2*i-1]
#
  jpislkm<-floor((n+1)*(valipit/supplen))-1
#  jpislkm<-obslkm-1
#

  if (jpislkm>=1){   #jos voidaan jakaa
    ssrvec2<-matrix(1,jpislkm,1) #ssrvec2:een talletet. kuhunkin mahd.
                  #jakopisteeseen liittyva ssr arvo i:nnelle muuttujalle 

    ordobs<-makeorder(obs,x,i) #order pointers acc. to the i:th coord.
    dordobs[i,]<-ordobs

    lefends<-matrix(1,jpislkm,1) 
    leftend<-0  #kullekin jakopisteelle pointer to leftwing observations
    for (j in 1:jpislkm){ #kayd i:nnelle muuttujalle mahd. jakopist. lapi
#
     jakopiste<-supplen*j/(n+1)
#      jakopiste<-valipit*j/(jpislkm+1) 
#
      #hila jaetaan vasempaan ja oikeaan hilaan
      leftrec<-rec 
      leftrec[2*i-1]<-rec[2*i-1]    
      leftrec[2*i]<-rec[2*i-1]+jakopiste 
      #
      leftend<-findobs(x,leftrec,ordobs,leftend,i)
      leftobslkm<-leftend
      lefends[j]<-leftend
      rightobslkm<-obslkm-leftobslkm
      #
      rightrec<-rec
      rightrec[2*i-1]<-rec[2*i-1]+jakopiste                  
      rightrec[2*i]<-rec[2*i]  
      #
#      volumeleft<-suppvolume*(leftrec[2*i]-leftrec[2*i-1])/
#                  (suppo[2*i]-suppo[2*i-1])
      volumeleft<-massone(leftrec)
      meanleft<-denmean(volumeleft,leftobslkm,n) 
                         #vas hilan estim arvo
#      volumeright<-suppvolume*(rightrec[2*i]-rightrec[2*i-1])/
#                  (suppo[2*i]-suppo[2*i-1])
      volumeright<-massone(rightrec)
      meanright<-denmean(volumeright,rightobslkm,n) 
                         #oik hilan estim arvo
      #
      ssrvec2[j]<-denssr(volumeleft,leftobslkm,n,method)+denssr(volumeright,rightobslkm,n,method)
    }

    minvali<-omaind(-ssrvec2) #indeksi, jossa ssr:n suurin arvo i. muuttuj.
#
    valvec[i]<-rec[2*i-1]+supplen*minvali/(n+1)
#    valvec[i]<-rec[2*i-1]+valipit*minvali/(jpislkm+1)
#
    ssrvec1[i]<-ssrvec2[minvali]        #min(ssrvec2)
    dleftend[i]<-lefends[minvali]
  }

  else ssrvec1[i]<-NA

  i<-i+1  
}
vec<-omaind(-ssrvec1)     #sen muuttujan numero joka halkaistu
val<-valvec[vec]          #halkaisupiste
leftrec<-rec
leftrec[2*vec]<-val   #vas rec:n loppupiste
rightrec<-rec
rightrec[2*vec-1]<-val  #oik rec:n alkupiste

if (dleftend[vec]==0){
  leftbeg<-NA
  leftend<-NA
  rightbeg<-beg
  rightend<-end
  }
else if (dleftend[vec]==obslkm){
  leftbeg<-beg
  leftend<-end
  rightbeg<-NA
  rightend<-NA
  }
else{
  leftbeg<-beg
  leftend<-beg+dleftend[vec]-1
  rightbeg<-leftend+1
  rightend<-end  
}
obspoint[beg:end]<-dordobs[vec,]

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
    










