findsplitlev<-function(x,rec,beg,end,obspoint,suppo,n,lambda)
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
leftorright<-matrix(0,d,1)    #                     puolen valinta

suppvolume<-massone(suppo)
i<-1
while (i<=d){     #kaydaan muuttujat lapi
  supplen<-suppo[2*i]-suppo[2*i-1]    #alkup. jaettavan valin pituus
  valipit<-rec[2*i]-rec[2*i-1]

  jpislkm<-floor((n+1)*(valipit/supplen))-1

  if (jpislkm>=1){   #jos voidaan jakaa
    ssrvec2<-matrix(1,jpislkm,1) #ssrvec2:een talletet. kuhunkin mahd.
                  #jakopisteeseen liittyva ssr arvo i:nnelle muuttujalle 
    lefends<-matrix(1,jpislkm,1) 
    lrchoices<-matrix(0,jpislkm,1)

    ordobs<-makeorder(obs,x,i) #order pointers acc. to the i:th coord.
    dordobs[i,]<-ordobs
    leftend<-0  #kullekin jakopisteelle pointer to leftwing observations

    for (j in 1:jpislkm){ #kayd i:nnelle muuttujalle mahd. jakopist. lapi

      jakopiste<-supplen*j/(n+1)

      #hila jaetaan vasempaan ja oikeaan hilaan
      leftrec<-rec 
      leftrec[2*i-1]<-rec[2*i-1]    
      leftrec[2*i]<-rec[2*i-1]+jakopiste 
 
      leftend<-findobs(x,leftrec,ordobs,leftend,i)
      leftobslkm<-leftend
      lefends[j]<-leftend
      rightobslkm<-obslkm-leftobslkm
 
      rightrec<-rec
      rightrec[2*i-1]<-rec[2*i-1]+jakopiste                  
      rightrec[2*i]<-rec[2*i]  
 
      volumeleft<-massone(leftrec)
      volumeright<-massone(rightrec)

      leftchoice<-exma(volumeleft,leftobslkm,n,lambda)
      rightchoice<-exma(volumeright,rightobslkm,n,lambda)
  
      ssrvec2[j]<-min(leftchoice,rightchoice)

      if (ssrvec2[j]==leftchoice) lrchoices[j]<-0 else lrchoices[j]<-1

    }

    minvali<-omaind(-ssrvec2) #indeksi, jossa ssr:n suurin arvo i. muuttuj.

    valvec[i]<-rec[2*i-1]+supplen*minvali/(n+1)

    ssrvec1[i]<-ssrvec2[minvali]        #min(ssrvec2)
    dleftend[i]<-lefends[minvali]
    leftorright[i]<-lrchoices[minvali]
  }

  else ssrvec1[i]<-NA

  i<-i+1  
}

vec<-omaind(-ssrvec1)     #sen muuttujan numero joka halkaistu
val<-valvec[vec]          #halkaisupiste
leftrec<-rec
leftrec[2*vec]<-val       #vas rec:n loppupiste
rightrec<-rec
rightrec[2*vec-1]<-val    #oik rec:n alkupiste
lorr<-leftorright[vec]

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


resu<-list(val=val,vec=vec,leftrec=leftrec,rightrec=rightrec,leftbeg=leftbeg,
leftend=leftend,rightbeg=rightbeg,rightend=rightend,obspoint=obspoint,
lorr=lorr)
return(resu)
}
    










