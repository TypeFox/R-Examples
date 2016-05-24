densplitF<-function(x,leaf,minlkm,suppo,
method="loglik",blokki=50,splitscan=0,seedf=1)
{
d<-length(x[1,])  #muuttujien lkm
n<-length(x[,1])  #havaintojen lkm

suppvol<-massone(suppo)
minvolume<-suppvol/(n+1)^d

maxnoden<-2*n   #suurin arviotu mahd silmujen lkm, n(1+1/2+1/4+...+1/n) 
maxpinen<-2*n

val<-matrix(1,maxnoden)       #huom haluttaisiin vektoreita
vec<-matrix(1,maxnoden)
mean<-matrix(1,maxnoden)
ssr<-matrix(1,maxnoden)
nelem<-matrix(1,maxnoden)
volume<-matrix(1,maxnoden)
left<-matrix(0,maxnoden)
right<-matrix(0,maxnoden)
low<-matrix(0,maxnoden,d)
upp<-matrix(0,maxnoden,d)

obspoint<-seq(1:n)               #pointers to the data (rows of x)
pinoparen<-matrix(0,maxpinen,1)
pinorecs<-matrix(0,maxpinen,d*2) #osioiden maaritelmat
pinopoint<-matrix(0,maxpinen,2)  #pointers to the pointers: location where the
                          #pointers to the datapoints are, for each rectangle
pinin<-1
pinoparen[pinin]<-0
pinorecs[pinin,]<-suppo
pinopoint[pinin,]<-c(1,n) 

# paaluuppi    do-until (pinin=0)

curleafnum<-1
curin<-0
while (pinin>=1){

  #  otetaan pinon paallimmainen tulokseen
  curin<-curin+1
  if (curin>maxnoden){
     val<-blokitus(val,blokki)
     vec<-blokitus(vec,blokki)
     nelem<-blokitus(nelem,blokki)  
     volume<-blokitus(volume,blokki)
     mean<-blokitus(mean,blokki)
     ssr<-blokitus(ssr,blokki)  
     left<-blokitus(left,blokki)
     right<-blokitus(right,blokki)
     low<-blokitus(low,blokki)
     upp<-blokitus(upp,blokki)
     maxnoden<-maxnoden+blokki
  }
  #
  curparent<-pinoparen[pinin]
  currec<-pinorecs[pinin,]    
  curpoint<-pinopoint[pinin,]
  curbeg<-curpoint[1]
  curend<-curpoint[2]
  pinin<-pinin-1
  #
  if (curparent>0) right[curparent]<-curin
  val[curin]<-NA                #aluksi ei puolitettu missaan kohdassa
  vec[curin]<-NA                #aluksi ei puolitettu mitaan muuttujaa
  nelem[curin]<-count(curbeg,curend)                   #havaintojen lkm
  volume[curin]<-massone(currec)
  mean[curin]<-denmean(volume[curin],nelem[curin],n)#,suppvol)#estim arv osiossa
  ssr[curin]<-denssr(volume[curin],nelem[curin],n,method)   #log likeli
  for (dimi in 1:d){
    low[curin,dimi]<-currec[2*dimi-1]
    upp[curin,dimi]<-currec[2*dimi]
  }

  #  jatketaan vas. alipuuhun
 
  jpistesti<-FALSE
  for (trun in 1:d){
     supplen<-suppo[2*trun]-suppo[2*trun-1]    #alkup. jaettavan valin pituus
     valipit<-currec[2*trun]-currec[2*trun-1]
     jpislkm<-floor((n+1)*(valipit/supplen))-1
     if (jpislkm>=1) jpistesti<-TRUE
  }
  #jpistesti<-(volume[curin]>=minvolume)  #????

  if (leaf==0){
    test<-((nelem[curin]>minlkm) && (jpistesti))
  }
  else{
    test<-((curleafnum<leaf) && (nelem[curin]>minlkm) && (jpistesti)) 
  }
  while (test){
    #  koska solmu jaettava, tehdaan jako
    #
    curleafnum<-curleafnum+1
    if (splitscan==0){
      jako<-findsplit(x,currec,curbeg,curend,obspoint,suppo,n,method)  
    }
    else{
       jako<-findsplit(x,currec,curbeg,curend,obspoint,suppo,n,method)
             #splitscan,seedf)  
    }

    left[curin]<-curin+1
    val[curin]<-jako$val
    vec[curin]<-jako$vec
    #
    rightrec<-jako$rightrec
    leftrec<-jako$leftrec
    leftbeg<-jako$leftbeg
    leftend<-jako$leftend
    rightbeg<-jako$rightbeg
    rightend<-jako$rightend
    obspoint<-jako$obspoint

    #    oikea lapsi paivitetaan pinoon

    pinin<-pinin+1
    if (pinin>maxpinen){
     pinoparen<-blokitus(pinoparen,blokki)
     pinorecs<-blokitus(pinorecs,blokki)
     pinopoint<-blokitus(pinopoint,blokki)
     maxpinen<-maxpinen+blokki
    }
    pinoparen[pinin]<-curin
    pinorecs[pinin,]<-rightrec
    pinopoint[pinin,]<-c(rightbeg,rightend)
    #    vasen lapsi paivitetaan tulokseen
    curin<-curin+1
    if (curin>maxnoden){
     val<-blokitus(val,blokki)
     vec<-blokitus(vec,blokki)
     nelem<-blokitus(nelem,blokki)  
     volume<-blokitus(volume,blokki)
     mean<-blokitus(mean,blokki)
     ssr<-blokitus(ssr,blokki)  
     left<-blokitus(left,blokki)
     right<-blokitus(right,blokki)
     low<-blokitus(low,blokki)
     upp<-blokitus(upp,blokki)
     maxnoden<-maxnoden+blokki
    }
    currec<-leftrec
    curbeg<-leftbeg
    curend<-leftend
    val[curin]<-NA
    vec[curin]<-NA
    nelem[curin]<-count(curbeg,curend)
    volume[curin]<-massone(currec)
    mean[curin]<-denmean(volume[curin],nelem[curin],n)#,suppvol)
    ssr[curin]<-denssr(volume[curin],nelem[curin],n,method)
    for (dimi in 1:d){
      low[curin,dimi]<-leftrec[2*dimi-1]
      upp[curin,dimi]<-leftrec[2*dimi]
    }

    if (leaf==0){
      test<-((nelem[curin]>minlkm) && (volume[curin]>=minvolume))
    } 
    else{
       test<-((curleafnum<leaf) &&
              (nelem[curin]>minlkm) && (volume[curin]>=minvolume)) 
    }

  }

}
val<-val[1:curin]  #tassa matriisi muuntuu vektoriksi!!!
vec<-vec[1:curin]
volume<-volume[1:curin]
mean<-mean[1:curin]
ssr<-ssr[1:curin]
nelem<-nelem[1:curin]
left<-left[1:curin]
right<-right[1:curin]
low<-low[1:curin,]
upp<-upp[1:curin,]
puu<-list(val=val,vec=vec,mean=mean,nelem=nelem,ssr=ssr,volume=volume,
left=left,right=right,
low=low,upp=upp)
return(puu)
}









