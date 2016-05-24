allocolo<-function(mlkm,pg,mlkmpre,pgpre,colopre,coloind,colofall,paletti)
{

lkmpre<-mlkmpre$lkm   #previous number of modes
lkm<-mlkm$lkm         #current number of modes
modecolo<-matrix("",lkm,1)

# calculate distances

if (!is.null(colopre)){
dist<-matrix(NA,lkm,lkmpre)   #NA is infty
for (i in 1:lkm){
   for (j in 1:lkmpre){
      cent<-pg$center[,mlkm$modloc[i]]
      centpre<-pgpre$center[,mlkmpre$modloc[j]]
      dist[i,j]<-etais(cent[1:2],centpre[1:2])
   }
}
}

# allocate colors

if (is.null(colopre)){
for (k in 1:lkm){
  modecolo[k]<-paletti[k]
}
coloind<-lkm
}
else if (lkm>=lkmpre){
for (k in 1:lkmpre){
  minimi<-min(dist,na.rm=TRUE)
  argmin<-which(minimi==dist)[1]
  yind<-ceiling(argmin/lkm)      
  xind<-argmin-(yind-1)*lkm       #index for first color
  dist[xind,]<-NA
  modecolo[k]<-colopre[yind]
}
k<-lkmpre+1
while (k<=lkm){
  modecolo[k]<-paletti[coloind+k-lkmpre]
  k<-k+1
}
coloind<-lkm
}
else{                 #lkm<lkmpre
for (k in 1:lkm){
   minimi<-min(dist,na.rm=TRUE)
   argmin<-which(minimi==dist)[1]
   #cur<-omaindmat(dist)
   yind<-ceiling(argmin/lkm)      
   xind<-argmin-(yind-1)*lkm       #index for first color
   dist[xind,]<-NA
   modecolo[k]<-colopre[yind]
}
colofall<-lkm
}

return(list(modecolo=modecolo,coloind=coloind,colofall=colofall))
}




allosimp<-function(vec,indvec,colot){
#
#minimi<-min(vec,na.rm=TRUE)
#cur<-which(vec==min(vec,na.rm=TRUE)) 
#
epsi<-0.0000001
reslen<-length(vec)
res<-matrix("",reslen,1)
#
len<-length(indvec)
for (i in 1:len){
   inde<-indvec[i]
   #vektori<-which(vec==vec[inde])
   #vektori<-which(((vec<=vec[inde]+epsi)&&(vec>=vec[inde]-epsi)))
   apu<-matrix(FALSE,reslen,1)
   for (l in 1:reslen){
     if ((vec[l]<=vec[inde]+epsi)&&(vec[l]>=vec[inde]-epsi)) apu[l]<-TRUE
   }
   vektori<-which(apu)
   len2<-length(vektori)
   for (j in 1:len2){
      sijo<-vektori[j]
      res[sijo]<-colot[i]
   }
}
#
return(res)
}
belongs<-function(obs,rec,coordi){
#Finds whether observation belongs to the rectangle
#
#obs is d-vector
#rec is 2*d-vector
#coordi is in 1:d
#
#Returns TRUE is obs is in rec, otherwise FALSE
#
#We need to check only whether coordi:s coordinate is in the
#interval
#
ans<-TRUE
if ((obs[coordi]<rec[2*coordi-1]) || (obs[coordi]>rec[2*coordi])) ans<-FALSE
return(ans)
}
blokitus.delt<-function(obj,blokki){
#
if (dim(t(obj))[1]==1) k<-1 else k<-length(obj[,1]) #rivien maara 
if (k==1){
  len<-length(obj)
  uusobj<-matrix(0,len+blokki,1)
  uusobj[1:len]<-obj
}
else{
  lev<-length(obj[1,])
  uusobj<-matrix(0,k+blokki,lev)
  uusobj[1:k,]<-obj
}
return(uusobj)
}
bootbagg<-function(dendat,seed,scatter=0)
{
set.seed(seed)
n<-dim(dendat)[1]
d<-dim(dendat)[2]
dendatout<-matrix(0,n,d)

for (i in 1:n){
   res<-ceiling(runif(1)*n)
   addi<-scatter*2*(runif(d)-0.5)
   dendatout[i,]<-dendat[res,]+addi
}

return(dendatout)
}
bootworpl<-function(dendat,seed,scatter=0)
{
set.seed(seed)
n<-dim(dendat)[1]
d<-dim(dendat)[2]
nout<-ceiling(n/2)
dendatout<-matrix(0,nout,d)

blacklist<-matrix(0,n,1)
found<-0

while (found<nout){
   res<-ceiling(runif(1)*n)
   if (blacklist[res]==0){
       found<-found+1
       addi<-scatter*2*(runif(d)-0.5)
       dendatout[found,]<-dendat[res,]+addi
       blacklist[res]<-1
   }

}

return(dendatout)
}
childcode<-function(level,runnumb){
#
#input: level and ordinal number in the listing of the full binary tree
#output: ordinal number of the left and right child
#
nodeNum<-0
i<-0
while (i<=(level-2)){
   nodeNum<-nodeNum+2^i
   i<-i+1
}
levBeg<-nodeNum+1
whichInLevel<-runnumb-levBeg+1
levBegNex<-levBeg+2^(level-1)
#
left<-levBegNex+(whichInLevel-1)*2
right<-left+1
#
return(list(left=left,right=right))
}
cluster.lst<-function(dendat,h,N=NULL,cut=NULL,lambda=NULL,complete=FALSE,
type="grid",labels="number",nodes=NULL,minobs=1)
{
# cut is in (0,1)
n<-dim(dendat)[1]

if (type=="grid"){ 
   pcf<-pcf.kern(dendat,h,N,radi=0)
   lst<-leafsfirst(pcf)
}
else{
    pcf<-pcf.greedy.kernel(dendat,h,minobs=minobs)
    lst<-leafsfirst.adagrid(pcf)
}
if (!is.null(cut)) clusterlevel<-cut*(max(lst$level)-min(lst$level))
if (!is.null(lambda)) clusterlevel<-lambda
if (!is.null(nodes)) clusterlevel<-NULL

if (type=="grid")
cd<-colors2data(dendat,pcf,lst,clusterlevel=clusterlevel,nodes=nodes,type="regular")
else
cd<-colors2data(dendat,pcf,lst,clusterlevel=clusterlevel,nodes=nodes,type="ada")

cls0<-cd$datacolo
mita0<-unique(cd$datacolo)
ind<-(mita0=="grey")
mita<-mita0
mita[length(mita0)]<-"grey"
mita[ind]<-mita0[length(mita0)]

clnum<-length(mita)
cnum<-clnum-1
nums<-matrix(0,cnum,1)
cls<-matrix(0,n,1)
for (i in 1:n){
    if (cls0[i]=="grey") cls[i]<-0
    else{ 
       for (j in 1:cnum) 
       if (mita[j]==cls0[i]){ 
          cls[i]<-j
          nums[j]<-nums[j]+1 
       }
    }
}

if (complete){
   indi<-(cls==0)
   data<-dendat[indi,]
   n0<-dim(data)[1]
   part<-matrix(0,n0,1)
   part0<-matrix("",n0,1)
   mito<-mita[1:cnum]
   for (i in 1:n0){
       arg<-data[i,]
       arvot<-matrix(0,cnum,1)
       for (j in 1:cnum){
           ota<-(cls==j)
           x<-dendat[ota,]
           arvot[j]<-(nums[j]/(n-n0))*kernesti.dens(arg,x,h=h)
       }
       makind<-(arvot==max(arvot))
       part[i]<-seq(1,cnum)[makind]
       part0[i]<-mito[makind]
   }
   cls[indi]<-part
   cls0[indi]<-part0
}

if (labels=="number") labels<-cls else labels<-cls0
return(labels)
}

colo2eprof<-function(ep,mt,as){
#
#ep result of scaletree
#mt result of modetree
#as result of allosimp: vector of colors
#
len<-length(ep$bigdepths)
mtlen<-length(as)
colors<-matrix("black",len,1)
#
for (i in 1:len){
  label<-ep$mlabel[i]       #label for mode
  if (label>0){ 
      smoot<-ep$smoot[i]    #smoothing paramter value/leafnum
      # we find the corresponding slot from "as" where
      # label corresponds and smmothing parameter value corresponds
      run<-1
      koesmoot<-mt$ycoor[run]
      koelabel<-mt$mlabel[run] 
      while (((koesmoot!=smoot) || (koelabel!=label)) && (run<=mtlen)){
         run<-run+1
         koesmoot<-mt$ycoor[run]
         koelabel<-mt$mlabel[run] 
      }
      # we have found the slot
      colors[i]<-as[run]
  }
}
#
return(colors)
}
count<-function(beg,end){
#Laskee hilaan kuuluvien pisteiden lkm:n
#
#hila on 2-vector
#
#returns integer >1
#
if (is.na(beg)) ans<-0              #or end==NA
else if ((beg==0) && (end==0)) ans<-0
else ans<-end-beg+1
return(ans)
}
denmean<-function(volume,obslkm,n)
{
#Calculates the value of the estimate at rectangle

#volum is >0
#obslkm is lkm-vector, pointers to the rows of x
#n is the sample size

#returns non-negative real number

#Value of the estimate is the number of observations
#in rec divided by the total number of obs, 
#and divided by the volume of the rectangle
#value=n_rec/(n*volume(rec))

ans<-obslkm/(n*volume)
return(ans)
}
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









densplit<-function(dendat,minobs=NULL,
leaf=0,method="loglik",splitscan=0,seedf=1,suppo=NULL)
{
n<-dim(dendat)[1]    #havaintojen lkm
d<-dim(dendat)[2]    #muuttujien lkm

if (is.null(minobs)) minobs<-ceiling(sqrt(n)/2)
if (is.null(suppo)) suppo<-supp(dendat)

indendat<-matrix(0,n*d+1,1)
for (i in 1:n){
   for (j in 1:d){
       indendat[1+(i-1)*d+j]=dendat[i,j]
   }
}

insuppo<-matrix(0,2*d+1,1)
insuppo[2:(2*d+1)]<-suppo

#step<-stepcalc(suppo,c(n,n))      n+1 <-> n
step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(suppo[2*i]-suppo[2*i-1])/(n+1)
}
instep<-matrix(0,2*d+1,1)
instep[2:(2*d+1)]<-step

if (method=="loglik") inmethod<-1 else inmethod<-2

suppvol<-massone(suppo)
minvolume<-suppvol/(n+1)^d

maxnodnumR<-100000

ds<-.C("densplitC",as.double(indendat),
                   as.integer(leaf),
                   as.integer(minobs),
                   as.double(insuppo),
                   as.integer(inmethod),
                   as.integer(splitscan),
                   as.integer(seedf),
                   as.integer(n),
                   as.integer(d),
                   as.double(suppvol),
                   as.double(minvolume),
                   as.double(instep),
                   val = integer(maxnodnumR+1),
                   vec = integer(maxnodnumR+1),
                   mean = double(maxnodnumR+1),
                   nelem = integer(maxnodnumR+1),
                   ssr = double(maxnodnumR+1),
                   volume =  double(maxnodnumR+1),
                   left = integer(maxnodnumR+1),
                   right = integer(maxnodnumR+1),
                   glow = integer(d*maxnodnumR+1),
                   gupp = integer(d*maxnodnumR+1),
                   nodenum = integer(1),
                   obspointout = integer(n+1),
                   obslow = integer(maxnodnumR+1),
                   obsupp = integer(maxnodnumR+1))

nodenum<-ds$nodenum
gval<-ds$val[2:(nodenum+1)]
vec<-ds$vec[2:(nodenum+1)]
mean<-ds$mean[2:(nodenum+1)]
nelem<-ds$nelem[2:(nodenum+1)]
ssr<-ds$ssr[2:(nodenum+1)]
volume<-ds$volume[2:(nodenum+1)]
left<-ds$left[2:(nodenum+1)]
right<-ds$right[2:(nodenum+1)]

obspoint<-ds$obspointout[2:(n+1)]
obslow<-ds$obslow[2:(nodenum+1)]
obsupp<-ds$obsupp[2:(nodenum+1)]

val<-matrix(0,nodenum,1)
for (i in 1:nodenum){
  if (vec[i]!=0) val[i]<-suppo[2*vec[i]-1]+gval[i]*step[vec[i]]
}
 
low<-matrix(0,nodenum,d)
upp<-matrix(0,nodenum,d)
for (i in 1:nodenum){
  for (j in 1:d){
      low[i,j]<-ds$glow[1+(i-1)*d+j]
      upp[i,j]<-ds$gupp[1+(i-1)*d+j]
  }
}

#low<-matrix(0,nodenum,d)
#upp<-matrix(0,nodenum,d)
#for (i in 1:nodenum){
#  for (j in 1:d){
#      low[i,j]<-suppo[2*j-1]+ds$glow[1+(i-1)*d+j]*step[j]
#      upp[i,j]<-suppo[2*j-1]+ds$gupp[1+(i-1)*d+j]*step[j]
#  }
#}

puu<-list(split=gval,    #gval=gval,val=t(val),
direc=vec,mean=mean,nelem=nelem,ssr=ssr,volume=volume,
left=left,right=right,low=low,upp=upp,
N=rep(n,d),support=suppo,step=step,
obspoint=obspoint,obslow=obslow,obsupp=obsupp,
minlkm=minobs)
return(puu)
}








densplitR<-function(x,minlkm,suppo,method="loglik",blokki=50)
{
#Makes a density tree.

#x is n*d data-matrix
#minlkm integer >=1
#  Kun solmun elem lkm on minlkm tai pienempi, sita ei enaa jaeta
#suppo on 2*xlkm-vektori, ilmoitaa kullekin muuttujalle kantajan

#tulos on lista(val,vec,mean,nelem,ssr,volum)
#val,vec,mean,nelem,ssr ovat vektoreita joiden pituus on puun silmujen lkm,
#vektorit sis. tiedot puun silmuista.
# val on kohta jossa muuttuja on puolitettu
# vec on puolitetun muuttujan symboli (rivin numero hila:ssa)
# mean on tiheysfunktion estimoitu arvo (silmun havaintojen normeerattu lkm)
# nelem on silmuun kuuluvien havaintojen lkm
# ssr on silmun hav lkm kertaa logaritmi estimaatin arvosta silmussa
# left
# right

#minlkm ohella voitaisiin kaytaa kynnysarvoa minssr;
#silmuja jaetaan siihen asti etta kaikkien lehtien ssr on yli minssr, mutta 
#seuraava jako johtaisi ssr:aan joka on pienempi tai yhtasuuri kuin minssr.
#Xploressa: opt=list(minsize,mincut,mindev)

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
     maxnoden<-maxnoden+blokki
  }
  
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
  #  jatketaan vas. alipuuhun
  while ((nelem[curin]>minlkm) && (volume[curin]>=minvolume)){
  
    # lisaa varmempi testi ks densplitF ??????
 
    #  koska solmu jaettava, tehdaan jako
    
    jako<-findsplit(x,currec,curbeg,curend,obspoint,suppo,n,method)  
    
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
puu<-list(val=val,vec=vec,mean=mean,nelem=nelem,ssr=ssr,volume=volume,
left=left,right=right)
return(puu)
}







densplitter22<-function(dendat,minobs=1,method="loglik")
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

suppo<-matrix(0,2*d,1)  
indet<-matrix(0,2*d,1)  
for (i in 1:d){
    suppo[2*i-1]<-min(dendat[,i])   
    suppo[2*i]<-max(dendat[,i])
    indet[2*i-1]<-seq(1,n)[(suppo[2*i-1]==dendat[,i])]
    indet[2*i]<-seq(1,n)[(suppo[2*i]==dendat[,i])]
}
notindet<-setdiff(seq(1,n),indet)
inside<-dendat[notindet,]
m<-dim(inside)[1]

x<-matrix(0,m+1,d)  # x contains split points and boundaries
for (i in 1:d){
    ordi<-order(inside[,i])
    x[1,i]<-min(dendat[,i])   
    x[m+1,i]<-max(dendat[,i])
    ala<-inside[ordi[1:(m-1)],i]
    yla<-inside[ordi[2:m],i]
    x[2:m,i]<-(ala+yla)/2
}

maksi<-2*n  #n^2
currecs<-matrix(0,maksi,2*d)
for (i in 1:d){
   currecs[1,2*i-1]<-1
   currecs[1,2*i]<-m+1
}
pinin<-1
finrecs<-matrix(0,maksi,2*d)
saatu<-0

while (pinin>0){
   rec<-currecs[pinin,]
   pinin<-pinin-1
   
   fs<-findsplitter(x,rec,n,method)     
   direc<-fs$vec
   point<-fs$val
   leftrec<-rec
   leftrec[2*direc]<-point
   rightrec<-rec
   rightrec[2*direc-1]<-point 

   leftdiffes<-matrix(0,d,1)
   rightdiffes<-matrix(0,d,1)
   for (dd in 1:d){
       leftdiffes[dd]<-leftrec[2*dd]-leftrec[2*dd-1]
       rightdiffes[dd]<-rightrec[2*dd]-rightrec[2*dd-1]
   }
   lkmleft<-min(leftdiffes)
   lkmright<-min(rightdiffes)
   if ((lkmleft>minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-leftrec
      currecs[pinin+2,]<-rightrec
      pinin<-pinin+2
   }
   if ((lkmleft>minobs)&&(lkmright<=minobs)){
      currecs[pinin+1,]<-leftrec
      pinin<-pinin+1
      finrecs[saatu+1,]<-rightrec
      saatu<-saatu+1
   }
   if ((lkmleft<=minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-rightrec
      pinin<-pinin+1
      finrecs[saatu+1,]<-leftrec
      saatu<-saatu+1
   }
}

recs<-finrecs[1:saatu,]
if (saatu==1) recs<-matrix(recs,1,2*d)
truerecs<-recs
for (dd in 1:d){ 
    truerecs[,2*dd-1]<-x[recs[,2*dd-1],dd]
    truerecs[,2*dd]<-x[recs[,2*dd],dd]
}

return(list(recs=truerecs,support=suppo))
}





densplitter33<-function(dendat,minobs=1,method="loglik")
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

suppo<-matrix(0,2*d,1)  
indet<-matrix(0,2*d,1)  
for (i in 1:d){
    suppo[2*i-1]<-min(dendat[,i])   
    suppo[2*i]<-max(dendat[,i])
    indet[2*i-1]<-seq(1,n)[(suppo[2*i-1]==dendat[,i])]
    indet[2*i]<-seq(1,n)[(suppo[2*i]==dendat[,i])]
}
notindet<-setdiff(seq(1,n),indet)
inside<-dendat[notindet,]
m<-dim(inside)[1]

x<-matrix(0,m+1,d)  # x contains split points and boundaries
for (i in 1:d){
    ordi<-order(inside[,i])
    x[1,i]<-min(dendat[,i])   
    x[m+1,i]<-max(dendat[,i])
    ala<-inside[ordi[1:(m-1)],i]
    yla<-inside[ordi[2:m],i]
    x[2:m,i]<-(ala+yla)/2
}

maksi<-2*n  #n^2
currecs<-matrix(0,maksi,2*d)
for (i in 1:d){
   currecs[1,2*i-1]<-1
   currecs[1,2*i]<-m+1
}
pinin<-1
finrecs<-matrix(0,maksi,2*d)
saatu<-0

while (pinin>0){
   rec<-currecs[pinin,]
   pinin<-pinin-1
   
   fs<-findsplitter(x,rec,n,method)     
   direc<-fs$vec
   point<-fs$val
   leftrec<-rec
   leftrec[2*direc]<-point
   rightrec<-rec
   rightrec[2*direc-1]<-point 

   leftdiffes<-matrix(0,d,1)
   rightdiffes<-matrix(0,d,1)
   for (dd in 1:d){
       leftdiffes[dd]<-leftrec[2*dd]-leftrec[2*dd-1]
       rightdiffes[dd]<-rightrec[2*dd]-rightrec[2*dd-1]
   }
   lkmleft<-min(leftdiffes)
   lkmright<-min(rightdiffes)
   if ((lkmleft>minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-leftrec
      currecs[pinin+2,]<-rightrec
      pinin<-pinin+2
   }
   if ((lkmleft>minobs)&&(lkmright<=minobs)){
      currecs[pinin+1,]<-leftrec
      pinin<-pinin+1
      finrecs[saatu+1,]<-rightrec
      saatu<-saatu+1
   }
   if ((lkmleft<=minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-rightrec
      pinin<-pinin+1
      finrecs[saatu+1,]<-leftrec
      saatu<-saatu+1
   }
}

recs<-finrecs[1:saatu,]
if (saatu==1) recs<-matrix(recs,1,2*d)
truerecs<-recs
for (dd in 1:d){ 
    truerecs[,2*dd-1]<-x[recs[,2*dd-1],dd]
    truerecs[,2*dd]<-x[recs[,2*dd],dd]
}

return(list(recs=truerecs,support=suppo))
}





densplitter44<-function(dendat,minobs=1,method="loglik")
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

suppo<-matrix(0,2*d,1)  
indet<-matrix(0,2*d,1)  
for (i in 1:d){
    suppo[2*i-1]<-min(dendat[,i])   
    suppo[2*i]<-max(dendat[,i])
    indet[2*i-1]<-seq(1,n)[(suppo[2*i-1]==dendat[,i])]
    indet[2*i]<-seq(1,n)[(suppo[2*i]==dendat[,i])]
}
notindet<-setdiff(seq(1,n),indet)
inside<-dendat[notindet,]
m<-dim(inside)[1]

maksi<-2*n  #n^2
currecs<-matrix(0,maksi,2*d)
currecs[1,]<-suppo
pinin<-1
finrecs<-matrix(0,maksi,2*d)
saatu<-0
curobs<-matrix(FALSE,maksi,m)
curobs[1,]<-TRUE
curdown<-matrix(0,maksi,d)
curhigh<-matrix(0,maksi,d)
curdown[1,]<-c(1,1)
curhigh[1,]<-c(m+1,m+1)
findown<-matrix(0,maksi,d)
finhigh<-matrix(0,maksi,d)

while (pinin>0){
   rec<-currecs[pinin,]   
   obs<-curobs[pinin,]
   recdown<-curdown[pinin,]
   rechigh<-curhigh[pinin,]
   pinin<-pinin-1
     
   x<-inside[obs,]
   fs<-findsplitter(x,rec,n,method,recdown)  #,rechigh)     
   direc<-fs$vec
   point<-fs$val
   gridpoint<-fs$valio

   leftobs<-(obs&(inside[,direc]<point))
   rightobs<-(obs&(inside[,direc]>point))

   leftrec<-rec
   leftrec[2*direc]<-point
   rightrec<-rec
   rightrec[2*direc-1]<-point 

   leftdown<-recdown
   lefthigh<-rechigh
   lefthigh[direc]<-gridpoint

   rightdown<-recdown
   righthigh<-rechigh
   rightdown[direc]<-gridpoint

   lkmleft<-sum(leftobs)
   lkmright<-sum(rightobs)
   if ((lkmleft>minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-leftrec
      currecs[pinin+2,]<-rightrec
      curobs[pinin+1,]<-leftobs
      curobs[pinin+2,]<-rightobs
      curdown[pinin+1,]<-leftdown
      curhigh[pinin+1,]<-lefthigh
      curdown[pinin+2,]<-rightdown
      curhigh[pinin+2,]<-righthigh
      pinin<-pinin+2
   }
   if ((lkmleft>minobs)&&(lkmright<=minobs)){
      currecs[pinin+1,]<-leftrec
      curobs[pinin+1,]<-leftobs
      curdown[pinin+1,]<-leftdown
      curhigh[pinin+1,]<-lefthigh
      pinin<-pinin+1
      finrecs[saatu+1,]<-rightrec
      findown[saatu+1,]<-rightdown
      finhigh[saatu+1,]<-righthigh
      saatu<-saatu+1
   }
   if ((lkmleft<=minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-rightrec
      curobs[pinin+1,]<-rightobs
      curdown[pinin+1,]<-rightdown
      curhigh[pinin+1,]<-righthigh
      pinin<-pinin+1
      finrecs[saatu+1,]<-leftrec
      findown[saatu+1,]<-leftdown
      finhigh[saatu+1,]<-lefthigh
      saatu<-saatu+1
   }
   if ((lkmleft<=minobs)&&(lkmright<=minobs)){
      finrecs[saatu+1,]<-leftrec
      finrecs[saatu+2,]<-rightrec
      findown[saatu+1,]<-leftdown
      findown[saatu+2,]<-rightdown
      finhigh[saatu+1,]<-lefthigh
      finhigh[saatu+2,]<-righthigh
      saatu<-saatu+2
   }
}

recs<-finrecs[1:saatu,]
down<-findown[1:saatu,]
high<-finhigh[1:saatu,]
if (saatu==1){
   recs<-matrix(recs,1,2*d)
   down<-matrix(down,1,d)
   high<-matrix(high,1,d)
}

grid<-matrix(0,m+1,d)  # grid contains split points and boundaries
for (i in 1:d){
    ordi<-order(inside[,i])
    grid[1,i]<-min(dendat[,i])   
    grid[m+1,i]<-max(dendat[,i])
    ala<-inside[ordi[1:(m-1)],i]
    yla<-inside[ordi[2:m],i]
    grid[2:m,i]<-(ala+yla)/2
}

return(list(recs=recs,support=suppo,grid=grid,down=down,high=high))
}





densplitter<-function(dendat,minobs=1,method="loglik",dyadic=FALSE)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

suppo<-matrix(0,2*d,1)  
indet<-matrix(0,2*d,1)  
for (i in 1:d){
    suppo[2*i-1]<-min(dendat[,i])   
    suppo[2*i]<-max(dendat[,i])
    indet[2*i-1]<-seq(1,n)[(suppo[2*i-1]==dendat[,i])]
    indet[2*i]<-seq(1,n)[(suppo[2*i]==dendat[,i])]
}
notindet<-setdiff(seq(1,n),indet)
inside<-dendat[notindet,]
m<-dim(inside)[1]

grid<-matrix(0,m+1,d)  # grid contains split points and boundaries
for (i in 1:d){
    ordi<-order(inside[,i])
    grid[1,i]<-min(dendat[,i])   
    grid[m+1,i]<-max(dendat[,i])
    ala<-inside[ordi[1:(m-1)],i]
    yla<-inside[ordi[2:m],i]
    grid[2:m,i]<-(ala+yla)/2
}

maksi<-2*n  #n^2
currecs<-matrix(0,maksi,2*d)
currecs[1,]<-suppo
pinin<-1
finrecs<-matrix(0,maksi,2*d)
saatu<-0
curobs<-matrix(FALSE,maksi,m)
curobs[1,]<-TRUE
curdown<-matrix(0,maksi,d)
curhigh<-matrix(0,maksi,d)
curdown[1,]<-rep(1,d)
curhigh[1,]<-rep(m+1,d)
findown<-matrix(0,maksi,d)
finhigh<-matrix(0,maksi,d)

while (pinin>0){
   rec<-currecs[pinin,]   
   obs<-curobs[pinin,]
   recdown<-curdown[pinin,]
   rechigh<-curhigh[pinin,]
   pinin<-pinin-1
     
   x<-inside[obs,]
   if (dyadic) fs<-findsplitter.dyadic(grid,x,rec,n,method,minobs,recdown,rechigh)     
   else fs<-findsplitter(grid,x,rec,n,method,minobs,recdown,rechigh)     
   direc<-fs$vec
   point<-fs$val
   gridpoint<-fs$valio

   leftobs<-(obs&(inside[,direc]<point))
   rightobs<-(obs&(inside[,direc]>point))

   leftrec<-rec
   leftrec[2*direc]<-point
   rightrec<-rec
   rightrec[2*direc-1]<-point 

   leftdown<-recdown
   lefthigh<-rechigh
   lefthigh[direc]<-gridpoint

   rightdown<-recdown
   righthigh<-rechigh
   rightdown[direc]<-gridpoint

   lkmleft<-sum(leftobs)
   lkmright<-sum(rightobs)
   if ((lkmleft>minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-leftrec
      currecs[pinin+2,]<-rightrec
      curobs[pinin+1,]<-leftobs
      curobs[pinin+2,]<-rightobs
      curdown[pinin+1,]<-leftdown
      curhigh[pinin+1,]<-lefthigh
      curdown[pinin+2,]<-rightdown
      curhigh[pinin+2,]<-righthigh
      pinin<-pinin+2
   }
   if ((lkmleft>minobs)&&(lkmright<=minobs)){
      currecs[pinin+1,]<-leftrec
      curobs[pinin+1,]<-leftobs
      curdown[pinin+1,]<-leftdown
      curhigh[pinin+1,]<-lefthigh
      pinin<-pinin+1
      finrecs[saatu+1,]<-rightrec
      findown[saatu+1,]<-rightdown
      finhigh[saatu+1,]<-righthigh
      saatu<-saatu+1
   }
   if ((lkmleft<=minobs)&&(lkmright>minobs)){
      currecs[pinin+1,]<-rightrec
      curobs[pinin+1,]<-rightobs
      curdown[pinin+1,]<-rightdown
      curhigh[pinin+1,]<-righthigh
      pinin<-pinin+1
      finrecs[saatu+1,]<-leftrec
      findown[saatu+1,]<-leftdown
      finhigh[saatu+1,]<-lefthigh
      saatu<-saatu+1
   }
   if ((lkmleft<=minobs)&&(lkmright<=minobs)){
      finrecs[saatu+1,]<-leftrec
      finrecs[saatu+2,]<-rightrec
      findown[saatu+1,]<-leftdown
      findown[saatu+2,]<-rightdown
      finhigh[saatu+1,]<-lefthigh
      finhigh[saatu+2,]<-righthigh
      saatu<-saatu+2
   }
}

recs<-finrecs[1:saatu,]
down<-findown[1:saatu,]
high<-finhigh[1:saatu,]
if (saatu==1){
   recs<-matrix(recs,1,2*d)
   down<-matrix(down,1,d)
   high<-matrix(high,1,d)
}

return(list(recs=recs,support=suppo,grid=grid,down=down,high=high))
}





denssr<-function(volume,havlkm,n,method="loglik",mix=NULL){
#Calculates the log-likelihood 
#
if (havlkm==0) vastaus<-0
else if (method=="loglik")
                      vastaus<-havlkm*log(havlkm/(n*volume))
else if (method=="projec"){
   if (is.null(mix))  vastaus<-havlkm^2/(n*volume)
   else vastaus<-mix*(2-mix)*havlkm^2/(n^2*volume)
}
#
return(vastaus)
}



dentree<-function(treeseq,dendat,leafnum=NULL){
#Returns a tree from a sequence of trees
#
#if leafnum is NULL, we use empirical smoothing parameter selection
#
n<-length(dendat[,1])
#poistettavat<-treeseq$leafs
#if (dim(t(poistettavat))[1]==1) maxrem<-1
#else maxrem<-length(poistettavat[,1])      
#
if (is.null(leafnum)){   #empirical smoothing param. selection
  ri<-riskesti(treeseq,n)
  indeksi<-omaind(ri)
#  indeksit<-detsikko(treeseq$leafs,indeksi) #Haet haluttu puu puit jonosta
#  alipuu<-poistamon(treeseq$tree,indeksit)
}
else{
  indeksi<-detsi(treeseq$info[,1],leafnum)
#  indeksit<-detsikko(treeseq$leafs,indeksi)
  #indeksit<-t(t(treeseq$leafs))[indeksi,] 
#  alipuu<-poistamon(treeseq$tree,indeksit)  #alipuu<-dpoimi(puuseq,indeksi)
}
#
# Tehdaan binaaripuusta paloittain vakio
#
xlkm<-length(dendat[1,])
epsi<-0.01
#kanta<-kantaja(dendat,epsi)
#palvak<-teeositus(alipuu,xlkm,kanta)
#
# Tehdaan paloittain vakiosta tiheyspuu
#
#return(list(binpuu=alipuu,palvak=palvak))
alipuu<-NULL
return(alipuu)
}



dernor<-function(x)
{
ans<--x*evanor(x)
return(ans)
}

detsi<-function(pysty,lehtilkm){
#Hakee pysty-vektorista "pysty" sen indeksin, jonka kohdalla
#esiintyy luku "lehtilkm"
#
lkm<-length(pysty)
i<-1
while ((pysty[i]!=lehtilkm) && (i<=lkm)) i<-i+1
return(i)
}
downhigh<-function(et)
{
leafnum<-length(et$left)
d<-length(et$N)
value<-matrix(0,leafnum,1)
down<-matrix(0,leafnum,d)
high<-matrix(0,leafnum,d)
infopointer<-matrix(0,leafnum,1)

leafloc<-findleafs(et$left,et$right)

efek<-0
i<-1
while (i<=leafnum){  

   if (!is.na(leafloc[i])) if (leafloc[i]==1){   #if (mean[node]>0){
       efek<-efek+1

       infopointer[i]<-efek
       value[efek]<-et$mean[i]
 
       for (j in 1:d){
           down[efek,j]<-et$low[i,j]
           high[efek,j]<-et$upp[i,j]
       }
   }
   i<-i+1
}

value<-value[1:efek]

if (efek>1){
   down<-down[1:efek,]
   high<-high[1:efek,]
}
else{
   apudown<-matrix(0,1,d)
   apuhigh<-matrix(0,1,d)
   for (ddd in 1:d){
        apudown[1,ddd]<-down[1,ddd]
        apuhigh[1,ddd]<-high[1,ddd]
   }
   down<-apudown
   high<-apuhigh
}

return(list(down=down,high=high,value=value,infopointer=infopointer))
}

drawcart<-function(treeseq,lehtilkm,suppo,plkm){
#Makes data for drawing a perspective plot.
#
#plkm on kuvaajan hilan pisteiden lkm
#
#koe<-drawcart(dendat,lehtilkm=3,treeseq,epsi=0,plkm=30)
#persp(koe$x,koe$y,koe$z,phi=30,theta=60) 
#
# Haetaan ensin haluttu puu puitten jonosta
tree<-treeseq$tree
delnodes<-treeseq$delnodes
delend<-treeseq$delend
leafs<-treeseq$leafs
indeksi<-detsi(leafs,lehtilkm)
endi<-delend[indeksi]
if (endi>0){  #if there is something to remove
  indeksit<-delnodes[1:endi]
  re<-remnodes(tree$left,tree$right,indeksit)  
  tree$left<-re$left
  tree$right<-re$right
}
# Tehdaan binaaripuusta paloittain vakio
pv<-partition(tree,suppo)       
recs<-pv$recs
values<-pv$values
#
ans<-drawgene(values,recs,plkm)
return(list(x=ans$x,y=ans$y,z=ans$z))
}










drawgene<-function(values,recs,plkm=60,ep1=0.5){
#Makes data for drawing a perspective plot.

#plkm on kuvaajan hilan pisteiden lkm
#ep1 makes 0-corona around the support (useful for densities)

#koe<-drawgene(values,recs,plkm=30)
#persp(koe$x,koe$y,koe$z,phi=30,theta=60) 

alkux<-min(recs[,1])-ep1
alkuy<-min(recs[,3])-ep1
loppux<-max(recs[,2])+ep1
loppuy<-max(recs[,4])+ep1
pitx<-(loppux-alkux)/plkm
pity<-(loppuy-alkuy)/plkm
x<-alkux+c(0:plkm)*pitx
y<-alkuy+c(0:plkm)*pity

reclkm<-length(values)
xdim<-length(x)
ydim<-length(y)
arvot<-matrix(0,xdim,ydim)

l<-1
while (l<=reclkm){
   begx<-recs[l,1]
   endx<-recs[l,2]
   begy<-recs[l,3]
   endy<-recs[l,4]

   begxind<-round(plkm*(begx-alkux)/(loppux-alkux))
   endxind<-round(plkm*(endx-alkux)/(loppux-alkux))
   begyind<-round(plkm*(begy-alkuy)/(loppuy-alkuy))
   endyind<-round(plkm*(endy-alkuy)/(loppuy-alkuy))

   arvot[begxind:endxind,begyind:endyind]<-values[l]

   l<-l+1
}

return(list(x=x,y=y,z=arvot))
#persp(x,y,arvot)
}










eval.bagg<-function(dendat,B,leaf,minobs=NULL,seed=1,
sample="bagg",prune="off",
splitscan=0,seedf=1,
scatter=0,
src="c",method="loglik")
{
#B number of bootstrap samples
#leaf number of leaves in the trees to be grown

n<-dim(dendat)[1]
d<-dim(dendat)[2]
if (is.null(minobs)) minobs<-ceiling(sqrt(n)/2)

suppo<-supp(dendat)+scatter*rep(c(-1,1),d)
step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(suppo[2*i]-suppo[2*i-1])/(n+1)
}
{
if (sample=="bagg"){
  dendatB<-bootbagg(dendat,seed,scatter) 
}
else{  # scheme=="baggworpl"
  dendatB<-bootworpl(dendat,seed,scatter)
}
}
{
if (prune=="off"){
   if (leaf==0){
        tr<-densplit(dendatB,minobs,leaf=0,
        method=method,splitscan=splitscan,seedf=seedf,suppo=suppo)
   }
   else{
        tr<-eval.greedy(dendatB,leaf,method=method,suppo=suppo)
   }
}
else{  #prune == on
  if (src=="c"){
    tr<-densplit(dendatB,minobs,leaf=0,
        method=method,splitscan=splitscan,seedf=seedf,suppo=suppo)
    treeseq<-prune(tr)
    approleaf<-roundlnum(treeseq$leafs,leaf)
    tr<-eval.pick(treeseq,approleaf)
  }
  else{
    tr<-densplitF(dendatB,0,minobs,suppo)
    treeseq<-prune(tr)
    approleaf<-roundlnum(treeseq$leafs,leaf)
    tr<-eval.pick(treeseq,approleaf)
  }
}
}

bi<-2
while (bi<=B){
   {
   if (sample=="bagg"){
      dendatB<-bootbagg(dendat,seed+bi-1)  
   }
   else{
      dendatB<-bootworpl(dendat,seed+bi-1)
   }
   }
   if (prune=="off"){
       if (leaf==0){
          trnew<-densplit(dendatB,minobs,leaf=0,
                 method=method,splitscan=splitscan,seedf=seedf,suppo=suppo)
       }
       else{
          trnew<-eval.greedy(dendatB,leaf,method=method,suppo=suppo)
       }
   }
   else{  #prune == on
      if (src=="c"){
        trnew<- densplit(dendatB,minobs,leaf=0,
                method=method,splitscan=splitscan,seedf=seedf,suppo=suppo)
      treeseq<-prune(trnew)
      approleaf<-roundlnum(treeseq$leafs,leaf)
      trnew<-eval.pick(treeseq,approleaf)
      }
      else{
      trnew<-densplitF(dendatB,0,minobs,suppo,splitscan=splitscan,seedf=seedf)
      treeseq<-prune(trnew)
      approleaf<-roundlnum(treeseq$leafs,leaf)
      trnew<-eval.pick(treeseq,approleaf)
      }
   }

   tr<-treeadd(tr,trnew,bi-1)
   bi<-bi+1
   
}

tr<-c(tr,list(support=suppo,step=step))
tr$N<-rep(dim(dendatB)[1],d)

tayd<-partigen.disc(tr)
tr$value<-tayd$value
tr$down<-tayd$down
tr$high<-tayd$high

return(tr)
}







eval.cart<-function(dendat,leaf,minobs=NULL)
{
n<-dim(dendat)[1]
if (is.null(minobs)) minobs<-ceiling(sqrt(n)/2)

bt<-densplit(dendat,minobs)
treeseq<-prune(bt)
aplnum<-roundlnum(treeseq$leafs,leaf)
eval<-eval.pick(treeseq,aplnum)  
eval$lnum<-aplnum

return(eval)
}



eval.greedy<-function(dendat,leaf,method="loglik",minobs=NULL,bound=0,
suppo=NULL)
{

# splitinte<-floor(log(leaf,base=2))
# make first splitinte splits, then leaf-2^splitinte additional splits

n<-length(dendat[,1])  #havaintojen lkm
d<-length(dendat[1,])  #muuttujien lkm

if (is.null(minobs)) minobs<-ceiling(sqrt(n)/2)

if (is.null(suppo)) suppo<-supp(dendat,blown=TRUE)


suppvol<-massone(suppo)

maxnoden<-2*leaf
val<-matrix(0,maxnoden,1) 
vec<-matrix(0,maxnoden,1)
mean<-matrix(0,maxnoden,1)
loglik<-matrix(0,maxnoden,1)
nelem<-matrix(0,maxnoden,1)
volume<-matrix(0,maxnoden,1)
left<-matrix(0,maxnoden,1)
right<-matrix(0,maxnoden,1)
low<-matrix(0,maxnoden,d)
upp<-matrix(0,maxnoden,d)

recs<-matrix(0,maxnoden,2*d)
begs<-matrix(0,maxnoden,1)
ends<-matrix(0,maxnoden,1)

step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(suppo[2*i]-suppo[2*i-1])/n
}

# for the C-function 

xdendat<-matrix(0,d*n+1,1)
obsoso<-matrix(0,n+1,1)
insuppo<-matrix(0,2*d+1,1)
insuppo[2:(2*d+1)]<-suppo
instep<-matrix(0,d+1,1)
instep[2:(d+1)]<-step
ingrec<-matrix(0,2*d+1,1)
if (method=="loglik") inmethod<-1 else inmethod<-2

# info for root node

val[1]<-0
vec[1]<-0
curvol<-massone(suppo)
mean[1]<-denmean(curvol,n,n)
loglik[1]<-denssr(curvol,n,n,method)
nelem[1]<-n
volume[1]<-curvol
left[1]<-0
right[1]<-0
for (k in 1:d){
  low[1,k]<-0           #suppo[2*k-1]
  upp[1,k]<-n           #suppo[2*k]
  recs[1,2*k-1]<-0      #suppo[2*k-1]
  recs[1,2*k]<-n        #suppo[2*k]
}
begs[1]<-1
ends[1]<-n

# initialize

obspoint<-seq(1:n)
curleafnum<-1
curnodenum<-1

while (curleafnum<leaf){

  locs<-leaflocs(left,right)
  curleafnum<-locs$leafnum

  incre<-matrix(0,curleafnum,1)  
     #for each leaf find the increase in loglik  
     #we choose to make split which increases most the total loglik
  valpool<-matrix(0,curleafnum,1)
  vecpool<-matrix(0,curleafnum,1)
  leftrecpool<-matrix(0,curleafnum,2*d)
  rightrecpool<-matrix(0,curleafnum,2*d)
  leftbegpool<-matrix(0,curleafnum,1)
  leftendpool<-matrix(0,curleafnum,1)
  rightbegpool<-matrix(0,curleafnum,1)
  rightendpool<-matrix(0,curleafnum,1)
  obspointpool<-matrix(0,curleafnum,n)

  failnum<-0  # count the number of nodes where we are not able to make split

  for (i in 1:curleafnum){

     loca<-locs$leafloc[i]
     currec<-recs[loca,]
     curbeg<-begs[loca]
     curend<-ends[loca]
     curloglik<-loglik[loca]
     {
     if ((curbeg==0) || (curend==0)) maara<-0
     else maara<-count(curbeg,curend)
     }

     smallest<-1       #smallest==1, when "currec" is the minimal bin
     for (j in 1:d){
         valli<-currec[2*j]-currec[2*j-1]
         if (valli>1) smallest<-0
     }   

     #volli<-massone(currec)*prod(step)
     #suhde<-volli/suppvol

if ((maara<=minobs) || (smallest==1)){  # || (suhde<bound)){
        incre[i]<-NA
        failnum<-failnum+1
}

else{

     for (j in curbeg:curend){
           obspointer=obspoint[j]
           obsoso[1+j-curbeg+1]=obspointer
           jin=j-curbeg+1
           for (l in 1:d){
               xdendat[1+(jin-1)*d+l]<-dendat[obspointer,l]
           }
     }

     ingrec[2:(2*d+1)]<-currec
 
jako<-.C("findsplitCC",as.double(xdendat),
                       as.integer(obsoso),
                       as.integer(maara),
                       as.integer(curbeg),
                       as.integer(curend),
                       as.double(insuppo),
                       as.double(instep),
                       as.integer(ingrec),
                       #as.double(inrec),
                       as.integer(n),
                       as.integer(d),
                       as.integer(inmethod),
                       as.double(bound),
                       val = integer(1),
                       vec = integer(1),
                       #leftrec = integer(2*d+1),
                       #rightrec = integer(2*d+1),
                       leftbeg = integer(1),
                       leftend = integer(1),
                       rightbeg = integer(1),
                       rightend = integer(1),
                       obspoint = integer(maara+1))

     vecu<-jako$vec
     valu<-jako$val   #suppo[2*vecu-1]+jako$val*step[vecu] 
     leftrec<-currec
     leftrec[2*vecu]<-valu
     rightrec<-currec
     rightrec[2*vecu-1]<-valu   
     #leftrec<-jako$leftrec[2:(2*d+1)]
     #rightrec<-jako$rightrec[2:(2*d+1)]

     leftbeg<-jako$leftbeg
     leftend<-jako$leftend
     rightbeg<-jako$rightbeg
     rightend<-jako$rightend
     for (li in 1:maara){
        obspoint[curbeg+li-1]<-jako$obspoint[li+1]
     }

     #lvolume<-massone(leftrec)
     #rvolume<-massone(rightrec)
     lvolume<-1
     rvolume<-1
     for (ji in 1:d){
          lvolume<-lvolume*(leftrec[2*ji]-leftrec[2*ji-1])*step[ji]
          rvolume<-rvolume*(rightrec[2*ji]-rightrec[2*ji-1])*step[ji]
     }

     lnelem<-count(leftbeg,leftend)   #leftend-leftbeg+1
     rnelem<-count(rightbeg,rightend)   #rightend-rightbeg+1
     newloglik<-denssr(lvolume,lnelem,n,method)+
                denssr(rvolume,rnelem,n,method)
     incre[i]<-newloglik-curloglik

     if ((lvolume/suppvol < bound) && (lnelem>0)) incre[i]<-NA
     if ((rvolume/suppvol < bound) && (rnelem>0)) incre[i]<-NA
 
     valpool[i]<-valu
     vecpool[i]<-vecu
     leftrecpool[i,]<-leftrec
     rightrecpool[i,]<-rightrec
     leftbegpool[i]<-leftbeg
     leftendpool[i]<-leftend
     rightbegpool[i]<-rightbeg
     rightendpool[i]<-rightend
     obspointpool[i,]<-obspoint

} #else (we may split because there are observations in the rec)

  }   #for (i in 1:curleafnum){

####################################################

# now "incre" is the optimization vector, we find the
# maximum of this vector: this gives the best split

allfail<-1
for (ii in 1:d){
  if (!is.na(incre[ii]))  allfail<-0
}

sd<-omaind(-incre)   #omaind minimizes, we want to maximize

if ((failnum==curleafnum)){  # || (allfail==1)){  
             # we have to finish (no potential splits exist)

cl<-curnodenum

return(list(split=val[1:cl],direc=vec[1:cl],mean=mean[1:cl],nelem=nelem[1:cl],
ssr=loglik[1:cl],volume=volume[1:cl],
left=left[1:cl],right=right[1:cl],low=low[1:cl,],upp=upp[1:cl,],
suppo=suppo,step=step))

}

else{

  # make the split sd

  locloc<-locs$leafloc[sd]
  val[locloc]<-valpool[sd] 
  vec[locloc]<-vecpool[sd]

  # create left child

  leftpoint<-curnodenum+1
  left[locloc]<-leftpoint

  recu<-leftrecpool[sd,]
  #volu<-massone(recu)
  volu<-1
  for (ji in 1:d){
     volu<-volu*(recu[2*ji]-recu[2*ji-1])*step[ji]
  }
  nelemu<-count(leftbegpool[sd],leftendpool[sd])  
        #leftendpool[sd]-leftbegpool[sd]+1

  val[leftpoint]<-0
  vec[leftpoint]<-0
  mean[leftpoint]<-denmean(volu,nelemu,n)
  loglik[leftpoint]<-denssr(volu,nelemu,n,method)
  nelem[leftpoint]<-nelemu
  volume[leftpoint]<-volu
  for (k in 1:d){
    low[leftpoint,k]<-recu[2*k-1]
    upp[leftpoint,k]<-recu[2*k]
  }
  upp[leftpoint,vec[locloc]]<-val[locloc]

  recs[leftpoint,]<-recu
  begs[leftpoint]<-leftbegpool[sd]
  ends[leftpoint]<-leftendpool[sd]

  # create right child

  rightpoint<-curnodenum+2
  right[locloc]<-rightpoint

  recu<-rightrecpool[sd,] 
  #volu<-massone(recu)
  volu<-1
  for (ji in 1:d){
     volu<-volu*(recu[2*ji]-recu[2*ji-1])*step[ji]
  }
  nelemu<-count(rightbegpool[sd],rightendpool[sd]) 
          #rightendpool[sd]-rightbegpool[sd]+1

  val[rightpoint]<-0
  vec[rightpoint]<-0
  mean[rightpoint]<-denmean(volu,nelemu,n)
  loglik[rightpoint]<-denssr(volu,nelemu,n,method)
  nelem[rightpoint]<-nelemu
  volume[rightpoint]<-volu
  for (k in 1:d){
    low[rightpoint,k]<-recu[2*k-1]
    upp[rightpoint,k]<-recu[2*k]
  }
  low[rightpoint,vec[locloc]]<-val[locloc]

  recs[rightpoint,]<-recu
  begs[rightpoint]<-rightbegpool[sd]
  ends[rightpoint]<-rightendpool[sd]

  # final updates

  curleafnum<-curleafnum+1
  curnodenum<-curnodenum+2
  obspoint<-obspointpool[sd,]

}  #end split making


}  #while (curleafnum<leaf)

cl<-curnodenum

##################################################
ll<-leaflocs(left[1:cl],right[1:cl])
leafloc<-ll$leafloc
leafnum<-ll$leafnum

value<-matrix(0,leafnum,1)
down<-matrix(0,leafnum,d)
high<-matrix(0,leafnum,d)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if (mean[node]>0){
     efek<-efek+1

     value[efek]<-mean[node]
 
     for (j in 1:d){
         down[efek,j]<-low[node,j]
         high[efek,j]<-upp[node,j]
     }
   }
   i<-i+1
}
value<-value[1:efek]
if (efek>1){
   down<-down[1:efek,]
   high<-high[1:efek,]
}
else{
   apudown<-matrix(0,1,d)
   apuhigh<-matrix(0,1,d)
   for (ddd in 1:d){
        apudown[1,ddd]<-down[1,ddd]
        apuhigh[1,ddd]<-high[1,ddd]
   }
   down<-apudown
   high<-apuhigh
}

###################################################

return(list(split=val[1:cl],direc=vec[1:cl],mean=mean[1:cl],nelem=nelem[1:cl],
ssr=loglik[1:cl],volume=volume[1:cl],
left=left[1:cl],right=right[1:cl],low=low[1:cl,],upp=upp[1:cl,],
support=suppo,step=step,
value=value,down=down,high=high,N=rep(n,d)))

}















eval.pick<-function(treeseq,leaf){
#
tree<-treeseq$tree
delnodes<-treeseq$delnodes
delend<-treeseq$delend
leafs<-treeseq$leafs
indeksi<-detsi(leafs,leaf)
endi<-delend[indeksi]
if (endi>0){       #if there is something to remove
  indeksit<-delnodes[1:endi]
  re<-remnodes(tree$left,tree$right,indeksit)
  tree$left<-re$left
  tree$right<-re$right
}  

####################################################
cl<-length(tree$left)
d<-length(tree$N)

down<-matrix(0,cl,d)
high<-matrix(0,cl,d)
for (i in 1:cl){
   for (j in 1:d){
      down[i,j]<-tree$low[i,j]+1
      high[i,j]<-tree$upp[i,j]
    }
}

ll<-leaflocs(tree$left[1:cl],tree$right[1:cl])
leafloc<-ll$leafloc
leafnum<-ll$leafnum

value<-matrix(0,leafnum,1)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if (tree$mean[node]>0){
     efek<-efek+1

     value[efek]<-tree$mean[node]
 
     for (j in 1:d){
         down[efek,j]<-tree$low[node,j]
         high[efek,j]<-tree$upp[node,j]
     }
   }
   i<-i+1
}
value<-value[1:efek]
if (efek>1){
   down<-down[1:efek,]
   high<-high[1:efek,]
}
else{
   apudown<-matrix(0,1,d)
   apuhigh<-matrix(0,1,d)
   for (ddd in 1:d){
        apudown[1,ddd]<-down[1,ddd]
        apuhigh[1,ddd]<-high[1,ddd]
   }
   down<-apudown
   high<-apuhigh
}

###################################################
tree$value<-value
tree$down<-down
tree$high<-high
   
return(tree)
}
eval.stage.gauss<-function(dendat,M,mugrid,siggrid=1,sigeka=TRUE,src="c",
sampstart=FALSE,boost=FALSE,N=60)
{
n<-length(dendat)

  if (src=="R"){
     resu<-stage.gaussR(dendat,M,mugrid,siggrid,sigeka,sampstart)
     return(resu)
  }
  else{
            if (boost) inboost<-1 else inboost<-0
            mu0<-mean(dendat)
            sig0<-sqrt(var(dendat))

            indendat<-c(0,dendat)
            inM<-M
            inmugrid<-c(0,mugrid)
            insiggrid<-c(0,siggrid)
            insigeka<-1
            if (sampstart) insampstart<-1 else insampstart<-0
            dictCard<-length(mugrid)
            dictCardSig<-length(siggrid)
            kg<-.C("stageGauss",
               as.double(mu0),
               as.double(sig0),
               as.double(indendat),
               as.integer(inM),
               as.double(inmugrid),
               as.double(insiggrid),
               as.integer(insigeka),
               as.integer(insampstart),
               as.integer(n),
               as.integer(dictCard),
               as.integer(dictCardSig),
               as.integer(inboost),
               muut = double(inM+1),
               sigit = double(inM+1),
               curmix = double(inM+1))
            sgmuut<-kg$muut[2:(inM+1)]
            sgsigit<-kg$sigit[2:(inM+1)]
            sgcurmix<-kg$curmix[2:(inM+1)]

            minu<-min(sgmuut)-2*max(sgsigit)
            maxi<-max(sgmuut)+2*max(sgsigit)
            support<-c(minu,maxi)
            pcf<-pcf.func("mixt",N,sig=sgsigit,M=sgmuut,p=sgcurmix,
            support=support)

            return(list(value=pcf$value,down=pcf$down,high=pcf$high,
                        support=pcf$support,N=pcf$N,
                        muut=sgmuut,sigit=sgsigit,curmix=sgcurmix))
  }

}

  

eval.stage<-function(dendat,leaf,M,pis=NULL,mcn=dim(dendat)[1],
minobs=NULL,seedi=1,
method="projec",bound=0)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

if (is.null(minobs)) minobs<-ceiling(sqrt(n)/2)

if (is.null(pis)){
   pis<-matrix(0,(M-1),1)
   for (k in 1:(M-1)) pis[k]<-2/(k+2)
}

tr<-eval.greedy(dendat,leaf,method,minobs)

suppo<-supp(dendat,blown=TRUE)
N<-rep(n,d)
step<-stepcalc(suppo,N)

i<-1
while (i<=(M-1)){

   seedi<-seedi+1
   mcdendat<-simutree(tr,mcn,seedi)

   mix<-pis[i]
   trnew<-myosplitpena(dendat,leaf,mcdendat,mix,suppo,step,minobs,method)
   #trnew<-myosplitpenaR(dendat,leaf,mcdendat,mix,suppo,step,minobs,method)

   tr<-treeadd(tr,trnew,mix=mix)

   i<-i+1
}

##################################################
ll<-leaflocs(tr$left,tr$right)
leafloc<-ll$leafloc
leafnum<-ll$leafnum

value<-matrix(0,leafnum,1)
down<-matrix(0,leafnum,d)
high<-matrix(0,leafnum,d)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if (tr$mean[node]>0){
     efek<-efek+1

     value[efek]<-tr$mean[node]
 
     for (j in 1:d){
         down[efek,j]<-tr$low[node,j]
         high[efek,j]<-tr$upp[node,j]
     }
   }
   i<-i+1
}
tr$value<-value[1:efek]
tr$down<-down[1:efek,]
tr$high<-high[1:efek,]
###################################################

tr$N<-N

return(tr)
}








exma<-function(volume,obslkm,n,lambda)
{

return(-obslkm/n+lambda*volume)

}
exmavec<-function(vol,obsnum,n,lambda)
{

len<-length(vol)
resu<-matrix(0,len,1)

for (i in 1:len){
  
  resu[i]<-exma(vol[i],obsnum[i],n,lambda)

}

return(resu)
}



findleafs<-function(left,right)
{
# Finds location of leafs in binary tree, living in vector
# left, right are itemnum-vectors

# returns itemnum-vector, 1 in the location of nodes 0 non-terminal
# NA in positions not belonging to tree

# vector where binary tree is living may be larger than cardinality
# of nodes of the tree

itemnum<-length(left)
leafloc<-matrix(NA,itemnum,1)
pino<-matrix(0,itemnum,1)
pino[1]<-1     #pino[1]<-root
pinin<-1
while (pinin>0){
    cur<-pino[pinin]      #take from stack
    pinin<-pinin-1
    if (left[cur]==0){    #if we are in leaf
       leafloc[cur]<-1
    }
    else{
       while (left[cur]>0){  #go to leaf and put right nodes to stack
           leafloc[cur]<-0
           pinin<-pinin+1
           pino[pinin]<-right[cur]
           cur<-left[cur]
       }
       leafloc[cur]<-1  #now we know we are in leaf
    }
}
return(leafloc)
} 







findobs<-function(x,rec,ordobs,leftend,coordi){
#Finds which observations belong to given rectangle
#
#x on n*d havaintomatriisi
#rec on d*2-vector, sis kuvauksen uudesta osiosta
#ordobs is lkm-vector, points to the rows of x, ordered with respect
#  to i:th coordinate of x
#leftend >=0, <=lkm, integer, which pointers in ordobs point to the 
#  observations on the left rectangle
#coordi in 1:d, according to this coordinate pointers are ordered 
#
#Returns leftend in 0:lkm, endpoint of those belonging to rec
#
#We can start from leftend and proceed to the right hand side,
#because it was already checked that the observations on the left  
#side do not belong to rec.
#
lkm<-length(ordobs)
#
if (leftend<lkm){ #otherwise all obs already are on the left
   lcount<-0
   for (i in (leftend+1):lkm){
       ind<-ordobs[i]
       obs<-x[ind,]
       if (belongs(obs,rec,coordi))  lcount<-lcount+1  
   }
   leftend<-leftend+lcount
}
return(leftend)
}
findobsR<-function(dendat, ordobs, leftrecend, coordi)
{
  maara<-length(ordobs)
  j<-1
  ind<-ordobs[j]
  havakoor<-dendat[ind,coordi]

  lcount<-0
  while ((havakoor<=leftrecend) && (j<maara)){
          lcount<-lcount+1
          j<-j+1
          ind<-ordobs[j]
          havakoor<-dendat[ind,coordi]
  }
  ind<-ordobs[maara]
  havakoor<-dendat[ind,coordi] 
  if (havakoor<=leftrecend) lcount<-lcount+1

  leftend<-lcount  #/* could be zero */

return(leftend)
}

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
    






gaussprod<-function(mu1,mu2,sig1=1,sig2=1)
{
resp<-(2*pi)^(-1/2)*(sig1^2+sig2^2)^(-1/2)*exp(-(mu1-mu2)^2/(2*(sig1^2+sig2^2)))
return(resp)
}







initial<-function(ssr,left,right){
#Prunes a tree by removing useless splits.
#
#ssr, left, right are nodelkm-vectors
#tree is a list(ssr,left,right,...)
#
#Returnes left,right
#
#densplit saattaa muodostaa puun, jossa jokin alipuu
#on huonompi kuin ko. alipuun juuri. 
#Muokataan puuta siten, etta alipuut ovat vahintaan yhta hyvia.
#Algoritmi sama kuin Donoho 95, ts. 
#lahdetaan yhta tasoa lehtia ylempaa, edetaan taso kerrallan ylospain, 
#kustakin jaosta tarkistetaan, onko se hyodyllinen.
#Huom. oikea lapsi loydetaan hakemalla vasemman alipuun loppupiste
#ja siirtymalla yksi eteenpain: i:n oikea lapsi onn endpoint(tree,i+1)+1
#huom pyoristysvirhe vertailussa ??????????
#Kutsuu: poistamon, endpoint, sort.
#
ep<-0.0000001  #pyoristysvirheen huomioon ottaminen
#
nodlkm<-length(ssr)
#
# Muodostetaan vektori, jossa kunkin solmun korkeus
kork<-matrix(0,nodlkm,1)
kork[1]<-1    #juuren korkeus on 1
i<-1
while (i<=nodlkm){
  if (right[i]>0){   #jos ei olla lehdessa
    kork[left[i]]<-kork[i]+1       #vasen lapsi yhta korkeammalla
    kork[right[i]]<-kork[i]+1      #oikea lapsi yhta korkeammalla
  }
  i<-i+1
}
# Muodostetaan vektori, jossa kustakin solmusta alkavan puun log-uskottavuus
#ssralip<-matrix(0,nodlkm,1)
#i<-nodlkm
#while (i>=1){
#  if (right[i]==0)    #jos ollaan lehdessa
#     ssralip[i]<-ssr[i] #alipuun ssr=solmun itsensa ssr
#  else ssralip[i]<-ssralip[left[i]]+ssralip[right[i]]
#            #muuten alipuun ssr=vasemman alipuun ssr + oik alip ssr
#  i<-i-1
#}
# Kaydaan puu lapi taso kerrallaan (miksi?) ja merkitaan muistiin poistettavat
#poistot<-matrix(0,nodlkm,1) #viite<-0   #poistojen maara
ssralip<-matrix(0,nodlkm,1)
tasolkm<-max(kork)     #tasojen lkm
k<-tasolkm             #aloitetaan korkeimmalta tasolta
while (k>0){
  i<-1
  while (i<=nodlkm){
    if (kork[i]==k){
      if (right[i]==0)    #jos ollaan lehdessa
         ssralip[i]<-ssr[i] #alipuun ssr=solmun itsensa ssr
      else{
          ssralip[i]<-ssralip[left[i]]+ssralip[right[i]]
            #muuten alipuun ssr=vasemman alipuun ssr + oik alip ssr
          if (ssralip[i]<=ssr[i]+ep){  
            #jos alipuu ei paranna, se poistetaan, huom pyoristysvirhe
             left[i]<-0 
             right[i]<-0
             #viite<-viite+1  #poistot[viite]<-i
        }
      }
    }
    i<-i+1
  }
  k<-k-1
}
#if (viite>=1){    #jos poistoja on tehty
#  poistot<-poistot[1:viite]
#}
#else{
#  poistot<-NULL
#}
return(list(right=right,left=left))
#return(list(dels=poistot,delnum=viite))
}












intpcf<-function(pcf)
{
value<-pcf$value
down<-pcf$down
high<-pcf$high
support<-pcf$support
N<-pcf$N

d<-length(N)  #dim(down)[2]
step<-stepcalc(support,N)
recnum<-length(value)
int<-0
rr<-1
while (rr<=recnum){
     recint<-matrix(0,2*d,1)
     for (dd in 1:d){
          recint[2*dd-1]<-down[rr,dd]
          recint[2*dd]<-high[rr,dd]
     }
     volint<-massone(recint)
     vol<-prod(step)*volint
     int<-int+value[rr]*vol
     rr<-rr+1
}

return(int)
}




leaflocs<-function(left,right){
#Finds location of leafs in binary tree, living in vector

itemnum<-length(left)
leafloc<-matrix(0,itemnum,1)
pino<-matrix(0,itemnum,1)
pino[1]<-1     #pino[1]<-root
pinin<-1
leafnum<-0
while (pinin>0){
    cur<-pino[pinin]      #take from stack
    pinin<-pinin-1
    if (left[cur]==0){    #if we are in leaf
       leafnum<-leafnum+1
       leafloc[leafnum]<-cur
    }
    else{
       while (left[cur]>0){  #go to leaf and put right nodes to stack
           pinin<-pinin+1
           pino[pinin]<-right[cur]
           cur<-left[cur]
       }
       #now we know we are in leaf
       leafnum<-leafnum+1
       leafloc[leafnum]<-cur  
    }
}
leafloc<-leafloc[1:leafnum]
return(list(leafloc=leafloc,leafnum=leafnum))
} 




leafsum<-function(info,root,left,right){
#Calculates the sum of info over leaves of the subtree
#
#info is itemnum-vector
#root is in 1:itemnum
#left, right are itemnum-vectors, links to childs, 0 if leaf
#
itemnum<-length(info)
sum<-0
pino<-matrix(0,itemnum,1)
pino[1]<-root
pinin<-1
while (pinin>0){
  cur<-pino[pinin]
  pinin<-pinin-1
  if (left[cur]==0){  #if we are in leaf, sum
    sum<-sum+info[cur]
  }
  else{  
    while (left[cur]>0){  #put right on the stack and go to left 
       pinin<-pinin+1
       pino[pinin]<-right[cur]
       cur<-left[cur]
    }
    sum<-sum+info[cur]   #now we know we are in leaf      
  }
}
return(sum)
} 
lefrig2par<-function(et)
{
#from left,right representation to parent representation

d<-length(et$N)
len<-length(et$mean)

parent<-matrix(NA,len,1)
level<-matrix(0,len,1)
volume<-matrix(0,len,1)
center<-matrix(0,d,len)

parent[1]<-0              #root has no children
level[1]<-1
volume[1]<-1

pino<-matrix(0,len,1)
pinoind<-1
pino[1]<-1
while (pinoind>0){
   cur<-pino[pinoind]
   pinoind<-pinoind-1
    
   if (et$left[cur]!=0){
      parent[et$left[cur]]<-cur
      parent[et$right[cur]]<-cur
      level[et$left[cur]]<-level[cur]+1
      level[et$right[cur]]<-level[cur]+1
      volume[et$left[cur]]<-volume[cur]/2
      volume[et$right[cur]]<-volume[cur]/2
      center[,et$left[cur]]<-rep(1,d)          #left is furhtest from origo
      center[,et$right[cur]]<-rep(0,d)
   }
   
   while (et$left[cur]>0){
       #
       # laita oikea pinoon, 
       #
       oikea<-et$right[cur]
       pinoind<-pinoind+1
       pino[pinoind]<-oikea
       #
       # go to left
       #
       cur<-et$left[cur]
       #
       if (et$left[cur]!=0){
          parent[et$left[cur]]<-cur
          parent[et$right[cur]]<-cur
          level[et$left[cur]]<-level[cur]+1
          level[et$right[cur]]<-level[cur]+1
          volume[et$left[cur]]<-volume[cur]/2
          volume[et$right[cur]]<-volume[cur]/2
          center[,et$left[cur]]<-rep(1,d)       #left is furthest from origo
          center[,et$right[cur]]<-rep(0,d)
       }
   }
}

parent2<-matrix(0,len,1)
level2<-matrix(0,len,1)
volume2<-matrix(0,len,1)
center2<-matrix(0,d,len)
codeb<-matrix(0,d,len)
laskuri<-0
for (i in 1:len){
  if (!is.na(parent[i])){
     laskuri<-laskuri+1
     parent2[laskuri]<-parent[i]
     level2[laskuri]<-level[i]
     volume2[laskuri]<-volume[i]
     center2[,laskuri]<-center[,i]
     codeb[laskuri]<-i
  }
}
parent2<-parent2[1:laskuri]
level2<-level2[1:laskuri]
volume2<-volume2[1:laskuri]
center2<-center[,1:laskuri]
codeb<-codeb[1:laskuri] 

i<-2
while (i<=laskuri){
   cod<-parent2[i]
   j<-1
   while ((j<=laskuri) && (cod!=codeb[j])){
      j<-j+1
   }
   parent2[i]<-j
   i<-i+1
}

return(list(parent=parent2,level=level2,center=center2,volume=volume2))
}









levefore<-function(dendat,B,leaf,minlkm=5,seed=1,lambda=0.01,thres=0.5,
sample="bagg",prune="off",
splitscan=0,seedf=1,
scatter=0,
src="c",method="loglik")
{
#B number of bootstrap samples
#leaf number of leaves in the trees to be grown

n<-dim(dendat)[1]
d<-dim(dendat)[2]
suppo<-supp(dendat)+scatter*rep(c(-1,1),d)
step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(suppo[2*i]-suppo[2*i-1])/(n+1)
}

if (sample=="bagg"){
  dendatB<-bootbagg(dendat,seed,scatter) 
}
else{  # scheme=="baggworpl"
  dendatB<-bootworpl(dendat,seed,scatter)
}

tr<-densplit(dendatB,minlkm,suppo,leaf=0,
             method=method,splitscan=splitscan,seedf=seedf)
treeseq<-prunelev(tr,lambda=lambda,n=n)
approleaf<-roundlnum(treeseq$leafs,leaf)
tr<-picktreelev(treeseq,approleaf)

bi<-2
while (bi<=B){
   
   if (sample=="bagg"){
      dendatB<-bootbagg(dendat,seed+bi-1)  
   }
   else{
      dendatB<-bootworpl(dendat,seed+bi-1)
   }

   trnew<-densplit(dendatB,minlkm,suppo,leaf=0,
                   method=method,splitscan=splitscan,seedf=seedf)
   treeseq<-prunelev(trnew,lambda=lambda,n=n)
   approleaf<-roundlnum(treeseq$leafs,leaf)
   trnew<-picktreelev(treeseq,approleaf)

   tr<-treeadd(tr,trnew,bi-1)
   bi<-bi+1
   
}

tr<-c(tr,list(suppo=suppo,step=step))

for (i in 1:length(tr$mean)){
    if (tr$mean[i]>=thres) tr$mean[i]<-1
    else tr$mean[i]<-0
}
#tr$mean<-round(tr$mean)

return(tr)
}







levsplitR<-function(x,minlkm,suppo,lambda=0.1,blokki=50)
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
     maxnoden<-maxnoden+blokki
  }
  
  curparent<-pinoparen[pinin]
  currec<-pinorecs[pinin,]    
  curpoint<-pinopoint[pinin,]
  curbeg<-curpoint[1]
  curend<-curpoint[2]
  pinin<-pinin-1
  
  if (curparent>0) right[curparent]<-curin
  val[curin]<-NA                #aluksi ei puolitettu missaan kohdassa
  vec[curin]<-NA                #aluksi ei puolitettu mitaan muuttujaa
  nelem[curin]<-count(curbeg,curend)                   #havaintojen lkm
  volume[curin]<-massone(currec)
  mean[curin]<-1
  ssr[curin]<-exma(volume[curin],nelem[curin],n,lambda)   #log likeli
  #  jatketaan vas. alipuuhun
  while ((nelem[curin]>minlkm) && (volume[curin]>=minvolume)){
  
    # lisaa varmempi testi ks densplitF ??????
 
    #  koska solmu jaettava, tehdaan jako
    
    jako<-findsplitlev(x,currec,curbeg,curend,obspoint,suppo,n,lambda)  
    
    left[curin]<-curin+1
    val[curin]<-jako$val
    vec[curin]<-jako$vec
    
    rightrec<-jako$rightrec
    leftrec<-jako$leftrec
    leftbeg<-jako$leftbeg
    leftend<-jako$leftend
    rightbeg<-jako$rightbeg
    rightend<-jako$rightend
    obspoint<-jako$obspoint
    #lrindi<-jako$lorr

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
     maxnoden<-maxnoden+blokki
    }
    currec<-leftrec
    curbeg<-leftbeg
    curend<-leftend
    val[curin]<-NA
    vec[curin]<-NA
    nelem[curin]<-count(curbeg,curend)
    volume[curin]<-massone(currec)
    mean[curin]<-1    #lrindi
    ssr[curin]<-exma(volume[curin],nelem[curin],n,lambda)
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
puu<-list(val=val,vec=vec,mean=mean,nelem=nelem,ssr=ssr,volume=volume,
left=left,right=right)
return(puu)
}







lokeroi<-function(A,b){
#
#a is a vector 0<=A_1<...<A_m
#b is number >0
#
#we want to find index i in 1,...,m so that A_{i} < b <= A_i  
#huom jos i=m, niin b>A_{m}
#
m<-length(A)
res<-m
i<-m-1
while (i>=1){
  if (b<=A[i]) res<-i
  i<-i-1
}
#
return(res)
}










lowupp<-function(tr)
{
   nodenum<-length(tr$left)
   d<-length(tr$N)
   low<-matrix(0,nodenum,d)
   upp<-matrix(0,nodenum,d)

   upp[1,]<-tr$N  

   pino<-matrix(0,nodenum,1)
   pinoin<-1
   pino[pinoin]<-1
   while (pinoin>0){      # go through the nodes of tr2
       node<-pino[pinoin]
       pinoin<-pinoin-1

       while (tr$left[node]>0){   
       
          direk<-tr$direc[node]
          split<-tr$split[node]

          leftnode<-tr$left[node]
          rightnode<-tr$right[node]
  
          low[leftnode,]<-low[node,]
          upp[leftnode,]<-upp[node,]
          upp[leftnode,direk]<-tr$split[node]

          low[rightnode,]<-low[node,]
          low[rightnode,direk]<-tr$split[node]
          upp[rightnode,]<-upp[node,]

          # put right node to the stack (if exists)

          pinoin<-pinoin+1
          pino[pinoin]<-rightnode

          # go to left

          node<-leftnode
       }
   }

return(list(low=low,upp=upp))
}

lstseq.bagg<-function(dendat,B,
lstree=NULL,level=NULL,
maxleaf=NULL,leafseq=NULL,
minobs=NULL,seed=1,sample="bagg",prune="off",
splitscan=0,seedf=1,scatter=0,src="c",method="loglik")
{
n<-dim(dendat)[1]
if (is.null(minobs)) minobs<-ceiling(sqrt(n)/2)

if (!is.null(maxleaf)){  
  leafseq<-seq(maxleaf,2)
  hnum<-maxleaf-1
}
else{
  hnum<-length(leafseq)
  if ((hnum>1) && (leafseq[1]>leafseq[2])) leafseq<-leafseq[seq(hnum,1)]
}

for (i in 1:hnum){   
      leaf<-leafseq[i]
      pcf<-eval.bagg(dendat,B,leaf,
           minobs=minobs,seed=seed,
           sample=sample,prune=prune,
           splitscan=splitscan,seedf=seedf,
           scatter=scatter,src=src,method=method)


      if (!is.null(lstree)) lf<-leafsfirst(pcf)
      if (!is.null(level)){ 
           lev<-level*max(pcf$value)  
           refe<-locofmax(pcf)
           st<-leafsfirst(pcf,lev=lev,refe=refe)
      }

      if (i==1){
           if (hnum==1){ 
               pcfseq<-pcf
               if (!is.null(lstree)) lstseq<-lf
               if (!is.null(level)) stseq<-st
           }
           else{ 
               pcfseq<-list(pcf)
               if (!is.null(lstree)) lstseq<-list(lf)
               if (!is.null(level)) stseq<-list(st)
          }
      }
      else{ 
          pcfseq<-c(pcfseq,list(pcf))
          if (!is.null(lstree)) lstseq<-c(lstseq,list(lf))
          if (!is.null(level)) stseq<-c(stseq,list(st))
      }

}

if (is.null(lstree))  lstseq<-NULL
if (is.null(level)) stseq<-NULL
return(list(lstseq=lstseq,pcfseq=pcfseq,stseq=stseq,
hseq=leafseq,type="bagghisto"))
}

lstseq.cart<-function(treeseq,maxleaf=NULL,lstree=NULL,level=NULL,indvec=NULL)
{

if ((is.null(indvec)) && (is.null(maxleaf))){ 
    leaf<-treeseq$leafs
    alfa<-treeseq$alfa
    alkm<-length(leaf)
}
else if (!is.null(maxleaf)){
    aplnum<-roundlnum(treeseq$leafs,maxleaf)
    indeksi<-detsi(treeseq$leafs,aplnum)
    leaf<-treeseq$leafs[indeksi:length(treeseq$leafs)]
    alfa<-treeseq$alfa[indeksi:length(treeseq$leafs)]
    alkm<-length(leaf)
}
else if (!is.null(indvec)){ 
  leaf<-treeseq$leafs[indvec]
  alfa<-treeseq$alfa[indvec]
  alkm<-length(indvec)
}

tuloleaf<-matrix(0,alkm,1)
tuloalfa<-matrix(0,alkm,1)
laskuri<-1

for (inds in alkm:1){  # start with the oversmoothed estimate 
     leafnum<-leaf[inds]
   
     pcf<-eval.pick(treeseq,leafnum)
     #pv<-partition(pcf,suppo)
     #lst<-profgene(pv$values,pv$recs,frekv=F,cvol=T,ccen=T,cfre=F)
     #lst<-proftree(pcf)
     if (!is.null(lstree)) lst<-leafsfirst(pcf)
     if (!is.null(level)){ 
           lev<-level*max(pcf$value)  
           refe<-locofmax(pcf)
           st<-leafsfirst(pcf,lev=lev,refe=refe)
     }

     if (inds==alkm){
        if (alkm==1){
               pcfseq<-pcf
               if (!is.null(lstree)) lstseq<-lst
               if (!is.null(level)) stseq<-st
       }
        else{
               pcfseq<-list(pcf)
               if (!is.null(lstree)) lstseq<-list(lst)
               if (!is.null(level)) stseq<-list(st)  
        }
     }
     else{
          pcfseq<-c(pcfseq,list(pcf))
          if (!is.null(lstree)) lstseq<-c(lstseq,list(lst))
          if (!is.null(level)) stseq<-c(stseq,list(st))
     }
     
    tuloleaf[laskuri]<-leaf[inds]
    tuloalfa[laskuri]<-alfa[inds]
    laskuri<-laskuri+1
}

if (is.null(lstree))  lstseq<-NULL
if (is.null(level)) stseq<-NULL
return(list(lstseq=lstseq,pcfseq=pcfseq,stseq=stseq,hseq=tuloalfa,
leaf=tuloleaf,type="carthisto"))
}

lstseq.greedy<-function(dendat,maxleaf,lstree=NULL,level=NULL)
#treeseq,lsets=FALSE,invalue=FALSE,indvec=NULL,
{
hseq<-seq(maxleaf,1)
hnum<-length(hseq)

for (i in 1:hnum){

     leaf<-hseq[i]
     pcf<-eval.greedy(dendat,leaf=leaf)

     if (!is.null(lstree)) lf<-leafsfirst(pcf)
     if (!is.null(level)){ 
           lev<-level*max(pcf$value)  
           refe<-locofmax(pcf)
           st<-leafsfirst(pcf,lev=lev,refe=refe)
     }
     if (i==1){
           if (hnum==1){ 
               pcfseq<-pcf
               if (!is.null(lstree)) lstseq<-lf
               if (!is.null(level)) stseq<-st
           }
           else{
               pcfseq<-list(pcf)
               if (!is.null(lstree)) lstseq<-list(lf)
               if (!is.null(level)) stseq<-list(st)
           }
      }
      else{
          pcfseq<-c(pcfseq,list(pcf))
          if (!is.null(lstree)) lstseq<-c(lstseq,list(lf))
          if (!is.null(level)) stseq<-c(stseq,list(st))
      }
}

if (is.null(lstree))  lstseq<-NULL
if (is.null(level)) stseq<-NULL
return(list(lstseq=lstseq,pcfseq=pcfseq,stseq=stseq,hseq=hseq,type="greedy"))
}




luo<-function(tree)
{
S<-tree$S
R<-tree$ssr                #R(i) on noden i ssr eli -log likeli

alknodlkm<-length(tree$left)
#leafloc<-findleafs(tree$left,tree$right)

N<-matrix(0,alknodlkm,1)   #number of leaves in the tree whose root is i
p<-matrix(0,alknodlkm,1)   #parent
g<-matrix(0,alknodlkm,1)   #(R(i)-S(i))/(N(i)-1), R(i) on noden i ssr
G<-matrix(0,alknodlkm,1)   #min{g(t),G(l(t)),G(r(t))}

t<-alknodlkm

while (t>=1){
  if (tree$left[t]==0){  #l(t)=0 eli ollaan lehdessa 
     N[t]<-1
     G[t]<-NA                 #\infty
     g[t]<-NA
  }
  else{     #if (!is.na(leafloc[t])){
     p[tree$left[t]]<-t       #p[t+1]<-t
     p[tree$right[t]]<-t      #p[endpoint(tree,t)+1]<-t
     N[t]<-N[tree$left[t]]+N[tree$right[t]]
     g[t]<-(R[t]-S[t])/(N[t]-1)
     G[t]<-omamindelt(g[t],omamindelt(G[tree$left[t]],G[tree$right[t]]))
  }
  t<-t-1
}

return(p=p,G=G,g=g,N=N)
}



makebina<-function(et){

#source("~/delt/R/makebina.R")
#list(frame)
#frame:
#var  splitting variable  or "<leaf>"
#n    number of cases
#dev  deviance of the node
#yval fitted value
#split/splits two-column matrix of the label for the left and right splits
#splits.cutleft  splits.cutright

len<-length(et$split)

var<-matrix("",len,1)
split<-matrix("",len,2)
n<-matrix(0,len,1)
dev<-matrix(0,len,1)
yval<-matrix(0,len,1)

leve<-matrix(0,len,1)    #level for each node
runnum<-matrix(0,len,1)  #ordinal number in the (imaginary) full bin tree 
runredu<-matrix(0,len,1)

leve[1]<-1
runnum[1]<-1

pino<-matrix(0,len,1)
pinoind<-1
pino[1]<-1
laskuri<-0
while (pinoind>0){
   cur<-pino[pinoind]
   pinoind<-pinoind-1
      
   laskuri<-laskuri+1
   n[laskuri]<-et$nelem[cur]                #+1
   dev[laskuri]<-et$ssr[cur]                #abs(et$ssr)
   yval[laskuri]<-et$mean[cur]
   if (et$left[cur]==0){ 
          var[laskuri]<-"<leaf>"
   } 
   else{ 
          var[laskuri]<-paste("x",as.character(et$direc[cur]))
   }
   #var[i]<-as.character(et$direc[i])
   if (et$left[cur]!=0){        #(!is.na(et$split[cur])){
     #split[laskuri,1]<-as.character(et$split[cur])
     #split[laskuri,2]<-as.character(et$split[cur])
     split[laskuri,1]<-format(et$split[cur],digits=2,nsmall=1)
     split[laskuri,2]<-format(et$split[cur],digits=2,nsmall=1)
   }
   runredu[laskuri]<-runnum[cur]
   
   while (et$left[cur]>0){
       #
       # laita oikea pinoon, laske "leve" ja "runnum"
       #
       oikea<-et$right[cur]
       #
       pinoind<-pinoind+1
       pino[pinoind]<-oikea
       #
       leve[oikea]<-leve[cur]+1
       runnum[oikea]<-childcode(leve[cur],runnum[cur])$right
       #
       # count "leve" ja "runnum", go to et$left, update
       #
       vasen<-et$left[cur]
       #
       leve[vasen]<-leve[cur]+1
       runnum[vasen]<-childcode(leve[cur],runnum[cur])$left
       #
       cur<-vasen
       #
       laskuri<-laskuri+1
       n[laskuri]<-et$nelem[cur]                #+1
       dev[laskuri]<-et$ssr[cur]                #abs(et$ssr)
       yval[laskuri]<-et$mean[cur]
       if (et$left[cur]==0){ 
           var[laskuri]<-"<leaf>" 
       } 
       else{ 
           var[laskuri]<-paste("x",as.character(et$direc[cur])) 
       }
       #var[i]<-as.character(et$direc[i])
       if (et$left[cur]!=0){    #if not child
         #split[laskuri,1]<-as.character(et$split[cur])
         #split[laskuri,2]<-as.character(et$split[cur])
         split[laskuri,1]<-format(et$split[cur],digits=2,nsmall=1)
         split[laskuri,2]<-format(et$split[cur],digits=2,nsmall=1)
       }
       runredu[laskuri]<-runnum[cur]
       #
   }
}
var<-var[1:laskuri]
split<-split[1:laskuri,]
n<-n[1:laskuri]
dev<-dev[1:laskuri]
yval<-yval[1:laskuri]
#
runredu<-runredu[1:laskuri]
#
#
frame<-data.frame(var=var,n=n,dev=dev,yval=yval) 
row.names(frame)<-runredu
frame$splits<-array(split,c(laskuri,2),
                    list(character(0), c("cutleft", "cutright")))
#
where<-matrix(1,len,1)
terms<-""
call<-""
y<-FALSE
weigths<-FALSE
#
bintree<-list(frame=frame,where=where,terms=terms,call=call,y=y,weigths=weigths)
attr(bintree,"class")<-"tree" 
return(bintree)
}





makeorder<-function(obspoint,x,coordi){
#Orders obspoint according to coordi:th coordinate of x
#
#obspoint is lkm-vector: pointers to rows of x
#x is n*d-matrix
#coordi is in 1:d
#
lkm<-length(obspoint)
#
redu<-matrix(0,lkm,1)
for (i in 1:lkm){
  obsind<-obspoint[i]
  redu[i]<-x[obsind,coordi]
}
ordobs<-matrix(0,lkm,1)
for (i in 1:lkm){
   pienin<-omaind(redu)
   ordobs[i]<-obspoint[pienin]
   redu[pienin]<-NA  #NA on plus aareton
}
return(ordobs)
}

makeparent<-function(left,right)
{
parent<-matrix(0,length(left),1)

pino<-matrix(0,length(left),1)
pinin<-1
pino[1]<-1

while (pinin>0){

    node<-pino[pinin]
    pinin<-pinin-1

    if (left[node]>0){
       parent[left[node]]<-node
       parent[right[node]]<-node

       pinin<-pinin+1
       pino[pinin]<-right[node]
    }

    while (left[node]>0){
       
        node<-left[node]

        if (left[node]>0){
           parent[left[node]]<-node
           parent[right[node]]<-node

           pinin<-pinin+1
           pino[pinin]<-right[node]
        }
    }

}

return(t(parent))
}
massone.delt<-function(rec){
#Calculates the mass of a rectangle.
#
#rec is (2*d)-vector, represents rectangle in d-space
#Returns a real number >0.
#
d<-length(rec)/2
vol<-1
for (j in 1:d){
  vol<-vol*(rec[2*j]-rec[2*j-1])
}
return(vol) 
}

massone<-function(rec){
#Calculates the mass of a rectangle.
#
#rec is (2*d)-vector, represents rectangle in d-space
#Returns a real number >0.
#
d<-length(rec)/2
vol<-1
for (j in 1:d){
  vol<-vol*(rec[2*j]-rec[2*j-1])
}
return(vol) 
}

meanstd<-function(apu){
#
wv<-dim(apu)[1]
tuldim<-dim(apu)[2]
cv<-matrix(0,tuldim,1)
cvstd<-matrix(0,tuldim,1)
#
for (lo in 1:tuldim){
   for (ek in 1:wv){
      cv[lo]<-cv[lo]+apu[ek,lo]
   }
   cv[lo]<-cv[lo]/wv
   for (ek in 1:wv){
      cvstd[lo]<-cvstd[lo]+(apu[ek,lo]-cv[lo])^2
   }
   cvstd[lo]<-sqrt(cvstd[lo]/wv)
   #
   #cv[lo]<-mean(apu[,lo])              
   #cvstd[lo]<-sd(apu[,lo])
}
#
return(list(cv=cv,cvstd=cvstd))
}
myosplitpena<-function(dendat,leafs,mcdendat,mix,suppo,step,minobs=0,
                       method="projec",bound=0)
{
#suppo<-supp(dendat,blown=TRUE)

n<-length(dendat[,1])  #havaintojen lkm
d<-length(dendat[1,])  #muuttujien lkm

suppvol<-massone(suppo)
mcn<-length(mcdendat[,1])

maxnoden<-2*leafs
val<-matrix(0,maxnoden,1) 
vec<-matrix(0,maxnoden,1)
mean<-matrix(0,maxnoden,1)
loglik<-matrix(0,maxnoden,1)
nelem<-matrix(0,maxnoden,1)
volume<-matrix(0,maxnoden,1)
left<-matrix(0,maxnoden,1)
right<-matrix(0,maxnoden,1)
low<-matrix(0,maxnoden,d)
upp<-matrix(0,maxnoden,d)

recs<-matrix(0,maxnoden,2*d)
begs<-matrix(0,maxnoden,1)
ends<-matrix(0,maxnoden,1)

mcbegs<-matrix(0,maxnoden,1)
mcends<-matrix(0,maxnoden,1)

# for the C-function 

xdendat<-matrix(0,d*n+1,1)
obsoso<-matrix(0,n+1,1)
mcxdendat<-matrix(0,d*mcn+1,1)
mcobsoso<-matrix(0,mcn+1,1)
insuppo<-matrix(0,2*d+1,1)
insuppo[2:(2*d+1)]<-suppo
instep<-matrix(0,d+1,1)
instep[2:(d+1)]<-step
ingrec<-matrix(0,2*d+1,1)
if (method=="loglik") inmethod<-1 else inmethod<-2

# info for root node

val[1]<-0
vec[1]<-0
curvol<-massone(suppo)
mean[1]<-denmean(curvol,n,n)
loglik[1]<-denssr(curvol,n,n,method)
nelem[1]<-n
volume[1]<-curvol
left[1]<-0
right[1]<-0
for (k in 1:d){
  low[1,k]<-0             #suppo[2*k-1]
  upp[1,k]<-n             #suppo[2*k]
  recs[1,2*k-1]<-0        #suppo[2*k-1]
  recs[1,2*k]<-n          #suppo[2*k]
}

begs[1]<-1
ends[1]<-n
mcbegs[1]<-1
mcends[1]<-mcn

# initialize

obspoint<-seq(1:n)
mcobspoint<-seq(1:mcn)

curleafnum<-1
curnodenum<-1

while (curleafnum<leafs){

  locs<-leaflocs(left,right)
  curleafnum<-locs$leafnum

  incre<-matrix(0,curleafnum,1)  
     #for each leaf find the increase in loglik  
     #we choose to make split which increases most the total loglik
  valpool<-matrix(0,curleafnum,1)
  vecpool<-matrix(0,curleafnum,1)
  leftrecpool<-matrix(0,curleafnum,2*d)
  rightrecpool<-matrix(0,curleafnum,2*d)

  leftbegpool<-matrix(0,curleafnum,1)
  leftendpool<-matrix(0,curleafnum,1)
  rightbegpool<-matrix(0,curleafnum,1)
  rightendpool<-matrix(0,curleafnum,1)
  obspointpool<-matrix(0,curleafnum,n)

  mcleftbegpool<-matrix(0,curleafnum,1)
  mcleftendpool<-matrix(0,curleafnum,1)
  mcrightbegpool<-matrix(0,curleafnum,1)
  mcrightendpool<-matrix(0,curleafnum,1)
  mcobspointpool<-matrix(0,curleafnum,mcn)

  failnum<-0  #count the number of nodes where we are not able to make split

  for (i in 1:curleafnum){

     loca<-locs$leafloc[i]
     curloglik<-loglik[loca]
     currec<-recs[loca,]

     curbeg<-begs[loca]
     curend<-ends[loca]
     mccurbeg<-mcbegs[loca]
     mccurend<-mcends[loca]

     {
     if ((curbeg==0) || (curend==0)) maara<-0
     else maara<-count(curbeg,curend)
     }
     {
     if ((mccurbeg==0) || (mccurend==0)) mcmaara<-0
     else mcmaara<-count(mccurbeg,mccurend)
     }

     smallest<-1
     for (j in 1:d){
         valli<-currec[2*j]-currec[2*j-1]
         if (valli>1) smallest<-0
     } 

     #volli<-massone(currec)*prod(step)
     #suhde<-volli/suppvol
  
if ((maara<=minobs) || (smallest==1) || (mcmaara<=minobs)){  # || (suhde<bound)){
        incre[i]<-NA
        failnum<-failnum+1
}

else{

     for (j in curbeg:curend){
           obspointer=obspoint[j]
           obsoso[1+j-curbeg+1]=obspointer
           jin=j-curbeg+1
           for (l in 1:d){
               xdendat[1+(jin-1)*d+l]<-dendat[obspointer,l]
           }
     }
     for (j in mccurbeg:mccurend){
           obspointer=mcobspoint[j]
           mcobsoso[1+j-mccurbeg+1]=obspointer
           jin=j-mccurbeg+1
           for (l in 1:d){
               mcxdendat[1+(jin-1)*d+l]<-mcdendat[obspointer,l]
           }
     }

     ingrec[2:(2*d+1)]<-currec
 
jako<-.C("findsplitpenaC",as.double(xdendat),
                          as.integer(obsoso),
                          as.integer(maara),
                          as.integer(curbeg),
                          as.integer(curend),
                          as.integer(n),
                          # mcdata
                          as.double(mcxdendat),
                          as.integer(mcobsoso),
                          as.integer(mcmaara),
                          as.integer(mccurbeg),
                          as.integer(mccurend),
                          as.integer(mcn),
                          # general
                          as.double(insuppo),
                          as.double(instep),
                          as.integer(ingrec),
                          as.integer(d),
                          as.integer(inmethod),
                          as.double(bound),
                          as.double(mix),
                          # output
                          val = integer(1),
                          vec = integer(1),
                          #
                          leftbeg = integer(1),
                          leftend = integer(1),
                          rightbeg = integer(1),
                          rightend = integer(1),
                          obspoint = integer(maara+1),
                          # same for mc-sample 
                          mcleftbeg = integer(1),
                          mcleftend = integer(1),
                          mcrightbeg = integer(1),
                          mcrightend = integer(1),
                          mcobspoint = integer(mcmaara+1))

     vecu<-jako$vec
     valu<-jako$val   #suppo[2*vecu-1]+jako$val*step[vecu] 
     leftrec<-currec
     leftrec[2*vecu]<-valu
     rightrec<-currec
     rightrec[2*vecu-1]<-valu   

     leftbeg<-jako$leftbeg
     leftend<-jako$leftend
     rightbeg<-jako$rightbeg
     rightend<-jako$rightend
     for (li in 1:maara){
        obspoint[curbeg+li-1]<-jako$obspoint[li+1]
     }

     # same for mc-sample
     mcleftbeg<-jako$mcleftbeg
     mcleftend<-jako$mcleftend
     mcrightbeg<-jako$mcrightbeg
     mcrightend<-jako$mcrightend
     for (li in 1:mcmaara){
        mcobspoint[mccurbeg+li-1]<-jako$mcobspoint[li+1]
     }

     lvolume<-1
     rvolume<-1
     for (ji in 1:d){
          lvolume<-lvolume*(leftrec[2*ji]-leftrec[2*ji-1])*step[ji]
          rvolume<-rvolume*(rightrec[2*ji]-rightrec[2*ji-1])*step[ji]
     }

     lnelem<-count(leftbeg,leftend)     #leftend-leftbeg+1
     rnelem<-count(rightbeg,rightend)   #rightend-rightbeg+1

     mclnelem<-count(mcleftbeg,mcleftend)     #leftend-leftbeg+1
     mcrnelem<-count(mcrightbeg,mcrightend)   #rightend-rightbeg+1

     newloglik<-denssr(lvolume,lnelem,n,method,mix)+denssr(rvolume,rnelem,n,method,mix)+2*(mix-1)*mix*mclnelem*denmean(lvolume,lnelem,n)/mcn+2*(mix-1)*mix*mcrnelem*denmean(rvolume,rnelem,n)/mcn

     incre[i]<-newloglik-curloglik
  
     if ((lvolume/suppvol < bound) && (lnelem>0)) incre[i]<-NA
     if ((rvolume/suppvol < bound) && (rnelem>0)) incre[i]<-NA
 
     valpool[i]<-valu
     vecpool[i]<-vecu
     leftrecpool[i,]<-leftrec
     rightrecpool[i,]<-rightrec

     leftbegpool[i]<-leftbeg
     leftendpool[i]<-leftend
     rightbegpool[i]<-rightbeg
     rightendpool[i]<-rightend
     obspointpool[i,]<-obspoint

     mcleftbegpool[i]<-mcleftbeg
     mcleftendpool[i]<-mcleftend
     mcrightbegpool[i]<-mcrightbeg
     mcrightendpool[i]<-mcrightend
     mcobspointpool[i,]<-mcobspoint

} #else (we may split because there are observations in the rec)

  }   #for (i in 1:curleafnum){

####################################################

# now "incre" is the optimization vector, we find the
# maximum of this vector: this gives the best split

allfail<-1
for (ii in 1:d){
  if (!is.na(incre[ii]))  allfail<-0
}

sd<-omaind(-incre)   #omaind minimizes, we want to maximize

if ((failnum==curleafnum)){  # || (allfail==1)){  # we have to finish

cl<-curnodenum

return(list(split=val[1:cl],direc=vec[1:cl],mean=mean[1:cl],nelem=nelem[1:cl],
ssr=loglik[1:cl],volume=volume[1:cl],
left=left[1:cl],right=right[1:cl],low=low[1:cl,],upp=upp[1:cl,],
support=suppo,step=step))

}

else{

  # make the split sd

  locloc<-locs$leafloc[sd]
  val[locloc]<-valpool[sd] 
  vec[locloc]<-vecpool[sd]

  # create left child

  leftpoint<-curnodenum+1
  left[locloc]<-leftpoint

  recu<-leftrecpool[sd,]
  volu<-1           #volu<-massone(recu)
  for (ji in 1:d){
     volu<-volu*(recu[2*ji]-recu[2*ji-1])*step[ji]
  }
  nelemu<-count(leftbegpool[sd],leftendpool[sd])  
        #leftendpool[sd]-leftbegpool[sd]+1

  val[leftpoint]<-0
  vec[leftpoint]<-0
  mean[leftpoint]<-denmean(volu,nelemu,n)
  loglik[leftpoint]<-denssr(volu,nelemu,n,method)
  nelem[leftpoint]<-nelemu
  volume[leftpoint]<-volu
  for (k in 1:d){
    low[leftpoint,k]<-recu[2*k-1]
    upp[leftpoint,k]<-recu[2*k]
  }
  upp[leftpoint,vec[locloc]]<-val[locloc]

  recs[leftpoint,]<-recu

  begs[leftpoint]<-leftbegpool[sd]
  ends[leftpoint]<-leftendpool[sd]
  mcbegs[leftpoint]<-mcleftbegpool[sd]
  mcends[leftpoint]<-mcleftendpool[sd]

  # create right child

  rightpoint<-curnodenum+2
  right[locloc]<-rightpoint

  recu<-rightrecpool[sd,] 
  volu<-1
  for (ji in 1:d){
     volu<-volu*(recu[2*ji]-recu[2*ji-1])*step[ji]
  }
  nelemu<-count(rightbegpool[sd],rightendpool[sd]) 
          #rightendpool[sd]-rightbegpool[sd]+1

  val[rightpoint]<-0
  vec[rightpoint]<-0
  mean[rightpoint]<-denmean(volu,nelemu,n)
  loglik[rightpoint]<-denssr(volu,nelemu,n,method)
  nelem[rightpoint]<-nelemu
  volume[rightpoint]<-volu
  for (k in 1:d){
    low[rightpoint,k]<-recu[2*k-1]
    upp[rightpoint,k]<-recu[2*k]
  }
  low[rightpoint,vec[locloc]]<-val[locloc]

  recs[rightpoint,]<-recu
  begs[rightpoint]<-rightbegpool[sd]
  ends[rightpoint]<-rightendpool[sd]
  mcbegs[rightpoint]<-mcrightbegpool[sd]
  mcends[rightpoint]<-mcrightendpool[sd]

  # final updates

  curleafnum<-curleafnum+1
  curnodenum<-curnodenum+2
  obspoint<-obspointpool[sd,]
  mcobspoint<-mcobspointpool[sd,]

}  #end split making


}  #while (curleafnum<leafs)

cl<-curnodenum

return(list(split=val[1:cl],direc=vec[1:cl],mean=mean[1:cl],nelem=nelem[1:cl],
ssr=loglik[1:cl],volume=volume[1:cl],
left=left[1:cl],right=right[1:cl],low=low[1:cl,],upp=upp[1:cl,],
support=suppo,step=step))

}















myosplitpenaR<-function(dendat,leafs,mcdendat,mix,suppo,step,minobs=0,
                        method="projec",bound=0)
{
#suppo<-supp(dendat,blown=TRUE)

n<-length(dendat[,1])  #havaintojen lkm
d<-length(dendat[1,])  #muuttujien lkm

suppvol<-massone(suppo)
mcn<-length(mcdendat[,1])

maxnoden<-2*leafs
val<-matrix(0,maxnoden,1) 
vec<-matrix(0,maxnoden,1)
mean<-matrix(0,maxnoden,1)
loglik<-matrix(0,maxnoden,1)
nelem<-matrix(0,maxnoden,1)
volume<-matrix(0,maxnoden,1)
left<-matrix(0,maxnoden,1)
right<-matrix(0,maxnoden,1)
low<-matrix(0,maxnoden,d)
upp<-matrix(0,maxnoden,d)

recs<-matrix(0,maxnoden,2*d)
begs<-matrix(0,maxnoden,1)
ends<-matrix(0,maxnoden,1)

mcbegs<-matrix(0,maxnoden,1)
mcends<-matrix(0,maxnoden,1)

# for the C-function 

xdendat<-matrix(0,d*n+1,1)
obsoso<-matrix(0,n+1,1)
mcxdendat<-matrix(0,d*mcn+1,1)
mcobsoso<-matrix(0,mcn+1,1)
insuppo<-matrix(0,2*d+1,1)
insuppo[2:(2*d+1)]<-suppo
instep<-matrix(0,d+1,1)
instep[2:(d+1)]<-step
ingrec<-matrix(0,2*d+1,1)
if (method=="loglik") inmethod<-1 else inmethod<-2

# info for root node

val[1]<-0
vec[1]<-0
curvol<-massone(suppo)
mean[1]<-denmean(curvol,n,n)
loglik[1]<-denssr(curvol,n,n,method)
nelem[1]<-n
volume[1]<-curvol
left[1]<-0
right[1]<-0
for (k in 1:d){
  low[1,k]<-0             #suppo[2*k-1]
  upp[1,k]<-n             #suppo[2*k]
  recs[1,2*k-1]<-0        #suppo[2*k-1]
  recs[1,2*k]<-n          #suppo[2*k]
}

begs[1]<-1
ends[1]<-n
mcbegs[1]<-1
mcends[1]<-mcn

# initialize

obspoint<-seq(1:n)
mcobspoint<-seq(1:mcn)

curleafnum<-1
curnodenum<-1

while (curleafnum<leafs){

  locs<-leaflocs(left,right)
  curleafnum<-locs$leafnum

  incre<-matrix(0,curleafnum,1)  
     #for each leaf find the increase in loglik  
     #we choose to make split which increases most the total loglik
  valpool<-matrix(0,curleafnum,1)
  vecpool<-matrix(0,curleafnum,1)
  leftrecpool<-matrix(0,curleafnum,2*d)
  rightrecpool<-matrix(0,curleafnum,2*d)

  leftbegpool<-matrix(0,curleafnum,1)
  leftendpool<-matrix(0,curleafnum,1)
  rightbegpool<-matrix(0,curleafnum,1)
  rightendpool<-matrix(0,curleafnum,1)
  obspointpool<-matrix(0,curleafnum,n)

  mcleftbegpool<-matrix(0,curleafnum,1)
  mcleftendpool<-matrix(0,curleafnum,1)
  mcrightbegpool<-matrix(0,curleafnum,1)
  mcrightendpool<-matrix(0,curleafnum,1)
  mcobspointpool<-matrix(0,curleafnum,mcn)

  failnum<-0  #count the number of nodes where we are not able to make split

  for (i in 1:curleafnum){

     loca<-locs$leafloc[i]
     curloglik<-loglik[loca]
     currec<-recs[loca,]

     curbeg<-begs[loca]
     curend<-ends[loca]
     mccurbeg<-mcbegs[loca]
     mccurend<-mcends[loca]

     {
     if ((curbeg==0) || (curend==0)) maara<-0
     else maara<-count(curbeg,curend)
     }
     {
     if ((mccurbeg==0) || (mccurend==0)) mcmaara<-0
     else mcmaara<-count(mccurbeg,mccurend)
     }

     smallest<-1
     for (j in 1:d){
         valli<-currec[2*j]-currec[2*j-1]
         if (valli>1) smallest<-0
     } 

     #volli<-massone(currec)*prod(step)
     #suhde<-volli/suppvol
  
if ((maara<=minobs) || (smallest==1) || (mcmaara<=minobs)){  # || (suhde<bound)){
        incre[i]<-NA
        failnum<-failnum+1
}

else{

     for (j in curbeg:curend){
           obspointer=obspoint[j]
           obsoso[1+j-curbeg+1]=obspointer
           jin=j-curbeg+1
           for (l in 1:d){
               xdendat[1+(jin-1)*d+l]<-dendat[obspointer,l]
           }
     }
     for (j in mccurbeg:mccurend){
           obspointer=mcobspoint[j]
           mcobsoso[1+j-mccurbeg+1]=obspointer
           jin=j-mccurbeg+1
           for (l in 1:d){
               mcxdendat[1+(jin-1)*d+l]<-mcdendat[obspointer,l]
           }
     }

     ingrec[2:(2*d+1)]<-currec
 
jako<-findsplitpenaR(dendat,currec,curbeg,curend,obspoint,
                     mcdendat,mccurbeg,mccurend,mcobspoint,mix,suppo,step)

     vecu<-jako$vec
     valu<-jako$val   #suppo[2*vecu-1]+jako$val*step[vecu] 
     leftrec<-currec
     leftrec[2*vecu]<-valu
     rightrec<-currec
     rightrec[2*vecu-1]<-valu   

     leftbeg<-jako$leftbeg
     leftend<-jako$leftend
     rightbeg<-jako$rightbeg
     rightend<-jako$rightend
     for (li in 1:maara){
        obspoint[curbeg+li-1]<-jako$obspoint[li]
     }

     # same for mc-sample
     mcleftbeg<-jako$mcleftbeg
     mcleftend<-jako$mcleftend
     mcrightbeg<-jako$mcrightbeg
     mcrightend<-jako$mcrightend
     for (li in 1:mcmaara){
        mcobspoint[mccurbeg+li-1]<-jako$mcobspoint[li]
     }

     lvolume<-1
     rvolume<-1
     for (ji in 1:d){
          lvolume<-lvolume*(leftrec[2*ji]-leftrec[2*ji-1])*step[ji]
          rvolume<-rvolume*(rightrec[2*ji]-rightrec[2*ji-1])*step[ji]
     }

     lnelem<-count(leftbeg,leftend)     #leftend-leftbeg+1
     rnelem<-count(rightbeg,rightend)   #rightend-rightbeg+1

     mclnelem<-count(mcleftbeg,mcleftend)     #leftend-leftbeg+1
     mcrnelem<-count(mcrightbeg,mcrightend)   #rightend-rightbeg+1

     newloglik<-denssr(lvolume,lnelem,n,method,mix)+denssr(rvolume,rnelem,n,method,mix)+2*(mix-1)*mix*mclnelem*denmean(lvolume,lnelem,n)/mcn+2*(mix-1)*mix*mcrnelem*denmean(rvolume,rnelem,n)/mcn

     incre[i]<-newloglik-curloglik
  
     if ((lvolume/suppvol < bound) && (lnelem>0)) incre[i]<-NA
     if ((rvolume/suppvol < bound) && (rnelem>0)) incre[i]<-NA
 
     valpool[i]<-valu
     vecpool[i]<-vecu
     leftrecpool[i,]<-leftrec
     rightrecpool[i,]<-rightrec

     leftbegpool[i]<-leftbeg
     leftendpool[i]<-leftend
     rightbegpool[i]<-rightbeg
     rightendpool[i]<-rightend
     obspointpool[i,]<-obspoint

     mcleftbegpool[i]<-mcleftbeg
     mcleftendpool[i]<-mcleftend
     mcrightbegpool[i]<-mcrightbeg
     mcrightendpool[i]<-mcrightend
     mcobspointpool[i,]<-mcobspoint

} #else (we may split because there are observations in the rec)

  }   #for (i in 1:curleafnum){

####################################################

# now "incre" is the optimization vector, we find the
# maximum of this vector: this gives the best split

allfail<-1
for (ii in 1:d){
  if (!is.na(incre[ii]))  allfail<-0
}

sd<-omaind(-incre)   #omaind minimizes, we want to maximize

if ((failnum==curleafnum)){  # || (allfail==1)){  # we have to finish

cl<-curnodenum

return(list(val=val[1:cl],vec=vec[1:cl],mean=mean[1:cl],nelem=nelem[1:cl],
ssr=loglik[1:cl],volume=volume[1:cl],
left=left[1:cl],right=right[1:cl],low=low[1:cl,],upp=upp[1:cl,],
suppo=suppo,step=step))

}

else{

  # make the split sd

  locloc<-locs$leafloc[sd]
  val[locloc]<-valpool[sd] 
  vec[locloc]<-vecpool[sd]

  # create left child

  leftpoint<-curnodenum+1
  left[locloc]<-leftpoint

  recu<-leftrecpool[sd,]
  volu<-1           #volu<-massone(recu)
  for (ji in 1:d){
     volu<-volu*(recu[2*ji]-recu[2*ji-1])*step[ji]
  }
  nelemu<-count(leftbegpool[sd],leftendpool[sd])  
        #leftendpool[sd]-leftbegpool[sd]+1

  val[leftpoint]<-0
  vec[leftpoint]<-0
  mean[leftpoint]<-denmean(volu,nelemu,n)
  loglik[leftpoint]<-denssr(volu,nelemu,n,method)
  nelem[leftpoint]<-nelemu
  volume[leftpoint]<-volu
  for (k in 1:d){
    low[leftpoint,k]<-recu[2*k-1]
    upp[leftpoint,k]<-recu[2*k]
  }
  upp[leftpoint,vec[locloc]]<-val[locloc]

  recs[leftpoint,]<-recu

  begs[leftpoint]<-leftbegpool[sd]
  ends[leftpoint]<-leftendpool[sd]
  mcbegs[leftpoint]<-mcleftbegpool[sd]
  mcends[leftpoint]<-mcleftendpool[sd]

  # create right child

  rightpoint<-curnodenum+2
  right[locloc]<-rightpoint

  recu<-rightrecpool[sd,] 
  volu<-1
  for (ji in 1:d){
     volu<-volu*(recu[2*ji]-recu[2*ji-1])*step[ji]
  }
  nelemu<-count(rightbegpool[sd],rightendpool[sd]) 
          #rightendpool[sd]-rightbegpool[sd]+1

  val[rightpoint]<-0
  vec[rightpoint]<-0
  mean[rightpoint]<-denmean(volu,nelemu,n)
  loglik[rightpoint]<-denssr(volu,nelemu,n,method)
  nelem[rightpoint]<-nelemu
  volume[rightpoint]<-volu
  for (k in 1:d){
    low[rightpoint,k]<-recu[2*k-1]
    upp[rightpoint,k]<-recu[2*k]
  }
  low[rightpoint,vec[locloc]]<-val[locloc]

  recs[rightpoint,]<-recu
  begs[rightpoint]<-rightbegpool[sd]
  ends[rightpoint]<-rightendpool[sd]
  mcbegs[rightpoint]<-mcrightbegpool[sd]
  mcends[rightpoint]<-mcrightendpool[sd]

  # final updates

  curleafnum<-curleafnum+1
  curnodenum<-curnodenum+2
  obspoint<-obspointpool[sd,]
  mcobspoint<-mcobspointpool[sd,]

}  #end split making


}  #while (curleafnum<leafs)

cl<-curnodenum

return(list(val=val[1:cl],vec=vec[1:cl],mean=mean[1:cl],nelem=nelem[1:cl],
ssr=loglik[1:cl],volume=volume[1:cl],
left=left[1:cl],right=right[1:cl],low=low[1:cl,],upp=upp[1:cl,],
suppo=suppo,step=step))

}















myosplitR<-function(dendat,leafs,method="loglik",minobs=0)
{
# split points: 1,...,n-1
# suppo: [0,n]  =  [min-step/2,max+step/2]

suppo<-supp(dendat,blown=TRUE)

n<-length(dendat[,1])  #havaintojen lkm
d<-length(dendat[1,])  #muuttujien lkm

maxnoden<-2*leafs
val<-matrix(0,maxnoden,1) 
vec<-matrix(0,maxnoden,1)
mean<-matrix(0,maxnoden,1)
loglik<-matrix(0,maxnoden,1)
nelem<-matrix(0,maxnoden,1)
volume<-matrix(0,maxnoden,1)
left<-matrix(0,maxnoden,1)
right<-matrix(0,maxnoden,1)
low<-matrix(0,maxnoden,d)
upp<-matrix(0,maxnoden,d)

recs<-matrix(0,maxnoden,2*d)
begs<-matrix(0,maxnoden,1)
ends<-matrix(0,maxnoden,1)

step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(suppo[2*i]-suppo[2*i-1])/n  #(n+1)
}

# info for root node

val[1]<-0
vec[1]<-0
curvol<-massone(suppo)
mean[1]<-denmean(curvol,n,n)
loglik[1]<-denssr(curvol,n,n,method)
nelem[1]<-n
volume[1]<-curvol
left[1]<-0
right[1]<-0
for (k in 1:d){
  low[1,k]<-0            #suppo[2*k-1]
  upp[1,k]<-n            #suppo[2*k]
  recs[1,2*k-1]<-0       #suppo[2*k-1]
  recs[1,2*k]<-n         #suppo[2*k]
}
begs[1]<-1
ends[1]<-n

# initialize

obspoint<-seq(1:n)
curleafnum<-1
curnodenum<-1

while (curleafnum<leafs){

  locs<-leaflocs(left,right)
  curleafnum<-locs$leafnum

  incre<-matrix(0,curleafnum,1)  
     #for each leaf find the increase in loglik  
     #we choose to make split which increases most the total loglik
  valpool<-matrix(0,curleafnum,1)
  vecpool<-matrix(0,curleafnum,1)
  leftrecpool<-matrix(0,curleafnum,2*d)
  rightrecpool<-matrix(0,curleafnum,2*d)
  leftbegpool<-matrix(0,curleafnum,1)
  leftendpool<-matrix(0,curleafnum,1)
  rightbegpool<-matrix(0,curleafnum,1)
  rightendpool<-matrix(0,curleafnum,1)
  obspointpool<-matrix(0,curleafnum,n)

  failnum<-0  # count the number of nodes where we are not able to make split

  for (i in 1:curleafnum){

     loca<-locs$leafloc[i]
     currec<-recs[loca,]
     curbeg<-begs[loca]
     curend<-ends[loca]
     curloglik<-loglik[loca]
     {
     if ((curbeg==0) || (curend==0)) maara<-0
     else maara<-count(curbeg,curend)
     }

     smallest<-1
     for (j in 1:d){
         valli<-currec[2*j]-currec[2*j-1]
         if (valli>1) smallest<-0
     }   

if ((maara<=minobs) || (smallest==1)){
        incre[i]<-NA
        failnum<-failnum+1
}

else{

     #for (j in curbeg:curend){
     #      obspointer=obspoint[j]
     #      obsoso[j-curbeg+1]=obspointer
     #}

     jako<-findsplitG(dendat,currec,  #inrec
                      curbeg,curend,obspoint,  #obsoso,
                      suppo,n,method)                

     vecu<-jako$vec
     valu<-jako$val   #+currec[2*vecu-1]  
     leftrec<-currec
     leftrec[2*vecu]<-valu
     rightrec<-currec
     rightrec[2*vecu-1]<-valu   

     leftbeg<-jako$leftbeg
     leftend<-jako$leftend
     rightbeg<-jako$rightbeg
     rightend<-jako$rightend
     #for (li in 1:maara){
     #    obspoint[curbeg+li-1]<-jako$obspoint[li]
     #}
     obspoint<-jako$obspoint

     lvolume<-1
     rvolume<-1
     for (ji in 1:d){
          lvolume<-lvolume*(leftrec[2*ji]-leftrec[2*ji-1])*step[ji]
          rvolume<-rvolume*(rightrec[2*ji]-rightrec[2*ji-1])*step[ji]
     }

     lnelem<-count(leftbeg,leftend)     #leftend-leftbeg+1
     rnelem<-count(rightbeg,rightend)   #rightend-rightbeg+1
     newloglik<-denssr(lvolume,lnelem,n,method)+
                denssr(rvolume,rnelem,n,method)
     incre[i]<-newloglik-curloglik
  
     valpool[i]<-valu
     vecpool[i]<-vecu
     leftrecpool[i,]<-leftrec
     rightrecpool[i,]<-rightrec
     leftbegpool[i]<-leftbeg
     leftendpool[i]<-leftend
     rightbegpool[i]<-rightbeg
     rightendpool[i]<-rightend
     obspointpool[i,]<-obspoint

} #else (we may split because there are observations in the rec)

  }   #for (i in 1:curleafnum){

sd<-omaind(-incre)   #omaind minimizes, we want to maximize

if (failnum==curleafnum){  # we have to finish

cl<-curnodenum

return(list(val=val[1:cl],vec=vec[1:cl],mean=mean[1:cl],nelem=nelem[1:cl],
ssr=loglik[1:cl],volume=volume[1:cl],
left=left[1:cl],right=right[1:cl],low=low[1:cl,],upp=upp[1:cl,],
suppo=suppo,step=step))

}

else{

  # make the split sd

  locloc<-locs$leafloc[sd]
  val[locloc]<-valpool[sd] 
  vec[locloc]<-vecpool[sd]

  # create left child

  leftpoint<-curnodenum+1
  left[locloc]<-leftpoint

  recu<-leftrecpool[sd,]
  volu<-1
  for (ji in 1:d){
     volu<-volu*(recu[2*ji]-recu[2*ji-1])*step[ji]
  }
  nelemu<-count(leftbegpool[sd],leftendpool[sd]) #leftendpool[sd]-leftbegpool[sd]+1

  val[leftpoint]<-0
  vec[leftpoint]<-0
  mean[leftpoint]<-denmean(volu,nelemu,n)
  loglik[leftpoint]<-denssr(volu,nelemu,n,method)
  nelem[leftpoint]<-nelemu
  volume[leftpoint]<-volu
  for (k in 1:d){
    low[leftpoint,k]<-recu[2*k-1]
    upp[leftpoint,k]<-recu[2*k]
  }
  upp[leftpoint,vec[locloc]]<-val[locloc]

  recs[leftpoint,]<-recu
  begs[leftpoint]<-leftbegpool[sd]
  ends[leftpoint]<-leftendpool[sd]

  # create right child

  rightpoint<-curnodenum+2
  right[locloc]<-rightpoint

  recu<-rightrecpool[sd,] 
  #volu<-massone(recu)
  volu<-1
  for (ji in 1:d){
     volu<-volu*(recu[2*ji]-recu[2*ji-1])*step[ji]
  }
  nelemu<-count(rightbegpool[sd],rightendpool[sd]) #rightendpool[sd]-rightbegpool[sd]+1

  val[rightpoint]<-0
  vec[rightpoint]<-0
  mean[rightpoint]<-denmean(volu,nelemu,n)
  loglik[rightpoint]<-denssr(volu,nelemu,n,method)
  nelem[rightpoint]<-nelemu
  volume[rightpoint]<-volu
  for (k in 1:d){
    low[rightpoint,k]<-recu[2*k-1]
    upp[rightpoint,k]<-recu[2*k]
  }
  low[rightpoint,vec[locloc]]<-val[locloc]

  recs[rightpoint,]<-recu
  begs[rightpoint]<-rightbegpool[sd]
  ends[rightpoint]<-rightendpool[sd]

  # final updates

  curleafnum<-curleafnum+1
  curnodenum<-curnodenum+2
  obspoint<-obspointpool[sd,]

}  #end split making

}  #while (curleafnum<leafs)

cl<-curnodenum

return(list(val=val[1:cl],vec=vec[1:cl],mean=mean[1:cl],nelem=nelem[1:cl],
ssr=loglik[1:cl],volume=volume[1:cl],
left=left[1:cl],right=right[1:cl],low=low[1:cl,],upp=upp[1:cl,],
suppo=suppo,step=step))

}















omaind.delt<-function(v){
#v on vektori, palautetaan indeksi jossa vektorin pienin arvo 
#
lkm<-length(v)
i<-1
while ((i<lkm) && (is.na(v[i]))) i<-i+1
if ((i==lkm) && (is.na(v[lkm]))) y<-1
 else
 if ((i==lkm) && (!is.na(v[lkm]))) y<-lkm
  else{
  apuu<-i
  valapu<-v[apuu]
  while (i<lkm){
    i<-i+1
    if ((!is.na(v[i])) && (v[i] < valapu)){
      apuu<-i
      valapu<-v[i]
    }
  }
y<-apuu
  }
return(y)
}
omaind<-function(v){
#v on vektori, palautetaan indeksi jossa vektorin pienin arvo 
#
lkm<-length(v)
i<-1
while ((i<lkm) && (is.na(v[i]))) i<-i+1
if ((i==lkm) && (is.na(v[lkm]))) y<-1
 else
 if ((i==lkm) && (!is.na(v[lkm]))) y<-lkm
  else{
  apuu<-i
  valapu<-v[apuu]
  while (i<lkm){
    i<-i+1
    if ((!is.na(v[i])) && (v[i] < valapu)){
      apuu<-i
      valapu<-v[i]
    }
  }
y<-apuu
  }
return(y)
}
omamindelt<-function(a,b){
#min(a,b), NA=infty, 
tulos<-F
if (is.na(a)) tulos<-b 
else{ if (is.na(b)) tulos<-a
      else{ if (a<=b) tulos<-a else tulos<-b}
}
return(tulos)
}
omaord2.delt<-function(a,b){
#Jarjestaa vektorin a vektorin b mukaiseen jarjestykseen
#
#a and b are lnum-vectors
#
lnum<-length(a)  
orda<-a               #tahan oikea jarjestys
ordb<-b
i<-1 
while (i<=lnum){
   pienin<-omaind(b)
   ordb[i]<-b[pienin]
   orda[i]<-a[pienin]
   b[pienin]<-NA         #NA on plus aareton
   i<-i+1
}
return(orda)
}





omaver<-function(a,b){
#(a<b), NA=infty
tulos<-F
if (is.na(a)) tulos<-F 
else{ if (is.na(b)) tulos<-T
      else{ if (a<b) tulos<-T}
}
return(tulos)
}
ositalog<-function(n,osimat,sara){
#
rownum<-dim(osimat)[1]
#rividim<-dim(osimat)[2]
#
#Lasketaan otoskoko
#las<-0
#j<-1
#while (j<=rividim)
#  if (osimat[saradim,j]!=NaN) las<-las+1
#  endif
#  j<-j+1
#endo
#n<-(saradim-1)*rividim+las
#Ryhdytaan paatehtavaan
#
tulos<-matrix(0,n,1)
fal<-0
tru<-1
#
#i=1
#while (i<=n)
#   tulos(i,1)<-false
#endo
#
i<-1
while (i<=rownum){
  if (!is.na(osimat[i,sara])) tulos[osimat[i,sara]]<-tru 
  i<-i+1
}
return(tulos)
}
osita<-function(n,wv,seed){
#
#if (wv>n) error
#
set.seed(seed)
#
koko<-floor(n/wv)
tulos<-matrix(0,koko+1,wv)
arpavec<-matrix(1,n,1)
i<-1
while (i<=n){
  arpavec[i]<-i
  i<-i+1
}
i<-1
while (i<=koko){
j<-1
while (j<=wv){
   uusidim<-n-((i-1)*wv+j)+1  
   arpa<-unidis(uusidim)
   tulos[i,j]<-arpavec[arpa]
   if (arpa==1){ 
       arpavec<-arpavec[2:uusidim]
   }
   else{
      if (arpa==uusidim){ 
        arpavec<-arpavec[1:(uusidim-1)]
      }
      else{ 
         arpavecnew<-matrix(0,uusidim-1,1)
         arpavecnew[1:(arpa-1)]<-arpavec[1:(arpa-1)]
         arpavecnew[arpa:(uusidim-1)]<-arpavec[(arpa+1):uusidim]
         arpavec<-arpavecnew
      }
   }
   j<-j+1
}
i<-i+1
} 
ylipitlkm<-n-wv*koko
j<-1
while (j<=ylipitlkm){
   uusidim<-n-(koko*wv+j)+1  
   arpa<-unidis(uusidim)
   tulos[koko+1,j]<-arpavec[arpa]
   if (arpa==1){ 
      arpavec=arpavec[2:uusidim]
   }
   else{
     if (arpa==uusidim){ 
        arpavec=arpavec[1:(uusidim-1)]
     }
     else{ 
         arpavecnew<-matrix(0,uusidim-1,1)
         arpavecnew[1:(arpa-1)]<-arpavec[1:(arpa-1)]
         arpavecnew[arpa:(uusidim-1)]<-arpavec[(arpa+1):uusidim]
         arpavec<-arpavecnew
     }
   } 
   j<-j+1
}
j<-ylipitlkm+1
while (j<=wv){
   tulos[koko+1,j]<-NA
   j<-j+1
}
#
return(tulos)
}





paf<-function(a,b){
#
rown<-dim(a)[1]
coln<-dim(a)[2]
res<-matrix(0,rown,coln)
#
counter<-0
for (i in 1:rown){
  if (b[i]==1){
     counter<-counter+1 
     res[counter,]<-a[i,]
  }
}
#
res<-res[1:counter,]
return(res)
}
partigen.disc<-function(tr)
{
cl<-length(tr$left)
d<-length(tr$N)

down<-matrix(0,cl,d)
high<-matrix(0,cl,d)

#for (i in 1:cl){
#   for (j in 1:d){
#      down[i,j]<-tr$low[i,j]+1
#      high[i,j]<-tr$upp[i,j]
#    }
#}

ll<-leaflocs(tr$left[1:cl],tr$right[1:cl])
leafloc<-ll$leafloc
leafnum<-ll$leafnum

value<-matrix(0,leafnum,1)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if (tr$mean[node]>0){
     efek<-efek+1

     value[efek]<-tr$mean[node]
 
     for (j in 1:d){
         down[efek,j]<-tr$low[node,j]
         high[efek,j]<-tr$upp[node,j]
     }
   }
   i<-i+1
}
value<-value[1:efek]
down<-down[1:efek,]
high<-high[1:efek,]

return(list(value=value,down=down,high=high))
}
partigenD<-function(tr,grid=TRUE,zerorecs=FALSE)
{
d<-length(tr$N)
step<-stepcalc(tr$support,tr$N)

nodenum<-length(tr$left)
left<-tr$left
right<-tr$right
{
if (grid){
  low<-matrix(0,nodenum,d)
  upp<-matrix(0,nodenum,d)
  for (i in 1:nodenum){
    for (j in 1:d){
      low[i,j]<-tr$low[i,j]*step[j]+tr$support[2*j-1]
      upp[i,j]<-tr$upp[i,j]*step[j]+tr$support[2*j-1]
    }
  }
}
else{
  low<-tr$low
  upp<-tr$upp
}
}
mean<-tr$mean

ll<-leaflocs(left,right)
leafloc<-ll$leafloc
leafnum<-ll$leafnum

values<-matrix(0,leafnum,1)
recs<-matrix(0,leafnum,2*d)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if ((mean[node]>0) || (zerorecs)){

     efek<-efek+1

     values[efek]<-mean[node]
 
     for (k in 1:d){

       recs[efek,2*k-1]<-low[node,k]
       recs[efek,2*k]<-upp[node,k]
     }

   }

   i<-i+1
}
values<-values[1:efek]
recs<-recs[1:efek,]
return(list(values=values,recs=recs))
}













partigen<-function(tr,grid=TRUE,zerorecs=FALSE)
{

d<-2
nodenum<-length(tr$left)
left<-tr$left
right<-tr$right
if (is.null(tr$step)) tr$step<-stepcalc(tr$support,tr$N)

if (grid){
  low<-matrix(0,nodenum,d)
  upp<-matrix(0,nodenum,d)
  for (i in 1:nodenum){
    for (j in 1:d){
      low[i,j]<-tr$low[i,j]*tr$step[j]+tr$support[2*j-1]
      upp[i,j]<-tr$upp[i,j]*tr$step[j]+tr$support[2*j-1]
    }
  }
}
else{
  low<-tr$low
  upp<-tr$upp
}

mean<-tr$mean

ll<-leaflocs(left,right)
leafloc<-ll$leafloc
leafnum<-ll$leafnum

values<-matrix(0,leafnum,1)
recs<-matrix(0,leafnum,2*d)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if ((mean[node]>0) || (zerorecs)){

     efek<-efek+1

     values[efek]<-mean[node]
 
     recs[efek,1]<-low[node,1]
     recs[efek,2]<-upp[node,1]
  
     recs[efek,3]<-low[node,2]
     recs[efek,4]<-upp[node,2]
   }

   i<-i+1
}
values<-values[1:efek]
recs<-recs[1:efek,]
return(list(values=values,recs=recs))
}













partitionlev<-function(tree,suppo){

#if (is.null(tree$label)) tree$label<-tree$mean

xlkm<-length(suppo)/2
sollkm<-length(tree$val)           #solmujen lkm

#Find parents and number of leaves:

leafloc<-findleafs(tree$left,tree$right)

N<-matrix(0,sollkm,1)   #number of leaves in the tree whose root is i
p<-matrix(0,sollkm,1)   #parent
t<-sollkm
while (t>=1){
  if ((!is.na(leafloc[t])) && (leafloc[t]==1)){  #l(t)=0 eli ollaan lehdessa
   N[t]<-1
  }
  else if ((!is.na(leafloc[t])) && (leafloc[t]==0)){ #non-terminal node
     p[tree$left[t]]<-t       #p[t+1]<-t
     p[tree$right[t]]<-t      #p[endpoint(tree,t)+1]<-t
     N[t]<-N[tree$left[t]]+N[tree$right[t]]
  }
  t<-t-1
}

leafnum<-N[1]
recs<-matrix(0,leafnum,2*xlkm)
values<-matrix(0,leafnum,1)
frekv<-matrix(0,leafnum,1)

if (leafnum==1){ 
  recs<-suppo
  values<-tree$val[1]
}
else{
 ind<-1
 for (i in 1:sollkm){
   if ((!is.na(leafloc[i])) && (leafloc[i]==1)){  #i is leaf
     values[ind]<-tree$mean[i]
     recs[ind,]<-suppo
     frekv[ind]<-tree$nelem[i]
     j<-i
     while (p[j]>0){  #we are not in the root
         pare<-p[j]
         vari<-tree$vec[pare]
         split<-tree$val[pare]
         if (tree$left[pare]==j){  #i is left child
           if (split<recs[ind,2*vari]){ #if we have new restriction 
               recs[ind,2*vari]<-split
           }
         }
         else{  #i is right child
           if (split>recs[ind,2*vari-1]){ #if we have new restriction 
               recs[ind,2*vari-1]<-split
           }
         }
         j<-pare
     }
     ind<-ind+1
   }
 }
}

return(list(values=values,recs=recs,frekv=frekv))
}













partition.old<-function(tree,suppo=NULL)
{
#Finds the partition corresponding to an evaluation tree.

#tree is a binary tree, list(val,vec,left,right,...)
#supp is 2*xlkm-vector, support of density

#Result is list(values,recs)
#  values is recnum-vector, values of the estimate
#  recs is recnum*(2*xlkm)-matrix

if (is.null(suppo)) suppo<-tree$support
xlkm<-length(suppo)/2
sollkm<-length(tree$left)           #solmujen lkm

# Find parents and number of leaves:

leafloc<-findleafs(tree$left,tree$right)

N<-matrix(0,sollkm,1)   #number of leaves in the tree whose root is i
p<-matrix(0,sollkm,1)   #parent
t<-sollkm
while (t>=1){
  if ((!is.na(leafloc[t])) && (leafloc[t]==1)){  #l(t)=0 eli ollaan lehdessa
   N[t]<-1
  }
  else if ((!is.na(leafloc[t])) && (leafloc[t]==0)){ #non-terminal node
     p[tree$left[t]]<-t       #p[t+1]<-t
     p[tree$right[t]]<-t      #p[endpoint(tree,t)+1]<-t
     N[t]<-N[tree$left[t]]+N[tree$right[t]]
  }
  t<-t-1
}
                     
leafnum<-N[1]
recs<-matrix(0,leafnum,2*xlkm)
values<-matrix(0,leafnum,1)
frekv<-matrix(0,leafnum,1)

if (leafnum==1){ 
  recs<-suppo
  values<-tree$mean[1]
}
else{
 ind<-1
 for (i in 1:sollkm){
   if ((!is.na(leafloc[i])) && (leafloc[i]==1)){  #i is leaf
     values[ind]<-tree$mean[i]
     recs[ind,]<-suppo
     frekv[ind]<-tree$nelem[i]
     j<-i
     while (p[j]>0){  #we are not in the root
         pare<-p[j]
         vari<-tree$direc[pare]
         split<-tree$split[pare]
         if (tree$left[pare]==j){  #i is left child
           if (split<recs[ind,2*vari]){ #if we have new restriction 
               recs[ind,2*vari]<-split
           }
         }
         else{  #i is right child
           if (split>recs[ind,2*vari-1]){ #if we have new restriction 
               recs[ind,2*vari-1]<-split
           }
         }
         j<-pare
     }
     ind<-ind+1
   }
 }
}
# clean zeros out
finrecs<-matrix(0,leafnum,2*xlkm)
finvalues<-matrix(0,leafnum,1)
finfrekv<-matrix(0,leafnum,1)
notzero<-0
for (i in 1:leafnum){
  if (values[i]>0){
      notzero<-notzero+1
      finvalues[notzero]<-values[i]
      finrecs[notzero,]<-recs[i,]
      finfrekv[notzero]<-frekv[i]
  }
}
finvalues<-finvalues[1:notzero]
finrecs<-finrecs[1:notzero,]
finfrekv<-finfrekv[1:notzero]

return(list(values=finvalues,recs=finrecs,frekv=finfrekv))
}













partition<-function(et,grid=TRUE,zerorecs=FALSE)
{

d<-length(et$support)/2
nodenum<-length(et$left)
left<-et$left
right<-et$right
if (is.null(et$step)) et$step<-stepcalc(et$support,et$N)

if (grid){
  low<-matrix(0,nodenum,d)
  upp<-matrix(0,nodenum,d)
  for (i in 1:nodenum){
    for (j in 1:d){
      low[i,j]<-et$low[i,j]*et$step[j]+et$support[2*j-1]
      upp[i,j]<-et$upp[i,j]*et$step[j]+et$support[2*j-1]
    }
  }
}
else{
  low<-et$low
  upp<-et$upp
}

mean<-et$mean

ll<-leaflocs(left,right)
leafloc<-ll$leafloc
leafnum<-ll$leafnum

values<-matrix(0,leafnum,1)
recs<-matrix(0,leafnum,2*d)

efek<-0
i<-1
while (i<=leafnum){  
   node<-leafloc[i]

   if ((mean[node]>0) || (zerorecs)){

     efek<-efek+1

     values[efek]<-mean[node]
     
     for (dd in 1:d){
         recs[efek,2*dd-1]<-low[node,dd]
         recs[efek,2*dd]<-upp[node,dd]
     }
   }

   i<-i+1
}

values<-values[1:efek]
if (efek==1) recs<-matrix(recs[1:efek,],1,2*d)
else recs<-recs[1:efek,]

return(list(values=values,recs=recs,support=et$support))
}













pcf.greedy.kernel.new<-function(dendat,h,leaf=round(dim(dendat)[1]/2),minobs=NULL, 
itermax = 200, type="cpp"){
if (type=="old"){
  eva<-eval.greedy(dendat,leaf)
  pa<-partition(eva)
  fs<-fs.calc.parti(pa,dendat,h)
  pcf<-eva
  pcf$value<-fs
}

if (type=="prune"){
    bt<-densplit(dendat,minobs=minobs)
    treeseq<-prune(bt)
    eva<-eval.pick(treeseq,leaf=treeseq$leafs[1])  
    pa<-partition(eva)
    pcf<-eva
    fs<-fs.calc.parti(pa,dendat,h)
    pcf$value<-fs
}

if (type=="greedy"){
   pa<-densplitter(dendat,minobs=minobs)
   pcf<-list(down=pa$down,high=pa$high,grid=pa$grid,support=pa$support,recs=pa$recs)
   fs<-fs.calc.parti(pa,dendat,h)
   pcf$value<-fs
}

if (type=="cpp"){
   pcf<-densplitter2(dendat, minobs=minobs, neld = TRUE, itermax = 500)
   pcf$value <-fsCalcParti(pcf$recs, dendat, h)
}

return(pcf)
}

pcf.greedy.kernel<-function(dendat,h,leaf=round(dim(dendat)[1]/2),
minobs=NULL,type="greedy")
{
if (type=="old"){
  eva<-eval.greedy(dendat,leaf)
  pa<-partition(eva)
  fs<-fs.calc.parti(pa,dendat,h)
  pcf<-eva
  pcf$value<-fs
}

if (type=="prune"){
    bt<-densplit(dendat,minobs=minobs)
    treeseq<-prune(bt)
    eva<-eval.pick(treeseq,leaf=treeseq$leafs[1])  
    pa<-partition(eva)
    pcf<-eva
    fs<-fs.calc.parti(pa,dendat,h)
    pcf$value<-fs
}

if (type=="greedy"){
   pa<-densplitter(dendat,minobs=minobs)
   pcf<-list(down=pa$down,high=pa$high,grid=pa$grid,support=pa$support,recs=pa$recs)
   fs<-fs.calc.parti(pa,dendat,h)
   pcf$value<-fs
}

if (type=="dyadic"){
   pa<-densplitter(dendat,minobs=minobs,dyadic=TRUE)
   pcf<-list(down=pa$down,high=pa$high,grid=pa$grid,support=pa$support,recs=pa$recs)
   fs<-fs.calc.parti(pa,dendat,h)
   pcf$value<-fs
}

if (type=="cpp"){
   pa<-densplitter2(dendat,minobs=minobs)
   pcf<-list(down=pa$down,high=pa$high,grid=pa$grid,support=pa$support,recs=pa$recs)
   fs<-fs.calc.parti(pa,dendat,h)
   pcf$value<-fs
}

return(pcf)
}



pickseq<-function(treeseq,suppo,lsets=FALSE,invalue=FALSE,parvec=NULL)
{

if (is.null(parvec)) leaf<-treeseq$leafs
else leaf<-treeseq$leafs[parvec]
alfa<-treeseq$alfa[parvec]

alkm<-length(parvec)
for (inds in alkm:1){  # start with the oversmoothed estimate 
     leafnum<-leaf[inds]
   
     tree<-eval.pick(treeseq,leafnum)
     #pv<-partition(tree,suppo)
     #curtree<-profgene(pv$values,pv$recs,frekv=F,cvol=T,ccen=T,cfre=F)
     curtree<-proftree(tree)

     if (inds==alkm){
        if (alkm==1){
            treelist<-curtree
        }
        else{
           treelist=list(curtree)
        }
     }
     else{
        treelist=c(treelist,list(curtree))
     }
}

return(list(treelist=treelist,alfa=alfa,leaf=leaf))
}
picktreelev<-function(treeseq,leafnum){

tree<-treeseq$tree
#delnodes<-treeseq$delnodes
#delend<-treeseq$delend
leafs<-treeseq$leafs
indeksi<-detsi(leafs,leafnum)
endi<-treeseq$delnodeend[indeksi]
if (endi>0){       #if there is something to remove
  indeksit<-treeseq$delnodes[1:endi]
  re<-remnodes(tree$left,tree$right,indeksit)
  tree$left<-re$left
  tree$right<-re$right
}                                  
return(tree)
}
plotpartilev<-function(pa,dendat=NULL,restri=NULL,pch=21,col="blue")
{
recs<-pa$recs

if (!is.null(dendat)) plot(dendat,xlab="",ylab="",pch=pch)
else{
  xmin<-min(recs[,1])
  xmax<-max(recs[,2])
  ymin<-min(recs[,3])
  ymax<-max(recs[,4])
  plot(0,0,type="n",ylim=c(ymin,ymax),xlab="",ylab="",xlim=c(xmin,xmax))
}

if (is.null(dim(recs))) len<-1 else len<-dim(recs)[1]

if (is.null(restri)) restric<-matrix(1,len,1)
else{
  restric<-matrix(0,len,1)
  restrilen<-length(restri)
  for (i in 1:restrilen){
     restric[restri[i]]<-1
  }
}

i<-1
while (i<=len){

 if (restric[i]==1){

    if (len==1){
      x<-c(recs[1],recs[1],recs[2],recs[2])
      y<-c(recs[3],recs[4],recs[4],recs[3])
    }
    else{
      x<-c(recs[i,1],recs[i,1],recs[i,2],recs[i,2])
      y<-c(recs[i,3],recs[i,4],recs[i,4],recs[i,3])
    }

    if (pa$values[i]==0) colo<-NA else colo=col  
    polygon(x,y,col=colo)

    #lines(c(recs[i,1],recs[i,1]),c(recs[i,3],recs[i,4]))
    #lines(c(recs[i,1],recs[i,2]),c(recs[i,4],recs[i,4]))
    #lines(c(recs[i,2],recs[i,2]),c(recs[i,4],recs[i,3]))
    #lines(c(recs[i,2],recs[i,1]),c(recs[i,3],recs[i,3]))
 }

 i<-i+1
}

if (!is.null(dendat)){

  points(dendat,pch=pch)

  suppor<-supp(dendat)

  x<-c(suppor[1],suppor[1],suppor[2],suppor[2])
  y<-c(suppor[3],suppor[4],suppor[4],suppor[3])
  polygon(x,y)
  #lines(c(suppor[1],suppor[1]),c(suppor[3],suppor[4]))
  #lines(c(suppor[1],suppor[2]),c(suppor[4],suppor[4]))
  #lines(c(suppor[2],suppor[2]),c(suppor[4],suppor[3]))
  #lines(c(suppor[2],suppor[1]),c(suppor[3],suppor[3]))
}

}






plotparti<-function(pa,d1=NULL,d2=NULL,
dendat=NULL,restri=NULL,pch=21,support=pa$support,col="black",cex.axis=1)
{
# support=NULL

recs<-pa$recs
if (!is.null(d1)) recs<-recs[,c(2*d1-1,2*d1,2*d2-1,2*d2)]

xmin<-min(recs[,1])
xmax<-max(recs[,2])
ymin<-min(recs[,3])
ymax<-max(recs[,4])

if (!is.null(dendat)) plot(dendat,xlab="",ylab="",pch=pch,cex.axis=cex.axis)
else plot(0,0,type="n",ylim=c(ymin,ymax),xlab="",ylab="",xlim=c(xmin,xmax),
     cex.axis=cex.axis)

len<-dim(recs)[1]

if (is.null(restri)) restric<-matrix(1,len,1)
else{
  restric<-matrix(0,len,1)
  restrilen<-length(restri)
  for (i in 1:restrilen){
     restric[restri[i]]<-1
  }
}

i<-1
while (i<=len){

 if (restric[i]==1){
    lines(c(recs[i,1],recs[i,1]),c(recs[i,3],recs[i,4]),col=col)
    lines(c(recs[i,1],recs[i,2]),c(recs[i,4],recs[i,4]),col=col)
    lines(c(recs[i,2],recs[i,2]),c(recs[i,4],recs[i,3]),col=col)
    lines(c(recs[i,2],recs[i,1]),c(recs[i,3],recs[i,3]),col=col)
 }

 i<-i+1
}

if (!is.null(support)){
  lines(c(support[1],support[1]),c(support[3],support[4]),col=col)
  lines(c(support[1],support[2]),c(support[4],support[4]),col=col)
  lines(c(support[2],support[2]),c(support[4],support[3]),col=col)
  lines(c(support[2],support[1]),c(support[3],support[3]),col=col)
}

}







preprocess<-function(ssr,left,right,mean)
{
#ssr=excess mass

nodlkm<-length(ssr)
ssralip<-matrix(0,nodlkm,1)
#label<-matrix(1,nodlkm,1)

# Muodostetaan vektori, jossa kunkin solmun korkeus
kork<-matrix(0,nodlkm,1)
kork[1]<-1    #juuren korkeus on 1
i<-1
while (i<=nodlkm){
  if (right[i]>0){   #jos ei olla lehdessa
    kork[left[i]]<-kork[i]+1       #vasen lapsi yhta korkeammalla
    kork[right[i]]<-kork[i]+1      #oikea lapsi yhta korkeammalla
  }
  i<-i+1
}
tasolkm<-max(kork)     #tasojen lkm

k<-tasolkm             #aloitetaan korkeimmalta tasolta
while (k>0){
  i<-1
  while (i<=nodlkm){  #nodlkm
    if (kork[i]==k){
      if (right[i]==0){    #jos ollaan lehdessa
         ssralip[i]<-ssr[i] #alipuun ssr=solmun itsensa ssr
      }
      else if ((left[left[i]]==0) && (left[right[i]]==0)){
         minnu<-min(ssralip[left[i]],ssralip[right[i]])
         minnu<-min(minnu,ssr[i]) 
         if (ssralip[left[i]]==minnu){  #remove right
              ssralip[i]<-ssralip[left[i]]
              mean[right[i]]<-0
              left[right[i]]<-0
              right[right[i]]<-0
          }
          else 
             if (ssralip[right[i]]==minnu){  #remove left
                ssralip[i]<-ssralip[right[i]]
                mean[left[i]]<-0
                left[left[i]]<-0
                right[left[i]]<-0
             }
             else{
                left[i]<-0 
                right[i]<-0
                ssralip[i]<-ssr[i]
             }
      }  
      else{
         minnu<-min(ssralip[left[i]],ssralip[right[i]])
         minnu<-min(minnu,ssralip[left[i]]+ssralip[right[i]])
         minnu<-min(minnu,ssr[i])  
         if (minnu==ssralip[left[i]]+ssralip[right[i]]){
            ssralip[i]<-ssralip[left[i]]+ssralip[right[i]]
         }
         else
            if (ssralip[left[i]]==minnu){  #remove right
              ssralip[i]<-ssralip[left[i]]
              mean[right[i]]<-0
              left[right[i]]<-0
              right[right[i]]<-0
            }
            else 
              if (ssralip[right[i]]==minnu){  #remove left
                ssralip[i]<-ssralip[right[i]]
                mean[left[i]]<-0
                left[left[i]]<-0
                right[left[i]]<-0
              }
              else{
                left[i]<-0 
                right[i]<-0
                ssralip[i]<-ssr[i]
              }
      #
      }
    }
    i<-i+1
  }
  k<-k-1
}

return(list(right=right,left=left,S=ssralip,mean=mean))
}



profdelt<-function(treeseq,leafnum,suppo,frekv=NULL,cvol=FALSE,ccen=FALSE,cfre=FALSE){
#Profiles an adaptive histogram

tree<-NULL
#tree<-picktree(treeseq,leafnum)

# Tehdaan binaaripuusta paloittain vakio
pv<-partition(tree,suppo)
recs<-pv$recs
values<-pv$values
#                               
pg<-profgene(values,recs,frekv,cvol,ccen,cfre)
#
return(pg)
}


prunelev<-function(bt,lambda=NULL,n=NULL){

len<-length(bt$mean)
bt$mean<-matrix(1,len,1)

if (!is.null(lambda)) bt$ssr<-exmavec(bt$volume,bt$nelem,n,lambda)

ini1<-preprocess(bt$ssr,bt$left,bt$right,bt$mean)
bt$S<-t(ini1$S)
bt$mean<-t(ini1$mean)
bt$left<-ini1$left
bt$right<-ini1$right


#ini<-initial(bt$ssr,bt$left,bt$right)
#bt$left<-ini$left
#bt$right<-ini$right

treeseq<-pruseqlev(bt)

return(treeseq)
}


prune<-function(et)
{

ini<-initial(et$ssr,et$left,et$right)
et$left<-ini$left
et$right<-ini$right
treeseq<-pruseq(et)

treeseq$support<-et$support
return(treeseq)
}
pruseqlev<-function(tree){

S<-tree$S
lu<-luo(tree)
p<-lu$p
G<-lu$G
g<-lu$g
N<-lu$N

left<-tree$left
right<-tree$right

alknodlkm<-length(tree$vec)  #solmujen lkm
alfalkm<-alknodlkm           #alfojen lkm <= lehtien lkm <= solmujen lkm

delnodebeg<-matrix(0,alfalkm,1)
delnodeend<-matrix(0,alfalkm,1)
delWbeg<-matrix(0,alfalkm,1)
delWend<-matrix(0,alfalkm,1)
leafs<-matrix(0,alfalkm,1)
alfa<-matrix(0,alfalkm,1)
loglik<-matrix(0,alfalkm,1)

delnodes<-matrix(0,alknodlkm,1)
delW<-matrix(0,alknodlkm,1)

alphamin<-0

if (N[1]==1){    #tree on triviaali 

  delnodes<-NULL
  delW<-NULL
  delnodebeg<-c(0)
  delnodeend<-c(0)
  delWbeg<-c(0)
  delWend<-c(0)
  leafs<-c(1)
  alfa<-c(0)
  loglik<-c(tree$ssr[1])
  tulos<-list(tree=tree,delnodes=delnodes,delW=delW,
              delnodebeg=delnodebeg,delnodeend=delnodeend,
              delWbeg=delWbeg,delWend=delWend,
              leafs=leafs,alfa=alfa,loglik=loglik)
}
else{ 
 k<-1           #tulee kertomaan alfojen lkm:n +1
 #remnodenum<-0  #poistettavien haarojen lkm
 #j<-0           #poistolaskuri

 delnodebeg[k]<-1
 delWbeg[k]<-1

 alpha<-alphamin

 while (N[1]>1){     #jos puu ei triviaali

   if (omaver(alpha,G[1])){    #(G[1]>alpha){
      leafs[k]<-N[1]
      alfa[k]<-alpha
      loglik[k]<-S[1] 
      alpha<-G[1]

      k<-k+1   #siirrytaan uuteen alfaan
      #j<-0     #poistolaskuri nollataan
 
      delnodebeg[k]<-delnodeend[k-1]+1 
      delnodeend[k]<-delnodeend[k-1]
      delWbeg[k]<-delWend[k-1]+1 
      delWend[k]<-delWend[k-1]


   }

   # Haetaan seuraavaksi solmu t, jonka alapuolelta voidaan katkaista.
   t<-1                        #aloitetaan juuresta
   while (omaver(G[t],g[t])){ #(G[t]<g[t]){  #jatk. kunnes g saav. miniminsa 
      if (omaver(G[left[t]],G[right[t]]))
                         #(omasam(G[t],G[left[t]])) 
           t<-left[t]                 #mennaan vasemmalle
      else t<-right[t]                #mennaan oikealle
   }

   delnodeend[k]<-delnodeend[k]+1
   #remnodenum<-remnodenum+1
   #j<-j+1
   delnodes[delnodeend[k]]<-t   

   # Tehdaan t:sta lehti
   N[t]<-1
   S[t]<-tree$ssr[t]
   G[t]<-NA                    #infty
   g[t]<-NA
   # Palataan juureen paivittaen N, S, g ja G.
   while (t>1){
      t<-p[t]
      ######################
      #S[t]<-S[left[t]]+S[right[t]]
         minnu<-min(S[left[t]],S[right[t]])
         minnu<-min(minnu,S[left[t]]+S[right[t]])
         minnu<-min(minnu,tree$ssr[t])  
         if (minnu==S[left[t]]+S[right[t]]){
            S[t]<-S[left[t]]+S[right[t]]
         }
         else
            if (S[left[t]]==minnu){  #remove right
              S[t]<-S[left[t]]
              delWend[k]<-delWend[k]+1
              delW[delWend[k]]<-right[t]
              #
              delnodeend[k]<-delnodeend[k]+1
              delnodes[delnodeend[k]]<-right[t] 
              #
              left[right[t]]<-0
              right[right[t]]<-0
              N[right[t]]<-1
              N[right[t]]<-1
              G[right[t]]<-NA
              G[right[t]]<-NA
              g[right[t]]<-NA
              g[right[t]]<-NA
            }
            else 
              if (S[right[t]]==minnu){  #remove left
                S[t]<-S[right[t]]
                delWend[k]<-delWend[k]+1
                delW[delWend[k]]<-left[t]
                #
                delnodeend[k]<-delnodeend[k]+1
                delnodes[delnodeend[k]]<-left[t] 
                #
                left[left[t]]<-0
                right[left[t]]<-0
                N[left[t]]<-1
                N[left[t]]<-1
                G[left[t]]<-NA
                G[left[t]]<-NA
                g[left[t]]<-NA
                g[left[t]]<-NA
              }
              else{
                S[t]<-tree$ssr[t]
                left[t]<-0 
                right[t]<-0
              }
      ###############################################
      if (left[t]==0){
                N[t]<-1
                G[t]<-NA
                g[t]<-NA
      }
      else{
         N[t]<-N[left[t]]+N[right[t]]
         g[t]<-(tree$ssr[t]-S[t])/(N[t]-1)
         G[t]<-omamindelt(g[t],omamindelt(G[left[t]],G[right[t]]))
      }

   }  #while t>1
 } #while N[t]>1

leafs[k]<-N[1]
alfa[k]<-alpha
loglik[k]<-S[1]  #palataan -log-likelista log-likeliin
alpha<-G[1]

leafs<-leafs[1:k]      #(k-1) kertoo alfojen maaran
alfa<-alfa[1:k]
loglik<-loglik[1:k]

delnodebeg<-delnodebeg[1:k]
delWbeg<-delWbeg[1:k]
delnodeend<-delnodeend[1:k]
delWend<-delWend[1:k]

delnodes<-delnodes[1:delnodeend[k]]
delW<-delW[1:delWend[k]]

tulos<-list(tree=tree,delnodes=delnodes,delW=delW,
            delnodebeg=delnodebeg,delnodeend=delnodeend,
            delWbeg=delWbeg,delWend=delWend,
            leafs=leafs,alfa=alfa,loglik=loglik)

}

return(tulos)
}







pruseq<-function(tree){
#Forms a sequence of trees, which minimize likelihood-complexity
#criterion
#
#tree on list(val,vec,mean,nelem,ssr,left,right), ks densplit
#
#Result on list(tree,subtree,leafs,alfa,loglik)
#-tree on alkup. puu,
#-subtree on alfalkm*nodelkm-matriisi: osapuu annettu niitten solmujen 
#  indekseina, joista alkavat "tree":n alipuut eivat kuulu ko. osapuuhun.
#  nodelkm on yli alipuitten otettu maksimi poistettavien indeksien
#  maarasta
#-leafs,alfa,loglik are alfalkm-vectors: 
#    sis alipuuun lehtien lkm:n, alfan ja alipuun uskottavuuden

#alklehlkm<-leafnum(tree,1)     #lehtien lkm
alknodlkm<-length(tree$left)    #solmujen lkm
#nodelkm<-alknodlkm-alklehlkm
alfalkm<-alknodlkm             #alfojen lkm <= lehtien lkm <= solmujen lkm

delnodes<-matrix(0,alknodlkm,1)
delend<-matrix(0,alfalkm,1)      
leafs<-matrix(0,alfalkm,1)
alfa<-matrix(0,alfalkm,1)
loglik<-matrix(0,alfalkm,1)

# Initializing:

leafloc<-findleafs(tree$left,tree$right)

N<-matrix(0,alknodlkm,1)   #number of leaves in the tree whose root is i
R<--tree$ssr               #R(i) on noden i ssr eli -log likeli
p<-matrix(0,alknodlkm,1)   #parent
S<-matrix(0,alknodlkm,1)   #sen alipuun ssr -(log-likeli), jonka juuri i
g<-matrix(0,alknodlkm,1)   #(R(i)-S(i))/(N(i)-1), R(i) on noden i ssr
G<-matrix(0,alknodlkm,1)   #min{g(t),G(l(t)),G(r(t))}
t<-alknodlkm
while (t>=1){
  if (!is.na(leafloc[t]) && (leafloc[t]==1)){  #l(t)=0 eli ollaan lehdessa 
   N[t]<-1
   S[t]<--tree$ssr[t]                        
   G[t]<-NA                 #\infty
  }
  else  if (!is.na(leafloc[t])){
     p[tree$left[t]]<-t       #p[t+1]<-t
     p[tree$right[t]]<-t      #p[endpoint(tree,t)+1]<-t
     S[t]<-S[tree$left[t]]+S[tree$right[t]]  #S[t+1]+S[endpoint(tree,t+1)+1]    
     N[t]<-N[tree$left[t]]+N[tree$right[t]]
     g[t]<-(R[t]-S[t])/(N[t]-1)
     G[t]<-omamindelt(g[t],omamindelt(G[tree$left[t]],G[tree$right[t]]))
  }
  t<-t-1
}

alphamin<-0

if (N[1]==1){    #tree on triviaali 

  delnodes<-NULL
  delend<-c(0)
  leafs<-c(1)
  alfa<-c(0)
  tulos<-list(tree=tree,delnodes=delnodes,delend=delend,leafs=leafs,alfa=alfa,
              loglik=loglik)
}
else{ 
 k<-1  #tulee kertomaan alfojen lkm:n +1
 delend[k]<-0
 j<-0  #poistolaskuri
 remnodenum<-0    #poistettavien haarojen lkm
 alpha<-alphamin
 while (N[1]>1){     #jos puu ei triviaali
   if (omaver(alpha,G[1])){    #(G[1]>alpha){
      leafs[k]<-N[1]
      alfa[k]<-alpha
      loglik[k]<--S[1]  #palataan -log-likelista log-likeliin
      alpha<-G[1]
      if (k>=2) delend[k]<-delend[k-1]+j
      k<-k+1   #siirrytaan uuteen alfaan
      j<-0     #poistolaskuri nollataan
   }
   #Haetaan seuraavaksi solmu t, jonka alapuolelta voidaan katkaista.
   t<-1                        #aloitetaan juuresta
   while (omaver(G[t],g[t])){ #(G[t]<g[t]){  #jatk. kunnes g saav. miniminsa 
      if (omaver(G[tree$left[t]],G[tree$right[t]]))
                         #(omasam(G[t],G[tree$left[t]])) 
         t<-tree$left[t]                 #mennaan vasemmalle
      else t<-tree$right[t]              #mennaan oikealle
   }
   remnodenum<-remnodenum+1
   delnodes[remnodenum]<-t   
   j<-j+1
   #Tehdaan t:sta lehti
   N[t]<-1
   S[t]<-R[t]
   G[t]<-NA                    #infty
   #Palataan juureen paivittaen N, S, g ja G.
   while (t>1){
      t<-p[t]
      S[t]<-S[tree$left[t]]+S[tree$right[t]]
      N[t]<-N[tree$left[t]]+N[tree$right[t]]
      g[t]<-(R[t]-S[t])/(N[t]-1)
      G[t]<-omamindelt(g[t],omamindelt(G[tree$left[t]],G[tree$right[t]]))
   }
 }

leafs[k]<-N[1]
alfa[k]<-alpha
loglik[k]<--S[1]  #palataan -log-likelista log-likeliin
alpha<-G[1]
if (k>=2) delend[k]<-delend[k-1]+j
k<-k+1

leafs<-leafs[1:(k-1)]      #(k-1) kertoo alfojen maaran
alfa<-alfa[1:(k-1)]
loglik<-loglik[1:(k-1)]
delnodes<-delnodes[1:remnodenum]
delend<-delend[1:(k-1)]
tulos<-list(tree=tree,delnodes=delnodes,delend=delend,leafs=leafs,alfa=alfa,
loglik=loglik)
}

return(tulos)
}







remnodes<-function(left,right,list){
#Removes branches from a tree
#
#list is a vector, branches whose root is mentioned in the list
#  will be removed
#
num<-length(list)
for (i in 1:num){
  left[list[i]]<-0
  right[list[i]]<-0
}
return(list(left=left,right=right))
}


rf2tree.old<-function(forest,suppo)
{
d<-length(suppo)/2

nr<-length(forest$ndbigtree)     #number of trees in the forest

nrtreemap<-length(forest$treemap)
map<-matrix(0,nrtreemap,1)
infopointer<-matrix(0,nrtreemap,1)
rootinfo<-matrix(0,nr,1)

# create infopointer
laskuri<-1
for (ij in 1:nrtreemap){
     if (forest$treemap[ij]==2) laskuri<-laskuri+1
     if (forest$treemap[ij]!=0){
            infopointer[ij]<-laskuri
            laskuri<-laskuri+1
     }
}
# create rootinfo
rootinfo[1]<-1
cusu<-0
ii<-1
while (ii<=nr){
   rootinfo[ii]<-cusu+1
   cusu<-cusu+forest$ndbigtree[ii]
   ii<-ii+1
}

totalrunner<-0
glob<-1

#####################################################################
while (glob <= nr){

# build a small tree

maxnrnodes<-2*forest$ndbigtree[glob]

left<-matrix(0,maxnrnodes,1)
right<-matrix(0,maxnrnodes,1)
val<-matrix(0,maxnrnodes,1)
vec<-matrix(0,maxnrnodes,1)
mean<-matrix(0,maxnrnodes,1)
nelem<-matrix(0,maxnrnodes,1)
low<-matrix(0,maxnrnodes,d)
upp<-matrix(0,maxnrnodes,d)

# create root

{
if (forest$treemap[1]!=0){
   #left[1]<-2
   #right[1]<-3

   locu<-rootinfo[glob]

   val[1]<-forest$upper[locu]
   vec[1]<-forest$mbest[locu]
   mean[1]<-forest$avnode[locu]
   for (si in 1:d){
      low[1,si]<-suppo[2*si-1]
      upp[1,si]<-suppo[2*si]
   }

   #map[1]<-2
   #map[2]<-3

   nodesleft<-1
   treeind<-1    #3
   
   prevbeg<-totalrunner #0 #1  #beg of previous level
   prevend<-totalrunner #0 #2  #end of previous level
   curbeg<-totalrunner+1  #1  #3
   curend<-totalrunner+2  #2  #6

   #totalrunner<-2
}
else{
   nodesleft<-0
   treeind<-1
}
}

while (nodesleft==1){

   nodesleft<-0
   prevlkm<-0
   curlkm<-curend-curbeg+1
   muisti<-matrix(0,curlkm,1)

   i<-curbeg

   while (i <= (curbeg+(curend-curbeg+1)/2-1)){
       indl<-curbeg+(2*(i-curbeg+1)-1)-1
       indr<-curbeg+2*(i-curbeg+1)-1
       #ind<-i-curbeg+1
       #parind<-prevbeg+ind-1

       if (forest$treemap[indl]!=0){
           
           if (prevbeg==prevend) curparent<-1 
           else{ 
                  parind<-vanh[indl-curbeg+1]    
                  parindmap<-prevbeg+parind-1           
                  curparent<-map[parindmap]
           }

           node<-treeind+1
           left[curparent]<-node

           locu<-infopointer[indl]

           val[node]<-forest$upper[locu]
           vec[node]<-forest$mbest[locu]
           mean[node]<-forest$avnode[locu]
           for (si in 1:d){
              low[node,si]<-low[curparent,si]
              upp[node,si]<-upp[curparent,si]
           }
           split<-vec[curparent]
           upp[node,split]<-val[curparent]

           map[indl]<-node
           treeind<-treeind+1
           prevlkm<-prevlkm+1
           nodesleft<-1
       
        #if (forest$treemap[indr]!=0)

           node<-treeind+1
           right[curparent]<-node

           locu<-infopointer[indr]

           val[node]<-forest$upper[locu]
           vec[node]<-forest$mbest[locu]
           mean[node]<-forest$avnode[locu]
           for (si in 1:d){
              low[node,si]<-low[curparent,si]
              upp[node,si]<-upp[curparent,si]
           }
           split<-vec[curparent]
           low[node,split]<-val[curparent]

           map[indr]<-node 
           treeind<-treeind+1
           prevlkm<-prevlkm+1
           nodesleft<-1
       }
       i<-i+1
   }

   prevbeg<-curbeg
   prevend<-curend
   curbeg<-curend+1
   curend<-curbeg+2*prevlkm-1

   vanh<-matrix(0,curend-curbeg+1,1)
   liuk<-0
   sep<-1
   while (sep <= (prevend-prevbeg+1)){
        if (forest$treemap[prevbeg+sep-1]!=0){
             liuk<-liuk+1
             vanh[2*liuk-1]<-sep
             vanh[2*liuk]<-sep
        }
        sep<-sep+1
   }

}

while ((forest$treemap[curend]==0) && (glob<nr)) curend<-curend+1
totalrunner<-curend-1

left<-left[1:treeind]
right<-right[1:treeind]
val<-val[1:treeind]
vec<-vec[1:treeind]
mean<-mean[1:treeind]
nelem<-nelem[1:treeind]
low<-low[1:treeind,]
upp<-upp[1:treeind,]

trnew<-list(val=val,vec=vec,mean=mean,nelem=nelem,
left=left,right=right,low=low,upp=upp)

# small tree ready
{
if (glob>1) tr<-treeadd(tr,trnew,d)
else tr<-trnew
}

glob<-glob+1

}

return(tr)
}









rf2tree<-function(rf,support,N=NULL)
{
forest<-rf$forest
d<-length(support)/2

if (is.null(N)) N<-rep(60,d)

nr<-length(forest$ndbigtree)     #number of trees in the forest
glob<-1
while (glob <= nr){

   # build trnew

   left<-forest$leftDaughter[,glob]
   right<-forest$rightDaughter[,glob]
   direc<-forest$bestvar[,glob]
   direc[(forest$bestvar[,glob]==0)]<-NA
   mean<-forest$nodepred[,glob]
   nodenum<-length(forest$xbestsplit[,glob])
   nelem<-rep(1,nodenum)
   split<-matrix(NA,nodenum,1)
   for (i in 1:nodenum){
       vec<-direc[i]
       ala<-support[2*vec-1]
       yla<-support[2*vec]
       splitti<-round(N[vec]*(forest$xbestsplit[i,glob]-ala)/(yla-ala))
       split[i]<-min(max(splitti,1),(N[vec]-1))
   }
   split[(forest$bestvar[,glob]==0)]<-NA

   trnew<-list(split=split,direc=direc,mean=mean,nelem=nelem,
   left=left,right=right,#low=low,upp=upp,
   N=N,support=support)

   lu<-lowupp(trnew)
   trnew$low<-lu$low
   trnew$upp<-lu$upp
   trnew$volume<-rep(1,nodenum)

   if (glob>1) tr<-treeadd(tr,trnew) else tr<-trnew

   glob<-glob+1
}

dh<-downhigh(tr)
tr$down<-dh$down
tr$high<-dh$high   
tr$value<-dh$value
tr$infopointer<-dh$infopointer

return(tr)
}
riskesti<-function(treeseq,n){
#Estimates risk for every alpha
#
#treeseq is list(tree,leafs,alfa,...)
#  tree is list(volum,nelem,...)
#n is the sample size
#
#Returns an alphalkm-vector
#
tr<-treeseq$tree
left<-tr$left
right<-tr$right
alfalkm<-length(treeseq$alfa)
toremove<-treeseq$delnodes
if (dim(t(toremove))[1]==1) maxrem<-1 else maxrem<-length(toremove[1,])
#mita jos toremove on skalaari?    
#
inum<-length(tr$vec)
ykk<-rep(1,inum)
nelem<-tr$nelem
volum<-tr$volum
#  estimated risk is sum of the info over leafs, info is vector which
#  we have to sum over leafs 
info<-nelem*(ykk-nelem*(1+1/n))/volum  #/n^2
#
risks<-matrix(0,alfalkm,1)
risks[1]<-leafsum(info,root=1,left,right) #kun alpha=0, ei ole poist mitaan
cursum<-risks[1]
for (i in 1:alfalkm){
    if (maxrem==1){ 
       poista<-toremove[i]
       sumsubtree<-leafsum(info,root=poista,left,right)
       cursum<-cursum-sumsubtree+info[poista]
       left[poista]<-0
       right[poista]<-0
    }
    else{
      j<-1
      while ((j<=maxrem) && (toremove[i,j]>0)){
         poista<-toremove[i,j]
         sumsubtree<-leafsum(info,root=poista,left,right)
         cursum<-cursum-sumsubtree+info[poista]
         left[poista]<-0
         right[poista]<-0  
         j<-j+1
      }
    }
    risks[i]<-cursum
}
return(t(risks))
}




roundlnum<-function(lseq,wanted)
{
len<-length(lseq)
runner<-1
cur<-lseq[runner]
while ((cur>wanted) && (runner<len)){
    runner<-runner+1
    cur<-lseq[runner]
}
if (runner==1){
   approwanted<-cur
}
else{
  if ((wanted-cur)<=(lseq[runner-1]-wanted)){
     approwanted<-cur
  }
  else{
     approwanted<-lseq[runner-1]
  }
}

return(approwanted)
}    
scaspa<-function(treeseq,bind,eind)
{

alkm<-eind-bind+1
modelkm<-matrix(0,alkm,1)

j<-1
for (i in bind:eind){
     leafnum<-treeseq$leafs[i]
     tree<-eval.pick(treeseq,leafnum)
     pv<-partition(tree)
     values<-pv$values  
     recs<-pv$recs
     if (length(values==1)) modelkm[j]<-1
     else{
       pg<-profgene(values,recs)
       parents<-pg$parent
       mlkm<-moodilkm(parents)
       modelkm[j]<-mlkm$lkm
     }
     j<-j+1
}

leafnums<-treeseq$leafs[bind:eind]
alfas<-treeseq$alfa[bind:eind]
return(list(moodilkm=t(modelkm),alfas=alfas,leafnums=leafnums))
}
simesi3d<-function(n,p1,p2,p3,s1,s2,s3,siemen,noisedim=0){
#muodostaa n*2 havaintomatriisin
#
#p1+p2+p3+p4=1
#
#tiheysfunktio on paloittain vakio 4 palassa:
#I: [0,split1]*[0,1]*[0,1]
#II: [split1,1]*[0,split2]*[0,1] 
#III: [split1,1]*[split2,1]*[0,split3]
#IV: [split1,1]*[split2,1]*[split3,1]
#
#p1 on palan I tn. ja p2 on palan II tn.,...
#p1=f1*s1, p2=(1-s1)*s2*f2, p3=(1-s1)*(1-s2)*s3*f3, 
#p4=(1-s1)*(1-s2)*(1-s3)*f4
#kokeillaan esim. p1=0.1, p2=0.2, p3=0.3, p4=0.4 
#
set.seed(siemen)
d<-noisedim+3
x<-matrix(0,n,d)     #havaintomatriisi
i<-1
while (i<=n){         
  U<-runif(1)        #arpoo mihin palaan havainto tulee
  U1<-runif(1)
  U2<-runif(1)
  U3<-runif(1)
  if (U<p1){         #tn:lla p1 palaan I
    x[i,1]<-s1*U1    #x-koord skaalataan valiin [0,split1]
    x[i,2]<-U2       #y-koordinaatti valiin [0,1]
    x[i,3]<-U3       #z-koordinaatti valiin [0,1]
  }
  else{ 
    if ((U>=p1) && (U<p1+p2)){     #pala II
       x[i,1]<-s1+(1-s1)*U1  #x-koord skaalataan valiin [split1,1]
       x[i,2]<-s2*U2         #y-koord skaalataan valiin [0,split2]
       x[i,3]<-U3            #z-koordinaatti valiin [0,1] 
    }
    else{
       if ((U>=p1+p2) && (U<p1+p2+p3)){    #pala III
          x[i,1]<-s1+(1-s1)*U1   #x-koord skaalataan valiin [split1,1]  
          x[i,2]<-s2+(1-s2)*U2   #y-koord skaalataan valiin [split2,1]
          x[i,3]<-s3*U3          #z-koord skaalataan valiin [0,split3]
       }
       else{
          x[i,1]<-s1+(1-s1)*U1   #x-koord skaalataan valiin [split1,1]  
          x[i,2]<-s2+(1-s2)*U2   #y-koord skaalataan valiin [split2,1]
          x[i,3]<-s3+(1-s3)*U3   #z-koord skaalataan valiin [split3,1]
       }
    }
  }
  if (noisedim>0){
     nd<-3
     while (nd<=(3+noisedim)){
         x[i,nd]<-runif(1)
         nd<-nd+1
     }
  }
i<-i+1
}
return(x)
}




simesi<-function(n,p1,p2,s1,s2,siemen,noisedim=0){
#muodostaa n*2 havaintomatriisin
#
#p1+p2+p3=1
#split1, split2 in (0,1)
#
#tiheysfunktio on paloittain vakio 3 palassa:
#I: [0,split1]*[0,1], II: [split1,1]*[0,split2], III: [split1,1]*[split2,1]
#p1 on palan I tn. ja p2 on palan II tn., 
#p1=f1*s1, p2=(1-s1)*s2*f2, p3=(1-s1)*(1-s2)*f3 
#kokeillaan esim. p1=0.1, p2=0.3, p3=0.6 
#
#tulos: 
#val: 0.5 NA 0.5 NA NA, 
#vec: 1   NA 2   NA NA,
#mean: havlkm/(n*tilavuus), likimain: 1  2*p1 2*(p2+p3) 4*p2 4*p3 
#nelem:  n  p1*n  (p1+p2)*n  p2*n  p3*n
#ssr: havlkm*log(mean) = havlkm*log(havlkm/(n*tilavuus)), 
#    likimain: n*log(1)=0      p1*n*log(2*p1)   (p1+p2)*n*log(2*(p2+p3))
#              p2*n*log(4*p2)  p3*n*log(4*p3) 
#
#p1<-0.1
#p2<-0.3
#p3<-0.6
#s1<-0.75
#s2<-0.5
#values<-c(p1/s1,p2/((1-s1)*s2),p3/((1-s1)*(1-s2)))
#recs<-matrix(0,3,4)
#recs[1,]<-c(0,s1,0,1)
#recs[2,]<-c(s1,1,0,s2)
#recs[3,]<-c(s1,1,s2,1)
#koe<-drawgene(values,recs,plkm=30)
#persp(koe$x,koe$y,koe$z,phi=30,theta=60) 
#
set.seed(siemen)
d<-noisedim+2
x<-matrix(0,n,d)     #havaintomatriisi
i<-1
while (i<=n){         
  U<-runif(1)        #arpoo mihin palaan havainto tulee
  U1<-runif(1)
  U2<-runif(1)
  if (U<p1){         #tn:lla p1 palaan I
    x[i,1]<-s1*U1   #x-koord skaalataan valiin [0,split1]
    x[i,2]<-U2       #y-koordinaatti valiin [0,1]
  }
  else{ if (U>p1+p2){          #tn:lla p3 palaan III
            x[i,1]<-s1+(1-s1)*U1   #x-koord skaalataan valiin [split1,1]  
            x[i,2]<-s2+(1-s2)*U2   #y-koord skaalataan valiin [split2,1]

          }
           else{                 #tn:lla p2 palaan II
             x[i,1]<-s1+(1-s1)*U1  #x-koord skaalataan valiin [split1,1]
             x[i,2]<-s2*U2      #y-koord skaalataan valiin [0,split2]
           }
  }
  if (noisedim>0){
     nd<-3
     while (nd<=(2+noisedim)){
         x[i,nd]<-runif(1)
         nd<-nd+1
     }
  }
i<-i+1
}
return(x)
}




simfssk<-function(n,noisedim,siemen){
#tekee n*d, d=2+noisedim, datamatriisin
#3 moodia, (c,0), (-c,3), (-c,-3)
#
d<-2+noisedim
hajo<-1
noisehajo<-sqrt(7)
c<-3^(3/2)/2
set.seed(siemen)
data<-matrix(rnorm(d*n),,d)   #n*d matriisi, valkoista kohinaa
data[,1:2]<-hajo*data[,1:2]
if (noisedim>0) data[,3:d]<-noisehajo*data[,3:d]
i<-1
while (i<=n){
  mu<-matrix(0,1,d)       #moodin keskipiste
  ehto<-runif(1)
  if (ehto<1/3){          #sekoitteiden painot samat
         mu[1,1]<-0
         mu[1,2]<-c
  } 
  else if (ehto>2/3){
         mu[1,1]<-3
         mu[1,2]<--c
  }
  else{
         mu[1,1]<--3
         mu[1,2]<--c
  }
  data[i,]<-data[i,]+mu
  i<-i+1
}
return(data)
}





















des12<-function(n,lkm,m,siemen){
#tekee n*6 data-matriisin, eli d=6
#2 moodia, m maarittaa moodien etaisyydet
#moodit tyyppia (m,m,0,...,0) ja (0,...,0) missa "lkm" kpl:tta m:ia. 
#
d<-6
set.seed(siemen)
data<-matrix(rnorm(d*n),,d)   #n*d matriisi, valkoista kohinaa
i<-1
while (i<=n){
  mu<-matrix(0,1,d)       #moodin keskipiste
  ehto<-runif(1)
  if (ehto>1/2){          #sekoitteiden painot samat
     j<-1
     while (j<=lkm){
         mu[1,j]<-m
         j<-j+1
     } 
  }
  data[i,]<-data[i,]+mu
  i<-i+1
}
return(data)
}





















simmix.delt<-function(n,M,sig,p,seed,dime=NULL){
#Simulates a mixture of l normal distributions in R^d,
#with diagonal cov matrices
#
#n is the sample size
#M is l*d-matrix, rows are the means
#sig is l*d-matrix, for l:th mixture d covariances
#p is l-vector, proportion for each mixture
#
#returns n*d-matrix
#
set.seed(seed) 

if (is.null(dime)){

if (dim(t(M))[1]==1) l<-1 else l<-length(M[,1])
d<-length(M[1,])
data<-matrix(rnorm(d*n),,d) #n*d matriisi, valkoista kohinaa 
for (i in 1:n){
   ehto<-runif(1)
   alku<-0
   loppu<-p[1]
   lippu<-0
   for (j in 1:(l-1)){
      if ((alku<=ehto) && (ehto<loppu)){
         data[i,]<-sig[j,]*data[i,]+M[j,]
         lippu<-1
      }
      alku<-alku+p[j]
      loppu<-loppu+p[j+1]
   }      
   if (lippu==0) data[i,]<-sig[l,]*data[i,]+M[l,]
}
}

if (!is.null(dime) && (dime==1)){
d<-1
l<-length(M)
data<-matrix(rnorm(d*n),,d) #n*d matriisi, valkoista kohinaa 
for (i in 1:n){
   ehto<-runif(1)
   alku<-0
   loppu<-p[1]
   lippu<-0
   for (j in 1:(l-1)){
      if ((alku<=ehto) && (ehto<loppu)){
         data[i]<-sig[j]*data[i]+M[j]
         lippu<-1
      }
      alku<-alku+p[j]
      loppu<-loppu+p[j+1]
   }      
   if (lippu==0) data[i]<-sig[l]*data[i]+M[l]
}
}

return(data)
}
simpp<-function(n,noisedim,D,siemen){
# tekee n*d, d=2+noisedim, datamatriisin
# 3 moodia, (0,0), (D,0), (D/2,h), h=D*sqrt(3)/2
# variance of noise dimensions is 1+D^2/6

d<-2+noisedim
hajo<-1
noisehajo<-sqrt(1+D^2/6)
h<-D*sqrt(3)/2
set.seed(siemen)
data<-matrix(rnorm(d*n),,d)   #n*d matriisi, valkoista kohinaa
data[,1:2]<-hajo*data[,1:2]
if (noisedim>0) data[,3:d]<-noisehajo*data[,3:d]
i<-1
while (i<=n){
  mu<-matrix(0,1,d)       #moodin keskipiste
  ehto<-runif(1)
  if (ehto<1/3){          #sekoitteiden painot samat
         mu[1,1]<-0
         mu[1,2]<-0
  } 
  else if (ehto>2/3){
         mu[1,1]<-D
         mu[1,2]<-0
  }
  else{
         mu[1,1]<-D/2
         mu[1,2]<-h
  }
  data[i,]<-data[i,]+mu
  i<-i+1
}
return(data)
}





















des12<-function(n,lkm,m,siemen){
#tekee n*6 data-matriisin, eli d=6
#2 moodia, m maarittaa moodien etaisyydet
#moodit tyyppia (m,m,0,...,0) ja (0,...,0) missa "lkm" kpl:tta m:ia. 
#
d<-6
set.seed(siemen)
data<-matrix(rnorm(d*n),,d)   #n*d matriisi, valkoista kohinaa
i<-1
while (i<=n){
  mu<-matrix(0,1,d)       #moodin keskipiste
  ehto<-runif(1)
  if (ehto>1/2){          #sekoitteiden painot samat
     j<-1
     while (j<=lkm){
         mu[1,j]<-m
         j<-j+1
     } 
  }
  data[i,]<-data[i,]+mu
  i<-i+1
}
return(data)
}





















simutree<-function(tr,mcn,seedi)
{

set.seed(seedi)

step<-stepcalc(tr$support,tr$N)
d<-length(tr$N)
mcdendat<-matrix(0,mcn,d)

ll<-leaflocs(tr$left,tr$right)
leafnum<-ll$leafnum
leafloc<-ll$leafloc

p<-matrix(0,leafnum,1)
for (i in 1:leafnum){
  loc<-leafloc[i]
  p[i]<-tr$volume[loc]*tr$mean[loc]
}
p<-p/sum(p)

for (i in 1:mcn){
   ehto<-runif(1)
   alku<-0
   loppu<-p[1]
   lippu<-0
   for (j in 1:(leafnum-1)){
      if ((alku<=ehto) && (ehto<loppu)){
         loc<-leafloc[j]
         ran<-runif(d)
         for (k in 1:d){
            ala<-tr$suppo[2*k-1]+step[k]*tr$low[loc,k]
            yla<-tr$suppo[2*k-1]+step[k]*tr$upp[loc,k]
            ran[k]<-ala+(yla-ala)*ran[k]   
         }
         mcdendat[i,]<-ran
         lippu<-1
      }
      alku<-alku+p[j]
      loppu<-loppu+p[j+1]
   }      
   if (lippu==0){
         loc<-leafloc[leafnum]
         ran<-runif(d)
         for (k in 1:d){
            ala<-tr$suppo[2*k-1]+step[k]*tr$low[loc,k]
            yla<-tr$suppo[2*k-1]+step[k]*tr$upp[loc,k]
            ran[k]<-ala+(yla-ala)*ran[k]   
         }
         mcdendat[i,]<-ran
   }
}

return(mcdendat)

}

slicing.recs<-function(pa,vecci,d1=1,d2=2){

lenni<-length(pa$values)
d<-length(pa$recs[1,])/2

values<-matrix(0,lenni,1)
recs<-matrix(0,lenni,4)

efek<-0

for (i in 1:lenni){

  currec<-pa$recs[i,]
  
  dimcal<-0
  onvalissa<-T
  j<-1
  while (j<=d){

     if ((j!=d1) && (j!=d2)){    
         ala<-currec[2*j-1]
         yla<-currec[2*j]
         if ((ala>vecci[j-dimcal]) || (yla<vecci[j-dimcal])) onvalissa<-F
     }
     else dimcal<-dimcal+1
     j<-j+1
  }
  if (onvalissa){
     efek<-efek+1
     values[efek]<-pa$values[i]
     recs[efek,1:2]<-pa$recs[i,(2*d1-1):(2*d1)]
     recs[efek,3:4]<-pa$recs[i,(2*d2-1):(2*d2)]

  }
}

values<-values[1:efek]
recs<-recs[1:efek,]

return(list(values=values,recs=recs))
}



stage.gaussR<-function(dendat,M,mugrid,siggrid=1,sigeka=TRUE,sampstart=TRUE)
{
n<-length(dendat)
dict.card<-length(mugrid)
dict.card.sig<-length(siggrid)

muut<-matrix(0,M,1)      #gives the means of the final mixture
sigit<-matrix(0,M,1)     #gives the std:s of the final mixture
piit<-matrix(0,M-1,1)
for (i in 1:(M-1)) piit[i]<-2/(i+2)
riskit<-matrix(0,dict.card,dict.card.sig)

if (sampstart){
   muut[1]<-mean(dendat)
   sigit[1]<-sqrt(var(dendat))
}
else{
   # haetaan 1. termi
   for (i in 1:dict.card){
      for (ii in 1:dict.card.sig){
         sqint<-gaussprod(0,0,siggrid[ii],siggrid[ii])
         val<-0
         for (j in 1:n){   
             point<-dendat[j]
             evapoint<-(point-mugrid[i])/siggrid[ii]
             val<-val+evanor(evapoint)/siggrid[ii]
         }
         riskit[i,ii]<--2*val+sqint
      }
   }
   mind<-which.min(t(riskit))
   sarat<-dict.card.sig
   imin<-ceiling(mind/sarat)
   jmin<-mind-(imin-1)*sarat

   muut[1]<-mugrid[imin]
   if (sigeka) sigit[1]<-1 else sigit[1]<-siggrid[jmin]
}

# haetaan termit 2-M
curmix<-matrix(0,M,1)   #estimaatin painot
curmix[1]<-1            #alussa vain yksi simppeli funktio
k<-1
while (k <= (M-1)){
   for (i in 1:dict.card){
      for (ii in 1:dict.card.sig){
          # calculate the -2*average of evaluations
          val<-0
          for (j in 1:n){   
              point<-dendat[j]
              evapoint<-(point-mugrid[i])/siggrid[ii]
              val<-val+evanor(evapoint)/siggrid[ii]
          }
          # calculate the inner product of the candidate with the k-1 estimate
          prodint<-0
          jj<-1
          while (jj<=k){
             prodint<-prodint+curmix[jj]*gaussprod(muut[jj],mugrid[i],
                                         sigit[jj],siggrid[ii])
             jj<-jj+1
          }
          # calculate the risk at the k:th step
          gammanpik<--2*piit[k]*val/n+piit[k]^2*gaussprod(0,0,siggrid[ii],
                                                siggrid[ii])
          riskit[i,ii]<-gammanpik+2*(1-piit[k])*piit[k]*prodint
      }
   }  
   mind<-which.min(t(riskit))
   sarat<-dict.card.sig
   imin<-ceiling(mind/sarat)
   jmin<-mind-(imin-1)*sarat

   muut[k+1]<-mugrid[imin]
   sigit[k+1]<-siggrid[jmin]

   curmix[1:k]<-(1-piit[k])*curmix[1:k]
   curmix[k+1]<-piit[k]
   k<-k+1
}

#sig<-matrix(1,M,1)
#et<-eval.func("mixt",N,sig=sigit,M=muut,p=curmix)  

return(list(muut=muut,sigit=sigit,curmix=curmix))
}  


                            

stepcalc<-function(suppo,N)
{
d<-length(N)
step<-matrix(0,d,1)
for (i in 1:d){
    step[i]<-(suppo[2*i]-suppo[2*i-1])/N[i]
}

return(step)
}




stepwiseR<-function(dendat,leafs,M,pis,mcn,minobs=0,seedi=1,
method="projec",bound=0)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

tr<-myosplitR(dendat,leafs,method,minobs)

suppo<-supp(dendat,blown=TRUE)
step<-matrix(0,d,1)
for (i in 1:d){
    step[i]<-(suppo[2*i]-suppo[2*i-1])/n
}

i<-1
while (i<=(M-1)){

   seedi<-seedi+1
   mcdendat<-simutree(tr,mcn,seedi)

   mix<-pis[i]
   trnew<-myosplitpenaR(dendat,leafs,mcdendat,mix,suppo,step,minobs)

   tr<-treeadd(tr,trnew,mix=mix)

   i<-i+1
}

return(tr)

}







supp<-function(dendat,epsi=0,blown=FALSE){
#Estimates the support of density
#
#dendat on n*xlkm matriisi
#epsi on tekn parametri
#
#Returns xlkm*2-matriisin
#
#kantajaksi estimoidaan [min-epsi,max+epsi]

n<-dim(dendat)[1]
xlkm<-length(dendat[1,])    #dendat matr sarakk lkm on muuttujien lkm
vast<-matrix(0,2*xlkm,1)  

for (i in 1:xlkm){

    minni<-min(dendat[,i])   
    maxxi<-max(dendat[,i])
    if (blown) epsi<-(maxxi-minni)/(2*(n-1))
   
    vast[2*i-1]<-minni-epsi     #sis valien alkupisteet
    vast[2*i]<-maxxi+epsi       #sis valien paatepisteet
}
return(vast)
}

treeadd<-function(tr1,tr2,cumnum=1,epsi=0,mix=NULL)
{

if (is.null(mix)) mix<-1/(cumnum+1)

d<-dim(tr1$low)[2]  #d<-length(tr1$step)     

ls<-leaflocs(tr1$left,tr1$right)

leafloc<-ls$leafloc
leafnum1<-ls$leafnum

nnum1<-length(tr1$left)
nnum2<-length(tr2$left)
maxnoden<-nnum1+leafnum1*nnum2

val<-matrix(NA,maxnoden,1)  
vec<-matrix(NA,maxnoden,1)
mean<-matrix(0,maxnoden,1)
#ssr<-matrix(0,maxnoden,1)
nelem<-matrix(0,maxnoden,1)
volume<-matrix(0,maxnoden,1)
left<-matrix(0,maxnoden,1)
right<-matrix(0,maxnoden,1)
low<-matrix(0,maxnoden,d)
upp<-matrix(0,maxnoden,d)

val[1:nnum1]<-tr1$split
vec[1:nnum1]<-tr1$direc
mean[1:nnum1]<-tr1$mean
#ssr[1:nnum1]<-tr1$ssr
nelem[1:nnum1]<-tr1$nelem
volume[1:nnum1]<-tr1$volume
left[1:nnum1]<-tr1$left
right[1:nnum1]<-tr1$right
low[1:nnum1,]<-tr1$low
upp[1:nnum1,]<-tr1$upp

curnum<-nnum1
pinotr2<-matrix(0,nnum2,1)
pinotr<-matrix(0,nnum2,1)   #upp and low will require own stack for tr
  # pinotr2 is containing pointers to tr2
  # pinotr is containing pointers to new tree tr

i<-1
while (i<=leafnum1){          # go through the leafs of tr1
 
    curleaf<-leafloc[i]

    pinoin<-1
    pinotr2[pinoin]<-1        # root
    pinotr[pinoin]<-curleaf

    while (pinoin>0){      # go through the nodes of tr2
       node<-pinotr2[pinoin]
       newleaf<-pinotr[pinoin]
       pinoin<-pinoin-1

       while (tr2$left[node]>0){   # then (!is.na(direk))

        direk<-tr2$direc[node]
        split<-tr2$split[node]

        if ((low[newleaf,direk]<split-epsi) &&
            (split+epsi<upp[newleaf,direk])){
             # make left and right children

             val[newleaf]<-split
             vec[newleaf]<-direk

             left[newleaf]<-curnum+1
             
             mean[curnum+1]<-(1-mix)*tr1$mean[curleaf]+
                             mix*tr2$mean[tr2$left[node]]
             low[curnum+1,]<-low[newleaf,]
             upp[curnum+1,]<-upp[newleaf,]
             upp[curnum+1,direk]<-split

             currec<-matrix(0,2*d,1)
             for (ii in 1:d){
                currec[2*ii-1]<-low[curnum+1,ii]
                currec[2*ii]<-upp[curnum+1,ii]
             }
             volume[curnum+1]<-massone(currec)*prod(tr1$step)

             right[newleaf]<-curnum+2

             mean[curnum+2]<-(1-mix)*tr1$mean[curleaf]+
                             mix*tr2$mean[tr2$right[node]]
             low[curnum+2,]<-low[newleaf,]
             low[curnum+2,direk]<-split
             upp[curnum+2,]<-upp[newleaf,]

             for (ii in 1:d){
                currec[2*ii-1]<-low[curnum+2,ii]
                currec[2*ii]<-upp[curnum+2,ii]
             }
             volume[curnum+2]<-massone(currec)*prod(tr1$step)

             # put right node to the stack (if exists)

             pinoin<-pinoin+1
             pinotr2[pinoin]<-tr2$right[node]
             pinotr[pinoin]<-curnum+2

             # go to left

             node<-tr2$left[node]
             newleaf<-curnum+1
                           
             # update curnum

             curnum<-curnum+2

        }
        else{
           if (split<=tr1$low[curleaf,direk]){
              #do not make children, go to right in tr2
              if (tr2$right[node]>0){
                 node<-tr2$right[node]
              }
           }
           else{ #(split>=tr1$upp[curleaf,direk])
              #do not make children, go to left in tr2
              if (tr2$left[node]>0){
                  node<-tr2$left[node]
              }
           }
       }

    } # while node>0

    } # loop for tr2
    
    i<-i+1
} # loop for tr1

val<-val[1:curnum]
vec<-vec[1:curnum]
mean<-mean[1:curnum]
nelem<-nelem[1:curnum]
volume<-volume[1:curnum]
left<-left[1:curnum]
right<-right[1:curnum]
low<-low[1:curnum,]
upp<-upp[1:curnum,]

tr<-list(split=val,direc=vec,mean=mean,nelem=nelem,volume=volume,
left=left,right=right,low=low,upp=upp,
N=tr1$N,support=tr1$support)   #step=tr1$step)
return(tr)
}














unidis<-function(ulot){
#tasainen jakauma vec:n elementeille
#ulot=dim(vec)
#tasainen jakauma joukossa {1,2,...,ulot}
#
arpa<-runif(1)
ele<-ceiling(ulot*arpa)
return(ele)
}
