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















