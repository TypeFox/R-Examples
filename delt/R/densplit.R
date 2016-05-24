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








