addnode<-function(inde,curre,curdep,left,right,parent,low,upp,N,numnode){
#(inde,curre,curdep,left,right,deplink,low,upp,enofatdep,N,numnode){
#
#inde is d-vector: index (gridpoint) to be added
#curre is pointer to vectors left,right,...
#
d<-length(inde)
apu<-depth2com(curdep,N)
curdir<-apu$direc
depatd<-apu$depind
depit<-log(N,base=2)
#depit[d]<-depit[d]+1
#
while (curdir<=(d-1)){
    ind<-inde[curdir]
    while (depatd<=depit[curdir]){
        mid<-(low[curre]+upp[curre])/2
        if (ind<=mid){
           left[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-low[curre]
           upp[numnode+1]<-floor(mid)
        }
        else{
           right[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-ceiling(mid)
           upp[numnode+1]<-upp[curre]
        }
        numnode<-numnode+1
        curre<-numnode
        depatd<-depatd+1
        curdep<-curdep+1
#        deplink[endofatdep[curdep]]<-numnode
#        deplink[numnode]<-0
#        endofatdep[curdep]<-numnode
    }
    #
    # Last node of this dimension (first node of next dimension)
    #
    curdir<-curdir+1
    ind<-inde[curdir]
    low[curre]<-1
    upp[curre]<-N[curdir]
    mid<-(low[curre]+upp[curre])/2
    if (ind<=mid){
           left[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-low[curre]
           upp[numnode+1]<-floor(mid)
    }
    else{
           right[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-ceiling(mid)
           upp[numnode+1]<-upp[curre]
    }
    depatd<-2
    numnode<-numnode+1
    curre<-numnode
    curdep<-curdep+1
#    deplink[endofatdep[curdep]]<-numnode
#    deplink[curre]<-0
#    endofatdep[curdep]<-numnode
}
#
# Last dimension 
#
ind<-inde[curdir]
while (depatd<=depit[curdir]){
        mid<-(low[curre]+upp[curre])/2
        if (ind<=mid){
           left[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-low[curre]
           upp[numnode+1]<-floor(mid)
        }
        else{
           right[curre]<-numnode+1
           parent[numnode+1]<-curre
           low[numnode+1]<-ceiling(mid)
           upp[numnode+1]<-upp[curre]
        }
        numnode<-numnode+1
        curre<-numnode
        depatd<-depatd+1
        curdep<-curdep+1
#        deplink[endofatdep[curdep]]<-numnode
#        deplink[curre]<-0
#        endofatdep[curdep]<-numnode
}
#
# Last node of last dimension
#
#left[curre]<-0
#right[curre]<-0
#
#return(list(numnode=numnode,left=left,right=right,deplink=deplink,low=low,
#upp=upp,endofatdep=endofatdep))
return(list(numnode=numnode,left=left,right=right,parent=parent,low=low,
upp=upp,nodeloc=numnode))
}










allokoi.new<-function(cur,vecs,lst,left,right,sibord)
{
# allocates space for all children of "cur"

# Calculate the number of childs and sum of volumes of childs
now<-left[cur]
childnum<-1
childvolume<-lst$volume[now]
while (right[now]>0){
  now<-right[now]
  childnum<-childnum+1
  childvolume<-childvolume+lst$volume[now]
}
 
gaplen<-(lst$volume[cur]-childvolume)/(childnum+1)

if (childnum==1){
   now<-left[cur]
   xbeg<-gaplen+vecs[cur,1]
   xend<-xbeg+lst$volume[now]
   ycoo<-lst$level[now]
   vecs[now,]<-c(xbeg,xend,ycoo)
}
else{
  siblinks<-matrix(0,childnum,1)  #make siblinks in right order
  now<-left[cur] 
  sior<-sibord[now]
  siblinks[sior]<-now
  while (right[now]>0){
    now<-right[now]
    sior<-sibord[now]
    siblinks[sior]<-now
  }
  xend<-vecs[cur,1]      #initialize xend 
  for (i in 1:childnum){
     now<-siblinks[i]
     xbeg<-gaplen+xend
     xend<-xbeg+lst$volume[now]
     ycoo<-lst$level[now]
     vecs[now,]<-c(xbeg,xend,ycoo)
  }
} 


return(vecs)
}


allokoi<-function(vecs,cur,child,sibling,sibord,levels,volumes)
{
#Finds coordinates of a node
#sibord,levels,volumes are nodenum-vector

# Calculate the number of childs and sum of volumes of childs
now<-child[cur]
childnum<-1
childvolume<-volumes[now]
while (sibling[now]>0){
  now<-sibling[now]
  childnum<-childnum+1
  childvolume<-childvolume+volumes[now]
}
 
gaplen<-(volumes[cur]-childvolume)/(childnum+1)

if (childnum==1){
   now<-child[cur]
   xbeg<-gaplen+vecs[cur,1]
   xend<-xbeg+volumes[now]
   ycoo<-levels[now]
   vecs[now,]<-c(xbeg,ycoo,xend,ycoo)
}
else{
  siblinks<-matrix(0,childnum,1)  #make siblinks in right order
  now<-child[cur] 
  sior<-sibord[now]
  siblinks[sior]<-now
  while (sibling[now]>0){
    now<-sibling[now]
    sior<-sibord[now]
    siblinks[sior]<-now
  }
  xend<-vecs[cur,1]      #initialize xend 
  for (i in 1:childnum){
     now<-siblinks[i]
     xbeg<-gaplen+xend
     xend<-xbeg+volumes[now]
     ycoo<-levels[now]
     vecs[now,]<-c(xbeg,ycoo,xend,ycoo)
  }
} 
return(vecs)
}



alloroot<-function(vecs,roots,sibord,levels,volumes)
{
rootnum<-length(roots)

# Calculate sum of volumes of roots
rootsvolume<-0
for (i in 1:rootnum){
  now<-roots[i]
  rootsvolume<-rootsvolume+volumes[now]
}

basis<-rootsvolume+rootsvolume/4
 
gaplen<-(basis-rootsvolume)/(rootnum+1)

rootlinks<-matrix(0,rootnum,1)  #make links in right order

if (rootnum==1) rootlinks[1]<-roots[1]  #1
else{
for (i in 1:rootnum){
  now<-roots[i]
  roor<-sibord[now]
  rootlinks[roor]<-now
}
}
xbeg<-0
xend<-0
for (i in 1:rootnum){
  now<-rootlinks[i]
  xbeg<-gaplen+xend
  xend<-xbeg+volumes[now]
  ycoo<-levels[now]
  vecs[now,]<-c(xbeg,ycoo,xend,ycoo)
}
return(vecs)
}
alpha.complex<-function(complex,dendat,alpha)
{
M<-dim(complex)[1]
n<-dim(dendat)[1]
d<-dim(dendat)[2]  # d<-dim(complex)[2]-1  

acomplex<-matrix(0,M,d+1)
lkm<-0
for (m in 1:M){
    simindex<-complex[m,]
    simplex<-dendat[simindex,]

    tulos<-0
    i<-1
    while ((i<=d) && (tulos==0)){
       v1<-simplex[i,]
       j<-i+1
       while ((j<=(d+1)) && (tulos==0)){
         v2<-simplex[j,]
         etais2<-sum((v1-v2)^2)
         if (etais2>alpha^2) tulos<-1
         j<-j+1
       }
       i<-i+1
    }
    if (tulos==0){ 
       lkm<-lkm+1
       acomplex[lkm,]<-complex[m,]
    }
}
acomplex<-acomplex[1:lkm,]

return(acomplex)
}

blokitus2<-function(obj,blokki){
#
sar<-length(obj[1,]) #sarakkeiden maara 
riv<-length(obj[,1]) #rivien maara 
#
uusobj<-matrix(0,riv,sar+blokki)
uusobj[,1:sar]<-obj
#
return(uusobj)
}
blokitus<-function(obj,blokki){
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
boundbox<-function(rec1,rec2)
{
# rec:s are 2*d-vectors

d<-length(rec1)/2
rec<-matrix(0,2*d,1)

for (i in 1:d){
    rec[2*i-1]<-min(rec1[2*i-1],rec2[2*i-1])
    rec[2*i]<-max(rec1[2*i],rec2[2*i])
}

return(rec)
}

branchmap<-function(estiseq,hseq=NULL,levnum=80,paletti=NULL,rootpaletti=NULL,
type="jump")
{
#type= "smooth", "jump", "diffe" 

if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:100])
if (is.null(rootpaletti)) rootpaletti<-colors()[102:110]

lstseq<-estiseq$lstseq
if (is.null(hseq))
   if (!is.null(estiseq$type)){
       if (estiseq$type=="bagghisto") hseq<--estiseq$hseq
       if (estiseq$type=="carthisto")  hseq<--estiseq$leaf
       if (estiseq$type=="kernel")  hseq<-estiseq$hseq
   }
   else hseq<-estiseq$hseq
hnum<-length(hseq)

if (hseq[1]>hseq[2]){  
    hseq<-hseq[seq(hnum,1)]
    apuseq<-list(lstseq[[hnum]])
    i<-2
    while (i <= hnum){
         apuseq<-c(apuseq,list(lstseq[[hnum-i+1]]))
         i<-i+1 
   }
   lstseq<-apuseq
}

maxlevel<-0
i<-1
while (i<=hnum){
   lst<-lstseq[[i]]
   maxlevel<-max(max(lst$level),maxlevel)
   i<-i+1
}

levstep<-maxlevel/(levnum-1)
level<-seq(0,maxlevel,levstep)
z<-matrix(0,length(level)+1,length(hseq)+1)
#col<-matrix("white",length(level),length(hseq))
colot<-matrix("white",(length(level))*(length(hseq)),1)

i<-1
while (i<=hnum){
    lst<-lstseq[[i]]  #[[hnum-i+1]]

    if ((type=="smooth") || (type=="diffe")){
          eb<-excmas.bylevel(lst,length(level)+1)
          if (type=="smooth")  z[,i]<-eb$levexc
          else z[,i]<-eb$diffe
    }
    mut<-multitree(lst$parent)
    ex<-excmas(lst)
    fb<-findbranch.pare(lst$parent)
    if (is.null(fb)) branchnum<-0 else branchnum<-length(fb)

    if (branchnum==0) toplevel<-max(lst$level) else toplevel<-min(lst$level[fb])
    # toplevel is the level of the next branch
    rootnum<-length(mut$roots)
    rootstep<-toplevel/rootnum
    ordroots<-order(ex[mut$roots])  #order(lst$level[mut$roots])

    exmassa<-0     
    k<-1
    while (k<=rootnum){
        ind<-mut$roots[ordroots[k]]
        exmassa<-exmassa+ex[ind]  
        k<-k+1
    }

    leveka<-1
    levend<-max(leveka,min(round(levnum*toplevel/maxlevel),levnum))
    curleveka<-leveka
    k<-1
    while (k<=rootnum){
        ind<-mut$roots[ordroots[k]]
        curexma<-ex[ind] 
curlevend<-max(curleveka,min(round(curleveka+(levend-leveka)*curexma/exmassa),levend))
        if (type=="jump") z[curleveka:curlevend,i]<-exmassa
        ##col[curleveka:curlevend,i]<-rootpaletti[k]
        aa<-(i-1)*(levnum)+min(curleveka,levnum)
        bb<-(i-1)*(levnum)+min(curlevend,levnum)
        colot[aa:bb]<-rootpaletti[k]
        curleveka<-curlevend+1
        k<-k+1
    }

    curlevel<-toplevel   # curlevel is the level of the previous branching
    ordbranches<-order(lst$level[fb])
    k<-1
    while (k<=branchnum){

        branchind<-ordbranches[k]
        branch<-fb[branchind]

        if (k==branchnum) toplevel<-max(lst$level) 
        else{
            nextbranch<-fb[ordbranches[k+1]]
            toplevel<-lst$level[nextbranch]
        }
        childnum<-2
        children<-c(mut$child[branch],mut$sibling[mut$child[branch]])
        ordchild<-order(ex[children])  #order(lst$level[children])

        exmassa<-0     
        l<-1
        while (l<=childnum){
             ind<-children[ordchild[l]]
             exmassa<-exmassa+ex[ind]  
             l<-l+1
        }

        leveka<-curlevend+1
        levend<-max(leveka,min(leveka+round(levnum*(toplevel-curlevel)/maxlevel),levnum))
        curleveka<-leveka
        l<-1
        while (l<=childnum){
           ind<-children[ordchild[l]]
           curexma<-ex[ind]
  curlevend<-max(curleveka,min(round(curleveka+(levend-leveka)*curexma/exmassa),levend)) 
           if (type=="jump") z[curleveka:curlevend,i]<-exmassa
           ##col[curleveka:curlevend,i]<-paletti[l]
           aa<-(i-1)*(levnum)+min(curleveka,levnum)
           bb<-(i-1)*(levnum)+min(curlevend,levnum)
           colot[aa:bb]<-paletti[l]
           curleveka<-curlevend+1
           l<-l+1
        }
        curlevel<-toplevel
        k<-k+1
    }
    i<-i+1
}

z[,dim(z)[2]]<-z[,dim(z)[2]-1]
z[dim(z)[1],]<-z[dim(z)[1]-1,]
hseq[length(hseq)+1]<-hseq[length(hseq)]+hseq[length(hseq)]-hseq[length(hseq)-1]
level[length(level)+1]<-level[length(level)]+level[length(level)]-level[length(level)-1]

z<-z/max(z)

# add one column to the matrix: a new first column

lisa<-1
zapu<-matrix(0,dim(z)[1],dim(z)[2]+lisa)
zapu[,1:lisa]<-0
zapu[,(lisa+1):(lisa+dim(z)[2])]<-z

yapu<-matrix(0,length(hseq)+lisa,1)
ystep<-hseq[2]-hseq[1]
yapu[lisa:1]<-seq(hseq[1]-ystep,hseq[1]-ystep*lisa,-ystep)
yapu[(lisa+1):(length(hseq)+lisa)]<-hseq

# add colors to the end 
levelo<-lisa*(dim(zapu)[1]-1)
colapu<-matrix("",length(colot)+levelo,1)
colapu[1:length(colot)]<-colot
colapu[(length(colot)+1):length(colapu)]<-colot[(length(colot)-levelo+1):length(colot)]

return(list(level=level,h=yapu,z=zapu,col=colapu))
}





ccentebag<-function(component,AtomlistAtom,AtomlistNext,low,upp,volume,
step,suppo)
{
d<-dim(low)[2]

componum<-length(component)
center<-matrix(0,componum,d)

for (i in 1:componum){
   curcente<-matrix(0,d,1)
   pointer<-component[i]
   while (pointer>0){
        atompointer<-AtomlistAtom[pointer]
        
        newcente<-matrix(0,d,1)
        for (j in 1:d){
            # calculate 1st volume of d-1 dimensional rectangle where
            # we have removed j:th dimension

            vol<-1
            k<-1
            while (k<=d){
               if (k!=j){
                  vol<-vol*(upp[atompointer,k]-low[atompointer,k])*step[k]
               }
               k<-k+1
            }

            ala<-suppo[2*j-1]+step[j]*low[atompointer,j]
            yla<-suppo[2*j-1]+step[j]*upp[atompointer,j]
            newcente[j]<-vol*(yla^2-ala^2)/2
        }

        curcente<-curcente+newcente
        pointer<-AtomlistNext[pointer]
   }
   center[i,]<-curcente/volume[i]
}
return(t(center))
}

ccentedya<-function(volofatom,component,AtomlistNext,AtomlistAtom,
volume,minim,h,delta,index,d){
#
componum<-length(component)
center<-matrix(0,componum,d)
#
for (i in 1:componum){
   curcente<-0
   pointer<-component[i]
   while (pointer>0){
        atompointer<-AtomlistAtom[pointer]
        inde<-index[atompointer,]
        newcente<-minim-h+delta*inde
        curcente<-curcente+newcente
        pointer<-AtomlistNext[pointer]
   }
   center[i,]<-volofatom*curcente/volume[i]
}
return(t(center))
}
ccente<-function(levels,items,mass){
#Calculates centers from a collection of level sets.
#center is 1st moment didided by volume.
#
#levels is tasolkm*N-matrix of 1:s and 0:s
#items is N*(2*d)-matrix
#mass is tasolkm-vector
#
#returns N*d-matrix of 1st moments.
#
N<-length(levels[,1])
d<-length(items[1,])/2
res<-matrix(0,N,d)
if (dim(t(levels))[1]==1) tasolkm<-1 else tasolkm<-length(levels[,1]) 
for (i in 1:tasolkm){
  lev2<-change(levels[i,])
  m<-length(lev2)
  vol<-matrix(0,d,1)
  for (j in 1:m){
    ind<-lev2[j]
    rec<-items[ind,]
    vol<-vol+cenone(rec)
  }
  res[i,]<-vol/mass[i]
}
return(t(res))
}

 

cenone<-function(rec){
#Calculates the 1st moment of a rectangle.
#
#rec is (2*d)-vector, represents rectangle in d-space
#Returns a d-vector.
#
d<-length(rec)/2
res<-matrix(0,d,1)
for (j in 1:d){
  apurec<-rec      #apurec such that is volume is equal to
  apurec[2*j-1]<-0 #volume of d-1 dimensional rectangle where
  apurec[2*j]<-1   #we have removed j:th dimension
  vajmas<-massone(apurec) 
  res[j]<-vajmas*(rec[2*j]^2-rec[2*j-1]^2)/2  
}
return(res) 
}

cfrekv<-function(levels,arvot){
#laskee tasojoukon osien frekvenssit
#arvo on reaaliluku
#kumu on kork*n-matriisi, n saraketta, kuvaa kork kpl:tta tasojoukon osia
#1 jos vastaava data-matriisin rivin indikoima pallo kuuluu tasojouon osaan
#muodostetaan matriisi, jonka 1. sarakkeessa "arvo", 
#2. sarakkeessa kunkin tasojoukon osan frekvenssi  
#ts. laskettu kuinka monesta pallosta tasojoukko on yhdistetty
#
tasolkm<-length(levels[,1])     #levels:n rivien maara
frek<-matrix(0,tasolkm,1)
a<-1
while (a<=tasolkm){
   frek[a]<-sum(levels[a,]*arvot)
   a<-a+1 
}
return(t(frek))
}







change<-function(levset){
#
#
len<-length(levset)
m<-sum(levset)
rindeksit<-matrix(0,m,1)
j<-1
for (i in 1:len){
    if (levset[i]==1){
       rindeksit[j]<-i
       j<-j+1
    }
}
return(rindeksit)
}

cintemul<-function(roots,child,sibling,volume,level){
#
#integrate function, over the level of roots, in the region of roots
#
itemnum<-length(child)
rootnum<-length(roots)
inte<-0
for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    valpino<-matrix(0,itemnum,1)  #level of parent
    pino[1]<-roots[i]
    valpino[1]<-0
    sibling[roots[i]]<-0
    #    
    pinin<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        valcur<-valpino[pinin] 
        pinin<-pinin-1
        #
        if (level[cur]>0){
           inte<-inte+(level[cur]-max(valcur,0))*volume[cur]
        }
        #
        if (sibling[cur]>0){
              pinin<-pinin+1
              pino[pinin]<-sibling[cur]
              valpino[pinin]<-valcur
        }
        while (child[cur]>0){    #go to left and put right nodes to stack
              valcur<-level[cur]
              cur<-child[cur]
              #
              if (level[cur]>0){
                 inte<-inte+(level[cur]-max(valcur,0))*volume[cur]
              }
              #
              if (sibling[cur]>0){  #if cur has siblings
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
                 valpino[pinin]<-valcur
             }
        }
    }
}
#
return(inte)
}

cinte<-function(values,volumes,parents){
#Calculates the integral of a piecewise continuous function.
#
len<-length(values)
int<-0
for (i in len:1){
  par<-parents[i]
  if (par==0) valpar<-0 
  else valpar<-values[par]
  int<-int+volumes[i]*(values[i]-valpar)
}
return(int)
}

colo2scem<-function(sp,mt,ca)
{
#sp result of scemprof
#mt result of modegraph
#ca result of coloallo
#origlis translates h-values from sp terminology to mt terminology

len<-length(sp$bigdepths)
mtlen<-length(ca)
colors<-matrix("black",len,1)

for (i in 1:len){
  label<-sp$mlabel[i]       #label for mode
  if (label>0){ 
      smoot<-sp$smoot[i]    #smoothing paramter value/leafnum
      # we find the corresponding slot from "ca" where
      # label corresponds and smoothing parameter value corresponds
      run<-1
      koesmoot<-mt$ycoor[run]
      koelabel<-mt$mlabel[run] 
      while (((koesmoot!=smoot) || (koelabel!=label)) && (run<=mtlen)){
         run<-run+1
         koesmoot<-mt$ycoor[run]
         koelabel<-mt$mlabel[run] 
      }
      # we have found the slot
      colors[i]<-ca[run]
  }
}

return(colors)
}








coloallo<-function(mt,paletti=NULL)
{
# fast allocation of colors (matching of modes)
# mt is mode tree
# paletti gives a list of colors

if (is.null(paletti))
paletti<-c("red","blue","green","turquoise","orange","navy",
"darkgreen","orchid",colors()[50:100])

d<-dim(mt$xcoor)[2]

snum<-0
for (i in 1:length(mt$mlabel)){
  if (mt$mlabel[i]==1) snum<-snum+1
}

xcoor<-mt$xcoor
ycoor<-mt$ycoor
mlabel<-mt$mlabel
lenni<-length(ycoor)

colot<-matrix("",lenni,1)

# find the locations for the information for each h

low<-matrix(0,snum,1)
upp<-matrix(0,snum,1)
low[1]<-1
glob<-2
while ((glob<=lenni) && (mlabel[glob]!=1)){
       glob<-glob+1
}
upp[1]<-glob-1
# now glob is at the start of new block
i<-2
while (i<=snum){
   low[i]<-glob
   glob<-glob+1
   while ((glob<=lenni) && (mlabel[glob]!=1)){
       glob<-glob+1
   }
   upp[i]<-glob-1
   i<-i+1
}

# first we allocate colors for the largest h

run<-1  #low[1]
while (run<=upp[1]){
   colot[run]<-paletti[run]
   run<-run+1
}

firstnewcolo<-run

i<-2

while (i<=snum){
   prenum<-upp[i-1]-low[i-1]+1
   curnum<-upp[i]-low[i]+1

   smallernum<-min(prenum,curnum)
   greaternum<-max(prenum,curnum)

   if (prenum==smallernum){   
        bases<-i
        compa<-i-1
   }
   else{
        bases<-i-1
        compa<-i
   }

      dista<-matrix(NA,smallernum,greaternum)
      for (ap in low[bases]:upp[bases]){
         for (be in low[compa]:upp[compa]){
           if (d==1){
               curcenter<-xcoor[ap]
               precenter<-xcoor[be]
           }
           else{
               curcenter<-xcoor[ap,]
               precenter<-xcoor[be,]
           }
           dista[be-low[compa]+1,ap-low[bases]+1]<-etais(curcenter,precenter)
         }
      }

      match<-matrix(0,smallernum,1)  #for each mode the best match
      findtie<-TRUE

      # find the best match for all and check whether there are ties
      match<-matrix(0,smallernum,1)
      for (bm in 1:smallernum){
          minimi<-min(dista[bm,],na.rm=TRUE)
          match[bm]<-which(minimi==dista[bm,])[1]
      }
      findtie<-FALSE
      bm<-1
      while ((bm<=smallernum) && (findtie==FALSE)){
         koe<-match[bm]
         bm2<-bm+1
         while (bm2<=smallernum){
            if (koe==match[bm2]){
                  findtie<-TRUE
            }
            bm2<-bm2+1
         }
         bm<-bm+1
      }
    
      onkayty<-FALSE

      while (findtie){

      onkayty<-TRUE
      tiematch<-matrix(0,smallernum,1)
      
      # find the best match for all
      bestmatch<-matrix(0,smallernum,1)
      for (bm in 1:smallernum){
          allna<-TRUE
          am<-1
          while ((am<=greaternum) && (allna)){
             if (!is.na(dista[bm,am])) allna<-FALSE
             am<-am+1
          }
          if (!(allna)){
             minimi<-min(dista[bm,],na.rm=TRUE)
             bestmatch[bm]<-which(minimi==dista[bm,])[1]
          }
          else bestmatch[bm]<-match[bm]
      }

      # find the first tie
      findtie<-FALSE

      tieset<-matrix(0,smallernum,1)
      bm<-1
      while ((bm<=smallernum) && (findtie==FALSE)){
         koe<-bestmatch[bm]
         bm2<-bm+1
         while (bm2<=smallernum){
            if (koe==bestmatch[bm2]){
                  findtie<-TRUE
                  tieset[bm]<-1
                  tieset[bm2]<-1
            }
            bm2<-bm2+1
         }
         bm<-bm+1
      }

      # solve the first tie
      if (findtie==TRUE){
         numofties<-sum(tieset)
         kavelija<-0
         tiepointer<-matrix(0,numofties,1) 
         # find the second best
         secondbest<-matrix(0,smallernum,1)
         for (bm in 1:smallernum){
            if (tieset[bm]==1){
               redudista<-dista[bm,]
               redudista[bestmatch[bm]]<-NA
               minimi<-min(redudista,na.rm=TRUE)
               secondbest[bm]<-which(minimi==redudista)[1]

               kavelija<-kavelija+1
               tiepointer[kavelija]<-bm
            }
         }
         # try different combinations       
         # try all subsets of size 2 from the set of ties
         numofsubsets<-choose(numofties,2)
            #gamma(numofties+1)/gamma(numofties-2+1)
         valuelist<-matrix(0,numofsubsets,1)
         vinnerlist<-matrix(0,numofsubsets,1)
         matchlist<-matrix(0,numofsubsets,1)
         runneri<-1
         eka<-1
         while (eka<=numofties){
            ekapo<-tiepointer[eka]
            toka<-eka+1
            while (toka<=numofties){
               tokapo<-tiepointer[toka]
               # try combinations for this subset (there are 2)
               # 1st combination
               fvinner<-ekapo
               fvinnermatch<-bestmatch[fvinner]
               floser<-tokapo
               flosermatch<-secondbest[floser]
               fvalue<-dista[fvinner,fvinnermatch]+dista[floser,flosermatch]
                # 2nd combination
               svinner<-tokapo
               svinnermatch<-bestmatch[svinner]
               sloser<-ekapo
               slosermatch<-secondbest[sloser]
               svalue<-dista[svinner,svinnermatch]+dista[sloser,slosermatch]
               # tournament
               if (fvalue<svalue){
                   valuelist[runneri]<-fvalue
                   vinnerlist[runneri]<-fvinner
                   matchlist[runneri]<-fvinnermatch
               }
               else{ 
                   valuelist[runneri]<-svalue
                   vinnerlist[runneri]<-svinner
                   matchlist[runneri]<-svinnermatch
               }
               runneri<-runneri+1 
               # 
               toka<-toka+1
            }
            eka<-eka+1
         }
         minimi<-min(valuelist,na.rm=TRUE)
         bestsub<-which(minimi==valuelist)[1]
         vinnerson<-vinnerlist[bestsub]
         matcherson<-matchlist[bestsub]

         tiematch[vinnerson]<-matcherson
         dista[vinnerson,]<-NA
         dista[,matcherson]<-NA

      }

      }  #while (findtie)

      if (onkayty){  #there was one tie
          
          for (sepo in 1:smallernum){
               if (tiematch[sepo]!=0) match[sepo]<-tiematch[sepo]
               else match[sepo]<-bestmatch[sepo]
          }
      }

      # finally allocate colors
      run<-1
      while (run<=smallernum){
          
          if (prenum==smallernum){
             xind<-run
             yind<-match[xind]
          }
          else{
             yind<-run
             xind<-match[yind]
          }

          colot[low[i]+yind-1]<-colot[low[i-1]+xind-1]    
          run<-run+1
      }
                    
      if (prenum<greaternum){

        run<-low[bases]
        while (run<=upp[bases]){
            if (colot[run]==""){
               colot[run]<-paletti[firstnewcolo]
               firstnewcolo<-firstnewcolo+1   
            }
            run<-run+1
        }

     }

     i<-i+1
}

return(colot)
}


















colobary.merge<-function(parent,level,colothre=min(level),paletti=NULL)
{
mt<-multitree(parent) #roots<-mt$roots child<-mt$child sibling<-mt$sibling

itemnum<-length(mt$child)
rootnum<-length(mt$roots)
if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])
hep<-1

colot<-matrix("",itemnum,1)
col<-colobary(parent,paletti)

for (i in 1:rootnum){
   curroot<-mt$roots[i]
   colot[curroot]<-col[curroot]  #"grey"  #paletti[hep]
   hep<-hep+1
   if (mt$child[curroot]>0){
      pino<-matrix(0,itemnum,1)
      pino[1]<-mt$child[curroot]
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          #if (level[mt$child[cur]]>colothre)
          #if (level[parent[cur]]>colothre)
          if (level[cur]>colothre)
              colot[cur]<-colot[parent[cur]]
          else{ 
                colot[cur]<-col[cur]  #"grey"  #paletti[hep]
                hep<-hep+1
          } 
          # put to the stack 
          if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
          }
          # go to left and put right nodes to the stack
          while (mt$child[cur]>0){   
             cur<-mt$child[cur]
             #if (level[mt$child[cur]]>colothre)
             #if (level[parent[cur]]>colothre)
             if (level[cur]>colothre) 
                       colot[cur]<-colot[parent[cur]]
             else{ 
                    colot[cur]<-col[cur]  #"grey"  #paletti[hep]
                    hep<-hep+1
             } 
             if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
             }
           }
       }#while (pinin>0)
   }
}       

ind<-(level<colothre)
colot[ind]<-"grey"
                    
return(colot)
}


colobary.nodes<-function(parent,nodes,paletti=NULL)
{
mt<-multitree(parent) #roots<-mt$roots child<-mt$child sibling<-mt$sibling

itemnum<-length(mt$child)
rootnum<-length(mt$roots)
if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])

colot<-matrix("",itemnum,1)
col<-colobary(parent,paletti)

nodenum<-length(parent)
allnodes<-matrix(0,nodenum,1)
#allnodes[1:length(nodes)]<-nodes
counter<-0  #length(nodes)
for (i in 1:length(nodes)){
    node<-nodes[i]
    tt<-travel.tree(parent,node)
    allnodes[(counter+1):(counter+length(tt))]<-tt
    counter<-counter+length(tt)
}
allnodes<-allnodes[1:counter]

for (i in 1:rootnum){
   curroot<-mt$roots[i]
   colot[curroot]<-col[curroot]  #"grey"  #paletti[hep]
   if (mt$child[curroot]>0){
      pino<-matrix(0,itemnum,1)
      pino[1]<-mt$child[curroot]
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          if (sum(cur==allnodes)>0)
              colot[cur]<-colot[parent[cur]]
          else{ 
                colot[cur]<-col[cur]  #"grey"  #paletti[hep]
          } 
          # put to the stack 
          if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
          }
          # go to left and put right nodes to the stack
          while (mt$child[cur]>0){   
             cur<-mt$child[cur]
             if (sum(cur==allnodes)>0)
                       colot[cur]<-colot[parent[cur]]
             else{ 
                    colot[cur]<-col[cur]  #"grey"  #paletti[hep]
             } 
             if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
             }
           }
       }#while (pinin>0)
   }
}       

ind<-setdiff(seq(1:itemnum),allnodes)
colot[ind]<-"grey"
                    
return(colot)
}


colobary<-function(parent,paletti,roots=NULL,
modecolo=NULL,modepointer=NULL #,segtype="char"
)
{
nodenum<-length(parent)
#if (segtype=="char") colot<-matrix("",nodenum,1) 
#else 
colot<-matrix(0,nodenum,1)

fb<-findbranch(parent)$indicator
modloc<-moodilkm(parent)$modloc
#if (repretype=="B"){
#   fb<-findbranchB(parent,roots)$indicator
#   modloc<-moodilkmB(parent)$modloc
#}

moodilkm<-length(modloc)
palerun<-0

# first allocate colors for modes
if (is.null(modecolo)){
   i<-1
   while (i<=moodilkm){
       cur<-modloc[i]
       palerun<-palerun+1
       colot[cur]<-paletti[palerun]
       i<-i+1
   }
}
else{
   # remove modecolo:s from paletti
   indu<-0
   for (pp in 1:length(paletti)) 
       for (ppp in 1:length(modecolo))
          if (paletti[pp]==modecolo[ppp]){ 
                 indu<-indu+1
                 paletti[pp]<-colors()[100+indu] 
          }
   
   i<-1
   while (i<=moodilkm){
       cur<-modepointer[i]
       colot[cur]<-modecolo[i]
       i<-i+1
   } 
}

# then allocate for others
i<-1
while (i<=moodilkm){
 
  cur<-modloc[i]
  while (parent[cur]>0){

     child<-parent[cur]

     if ((fb[cur]==1) && (colot[child]==0)){ #cur is a result of a branch 
           palerun<-palerun+1
           colot[child]<-paletti[palerun]
     }      
     else if (colot[child]==0) colot[child]<-colot[cur]

     cur<-child
  }
  i<-i+1
}

return(colot)
}    


colobary.roots<-function(parent,level,colothre=min(level),paletti=NULL)
{
mt<-multitree(parent) #roots<-mt$roots child<-mt$child sibling<-mt$sibling

itemnum<-length(mt$child)
colot<-matrix("",itemnum,1)
rootnum<-length(mt$roots)
if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])
hep<-1

for (i in 1:rootnum){
   curroot<-mt$roots[i]
   colot[curroot]<-paletti[hep]
   hep<-hep+1
   if (mt$child[curroot]>0){
      pino<-matrix(0,itemnum,1)
      pino[1]<-mt$child[curroot]
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          if ((mt$sibling[cur]==0)||(level[parent[cur]]<colothre)) 
              colot[cur]<-colot[parent[cur]]
          else{ 
                colot[cur]<-paletti[hep]
                hep<-hep+1
          } 
          # put to the stack 
          if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
          }
          # go to left and put right nodes to the stack
          while (mt$child[cur]>0){   
             cur<-mt$child[cur]
             if ((mt$sibling[cur]==0)||(level[parent[cur]]<colothre)) 
                       colot[cur]<-colot[parent[cur]]
             else{ 
                    colot[cur]<-paletti[hep]
                    hep<-hep+1
             } 
             if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
             }
           }
       }#while (pinin>0)
   }
}                           
return(colot)
}




 colors2data<-function(dendat,pcf,lst,paletti=NULL,clusterlevel=NULL,nodes=NULL,
type="regular")
{
if (type=="regular")
return( colorsofdata(dendat,pcf,lst,paletti=paletti,clusterlevel=clusterlevel,nodes=nodes) )
else
return( colorsofdata.adagrid(dendat,pcf,lst,paletti=paletti,clusterlevel=clusterlevel,nodes=nodes) )

}

colorsofdata2<-function(dendat,pcf,lst,paletti=NULL,
clusterlevel=NULL,nodes=NULL)
{
# links from dendat to rec to node to color
# "lst$infopointer" gives links from nodes to recs

n<-dim(dendat)[1]
d<-dim(dendat)[2]
rnum<-length(pcf$value)

step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
    
if (is.null(paletti))
paletti<-c("red","blue","green",
    "orange","navy","darkgreen",
    "orchid","aquamarine","turquoise",
    "pink","violet","magenta","chocolate","cyan",
    colors()[50:657],colors()[50:657])

# links from node to color
if ((is.null(clusterlevel))&&(is.null(nodes))) col<-colobary(lst$parent,paletti)
if (!is.null(clusterlevel))
col<-colobary.merge(lst$parent,lst$level,colothre=clusterlevel,paletti)
if (!is.null(nodes))
col<-colobary.nodes(lst$parent,nodes,paletti)

# links from rec to node (invert the links in infopointer)
nodefinder<-matrix(0,rnum,1)
for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i

# find links from dendat to rec
den2pcf<-matrix(0,n,1)
pcf2den<-matrix(0,rnum,1)
value<-matrix(0,n,1)
for (i in 1:n){
    j<-1
    while (j<=rnum){
         inside<-TRUE
         coordi<-1
         while ((inside) && (coordi<=d)){
             ala<-pcf$down[j,coordi]
             yla<-pcf$high[j,coordi]
             ala<-pcf$support[2*coordi-1]+ala*step[coordi]
             yla<-pcf$support[2*coordi-1]+yla*step[coordi]
             if ((dendat[i,coordi]<ala) || (dendat[i,coordi]>yla)) 
                         inside<-FALSE
             coordi<-coordi+1
         }
         if (inside){
            den2pcf[i]<-j
            pcf2den[j]<-i
            value[i]<-pcf$value[j]
         }
         j<-j+1
    }
}

datcol<-matrix("white",n,1)
for (i in 1:n){
    eka<-den2pcf[i]
    if (eka>0) tok<-nodefinder[eka]
    if ((eka>0)&&(tok>0)) datcol[i]<-col[tok]
}

or<-order(value,decreasing=FALSE)
return(list(datacolo=datcol,ord=or))
}



colorsofdata3<-function(dendat,pcf,lst,paletti=NULL, clusterlevel=NULL,nodes=NULL){
# this version written made by Sauli Herrala
# links from dendat to rec to node to color
# "lst$infopointer" gives links from nodes to recs

  n<-dim(dendat)[1]
  d<-dim(dendat)[2]
  rnum<-length(pcf$value)
  
  i <- 1:d
  step <-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
	  
  if (is.null(paletti))
  paletti<-c("red","blue","green",
      "orange","navy","darkgreen",
      "orchid","aquamarine","turquoise",
      "pink","violet","magenta","chocolate","cyan",
      colors()[50:657],colors()[50:657])
  
  # links from node to color
  if ((is.null(clusterlevel))&&(is.null(nodes))) col<-colobary(lst$parent,paletti)
  if (!is.null(clusterlevel)) col<-colobary.merge(lst$parent,lst$level,colothre=clusterlevel,paletti)
  if (!is.null(nodes)) col<-colobary.nodes(lst$parent,nodes,paletti)
  # links from rec to node (invert the links in infopointer)
 
  nodefinder<-matrix(0,rnum,1)
  for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i
  
  den2pcf<-matrix(0,n,1)
  pcf2den<-matrix(0,rnum,1)
  value<-matrix(0,n,1)
  ala <- pcf$down
  yla <- pcf$high
  alaTesti <- t(ala) * step + pcf$support[2*1:ncol(ala) -1]
  ylaTesti <- t(yla) * step + pcf$support[2*1:ncol(ala) -1]
  
  for (i in 1:n){
    bol <- (c(dendat[i, ]) < alaTesti) | (unlist(dendat[i,]) >  ylaTesti)
	j <- which.min(colSums(bol))
	den2pcf[i] <- j
    pcf2den[j] <- i
    value[i] <- pcf$value[j]	
  } 	
  
  datcol<-matrix("white",n,1)
  tok <- 0
  for (i in 1:n){
      eka<-den2pcf[i]
      if (eka>0) tok<-nodefinder[eka]
      if (tok>0) datcol[i]<-col[tok]
  }
  
  or<-order(value,decreasing=FALSE)
  return(list(datacolo=datcol,ord=or))
}


colorsofdata.adagrid2<-function(dendat,pcf,lst,paletti=NULL,
clusterlevel=NULL,nodes=NULL)
{
# links from dendat to rec to node to color
# "lst$infopointer" gives links from nodes to recs

n<-dim(dendat)[1]
d<-dim(dendat)[2]
rnum<-length(pcf$value)

if (is.null(paletti))
paletti<-c("red","blue","green",
    "orange","navy","darkgreen",
    "orchid","aquamarine","turquoise",
    "pink","violet","magenta","chocolate","cyan",
    colors()[50:657],colors()[50:657])

# links from node to color
if ((is.null(clusterlevel))&&(is.null(nodes))) col<-colobary(lst$parent,paletti)
if (!is.null(clusterlevel))
col<-colobary.merge(lst$parent,lst$level,colothre=clusterlevel,paletti)
if (!is.null(nodes))
col<-colobary.nodes(lst$parent,nodes,paletti)

# links from rec to node (invert the links in infopointer)
nodefinder<-matrix(0,rnum,1)
for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i

# find links from dendat to rec
den2pcf<-matrix(0,n,1)
pcf2den<-matrix(0,rnum,1)
value<-matrix(0,n,1)
for (i in 1:n){
    j<-1
    while (j<=rnum){
         inside<-TRUE
         coordi<-1
         while ((inside) && (coordi<=d)){
             ala<-pcf$down[j,coordi]
             yla<-pcf$high[j,coordi]
             ala<-pcf$grid[ala,coordi]   #pcf$support[2*coordi-1]+ala*step[coordi]
             yla<-pcf$grid[yla,coordi]   #pcf$support[2*coordi-1]+yla*step[coordi]
             if ((dendat[i,coordi]<ala) || (dendat[i,coordi]>yla)) 
                         inside<-FALSE
             coordi<-coordi+1
         }
         if (inside){
            den2pcf[i]<-j
            pcf2den[j]<-i
            value[i]<-pcf$value[j]
         }
         j<-j+1
    }
}

datcol<-matrix("white",n,1)
for (i in 1:n){
    eka<-den2pcf[i]
    if (eka>0) tok<-nodefinder[eka]
    if ((eka>0)&&(tok>0)) datcol[i]<-col[tok]
}

or<-order(value,decreasing=FALSE)
return(list(datacolo=datcol,ord=or))
}



colorsofdata.adagrid.new<-function(dendat, pcf, lst, paletti = NULL, clusterlevel=NULL, 
nodes=NULL){
  # links from dendat to rec to node to color
  # "lst$infopointer" gives links from nodes to recs

  n<-dim(dendat)[1]
  d<-dim(dendat)[2]
  rnum<-length(pcf$value)
  
  	  
  if (is.null(paletti)){
    paletti<-c("red","blue","green","orange","navy","darkgreen",
      "orchid","aquamarine","turquoise", "pink","violet","magenta",
      "chocolate","cyan", colors()[50:657],colors()[50:657])
  }
  # links from node to color
  if ((is.null(clusterlevel))&&(is.null(nodes))) col<-colobary(lst$parent,paletti)
  if (!is.null(clusterlevel)) col<-colobary.merge(lst$parent,lst$level,colothre=clusterlevel,paletti)
  if (!is.null(nodes)) col<-colobary.nodes(lst$parent,nodes,paletti)

  # links from rec to node (invert the links in infopointer)
  nodefinder<-matrix(0,rnum,1)
  for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i

  # find links from dendat to rec
  den2pcf<-matrix(0, n, 1)
  pcf2den<-matrix(0, rnum, 1)
  value<-matrix(0, n, 1)
  ala <- pcf$down
  yla <- pcf$high
  # a bit complex
  ala <- sapply(1:ncol(ala), function(x, pcf, ala) pcf$grid[ala[, x], x], pcf = pcf, ala = ala)
  yla <- sapply(1:ncol(yla), function(x, pcf, yla) pcf$grid[yla[, x], x], pcf = pcf, yla = yla)
  ala <- t(ala)
  yla <- t(yla)
  
  prc <- proc.time()
  for (i in 1:n){
    bol <- (c(dendat[i, ]) < ala) | (c(dendat[i,]) >  yla)
	j <- which.min(colSums(bol))
	den2pcf[i] <- j
    pcf2den[j] <- i
    value[i] <- pcf$value[j]	
  } 	
  datcol <- matrix("white",n,1)
  for (i in 1:n){
    eka<-den2pcf[i]
    if (eka>0) tok<-nodefinder[eka]
    if ((eka>0)&&(tok>0)) datcol[i]<-col[tok]
  } 
  or<-order(value,decreasing=FALSE)
  return(list(datacolo=datcol,ord=or))
}


colorsofdata.adagrid <- function(dendat, pcf, lst, paletti = NULL, clusterlevel=NULL, nodes=NULL){
  # links from dendat to rec to node to color
  # "lst$infopointer" gives links from nodes to recs
# this version written made by Sauli Herrala

  n<-dim(dendat)[1]
  d<-dim(dendat)[2]
  rnum<-length(pcf$value)
  
  	  
  if (is.null(paletti)){
    paletti<-c("red","blue","green","orange","navy","darkgreen",
      "orchid","aquamarine","turquoise", "pink","violet","magenta",
      "chocolate","cyan", colors()[50:657],colors()[50:657])
  }
  # links from node to color
  if ((is.null(clusterlevel))&&(is.null(nodes))) col<-colobary(lst$parent,paletti)
  if (!is.null(clusterlevel)) col<-colobary.merge(lst$parent,lst$level,colothre=clusterlevel,paletti)
  if (!is.null(nodes)) col<-colobary.nodes(lst$parent,nodes,paletti)

  # links from rec to node (invert the links in infopointer)
  nodefinder<-matrix(0,rnum,1)
  for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i

  # find links from dendat to rec
  den2pcf<-matrix(0, n, 1)
  pcf2den<-matrix(0, rnum, 1)
  value<-matrix(0, n, 1)
  ala <- pcf$down
  yla <- pcf$high
  # a bit complex
  ala <- sapply(1:ncol(ala), function(x, pcf, ala) pcf$grid[ala[, x], x], pcf = pcf, ala = ala)
  yla <- sapply(1:ncol(yla), function(x, pcf, yla) pcf$grid[yla[, x], x], pcf = pcf, yla = yla)
  ala <- t(ala)
  yla <- t(yla)
  
  prc <- proc.time()
  for (i in 1:n){
    bol <- (c(dendat[i, ]) < ala) | (c(dendat[i,]) >  yla)
	j <- which.min(colSums(bol))
	den2pcf[i] <- j
    pcf2den[j] <- i
    value[i] <- pcf$value[j]	
  } 	
  datcol <- matrix("white",n,1)
  for (i in 1:n){
    eka<-den2pcf[i]
    if (eka>0) tok<-nodefinder[eka]
    if ((eka>0)&&(tok>0)) datcol[i]<-col[tok]
  } 
  or<-order(value,decreasing=FALSE)
  return(list(datacolo=datcol,ord=or))
}

colorsofdata.new<-function(dendat, pcf, lst, paletti=NULL, clusterlevel=NULL, nodes=NULL){
# links from dendat to rec to node to color
# "lst$infopointer" gives links from nodes to recs
  n<-dim(dendat)[1]
  d<-dim(dendat)[2]
  rnum<-length(pcf$value)
  
  i <- 1:d
  step <-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
	  
  if (is.null(paletti))
  paletti<-c("red","blue","green",
      "orange","navy","darkgreen",
      "orchid","aquamarine","turquoise",
      "pink","violet","magenta","chocolate","cyan",
      colors()[50:657],colors()[50:657])
  
  # links from node to color
  if ((is.null(clusterlevel))&&(is.null(nodes))) col<-colobary(lst$parent,paletti)
  if (!is.null(clusterlevel)) col<-colobary.merge(lst$parent,lst$level,colothre=clusterlevel,paletti)
  if (!is.null(nodes)) col<-colobary.nodes(lst$parent,nodes,paletti)
  # links from rec to node (invert the links in infopointer)
 
  nodefinder<-matrix(0,rnum,1)
  for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i
  
  den2pcf<-matrix(0,n,1)
  pcf2den<-matrix(0,rnum,1)
  value<-matrix(0,n,1)
  ala <- pcf$down
  yla <- pcf$high
  alaTesti <- t(ala) * step + pcf$support[2*1:ncol(ala) -1]
  ylaTesti <- t(yla) * step + pcf$support[2*1:ncol(ala) -1]
  
  for (i in 1:n){
    bol <- (c(dendat[i, ]) < alaTesti) | (c(dendat[i,]) >  ylaTesti)
	j <- which.min(colSums(bol))
	den2pcf[i] <- j
    pcf2den[j] <- i
    value[i] <- pcf$value[j]	
  } 	
  
  datcol<-matrix("white",n,1)
  tok <- 0
  for (i in 1:n){
      eka<-den2pcf[i]
      if (eka>0) tok<-nodefinder[eka]
      if (tok>0) datcol[i]<-col[tok]
  }
  
  or<-order(value,decreasing=FALSE)
  return(list(datacolo=datcol,ord=or))
}





colorsofdata<-function(dendat, pcf, lst, paletti=NULL, clusterlevel=NULL, nodes=NULL)
{
# links from dendat to rec to node to color
# "lst$infopointer" gives links from nodes to recs
# this version written made by Sauli Herrala

  n<-dim(dendat)[1]
  d<-dim(dendat)[2]
  rnum<-length(pcf$value)
  
  i <- 1:d
  step <-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
	  
  if (is.null(paletti))
  paletti<-c("red","blue","green",
      "orange","navy","darkgreen",
      "orchid","aquamarine","turquoise",
      "pink","violet","magenta","chocolate","cyan",
      colors()[50:657],colors()[50:657])
  
  # links from node to color
  if ((is.null(clusterlevel))&&(is.null(nodes))) col<-colobary(lst$parent,paletti)
  if (!is.null(clusterlevel)) col<-colobary.merge(lst$parent,lst$level,colothre=clusterlevel,paletti)
  if (!is.null(nodes)) col<-colobary.nodes(lst$parent,nodes,paletti)
  # links from rec to node (invert the links in infopointer)
 
  nodefinder<-matrix(0,rnum,1)
  for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i
  
  den2pcf<-matrix(0,n,1)
  pcf2den<-matrix(0,rnum,1)
  value<-matrix(0,n,1)
  ala <- pcf$down
  yla <- pcf$high
  alaTesti <- t(ala) * step + pcf$support[2*1:ncol(ala) -1]
  ylaTesti <- t(yla) * step + pcf$support[2*1:ncol(ala) -1]
  
  for (i in 1:n){
    bol <- (c(dendat[i, ]) < alaTesti) | (c(dendat[i,]) >  ylaTesti)
	j <- which.min(colSums(bol))
	den2pcf[i] <- j
    pcf2den[j] <- i
    value[i] <- pcf$value[j]	
  } 	
  
  datcol<-matrix("white",n,1)
  tok <- 0
  for (i in 1:n){
      eka<-den2pcf[i]
      if (eka>0) tok<-nodefinder[eka]
      if (tok>0) datcol[i]<-col[tok]
  }
  
  or<-order(value,decreasing=FALSE)
  return(list(datacolo=datcol,ord=or))
}



colorsofdata.tail<-function(dendat,lst,paletti=NULL)
{
# links from dendat to node to color
# "lst$infopointer" gives links from nodes to data

n<-dim(dendat)[1]
d<-dim(dendat)[2]
    
if (is.null(paletti))
paletti<-c("red","blue","green",
    "orange","navy","darkgreen",
    "orchid","aquamarine","turquoise",
    "pink","violet","magenta","chocolate","cyan",
    colors()[50:657],colors()[50:657])

# links from node to color
col<-colobary(lst$parent,paletti)

# links from dendat to node (invert the links in infopointer)
nodefinder<-matrix(0,n,1)
for (i in 1:n) nodefinder[lst$infopointer[i]]<-i

datcol<-matrix("white",n,1)
for (i in 1:n){
    tok<-nodefinder[i]
    datcol[i]<-col[tok]
}

return(datacolo=datcol)
}



complex.rips<-function(dendat,rho)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

lkm<-0
complex<-matrix(0,1,d+1)

nn<-nn.indit(dendat)
maxk<-n-1
nnr<-nn.radit(dendat,maxk)

for (i in 1:n){
   rcur<-nnr[i,1]
   j<-1
   while ( (rcur<=rho) && (j<=(n-1)) ){
        cind<-nn[i,j]

        if (cind>i){
        # find the connection
        rcur1<-nnr[i,j]
        rcur2<-nnr[cind,1]
        kk<-j+1
        found<-FALSE
        while ( (rcur1<=rho) && (kk<=(n-1)) && (!found) ){
           koe1<-nn[i,kk]
           ll<-1
           while ( (rcur2<=rho) && (ll<=(n-1)) && (!found) ){
                koe2<-nn[cind,ll]
                if (koe1==koe2){
                     found<-TRUE
                     addi<-matrix(c(i,cind,koe1),1,3)
                     if (lkm==0) complex<-matrix(c(i,cind,koe1),1,3)
                     else complex<-rbind(complex,addi)
                     lkm<-lkm+1 
                } 
                ll<-ll+1 
                if (ll<=(n-1)) rcur2<-nnr[cind,ll]     
           }
           kk<-kk+1
           if (kk<=(n-1)) rcur1<-nnr[i,kk]
        }        
        # connection search end
        }

        j<-j+1
        if (j<=(n-1)) rcur<-nnr[i,j]
   }
}

if (lkm==0) complex<-NULL
return(complex)
}

cumu<-function(values,recs,frekv=NULL){
#Finds level sets of a piecewise constant function 
#
#values is recnum-vector
#recs is recnum*(2*d)-matrix
#frekv is recnum-vector
#
#returns list(levels,lsets,recs)
#   levels is levnum-vector,
#   lsets is levnum*atomnum-matrix,
#   atoms is recs but rows in different order
#   frekv is also only ordered differently

jarj<-omaord(values,recs,frekv)
values<-jarj$values
recs<-jarj$recs
frekv<-jarj$frekv
recnum<-length(values) #=length(recs[,1])#numb of lev is in the worst case
levels<-matrix(0,recnum,1)      
lsets<-matrix(1,recnum,recnum)  #same as the number of recs
       #at the beginning we mark everything belonging to level sets
       #next we start removing recs from level sets
levels[1]<-values[1]    #smallest values are first, first row of levels
                       #contains already 1:s 
curval<-values[1]
curlev<-1
for (i in 1:recnum){
  if (values[i]<=curval) lsets[(curlev+1):recnum,i]<-0
    else{
      curlev<-curlev+1      
      curval<-values[i]
      levels[curlev]<-values[i]
      if ((curlev+1)<=recnum) lsets[(curlev+1):recnum,i]<-0
    }
}
levels<-levels[1:curlev]
lsets<-lsets[1:curlev,]
return(list(levels=levels,lsets=lsets,atoms=recs,frekv=frekv))
}
cutmut<-function(mut,level,levels){
#
roots<-mut$roots
child<-mut$child
sibling<-mut$sibling
#
itemnum<-length(child)
rootnum<-length(roots)       
#
newroots<-matrix(0,itemnum,1)
newsibling<-sibling
ind<-0
#
for (i in 1:rootnum){
      curroot<-roots[i]
      pino<-matrix(0,itemnum,1)
      pino[1]<-curroot
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          # 
          # if cur acrosses the level, make cur root
          if (levels[cur]>level){
             ind<-ind+1
             newroots[ind]<-cur     # add to list
             newsibling[cur]<-0     # remove siblings
          }               
          # put to the stack
          if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
          }
          # go to left and put right nodes to the stack
          # if we have not already crossed the level
          while ((child[cur]>0) && (levels[cur]<=level)){
             cur<-child[cur]
             # if cur acrosses the level, make cur root
             if (levels[cur]>level){
                ind<-ind+1
                newroots[ind]<-cur  # add to list
                newsibling[cur]<-0     # remove siblings
             }  
             if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
             }
           }
       }
}                                        
if (ind>0){
    newroots<-newroots[1:ind]
}
else{
    newroots<-NULL
}
return(list(roots=newroots,sibling=newsibling))
}





cutvalue<-function(roots,child,sibling,level,component,
          AtomlistAtom,AtomlistNext,valnum){
#
#from the cutted multitree, form a "newvalue",
#which gives quantized values for the kernel estimate,
#in addition the values are cutted, so that one mode is 
#removed (input is cutted multitree)
#
itemnum<-length(child)
rootnum<-length(roots)
newvalue<-matrix(0,valnum,1)

for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    pino[1]<-roots[i]
    #    
    pinin<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        #
        node<-cur
        compo<-component[node]
        ato<-compo                          #ato is pointer to "value"
        while (ato>0){
           newvalue[AtomlistAtom[ato]]<-level[node]
           ato<-AtomlistNext[ato]
        }
        #
        if (sibling[cur]>0){
              pinin<-pinin+1
              pino[pinin]<-sibling[cur]
        }
        while (child[cur]>0){    #go to left and put right nodes to stack
              cur<-child[cur]
              #
              node<-cur
              compo<-component[node]
              ato<-compo                    #ato is pointer to "value"
              while (ato>0){
                  newvalue[AtomlistAtom[ato]]<-level[node]
                  ato<-AtomlistNext[ato]
              }
              #
              if (sibling[cur]>0){  #if cur has siblings
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
             }
        }
    }
}
#
return(newvalue)
}


cvolumbag<-function(component,AtomlistAtom,AtomlistNext,low,upp,steppi)
{
d<-dim(low)[2]

componum<-length(component)
volume<-matrix(0,componum,1)

for (i in 1:componum){
   curvolu<-0
   pointer<-component[i]
   while (pointer>0){
        atto<-AtomlistAtom[pointer]

        vol<-1
        for (j in 1:d){
            vol<-vol*(upp[atto,j]-low[atto,j])*steppi[j]
        }

        curvolu<-curvolu+vol
        pointer<-AtomlistNext[pointer]
   }
   volume[i]<-curvolu
}
return(volume)
}
cvolumdya<-function(volofatom,component,AtomlistNext){
#
componum<-length(component)
volume<-matrix(0,componum,1)
#
# it is enough to calculate the number of aoms in each component
for (i in 1:componum){
   numofatoms<-0
   pointer<-component[i]
   while (pointer>0){
        numofatoms<-numofatoms+1
        pointer<-AtomlistNext[pointer]
   }
   volume[i]<-numofatoms*volofatom
}
return(volume)
}
cvolum<-function(levels,items){
#Calculates volumes of set of level sets
#
#levels is tasolkm*N-matrix of 1:s and 0:s
#items is N*(2*d)-matrix
#
#returns N-vector of volumes
#
N<-length(levels[,1])
res<-matrix(0,N,1)
if (dim(t(levels))[1]==1) tasolkm<-1 else tasolkm<-length(levels[,1]) 
for (i in 1:tasolkm){
  lev2<-change(levels[i,])
  m<-length(lev2)
  vol<-0
  for (j in 1:m){
    ind<-lev2[j]
    rec<-items[ind,]
    vol<-vol+massone(rec)
  }
  res[i]<-vol
}
return(t(res))
}


 

declevdya<-function(beg,AtomlistNext,AtomlistAtom,kg,N,nodenumOfDyaker,
terminalnum,d){
#
#beg is pointer to AtomlistAtom
#nodenumOfDyaker is the num of nodes of the _original dyaker_
#terminalnum is the num of terminal nodes of the "current dyaker"
#
#kg=kernel estimate is represented as binary tree with "nodenumOfDyaker" nodes:
#   vectors whose length is "nodenum":
#     -left,right,parent
#     -infopointer: pointer to "value" and "index" (only for terminal nodes)
#   additional data structures  
#     -value, index
#
#to be created:
#-separy, vector with length "nodenumOfDyaker", points to begsSepaBegs
#-begsSepaBegs, begsSepaNext, begsLeftBoun, begsRighBoun
#   vectors of same length as "value" and "index",
#   list of separate sets, index gives starting point for set in atomsSepaNext
#-atomsSepaAtom, atomsSepaNext, atomsLBounNext, atomsRBounNext
#   vector of same length as "value" and "index",
#   list of atoms in separate sets, index gives the atom in value and index
#
#return: 
#begs: list of beginnings of lists
#AtomlistAtoms, AtomlistNext: list of lists of atoms
#
left<-kg$left
right<-kg$right
parent<-kg$parent
index<-kg$index
nodefinder<-kg$nodefinder
#infopointer<-kg$infopointer
#
separy<-matrix(0,nodenumOfDyaker,1)
#
begsSepaNext<-matrix(0,terminalnum,1)
begsSepaBegs<-matrix(0,terminalnum,1)   #pointers to begsSepaAtoms
begsLeftBoun<-matrix(0,terminalnum,1)
begsRighBoun<-matrix(0,terminalnum,1)
#
atomsSepaNext<-matrix(0,terminalnum,1)  #pointers to value,index
atomsSepaAtom<-matrix(0,terminalnum,1)
atomsLBounNext<-matrix(0,terminalnum,1)
atomsRBounNext<-matrix(0,terminalnum,1)
#
nextFloor<-matrix(0,terminalnum,1)
currFloor<-matrix(0,terminalnum,1)
already<-matrix(0,nodenumOfDyaker,1)
#
#############################
# INITIALIZING: "we go through the nodes at depth "depoftree""
# we make currFloor to be one over bottom floor and initialize 
# separy, boundary, atomsSepaAtom, atomsBounAtom
##############################
#
lkm<-0
curlkm<-0
curre<-beg
while(curre>0){
    lkm<-lkm+1
    atom<-AtomlistAtom[curre]
    node<-nodefinder[atom]
    #
    separy[node]<-lkm
    atomsSepaAtom[lkm]<-atom
    #
    exists<-parent[node]
    if (already[exists]==0){
         curlkm<-curlkm+1
         currFloor[curlkm]<-exists 
         already[exists]<-1
    }
    #
    curre<-AtomlistNext[curre]
}     # obs terminalnum=lkm
# initialize the rest
begsSepaBegs[1:terminalnum]<-seq(1,terminalnum)       
begsLeftBoun[1:terminalnum]<-seq(1,terminalnum)
begsRighBoun[1:terminalnum]<-seq(1,terminalnum)
#obs: we need not change 
#begsSepaNext, atomsSepaNext, atomsLBounNext, atomsRBounNext
#since at the beginning set consist only one member:
#pointer is always 0, since we do not have followers 
#
###########################
#
# START the MAIN LOOP
###########################
i<-d
while (i >= 2){
   j<-log(N[i],base=2)   #depth at direction d
   while (j>=1){
        nexlkm<-0
        k<-1
        while (k <= curlkm){
            node<-currFloor[k]  
            # we create simultaneously the upper floor
            exists<-parent[node]
            if (already[exists]==0){
                 nexlkm<-nexlkm+1
                 nextFloor[nexlkm]<-exists
                 already[exists]<-1
            }
####################################################
#            now we join childs            
####################################################
            leftbeg<-left[node]
            rightbeg<-right[node]
            direction<-i
            jg<-joingene(node,leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index)
#
separy<-jg$separy
begsSepaNext<-jg$begsSepaNext
begsSepaBegs<-jg$begsSepaBegs
begsLeftBoun<-jg$begsLeftBoun
begsRighBoun<-jg$begsRighBoun
atomsSepaNext<-jg$atomsSepaNext
atomsSepaAtom<-jg$atomsSepaAtom
atomsLBounNext<-jg$atomsLBounNext
atomsRBounNext<-jg$atomsRBounNext
           k<-k+1
        }
        j<-j-1
        curlkm<-nexlkm
        currFloor<-nextFloor
   }
   # now we move to the next direction, correct boundaries
   begsLeftBoun<-begsSepaBegs
   begsRighBoun<-begsSepaBegs
   #
   atomsLBounNext<-atomsSepaNext
   atomsRBounNext<-atomsSepaNext
   #
   i<-i-1
}
#########################
# ENO OF MAIN LOOP
#########################
#
#
###################
# LAST DIMENSION WILL BE handled, (because this contains root node
##################
   i<-1
   j<-log(N[i],base=2)   #depth at direction d
   while (j>=2){
        nexlkm<-0
        k<-1
        while (k <= curlkm){
            node<-currFloor[k]  
            # we create simultaneously the upper floor
            exists<-parent[node]
            if (already[exists]==0){
                 nexlkm<-nexlkm+1
                 nextFloor[nexlkm]<-exists
                 already[exists]<-1
            }
####################################################
#            now we join childs            
#if (right[parent[node]]==node)  #if node is right child 
####################################################
            leftbeg<-left[node]
            rightbeg<-right[node]
            direction<-1 
jg<-joingene(node,leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index)
#
separy<-jg$separy
begsSepaNext<-jg$begsSepaNext
begsSepaBegs<-jg$begsSepaBegs
begsLeftBoun<-jg$begsLeftBoun
begsRighBoun<-jg$begsRighBoun
#
atomsSepaNext<-jg$atomsSepaNext
atomsSepaAtom<-jg$atomsSepaAtom
atomsLBounNext<-jg$atomsLBounNext
atomsRBounNext<-jg$atomsRBounNext
           k<-k+1
        }
        j<-j-1
        curlkm<-nexlkm
        currFloor<-nextFloor
   }
#########################
# ROOT NODE, we do not anymore update boundaries
#########################
            k<-1
            node<-currFloor[k]  
  while (k <= 1){
####################################################
#            now we join childs            
#if (right[parent[node]]==node)  #if node is right child 
####################################################
            leftbeg<-left[node]
            rightbeg<-right[node]
            if ((leftbeg==0) || (separy[leftbeg]==0)){
                              #if left child does not exist    
                separy[node]<-separy[rightbeg]
            }
            else{   #eka else
                if ((rightbeg==0) || (separy[rightbeg]==0)){  
                              #right child does not exist  
                   separy[node]<-separy[leftbeg]
                } 
                else{   #toka else: both children exist
                    #check whether left boundary of right child is empty
                    Lempty<-TRUE
                    note<-separy[rightbeg]
                    while (note>0){
                        if (begsLeftBoun[note]>0){
                              Lempty<-FALSE
                        }
                        note<-begsSepaNext[note]
                     }
                     #check whether right bound of left child is empty     
                     Rempty<-TRUE
                     note<-separy[leftbeg]
                     while (note>0){
                          if (begsRighBoun[note]>0){
                                 Rempty<-FALSE
                          }
                          note<-begsSepaNext[note]
                     }
                     #check whether one of boundaries is empty
                     if (Lempty || Rempty){
                            #one of boundaries is empty         
############
#concatenating separate parts
#separy[node]<- concatenate separy[leftbeg],separy[rightbeg]
###########
akku<-separy[leftbeg]
while (begsSepaNext[akku]>0){
  akku<-begsSepaNext[akku]
}                           
begsSepaNext[akku]<-separy[rightbeg]
#
separy[node]<-separy[leftbeg]
####################
#end of concatenating, handle next boundaries
###################
                    }
                    else{  #both children exist, both boundaries non-empty  
direction<-i
jc<-joinconne(leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index)  #direction<-i
#
separy<-jc$separy
separy[node]<-jc$totbegSepary 
#
begsSepaNext<-jc$begsSepaNext
begsSepaBegs<-jc$begsSepaBegs
begsLeftBoun<-jc$begsLeftBoun
begsRighBoun<-jc$begsRighBoun
#
atomsSepaNext<-jc$atomsSepaNext
atomsSepaAtom<-jc$atomsSepaAtom
atomsLBounNext<-jc$atomsLBounNext
atomsRBounNext<-jc$atomsRBounNext
                        #
                    }
                } #toka else
            } #eka else
 k<-k+1
}
###########################################
#          end of child joining
###########################################
###################################
# END of ROOT
###################################
#
###########################
###########################
totbegSepary<-separy[node]
return(list(totbegSepary=totbegSepary,
begsSepaNext=begsSepaNext,begsSepaBegs=begsSepaBegs,
atomsSepaNext=atomsSepaNext,atomsSepaAtom=atomsSepaAtom))
}













declevgen<-function(tobehandled,curterminalnum,
left,right,val,vec,infopointer,parent,
low,upp,
d)
# source("~/denpro/R/declevgen.R")
# dl<-declevgen(tobehandled,curterminalnum,left,right,val,vec,
# infopointer,parent,low,upp,d)

{
nodelkm<-length(left)

separy<-matrix(0,nodelkm,1)          #=infopointer?

begsSepaNext<-matrix(0,curterminalnum,1)
begsSepaBegs<-matrix(0,curterminalnum,1)    #pointers to begsSepaAtom
begsLeftBoun<-matrix(0,curterminalnum,1)
begsRighBoun<-matrix(0,curterminalnum,1)

atomsSepaNext<-matrix(0,curterminalnum,1)   #pointers to value,index
atomsSepaAtom<-matrix(0,curterminalnum,1)
atomsLBounNext<-matrix(0,curterminalnum,1)
atomsRBounNext<-matrix(0,curterminalnum,1)

leafloc<-findleafs(left,right)

lkm<-0
node<-nodelkm
while (node>=1){ 
          # root is in position 1
    if ((leafloc[node]==1) && (tobehandled[node]==1)){   #we are in leaf 

          tobehandled[parent[node]]<-1

          lkm<-lkm+1
          separy[node]<-lkm
          atomsSepaAtom[lkm]<-infopointer[node]

          begsSepaBegs[lkm]<-lkm

          #begsLeftBoun[lkm]<-lkm
          #begsRighBoun[lkm]<-lkm

          #obs: we need not change
          #begsSepaNext, atomsSepaNext, atomsLBounNext, atomsRBounNext
          #since at the beginning set consist only one member:
          #pointer is always 0, since we do not have followers

    }
    else if (tobehandled[node]==1){   #not a leaf

        tobehandled[parent[node]]<-1

        leftbeg<-left[node]
        rightbeg<-right[node]

        if ((leftbeg==0) || (separy[leftbeg]==0)){
           #if left child does not exist

           #note that since we consider subsets of the
           #terminal nodes of the original tree, it may happen
           #that leftbeg>0 but left child does not exist
          
           separy[node]<-separy[rightbeg]

           #we need that all lists contain as many members
           #left boundary is empty, but we will make it a list
           #of empty lists

           #note<-separy[node]
           #while (note>0){
           #    begsLeftBoun[note]<-0
           #    note<-begsSepaNext[note]
           #}
           # right boundary stays same as for rightbeg
        }
        else{   # eka else
            if ((rightbeg==0) || (separy[rightbeg]==0)){
            #right child does not exist
    
                   separy[node]<-separy[leftbeg]
    
                   # left boundary stays same as for leftbeg
                   # right boundary is empty
    
                   #note<-separy[node]
                   #while (note>0){
                   #    begsRighBoun[note]<-0
                   #    note<-begsSepaNext[note]
                   #}
            }
            else{   #toka else: both children exist
                    #create boundaries 
                    
                    direktiooni<-vec[node]
                    splittiini<-val[node]     
               
                    #left boundary of right child :
                    #create/check whether empty 
                   
                    Lempty<-TRUE
                    note<-separy[rightbeg]
                    while (note>0){
                        thisnoteempty<-TRUE

                        poiju<-begsSepaBegs[note]
                        prevpoiju<-poiju 
                        while (poiju>0){
                           aatto<-atomsSepaAtom[poiju]
                           if (!(splittiini<low[aatto,direktiooni])){
                               #this atom belongs to boundary
                               if (thisnoteempty==TRUE){
                                   #poiju is the 1st non-empty
                                   begsLeftBoun[note]<-poiju
                               }                
                               Lempty<-FALSE   
                               atomsLBounNext[prevpoiju]<-poiju    
                               atomsLBounNext[poiju]<-0
                               prevpoiju<-poiju

                               thisnoteempty<-FALSE
                           }
                           poiju<-atomsSepaNext[poiju]
                        }
                        if (thisnoteempty) begsLeftBoun[note]<-0
                        note<-begsSepaNext[note]
                     }

                     #right boundary of left child

                     Rempty<-TRUE
                     note<-separy[leftbeg]
                     while (note>0){
                        thisnoteempty<-TRUE
                       
                        poiju<-begsSepaBegs[note]
                        prevpoiju<-poiju 
                        while (poiju>0){
                           aatto<-atomsSepaAtom[poiju]
                           if (!(splittiini>upp[aatto,direktiooni])){
                               #this atom belongs to boundary
                               if (thisnoteempty==TRUE){
                                   #poiju is the 1st non-empty
                                   begsRighBoun[note]<-poiju
                               }                
                               Rempty<-FALSE   
                               atomsRBounNext[prevpoiju]<-poiju    
                               atomsRBounNext[poiju]<-0
                               prevpoiju<-poiju 
                           
                               thisnoteempty<-FALSE
                           }
                           poiju<-atomsSepaNext[poiju]
                        }
                        if (thisnoteempty) begsRighBoun[note]<-0
                        note<-begsSepaNext[note]
                     }

                     #check whether one of boundaries is empty

                     if (Lempty || Rempty){
                        #one of boundaries is empty
                        ############
                        #concatenating separate parts
                        ############
                        akku<-separy[leftbeg]
                        #begsRighBoun[akku]<-0 
                        # right boundaries of sets in left child are empty
                        # begsLeftBoun[akku] does not change
                        
                        #while (begsSepaNext[akku]>0){
                        #    akku<-begsSepaNext[akku]
                        #    begsRighBoun[akku]<-0
                        #}
                        begsSepaNext[akku]<-separy[rightbeg] 
                        #concatenate list of separate sets

                        separy[node]<-separy[leftbeg]
                        akku<-separy[rightbeg]
                          
                        #begsLeftBoun[akku]<-0 
                        #left boundaries of sets in right child are empty

                        #while (begsSepaNext[akku]>0){
                        #   akku<-begsSepaNext[akku]
                        #   begsLeftBoun[akku]<-0
                        #}
                        ############# 
                        #end of concatenating
                        #############
                     }
                     else{  #both children exist, both boundaries non-empty
                     direction<-vec[node]
                     jc<-joincongen(leftbeg,rightbeg,separy,
                     begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
                     atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
                     direction,low,upp)   
                      
                     separy<-jc$separy
                     separy[node]<-jc$totbegSepary

                     begsSepaNext<-jc$begsSepaNext
                     begsSepaBegs<-jc$begsSepaBegs
                     #begsLeftBoun<-jc$begsLeftBoun
                     #begsRighBoun<-jc$begsRighBoun

                     atomsSepaNext<-jc$atomsSepaNext
                     atomsSepaAtom<-jc$atomsSepaAtom
                     #atomsLBounNext<-jc$atomsLBounNext
                     #atomsRBounNext<-jc$atomsRBounNext
                     }
            } #toka else
        } #eka else
    }  #else not a leaf
    node<-node-1
}

return(list(
totbegSepary=separy[1],
begsSepaNext=begsSepaNext,
begsSepaBegs=begsSepaBegs,
atomsSepaNext=atomsSepaNext,
atomsSepaAtom=atomsSepaAtom
))
}













declevnew<-function(rindeksit,linkit,n){
#Splits level set to disconnected subsets
#
#rindeksit is m-vector, links to observations, levset is union of m atoms
#linkit is n*n-matrix, n is the total number of atoms,
#  describes which atoms touch each other
#
#Returns sublevnum*n-matrix, describes disconnected parts of levset  
#
m<-length(rindeksit)
tulos<-matrix(0,m,n)      #in the worst case we split to m parts
#                         #blokit !!!!!!!!!!!!!!!!
merkatut<-matrix(0,m,1) #laitetaan 1 jos on jo sijoitettu johonkin tasojouk.
pino<-matrix(0,m,1) #pinoon laitetaan aina jos koskettaa, max kosketuksia m
pinind<-0           #pinossa viitataan rindeksin elementteihin
curleima<-1
i<-1   #i ja j viittavat rindeksit-vektoriin, jonka alkiot viittavat atomeihin
while (i<=m){
  if (merkatut[i]==0){  #jos ei viela merkattu niin 
    pinind<-pinind+1    #pannaan pinoon
    pino[pinind]<-i    
    while (pinind>0){
      curviite<-pino[pinind]  #otetaan pinosta viite rindeksit-vektoriin
                              #jossa puolestaan viitteet itse palloihin
      curpallo<-rindeksit[curviite]  #haetaan viite palloon
      pinind<-pinind-1  
      tulos[curleima,curpallo]<-1    #laitetaan pallo ko tasojoukkoon
#      merkatut[curviite]<-1          #merkataan kaytetyksi
      j<-1
      while (j<=m){          #pannnaan linkeista pinoon
        ehdokas<-rindeksit[j]         #kaydaan ko tasojoukon atomit lapi
        touch<-(linkit[curpallo,ehdokas]==1)
        if ((touch) && (merkatut[j]==0)){
          pinind<-pinind+1      
          pino[pinind]<-j  
          merkatut[j]<-1
        }
        j<-j+1
      }
    }
    curleima<-curleima+1   #uusi leima
  }
  i<-i+1
} 
tullos<-matrix(0,(curleima-1),n)
tullos[1:(curleima-1),]<-tulos[1:(curleima-1),] #poistetaan ylimaaraiset
return(tullos)
}


















declev<-function(rindeksit,linkit,n){
#Splits level set to disconnected subsets
#
#rindeksit is m-vector, links to observations, levset is union of m atoms
#linkit is n*n-matrix, n is the total number of atoms,
#  describes which atoms touch each other
#
#Returns sublevnum*n-matrix, describes disconnected parts of levset  
#
m<-length(rindeksit)
tulos<-matrix(0,m,n)      #in the worst case we split to m parts
#                         #blokit !!!!!!!!!!!!!!!!
merkatut<-matrix(0,m,1) #laitetaan 1 jos on jo sijoitettu johonkin tasojouk.
pino<-matrix(0,m,1) #pinoon laitetaan aina jos koskettaa, max kosketuksia m
pinind<-0           #pinossa viitataan rindeksin elementteihin
curleima<-1
i<-1   #i ja j viittavat rindeksit-vektoriin, jonka alkiot viittavat atomeihin
while (i<=m){
  if (merkatut[i]==0){  #jos ei viela merkattu niin 
    pinind<-pinind+1    #pannaan pinoon
    pino[pinind]<-i    
    while (pinind>0){
      curviite<-pino[pinind]  #otetaan pinosta viite rindeksit-vektoriin
                              #jossa puolestaan viitteet itse palloihin
      curpallo<-rindeksit[curviite]  #haetaan viite palloon
      pinind<-pinind-1  
      tulos[curleima,curpallo]<-1    #laitetaan pallo ko tasojoukkoon
#      merkatut[curviite]<-1          #merkataan kaytetyksi
      j<-1
      while (j<=m){          #pannnaan linkeista pinoon
        curkoske<-linkit[curpallo,]   #ne joihin curpallo koskettaa
        ehdokas<-rindeksit[j]         #kaydaan ko tasojoukon atomit lapi
        touch<-onko(curkoske,ehdokas) #onko ehdokas rivilla "curkoske"
            #touch<-(linkit[curpallo,rindeksit[j]]==1)
        if ((touch) && (merkatut[j]==0)){
          pinind<-pinind+1      
          pino[pinind]<-j  
          merkatut[j]<-1
        }
        j<-j+1
      }
    }
    curleima<-curleima+1   #uusi leima
  }
  i<-i+1
} 
tullos<-matrix(0,(curleima-1),n)
tullos[1:(curleima-1),]<-tulos[1:(curleima-1),] #poistetaan ylimaaraiset
return(tullos)
}


















decombag<-function(numofall,levseq,
left,right,val,vec,infopointer,parent,
nodenumOfTree,terminalnum,
value,low,upp,nodefinder,
d)
{
#Makes density tree
#returns list(parent,level,component)
#  -component is pointer to AtomlistAtom, AtomlistNext, ei tarvita

#apu2<-0

levnum<-length(levseq)

AtomlistNext<-matrix(0,numofall,1)
AtomlistAtom<-matrix(0,numofall,1) #point to value,..: values in 1,...,atomnum

componum<-numofall
parentLST<-matrix(0,componum,1)
level<-matrix(0,componum,1)
component<-matrix(0,componum,1)    

pinoComponent<-matrix(0,componum,1)   #pointer to component, level,...
pinoTaso<-matrix(0,componum,1)        #ordinal of level (pointer to levseq)

# Initilize the lists

AtomlistAtom[1:terminalnum]<-seq(1,terminalnum)
AtomlistNext[1:(terminalnum-1)]<-seq(2,terminalnum)
AtomlistNext[terminalnum]<-0
listEnd<-terminalnum
beg<-1

# Let us divide the lowest level set to disconnected parts

begi<-1
tobehandled<-matrix(0,nodenumOfTree,1)
atto<-AtomlistAtom[begi]
while (begi>0){
   if (value[atto]>0){
      node<-nodefinder[atto]
      tobehandled[node]<-1
   }    
   begi<-AtomlistNext[begi]
   atto<-AtomlistAtom[begi]
}
curterminalnum<-terminalnum

dld<-declevgen(tobehandled,curterminalnum,
left,right,val,vec,infopointer,parent,
low,upp,
d)
 
totbegSepary<-dld$totbegSepary
begsSepaNext<-dld$begsSepaNext
begsSepaBegs<-dld$begsSepaBegs
atomsSepaNext<-dld$atomsSepaNext
atomsSepaAtom<-dld$atomsSepaAtom

lc<-listchange(AtomlistAtom,AtomlistNext,totbegSepary,
begsSepaNext,begsSepaBegs,atomsSepaNext,atomsSepaAtom,
curterminalnum,beg)
begs<-lc$begs
AtomlistNext<-lc$AtomlistNext
AtomlistAtom<-lc$AtomlistAtom

koko<-length(begs)
# Talletetaan osat 
component[1:koko]<-begs
level[1:koko]<-levseq[1]       #arvo toistuu 
efek<-koko                     #kirjataan uusien osien lkm  ????? jos vain yksi
# Laitetaan kaikki osat pinoon
pinoComponent[1:koko]<-seq(1,koko)      #1,2,...,koko
pinoTaso[1:koko]<-1     #kaikki osat kuuluvat alimpaan tasojoukkoon
pinind<-koko            #indeksi pinoon

if (levnum>1){  while (pinind>=1){
  # Take from stack
  ind<-pinoComponent[pinind]      #indeksi tasoon
  levind<-pinoTaso[pinind]        #ko tason korkeus
  pinind<-pinind-1                #otettu pinosta
  partlevsetbeg<-component[ind]  
  higlev<-levseq[levind+1]
  
  # Make intersection with the curr. component and higher lev.set
  PrevlistEnd<-listEnd
  addnum<-0      #num of atoms in the intersection
  removenum<-0   #num of atoms which have to be removed to get intersec.
  runner<-partlevsetbeg
  origiListEnd<-listEnd
  while ((runner>0) && (runner<=origiListEnd)){
      atom<-AtomlistAtom[runner]
      arvo<-value[atom]
      if (arvo>=higlev){
          listEnd<-listEnd+1    
          AtomlistAtom[listEnd]<-atom
          AtomlistNext[listEnd]<-listEnd+1 
          addnum<-addnum+1     
      }
      else{           
          removenum<-removenum+1
      }                                
      runner<-AtomlistNext[runner] 
  }
  AtomlistNext[listEnd]<-0      # we have to correct the end to terminate
  if (addnum>0){
      AtomlistNext[PrevlistEnd]<-0
      beghigher<-PrevlistEnd+1 
  }
  if (removenum==0){  #jos leikkaus ei muuta, niin tasoj sailyy samana  
      level[ind]<-levseq[levind+1]  #remove lower part
          #component and parentLST stay same, it is enough to change level
      if (levind+1<levnum){ #jos ei olla korkeimmalla tasolla,laita pinoon
        pinoComponent[pinind+1]<-ind     
        pinoTaso[pinind+1]<-levind+1 #tasojouk taso on levind+1 
        pinind<-pinind+1       
      }
  }
  else if (addnum>0){     #leikkaus ei tyhja
      beg<-beghigher

      begi<-beghigher
      tobehandled<-matrix(0,nodenumOfTree,1)
      while (begi>0){
          atto<-AtomlistAtom[begi]
          node<-nodefinder[atto]
          tobehandled[node]<-1
          begi<-AtomlistNext[begi]
      }
      curterminalnum<-addnum

      dld<-declevgen(tobehandled,curterminalnum,
      left,right,val,vec,infopointer,parent,
      low,upp,
      d)
 
      totbegSepary<-dld$totbegSepary
      begsSepaNext<-dld$begsSepaNext
      begsSepaBegs<-dld$begsSepaBegs
      atomsSepaNext<-dld$atomsSepaNext
      atomsSepaAtom<-dld$atomsSepaAtom

#      apu2<-apu2+1
#if (apu2==1) an<-atomsSepaNext     

      lc<-listchange(AtomlistAtom,AtomlistNext,totbegSepary,
      begsSepaNext,begsSepaBegs,atomsSepaNext,atomsSepaAtom,
      curterminalnum,beg)
      begs<-lc$begs
      AtomlistNext<-lc$AtomlistNext
      AtomlistAtom<-lc$AtomlistAtom

      koko<-length(begs)    #jos vain yksi ?????????
      #
      # paivitetaan kumu tulokseen
      level[(efek+1):(efek+koko)]<-levseq[levind+1]   #arvo toistuu  
      parentLST[(efek+1):(efek+koko)]<-ind
      component[(efek+1):(efek+koko)]<-begs
      efek<-efek+koko
      if (levind+1<levnum){ #jos ei olla korkeimmalla tasolla,laita pinoon
        pinoComponent[(pinind+1):(pinind+koko)]<-seq(efek-koko+1,efek)
        pinoTaso[(pinind+1):(pinind+koko)]<-levind+1  #tasjouk tas on levind+1 
        pinind<-pinind+koko       
      }
  }
}} 
level<-t(level[1:efek])
parentLST<-t(parentLST[1:efek])
component<-t(component[1:efek])
#
return(list(level=level,parent=parentLST,
component=component,AtomlistAtom=AtomlistAtom,AtomlistNext=AtomlistNext))
}









decomdya<-function(numofall,atomnum,levseq,kg,N,nodenumOfDyaker){
#Makes density tree
#
#returns list(parent,level,component)
#  -component is pointer to AtomlistAtom, AtomlistNext, ei tarvita
#
d<-length(N)
levnum<-length(levseq)
#
AtomlistNext<-matrix(0,numofall,1)
AtomlistAtom<-matrix(0,numofall,1) #point to value,..: values in 1,...,atomnum
#
componum<-numofall
parent<-matrix(0,componum,1)
level<-matrix(0,componum,1)
component<-matrix(0,componum,1)    
#
pinoComponent<-matrix(0,componum,1)   #pointer to component, level,...
pinoTaso<-matrix(0,componum,1)        #ordinal of level (pointer to levseq)
#
# Initilize the lists
#
AtomlistAtom[1:atomnum]<-seq(1,atomnum)
AtomlistNext[1:(atomnum-1)]<-seq(2,atomnum)
AtomlistNext[atomnum]<-0
listEnd<-atomnum
#
# Let us divide the lowest level set to disconnected parts
#
beg<-1
terminalnum<-atomnum
dld<-declevdya(beg,AtomlistNext,AtomlistAtom,kg,N,nodenumOfDyaker,
terminalnum,d)  #terminalnum<-atomnum 
totbegSepary<-dld$totbegSepary
begsSepaNext<-dld$begsSepaNext
begsSepaBegs<-dld$begsSepaBegs
atomsSepaNext<-dld$atomsSepaNext
atomsSepaAtom<-dld$atomsSepaAtom
#
lc<-listchange(AtomlistAtom,AtomlistNext,totbegSepary,
begsSepaNext,begsSepaBegs,atomsSepaNext,atomsSepaAtom,
terminalnum,beg)
begs<-lc$begs
AtomlistNext<-lc$AtomlistNext
AtomlistAtom<-lc$AtomlistAtom
#
koko<-length(begs)
# Talletetaan osat 
component[1:koko]<-begs
level[1:koko]<-levseq[1]       #arvo toistuu 
efek<-koko                     #kirjataan uusien osien lkm  ????? jos vain yksi
# Laitetaan kaikki osat pinoon
pinoComponent[1:koko]<-seq(1,koko)      #1,2,...,koko
pinoTaso[1:koko]<-1     #kaikki osat kuuluvat alimpaan tasojoukkoon
pinind<-koko            #indeksi pinoon
# 
if (levnum>1){  while (pinind>=1){
  # Take from stack
  ind<-pinoComponent[pinind]      #indeksi tasoon
  levind<-pinoTaso[pinind]        #ko tason korkeus
  pinind<-pinind-1                #otettu pinosta
  partlevsetbeg<-component[ind]  
  higlev<-levseq[levind+1]
  #
  # Make intersection with the curr. component and higher lev.set
  PrevlistEnd<-listEnd
  addnum<-0      #num of atoms in the intersection
  removenum<-0   #num of atoms which have to be removed to get intersec.
  runner<-partlevsetbeg
  origiListEnd<-listEnd
  value<-kg$value
  while ((runner>0) && (runner<=origiListEnd)){
      atom<-AtomlistAtom[runner]
      arvo<-value[atom]
      if (arvo>=higlev){
          listEnd<-listEnd+1    
          AtomlistAtom[listEnd]<-atom
          AtomlistNext[listEnd]<-listEnd+1 
          addnum<-addnum+1     
      }
      else{           
          removenum<-removenum+1
      }                                
      runner<-AtomlistNext[runner] 
  }
  AtomlistNext[listEnd]<-0      # we have to correct the end to terminate
  if (addnum>0){
      AtomlistNext[PrevlistEnd]<-0
      beghigher<-PrevlistEnd+1 
  }
  if (removenum==0){  #jos leikkaus ei muuta, niin tasoj sailyy samana  
      level[ind]<-levseq[levind+1]  #remove lower part
          #component and parent stay same, it is enough to change level
      if (levind+1<levnum){ #jos ei olla korkeimmalla tasolla,laita pinoon
        pinoComponent[pinind+1]<-ind     
        pinoTaso[pinind+1]<-levind+1 #tasojouk taso on levind+1 
        pinind<-pinind+1       
      }
  }
  else if (addnum>0){     #leikkaus ei tyhja
      beg<-beghigher
      terminalnum<-addnum
      dld<-declevdya(beg,AtomlistNext,AtomlistAtom,kg,N,nodenumOfDyaker,
                     terminalnum,d) 
totbegSepary<-dld$totbegSepary
begsSepaNext<-dld$begsSepaNext
begsSepaBegs<-dld$begsSepaBegs
atomsSepaNext<-dld$atomsSepaNext
atomsSepaAtom<-dld$atomsSepaAtom
#
lc<-listchange(AtomlistAtom,AtomlistNext,totbegSepary,
begsSepaNext,begsSepaBegs,atomsSepaNext,atomsSepaAtom,
terminalnum,beg)
begs<-lc$begs
AtomlistNext<-lc$AtomlistNext
AtomlistAtom<-lc$AtomlistAtom
#
      koko<-length(begs)    #jos vain yksi ?????????
      #
      # paivitetaan kumu tulokseen
      level[(efek+1):(efek+koko)]<-levseq[levind+1]   #arvo toistuu  
      parent[(efek+1):(efek+koko)]<-ind
      component[(efek+1):(efek+koko)]<-begs
      efek<-efek+koko
      if (levind+1<levnum){ #jos ei olla korkeimmalla tasolla,laita pinoon
        pinoComponent[(pinind+1):(pinind+koko)]<-seq(efek-koko+1,efek)
        pinoTaso[(pinind+1):(pinind+koko)]<-levind+1  #tasjouk tas on levind+1 
        pinind<-pinind+koko       
      }
  }
}} 
level<-t(level[1:efek])
parent<-t(parent[1:efek])
component<-t(component[1:efek])
#
return(list(level=level,parent=parent,
component=component,AtomlistAtom=AtomlistAtom,AtomlistNext=AtomlistNext))
}









decom<-function(lsets,levels,links,alkublokki,blokki){
#Makes density tree,  edellyttaa jarjestyksen ???????
#
#lsets is levnum*atomnum-matrix
#levels is levnum-vector
#links is  atomnum*maxtouchnum-matrix, pointers to atoms,
#  for each atom we indicate which atoms it touches
#
#lsets matriisiin lisataan riveja, silla aikaisemmat rivit jakautuvat
#erillisiin osiinsa, lisataan levels:n vastaavat alkiot 
#
#Returns list(lsets,levels,parents)
# parents and levels are newlevnum-vectors
# lsets is newlevnum*atomnum
#
if (dim(t(lsets))[1]==1) levnum<-1 else levnum<-length(lsets[,1])
                                              #rows of lsets 
if (levnum==1) atomnum<-length(lsets) else
atomnum<-length(lsets[1,])     #maxalkio on n
newlevnum<-levnum*atomnum       #karkea arvio, jokainen tasojoukko 
#                               #voi periaatteessa jakautua n:aan osaan
newlsets<-matrix(0,alkublokki,atomnum) 
newlevels<-matrix(0,alkublokki,1)
parents<-matrix(0,alkublokki,1)
#
curblokki<-alkublokki
#
pino<-matrix(0,newlevnum,2) #1.col indeksi newlsets:iin, 2.col ind tasoon
a<-1
b<-2
#
if (levnum==1) levset<-lsets else
levset<-lsets[1,]                   #alin tasojoukko
rindeksit<-change(levset)            #change the representation
kumu<-declev(rindeksit,links,atomnum)  #jaetaan alin tasoj. osiin
if (dim(t(kumu))[1]==1) koko<-1 else koko<-length(kumu[,1]) #osien lkm
# Talletetaan osat 
newlsets[1:koko,]<-kumu   
newlevels[1:koko]<-levels[1]     #arvo toistuu 
efek<-koko                       #kirjataan uusien osien lkm
# Laitetaan kaikki osat pinoon
pino[1:koko,a]<-seq(1,koko)      #1,2,...,koko
pino[1:koko,b]<-1   #kaikki osat kuuluvat alimpaan tasojoukkoon
pinind<-koko        #indeksi pinoon
# 
if (levnum>1){  while (pinind>=1){
  ind<-pino[pinind,a]         #indeksi tasoon
  levind<-pino[pinind,b]      #ko tason korkeus
  pinind<-pinind-1            #otettu pinosta
  partlevset<-newlsets[ind,]
  higlevset<-lsets[levind+1,] #huom levind<levnum
  levset<-partlevset*higlevset
  if (sum(partlevset-levset)==0){  
          #jos leikkaus ei muuta, niin tasoj sailyy samana  
      newlevels[ind]<-levels[levind+1]  #poistetaan alempi osa      
      if (levind+1<levnum){ #jos ei olla korkeimmalla tasolla,laita pinoon
        pino[pinind+1,a]<-ind     
        pino[pinind+1,b]<-levind+1 #tasojouk taso on levind+1 
        pinind<-pinind+1       
      }
  }
  else if (sum(levset)>0){   #leikkaus ei tyhja
      rindeksit<-change(levset)
      kumu<-declev(rindeksit,links,atomnum) #jaet tasoj. osa osiin
      if (dim(t(kumu))[1]==1) koko<-1 else koko<-length(kumu[,1]) #osien lkm
      if ((efek+koko)>curblokki){   
        newlsets<-blokitus(newlsets,blokki)
        newlevels<-blokitus(newlevels,blokki)
        parents<-blokitus(parents,blokki)
        curblokki<-curblokki+blokki    
      }
      newlsets[(efek+1):(efek+koko),]<-kumu  #paivitetaan kumu tulokseen
      newlevels[(efek+1):(efek+koko),]<-levels[levind+1]   #arvo toistuu  
      parents[(efek+1):(efek+koko),]<-ind
      efek<-efek+koko
      if (levind+1<levnum){ #jos ei olla korkeimmalla tasolla,laita pinoon
        pino[(pinind+1):(pinind+koko),a]<-seq(efek-koko+1,efek)
        pino[(pinind+1):(pinind+koko),b]<-levind+1 #tasojouk taso on levind+1 
        pinind<-pinind+koko       
      }
  }
}} 
newlevels<-t(newlevels[1:efek])
newlsets<-newlsets[1:efek,]
parents<-t(parents[1:efek])
return(list(lsets=newlsets,levels=newlevels,parents=parents))
}








den2reg<-function(dendat,binlkm,kantaja)
{
#muuntaa tiheysdatan regressiodataksi

#dendat on tiheysdatan sisaltava n*xlkm matriisi
#binlkm on niitten lokeroitten maara joihin yksi muuttuja ositetaan
#otosavaruus ositetaan binlkm^p lokeroon
#kantaja on 2*xlkm vaakavektori, sisaltaa kantajan kullekin 
#muuttujalle, olet etta kaikki dendat:n hav todella sisaltyvat kantajaan.

#palauttaa: listan vast
#vast.dep on saatu*1 vektori, sisaltaa frekvenssit,
#(saatu on erillisten diskretoitujen havaintojen lkm) 
#vast.ind on saatu*xlkm matriisi, sisaltaa niitten lokeroitten 
#keskipisteet joita (diskretoidussa) datassa esiintyy vahintaan yksi, 
#vast.hila on xlkm*3 matriisi, 
#vast.hila:n i:s rivi sis. i:nnelle muuttujalle 
#pisteiden maaran
#ensimmaisen regressiodatan havaintopisteen,
#viimeisen regressiodatan havaintopisteen.

xlkm<-length(dendat[1,]) #dendat:in sarakkeitten lkm on muuttujien lkm
n<-length(dendat[,1])    #dendat:in rivien lkm on havaintojen maara
hila<-matrix(1,xlkm,3)   #hila on xlkm*3 matriisi
binhila<-hila            #binhila on xlkm*3 matriisi
valpit<-matrix(1,xlkm,1)      #hilan valien pituudet
hila[,1]<-binlkm*hila[,1]
binhila[,1]<-binlkm*binhila[,1]
i<-1
while (i<=xlkm){
     binhila[i,2]<-kantaja[i,1]   #min(dendat[,i])-epsi 
     binhila[i,3]<-kantaja[i,2]     #max(dendat[,i])+epsi 
     valpit[i]<-(binhila[i,3]-binhila[i,2])/binlkm
       #binhila:n i:s rivi sis. i:nnelle muuttujalle 
       #lokeroinnin alkupisteen, arvoalueen loppupisteen   
       #valpit: arvoalueen pituus jaettuna lokeroitten lkm:lla
       #eli yhden lokeron leveys. 
     i<-i+1
}
hila[,2]<-binhila[,2]+valpit/2
hila[,3]<-binhila[,3]-valpit/2
       #hila:n i:s rivi sis. i:nnelle muuttujalle 
       #ensimmaisen regressiodatan havaintopisteen,
       #viimeisen regressiodatan havaintopisteen
#if (valpit<=0) stop("in some variable there is no variation")
hiladat<-matrix(1,n,xlkm) #muunnetaan dendat hiladat:iksi
                          #ts. pyoristetaan havainnot 
                          #hiladat sis. diskretoidut havainnot
i<-1
while (i<=n){       #kaydaan lapi aineisto
    #pyoristetaan i:s havainto hilapisteeseen
    j<-1
    while (j<=xlkm){
       alavali<-floor((dendat[i,j]-binhila[j,2])/valpit[j])
           #alavali ilmaisee monennessako lokerossa hav. sijaitsee
       hiladat[i,j]<-binhila[j,2]+alavali*valpit[j]+valpit[j]/2
       j<-j+1
    }
i<-i+1
}
xtulos<-matrix(0,n,xlkm) #periaatteessa mahdollista etta kaikki n
                         #havaintoa ovat eri lokeroissa, siksi
                         #laitetaan xtulos matriisin rivien maaraksi n
ytulos<-matrix(0,n,1)
xtulos[1,]<-hiladat[1,]  #hiladat:in ensimmainen rivi esiintyy ainakin kerran
ytulos[1]<-1             #sen frekvenssi ainakin yksi
saatu<-1                 #toistaiseksi yksi erillinen havainto
i<-1
while (i<n){  #kaydaan lapi aineisto
   i<-i+1     #aloitetaan kakkosesta
   lippu<-0   #apriori kyseessa uusi lajityyppi
   j<-1
   while ((j<=saatu) && (lippu==0)){  #kaydaan lapi keratyt havinnot
       if (all(hiladat[i,]==xtulos[j,])){  #jos on jo saatu
            lippu<-1     #liputetaan etta havaittiin toisto
            jind<-j      #merkataan indeksi frekvenssin paivitysta varten
       }
       j<-j+1       
   }
   if (lippu==1) ytulos[jind]<-ytulos[jind]+1 
            #jos saatiin toisto, paivitetaan frekvenssi 
      else{ 
         saatu<-saatu+1      #jos saatiin uusi, lisataan saatu:un yksi ja
         xtulos[saatu,]<-hiladat[i,]  #merkitaan uusi lajityyppi muistiin
         ytulos[saatu]<-1    #uuden lajityypin frekvenssi on aluksi yksi
      }
}
xtulos<-xtulos[1:saatu,]
ytulos<-ytulos[1:saatu]
ytulos<-t(t(ytulos))
if (xlkm==1) xtulos<-t(t(xtulos))
return(list(dep=ytulos,ind=xtulos,hila=hila))
}

dend2parent<-function(hc,dendat)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]
nodenum<-length(hc$height)+n
parent<-matrix(0,nodenum,1)
volume<-matrix(0,nodenum,1)
center<-matrix(0,d,nodenum)
level<-matrix(0,nodenum,1)
level[(n+1):nodenum]<-hc$height
pointers<-matrix(0,n,1)      #pointers dendat to tree 

joinnum<-length(hc$height)
lkm<-matrix(0,joinnum,1)
vec1<-matrix(0,1,d)
vec2<-matrix(0,1,d)
v1<-0
v2<-0
vapaa<-1
for (i in 1:joinnum){
   node1<-hc$merge[i,1]
   node2<-hc$merge[i,2]
   if (node1<0){
        parent[vapaa]<-i+n
        for (j in 1:d) vec1[j]<-dendat[-node1,j]
        v1<-1
        center[,vapaa]<-vec1
        volume[vapaa]<-v1
        pointers[-node1]<-vapaa
        lkm1<-1
        vapaa<-vapaa+1
   }
   else{
        parent[node1+n]<-i+n
        vec1<-center[,node1+n]
        v1<-volume[node1+n]
        lkm1<-lkm[node1]
   }
   if (node2<0){
        parent[vapaa]<-i+n
        for (j in 1:d) vec2[j]<-dendat[-node2,j]
        v2<-1
        center[,vapaa]<-vec2
        volume[vapaa]<-v2
        pointers[-node2]<-vapaa
        lkm2<-1
        vapaa<-vapaa+1
   }
   else{
        parent[node2+n]<-i+n
        vec2<-center[,node2+n]
        v2<-volume[node2+n]
        lkm2<-lkm[node2]
   }
   volume[i+n]<-1.1*(v1+v2)
   center[,i+n]<-(lkm1*vec1+lkm2*vec2)/(lkm1+lkm2)
   lkm[i]<-lkm1+lkm2
}

apoin<-matrix(0,n,1)
for (i in 1:n) apoin[i]<-nodenum-pointers[i]+1
apar<-parent[nodenum:1]
apar2<-matrix(0,n,1)
for (i in 1:nodenum) if (apar[i]!=0) apar2[i]<-nodenum-apar[i]+1

return(list(parent=apar2,level=level[nodenum:1],
volume=volume[nodenum:1],center=center[,nodenum:1],pointers=apoin))
}


dendat2lst<-function(dendat,lst,pcf)
{
# compare liketree

rnum<-length(pcf$value)
nodefinder<-matrix(0,rnum,1)
for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i

n<-dim(dendat)[1]
d<-dim(dendat)[2]
den2lst<-matrix(0,n,1)

step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

# find links from dendat to pcf
den2pcf<-matrix(0,n,1)
pcf2den<-matrix(0,rnum,1)
for (i in 1:n){
    j<-1
    while (j<=rnum){
         inside<-TRUE
         coordi<-1
         while ((inside) && (coordi<=d)){
             ala<-pcf$down[j,coordi]
             yla<-pcf$high[j,coordi]
             ala<-pcf$support[2*coordi-1]+ala*step[coordi]
             yla<-pcf$support[2*coordi-1]+yla*step[coordi]
             if ((dendat[i,coordi]<ala) || (dendat[i,coordi]>yla)) 
                         inside<-FALSE
             coordi<-coordi+1
         }
         if (inside){
                  den2pcf[i]<-j
                  pcf2den[j]<-i
         }
         j<-j+1
    }
}


for (i in 1:n){
    pcfi<-den2pcf[i]
    den2lst[i]<-nodefinder[pcfi]
}

return(den2lst)
}

depth2com<-function(dep,N){
#
d<-length(N)
logn<-log(N,base=2)
cusu<-cumsum(logn)
ind<-1
while ((ind<=d) && ((dep-cusu[ind])>0)){
  ind<-ind+1
}
direc<-min(ind,d)
if (direc==1){
  depind<-dep
}
else{
  depind<-dep-cusu[direc-1]
}
return(list(direc=direc,depind=depind))
}
depth<-function(mt){
#finds the dephts of the nodes
#
#mt is a result from multitree or pruneprof
#
roots<-mt$roots
child<-mt$child
sibling<-mt$sibling
#
itemnum<-length(child)
depths<-matrix(0,itemnum,1)
#
rootnum<-length(roots)
#
for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    pino[1]<-roots[i]  
    pinin<-1
    depths[roots[i]]<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        if (sibling[cur]>0){    #put right to stack
             pinin<-pinin+1
             pino[pinin]<-sibling[cur]
             depths[sibling[cur]]<-depths[cur]
        }
        while (child[cur]>0){    #go to leaf and put right nodes to stack
             chi<-child[cur]
             depths[chi]<-depths[cur]+1
             if (sibling[chi]>0){ 
                   pinin<-pinin+1
                   pino[pinin]<-sibling[chi]
                   depths[sibling[chi]]<-depths[cur]+1
             }
             cur<-chi
        }
    }
}
return(depths)
}







digit<-function(luku,base){
#Gives representation of luku for system with base
#
#luku is a natural number >=0
#base is d-vector of integers >=2, d>=2, 
#base[d] tarvitaan vain tarkistamaan onko luku rajoissa
#
#Returns d-vector of integers.
#
#example: digit(52,c(10,10)), returns vector (2,5)
#
d<-length(base)
digi<-matrix(0,d,1)
jako<-matrix(0,d,1)
jako[d]<-base[1]
for (i in (d-1):1){
  jako[i]<-base[d-i+1]*jako[i+1]
}
vah<-0
for (i in 1:(d-1)){
  digi[i]<-floor((luku-vah)/jako[i+1]) #if digi[i]>base[i], then ERROR
  vah<-vah+digi[i]*jako[i+1]
}
digi[d]<-luku-vah
# annetaan vastaus kaanteisesti se 2354 annetaan c(4,5,3,2)
# talloin vastaavuus sailyy base:n kanssa 
#apu<-matrix(0,d,1)
#for (i in 1:d){
#  apu[i]<-digi[d-i+1]
#}
apu<-digi[d:1]
return(apu)
}
dist.func<-function(dendat,xepsi=0,yepsi=0,col="black",type="distr",
log="y",cex.axis=1,dendat2=NULL,dendat3=NULL,col2="red",col3="blue",
pch2=20,pch3=20,split=median(dendat),xlim=NULL,xaxt="s",yaxt="s")
{
n<-length(dendat)

if (type=="distr"){

plot(x="",y="",
xlim=c(min(dendat)-xepsi,max(dendat)+xepsi), ylim=c(0-yepsi,1+yepsi),
xlab="",ylab="",cex.axis=cex.axis)
ycur<-0
ordi<-order(dendat)
dendatord<-dendat[ordi]
for(i in 1:(n-1)){
    segments(dendatord[i],ycur,dendatord[i],ycur+1/n,col=col)
    segments(dendatord[i],ycur+1/n,dendatord[i+1],ycur+1/n,col=col)
    ycur<-ycur+1/n
}
segments(dendatord[n],ycur,dendatord[n],ycur+1/n,col=col)
segments(dendatord[n],ycur+1/n,max(dendat)+xepsi,ycur+1/n,col=col)

}
else if ((type=="right.tail") || (type=="left.tail")){

  if (type=="right.tail"){
         redu.ind<-(dendat>split) 
         dendat.redu<-dendat[redu.ind]
  }
  else{
         redu.ind<-(dendat<split)
         dendat.redu<--dendat[redu.ind]
  }
  ordi<-order(dendat.redu)
  dendat.ord<-dendat.redu[ordi]
  nredu<-length(dendat.redu)
  level<-seq(nredu,1)
  if (type=="right.tail")
  plot(dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,xlim=xlim)
  else
  plot(-dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,xlim=xlim,
  xaxt=xaxt,yaxt=yaxt)

  #ordi<-order(dendat)
  #dendat.ord<-dendat[ordi]
  #medi.ind<-floor(n/2)
  #dendat.redu<-dendat.ord[medi.ind:n]
  #nredu<-length(dendat.redu)
  #level<-seq(nredu,1)
  #plot(dendat.redu,level,log="y",xlab="",ylab="")

  if (!is.null(dendat2)){

     if (type=="right.tail"){
          redu.ind<-(dendat2>split) 
          dendat.redu<-dendat2[redu.ind]
     }
     else{
          redu.ind<-(dendat2<split)
          dendat.redu<--dendat2[redu.ind]
     }
     ordi<-order(dendat.redu)
     dendat.ord<-dendat.redu[ordi]
     nredu<-length(dendat.redu)
     level<-seq(nredu,1)
     if (type=="right.tail")
     matplot(dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,
     add=TRUE,col=col2,pch=pch2)
     else
     matplot(-dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,
     add=TRUE,col=col2,pch=pch2)

  }

  if (!is.null(dendat3)){

     if (type=="right.tail"){
          redu.ind<-(dendat3>split) 
          dendat.redu<-dendat3[redu.ind]
     }
     else{
          redu.ind<-(dendat3<split)
          dendat.redu<--dendat3[redu.ind]
     }
     ordi<-order(dendat.redu)
     dendat.ord<-dendat.redu[ordi]
     nredu<-length(dendat.redu)
     level<-seq(nredu,1)
      if (type=="right.tail")
     matplot(dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,
     add=TRUE,col=col3,pch=pch3)
     else
     matplot(-dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,
     add=TRUE,col=col3,pch=pch3)
  }

}

}

dotouchgen<-function(indelow1,indeupp1,indelow2,indeupp2,direction){
#
epsi<-0
d<-length(indelow1)
touch<-TRUE
i<-1
while (i<=d){
     if ((i != direction) &&
         ((indelow1[i]>indeupp2[i]+epsi) || (indeupp1[i]<indelow2[i]-epsi))){
              touch<-FALSE
     }
     i<-i+1
}
return(touch)
}
dotouch<-function(inde1,inde2,direction){
#
d<-length(inde1)
touch<-TRUE
for (i in direction:d){
     if ((inde1[i]>inde2[i]+1) || (inde1[i]<inde2[i]-1)){
           touch<-FALSE
     }
}
return(touch)
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










drawhist<-function(dendat,binlkm,epsi=0,plkm){
#piirtaa 2-ulotteisessa tapauksessa histogramma estimaattorin kuvaajan 
#
#plkm on kuvaajan hilan pisteiden lkm
#
#dendat<-matrix(rnorm(20),10) 
#koe<-drawhist(dendat,binlkm=3,plk=30)
#persp(koe$x,koe$y,koe$z,phi=30,theta=60)
#
hi<-histo(dendat,binlkm,epsi)
recs<-hi$recs
values<-hi$values   
#
ans<-drawgene(values,recs,plkm)
return(list(x=ans$x,y=ans$y,z=ans$z))
}                                      








draw.kern<-function(value,index,N,support,minval=0,dendat=NULL,h=NULL)
{

d<-length(N)

if (d==2){

x<-matrix(0,N[1]+2,1)
y<-matrix(0,N[2]+2,1)
z<-matrix(minval,N[1]+2,N[2]+2)
#col<-matrix("black",dim(z)[1]*dim(z)[2],1)

minim<-matrix(0,d,1)  #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-support[2*i-1]
  maxim[i]<-support[2*i]
}
delta<-(maxim-minim)/(N+1)

indenum<-dim(index)[1]

i<-1
while (i<=indenum){
   inde<-index[i,]
   z[1+inde[1],1+inde[2]]<-value[i]
   #col[1+inde[1]+dim(z)[1]*inde[2]]<-ts[i]
   i<-i+1
}

i<-1
while (i<=N[1]){
   x[1+i]<-support[1]+delta[1]*i
   i<-i+1
}

i<-1
while (i<=N[2]){
   y[1+i]<-support[3]+delta[2]*i
   i<-i+1
}

x[1]<-support[1]
x[N[1]+2]<-support[2]
y[1]<-support[3]
y[N[2]+2]<-support[4]

return(list(x=x,y=y,z=z)) #col=col[length(col):1]))

}

else{    #d=1

x<-matrix(0,N+2,1)
y<-matrix(0,N+2,1)

minim<-min(dendat)
maxim<-max(dendat)
hmax<-max(h)
delta<-(maxim-minim+2*hmax)/(N+1)

indenum<-dim(index)[1]

i<-1
while (i<=indenum){

   inde<-index[i]
   point<-minim-hmax+delta*inde
    
   y[1+inde]<-value[i]
   x[1+inde]<-point

   i<-i+1
}
x[1]<-minim-hmax
x[N+2]<-minim-hmax+delta*N+delta

return(list(x=x,y=y))

}

}











draw.levset<-function(pcf,lev=NULL,bary=NULL,propor=0.1,col=NULL,
bound=NULL,dendat=NULL,xaxt="s",yaxt="s",cex.axis=1)
{

if (is.null(lev)) lev<-propor*max(pcf$value)

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]=(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i];

if (is.null(bound)){
  xmin<-pcf$support[1]
  xmax<-pcf$support[2]
  ymin<-pcf$support[3]
  ymax<-pcf$support[4]
}
else{
  xmin<-bound[1]
  xmax<-bound[2]
  ymin<-bound[3]
  ymax<-bound[4]
}

if (is.null(bary))
   plot(xmin,ymin,type="n",xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
   pch=20,col="red",xaxt=xaxt,yaxt=yaxt,cex.axis=cex.axis)
else
  plot(x=bary[1],y=bary[2],
  xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
  pch=20,col="red",xaxt=xaxt,yaxt=yaxt,cex.axis=cex.axis)

lenni<-length(pcf$value)
for (i in 1:lenni){
  if (pcf$value[i]>=lev){

     x1<-pcf$support[1]+step[1]*pcf$down[i,1]
     x2<-pcf$support[1]+step[1]*pcf$high[i,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[i,2]
     y2<-pcf$support[3]+step[2]*pcf$high[i,2] 

     if (is.null(col)) polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2))  
     else  polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=col[i],lty="blank")
  }
  i<-i+1
}

if (!is.null(dendat)) points(dendat)

#points(x=bary[1],y=bary[2],pch=20,col="red")

}



draw.pcf<-function(pcf,pnum=rep(32,length(pcf$N)),corona=5,dens=FALSE,minval=0,
drawkern=TRUE)
{
#Makes data for drawing a perspective plot.
#pnum on kuvaajan hilan pisteiden lkm
#corona makes corona around the support (useful for densities)

d<-length(pcf$N)

if (d==2){

#col=matrix("black",length(pcf$value),1)

step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
}

if ((drawkern)&&(!is.null(pcf$index))){
   return(draw.kern(pcf$value,pcf$index,pcf$N,pcf$support,minval=minval))
}

else{
     pit<-matrix(0,d,1)
     for (i in 1:d){
         pit[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pnum[i]
     }
     alkux<-pcf$support[1]+pit[1]/2-corona*pit[1]
     alkuy<-pcf$support[3]+pit[2]/2-corona*pit[2]
     loppux<-pcf$support[2]-pit[1]/2+corona*pit[1]
     loppuy<-pcf$support[4]-pit[2]/2+corona*pit[2]

     pnum2<-pnum+2*corona
     x<-alkux+c(0:(pnum2[1]-1))*pit[1]
     y<-alkuy+c(0:(pnum2[2]-1))*pit[2]

     reclkm<-length(pcf$value)
     xdim<-length(x)
     ydim<-length(y)
     arvot<-matrix(minval,xdim,ydim)

     l<-1
     while (l<=reclkm){
         begx<-pcf$support[1]+step[1]*(pcf$down[l,1])   #recs[l,1]
         endx<-pcf$support[1]+step[1]*(pcf$high[l,1])     #recs[l,2]
         begy<-pcf$support[3]+step[2]*(pcf$down[l,2])   #recs[l,3]
         endy<-pcf$support[3]+step[2]*(pcf$high[l,2])     #recs[l,4]

         begxind<-round(pnum2[1]*(begx-alkux)/(loppux-alkux))
         endxind<-round(pnum2[1]*(endx-alkux)/(loppux-alkux))
         begyind<-round(pnum2[2]*(begy-alkuy)/(loppuy-alkuy))
         endyind<-round(pnum2[2]*(endy-alkuy)/(loppuy-alkuy))

         arvot[begxind:endxind,begyind:endyind]<-pcf$value[l]
         #col[(begxind+(ydim-1)*begxind):(endxind+(ydim-1)*endyind)]<-ts[l]

         l<-l+1
     }

     #return(apu)
     return(list(x=x,y=y,z=arvot))
}  # else

}  # if (d==2)

else{  #d==1

if (dens){ 
  N<-pcf$N+2 
  alku<-1
}
else{
  N<-pcf$N
  alku<-0
}

x<-matrix(0,N,1)
y<-matrix(0,N,1)
step<-(pcf$support[2]-pcf$support[1])/pcf$N
minim<-pcf$support[1]
maxim<-pcf$support[2]
indenum<-length(pcf$value)

i<-1
while (i<=indenum){

   inde<-pcf$high[i]
   point<-minim+step*inde-step/2
    
   y[alku+inde]<-pcf$value[i]
   x[alku+inde]<-point

   i<-i+1
}

if (dens){
  x[1]<-minim-step/2
  x[N]<-maxim+step/2
}

# remove zeros
if (dens){
  numposi<-0
  for (i in 1:length(y)){
    if (y[i]>0){
       numposi<-numposi+1
       y[numposi]<-y[i]
       x[numposi]<-x[i]
    }
  }
  x<-x[1:numposi]
  y<-y[1:numposi]
}

or<-order(x)
xor<-x[or]
yor<-y[or]

return(list(x=xor,y=yor))

}  # else d=1

}











epane<-function(x,h){
#
d<-length(x)
val<-1
for (i in 1:d){
   val<-val*3*(1-(x[i]/h)^2)/2
}
return(val)
}
epan<-function(x)
{
d<-length(x)
val<-1
for (i in 1:d){
   val<-val*3*(1-x[i]^2)/2
}
return(val)
}
etais<-function(x,y){
#laskee euklid etais nelion vektorien x ja y valilla
#
pit<-length(x)
vast<-0
i<-1
while (i<=pit){
  vast<-vast+(x[i]-y[i])^2
  i<-i+1
}
return(vast)
}
etaisrec<-function(point,rec)
{
# calculates the squared diatance of a point to a rectangle

d<-length(rec)/2

res<-0
for (i in 1:d){
   if (point[i]>rec[2*i]) res<-res+(point[i]-rec[2*i])^2
   else if (point[i]<rec[2*i-1]) res<-res+(point[i]-rec[2*i-1])^2
}

return(res)

}



eva.clayton<-function(x,t,marginal="unif",sig=c(1,1),df=1)
# t>0
{
u<-x[1]
v<-x[2]
marg1<-1
marg2<-1

if (marginal=="normal"){
   u<-pnorm(x[1]/sig[1])
   v<-pnorm(x[2]/sig[2])
   marg1<-evanor(x[1]/sig[1])/sig[1]
   marg2<-evanor(x[2]/sig[2])/sig[2]
}
if (marginal=="student"){
   u<-pt(x[1],df)
   v<-pt(x[2],df)
   marg1<-dt(x[1],df)
   marg2<-dt(x[2],df)
}

val<-(1+t)*(u*v)^(-1-t)*(u^(-t)+v^(-t)-1)^(-2-1/t)*marg1*marg2

#if (val<0) val<-0

return(val)
}

eva.cop6<-function(x,t,marginal="unif",sig=c(1,1))
# t in [1,\infty)
{
u<-x[1]
v<-x[2]
marg1<-1
marg2<-1

if (marginal=="normal"){
   u<-pnorm(x[1]/sig[1])
   v<-pnorm(x[2]/sig[2])
   marg1<-evanor(x[1]/sig[2])/sig[1]
   marg2<-evanor(x[2]/sig[2])/sig[2]
}

val<-t*(1-u)^(t-1)*(1-v)^(t-1)*
((1-u)^t+(1-v)^t-(1-u)^t*(1-v)^t)^(1/t-2)*
(-(1/t-1)*(1-(1-u)^t)*(1-(1-v)^t)+
 (1-u)^t+(1-v)^t-(1-u)^t*(1-v)^t)*marg1*marg2

#if (val<0) val<-0

return(val)
}

eva.copula<-function(x,type="gauss",marginal="unif",sig=rep(1,length(x)),r=0,
t=rep(4,length(x)),g=1)
{
# sig is std of marginals, r is the correlation coeff, 
# t is deg of freedom

d<-length(x)
marg<-matrix(0,d,1)
u<-matrix(0,d,1)

if (marginal=="unif"){
   for (i in 1:d){
      u[i]<-x[i]/sig[i]  #+1/2
      marg[i]<-1/sig[i]
   }
}
if ((marginal=="normal")||(marginal=="gauss")){
   for (i in 1:d){
      u[i]<-pnorm(x[i]/sig[i])
      marg[i]<-evanor(x[i]/sig[i])/sig[i]
   }
}
if (marginal=="student"){
   for (i in 1:d){
      u[i]<-pt(x[i]/sig[i],df=t[i])
      marg[i]<-dt(x[i]/sig[i],df=t[i])/sig[i]
   }
}
if (type=="gauss"){
   d<-2
   x1<-qnorm(u[1],sd=1)
   x2<-qnorm(u[2],sd=1)

#   produ<-dnorm(x1,sd=1)*dnorm(x2,sd=1)
#   nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
#   vakio<-(2*pi)^(-d/2) 
#   g<-vakio*(1-r^2)^(-1/2)*exp(-(1/2)*nelio)
#   val<-g/produ*marg[1]*marg[2]

  copuval<-(1-r^2)^(-1/2)*
  exp(-(x1^2+x2^2-2*r*x1*x2)/(2*(1-r^2)))/exp(-(x1^2+x2^2)/2)
  val<-copuval*marg[1]*marg[2]

}
if (type=="gumbel"){
  link<-function(y,g){ return ( (-log(y))^g ) }
  linkinv<-function(y,g){ return ( exp(-y^(1/g)) ) }
  der1<-function(y,g){ return ( -g*(-log(y))^(g-1)/y ) }
  der2<-function(y,g){ return ( g*y^(-2)*(-log(y))^(g-2)*(g-1-log(y)) ) }
  
  linky<-link(u,g)
  a<-sum(linky)
  b<-linkinv(a,g)
  der1b<-der1(b,g)
  der2b<-der2(b,g)
  psi<--der2b*der1b^(-3)
  deriy<-der1(u,g) 
  val<-psi*abs(prod(deriy))*prod(marg)
}
if (type=="frank"){
  link<-function(y,g){ return ( -log((exp(-g*y)-1)/(exp(-g)-1)) ) }
  linkinv<-function(y,g){ return ( -log(1+(exp(-g)-1)/exp(y))/g  ) }
  der1<-function(y,g){ return ( g*exp(-g*y)/(exp(-g*y)-1) ) }
  der2<-function(y,g){ return ( g^2*exp(-g*y)/(exp(-g*y)-1)^2 ) }
  
  linky<-link(u,g)
  a<-sum(linky)
  b<-linkinv(a,g)
  der1b<-der1(b,g)
  der2b<-der2(b,g)
  psi<--der2b*der1b^(-3)
  deriy<-der1(u,g) 
  val<-psi*abs(prod(deriy))*prod(marg)
}
if (type=="clayton"){
  link<-function(y,g){ return ( y^(-g)-1 ) }
  linkinv<-function(y,g){ return ( (y+1)^(-1/g) ) }
  der1<-function(y,g){ return ( -g*y^(-g-1) ) }
  der2<-function(y,g){ return ( g*(g+1)*y^(-g-2) ) }
  
  linky<-link(u,g)
  a<-sum(linky)
  b<-linkinv(a,g)
  der1b<-der1(b,g)
  der2b<-der2(b,g)
  psi<--der2b*der1b^(-3)
  deriy<-der1(u,g) 
  val<-psi*abs(prod(deriy))*prod(marg)
}

return(val)
}






eva.gauss<-function(x,t=1,marginal="unif",sig=c(1,1),r=0,tapa1=TRUE)
{
#  sig is std of marginals

if (marginal=="unif"){
   u<-x[1]/sig[1]+1/2
   v<-x[2]/sig[2]+1/2
   marg1<-1/sig[1]
   marg2<-1/sig[2]
}
if (marginal=="normal"){
   u<-pnorm(x[1]/sig[1])
   v<-pnorm(x[2]/sig[2])
   marg1<-evanor(x[1]/sig[1])/sig[1]
   marg2<-evanor(x[2]/sig[2])/sig[2]
}
if (marginal=="student"){
   u<-pt(x[1]/sig[1],df=t)
   v<-pt(x[2]/sig[2],df=t)
   marg1<-dt(x[1]/sig[1],df=t)/sig[1]
   marg2<-dt(x[2]/sig[2],df=t)/sig[2]
}

d<-2
x1<-qnorm(u,sd=1)
x2<-qnorm(v,sd=1)

if (tapa1){
  produ<-dnorm(x1,sd=1)*dnorm(x2,sd=1)
  nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
  vakio<-(2*pi)^(-d/2) 
  g<-vakio*(1-r^2)^(-1/2)*exp(-(1/2)*nelio)
  val<-g/produ*marg1*marg2
}
else{
  # x1,x2 -> copuval
  copuval<-(1-r^2)^(-1/2)*
  exp(-(x1^2+x2^2-2*r*x1*x2)/(2*(1-r^2)))/exp(-(x1^2+x2^2)/2)
  val<-copuval*marg1*marg2
}

return(val)
}

eva.hat<-function(x,a=0.5,b=0.5)
{
# 0<a<1, b<1
# if b<a then marginal is unimodal
# if a^2 < b < a then not star unimodal

d<-length(x)
eta<-sum(x^2)     #vektorin x pituuden nelio
normvakio<-((2*pi)^d*(a^(-d)-b))^(-1)
tulos<-normvakio*(exp(-a^2*eta/2)-b*exp(-eta/2))

return(tulos)
}
eval.func.1D<-function(func,N,support=NULL,g=1,std=1,distr=FALSE,
M=NULL,sig=NULL,p=NULL,a=0.5,b=0.5,d=2)
{
if (func=="gauss"){
   norma<-(2*pi)^(-1/2)
   funni<-function(t){ fu<-exp(-t^2/2); return( norma*fu ) }
}
if (func=="polynomial"){
   support<-c(-std,std)
   norma<-(2*(1-1/(g+1)))^(-1)
   funni<-function(t){ fu<-1-abs(t)^g; return( norma*fu ) }
}
if (func=="student"){
   norma<-gamma((g+1)/2)/((g*pi)^(1/2)*gamma(g/2))
   funni<-function(t){ fu<-(1+t^2/g)^(-(g+1)/2); return( norma*fu ) }
   #y<-dt(x,df=g)
}
if (func=="exponential"){
   norma<-1/2
   funni<-function(t){ fu<-exp(-abs(t)); return( norma*fu ) }
}
if (func=="exponential"){
   norma<-1/2
   funni<-function(t){ fu<-exp(-abs(t)); return( norma*fu ) }
}
if (func=="mixt"){
   funni<-function(t){ 
       mixnum<-length(p)
       val<-0
       for (mi in 1:mixnum){
               evapoint<-(t-M[mi])/sig[mi]
               val<-val+p[mi]*evanor(evapoint)/sig[mi]
        }
        return( val ) 
   }
}
if (func=="hat"){
   normavak<-((2*pi)^d*(a^(-d)-b))^(-1)
   norma<-normavak*(2*pi)^((d-1)/2)
   funni<-function(t){  #(t,a,b,d,...){ 
          fu<-a^(1-d)*exp(-a^2*t^2)-b*exp(-t^2/2); return( norma*fu ) }
}

if (is.null(support)) support<-c(-1,1)

value<-matrix(0,N,1)
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

if (!distr){
   for (i in 1:N){
       inde<-i
       t<-lowsuppo+step*inde-step/2
       value[i]<-funni(t/std)/std
   }
}
else{
   inde<-1
   t<-lowsuppo+step*inde-step/2
   value[1]<-step*funni(t/std)/std
   for (i in 2:N){
       inde<-i
       t<-lowsuppo+step*inde-step/2
       value[i]<-value[i-1]+step*funni(t/std)/std
       #funni(t/std,g=g,a=a,b=b,d=d)/std
   }
}

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

res<-list(
value=value,
down=down,high=high,
#down=index-1,high=index,  
support=support,N=N)

return(res)
}

                              

eval.func.dD<-function(func,N,
sig=rep(1,length(N)),support=NULL,theta=NULL,g=1,
M=NULL,p=NULL,mul=3,
t=rep(1,length(N)),marginal="normal",r=0,
mu=NULL,xi=NULL,Omega=NULL,alpha=NULL,df=NULL,a=0.5,b=0.5
)   
# func== "mixt", "epan", "cop1"
{
d<-length(N)
recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (is.null(support)){

   if (func=="mixt"){
     support<-matrix(0,2*d,1)
        for (i in 1:d){
           support[2*i-1]<-min(M[,i]-mul*sig[,i])
           support[2*i]<-max(M[,i]+mul*sig[,i])
        }
   }

   if (func=="epan"){
      if (is.null(sig)) sig<-c(1,1)
      support<-matrix(0,2*d,1)
      for (i in 1:d){
          support[2*i-1]<--sig[i]
          support[2*i]<-sig[i]
      }
   }
}

if ((marginal=="unif")) support<-c(0,sig[1],0,sig[2])
# && (is.null(support))) 
#support<-c(-sig[1]/2,sig[1]/2,-sig[2]/2,sig[2]/2)


lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

numpositive<-0
for (i in 1:recnum){
    inde<-digit(i-1,N)+1
    #if ((inde[1]==0) && (inde[2]==N[2])) inde<-c(0,0)
    point<-lowsuppo+step*inde-step/2

    if (!is.null(theta)){
         rotmat<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
         point<-rotmat%*%point
    }

    if (func=="prod") val<-eva.prod(point,marginal,g)
    if (func=="skewgauss") val<-eva.skewgauss(point,mu,sig,alpha)
    #if (func=="dmsn") val<-dmsn(point,xi,Omega,alpha)
    if (func=="student") val<-eva.student(point,t,marginal,sig,r,df)
    if (func=="gumbel") val<-eva.copula(point,
        type="gumbel",marginal=marginal,sig=sig,r=r,t=t,g=g)
    if (func=="frank") val<-eva.copula(point,
        type="frank",marginal=marginal,sig=sig,t=t,g=g)
    if (func=="plackett") val<-eva.plackett(point,t,marginal,sig)
    if (func=="clayton2") val<-eva.clayton(point,t,marginal,sig,df)
    if (func=="clayton") val<-eva.copula(point,
        type="clayton",marginal=marginal,sig=sig,r=r,t=t,g=g)
    if (func=="cop6") val<-eva.cop6(point,t,marginal,sig)
    if (func=="epan") val<-epan(point)
    if (func=="gauss") val<-eva.copula(point,
        type="gauss",marginal=marginal,sig=sig,r=r,t=t)
    if (func=="normal") val<-eva.gauss(point,t=t,marginal=marginal,sig=sig,r=r)
    if (func=="mixt"){
        val<-0
        mixnum<-length(p)
        for (mi in 1:mixnum){
            evapoint<-(point-M[mi,])/sig[mi,]
            val<-val+p[mi]*evanor(evapoint)/prod(sig[mi,])
        }
    }
   if (func=="hat") val<-eva.hat(point,a=a,b=b)

    if (val>0){
       numpositive<-numpositive+1
       value[numpositive]<-val
       index[numpositive,]<-inde
    }
}

value<-value[1:numpositive]
index<-index[1:numpositive,]
down<-index-1
high<-index

res<-list(
value=value,index=index,
down=down,high=high,  #step=delta,
support=support,N=N)

return(res)
}                              











eva.lognormal<-function(t)
{

ans<-(2*pi)^(-1/2)*t^(-1)*exp(-(log(t))^2/2)

return(ans)
}

evanor<-function(x){

d<-length(x) 
eta<-sum(x^2)     #vektorin x pituuden nelio
normvakio<-(sqrt(2*pi))^{-d}
tulos<-exp(-eta/2)*normvakio
return(tulos)
}
eva.plackett<-function(x,t,marginal="unif",sig=c(1,1))
# t>=0, t \neq 1
{
u<-x[1]
v<-x[2]
marg1<-1
marg2<-1

if (marginal=="normal"){
   u<-pnorm(x[1]/sig[1])
   v<-pnorm(x[2]/sig[2])
   marg1<-evanor(x[1]/sig[1])/sig[1]
   marg2<-evanor(x[2]/sig[2])/sig[2]
}

val<-t*(1+(t-1)*(u+v-2*u*v))*((1+(t-1)*(u+v))^2-4*t*(t-1)*u*v)^(-3/2)*marg1*marg2
#if (val<0) val<-0

return(val)
}

eva.prod<-function(x,marginal="student",g=1)
{
if (marginal=="student"){
    d<-1
    vakio<-gamma((g+d)/2)/((g*pi)^(d/2)*gamma(g/2))
    y<-vakio*(1+x^2/g)^(-(g+d)/2)
    val<-prod(y)
}
if (marginal=="student.old"){
    d<-1
    vakio<-gamma((g+d)/2)/((g*pi)^(d/2)*gamma(g/2))
    x1<-vakio*(1+x[1]^2/g)^(-(g+d)/2)
    x2<-vakio*(1+x[2]^2/g)^(-(g+d)/2)
    val<-x1*x2
}
if (marginal=="studentR"){
    #x1<-dt(x[1],df=g)
    #x2<-dt(x[2],df=g)
    y<-dt(x,df=g)
    val<-prod(y)  
}
if (marginal=="polyno.old"){
    vakio<-2*(1-1/(g+1))
    y<-vakio*abs(1-x)^g
    val<-prod(y)
}
if (marginal=="polyno"){ 
    vakio<-1/(2*(1-1/(g+1)))
    y<-vakio*abs(1-abs(x)^g)
    val<-prod(y)
}
if (marginal=="double"){
    vakio<-1/2
    y<-exp(-abs(x))
    val<-prod(y)
}
if (marginal=="gauss"){
    vakio<-(2*pi)^(-1/2)
    y<-exp(-x^2/2)
    val<-prod(y)
}

  
return(val)
}

eva.skewgauss<-function(x,mu,sig,alpha)
{

norvak<-prod(sig)^(-1)
point<-(x-mu)/sig
en<-evanor(point)     #dnorm(poi)

point2<-alpha%*%((x-mu)/sig)
pn<-pnorm(point2)

ans<-2*en*pn

return(ans)
}
eva.student<-function(x,t=rep(4,length(x)),
marginal="unif",sig=c(1,1),r=0,df=1)
# t>2 
#  sig is std of marginals
{
if (marginal=="unif"){
   u<-x[1]/sig[1]
   v<-x[2]/sig[2]
   marg1<-1/sig[1]
   marg2<-1/sig[2]
}
if (marginal=="normal"){
   u<-pnorm(x[1]/sig[1])
   v<-pnorm(x[2]/sig[2])
   marg1<-evanor(x[1]/sig[1])/sig[1]
   marg2<-evanor(x[2]/sig[2])/sig[2]
}
if (marginal=="student"){
   u<-pt(x[1]/sig[1],df=t[1])
   v<-pt(x[2]/sig[2],df=t[2])
   marg1<-dt(x[1]/sig[1],df=t[1])/sig[1]
   marg2<-dt(x[2]/sig[2],df=t[2])/sig[2]
}

d<-2
x1<-qt(u,df=df)
x2<-qt(v,df=df)
produ<-dt(x1,df=df)*dt(x2,df=df)

nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
vakio<-gamma((df+d)/2)/((df*pi)^(d/2)*gamma(df/2))
ga<-vakio*(1-r^2)^(-1/2)*(1+nelio/df)^(-(df+d)/2)
#ga<-(1-r^2)^(1/2)*(1+(x1^2+x2^2-2*r*x1*x2)/(t*(1-r^2)))^(-(t+d)/2)

val<-ga/produ*marg1*marg2

return(val)
}

eva.t<-function(x,df,mu,Sigma)
{
norma<-gamma((df+1)/2)/((df*pi)^(1/2)*gamma(df/2))
fu<-(1+x^2/df)^(-(df+1)/2)
val<-norma*fu

return(val)
}


excmas.bylevel<-function(lst,levnum)
{
#source("~/denpro/R/excmas.bylevel.R")
#excmas.bylevel(lst,20)

levexc<-matrix(0,levnum,1)

maxlev<-max(lst$level)
step<-maxlev/levnum
nodelkm<-length(lst$parent)

mlkm<-moodilkm(lst$parent)
modloc<-mlkm$modloc    #pointers to modes
lkm<-mlkm$lkm       

added<-matrix(0,nodelkm,1)  #1 if we have visited this node

i<-1
while (i<=lkm){
    node<-modloc[i]
    # calculate curexc
    par<-lst$parent[node]
    if (par==0) valpar<-0 else valpar<-lst$level[par] 
    curexc<-(lst$level[node]-valpar)*lst$volume[node]
    
    nodelevind<-min(max(round(lst$level[node]/step),1),levnum)    
    levexc[1:nodelevind]<-levexc[1:nodelevind]+curexc

    while (lst$parent[node]>0){
         node<-lst$parent[node]
         if (added[node]==0){   
           # calculate curexc
           par<-lst$parent[node]
           if (par==0) valpar<-0 else valpar<-lst$level[par] 
           curexc<-(lst$level[node]-valpar)*lst$volume[node] 
           
           nodelevind<-min(max(round(lst$level[node]/step),1),levnum)    
           levexc[1:nodelevind]<-levexc[1:nodelevind]+curexc

           added[node]<-1 
         }
    }
    i<-i+1
}

levexc<-levexc/levexc[1]

diffe<-matrix(0,length(levexc),1)
for (i in 1:(length(levexc)-1)) diffe[i]<-(levexc[i+1]-levexc[i])/step
diffe[length(diffe)]<-diffe[length(diffe)-1]

return(list(levexc=levexc,diffe=diffe))
}




excmas<-function(lst){
#
parents<-lst$parent
volumes<-lst$volume
levels<-lst$level
#
nodelkm<-length(parents)
excmasses<-matrix(1,nodelkm,1)
#
mlkm<-moodilkm(parents)
modloc<-mlkm$modloc    #pointers to modes
lkm<-mlkm$lkm       
#
added<-matrix(0,nodelkm,1)  #1 if we have visited this node
#
for (i in 1:lkm){
    node<-modloc[i]
    # calculate curexc
    par<-parents[node]
    if (par==0) valpar<-0 else valpar<-levels[par] 
    curexc<-(levels[node]-valpar)*volumes[node]
    #
    excmasses[node]<-curexc
    while (parents[node]>0){
         node<-parents[node]
         if (added[node]==0){   
           # calculate curexc
           par<-parents[node]
           if (par==0) valpar<-0 else valpar<-levels[par] 
           curexc<-curexc+(levels[node]-valpar)*volumes[node] 
           #
           excmasses[node]<-curexc 
           added[node]<-1 
         }
         else{   #add only previous cumulative 
            excmasses[node]<-excmasses[node]+curexc
         }
    }
}
return(t(excmasses))
}




exmap<-function(estiseq,mt,hind=NULL,hseq=NULL,
n=NULL,moteslist=NULL,ylist=NULL)
{
#moteslist is a list of alpha values for every node
#not just for the branch nodes, but it might be nonsense for others

pk<-estiseq$lstseq
if (is.null(hseq)) hseq<-mt$hseq
if (is.null(hind)) hind<-c(1:length(mt$hseq))
slis<-mt$hseq[hind]

d<-dim(pk[[1]]$center)[1]

if (is.null(ylist)) ylist<-c(length(slis):1)

hrun<-1
for (i in 1:length(slis)){
   while (hseq[hrun]!=slis[i]){
      hrun<-hrun+1
   }
   if (i==1) treelist<-list(pk[[hrun]])
   else  treelist=c(treelist,list(pk[[hrun]]))
}

parnum<-length(slis)
veclkm<-0

if (d==1){
  crit<-max(treelist[[1]]$center)
}
else{
  crit<-rep(0,d)
}

for (i in 1:parnum){
     scur<-slis[i]

     if (!is.null(ylist))  yhigh<-ylist[i]
     else yhigh<-scur     

     profile<-treelist[[i]]
     
     if (!is.null(moteslist))  motes<-moteslist[[i]]
     else motes<-NULL

     level<-scur
     levelplot<-yhigh

     vecplu<-prof2vecs(profile,levelplot,n,crit,motes=motes)
     vecs<-vecplu$vecs
     depths<-vecplu$depths
     motes<-vecplu$motes
     mlabel<-vecplu$mlabel
     vecnum<-length(depths)
     smoot<-matrix(level,vecnum,1)
     
     # concatenate to big's
     
     veclkmold<-veclkm
     veclkm<-veclkm+vecnum
     if (veclkmold==0){   
        bigvecs<-vecs
        bigdepths<-depths
        bigmotes<-motes
        bigmlabel<-mlabel
        bigsmoot<-smoot
     }
     else{
       temp<-matrix(0,veclkm,4)
       temp[1:veclkmold,]<-bigvecs
       temp[(veclkmold+1):veclkm,]<-vecs
       bigvecs<-temp
       
       tempdep<-matrix(0,veclkm,1)
       tempdep[1:veclkmold]<-bigdepths
       tempdep[(veclkmold+1):veclkm]<-depths
       bigdepths<-tempdep
       
       tempmoo<-matrix(0,veclkm,1)
       tempmoo[1:veclkmold]<-bigmotes
       tempmoo[(veclkmold+1):veclkm]<-motes
       bigmotes<-tempmoo
       
       templab<-matrix(0,veclkm,1)
       templab[1:veclkmold]<-bigmlabel
       templab[(veclkmold+1):veclkm]<-mlabel
       bigmlabel<-templab
       
       tempsmoo<-matrix(0,veclkm,1)
       tempsmoo[1:veclkmold]<-bigsmoot
       tempsmoo[(veclkmold+1):veclkm]<-smoot
       bigsmoot<-tempsmoo
     }
}    
#if (makeplot==T) plotvecs(bigvecs,segme=T,bigdepths) 

return(list(bigvecs=bigvecs,bigdepths=t(bigdepths),motes=t(bigmotes),mlabel=t(bigmlabel),smoot=t(bigsmoot)))
}







explo.compa<-function(dendat,seed=1)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

cova<-cov(dendat)
mu<-mean(data.frame(dendat))

eig<-eigen(cov(dendat),symmetric=TRUE)
sigsqm<-eig$vectors%*%diag(eig$values^{1/2})  #%*%t(eig$vectors)

set.seed(seed)
symmedata<-matrix(rnorm(d*n),n,d)
dendat.simu<-t(sigsqm%*%t(symmedata))

return(dendat.simu)
}

findbnodes<-function(lst,modenum=1,num=NULL)
{
# prunes from a level set tree "lst" the modes with "num" 
# smallest excess masses 
# or the modes with smaller excess mass than "exmalim"

if (is.null(num)){
    curmodenum<-moodilkm(lst$parent)$lkm
    num<-curmodenum-modenum
}

len<-length(lst$parent)
child.frekve<-matrix(0,len,1)
for (i in 1:len){
    if (lst$parent[i]>0) 
    child.frekve[lst$parent[i]]<-child.frekve[lst$parent[i]]+1
}

ml<-moodilkm(lst$parent)
mode.list<-ml$modloc
roots.of.modes<-matrix(0,length(mode.list),1)
for (aa in 1:length(mode.list)){
    node<-mode.list[aa]
    while ((lst$parent[node]>0) && (child.frekve[lst$parent[node]]==1)){ 
         node<-lst$parent[node]
    }
    roots.of.modes[aa]<-node
}

em<-excmas(lst)
or<-order(em[roots.of.modes],decreasing=TRUE)
#nodes<-ml$modloc[or[1:modenum]]
nodes<-roots.of.modes[or[1:modenum]]

return(nodes=nodes)
}


findbranch.pare<-function(parent)
{
# finds the nodes who have more than 1 child

len<-length(parent)
frekve<-matrix(0,len,1)

for (i in 1:len){
   if (parent[i]>0) frekve[parent[i]]<-frekve[parent[i]]+1
}

tulos<-matrix(0,len,1)

for (i in 1:len){
     #if (parent[i]==0) tulos[i]<-1
     #else 
     if ((parent[i]!=0) && (frekve[parent[i]]>1)){ #result of a branching
                 tulos[parent[i]]<-1
    }
}

if (sum(tulos)==0) ans<-NULL else ans<-which(tulos==1)

return(ans)
}    
findbranch<-function(parent,colo="red1",pch=22)
{
# finds the nodes which make the tree of the branches

#pch=19: solid circle, pch=20: bullet (smaller circle), 
#pch=21: circle, pch=22: square, 
#pch=23: diamond, pch=24: triangle point-up, 
#pch=25: triangle point down. 

len<-length(parent)
frekve<-matrix(0,len,1)

for (i in 1:len){
   if (parent[i]>0) frekve[parent[i]]<-frekve[parent[i]]+1
}

tulos<-matrix(0,len,1)
colovec<-matrix("black",len,1)
pchvec<-matrix(21,len,1)

for (i in 1:len){
    if (parent[i]==0){ #root node
             tulos[i]<-1  
             colovec[i]<-colo
             pchvec[i]<-pch
 
    }  
    else if (frekve[parent[i]]>1){ #result of a branching
                 tulos[i]<-1
                 colovec[i]<-colo
                 pchvec[i]<-pch  
    }
}

return(list(indicator=tulos,colovec=colovec,pchvec=pchvec))
}    
findend<-function(inde,left,right,low,upp,N){
#
lenn<-length(inde)
current<-1
dep<-1
if ((left[current]==0) && (right[current]==0)){
   exists<-FALSE
   }
else{
   exists<-TRUE
}
while ((exists) && ((left[current]>0) || (right[current]>0))){
     mid<-(low[current]+upp[current])/2 
     direc<-depth2com(dep,N)$direc 
     if (inde[direc]<=mid){
           if (left[current]>0){ 
                current<-left[current]
                dep<-dep+1
           }
           else{ 
                exists<-FALSE
           }
     }
     else{   
           if (right[current]>0){
                  current<-right[current]
                  dep<-dep+1
           }
           else{
               exists<-FALSE
           }
     }
}
return(list(exists=exists,location=current,dep=dep))
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







findneighbor<-function(lst,node)
{
mu<-multitree(lst$parent)

no<-lst$parent[node]
while ((no!=0) && (mu$sibling[mu$child[no]]==0)){
     no<-lst$parent[no]
}


return(no)
}

fs.calc.parti<-function(pa,dendat,h)
{
#type =  "barys", "means", "mins", "maxs"
lkm<-dim(pa$recs)[1]
n<-dim(dendat)[1]
d<-dim(dendat)[2]
fs<-matrix(0,lkm,1)

for (i in 1:lkm){
    recu<-pa$recs[i,]
    arg<-matrix(0,d,1)
    for (j in 1:d){
        arg[j]<-(recu[2*j-1]+recu[2*j])/2
    }
    fs[i]<-kernesti.dens(arg,dendat,h)
}

return(fs)
}

fs.calc<-function(complex,dendat,h,type="barys")
{
#type =  "barys", "means", "mins", "maxs"
lkm<-dim(complex)[1]
n<-dim(dendat)[1]
d<-dim(dendat)[2]
fs<-matrix(0,lkm,1)

if (type!="barys"){
   f<-matrix(0,n,1)
   for (i in 1:n){
      arg<-dendat[i,]
      f[i]<-kernesti.dens(arg,dendat,h)
   }
   for (i in 1:lkm){
       vs<-complex[i,]
       vals<-f[vs]
       if (type=="means") fs[i]<-mean(vals) 
       if (type=="maxs") fs[i]<-max(vals) 
       if (type=="mins") fs[i]<-min(vals) 
   }
}

if (type=="barys"){
   #barys<-matrix(0,lkm,d)
   for (i in 1:lkm){
      simple<-complex[i,]
      simp<-dendat[simple,]
      arg<-colSums(simp)/(d+1)  #arg<-barys[i,]
      fs[i]<-kernesti.dens(arg,dendat,h)
   }
}

return(fs)
}

graph.matrix.level<-function(dendat,tt=NULL,permu=seq(1:dim(dendat)[1]),
col=seq(1:2000),config="new",shift=0.1,segme=TRUE,poin=FALSE,epsi=0,
ystart=0.5,pch=21,cex=1, yaxt="s",cex.axis=1,texto=TRUE)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

origins<-matrix(0,d,1)
starts<-matrix(0,d,1)
range<-matrix(0,d,1)
minis<-matrix(0,d,1)
means<-matrix(0,d,1)
for (i in 1:d){
    minis[i]<-min(dendat[,i]) 
    means[i]<-mean(dendat[,i])
}

for (i in 1:d) range[i]<-max(dendat[,i])-min(dendat[,i]) 
starts[1]<-0 #min(dendat[,1])
i<-2
while (i<=d){ 
      starts[i]<-starts[i-1]+range[i-1]+epsi
      i<-i+1
}
if (config=="new")
  for (i in 1:d) origins[i]<-starts[i]+mean(dendat[,i])-min(dendat[,i])
else
  for (i in 1:d) origins[i]<-starts[i]+range[i]/2-min(dendat[,i])

  #starts[i]+range[i]/2

plot(x="",y="",xlim=c(starts[1],starts[d]+range[d]),ylim=c(ystart,n+0.5),
xlab="",ylab="",xaxt="n",yaxt=yaxt,cex.axis=cex.axis)

if (is.null(tt)){

   for (j in 1:d){
      for (i in 1:n){
           indo<-permu[i]
           x0<-dendat[indo,j]+starts[j]-min(dendat[,j]) 
                    #dendat[indo,j]+origins[j]
           y0<-i
           x1<-origins[j]  
           y1<-i
           if (segme)
           segments(x0, y0, x1, y1, col=col[indo])
           if (poin)          
           points(x0, y0, col=col[indo], pch=pch, cex=cex)


       }
       if (j>1){
            beg<-starts[j]-epsi/2 
            segments(beg,0.5,beg,n+0.5,lty=1,lwd=2)
       }
       segments(origins[j],1,origins[j],n,lty=2)
       if (texto){
       text(origins[j],ystart,as.character(round(means[j],digits=2)),cex=cex)
       text(starts[j]+shift,ystart,as.character(round(minis[j],digits=2)),
       cex=cex)
       }
   }

}

else{
   paletti<-seq(1:2000)
   coli<-colobary(tt$parent,paletti)  #,segtype="num")

   cente<-c(mean(dendat[,1]),mean(dendat[,2]))
   dendat3<-matrix(0,n,d)
   newcolo<-matrix(0,n,1)
   maxseg<-max(coli)
   curbeg<-1
   i<-maxseg
   while (i >= 1){
       inditree<-which(coli==i)
       indidendat<-tt$infopointer[inditree]
       curend<-curbeg+length(indidendat)-1
       curseg<-dendat[indidendat,]
       leveli<-sqrt(sum(curseg-cente)^2)    #pituus(curseg-cente))
       or<-order(leveli)
       if (length(or)>1) orseg<-curseg[or,] else orseg<-curseg
       dendat3[curbeg:curend,]<-orseg
       newcolo[curbeg:curend]<-i
       curbeg<-curend+1  #curbeg+length(indit)+1
       i<-i-1
   }

   for (j in 1:d){
       for (i in 1:n){
           x0<-dendat3[i,j]+starts[j]-min(dendat3[,j]) 
               #dendat3[i,j]+origins[j]
           y0<-i
           x1<-origins[j]
           y1<-i
           #segments(x0, y0, x1, y1,col=newcolo[i])
           if (segme)
           segments(x0, y0, x1, y1, col=newcolo[i])
           if (poin)          
           points(x0, y0, col=newcolo[i],pch=pch,cex=cex)

       }
       if (j>1){
            beg<-starts[j]-epsi/2 
            segments(beg,0.5,beg,n+0.5,lty=1,lwd=2)
       }
       segments(origins[j],1,origins[j],n,lty=2)
       if (texto){
       text(origins[j],ystart,as.character(round(means[j],digits=2)),cex=cex)
       text(starts[j]+shift,ystart,as.character(round(minis[j],digits=2)),
       cex=cex)
       } 
   }

}  # else

}






graph.matrix<-function(dendat,type="level",
tt=NULL,permu=seq(1:dim(dendat)[1]),col=seq(1:2000),
config="new",shift=0.1,segme=TRUE,poin=FALSE,epsi=0,ystart=0.5, 
pch=21, cex=1, cex.axis=1, yaxt="s",
# profile:
ylen=100,profcol=rep("black",n),texto=TRUE)
{
if (type=="level") 
 graph.matrix.level(dendat, tt=tt, permu=permu, col=col,
 config=config, shift=shift, segme=segme, poin=poin, epsi=epsi, ystart=ystart,
 pch=pch, cex=cex, yaxt=yaxt, cex.axis=cex.axis, texto=texto)

else{  # type="profile"
 n<-dim(dendat)[1]
 d<-dim(dendat)[2]
 x<-seq(1:n)
 y<-seq(1:ylen)
 z<-matrix(0,length(x),length(y))
 varit<-matrix("",length(x),length(y))
 ala<-matrix(0,d,1)
 for (i in 1:d) ala[i]<-min(dendat[,i])
 yla<-matrix(0,d,1)
 for (i in 1:d) yla[i]<-max(dendat[,i])
 range<-yla-ala
 alaind<-matrix(0,d,1)
 alaind[1]<-1
 for (i in 2:d) 
     alaind[i]<-min(alaind[i-1]+round(ylen*range[i]/sum(range))+1,ylen)
 ylaind<-matrix(0,d,1)
 ylaind[d]<-ylen
 for (i in 1:(d-1)) ylaind[i]<-alaind[i+1]+1
 plot(x="",y="",xlim=c(0,n),ylim=c(0,ylen),xlab="",ylab="",yaxt="n",xaxt="n")
 for (i in 1:n){
     for (j in 1:d){
         suht<-(dendat[i,j]-ala[j])/range[j]
         korkeus<-round(suht*(ylaind[j]-alaind[j]))
         alku<-alaind[j]
         loppu<-min(max(alku+korkeus,1),ylen)
         if (alku>=loppu) loppu<-loppu+1
         polygon(x=c(i-1,i-1,i,i),y=c(alku,loppu,loppu,alku),col=profcol[i],
                 lty="blank")
         #z[i,alku:loppu]<-1
     }
 }
 #image(x,y,z,col=c("white","black"),xlab="",ylab="",xaxt="n",yaxt="n")
}

}




hgrid<-function(h1,h2,lkm,base=10)
{
step<-(h2-h1)/(lkm-1)

if (is.null(base)){
   hseq<-seq(h2,h1,-step)
}
else{
   a<-(h2-h1)/(base^(h2)-base^(h1))
   b<-h1-a*base^(h1)
   un<-seq(h2,h1,-step)
   hseq<-a*base^(un)+b
}

return(hseq)
}

histo1d<-function(dendat,binlkm,ala=NULL,yla=NULL,
pic=TRUE,brush=NULL,brushcol=c("blue"),col=NULL,border=NULL,
xlab="",ylab="",cex.lab=1,cex.axis=1,data=FALSE,
weights=rep(1,length(dendat)),normalization=TRUE,
height=NULL,subweights=NULL,graphplot=FALSE)
{
if (is.null(ala)) ala<-min(dendat)
if (is.null(yla)) yla<-max(dendat)
step<-(yla-ala)/binlkm
frekv<-matrix(0,binlkm,1)
value<-matrix(0,binlkm,1)
if (!is.null(brush)){
   cnum<-max(brush)
   shade <-matrix(0,binlkm,cnum)
}
if (!is.null(subweights)) taint<-matrix(0,binlkm,1)
n<-length(dendat)
for (i in 1:n){
   hava<-dendat[i]
   weight<-weights[i]
   ind<-min(binlkm,floor((hava-ala)/step)+1)
   frekv[ind]<-frekv[ind]+weight
   if ((!is.null(brush)) && (brush[i]>0)) 
              shade[ind,brush[i]]<-shade[ind,brush[i]]+1
   if (!is.null(subweights)) taint[ind]<-taint[ind]+n*subweights[i]
}
if (normalization) value<-frekv/(n*step) else value<-frekv
if (!is.null(brush)) shade<-shade/(n*step)
if ((normalization) && (!is.null(subweights))) taint<-taint/(n*step)

if (pic){
   if (is.null(height)) height<-max(value)
   plot(x="",y="",xlab=xlab,ylab=ylab,xlim=c(ala,yla),ylim=c(0,height),
   cex.lab=cex.lab,cex.axis=cex.axis)
   for (i in 1:binlkm){
          xala<-ala+(i-1)*step
          xyla<-xala+step
          y<-value[i]

          if (graphplot){
               if (i==1) yeka<-0 else yeka<-value[i-1]
               if (i==binlkm) ytok<-0 else ytok<-value[i]
               segments(xala,yeka,xala,ytok)
               segments(xala,ytok,xyla,ytok)
          } 
          else
          polygon(c(xala,xala,xyla,xyla),c(0,y,y,0),col=col,border=border)

          if (!is.null(brush)){
              y0<-0
              for (kk in 1:cnum){
                  y<-y0+shade[i,kk]
                  polygon(c(xala,xala,xyla,xyla),c(y0,y,y,y0),col=brushcol[kk])
                  y0<-y
              }
          }
          if (!is.null(subweights)){
              if (graphplot){
                 if (i==1) yeka<-0 else yeka<-taint[i-1]
                 if (i==binlkm) ytok<-0 else ytok<-taint[i]
                 segments(xala,yeka,xala,ytok,col=brushcol)
                 segments(xala,ytok,xyla,ytok,col=brushcol)
              } 
              else{
                 y<-taint[i]
                 polygon(c(xala,xala,xyla,xyla),c(0,y,y,0),col=brushcol)
              } 
          }
   }
}
if (data){
     return(list(frekv=frekv,ala=ala,step=step,value=value))
}
}



histo2data<-function(pcf){

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]=(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i];
xmin<-pcf$support[1]
xmax<-pcf$support[2]
ymin<-pcf$support[3]
ymax<-pcf$support[4]
zmin<-pcf$support[5]
zmax<-pcf$support[6]

nnew<-length(pcf$value)
desdat<-matrix(0,nnew,3)
for (i in 1:nnew){
     x1<-pcf$support[1]+step[1]*pcf$down[i,1]
     x2<-pcf$support[1]+step[1]*pcf$high[i,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[i,2]
     y2<-pcf$support[3]+step[2]*pcf$high[i,2] 
     z1<-pcf$support[5]+step[3]*pcf$down[i,3]
     z2<-pcf$support[5]+step[3]*pcf$high[i,3] 
     desdat[i,]<-c((x1+x2)/2,(y1+y2)/2,(z1+z2)/2)
}

f0<-sqrt(pcf$value)
colo<-1-(f0-min(f0)+0.5)/(max(f0)-min(f0)+0.5)
col<-gray(colo)

return(list(dendat=desdat,col=col))
}


histo<-function(dendat,binlkm,epsi=0)
{
# Constructs a histogram estimate: result is given by giving level
# sets of the estimate

supp<-support(dendat,epsi)
regdat<-den2reg(dendat,binlkm,supp)
palvak<-makehis(regdat)
values<-palvak$values
recs<-palvak$recs

integ<-0
recnum<-length(values)
for (i in 1:recnum){
   integ<-integ+values[i]*massone(recs[i,])
}

values<-values/integ

return(list(values=values,recs=recs))
}




intersec.edges<-function(edge1,edge2)
{
# returns 1 if there is intersection otherwise 0
# edge1, edge2 are d*d matrices
# rows are points in R^d
# edge is d points in R^d

d<-2

x1<-edge1[1,]
x2<-edge1[2,]
y1<-edge2[1,]
y2<-edge2[2,]

A<-matrix(0,d,d)
A[1,1]<-x1[1]-x2[1]
A[2,1]<-x1[2]-x2[2]
A[1,2]<--(y1[1]-y2[1])
A[2,2]<--(y1[2]-y2[2])
tulos<-0
if (det(A)!=0){
    invA<-solve(A,diag(rep(1,d)))
    vec<-matrix(y2-x2,2,1)
    tu<-invA%*%vec
    if ( (tu[1]>=0) && (tu[1]<=1) && (tu[2]>=0) && (tu[2]<=1) ) tulos<-1
}

return(tulos)
}


intersec<-function(taso,endind,cur,uni){
#Makes from a set of rectangles "cur" a larger rectangle.
#For a given rectangle in cur, we make intersection with
#rectangles in uni, starting with rectangle after the point
#indicated in "endind". 
#Result is a set of k over "taso" rectangles, where k is the number
#of rectangles in "uni". 
#uni has the basic sets, we have in cur all the (taso-1)-fold
#intersections of uni, we want to form taso-fold intersections,
#in endind we have index of the last rectangle in (taso-1)-fold
#intersection: we have to form intersections with all the rest
#rectangles in uni. Thus result has the size: how many subsets of
#size taso, we can take from a set of size k. 
#
#taso is integer >=1
#endind is l-vector
#cur is l*(2*d)-matrix
#uni is k*(2*d)-matrix, huom oletetaan etta k>1!!!!!
#
#Return NA if there is no intersectio, otherwise
#list(set=tulos,endind=newendind)
#
k<-length(uni[,1])   #rows of uni
d2<-length(uni[1,])  #col of uni is the 2*d
tulrow<-choose(k,taso)
#tulrow<-gamma(k+1)/(gamma(k-taso+1)*gamma(taso+1))#k yli taso,gamma(k)=(k-1)!
newendind<-matrix(0,tulrow,1)
tulos<-matrix(0,tulrow,d2)
ind<-0        #indeksi to tulos and newendind
if (dim(t(cur))[1]==1) a<-1 else a<-length(cur[,1])  #rows of cur
for (i in 1:a){
  if (endind[i]<k){
    for (j in (endind[i]+1):k){
      if (a==1) apu<-leikkaa(cur,uni[j,])
        else apu<-leikkaa(cur[i,],uni[j,])
         #for (l in 1:d){
         #  tulos[ind,2*l-1]<-max(cur[i,2*l-1],uni[j,2*l-1])
         #  tulos[ind,2*l]<-min(cur[i,2*l],uni[j,2*l])
         #}      
      if (!is.na(apu)){  #if there is intersection, save the result
         ind<-ind+1
         tulos[ind,]<-apu
         newendind[ind]<-j
      }
    }
  }
}
if (ind==0) palauta<-NA
else{
  tulos<-tulos[1:ind,]
  newendind<-newendind[1:ind]
  palauta<-list(set=tulos,endind=newendind)
}
return(palauta) 
}







intersec.simpces2<-function(simp1,simp2)
{
# returns 1 if there is intersection otherwise 0
# simp1, simp2 are (d+1)*d matrices

d<-2
tulos<-0
for (i in 1:d){
for (j in (i+1):(d+1)){
       x1<-simp1[i,]
       x2<-simp1[j,]
       for (ii in 1:d){
       for (jj in (ii+1):(d+1)){
           y1<-simp2[ii,]
           y2<-simp2[jj,]

A<-matrix(0,d,d)
A[1,1]<-x1[1]-x2[1]
A[2,1]<-x1[2]-x2[2]
A[1,2]<--(y1[1]-y2[1])
A[2,2]<--(y1[2]-y2[2])
tulos<-0
if (det(A)!=0){
    invA<-solve(A,diag(rep(1,d)))
    vec<-matrix(y2-x2,2,1)
    tu<-invA%*%vec
    if ( (tu[1]>=0) && (tu[1]<=1) && (tu[2]>=0) && (tu[2]<=1) ) tulos<-1
}


       }
       }  
}
}

return(tulos)
}


intersec.simpces<-function(simp1,simp2)
{
# returns 1 if there is intersection otherwise 0
# simp1, simp2 are (d+1)*d matrices

d<-2

tulos<-is.inside(simp1,simp2)
if (tulos==0) tulos<-is.inside(simp2,simp1)

if (tulos==0)
for (i in 1:d){
for (j in (i+1):(d+1)){
       x1<-simp1[i,]
       x2<-simp1[j,]
       for (ii in 1:d){
       for (jj in (ii+1):(d+1)){
           y1<-simp2[ii,]
           y2<-simp2[jj,]

           edge1<-matrix(0,2,2)
           edge2<-matrix(0,2,2)

           edge1[1,]<-x1
           edge1[2,]<-x2
           edge2[1,]<-y1
           edge2[2,]<-y2

           tulos<-intersec.edges(edge1,edge2)

       }
       }  
}
}


return(tulos)
}

is.inside<-function(simp1,simp2)
{
# simp1, simp2 (d+1)*d matrices
# returns 1 if simp1 is inside simp2

d<-2  #dim(simp1)[2]

deet<-matrix(0,3,1)
lk<-1
for (ii in 1:(d+1)){
    v1<-simp2[ii,]
    jj<-ii+1
    while (jj<=(d+1)){
         v2<-simp2[jj,]
         deet[lk]<-sqrt( sum((v1-v2)^2) )
         jj<-jj+1
         lk<-lk+1
    }
}
rho<-max(deet)

tulos<-1
i<-1
while ( (i<=(d+1)) && (tulos==1) ){
    vertice1<-simp1[i,]
    j<-1
    while ( (j<=(d+1)) && (tulos==1) ){
        vertice2<-simp2[j,]
        eta<-sqrt( sum((vertice1-vertice2)^2) )
        if (eta>rho) tulos<-0
        j<-j+1
    }
    i<-i+1
}

return(tulos)
}


is.inside.simp.bary<-function(point,simple)
{
# point is d-vector
# simple is (d+1)*d matrix of vertices
# return 1 if is inside
# use barycentric coordinates

d<-2
v1<-simple[1,]
v2<-simple[2,]
v3<-simple[3,]

x<-point[1]
y<-point[2]
x1<-v1[1]
y1<-v1[2]
x2<-v2[1]
y2<-v2[2]
x3<-v3[1]
y3<-v3[2]


l1<-((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3))
l2<-((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3))
l3<-1-l1-l2

if ((0<l1) && (l1<1) && (0<l2) && (l2<1) && (0<l3) && (l3<1)) tulos<-1
else tulos<-0

return(tulos)
}

is.inside.simp.long<-function(point,simple,rho)
{
# point is d-vector
# simple is (d+1)*d matrix of vertices

eps<-rho/100000
d<-2
tulos<-1
i<-1
while ( (i<=(d+1)) && (tulos==1) ){
      y1<-point
      y2<-simple[i,]
      if (y1[1]<y2[1]) y21<-y2[1]-eps else y21<-y2[1]+eps
      if (y1[2]<y2[2]) y22<-y2[2]-eps else y22<-y2[2]+eps
      for (jj in 1:d){
          for (kk in (jj+1):(d+1)){
             x1<-simple[jj,]
             x2<-simple[kk,]

             edge1<-matrix(0,2,2)
             edge2<-matrix(0,2,2)

             edge1[1,]<-x1
             edge1[2,]<-x2
             edge2[1,]<-y1
             edge2[2,]<-y2  #c(y21,y22)

             ints<-intersec.edges(edge1,edge2)
             if (ints==1) tulos<-0
          }
      }  
      i<-i+1
}

return(tulos)
}


is.inside.simp<-function(point,simple,rho)
{
d<-length(point)
dd<-1
sisalla2<-1
while ( (dd<=(d+1)) && (sisalla2==1) ){
      karki<-simple[dd,]
      dista2<-sum((point-karki)^2)
      if (dista2>rho^2){ 
            sisalla2<-0
      }
      dd<-dd+1
}

return(sisalla2)
}






joincongen<-function(leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,low,upp)
{
#upp and low are terminalnum*d matrices
#vrt index: we have links from the tree to this data structure
#we know that coordinate "direction" touches
#
#
# 1. We find the startpoints for the sets
#
# we will make three vectors of startpoints:
# a. startpoints for the whole sets (concatenate left and right child)
#    (startpointsS)
# b. startpoints for the right boundary of the left child and 
#    left boundary of the right child
#    (startpointsB)
# c.1. startpoints for the left boundary of the left child a
#    (startpointsNewBleft)
# c.2. startpoints for the right boundary of the right child
#    (startpointsNewBright)
#
#
# startpoints are pointers to 
# a. atomsSepaAtoms/atomsSepaNext
# b. atomsSepaAtoms/atomsRBounNext, atomsSepaAtoms/atomsLBounNext
# c.1. atomsSepaAtoms/atomsLBounNext
# c.2.  atomsSepaAtoms/atomsRBounNext
#
# startpoints are found from 
# a. begsSepaBegs
# b. begsLeftBoun, begsRighBoun
# c.1 begsLeftBoun, 
# c.2 begsRighBoun
#
# b. is used to check which touch
# a., c.1. and c.2. are joined together
# Note that some sets of the boundary are empty 
# (we store 0 to the respective location in begsLeftBoun, begsRighBoun )

suppi<-length(begsSepaNext)
startpointsS<-matrix(0,suppi,1)
startpointsB<-matrix(0,suppi,1)
startpointsNewBleft<-matrix(0,suppi,1)
startpointsNewBright<-matrix(0,suppi,1)
#new boundary: left bound. of left child, right b. of right child

induksi<-1
anfang<-separy[leftbeg]
startpointsS[induksi]<-anfang
startpointsB[induksi]<-begsRighBoun[anfang]
startpointsNewBleft[induksi]<-begsLeftBoun[anfang]
while (begsSepaNext[anfang]>0){
  anfang<-begsSepaNext[anfang]
  induksi<-induksi+1
  startpointsS[induksi]<-begsSepaBegs[anfang]
  startpointsB[induksi]<-begsRighBoun[anfang]
  startpointsNewBleft[induksi]<-begsLeftBoun[anfang]  
}
mleft<-induksi
induksi<-induksi+1
anfang<-separy[rightbeg]
startpointsS[induksi]<-anfang
startpointsB[induksi]<-begsLeftBoun[anfang]
startpointsNewBright[induksi]<-begsRighBoun[anfang]
while (begsSepaNext[anfang]>0){
  anfang<-begsSepaNext[anfang]
  induksi<-induksi+1
  startpointsS[induksi]<-begsSepaBegs[anfang]
  startpointsB[induksi]<-begsLeftBoun[anfang]
  startpointsNewBright[induksi]<-begsRighBoun[anfang]  

}
startpointsS<-startpointsS[1:induksi]
startpointsB<-startpointsB[1:induksi]
startpointsNewBleft<-startpointsNewBleft[1:induksi]
startpointsNewBright<-startpointsNewBright[1:induksi]
m<-induksi
mright<-m-mleft  


# 2. We make "links" matrix and apply declev

# We utilize previous programs

linkit<-matrix(0,m,m)
do<-1
while (do <= mleft){
   beg1<-startpointsB[do]    #could be 0
   re<-mleft+1
   while (re <= m){
       beg2<-startpointsB[re]    #could be 0
       conne<-FALSE
       begbeg1<-beg1
       while (begbeg1>0){
            begbeg2<-beg2
            while (begbeg2>0){
                atom1<-atomsSepaAtom[begbeg1]
                indelow1<-low[atom1,]
                indeupp1<-upp[atom1,]
                atom2<-atomsSepaAtom[begbeg2]
                indelow2<-low[atom2,]
                indeupp2<-upp[atom2,]
                if (dotouchgen(indelow1,indeupp1,indelow2,indeupp2,direction)){
                   conne<-TRUE
                }
                begbeg2<-atomsLBounNext[begbeg2]
            }
            begbeg1<-atomsRBounNext[begbeg1]
       }                
       if (conne){
           linkit[do,re]<-1
       }
       re<-re+1
   }
   do<-do+1
}
for (do in (mleft+1):m){
   beg1<-startpointsB[do]
   for (re in 1:mleft){
       beg2<-startpointsB[re]
       conne<-FALSE
       begbeg1<-beg1
       while (begbeg1>0){
            begbeg2<-beg2
            while (begbeg2>0){
                atom1<-atomsSepaAtom[begbeg1]
                indelow1<-low[atom1,]
                indeupp1<-upp[atom1,]
                atom2<-atomsSepaAtom[begbeg2]
                indelow2<-low[atom2,]
                indeupp2<-upp[atom2,]
                if (dotouchgen(indelow1,indeupp1,indelow2,indeupp2,direction)){
                   conne<-TRUE
                }
                begbeg2<-atomsRBounNext[begbeg2]
            }
            begbeg1<-atomsLBounNext[begbeg1]
       }                
       if (conne){
           linkit[do,re]<-1
      }
   }
} 
# huom ylla on nopeutettu, koska tiedetaan, etta atomit
# 1,...,mleft eivat koske toisiaan ja samoin atomit mleft+1,...,m
#
# next we apply "declev" 
rindeksitB<-seq(1,m)
res<-declevnew(rindeksitB,linkit,m)   #res is sepnum*m-matrix of 0/1
sepnum<-dim(res)[1]
# 
# res is sepnum*m-matrix, 1 in some row indicates that set (atom)
# belongs to this component, 0 in other positions
#
# output could be also a vector which contains pointers
# to a list of elements (in one list we find those sets which
# belong to the same component
#
#compopointer<-matrix(0,sepnum,1) 
#compoSet<-matrix(0,m,1)
#compoNext<-matrix(0,m,1)
#
#
#3. We join the sets 
#
# We join the sets whose startpoints are in 
# startpointsS and startpointsNewBleft, startpointsNewBright
# We have pointers separy[leftbeg] and separy[rightbeg]
# which contain pointers to lists which we can utilize
# to make a new list (these two lists contain together at most as many 
# elements as we need)
# 
# cut first list or (join these two lists and cut second)
#
TotalBeg<-separy[leftbeg]
#
tavoite<-1
hiihtaja<-TotalBeg
while ((begsSepaNext[hiihtaja]>0) && (tavoite<sepnum)){
   hiihtaja<-begsSepaNext[hiihtaja]
   tavoite<-tavoite+1
}  
if (tavoite<sepnum){  #now hiihtaja points to the end of the first list
   #join the lists
   begsSepaNext[hiihtaja]<-separy[rightbeg]
   #we continue
   hiihtaja<-separy[rightbeg]
   tavoite<-tavoite+1
   while ((begsSepaNext[hiihtaja]>0) && (tavoite<sepnum)){
      hiihtaja<-begsSepaNext[hiihtaja]
      tavoite<-tavoite+1
   }    
   begsSepaNext[hiihtaja]<-0
}
else{  #we have reached goal, cut without joining
   begsSepaNext[hiihtaja]<-0
}
#
#
nykyinen<-TotalBeg
i<-1
while (i<= sepnum){
  len<-sum(res[i,])            # number of sets to be joined
  #
  # we find vectors which contain pointer to the beginnings
  # of lists of atoms
  #
  osoittajaS<-matrix(0,len,1)  #make vector of pointers to the begs of sets
  osoittajaNewBleft<-matrix(0,len,1)
  osoittajaNewBright<-matrix(0,len,1)
  laskuri<-1
  for (j in 1:m){
     if (res[i,j]==1){
          osoittajaS[laskuri]<-startpointsS[j]  
          osoittajaNewBleft[laskuri]<-startpointsNewBleft[j]    #could be 0
          osoittajaNewBright[laskuri]<-startpointsNewBright[j]  #could be 0
          laskuri<-laskuri+1
     }    
  }
  #
  # handle separy 
  #
  begsSepaBegs[nykyinen]<-osoittajaS[1]    #always non-zero
  #
  k<-1
  while (k<=(len-1)){    
      curre<-osoittajaS[k]
      while (atomsSepaNext[curre]>0){    #find the end
          curre<-atomsSepaNext[curre]
      }
      atomsSepaNext[curre]<-osoittajaS[k+1]
      k<-k+1
  }
  #
  # handle left boundary
  #
  # set kL=0 if all 0 , otherwise kL first nonzero
  #
  k<-1
  while ((k<=len) && (osoittajaNewBleft[k]==0)){
      k<-k+1
  }
  if (k>len){   # all zero
     kL<-0
     begsLeftBoun[nykyinen]<-0
  }
  else{         # kL is first non zero
     kL<-k
     begsLeftBoun[nykyinen]<-osoittajaNewBleft[kL]  
  #
  # update the list of left boundaries
  # concatenate the lists of atoms
  #
  k<-kL
  while (k<=(len-1)){    
      curre<-osoittajaNewBleft[k]         # curre is not zero
      while (atomsLBounNext[curre]>0){    #find the end
            curre<-atomsLBounNext[curre]
      }
      # find the next non zero
      k<-k+1
      while ((k<=len) && (osoittajaNewBleft[k]==0)){
           k<-k+1
      }
      if (k>len){
          atomsLBounNext[curre]<-0
      }
      else{  # found nonzero
          atomsLBounNext[curre]<-osoittajaNewBleft[k]
      }
  }
  }
  #
  # handle right boundary
  #
  # set kR=0 if all 0 , otherwise kR first nonzero
  #
  k<-1
  while ((k<=len) && (osoittajaNewBright[k]==0)){
      k<-k+1
  }
  if (k>len){
     kR<-0
     begsRighBoun[nykyinen]<-0
  }
  else{
     kR<-k
     begsRighBoun[nykyinen]<-osoittajaNewBright[kR]  
  #
  # update the list of right boundaries
  # concatenate the lists of atoms
  #
  k<-kR
  while (k<=(len-1)){    
      curre<-osoittajaNewBright[k]         # curre is not zero
      while (atomsRBounNext[curre]>0){    #find the end
            curre<-atomsRBounNext[curre]
      }
      # find the next non zero
      k<-k+1
      while ((k<=len) && (osoittajaNewBright[k]==0)){
           k<-k+1
      }
      if (k>len){
          atomsRBounNext[curre]<-0
      }
      else{  # found nonzero
          atomsRBounNext[curre]<-osoittajaNewBright[k]
      }
  }
  }
  #
  # we move to the next sepaset
  nykyinen<-begsSepaNext[nykyinen]
  i<-i+1
}
#
return(list(totbegSepary=TotalBeg,separy=separy,
begsSepaNext=begsSepaNext,begsSepaBegs=begsSepaBegs,
begsLeftBoun=begsLeftBoun,begsRighBoun=begsRighBoun,
atomsSepaNext=atomsSepaNext,atomsSepaAtom=atomsSepaAtom,
atomsLBounNext=atomsLBounNext,atomsRBounNext=atomsRBounNext))
}























joinconne<-function(leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index){
#
#
# 1. We find the startpoints for the sets
#
# we will make three vectors of startpoints:
# a. startpoints for the whole sets (concatenate left and right child)
#    (startpointsS)
# b. startpoints for the right boundary of the left child and 
#    left boundary of the right child
#    (startpointsB)
# c.1. startpoints for the left boundary of the left child a
#    (startpointsNewBleft)
# c.2. startpoints for the right boundary of the right child
#    (startpointsNewBright)
#
#
# startpoints are pointers to 
# a. atomsSepaAtoms/atomsSepaNext
# b. atomsSepaAtoms/atomsRBounNext, atomsSepaAtoms/atomsLBounNext
# c.1. atomsSepaAtoms/atomsLBounNext
# c.2.  atomsSepaAtoms/atomsRBounNext
#
# startpoints are found from 
# a. begsSepaBegs
# b. begsLeftBoun, begsRighBoun
# c.1 begsLeftBoun, 
# c.2 begsRighBoun
#
# b. is used to check which touch
# a., c.1. and c.2. are joined together
# Note that some sets of the boundary are empty 
# (we store 0 to the respective location in begsLeftBoun, begsRighBoun )
#
suppi<-length(begsSepaNext)
startpointsS<-matrix(0,suppi,1)
startpointsB<-matrix(0,suppi,1)
startpointsNewBleft<-matrix(0,suppi,1)
startpointsNewBright<-matrix(0,suppi,1)
#new boundary: left bound. of left child, right b. of right child
#
induksi<-1
anfang<-separy[leftbeg]
startpointsS[induksi]<-anfang
startpointsB[induksi]<-begsRighBoun[anfang]
startpointsNewBleft[induksi]<-begsLeftBoun[anfang]
while (begsSepaNext[anfang]>0){
  anfang<-begsSepaNext[anfang]
  induksi<-induksi+1
  startpointsS[induksi]<-begsSepaBegs[anfang]
  startpointsB[induksi]<-begsRighBoun[anfang]
  startpointsNewBleft[induksi]<-begsLeftBoun[anfang]  
}
mleft<-induksi
induksi<-induksi+1
anfang<-separy[rightbeg]
startpointsS[induksi]<-anfang
startpointsB[induksi]<-begsLeftBoun[anfang]
startpointsNewBright[induksi]<-begsRighBoun[anfang]
while (begsSepaNext[anfang]>0){
  anfang<-begsSepaNext[anfang]
  induksi<-induksi+1
  startpointsS[induksi]<-begsSepaBegs[anfang]
  startpointsB[induksi]<-begsLeftBoun[anfang]
  startpointsNewBright[induksi]<-begsRighBoun[anfang]  

}
startpointsS<-startpointsS[1:induksi]
startpointsB<-startpointsB[1:induksi]
startpointsNewBleft<-startpointsNewBleft[1:induksi]
startpointsNewBright<-startpointsNewBright[1:induksi]
m<-induksi
mright<-m-mleft  
#
#
# 2. We make "links" matrix and apply declev
#
# We utilize previous programs
#
linkit<-matrix(0,m,m)
do<-1
while (do <= mleft){
   beg1<-startpointsB[do]    #could be 0
   re<-mleft+1
   while (re <= m){
       beg2<-startpointsB[re]    #could be 0
       conne<-FALSE
       begbeg1<-beg1
       while (begbeg1>0){
            begbeg2<-beg2
            while (begbeg2>0){
                atom1<-atomsSepaAtom[begbeg1]
                inde1<-index[atom1,]
                atom2<-atomsSepaAtom[begbeg2]
                inde2<-index[atom2,]
                if (dotouch(inde1,inde2,direction)){
                   conne<-TRUE
                }
                begbeg2<-atomsLBounNext[begbeg2]
            }
            begbeg1<-atomsRBounNext[begbeg1]
       }                
       if (conne){
           linkit[do,re]<-1
       }
       re<-re+1
   }
   do<-do+1
}
for (do in (mleft+1):m){
   beg1<-startpointsB[do]
   for (re in 1:mleft){
       beg2<-startpointsB[re]
       conne<-FALSE
       begbeg1<-beg1
       while (begbeg1>0){
            begbeg2<-beg2
            while (begbeg2>0){
                atom1<-atomsSepaAtom[begbeg1]
                inde1<-index[atom1,]
                atom2<-atomsSepaAtom[begbeg2]
                inde2<-index[atom2,]
                if (dotouch(inde1,inde2,direction)){
                   conne<-TRUE
                }
                begbeg2<-atomsRBounNext[begbeg2]
            }
            begbeg1<-atomsLBounNext[begbeg1]
       }                
       if (conne){
           linkit[do,re]<-1
      }
   }
} 
# huom ylla on nopeutettu, koska tiedetaan, etta atomit
# 1,...,mleft eivat koske toisiaan ja samoin atomit mleft+1,...,m
#
# next we apply "declev" 
rindeksitB<-seq(1,m)
res<-declevnew(rindeksitB,linkit,m)   #res is sepnum*m-matrix of 0/1
sepnum<-dim(res)[1]
# 
# res is sepnum*m-matrix, 1 in some row indicates that set (atom)
# belongs to this component, 0 in other positions
#
# output could be also a vector which contains pointers
# to a list of elements (in one list we find those sets which
# belong to the same component
#
#compopointer<-matrix(0,sepnum,1) 
#compoSet<-matrix(0,m,1)
#compoNext<-matrix(0,m,1)
#
#
#3. We join the sets 
#
# We join the sets whose startpoints are in 
# startpointsS and startpointsNewBleft, startpointsNewBright
# We have pointers separy[leftbeg] and separy[rightbeg]
# which contain pointers to lists which we can utilize
# to make a new list (these two lists contain together at most as many 
# elements as we need)
# 
# cut first list or (join these two lists and cut second)
#
TotalBeg<-separy[leftbeg]
#
tavoite<-1
hiihtaja<-TotalBeg
while ((begsSepaNext[hiihtaja]>0) && (tavoite<sepnum)){
   hiihtaja<-begsSepaNext[hiihtaja]
   tavoite<-tavoite+1
}  
if (tavoite<sepnum){  #now hiihtaja points to the end of the first list
   #join the lists
   begsSepaNext[hiihtaja]<-separy[rightbeg]
   #we continue
   hiihtaja<-separy[rightbeg]
   tavoite<-tavoite+1
   while ((begsSepaNext[hiihtaja]>0) && (tavoite<sepnum)){
      hiihtaja<-begsSepaNext[hiihtaja]
      tavoite<-tavoite+1
   }    
   begsSepaNext[hiihtaja]<-0
}
else{  #we have reached goal, cut without joining
   begsSepaNext[hiihtaja]<-0
}
#
#
nykyinen<-TotalBeg
i<-1
while (i<= sepnum){
  len<-sum(res[i,])            # number of sets to be joined
  #
  # we find vectors which contain pointer to the beginnings
  # of lists of atoms
  #
  osoittajaS<-matrix(0,len,1)  #make vector of pointers to the begs of sets
  osoittajaNewBleft<-matrix(0,len,1)
  osoittajaNewBright<-matrix(0,len,1)
  laskuri<-1
  for (j in 1:m){
     if (res[i,j]==1){
          osoittajaS[laskuri]<-startpointsS[j]  
          osoittajaNewBleft[laskuri]<-startpointsNewBleft[j]    #could be 0
          osoittajaNewBright[laskuri]<-startpointsNewBright[j]  #could be 0
          laskuri<-laskuri+1
     }    
  }
  #
  # handle separy 
  #
  begsSepaBegs[nykyinen]<-osoittajaS[1]    #always non-zero
  #
  k<-1
  while (k<=(len-1)){    
      curre<-osoittajaS[k]
      while (atomsSepaNext[curre]>0){    #find the end
          curre<-atomsSepaNext[curre]
      }
      atomsSepaNext[curre]<-osoittajaS[k+1]
      k<-k+1
  }
  #
  # handle left boundary
  #
  # set kL=0 if all 0 , otherwise kL first nonzero
  #
  k<-1
  while ((k<=len) && (osoittajaNewBleft[k]==0)){
      k<-k+1
  }
  if (k>len){   # all zero
     kL<-0
     begsLeftBoun[nykyinen]<-0
  }
  else{         # kL is first non zero
     kL<-k
     begsLeftBoun[nykyinen]<-osoittajaNewBleft[kL]  
  #
  # update the list of left boundaries
  # concatenate the lists of atoms
  #
  k<-kL
  while (k<=(len-1)){    
      curre<-osoittajaNewBleft[k]         # curre is not zero
      while (atomsLBounNext[curre]>0){    #find the end
            curre<-atomsLBounNext[curre]
      }
      # find the next non zero
      k<-k+1
      while ((k<=len) && (osoittajaNewBleft[k]==0)){
           k<-k+1
      }
      if (k>len){
          atomsLBounNext[curre]<-0
      }
      else{  # found nonzero
          atomsLBounNext[curre]<-osoittajaNewBleft[k]
      }
  }
  }
  #
  # handle right boundary
  #
  # set kR=0 if all 0 , otherwise kR first nonzero
  #
  k<-1
  while ((k<=len) && (osoittajaNewBright[k]==0)){
      k<-k+1
  }
  if (k>len){
     kR<-0
     begsRighBoun[nykyinen]<-0
  }
  else{
     kR<-k
     begsRighBoun[nykyinen]<-osoittajaNewBright[kR]  
  #
  # update the list of right boundaries
  # concatenate the lists of atoms
  #
  k<-kR
  while (k<=(len-1)){    
      curre<-osoittajaNewBright[k]         # curre is not zero
      while (atomsRBounNext[curre]>0){    #find the end
            curre<-atomsRBounNext[curre]
      }
      # find the next non zero
      k<-k+1
      while ((k<=len) && (osoittajaNewBright[k]==0)){
           k<-k+1
      }
      if (k>len){
          atomsRBounNext[curre]<-0
      }
      else{  # found nonzero
          atomsRBounNext[curre]<-osoittajaNewBright[k]
      }
  }
  }
  #
  # we move to the next sepaset
  nykyinen<-begsSepaNext[nykyinen]
  i<-i+1
}
#
return(list(totbegSepary=TotalBeg,separy=separy,
begsSepaNext=begsSepaNext,begsSepaBegs=begsSepaBegs,
begsLeftBoun=begsLeftBoun,begsRighBoun=begsRighBoun,
atomsSepaNext=atomsSepaNext,atomsSepaAtom=atomsSepaAtom,
atomsLBounNext=atomsLBounNext,atomsRBounNext=atomsRBounNext))
}























joingene<-function(node,leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index){
#
            if ((leftbeg==0) || (separy[leftbeg]==0)){
#if left child does not exist    
#note that since we consider subsets of the
#terminal nodes of the original tree, it may happen
#that leftbeg>0 but left child does not exist
                separy[node]<-separy[rightbeg]
                #we need that all lists contain as many members
                #left boundary is empty, but we will make it a list
                #of empty lists
                note<-separy[node]
                while (note>0){
                      begsLeftBoun[note]<-0
                      note<-begsSepaNext[note]
                }
                # right boundary stays same as for rightbeg
            }
            else{   # eka else
                if ((rightbeg==0) || (separy[rightbeg]==0)){  
                              #right child does not exist
                   separy[node]<-separy[leftbeg]
                   # left boundary stays same as for leftbeg
                   # right boundary is empty
                   note<-separy[node]
                   while (note>0){
                       begsRighBoun[note]<-0
                       note<-begsSepaNext[note]
                   }
                } 
                else{   #toka else: both children exist
                    #check whether left boundary of right child is empty
                    Lempty<-TRUE
                    note<-separy[rightbeg]
                    while (note>0){
                        if (begsLeftBoun[note]>0){
                              Lempty<-FALSE
                        }
                        note<-begsSepaNext[note]
                     }
                     #check whether right bound of left child is empty     
                     Rempty<-TRUE
                     note<-separy[leftbeg]
                     while (note>0){
                          if (begsRighBoun[note]>0){
                                 Rempty<-FALSE
                          }
                          note<-begsSepaNext[note]
                     }
                     #check whether one of boundaries is empty
                     if (Lempty || Rempty){
                             #one of boundaries is empty
############
#concatenating separate parts
#and updating boundaries for the separate parts
#separy[node]<- concatenate separy[leftbeg],separy[rightbeg]
###########
akku<-separy[leftbeg]
begsRighBoun[akku]<-0 #right boundaries of sets in left child are empty
                      # begsLeftBoun[akku] does not change
while (begsSepaNext[akku]>0){
  akku<-begsSepaNext[akku]
  begsRighBoun[akku]<-0
}                           
begsSepaNext[akku]<-separy[rightbeg] #concatenate list of separate sets
separy[node]<-separy[leftbeg]
akku<-separy[rightbeg]
begsLeftBoun[akku]<-0 #left boundaries of sets in right child are empty
while (begsSepaNext[akku]>0){
  akku<-begsSepaNext[akku]
  begsLeftBoun[akku]<-0
}        
####################
#end of concatenating
###################
                    }
                    else{  #both children exist, both boundaries non-empty  
jc<-joinconne(leftbeg,rightbeg,separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
direction,index)   #direction<-i
#
separy<-jc$separy
separy[node]<-jc$totbegSepary 
#
begsSepaNext<-jc$begsSepaNext
begsSepaBegs<-jc$begsSepaBegs
begsLeftBoun<-jc$begsLeftBoun
begsRighBoun<-jc$begsRighBoun
#
atomsSepaNext<-jc$atomsSepaNext
atomsSepaAtom<-jc$atomsSepaAtom
atomsLBounNext<-jc$atomsLBounNext
atomsRBounNext<-jc$atomsRBounNext
                        #
                    }
                } #toka else
            } #eka else
###########################################
#          end of child joining
###########################################
return(list(separy=separy,
begsSepaNext=begsSepaNext,begsSepaBegs=begsSepaBegs,
begsLeftBoun=begsLeftBoun,begsRighBoun=begsRighBoun,
atomsSepaNext=atomsSepaNext,atomsSepaAtom=atomsSepaAtom,
atomsLBounNext=atomsLBounNext,atomsRBounNext=atomsRBounNext))
}










kereva<-function(dendat,h,N,kernel="epane",trunc=3,threshold=0.0000001,
hw=NULL,weig=NULL)
{
#weig=rep(1/dim(dendat)[1],dim(dendat)[1]))

#source("~/kerle/profkernCRC.R")
#dyn.load("/home/jsk/kerle/kerCeva")
#dyn.load("/home/jsk/kerle/kerleCversio")
#pk2<-profkernCRC(dendat,h,N,Q)

#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(8,8)
#Q<-3

n<-dim(dendat)[1]
d<-dim(dendat)[2]  #length(N)

if (kernel=="gauss") h<-h*trunc   #trunc<-3

if (is.null(weig)) weig<-rep(1/n,n) 

if (!is.null(hw)){
   weig<-weightsit(n,hw)

   dendatnew<-dendat
   weignew<-weig
   cumul<-0
   for (i in 1:n){
        if (weig[i]>0){
            cumul<-cumul+1
            dendatnew[cumul,]<-dendat[i,]
            weignew[cumul]<-weig[i] 
        }
   }
   dendat<-dendatnew[1:cumul,]
   weig<-weignew[1:cumul]
   n<-cumul
}

inweig<-matrix(0,n+1,1)
inweig[2:(n+1)]<-weig

hnum<-length(h)
mnn<-maxnodenum(dendat,h,N,n,d)
extMaxnode<-mnn$maxnode
extMaxvals<-mnn$maxpositive
{
if (hnum>1){
 inh<-matrix(0,hnum+1,1)
 inh[2:(hnum+1)]<-h
}
else{
 inh<-h
}
}
inN<-matrix(0,d+1,1)
inN[2:(d+1)]<-N

if (kernel=="radon") kertype<-3
else if (kernel=="epane") kertype<-1 
else kertype<-2  # gaussian

kg<-.C("kergrid",
               as.integer(extMaxnode),
               as.integer(extMaxvals),
               as.double(dendat),
               as.double(inh),
               as.integer(inN),
               as.integer(n),
               as.integer(hnum),
               as.integer(d),
               as.integer(kertype),
               as.double(trunc),
               as.double(threshold),
               as.double(inweig),
               ioleft = integer(extMaxnode+1),
               ioright = integer(extMaxnode+1),
               ioparent = integer(extMaxnode+1),
               infopointer = integer(extMaxnode+1),
               iolow = integer(extMaxnode+1),
               ioupp = integer(extMaxnode+1),
               value = double(hnum*extMaxvals),
               index = integer(d*extMaxvals),
               nodefinder = integer(extMaxvals),
               numpositive = integer(1),
               numnode = integer(1),
PACKAGE = "denpro")

#left<-kg$ioleft[2:(kg$numnode+1)]
#right<-kg$ioright[2:(kg$numnode+1)]
#parent<-kg$ioparent[2:(kg$numnode+1)]
#infopointer<-kg$infopointer[2:(kg$numnode+1)]
#iolow<-kg$iolow[2:(kg$numnode+1)]
#ioupp<-kg$ioupp[2:(kg$numnode+1)]

value<-kg$value[2:(kg$numpositive+1)]
#nodefinder<-kg$nodefinder[2:(kg$numpositive+1)]
vecindex<-kg$index[2:(d*kg$numpositive+1)]
index<-matrix(0,kg$numpositive,d)
for (i in 1:kg$numpositive){
  for (j in 1:d){
     index[i,j]<-vecindex[(i-1)*d+j]
  }
}

#return(list(left=left,right=right,parent=parent,infopointer=infopointer,
#low=low,upp=upp,value=value,index=index,nodefinder=nodefinder))

suppo<-matrix(0,2*d,1)
for (i in 1:d){
   suppo[2*i-1]<-min(dendat[,i])-h
   suppo[2*i]<-max(dendat[,i])+h
}

step<-matrix(0,d,1)
for (i in 1:d) step[i]=(suppo[2*i]-suppo[2*i-1])/N[i];

recnum<-dim(index)[1]
low<-matrix(0,recnum,d)
upp<-matrix(0,recnum,d)
for (i in 1:recnum){
     low[i,]<-index[i,]-1
     upp[i,]<-index[i,]
}

return(list(value=value,index=index,
down=low,high=upp,N=N,step=step,support=suppo,n=n))

}
kergrid<-function(dendat,h,N){
#
#dendat is n*d- matrix of observations, 
#h is vector of positive smoothing parameters
#N is d-vector of the (dyadic) number of grid points for each direction
#
#dendat<-matrix(rnorm(20),10)
#h<-c(0.8,1,1.2)
#N<-c(4,4)
#
hnum<-length(h)
n<-dim(dendat)[1]
d<-dim(dendat)[2]
depth<-log(N,base=2)   
depoftree<-sum(depth)+1
#
minim<-matrix(0,d,1)  #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-min(dendat[,i])  
  maxim[i]<-max(dendat[,i])
}
hmax<-max(h)
delta<-(maxim-minim+2*hmax)/(N+1)
#
mindelta<-min(delta)
maxpositive<-n*(2*hmax/mindelta)^d
bigd<-sum(log(N,base=2))
maxnode<-bigd*ceiling(maxpositive)
#
numnode<-1
left<-matrix(0,maxnode,1)
right<-matrix(0,maxnode,1)
parent<-matrix(0,maxnode,1)
infopointer<-matrix(0,maxnode,1)
low<-matrix(0,maxnode,1)
low[1]<-1
upp<-matrix(0,maxnode,1)
upp[1]<-N[1]
#
numpositive<-0
value<-matrix(0,maxpositive,hnum)
index<-matrix(0,maxpositive,d)
nodefinder<-matrix(0,maxpositive,1)
#
gridlow<-matrix(0,d,1)
gridupp<-matrix(0,d,1)
#
for (i in 1:n){
 for (hrun in 1:hnum){ 
   #find the grid points in the support 
   for (j in 1:d){  
      gridlow[j]<-floor(((dendat[i,j]-minim[j])/delta[j])+1)
      gridupp[j]<-ceiling(((dendat[i,j]-minim[j]+2*h[hrun])/delta[j])-1)
   }
   base<-gridupp-gridlow+1
   gridcard<-prod(base)
   k<-0
   while (k<=(gridcard-1)){
      if (d>1){  
          inde<-digit(k,base)   #inde is d-vector
          inde<-inde+gridlow
      }
      else{
          inde<-gridlow+k
      }
      point<-minim-h[hrun]+delta*inde     #point is d-vector  
      val<-epane(point-dendat[i,],h[hrun])
      #find whether gridpoint is already in tree
      fe<-findend(inde,left,right,low,upp,N)
      if (fe$exists){
           pointer<-infopointer[fe$location]
           curval<-value[pointer,hrun]
           value[pointer,hrun]<-curval+val/n
      }
      else{  #gridpoint was not yet in the tree
         curre<-fe$location
         curdep<-fe$dep
         #
         ad<-addnode(inde,curre,curdep,left,right,parent,low,upp,N,numnode)
         numnode<-ad$numnode
         left<-ad$left
         right<-ad$right
         parent<-ad$parent
         low<-ad$low
         upp<-ad$upp
         nodeloc<-ad$nodeloc
         #
         numpositive<-numpositive+1
         infopointer[numnode]<-numpositive
         value[numpositive,hrun]<-val/n
         index[numpositive,]<-inde
         nodefinder[numpositive]<-nodeloc
      }
      k<-k+1 
   }
 }
}
left<-left[1:numnode]
right<-right[1:numnode]
parent<-parent[1:numnode]
infopointer<-infopointer[1:numnode]
#deplink<-deplink[1:numnode]
low<-low[1:numnode]
upp<-upp[1:numnode]
#
value<-value[1:numpositive,]
index<-index[1:numpositive,]
nodefinder<-nodefinder[1:numpositive]
return(list(left=left,right=right,parent=parent,infopointer=infopointer,
low=low,upp=upp,value=value,index=index,nodefinder=nodefinder))
}                              








kernesti.dens<-function(arg,x,h=1,kernel="gauss",g=NULL,gernel="gauss")
{
d<-length(arg)

if (d>1){

if (length(h)==1) h<-rep(h,d)

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ ans<-(rowSums(xx^2) <= 1) 
                      return( ans ) }

argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
xxx<-sweep(x-argu,2,h,"/")
w<-ker(xxx)/prod(h)
est<-sum(w)/length(w)

if (!is.null(g)){

   n<-dim(x)[1]
   if (gernel=="bart") 
   ger<-function(xx){ return( (1-rowSums(xx^2))*(rowSums(xx^2)<= 1) ) }
   if (gernel=="gauss") 
   ger<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   if (gernel=="uniform") 
   ger<-function(xx){ ans<-(rowSums(xx^2)<= 1) 
                      return( ans ) }

   argui<-matrix(seq(n,1,-1),n,1)
   w<-ker((x-argu)/h)/prod(h)*ger((n-argui)/g)/g
   est<-sum(w)/length(w)
}
}
else{  # d==1  #########################################

if (kernel=="gauss") ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }

x<-matrix(x,length(x),1)
w<-ker((x-arg)/h)/h^d   #weights<-w/sum(w)
est<-sum(w)/length(w)

if (!is.null(g)){

   n<-length(x)
   if (gernel=="bart") 
   ger<-function(xx){ return( (1-rowSums(xx^2))*(rowSums(xx^2)<= 1) ) }
   if (gernel=="gauss") 
   ger<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   if (gernel=="uniform") 
   ger<-function(xx){ ans<-(rowSums(xx^2)<= 1) 
                      return( ans ) }

   argui<-matrix(seq(n,1,-1),n,1)
   w<-ker((x-arg)/h)/h^d*ger((n-argui)/g)/g
   est<-sum(w)/length(w)
}

}

return(est)
}



lambda2emass<-function(lambda,m,M,sig,p,support=NULL,seed=1,mul=2)
{
#m is the number of Monte Carlo samples
#M is l*d-matrix, rows are the means
#sig is l*d-matrix, for l:th mixture d covariances
#p is l-vector, proportion for each mixture

set.seed(seed)
l<-dim(M)[1]
d<-dim(M)[2]
if (is.null(support)){
   support<-matrix(0,2*d,1)
   for (i in 1:d){
       support[2*i-1]<-min(M[,i]-mul*sig[,i])
       support[2*i]<-max(M[,i]+mul*sig[,i])
   }
}

maksi<-0
for (i in 1:l){
    zig<-sig[i,]
    maksi<-maksi+p[i]*evanor(0)/prod(zig)
}

boxvol<-1
for (i in 1:d) boxvol<-boxvol*(support[2*i]-support[2*i-1])
boxvol<-boxvol*maksi

inside<-0
for (i in 1:m){
    x<-matrix(0,d,1)
    ran<-runif(d+1)
    for (j in 1:d){
        beg<-support[2*j-1]
        end<-support[2*j] 
        x[j]<-beg+(end-beg)*ran[j]
    }
    y<-0+(maksi-0)*ran[d+1]

    arvo<-0
    for (j in 1:l){
        zig<-sig[j,]
        mu<-M[j,]
        arvo<-arvo+p[j]*evanor((x-mu)/zig)/prod(zig)
    }
    
    if ((y<=arvo)&&(y>=lambda)) inside<-inside+1
}

emass<-boxvol*inside/m
return(emass)
}


leafsfirst.adagrid<-function(pcf)
{
down<-pcf$down
high<-pcf$high
support<-pcf$support
grid<-pcf$grid
value<-pcf$value

d<-dim(down)[2]
lkm<-dim(down)[1]

distat<-pcf$value
infopointer<-seq(1,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
number<-matrix(0,lkm,1)
atomlist<-matrix(0,lkm,lkm)
atomnumb<-matrix(0,lkm,1)
center<-matrix(0,d,lkm)
distcenter<-matrix(0,lkm,d)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]

downi<-down[infopointer[node],]
highi<-high[infopointer[node],]
simp<-matrix(0,2*d,1)
for (i in 1:d){
    simp[2*i-1]<-grid[downi[i],i]
    simp[2*i]<-grid[highi[i],i]
}
volume[node]<-massone(simp)
for (dd in 1:d){
   center[dd,node]<-(simp[2*dd-1]+simp[2*dd])/2
}
number[node]<-1
atomlist[node,1]<-infopointer[node]
atomnumb[node]<-1

beg<-node                 #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-down[note,i]   
  boundrec[node,2*i]<-high[note,i]  
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)    #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-down[note,i]  
         rec1[2*i]<-high[note,i] 
    }
    boundrec[node,]<-rec1

    # radius
    radius[node]<-distat[ord[node]]

    downi<-down[infopointer[node],]
    highi<-high[infopointer[node],]
    simp<-matrix(0,2*d,1)
    for (i in 1:d){
       simp[2*i-1]<-grid[downi[i],i]
       simp[2*i]<-grid[highi[i],i]
    }
    volume[node]<-massone(simp) 
    for (dd in 1:d){
       center[dd,node]<-(simp[2*dd-1]+simp[2*dd])/2
    } 
    number[node]<-1
    atomlist[node,1]<-infopointer[node]
    atomnumb[node]<-1

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
       istouch<-touchstep(node,curroot,boundrec,child,sibling,infopointer,down,high,)

        if (istouch==1){

           # paivita parent, child, sibling, volume 
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }

           number[node]<-number[node]+number[curroot]
           center[,node]<-(center[,node]*volume[node]+center[,curroot]*volume[curroot])/(volume[node]+volume[curroot])
           volume[node]<-volume[node]+volume[curroot]


           atomlist[node,(atomnumb[node]+1):(atomnumb[node]+atomnumb[curroot])]<-atomlist[curroot,1:atomnumb[curroot]]
           atomnumb[node]<-atomnumb[curroot]+atomnumb[node]

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]

        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter

maxdis<-distat[ord[length(ord)]]
#center<-matrix(0,d,lkm)
#for (i in 1:lkm){
#    for (j in 1:d){
#       ala<-grid[down[infopointer[i],j],j]
#       yla<-grid[high[infopointer[i],j],j]
#       center[j,i]<-(ala+yla)/2
#    }
#}

lf<-list(
parent=parent,volume=volume,center=center,level=radius,
root=root,
infopointer=infopointer,
maxdis=maxdis,
dendat=pcf$dendat,
atomlist=atomlist,atomnumb=atomnumb)

return(lf)
}


leafsfirst.bondary<-function(pcf=NULL,lev=NULL,refe=NULL,type="lst",levmet="radius",
ordmet="etaisrec",ngrid=NULL,
dendat=NULL,rho=0,propor=NULL,lowest="dens",f=NULL)
{
# pcf is a piecewise constant object
# type= "lst"/"shape"
# levmet= "radius"/"proba"

if (lowest=="dens") lowest<-0 else lowest<-min(pcf$value)

if ((!is.null(lev)) || (!is.null(propor))){
    type<-"shape"
    if (!is.null(propor)) lev<-propor*max(pcf$value)
    if (is.null(refe)) refe<-locofmax(pcf)
}
if (!is.null(dendat)) type<-"tail"

if (type=="tail"){
   d<-dim(dendat)[2]
   pcf$high<-dendat
   pcf$down<-dendat
   if (is.null(refe)){
       refe<-matrix(0,1,d)
       for (i in 1:d) refe[1,i]<-mean(dendat[,i])
       refe<-refe[1:d]
   }
}
else{
  d<-length(pcf$N)
  step<-matrix(0,d,1)
  for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
}

if (type=="lst"){
  lkm<-length(pcf$value)
  distat<-pcf$value-lowest
  infopointer<-seq(1,lkm)     # links from nodes to recs
}
else if (type=="shape"){
  lenni<-length(pcf$value)
  distat<-matrix(0,lenni,1)
  infopointer<-matrix(0,lenni,1)
  lkm<-0
  for (i in 1:lenni){
    if (pcf$value[i]>=lev){
       lkm<-lkm+1
       nod<-i  #nod<-pcf$nodefinder[i]
       if (ordmet=="etaisrec"){
           recci<-matrix(0,2*d,1)
           for (jj in 1:d){
              recci[2*jj-1]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
              recci[2*jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
           }
           distat[lkm]<-etaisrec(refe,recci)
       }
       else{
          lowi<-matrix(0,d,1)
          uppi<-matrix(0,d,1)
          for (jj in 1:d){
             lowi[jj]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
             uppi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
          }
          baryc<-lowi+(uppi-lowi)/2  
          distat[lkm]<-etais(baryc,refe)
       }
       infopointer[lkm]<-i
    }
  }
}
else{  #type=="tail"
   d<-dim(dendat)[2]
   n<-dim(dendat)[1]
   lkm<-dim(dendat)[1]
   distat<-sqrt(pituus(dendat-t(matrix(refe,d,n))))
   infopointer<-seq(1,lkm)
}

distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   
#if (length(rho)==1) rho<-rep(rho,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
proba<-matrix(0,lkm,1)
ekamome<-matrix(0,lkm,d)
distcenter<-matrix(0,lkm,d)
branchradius<-matrix(0,lkm,1)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]
branchradius[node]<-radius[node]

if (type!="tail"){
  # volume calculation
  vol<-1
  k<-1
  ip<-infopointer[node]  #pcf$nodefinder[infopointer[node]]
  while (k<=d){
      vol<-vol*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
      k<-k+1
  }
  volume[node]<-vol
  ip2<-infopointer[node]
  proba[node]<-pcf$value[ip2]*vol

  # ekamome calculation
  newcente<-matrix(0,d,1)
  for (j in 1:d){
    volmin<-1
    k<-1
    while (k<=d){
       if (k!=j){
          volmin<-volmin*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
       }
       k<-k+1
    }
    ala<-pcf$support[2*j-1]+step[j]*pcf$down[ip,j]
    yla<-pcf$support[2*j-1]+step[j]*pcf$high[ip,j]
    newcente[j]<-volmin*(yla^2-ala^2)/2
  }
  ekamome[node,]<-newcente
  distcenter[node,]<-newcente/vol
}
else{  # type=tail
  if (is.null(f)) volume[node]<-1
  else{ 
       ip<-infopointer[node] 
       volume[node]<-1/(f[ip]*length(f))
  }
}

beg<-node             #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcf$down[note,i]   
  boundrec[node,2*i]<-pcf$high[note,i]  
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)  #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-pcf$down[note,i]  
         rec1[2*i]<-pcf$high[note,i] 
    }
    boundrec[node,]<-rec1

    # radius
    radius[node]<-distat[ord[node]]
    branchradius[node]<-radius[node]
    if (type!="tail"){
       # volume calculation
       vol<-1
       k<-1
       ip<-infopointer[node]    #pcf$nodefinder[infopointer[node]]
       while (k<=d){
          vol<-vol*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
          k<-k+1
       }
       volume[node]<-vol
       ip2<-infopointer[node]
       proba[node]<-pcf$value[ip2]*vol

       # ekamome calculation
       newcente<-matrix(0,d,1)
       for (jj in 1:d){
            volmin<-1
            k<-1
            while (k<=d){
               if (k!=jj){
                   volmin<-volmin*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
               }
               k<-k+1
            }
            ala<-pcf$support[2*jj-1]+step[jj]*pcf$down[ip,jj]
            yla<-pcf$support[2*jj-1]+step[jj]*pcf$high[ip,jj]
            newcente[jj]<-volmin*(yla^2-ala^2)/2
       }
       ekamome[node,]<-newcente
       distcenter[node,]<-newcente/vol
    }
    else{     #type==tail
       if (is.null(f)) volume[node]<-1
       else{ 
          ip<-infopointer[node] 
          volume[node]<-1/(f[ip]*length(f))
       }
    }

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        rhocur<-rho   #rho[infopointer[node]]  
        istouch<-touchstep.boundary(node,curroot,boundrec,child,sibling,
                           infopointer,pcf$down,pcf$high,rhocur)
        if (istouch==1){
{
           # paivita parent, child, sibling, volume ekamome
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }
           if (type!="tail"){
              volume[node]<-volume[node]+volume[curroot]
              proba[node]<-proba[node]+proba[curroot]
              ekamome[node,]<-ekamome[node,]+ekamome[curroot,]
           }
           else{  # type == tail
              volume[node]<-volume[node]+volume[curroot]
           }

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])
           if (branchradius[node]<=branchradius[curroot]) 
                  distcenter[node,]<-distcenter[curroot,]
           branchradius[node]<-max(branchradius[node],branchradius[curroot])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]
}
        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter
# lf is the level set tree or the shape tree
if (type!="tail"){
   for (i in 1:lkm){
      for (j in 1:d){
          ekamome[i,j]<-ekamome[i,j]/volume[i]
      }
   }
   bary<-ekamome[root,]
}
if (type=="shape"){
  maxdis<-sqrt(distat[ord[length(ord)]])
  if (levmet=="proba")
     level<-taillevel(root,#child,sibling,
            parent,volume,proba)
  else 
     level<-sqrt(radius)
}
else{ #type="lst"
     level<-radius+lowest
     maxdis<-distat[ord[length(ord)]]
}
if (type=="tail"){
   center<-t(dendat[infopointer,])
}

if (type!="tail"){
  lf<-list(
  parent=parent,volume=volume,center=t(ekamome),level=level,
  root=root,
  #child=child,sibling=sibling,  #virhe??
  infopointer=infopointer,
  proba=proba,#radius=radius,
  #branchradius=sqrt(branchradius),
  distcenter=t(distcenter),
  refe=refe,maxdis=maxdis,bary=bary,lev=lev)
}
else{
  lf<-list(
  parent=parent,volume=volume,center=center,level=level,
  root=root,
  #child=child,sibling=sibling,  #virhe??
  infopointer=infopointer,
  #proba=proba,#radius=radius,
  #branchradius=sqrt(branchradius),
  #distcenter=t(distcenter),
  refe=refe,maxdis=maxdis,
  dendat=dendat)
}

# if ngrid given, reduce the lst
if (!is.null(ngrid)){
    stepsi<-maxdis/ngrid
    radii<-seq(0,maxdis,stepsi)
    lf<-treedisc(lf,pcf,r=radii,type=type)
}

return(lf)
}





leafsfirst.complex<-function(dendat,f,complex,rho=0)
{
# complex is lkm*(d+1) matrix: pointers to dendat
# lambdas is lkm vector of levels

d<-dim(dendat)[2]  #dim(complex)[2]-1
lkm<-dim(complex)[1]

lambdas<-matrix(0,lkm,1)
mids<-matrix(0,lkm,d)
for (i in 1:lkm){
    vs<-complex[i,]
    vals<-f[vs]
    lambdas[i]<-min(vals)
    mids[i,1]<-mean(dendat[vs,1])
    mids[i,2]<-mean(dendat[vs,2])
}

pcfhigh<-mids+rho/2
pcfdown<-mids-rho/2

distat<-lambdas
infopointer<-seq(1,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
number<-matrix(0,lkm,1)
atomlist<-matrix(0,lkm,lkm)
atomnumb<-matrix(0,lkm,1)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]

simple<-complex[infopointer[node],]
simp<-dendat[simple,]
volume[node]<-voltriangle(simp)   #kappa*pi*rho^2
number[node]<-1
atomlist[node,1]<-infopointer[node]
atomnumb[node]<-1

beg<-node                 #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcfdown[note,i]   
  boundrec[node,2*i]<-pcfhigh[note,i]  
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)    #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-pcfdown[note,i]  
         rec1[2*i]<-pcfhigh[note,i] 
    }
    boundrec[node,]<-rec1

    # radius
    radius[node]<-distat[ord[node]]

    simple<-complex[infopointer[node],]
    simp<-dendat[simple,]
    volume[node]<-voltriangle(simp)  #kappa*pi*rho[infopointer[node]]^2
    number[node]<-1
    atomlist[node,1]<-infopointer[node]
    atomnumb[node]<-1

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        istouch<-touchstep.complex(node,curroot,boundrec,child,sibling,
                                infopointer,pcfdown,pcfhigh,dendat,complex)
        if (istouch==1){

           # paivita parent, child, sibling, volume 
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }

           number[node]<-number[node]+number[curroot]
           volume[node]<-volume[node]+volume[curroot]
                         #kappa*number[node]*pi*rho[1]^2
           atomlist[node,(atomnumb[node]+1):(atomnumb[node]+atomnumb[curroot])]<-atomlist[curroot,1:atomnumb[curroot]]
           atomnumb[node]<-atomnumb[curroot]+atomnumb[node]

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]

        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter

maxdis<-distat[ord[length(ord)]]
center<-t(mids[infopointer,])

lf<-list(
parent=parent,volume=volume,center=center,level=radius,
root=root,
infopointer=infopointer,
maxdis=maxdis,
dendat=dendat,rho=rho,
atomlist=atomlist,atomnumb=atomnumb)

return(lf)
}


leafsfirst.complex.volu<-function(lst,dendat,complex,rho,vols,M=1000,grid=1,
seed=1)
{
itemnum<-length(lst$volume)
volume<-matrix(0,itemnum,1)
kapat<-matrix(0,itemnum,1)
d<-dim(dendat)[2]

if (grid==0){

for (note in 1:itemnum){
  atomit<-lst$atomlist[note,1:lst$atomnumb[note]]
  pisteet<-matrix(complex[atomit,],lst$atomnumb[note],d+1)
  voltti<-montecarlo.complex(dendat,pisteet,rho,M,seed=seed)
  #if (lst$parent[note]>0) voltti<-min(voltti,lst$volume[lst$parent[note]])
  volume[note]<-voltti
}
lst$volume<-volume

}
else{

volume.root<-montecarlo.complex(dendat,complex,rho,M,seed=seed)
volume.sum<-sum(vols)  #itemnum*pi*rho^2
kappa<-volume.root/volume.sum

# 1) lasketaan kapat grid kpl

kapat.lyhyt<-matrix(0,grid,1)
levet.lyhyt<-matrix(0,grid,1)
kapat.lyhyt[1]<-kappa
levet.lyhyt[1]<-1
if (grid>1){
   levstep<-floor(itemnum/grid)
   or<-order(lst$level)
   for (i in 2:grid){
         levlok<-(i-1)*levstep
         note<-or[levlok]
         atomit<-lst$atomlist[note,1:lst$atomnumb[note]]
         pisteet<-matrix(complex[atomit,],lst$atomnumb[note],d+1)
         volume.nyt<-montecarlo.complex(dendat,pisteet,rho,M,seed=seed)
         volume.sum<-sum(vols[atomit])    #lst$atomnumb[note]*rho^2/2
         kapat.lyhyt[i]<-volume.nyt/volume.sum
         kapat.lyhyt[i]<-min(kapat.lyhyt[i],kapat.lyhyt[i-1])
         levet.lyhyt[i]<-levlok    
   }
}

# 2) interpoloidaan muut kapat

ra<-rank(lst$level)
for (i in 1:itemnum){
    ranko<-ra[i]
    lohko<-ceiling(grid*ranko/itemnum)
    kapa.ala<-kapat.lyhyt[lohko]
    if (lohko<grid) kapa.yla<-kapat.lyhyt[lohko+1] 
    else            kapa.yla<-kapat.lyhyt[grid]
    leve.ala<-levet.lyhyt[lohko]
    if (lohko<grid) leve.yla<-levet.lyhyt[lohko+1] 
    else            leve.yla<-itemnum
    kappa<-kapa.ala+(ranko-leve.ala)*(kapa.yla-kapa.ala)/(leve.yla-leve.ala)
    kapat[i]<-kappa
    atomit<-lst$atomlist[i,1:lst$atomnumb[i]]
    volume.pot<-kappa*sum(vols[atomit]) 
    #volume.pot<-kappa*lst$atomnumb[i]*rho^2/2
    #if (lst$parent[i]>0) volume[i]<-min(volume.pot,lst$volume[lst$parent[i]])
    volume[i]<-volume.pot
}
lst$volume<-volume

}

lst$volume<-volume
return(lst)
}

leafsfirst.delaunay<-function(dendat,complex,fs,rho=0)
{
# complex is lkm*(d+1) matrix: pointers to dendat
# fs is lkm vector of values of the function 
# lambdas is lkm vector of levels

d<-dim(dendat)[2]  #dim(complex)[2]-1
lkm<-dim(complex)[1]
lambdas<-fs

mids<-matrix(0,lkm,d)
for (i in 1:lkm){
    vs<-complex[i,]
    mids[i,1]<-mean(dendat[vs,1])
    mids[i,2]<-mean(dendat[vs,2])
}
pcfhigh<-mids+rho/2
pcfdown<-mids-rho/2

distat<-lambdas
infopointer<-seq(1,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
ekamome<-matrix(0,lkm,d)
number<-matrix(0,lkm,1)
atomlist<-matrix(0,lkm,lkm)
atomnumb<-matrix(0,lkm,1)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]

# volume calculation
simple<-complex[infopointer[node],]
simp<-dendat[simple,]
volume[node]<-volsimplex(simp)   #kappa*pi*rho^2

number[node]<-1
atomlist[node,1]<-infopointer[node]
atomnumb[node]<-1

# ekamome calculation
ekamome[node,]<-colSums(simp)/(d+1)*volume[node]
#ekamome[node,]<-simp[1,]*volume[node]

beg<-node                 #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcfdown[note,i]   
  boundrec[node,2*i]<-pcfhigh[note,i]  
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)    #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-pcfdown[note,i]  
         rec1[2*i]<-pcfhigh[note,i] 
    }
    boundrec[node,]<-rec1

    # radius
    radius[node]<-distat[ord[node]]

    # volume
    simple<-complex[infopointer[node],]
    simp<-dendat[simple,]
    volume[node]<-volsimplex(simp)  #kappa*pi*rho[infopointer[node]]^2

    number[node]<-1
    atomlist[node,1]<-infopointer[node]
    atomnumb[node]<-1

    # ekamome calculation
    ekamome[node,]<-colSums(simp)/(d+1)*volume[node]
    #ekamome[node,]<-simp[1,]*volume[node]

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        istouch<-touchstep.delaunay(node,curroot,boundrec,child,sibling,
                                infopointer,pcfdown,pcfhigh,dendat,complex)
        if (istouch==1){

           # paivita parent, child, sibling, volume 
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }

           number[node]<-number[node]+number[curroot]
           volume[node]<-volume[node]+volume[curroot]#kappa*number[node]*pi*rho[1]^2
           ekamome[node,]<-ekamome[node,]+ekamome[curroot,]
           atomlist[node,(atomnumb[node]+1):(atomnumb[node]+atomnumb[curroot])]<-atomlist[curroot,1:atomnumb[curroot]]
           atomnumb[node]<-atomnumb[curroot]+atomnumb[node]

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]

        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter
# lf is the level set tree or the shape tree
for (i in 1:lkm){
    for (j in 1:d){
       ekamome[i,j]<-ekamome[i,j]/volume[i]
    }
}
bary<-ekamome[root,]
center=t(ekamome)
#center<-t(mids[infopointer,])

maxdis<-distat[ord[length(ord)]]

lf<-list(
parent=parent,volume=volume,center=center,level=radius,
root=root,
infopointer=infopointer,
maxdis=maxdis,
dendat=dendat,rho=rho,
atomlist=atomlist,atomnumb=atomnumb)

return(lf)
}


leafsfirst.intpol<-function(dendat, f, rho=0, dist.type="euclid")
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]
pcfhigh<-dendat+rho
pcfdown<-dendat-rho

distat<-f
lkm<-n
infopointer<-seq(1,lkm)
if (length(rho)==1) rho<-rep(rho,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
number<-matrix(0,lkm,1)
atomlist<-matrix(0,lkm,lkm)
atomnumb<-matrix(0,lkm,1)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]

volume[node]<-1  #kappa*pi*rho[1]^2
number[node]<-1
atomlist[node,1]<-infopointer[node]
atomnumb[node]<-1

beg<-node                 #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcfdown[note,i]   
  boundrec[node,2*i]<-pcfhigh[note,i]  
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)    #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-pcfdown[note,i]  
         rec1[2*i]<-pcfhigh[note,i] 
    }
    boundrec[node,]<-rec1

    # radius
    radius[node]<-distat[ord[node]]

    volume[node]<-1  #kappa*pi*rho[infopointer[node]]^2
    number[node]<-1
    atomlist[node,1]<-infopointer[node]
    atomnumb[node]<-1

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        #rhocur<-rho[infopointer[node]]  
        istouch<-touchstep.tail(node,curroot,boundrec,child,sibling,
                                infopointer,pcfdown,pcfhigh,rho,dendat,
                                dist.type=dist.type)
        if (istouch==1){

           # paivita parent, child, sibling, volume 
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }

           number[node]<-number[node]+number[curroot]
           volume[node]<-volume[node]+volume[curroot]
                         #kappa*number[node]*pi*rho[1]^2
           atomlist[node,(atomnumb[node]+1):(atomnumb[node]+atomnumb[curroot])]<-atomlist[curroot,1:atomnumb[curroot]]
           atomnumb[node]<-atomnumb[curroot]+atomnumb[node]

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]

        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter

maxdis<-distat[ord[length(ord)]]
center<-t(dendat[infopointer,])

lf<-list(
parent=parent,volume=volume,center=center,level=radius,
root=root,
infopointer=infopointer,
maxdis=maxdis,
dendat=dendat,rho=rho,
atomlist=atomlist,atomnumb=atomnumb)

return(lf)
}









leafsfirst.intpol.volu<-function(lst, dendat, rho, M=1000, grid=1)
{
itemnum<-length(lst$volume)
volume<-matrix(0,itemnum,1)
kapat<-matrix(0,itemnum,1)
d<-dim(dendat)[2]

volume.root<-montecarlo.ball(dendat,rho,M)
volume.sum<-itemnum*pi*rho^2
kappa<-volume.root/volume.sum

# 1) lasketaan kapat grid kpl

kapat.lyhyt<-matrix(0,grid,1)
levet.lyhyt<-matrix(0,grid,1)
kapat.lyhyt[1]<-kappa
levet.lyhyt[1]<-1
if (grid>1){
   levstep<-floor(itemnum/grid)
   or<-order(lst$level)
   for (i in 2:grid){
         levlok<-(i-1)*levstep
         note<-or[levlok]
         atomit<-lst$atomlist[note,1:lst$atomnumb[note]]
         pisteet<-matrix(dendat[atomit,],lst$atomnumb[note],d)
         volume.nyt<-montecarlo.ball(pisteet,rho,M)
         volume.sum<-lst$atomnumb[note]*pi*rho^2
         kapat.lyhyt[i]<-volume.nyt/volume.sum
         kapat.lyhyt[i]<-min(kapat.lyhyt[i],kapat.lyhyt[i-1])
         levet.lyhyt[i]<-levlok    
   }
}

# 2) interpoloidaan muut kapat

ra<-rank(lst$level)
for (i in 1:itemnum){
    ranko<-ra[i]
    lohko<-ceiling(grid*ranko/itemnum)
    kapa.ala<-kapat.lyhyt[lohko]
    if (lohko<grid) kapa.yla<-kapat.lyhyt[lohko+1] 
    else            kapa.yla<-kapat.lyhyt[grid]
    leve.ala<-levet.lyhyt[lohko]
    if (lohko<grid) leve.yla<-levet.lyhyt[lohko+1] 
    else            leve.yla<-itemnum
    kappa<-kapa.ala+(ranko-leve.ala)*(kapa.yla-kapa.ala)/(leve.yla-leve.ala)
    kapat[i]<-kappa
    volume.pot<-kappa*lst$atomnumb[i]*pi*rho^2
    #if (lst$parent[i]>0) volume[i]<-min(volume.pot,lst$volume[lst$parent[i]])
    volume[i]<-volume.pot
}

lst$volume<-volume
return(lst)
}


leafsfirst.lst<-function(pcf, ngrid=NULL, predictor=NULL, type=NULL)
{
rho<-0

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

lkm<-length(pcf$value)
distat<-pcf$value
infopointer<-seq(1,lkm)     # links from nodes to recs

distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   
if (length(rho)==1) rho<-rep(rho,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
ekamome<-matrix(0,lkm,d)
distcenter<-matrix(0,lkm,d)
branchradius<-matrix(0,lkm,1)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]
branchradius[node]<-radius[node]

volume[node]<-1

beg<-node             #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcf$down[note,i]   
  boundrec[node,2*i]<-pcf$high[note,i]  
}

found.predictor.node<-FALSE
if ((!is.null(predictor))&&(!found.predictor.node)){
   predictor.rec<-matrix(0,2*d,1)
   for (ii in 1:d){ 
     predictor.rec[2*ii-1]<-floor((predictor[ii]-pcf$support[2*ii-1])/step[ii])
     predictor.rec[2*ii]<-ceiling((predictor[ii]-pcf$support[2*ii-1])/step[ii])
   }
   if (touch(predictor.rec,boundrec[node,])) predictor.node<-node
}
else predictor.node<-NULL

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)  #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-pcf$down[note,i]  
         rec1[2*i]<-pcf$high[note,i] 
    }
    boundrec[node,]<-rec1

    if ((!is.null(predictor))&&(!found.predictor.node)){
       if (touch(predictor.rec,boundrec[node,])) predictor.node<-node
    }

    # radius
    radius[node]<-distat[ord[node]]
    branchradius[node]<-radius[node]
    
    volume[node]<-1

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        rhocur<-rho[infopointer[node]]  
        istouch<-touchstep(node,curroot,boundrec,child,sibling,
                           infopointer,pcf$down,pcf$high,rhocur)
        if (istouch==1){
{
           # paivita parent, child, sibling, volume ekamome
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }
           
           volume[node]<-volume[node]+volume[curroot]

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])
           if (branchradius[node]<=branchradius[curroot]) 
                  distcenter[node,]<-distcenter[curroot,]
           branchradius[node]<-max(branchradius[node],branchradius[curroot])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]
}
        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter

for (i in 1:lkm){
      for (j in 1:d){
          ekamome[i,j]<-ekamome[i,j]/volume[i]
      }
}
bary<-ekamome[root,]

level<-radius
maxdis<-distat[ord[length(ord)]]

lf<-list(
  parent=parent,volume=volume,center=t(ekamome),level=level,
  root=root,
  infopointer=infopointer,
  distcenter=t(distcenter),
  maxdis=maxdis,bary=bary,predictor.node=predictor.node)


# if ngrid given, reduce the lst
if (!is.null(ngrid)){
    stepsi<-maxdis/ngrid
    radii<-seq(0,maxdis,stepsi)
    lf<-treedisc(lf,pcf,r=radii,type=type)
}

return(lf)
}





leafsfirst.new<-function(pcf=NULL, lev=NULL, refe=NULL, type="lst",
levmet="radius", ordmet="etaisrec", ngrid=NULL,
dendat=NULL, rho=0, propor=NULL, dist.type="euclid")
{
# pcf is a piecewise constant object
# type= "lst"/"shape"
# levmet= "radius"/"proba"

if ((!is.null(lev)) || (!is.null(propor))) type<-"shape"
if (!is.null(dendat)) type<-"tail"

if (type=="tail") 

lst<-leafsfirst.tail(dendat=dendat, rho=rho, refe=refe, dist.type=dist.type)


return(lst)
}


leafsfirst.nn<-function(pcf=NULL,lev=NULL,refe=NULL,type="lst",levmet="radius",
ordmet="etaisrec",ngrid=NULL,
dendat=NULL,rho=0,propor=NULL)
{
# pcf is a piecewise constant object
# type= "lst"/"shape"
# levmet= "radius"/"proba"

if ((!is.null(lev)) || (!is.null(propor))){
    type<-"shape"
    if (!is.null(propor)) lev<-propor*max(pcf$value)
    if (is.null(refe)) refe<-locofmax(pcf)
}
if (!is.null(dendat)) type<-"tail"

if (type=="tail"){
   d<-dim(dendat)[2]
   pcf$high<-dendat+rho  #[infopointer[node]]  
   pcf$down<-dendat-rho 
   if (is.null(refe)){
       refe<-matrix(0,1,d)
       for (i in 1:d) refe[1,i]<-mean(dendat[,i])
       refe<-refe[1:d]
   }
}
else{
  d<-length(pcf$N)
  step<-matrix(0,d,1)
  for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
}

if (type=="lst"){
  lkm<-length(pcf$value)
  distat<-pcf$value
  infopointer<-seq(1,lkm)     # links from nodes to recs
}
else if (type=="shape"){
  lenni<-length(pcf$value)
  distat<-matrix(0,lenni,1)
  infopointer<-matrix(0,lenni,1)
  lkm<-0
  for (i in 1:lenni){
    if (pcf$value[i]>=lev){
       lkm<-lkm+1
       nod<-i  #nod<-pcf$nodefinder[i]
       if (ordmet=="etaisrec"){
           recci<-matrix(0,2*d,1)
           for (jj in 1:d){
              recci[2*jj-1]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
              recci[2*jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
           }
           distat[lkm]<-etaisrec(refe,recci)
       }
       else{
          lowi<-matrix(0,d,1)
          uppi<-matrix(0,d,1)
          for (jj in 1:d){
             lowi[jj]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
             uppi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
          }
          baryc<-lowi+(uppi-lowi)/2  
          distat[lkm]<-etais(baryc,refe)
       }
       infopointer[lkm]<-i
    }
  }
}
else{  #type=="tail"
   d<-dim(dendat)[2]
   n<-dim(dendat)[1]
   lkm<-dim(dendat)[1]
   distat<-sqrt(pituus(dendat-t(matrix(refe,d,n))))
   infopointer<-seq(1,lkm)
}

distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
proba<-matrix(0,lkm,1)
ekamome<-matrix(0,lkm,d)
distcenter<-matrix(0,lkm,d)
branchradius<-matrix(0,lkm,1)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]
branchradius[node]<-radius[node]

if (type!="tail"){
  # volume calculation
  vol<-1
  k<-1
  ip<-infopointer[node]  #pcf$nodefinder[infopointer[node]]
  while (k<=d){
      vol<-vol*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
      k<-k+1
  }
  volume[node]<-vol
  ip2<-infopointer[node]
  proba[node]<-pcf$value[ip2]*vol

  # ekamome calculation
  newcente<-matrix(0,d,1)
  for (j in 1:d){
    volmin<-1
    k<-1
    while (k<=d){
       if (k!=j){
          volmin<-volmin*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
       }
       k<-k+1
    }
    ala<-pcf$support[2*j-1]+step[j]*pcf$down[ip,j]
    yla<-pcf$support[2*j-1]+step[j]*pcf$high[ip,j]
    newcente[j]<-volmin*(yla^2-ala^2)/2
  }
  ekamome[node,]<-newcente
  distcenter[node,]<-newcente/vol
}
else{
  volume[node]<-1
}

beg<-node             #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcf$down[note,i]   
  boundrec[node,2*i]<-pcf$high[note,i]  
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)  #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-pcf$down[note,i]  
         rec1[2*i]<-pcf$high[note,i] 
    }
    boundrec[node,]<-rec1

    # radius
    radius[node]<-distat[ord[node]]
    branchradius[node]<-radius[node]
    if (type!="tail"){
       # volume calculation
       vol<-1
       k<-1
       ip<-infopointer[node]    #pcf$nodefinder[infopointer[node]]
       while (k<=d){
          vol<-vol*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
          k<-k+1
       }
       volume[node]<-vol
       ip2<-infopointer[node]
       proba[node]<-pcf$value[ip2]*vol

       # ekamome calculation
       newcente<-matrix(0,d,1)
       for (jj in 1:d){
            volmin<-1
            k<-1
            while (k<=d){
               if (k!=jj){
                   volmin<-volmin*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
               }
               k<-k+1
            }
            ala<-pcf$support[2*jj-1]+step[jj]*pcf$down[ip,jj]
            yla<-pcf$support[2*jj-1]+step[jj]*pcf$high[ip,jj]
            newcente[jj]<-volmin*(yla^2-ala^2)/2
       }
       ekamome[node,]<-newcente
       distcenter[node,]<-newcente/vol
    }
    else{
       volume[node]<-1
    }

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        if (type=="tail") 
            istouch<-touchstep.tail(node,curroot,boundrec,child,sibling,
                                    infopointer,pcf$down,pcf$high) 
        else istouch<-touchstep(node,curroot,boundrec,child,sibling,
                           infopointer,pcf$down,pcf$high)
        if (istouch==1){

           # paivita parent, child, sibling, volume ekamome
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }
           if (type!="tail"){
              volume[node]<-volume[node]+volume[curroot]
              proba[node]<-proba[node]+proba[curroot]
              ekamome[node,]<-ekamome[node,]+ekamome[curroot,]
           }
           else{
              volume[node]<-volume[node]+volume[curroot]
           }

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])
           if (branchradius[node]<=branchradius[curroot]) 
                  distcenter[node,]<-distcenter[curroot,]
           branchradius[node]<-max(branchradius[node],branchradius[curroot])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]

        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter
# lf is the level set tree or the shape tree
if (type!="tail"){
   for (i in 1:lkm){
      for (j in 1:d){
          ekamome[i,j]<-ekamome[i,j]/volume[i]
      }
   }
   bary<-ekamome[root,]
}
if (type=="shape"){
  maxdis<-sqrt(distat[ord[length(ord)]])
  if (levmet=="proba")
     level<-taillevel(root,#child,sibling,
            parent,volume,proba)
  else 
     level<-sqrt(radius)
}
else{ #type="lst"
     level<-radius
     maxdis<-distat[ord[length(ord)]]
}
if (type=="tail"){
   center<-t(dendat[infopointer,])
}

if (type!="tail"){
  lf<-list(
  parent=parent,volume=volume,center=t(ekamome),level=level,
  root=root,
  #child=child,sibling=sibling,  #virhe??
  infopointer=infopointer,
  proba=proba,#radius=radius,
  #branchradius=sqrt(branchradius),
  distcenter=t(distcenter),
  refe=refe,maxdis=maxdis,bary=bary,lev=lev)
}
else{
  lf<-list(
  parent=parent,volume=volume,center=center,level=level,
  root=root,
  #child=child,sibling=sibling,  #virhe??
  infopointer=infopointer,
  #proba=proba,#radius=radius,
  #branchradius=sqrt(branchradius),
  #distcenter=t(distcenter),
  refe=refe,maxdis=maxdis,
  dendat=dendat)
}

# if ngrid given, reduce the lst
if (!is.null(ngrid)){
    stepsi<-maxdis/ngrid
    radii<-seq(0,maxdis,stepsi)
    lf<-treedisc(lf,pcf,r=radii,type=type)
}

return(lf)
}





leafsfirst<-function(pcf=NULL,lev=NULL,refe=NULL,type="lst",levmet="radius",
ordmet="etaisrec",ngrid=NULL,
dendat=NULL,rho=0,propor=NULL,lowest="dens",f=NULL)
{
# pcf is a piecewise constant object
# type= "lst"/"shape"
# levmet= "radius"/"proba"

if (lowest=="dens") lowest<-0 else lowest<-min(pcf$value)

if ((!is.null(lev)) || (!is.null(propor))){
    type<-"shape"
    if (!is.null(propor)) lev<-propor*max(pcf$value)
    if (is.null(refe)) refe<-locofmax(pcf)
}
if (!is.null(dendat)) type<-"tail"

if (type=="tail"){
   d<-dim(dendat)[2]
   pcf$high<-dendat
   pcf$down<-dendat
   if (is.null(refe)){
       refe<-matrix(0,1,d)
       for (i in 1:d) refe[1,i]<-mean(dendat[,i])
       refe<-refe[1:d]
   }
}
else{
  d<-length(pcf$N)
  step<-matrix(0,d,1)
  for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
}

if (type=="lst"){
  lkm<-length(pcf$value)
  distat<-pcf$value-lowest
  infopointer<-seq(1,lkm)     # links from nodes to recs
}
else if (type=="shape"){
  lenni<-length(pcf$value)
  distat<-matrix(0,lenni,1)
  infopointer<-matrix(0,lenni,1)
  lkm<-0
  for (i in 1:lenni){
    if (pcf$value[i]>=lev){
       lkm<-lkm+1
       nod<-i  #nod<-pcf$nodefinder[i]
       if (ordmet=="etaisrec"){
           recci<-matrix(0,2*d,1)
           for (jj in 1:d){
              recci[2*jj-1]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
              recci[2*jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
           }
           distat[lkm]<-etaisrec(refe,recci)
       }
       else{
          lowi<-matrix(0,d,1)
          uppi<-matrix(0,d,1)
          for (jj in 1:d){
             lowi[jj]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
             uppi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
          }
          baryc<-lowi+(uppi-lowi)/2  
          distat[lkm]<-etais(baryc,refe)
       }
       infopointer[lkm]<-i
    }
  }
}
else{  #type=="tail"
   d<-dim(dendat)[2]
   n<-dim(dendat)[1]
   lkm<-dim(dendat)[1]
   distat<-sqrt(pituus(dendat-t(matrix(refe,d,n))))
   infopointer<-seq(1,lkm)
}

distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   
#if (length(rho)==1) rho<-rep(rho,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
proba<-matrix(0,lkm,1)
ekamome<-matrix(0,lkm,d)
distcenter<-matrix(0,lkm,d)
branchradius<-matrix(0,lkm,1)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]
branchradius[node]<-radius[node]

if (type!="tail"){
  # volume calculation
  vol<-1
  k<-1
  ip<-infopointer[node]  #pcf$nodefinder[infopointer[node]]
  while (k<=d){
      vol<-vol*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
      k<-k+1
  }
  volume[node]<-vol
  ip2<-infopointer[node]
  proba[node]<-pcf$value[ip2]*vol

  # ekamome calculation
  newcente<-matrix(0,d,1)
  for (j in 1:d){
    volmin<-1
    k<-1
    while (k<=d){
       if (k!=j){
          volmin<-volmin*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
       }
       k<-k+1
    }
    ala<-pcf$support[2*j-1]+step[j]*pcf$down[ip,j]
    yla<-pcf$support[2*j-1]+step[j]*pcf$high[ip,j]
    newcente[j]<-volmin*(yla^2-ala^2)/2
  }
  ekamome[node,]<-newcente
  distcenter[node,]<-newcente/vol
}
else{  # type=tail
  if (is.null(f)) volume[node]<-1
  else{ 
       ip<-infopointer[node] 
       volume[node]<-1/(f[ip]*length(f))
  }
}

beg<-node             #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcf$down[note,i]   
  boundrec[node,2*i]<-pcf$high[note,i]  
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)  #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-pcf$down[note,i]  
         rec1[2*i]<-pcf$high[note,i] 
    }
    boundrec[node,]<-rec1

    # radius
    radius[node]<-distat[ord[node]]
    branchradius[node]<-radius[node]
    if (type!="tail"){
       # volume calculation
       vol<-1
       k<-1
       ip<-infopointer[node]    #pcf$nodefinder[infopointer[node]]
       while (k<=d){
          vol<-vol*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
          k<-k+1
       }
       volume[node]<-vol
       ip2<-infopointer[node]
       proba[node]<-pcf$value[ip2]*vol

       # ekamome calculation
       newcente<-matrix(0,d,1)
       for (jj in 1:d){
            volmin<-1
            k<-1
            while (k<=d){
               if (k!=jj){
                   volmin<-volmin*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
               }
               k<-k+1
            }
            ala<-pcf$support[2*jj-1]+step[jj]*pcf$down[ip,jj]
            yla<-pcf$support[2*jj-1]+step[jj]*pcf$high[ip,jj]
            newcente[jj]<-volmin*(yla^2-ala^2)/2
       }
       ekamome[node,]<-newcente
       distcenter[node,]<-newcente/vol
    }
    else{     #type==tail
       if (is.null(f)) volume[node]<-1
       else{ 
          ip<-infopointer[node] 
          volume[node]<-1/(f[ip]*length(f))
       }
    }

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        rhocur<-rho   #rho[infopointer[node]]  
        istouch<-touchstep(node,curroot,boundrec,child,sibling,
                           infopointer,pcf$down,pcf$high,rhocur)
        if (istouch==1){
{
           # paivita parent, child, sibling, volume ekamome
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }
           if (type!="tail"){
              volume[node]<-volume[node]+volume[curroot]
              proba[node]<-proba[node]+proba[curroot]
              ekamome[node,]<-ekamome[node,]+ekamome[curroot,]
           }
           else{  # type == tail
              volume[node]<-volume[node]+volume[curroot]
           }

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])
           if (branchradius[node]<=branchradius[curroot]) 
                  distcenter[node,]<-distcenter[curroot,]
           branchradius[node]<-max(branchradius[node],branchradius[curroot])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]
}
        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter
# lf is the level set tree or the shape tree
if (type!="tail"){
   for (i in 1:lkm){
      for (j in 1:d){
          ekamome[i,j]<-ekamome[i,j]/volume[i]
      }
   }
   bary<-ekamome[root,]
}
if (type=="shape"){
  maxdis<-sqrt(distat[ord[length(ord)]])
  if (levmet=="proba")
     level<-taillevel(root,#child,sibling,
            parent,volume,proba)
  else 
     level<-sqrt(radius)
}
else{ #type="lst"
     level<-radius+lowest
     maxdis<-distat[ord[length(ord)]]
}
if (type=="tail"){
   center<-t(dendat[infopointer,])
}

if (type!="tail"){
  lf<-list(
  parent=parent,volume=volume,center=t(ekamome),level=level,
  root=root,
  #child=child,sibling=sibling,  #virhe??
  infopointer=infopointer,
  proba=proba,#radius=radius,
  #branchradius=sqrt(branchradius),
  distcenter=t(distcenter),
  refe=refe,maxdis=maxdis,bary=bary,lev=lev)
}
else{
  lf<-list(
  parent=parent,volume=volume,center=center,level=level,
  root=root,
  #child=child,sibling=sibling,  #virhe??
  infopointer=infopointer,
  #proba=proba,#radius=radius,
  #branchradius=sqrt(branchradius),
  #distcenter=t(distcenter),
  refe=refe,maxdis=maxdis,
  dendat=dendat)
}

# if ngrid given, reduce the lst
if (!is.null(ngrid)){
    stepsi<-maxdis/ngrid
    radii<-seq(0,maxdis,stepsi)
    lf<-treedisc(lf,pcf,r=radii,type=type)
}

return(lf)
}





leafsfirst.shape<-function(pcf=NULL, lev=NULL, refe=NULL, levmet="radius",
ordmet="etaisrec", propor=NULL)
{
rho<-0

if (!is.null(propor)) lev<-propor*max(pcf$value)
if (is.null(refe)) refe<-locofmax(pcf)

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

  lenni<-length(pcf$value)
  distat<-matrix(0,lenni,1)
  infopointer<-matrix(0,lenni,1)
  lkm<-0
  for (i in 1:lenni){
    if (pcf$value[i]>=lev){
       lkm<-lkm+1
       nod<-i  #nod<-pcf$nodefinder[i]
       if (ordmet=="etaisrec"){
           recci<-matrix(0,2*d,1)
           for (jj in 1:d){
              recci[2*jj-1]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
              recci[2*jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
           }
           distat[lkm]<-etaisrec(refe,recci)
       }
       else{
          lowi<-matrix(0,d,1)
          uppi<-matrix(0,d,1)
          for (jj in 1:d){
             lowi[jj]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
             uppi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
          }
          baryc<-lowi+(uppi-lowi)/2  
          distat[lkm]<-etais(baryc,refe)
       }
       infopointer[lkm]<-i
    }
  }

distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   
#if (length(rho)==1) rho<-rep(rho,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
proba<-matrix(0,lkm,1)
ekamome<-matrix(0,lkm,d)
distcenter<-matrix(0,lkm,d)
branchradius<-matrix(0,lkm,1)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]
branchradius[node]<-radius[node]


  # volume calculation
  vol<-1
  k<-1
  ip<-infopointer[node]  #pcf$nodefinder[infopointer[node]]
  while (k<=d){
      vol<-vol*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
      k<-k+1
  }
  volume[node]<-vol
  ip2<-infopointer[node]
  proba[node]<-pcf$value[ip2]*vol

  # ekamome calculation
  newcente<-matrix(0,d,1)
  for (j in 1:d){
    volmin<-1
    k<-1
    while (k<=d){
       if (k!=j){
          volmin<-volmin*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
       }
       k<-k+1
    }
    ala<-pcf$support[2*j-1]+step[j]*pcf$down[ip,j]
    yla<-pcf$support[2*j-1]+step[j]*pcf$high[ip,j]
    newcente[j]<-volmin*(yla^2-ala^2)/2
  }
  ekamome[node,]<-newcente
  distcenter[node,]<-newcente/vol

beg<-node             #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcf$down[note,i]   
  boundrec[node,2*i]<-pcf$high[note,i]  
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)  #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-pcf$down[note,i]  
         rec1[2*i]<-pcf$high[note,i] 
    }
    boundrec[node,]<-rec1

    # radius
    radius[node]<-distat[ord[node]]
    branchradius[node]<-radius[node]

       # volume calculation
       vol<-1
       k<-1
       ip<-infopointer[node]    #pcf$nodefinder[infopointer[node]]
       while (k<=d){
          vol<-vol*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
          k<-k+1
       }
       volume[node]<-vol
       ip2<-infopointer[node]
       proba[node]<-pcf$value[ip2]*vol

       # ekamome calculation
       newcente<-matrix(0,d,1)
       for (jj in 1:d){
            volmin<-1
            k<-1
            while (k<=d){
               if (k!=jj){
                   volmin<-volmin*(pcf$high[ip,k]-pcf$down[ip,k])*step[k]
               }
               k<-k+1
            }
            ala<-pcf$support[2*jj-1]+step[jj]*pcf$down[ip,jj]
            yla<-pcf$support[2*jj-1]+step[jj]*pcf$high[ip,jj]
            newcente[jj]<-volmin*(yla^2-ala^2)/2
       }
       ekamome[node,]<-newcente
       distcenter[node,]<-newcente/vol

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
      
        istouch<-touchstep(node,curroot,boundrec,child,sibling,
                           infopointer,pcf$down,pcf$high,rho)
        if (istouch==1){
{
           # paivita parent, child, sibling, volume ekamome
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }

              volume[node]<-volume[node]+volume[curroot]
              proba[node]<-proba[node]+proba[curroot]
              ekamome[node,]<-ekamome[node,]+ekamome[curroot,]

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])
           if (branchradius[node]<=branchradius[curroot]) 
                  distcenter[node,]<-distcenter[curroot,]
           branchradius[node]<-max(branchradius[node],branchradius[curroot])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]
}
        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter
# lf is the level set tree or the shape tree

   for (i in 1:lkm){
      for (j in 1:d){
          ekamome[i,j]<-ekamome[i,j]/volume[i]
      }
   }
   bary<-ekamome[root,]

  maxdis<-sqrt(distat[ord[length(ord)]])
  if (levmet=="proba")
     level<-taillevel(root,#child,sibling,
            parent,volume,proba)
  else 
     level<-sqrt(radius)

  lf<-list(
  parent=parent,volume=volume,center=t(ekamome),level=level,
  root=root,
  #child=child,sibling=sibling,  #virhe??
  infopointer=infopointer,
  proba=proba,#radius=radius,
  #branchradius=sqrt(branchradius),
  distcenter=t(distcenter),
  refe=refe,maxdis=maxdis,bary=bary,lev=lev)

return(lf)
}

leafsfirst.tail<-function(dendat, rho=0, refe=NULL, dist.type="euclid")
{

n<-dim(dendat)[1]
d<-dim(dendat)[2]
pcfhigh<-dendat+rho
pcfdown<-dendat-rho
if (is.null(refe)){
      refe<-matrix(0,1,d)
      for (i in 1:d) refe[1,i]<-mean(dendat[,i])
      refe<-refe[1:d]
}

distat<-sqrt(pituus(dendat-t(matrix(refe,d,n))))
lkm<-n
infopointer<-seq(1,lkm)
if (length(rho)==1) rho<-rep(rho,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]

volume[node]<-1

beg<-node                 #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-pcfdown[note,i]   
  boundrec[node,2*i]<-pcfhigh[note,i]  
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)  #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-pcfdown[note,i]  
         rec1[2*i]<-pcfhigh[note,i] 
    }
    boundrec[node,]<-rec1

    # radius
    radius[node]<-distat[ord[node]]

    volume[node]<-1

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        #rhocur<-rho[infopointer[node]]  
        istouch<-touchstep.tail(node,curroot,boundrec,child,sibling,
                                infopointer,pcfdown,pcfhigh,rho,dendat,
                                dist.type=dist.type)

        if (istouch==1){

           # paivita parent, child, sibling, volume 
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }

           volume[node]<-volume[node]+volume[curroot]

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]

        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter

maxdis<-distat[ord[length(ord)]]
center<-t(dendat[infopointer,])

lf<-list(
parent=parent,volume=volume,center=center,level=radius,
root=root,
infopointer=infopointer,
refe=refe,maxdis=maxdis,
dendat=dendat)

return(lf)
}






leafsfirst.visu<-function(tt,pcf,lev=NULL,refe=NULL,type="lst",
levmet="radius",ordmet="etaisrec",
lkmbound=NULL,radius=NULL,
orde="furthest",suppo=T,propor=NULL,lty=NULL,numbers=TRUE,
sigcol="lightblue",cex.axis=1,cex=1)
{

if ((!is.null(lev)) || (!is.null(propor))){
    type<-"shape"
    if (is.null(refe)) refe<-locofmax(pcf)
    if (!is.null(propor)) lev<-propor*max(pcf$value)
}
if (is.null(refe)) refe<-locofmax(pcf)

pp<-plotprof(tt,plot=FALSE,data=TRUE)
vecs<-pp$vecs

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

# order the atoms for the level set with level "lev"

lenni<-length(pcf$value)
distat<-matrix(0,lenni,1)
infopointer<-matrix(0,lenni,1)

if (type=="lst"){
  lkm<-lenni
  distat<-pcf$value
  infopointer<-seq(1,lkm)
}
else{

lkm<-0
for (i in 1:lenni){
  if (pcf$value[i]>=lev){
     lkm<-lkm+1
     nod<-i  #nod<-pcf$nodefinder[i]
     if (ordmet=="etaisrec"){
         recci<-matrix(0,2*d,1)
         for (jj in 1:d){
            recci[2*jj-1]<-pcf$support[2*jj-1]+step[jj]*pcf$down[nod,jj]
            recci[2*jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
         }
         distat[lkm]<-etaisrec(refe,recci)
     }
     else{
         lowi<-matrix(0,d,1)
         uppi<-matrix(0,d,1)
         for (jj in 1:d){
            lowi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$down[nod,jj]
            uppi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
         }
         baryc<-lowi+(uppi-lowi)/2
         distat[lkm]<-etais(baryc,refe)  #etais(baryc[lk m,],baryind)
     }
     infopointer[lkm]<-i
  }
}

}  #else

distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   #pointe->pcf$value,pcf$nodefinder

ord<-order(distat)
infopointer<-infopointer[ord]

if (suppo){
  xmin<-pcf$support[1]
  xmax<-pcf$support[2]
  ymin<-pcf$support[3]
  ymax<-pcf$support[4]
}
else{
  xmin<-tt$refe[1]-tt$maxdis  #pcf$support[1]
  xmax<-tt$refe[1]+tt$maxdis  #pcf$support[2]
  ymin<-tt$refe[1]-tt$maxdis  #pcf$support[3]
  ymax<-tt$refe[2]+tt$maxdis  #pcf$support[4]
}

plot(x=refe[1],y=refe[2],xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
pch=20,cex.axis=cex.axis) #,col="red")

i<-1
while (i<=lkm){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]   #ip<-pcf$nodefinder[infopointer[node]]

     x1<-pcf$support[1]+step[1]*pcf$down[ip,1]
     x2<-pcf$support[1]+step[1]*pcf$high[ip,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[ip,2]
     y2<-pcf$support[3]+step[2]*pcf$high[ip,2] 
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="gray",lty=lty)

     i<-i+1
}

if (!is.null(lkmbound)){
  i<-1
  while (i<=lkmbound){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]  #ip<-pcf$nodefinder[infopointer[node]]

     x1<-pcf$support[1]+step[1]*pcf$down[ip,1]
     x2<-pcf$support[1]+step[1]*pcf$high[ip,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[ip,2]
     y2<-pcf$support[3]+step[2]*pcf$high[ip,2] 
     dev.set(which = dev.next())
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=sigcol,lty=lty)
     #points(x=refe[1],y=refe[2],pch=20,col="red")
     if (numbers) text(x=x1+(x2-x1)/2,y=y1+(y2-y1)/2,paste(i),cex=cex)

     i<-i+1
  }
}
else{
  i<-1
  radu<-tt$level[lkm]  #tt$madxdis
  while (radu>=radius){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]  #ip<-pcf$nodefinder[infopointer[node]]

     x1<-pcf$support[1]+step[1]*pcf$down[ip,1]
     x2<-pcf$support[1]+step[1]*pcf$high[ip,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[ip,2]
     y2<-pcf$support[3]+step[2]*pcf$high[ip,2] 
     dev.set(which = dev.next())
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="blue",lty=lty)
     points(x=refe[1],y=refe[2],pch=20,col="red")

     i<-i+1
     radu<-tt$level[node]
  }
}

}

leikkaa<-function(rec1,rec2){
#Makes an intersection of rectangles rec1, rec2
#rec1,rec2 are 2*d vectors
#
#Returns 2*d-vector or NA if intersection is empty
#
d<-length(rec1)/2
tulos<-matrix(0,2*d,1)
i<-1
while ((i<=d) && (!is.na(tulos))){  
    tulos[2*i-1]<-max(rec1[2*i-1],rec2[2*i-1])
    tulos[2*i]<-min(rec1[2*i],rec2[2*i])
    if (tulos[2*i]<=tulos[2*i-1]) tulos<-NA
    i<-i+1
}
return(tulos)
}

levord<-function(beg,sibling,sibord,centers,crit){
#order at the given level
#
# find first
#
itemnum<-length(sibling)
diffe<-matrix(NA,itemnum,1)    #NA is infty
cur<-beg
curre<-centers[,cur]
diffe[cur]<-etais(curre,crit)
sibnum<-1   #if beg has no siblings, then sibnum=1 (beg is its own sibling)
while (sibling[cur]>0){
     cur<-sibling[cur] 
     curre<-centers[,cur]
     diffe[cur]<-etais(curre,crit)
     sibnum<-sibnum+1
}
first<-omaind(-diffe)
#sibord[first]<-1
#
# find distances to first
#
firstcenter<-centers[,first]
distofir<-matrix(NA,itemnum,1)
cur<-beg
curre<-centers[,cur]
distofir[cur]<-etais(curre,firstcenter)
while (sibling[cur]>0){
    cur<-sibling[cur]
    curre<-centers[,cur]
    distofir[cur]<-etais(curre,firstcenter)
}
#  
# fill sibord in the order of closest to first 
#
ind<-1
remain<-sibnum
while (remain>0){  
     nex<-omaind(distofir)
     sibord[nex]<-ind
     distofir[nex]<-NA 
     ind<-ind+1
     remain<-remain-1     
}
return(sibord)
}











liketree<-function(dendat,pcf,lst)
{

# "lst$infopointer" gives links from nodes to recs
# invert the links
rnum<-length(pcf$value)
nodefinder<-matrix(0,rnum,1)
for (i in 1:rnum) nodefinder[lst$infopointer[i]]<-i

n<-dim(dendat)[1]
d<-dim(dendat)[2]

step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

# find links from dendat to pcf
# (for simplicity we delete multiple observations in a bin
den2pcf<-matrix(0,n,1)
pcf2den<-matrix(0,rnum,1)
varauslista<-matrix(0,rnum,1)
dendat2<-matrix(0,n,d)
n2<-0
for (i in 1:n){
    j<-1
    notjet<-TRUE
    while ((j<=rnum) && (notjet)){
         inside<-TRUE
         coordi<-1
         while ((inside) && (coordi<=d)){
             ala<-pcf$down[j,coordi]
             yla<-pcf$high[j,coordi]
             ala<-pcf$support[2*coordi-1]+ala*step[coordi]
             yla<-pcf$support[2*coordi-1]+yla*step[coordi]
             if ((dendat[i,coordi]<ala) || (dendat[i,coordi]>yla)) 
                         inside<-FALSE
             coordi<-coordi+1
         }
         if (inside){
              notjet<-FALSE
              if (varauslista[j]==0){ 
                  varauslista[j]<-1
                  n2<-n2+1
                  dendat2[n2,]<-dendat[i,]
                  den2pcf[n2]<-j
                  pcf2den[j]<-n2
              }
         }
         j<-j+1
    }
}
dendat2<-dendat2[1:n2,]

# make tree
parent<-matrix(0,n2,1)
center<-matrix(0,d,n2)
level<-matrix(0,n2,1)
for (i in 1:n2){
   rec<-den2pcf[i]
   node<-nodefinder[rec]
   level[i]<-lst$level[node]

   obs<-0
   curnode<-node
   notfound<-TRUE
   while ((notfound) && (lst$parent[curnode]>0)){
         curnode<-lst$parent[curnode] 
         rec<-lst$infopointer[curnode]
         if (pcf2den[rec]>0){ 
              notfound<-FALSE
              obs<-pcf2den[rec]
         }
   }
   parent[i]<-obs
}
center<-t(dendat2)

return(list(parent=parent,center=center,level=level,
dendat=dendat2,infopoint=seq(1:n2)))
}

listchange<-function(AtomlistAtom,AtomlistNext,totbegSepary,
begsSepaNext,begsSepaBegs,atomsSepaNext,atomsSepaAtom,
terminalnum,beg){
#
#create begs: beginnings of lists of atoms
#beg is index to AtomlistAtom/Next
#totbegsepary is index to begsSepaBegs/Next
#
begs<-matrix(0,terminalnum,1)
#
runnerBegs<-totbegSepary  #total beginning of list is at the root of the tree
runnerOrigi<-beg
runnerOrigiprev<-beg
sepalkm<-0
while (runnerBegs>0){
   sepalkm<-sepalkm+1
   runnerAtoms<-begsSepaBegs[runnerBegs]
   begs[sepalkm]<-runnerOrigi
   #  first step (in order to get also runnerOrigiprev to play)
   AtomlistAtom[runnerOrigi]<-atomsSepaAtom[runnerAtoms]
   runnerOrigiprev<-runnerOrigi
   runnerOrigi<-AtomlistNext[runnerOrigi]
   runnerAtoms<-atomsSepaNext[runnerAtoms]
   while (runnerAtoms>0){
       AtomlistAtom[runnerOrigi]<-atomsSepaAtom[runnerAtoms]
       runnerOrigiprev<-runnerOrigi
       runnerOrigi<-AtomlistNext[runnerOrigi]
       runnerAtoms<-atomsSepaNext[runnerAtoms]
   }
   AtomlistNext[runnerOrigiprev]<-0      #mark the end of the list
   runnerBegs<-begsSepaNext[runnerBegs]
}
#
begs<-begs[1:sepalkm]
#
return(list(begs=begs,AtomlistAtom=AtomlistAtom,AtomlistNext=AtomlistNext))
}
locofmax<-function(pcf)
{
d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]
}

nod<-which.max(pcf$value)

lowi<-matrix(0,d,1)
uppi<-matrix(0,d,1)
for (jj in 1:d){
    lowi[jj]<-pcf$support[2*jj-1]+step[jj]*(pcf$down[nod,jj])
    uppi[jj]<-pcf$support[2*jj-1]+step[jj]*pcf$high[nod,jj]
}
baryc<-lowi+(uppi-lowi)/2  

return(baryc)
}

lst2xy<-function(lst,type="radius",gnum=1000)
{
# gives the x and y vectors of a volume transform

if (type=="radius") pv<-plotvolu(lst,data=T,toplot=F)
else{
   lst2<-lst
   lst2$volume<-lst$proba
   pv<-plotvolu(lst2,data=T,toplot=F)
}

lenni<-length(pv$xcoor)/2
xs<-t(matrix(pv$xcoor,2,lenni))
ys<-matrix(0,lenni,1)
for (i in 1:lenni) ys[i]<-pv$ycoor[2*i]

or<-order(ys)
xs<-xs[or,]
ys<-ys[or]

xlow<-min(xs)
xhig<-max(xs)
xstep<-(xhig-xlow)/gnum
x<-seq(xlow,xhig,xstep)
y<-matrix(0,length(x),1)

x[1]<-xs[1,1]
x[length(x)]<-xs[1,2]
i<-2
while (i <= lenni){
  lowind<-round(length(x)*(xs[i,1]-xlow)/(xhig-xlow))
  higind<-round(length(x)*(xs[i,2]-xlow)/(xhig-xlow))
  y[higind:lowind]<-ys[i]
  i<-i+1
}

return(list(x=x,y=y))
}



lstseq.kern<-function(dendat,hseq,N,lstree=NULL,level=NULL,
Q=NULL,kernel="gauss",hw=NULL,algo="leafsfirst",support=NULL)
{
hnum<-length(hseq)
if ((hnum>1) && (hseq[1]<hseq[2])) hseq<-hseq[seq(hnum,1)]

if (algo=="leafsfirst"){

  for (i in 1:hnum){   
      h<-hseq[i]
      pcf<-pcf.kern(dendat,h,N,kernel=kernel,support=support)
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

}
else{  #algo=="decomdyna"
  lstseq<-profkern(dendat,hseq,N,Q,kernel=kernel,hw=hw)
}

if (is.null(lstree)) lstseq<-NULL
if (is.null(level)) stseq<-NULL
return(list(lstseq=lstseq,pcfseq=pcfseq,stseq=stseq,hseq=hseq,type="kernel"))
}

makehis<-function(regdat)
{
xlkm<-length(regdat$hila[,1])  #muuttujien lkm
valipit<-matrix(0,1,xlkm)
i<-1
while (i<=xlkm){
  if (regdat$hila[i,1]>1) 
    valipit[i]<-(regdat$hila[i,3]-regdat$hila[i,2])/(regdat$hila[i,1]-1)
  i<-i+1
}
lnum<-length(regdat$ind[,1])   #length(regdat$dep)
items<-matrix(0,lnum,2*xlkm)
arvot<-matrix(0,lnum,1)
i<-1
while (i<=lnum){
  arvot[i]<-regdat$dep[i]  
  j<-1
  while (j<=xlkm){
    items[i,2*j-1]<-regdat$ind[i,j]-valipit[j]/2
    items[i,2*j]<-regdat$ind[i,j]+valipit[j]/2
    j<-j+1
  }
  i<-i+1
}
return(list(values=arvot,recs=items))
}


makeinfo<-function(left,right,mean,low,upp)
{
lehdet<-findleafs(left,right)

d<-dim(low)[2]
nodenum<-length(lehdet)         #length(left)

value<-matrix(0,nodenum,1)
infolow<-matrix(0,nodenum,d)
infoupp<-matrix(0,nodenum,d)
nodefinder<-matrix(0,nodenum,1)
infopointer<-matrix(0,nodenum,1)

runner<-1
leafnum<-0
while (runner<=nodenum){
  if ((!is.na(lehdet[runner])) && (lehdet[runner]==1) && (mean[runner]>0)){  
      # we are in leaf where the value is positive
      leafnum<-leafnum+1
      value[leafnum]<-mean[runner]
      nodefinder[leafnum]<-runner
      infolow[leafnum,]<-low[runner,]
      infoupp[leafnum,]<-upp[runner,]

      infopointer[runner]<-leafnum
  }
  runner<-runner+1
}
value<-value[1:leafnum]
nodefinder<-nodefinder[1:leafnum]
infolow<-infolow[1:leafnum,]
infoupp<-infoupp[1:leafnum,]

return(list(value=value,low=infolow,upp=infoupp,nodefinder=nodefinder,
infopointer=infopointer,
terminalnum=leafnum))
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
massat<-function(rec){
#Calculates a vector of masses of a set of rectangles
#
#rec is k*(2*d)-matrix, represents k rectangles in d-space
#Returns a k-vector
#
#if (dim(t(rec))[1]==1) k<-1 else k<-length(rec[,1])  #rows of rec
if (dim(t(rec))[1]==1){
 d<-length(rec)/2
 vol<-1
 j<-1
   while ((j<=d) && (vol>0)){
     if (rec[2*j]<=rec[2*j-1]) vol<-0
     else vol<-vol*(rec[2*j]-rec[2*j-1])
     j<-j+1
   }
  tulos<-vol
}
else{
 k<-length(rec[,1])
 d<-length(rec[1,])/2
 tulos<-matrix(0,k,1)
 for (i in 1:k){
   vol<-1
   j<-1
   while ((j<=d) && (vol>0)){
     if (rec[i,2*j]<=rec[i,2*j-1]) vol<-0
     else vol<-vol*(rec[i,2*j]-rec[i,2*j-1])
     j<-j+1
   }
   tulos[i]<-vol
 }
}
return(tulos)
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

maxnodenum<-function(dendat,h,N,n,d)
{
minim<-matrix(0,d,1)
maxim<-matrix(0,d,1)
i<-1
while (i<=d){
  minim[i]<-min(dendat[,i])  
  maxim[i]<-max(dendat[,i])
  i<-i+1;
}
hmax<-max(h)
delta<-(maxim-minim+2*hmax)/(N+1)
mindelta<-min(delta)
maxpositive<-ceiling(n*(2*hmax/mindelta)^d)
bigd<-sum(log(N,base=2))
maxnode<-ceiling(bigd*maxpositive)

return(list(maxnode=maxnode,maxpositive=maxpositive));
}
modecent<-function(lst){
#
parents<-lst$parent
levels<-lst$level
volumes<-lst$volume   
centers<-lst$center
d<-dim(centers)[1]            #d<-length(centers[,1])
#
mlkm<-moodilkm(parents)
modloc<-mlkm$modloc
lkm<-mlkm$lkm
#
mut<-multitree(parents)
roots<-mut$roots
child<-mut$child
sibling<-mut$sibling 
#
crit<-rep(0,d)               #order so that 1st closest to origo
sibord<-siborder(mut,crit,centers)   
#
itemnum<-length(parents)
vecs<-matrix(NA,itemnum,4)
vecs<-alloroot(vecs,roots,sibord,levels,volumes)
vecs<-plotdata(roots,child,sibling,sibord,levels,volumes,vecs) 
#
res<-matrix(0,lkm,d)
#
for (i in 1:lkm){
   sija<-modloc[i]
   res[i,]<-centers[,sija]
}
#
ord<-vecs[,1]   #in this order we want modes
ordpick<-matrix(0,lkm,1)
for (i in 1:lkm){
  sija<-modloc[i]
  ordpick[i]<-ord[sija]
}
#
pointer<-seq(1:lkm)
pointer<-omaord2(pointer,ordpick) #pointer on the order determined by ord
#
endres<-res
for (i in 1:lkm){
  sija<-pointer[i]
  endres[i,]<-res[sija,]
}
#
return(endres)
}

modegraph<-function(estiseq,hseq=NULL,paletti=NULL)  #,reverse=F)
{
# we want that the largest h is first (1 mode, oversmoothing)

if (is.null(hseq))
   if (!is.null(estiseq$type)){
       if (estiseq$type=="greedy") hseq<--estiseq$hseq
       if (estiseq$type=="bagghisto") hseq<--estiseq$hseq
       if (estiseq$type=="carthisto")  hseq<--estiseq$leaf
       if (estiseq$type=="kernel")  hseq<-estiseq$hseq    
   }
   else hseq<-estiseq$hseq

hnum<-length(hseq)

treelist<-estiseq$lstseq
d<-dim(treelist[[1]]$center)[1]

if (hseq[1]<hseq[2]){   #(reverse){  
    #if ((hnum>1) && (is.null(hseq))) 
    hseq<-hseq[seq(hnum,1)]
    apuseq<-list(treelist[[hnum]])
    i<-2
    while (i <= hnum){
         apuseq<-c(apuseq,list(treelist[[hnum-i+1]]))
         i<-i+1 
   }
   treelist<-apuseq
}

if (is.null(paletti))
paletti<-c("red","blue","green","turquoise","orange","navy",
"darkgreen","orchid",colors()[50:100])

low<-matrix(0,hnum,1)
upp<-matrix(0,hnum,1)
tot<-moodilkm(treelist[[1]]$parent)$lkm  #tot is the number of modes over all lst:s
low[1]<-1
upp[1]<-tot
i<-2
while (i <= hnum){
  lkmm<-moodilkm(treelist[[i]]$parent)$lkm
  tot<-tot+lkmm
  low[i]<-upp[i-1]+1
  upp[i]<-low[i]+lkmm-1
  i<-i+1
}

xcoor<-matrix(0,tot,d)
ycoor<-matrix(0,tot,1)
parent<-matrix(0,tot,1)
mlabel<-matrix(0,tot,1)
flucpoints<-matrix(0,hnum,1)
nodepointer<-matrix(0,tot,1)
colot<-matrix("",tot,1)

# first we allocate colors for the largest h
colrun<-1  #low[1]
while (colrun<=upp[1]){
   colot[colrun]<-paletti[colrun]
   colrun<-colrun+1
}

laskuri<-1
srun<-1
mlkmpre<-1
flucnum<-0
while (srun<=hnum){  
    mlkm<-moodilkm(treelist[[srun]]$parent)
    if (mlkmpre < mlkm$lkm){
          flucnum<-flucnum+1
          flucpoints[flucnum]<-srun
    }

    for (j in 1:mlkm$lkm){
        loca<-mlkm$modloc[j]
        if (d>1){
           for (dim in 1:d){
              xcoor[laskuri,dim]<-treelist[[srun]]$center[dim,loca]
           }
        }
        else{
              xcoor[laskuri]<-treelist[[srun]]$center[loca]
        }
        ycoor[laskuri]<-hseq[srun]
        mlabel[laskuri]<-j
        nodepointer[laskuri]<-loca

        laskuri<-laskuri+1
    }

    if (srun>1){

       vec1<-matrix(0,mlkmpre,d)
       vec2<-matrix(0,mlkm$lkm,d)
       vec1[1:mlkmpre,]<-xcoor[low[srun-1]:upp[srun-1],]
       vec2[1:mlkm$lkm,]<-xcoor[low[srun]:upp[srun],]
       vm<-vectomatch(vec1,vec2)
       for (jj in low[srun]:upp[srun]){
           parent[jj]<-vm$parent[jj-low[srun]+1]+low[srun-1]-1
           if (vm$newnode[jj-low[srun]+1]==1){ 
                colot[jj]<-paletti[colrun]
                colrun<-colrun+1
           }
           else colot[jj]<-colot[parent[jj]]
      }
    }

    mlkmpre<-mlkm$lkm
    srun<-srun+1 
}

xcoor<-xcoor[1:(laskuri-1),]
ycoor<-ycoor[1:(laskuri-1)]
parent<-parent[1:(laskuri-1)]
colot<-colot[1:(laskuri-1)]
mlabel<-mlabel[1:(laskuri-1)]
nodepointer<-nodepointer[1:(laskuri-1)]
flucpoints<-flucpoints[1:flucnum]

mt<-list(xcoor=xcoor,ycoor=t(ycoor),
parent=parent,colot=colot,hseq=hseq,type=estiseq$type,
upp=upp,low=low,
mlabel=t(mlabel),
flucpoints=t(flucpoints),
nodepointer=t(nodepointer)
)

return(mt)
}




modetestgauss<-function(lst,n)
{

len<-length(lst$parent)
testvec<-matrix(0,len,1)    #this is output

em<-excmas(lst)

for (i in 1:len){

   if (lst$parent[i]!=0) val<-lst$level[lst$parent[i]]
   else val<-0

   a<-sqrt(n)*em[i]/sqrt(val*lst$volume[i])
   testvec[i]<-2*(1-pnorm(a))
}

return(testvec)
}
modetest<-function(pk,pknum,
h=NULL,N=NULL,Q=NULL,bootnum=NULL,delta=NULL,nsimu=NULL,minim=NULL,
type="boots",kernel="gauss",
n=NULL)
{

#pk is a list of level set trees
#h is vector of smoothing parameter values
#M is the number of bootstrap samples to be generated

run<-1
while (run<=pknum){
   curlst<-pk[[run]]

   if (type=="boots"){
       curh<-h[run]
       curmotes<-modetestydin(curlst,curh,N,Q,bootnum,delta,nsimu,minim,kernel)
   }
   else{
       curmotes<-modetestgauss(curlst,n)
   }

   if (run==1){
      if (pknum==1){
          moteslist<-curmotes
      }
      else{
          moteslist=list(curmotes)
      }
   }
   else{
      moteslist=c(moteslist,list(curmotes))
   }
   run<-run+1
}
#
return(moteslist)
}



modetestydin<-function(lstree,h,N,Q,bootnum,delta,nsimu,minim,kernel){

#we will approximate the estimate with a function whose values are
#levels of level set tree, this estimate has already been 
#normalized to integrate to one

index<-lstree$index
etnum<-dim(index)[1]
d<-length(N)

len<-length(lstree$parent)
mut<-multitree(lstree$parent)
mt<-pruneprof(mut)

#branchnodes<-findbranchMT(mt,len)
branchn<-findbranch(lstree$parent)$indicator
bnumbeg<-length(branchn)
branchnodes<-matrix(0,bnumbeg,1)
bnum<-0
for (i in 1:bnumbeg){
  if (branchn[i]==1){
      bnum<-bnum+1
      branchnodes[bnum]<-i
  }
}
branchnodes<-branchnodes[1:bnum]

testvec<-matrix(0,len,1)    #this is output

i<-1
while (i<=bnum){
   brnode<-branchnodes[i]
   
   #first cut the level set tree
   #3 cases: 1. brnode is linked from parent
   #         2. brnode is linked from previous sibling

   newchild<-mut$child
   newsibling<-mut$sibling
   newroots<-mut$roots

   brpare<-lstree$parent[brnode]
  
   if (brpare>0){  # brnode is not a root
 
   etsi<-mut$child[brpare]
{  if (etsi==brnode){
       newchild[brpare]<-mut$sibling[etsi]
   }
   else{
      while (etsi!=brnode){
         prevetsi<-etsi
         etsi<-mut$sibling[etsi]
      }
      newsibling[prevetsi]<-mut$sibling[etsi]
   }
} 
  
   #normalize the estimate to integrate to one
   
   kerroin<-cintemul(newroots,newchild,newsibling,lstree$volume,lstree$level)
   newlevel<-lstree$level/kerroin
   
   #creat value-vector "newvalue" with cutted values
   #newvalue*volofatom is probablility vector
   
   newvalue<-cutvalue(newroots,newchild,newsibling,
             newlevel,lstree$component,
             lstree$AtomlistAtom,lstree$AtomlistNext,etnum)
   
   #calculate the test statistics

         tstat<-excmas(lstree)[brnode]
          # cuttedlevel<-lstree$level[brpare]
          #cintemul(brnode,mut$child,mut$sibling,
          #lstree$volume,lstree$level-cuttedlevel)
   
   #generate samples from the cutted estimate
   
   overfrekv<-0
   j<-1
   while (j<=bootnum){
      dendatj<-simukern(nsimu,d,seed=j,newvalue,index,delta,minim,h)
      pk<-profkern(dendatj,h,N,Q,compoinfo=T,kernel=kernel)

      #find modes which are in the set associated with node brnode

      mlkm<-moodilkm(pk$parent)
      setissa<-matrix(0,mlkm$lkm,1)
      setissalkm<-0
      for (mrun in 1:mlkm$lkm){
           mloc<-mlkm$modloc[mrun]
           kandi<-pk$center[,mloc]
           torf<-onsetissa(kandi,h,delta,minim,
                           brnode,lstree$component,
                           lstree$index,
                           lstree$AtomlistAtom,lstree$AtomlistNext)
           if (torf){
               setissalkm<-setissalkm+1
               setissa[setissalkm]<-mloc 
           }
      }
      if (setissalkm>0){ 
        
         setissa<-setissa[1:setissalkm]

      #calculate the excess mass over "cuttedlevel" for
      #the modes in setissa 
      #Note that there may be a branching after "cuttedlevel"
      #We find the multitree for pk, then we cut this tree
      #by choosing the root node to be those nodes which are
      #arrived from modes (stop when the cuttedlevel is passed)
      
      cuttedlevel<-lstree$level[brpare]
      pkmut<-multitree(pk$parent)
      pkroots<-matrix(0,setissalkm,1)
      pkrootslkm<-0
      for (mrun in 1:setissalkm){
         node<-setissa[mrun]
         while ((pk$level[node]>cuttedlevel) && (node>0)){
              node<-pk$parent[node]
         }
         if (node>0){
             pkrootslkm<-pkrootslkm+1
             pkroots[pkrootslkm]<-node
         }
      }
      if (pkrootslkm>0){
         pkroots<-pkroots[1:pkrootslkm]

         bootstat<-cintemul(pkroots,pkmut$child,pkmut$sibling,
                            pk$volume,pk$level-cuttedlevel)
      }
      else{    #setissalkm>0, pkrootslkm==0
         bootstat<-0
      }
      
      }
      else{    #setissalkm==0
          bootstat<-0
      }

      if (bootstat<tstat){
          overfrekv<-overfrekv+1
      }
      j<-j+1
   }
   testvec[brnode]<-overfrekv/bootnum    #p-values are returned

   }

   i<-i+1
}

return(testvec)
}













montecarlo.ball<-function(dendat,rho,M,seed=1,type="ball")
{
# dendat on n*d matriisi
n<-dim(dendat)[1]
d<-dim(dendat)[2]

if (type=="ball"){
   keski<-colMeans(dendat)
   etais<-matrix(0,n,1)
   for (i in 1:n) etais[i]<-sqrt(sum((keski-dendat[i,])^2))
   masi<-max(etais)
   sade<-masi+rho
   set.seed(seed)
   polap<-sade*sqrt(runif(M))
   polax<-matrix(rnorm(M*d),M,d)
   varia<-matrix(0,M,d)
   for (i in 1:M) varia[i,]<-keski+polap[i]*polax[i,]/sqrt(sum(polax[i,]^2))
}
else {   # type=rectangular
  lows<-matrix(0,d,1)
  higs<-matrix(0,d,1)
  for (i in 1:d){
      ma<-max(dendat[,i])
      mi<-min(dendat[,i])
      lows[i]<-mi-rho
      higs[i]<-ma+rho
  }
  set.seed(seed)
  varia<-matrix(runif(M*d),M,d)
  for (i in 1:d) varia[,i]<-varia[,i]*(higs[i]-lows[i])+lows[i]
}

count<-matrix(0,M,1)
for (i in 1:M){
    point<-varia[i,]
    sisalla<-0
    j<-1
    while ((j<=n)&&(sisalla==0)){
        dista2<-sum((point-dendat[j,])^2)
        if (dista2<=rho^2){ 
            sisalla<-1
            count[i]<-1
        }
        j<-j+1
    }
}

if (type=="ball"){
  voluball<-pi*sade^2
  volu<-voluball*sum(count)/M
}
else{
  volurec<-1
  for (i in 1:d) volurec<-volurec*(higs[i]-lows[i])
  volu<-volurec*sum(count)/M
}

return(volu)
}


montecarlo.complex<-function(dendat,complex,rho,M,seed=1)
{
# dendat on n*d matriisi
n<-dim(dendat)[1]
d<-dim(dendat)[2]
lkm<-dim(complex)[1]

# create Monte Carlo sample
  ota<-c(complex)
  dendat2<-dendat[ota,]
  lows<-matrix(0,d,1)
  higs<-matrix(0,d,1)
  for (i in 1:d){
      ma<-max(dendat2[,i])
      mi<-min(dendat2[,i])
      lows[i]<-mi
      higs[i]<-ma
  }
  set.seed(seed)
  varia<-matrix(runif(M*d),M,d)
  for (i in 1:d) varia[,i]<-varia[,i]*(higs[i]-lows[i])+lows[i]

# laske kuinka monta joukossa

count<-matrix(0,M,1)
for (i in 1:M){
    point<-varia[i,]
    sisalla<-0
    j<-1
    while ( (j<=lkm) && (sisalla==0) ){
        simpl<-complex[j,]
        simple<-dendat[simpl,]
        sisalla2<-is.inside.simp.bary(point,simple)
        #sisalla2<-is.inside.simp.long(point,simple,rho)
        #sisalla2<-is.inside.simp(point,simple,rho)
        if (sisalla2==1){ 
               sisalla<-1
               count[i]<-1
        }
        j<-j+1
    }
}

# lakse tilavuus

  volurec<-1
  for (i in 1:d) volurec<-volurec*(higs[i]-lows[i])
  volu<-volurec*sum(count)/M

return(volu)
}

moodilkm<-function(vanhat)
{
#Lasketaan moodien lukumaara tiheyspuusta.
#Tiheyspuusta kaytettavissa vektori vanhat.
#Mikali solmu ei ole minkaan solmun vanhempi, se on lehti.

pit<-length(vanhat)
leima<-matrix(0,pit,1)
i<-1
while (i<=pit){
  solmu<-vanhat[i]
  leima[solmu]<-1
  i<-i+1
}
eimoodi<-sum(leima)
lkm<-(pit-eimoodi)
ykk<-rep(1,pit)
modnodes<-ykk-leima
#
moodiloc<-matrix(0,lkm,1)
ind<-1
for (i in 1:pit){
  if (modnodes[i]==1){
       moodiloc[ind]<-i
       ind<-ind+1
  }
}
#
return(list(lkm=lkm,modnodes=t(modnodes),modloc=t(moodiloc)))
}


mtest<-function(profile,n){
#
parents<-profile$parent
volumes<-profile$volume
levels<-profile$level
#
nodelkm<-length(parents)
lowmasses<-matrix(1,nodelkm,1)
for (i in 1:nodelkm){
   par<-parents[i]
   if (par==0){
      lowmasses[i]<-levels[i]*volumes[i]
   }
   else{
      lowmasses[i]<-levels[par]*volumes[i]
   }
}
testcrit<-sqrt(lowmasses/n)
#
return(t(testcrit))
#return(t(lowmasses))
}






multitree<-function(parents)
{
#Makes sibling-links and child-links

itemnum<-length(parents)
sibling<-matrix(0,itemnum,1)
child<-matrix(0,itemnum,1)
roots<-matrix(0,itemnum,1)
siborder<-matrix(0,itemnum,1)

rootnum<-0
for (i in itemnum:1){
  par<-parents[i]
  if (par==0){   #i is root (does not have parent)
     rootnum<-rootnum+1
     roots[rootnum]<-i
     siborder[i]<-rootnum
  } 
  else{          #i has parent
      if (child[par]==0){  #no childs so far
        child[par]<-i
        siborder[i]<-1
      }
      else{    #go to the end of sibling list
        chi<-child[par]
        sibsi<-1
        while(sibling[chi]>0){
           chi<-sibling[chi]
           sibsi<-sibsi+1
        }
        sibling[chi]<-i    #put to the sibling list
        siborder[i]<-sibsi+1
      }
  }
}
roots<-roots[1:rootnum]
return(list(child=child,sibling=sibling,roots=roots,siborder=siborder))
}
negapart<-function(pcf)
{
pcf$value<--pmin(pcf$value,0)
return(pcf)
}
nn.indit<-function(dendat)
{
n<-dim(dendat)[1]
maxk<-n-1
indmat<-matrix(0,n,maxk)

eta<-dist(dendat)
#i<j eta[n*(i-1) - i*(i-1)/2 + j-i]

for (i in 2:(n-1)){
   i1<-seq(1,i-1)
   j1<-i
   irow1<-eta[n*(i1-1) - i1*(i1-1)/2 + j1-i1]
   j2<-seq(i+1,n)
   irow2<-eta[n*(i-1) - i*(i-1)/2 + j2-i]
   irow<-c(irow1,irow2)
   or<-order(irow)
   poisi<-c(seq(1,i-1),seq(i+1,n))
   indmat[i,]<-poisi[or]
}

i<-1
j<-seq(i+1,n)
irow<-eta[n*(i-1) - i*(i-1)/2 + j-i]
or<-order(irow)
poisi<-seq(2,n)
indmat[i,]<-poisi[or]

i<-n
i1<-seq(1,n-1)
j<-i
irow<-eta[n*(i1-1) - i1*(i1-1)/2 + j-i1]
or<-order(irow)
poisi<-seq(1,n-1)
indmat[i,]<-poisi[or]

return(indmat)
}
nn.likeset<-function(dendat,radmat,k,p=0.1,lambda=NULL)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

volunitball<-volball(1,d)

radit<-radmat[,k]
evat<-k/(n*radit^d*volunitball)
if (is.null(lambda)){
  maksi<-max(evat,na.rm=TRUE)
  lambda<-p*maksi
}
grt<-(evat>=lambda)

#dendatsub<-dendat[grt,]

return(grt)
}

nn.radit<-function(dendat,maxk)
{
n<-dim(dendat)[1]
radmat<-matrix(0,n,maxk)

eta<-dist(dendat)
#i<j eta[n*(i-1) - i*(i-1)/2 + j-i]

for (i in 2:(n-1)){
   i1<-seq(1,i-1)
   j1<-i
   irow1<-eta[n*(i1-1) - i1*(i1-1)/2 + j1-i1]
   j2<-seq(i+1,n)
   irow2<-eta[n*(i-1) - i*(i-1)/2 + j2-i]
   irow<-c(irow1,irow2)
   or<-order(irow)
   radmat[i,]<-irow[or[1:maxk]]
}

i<-1
j<-seq(i+1,n)
irow<-eta[n*(i-1) - i*(i-1)/2 + j-i]
or<-order(irow)
radmat[i,]<-irow[or[1:maxk]]

i<-n
i1<-seq(1,n-1)
j<-i
irow<-eta[n*(i1-1) - i1*(i1-1)/2 + j-i1]
or<-order(irow)
radmat[i,]<-irow[or[1:maxk]]

return(radmat)
}


nnt<-function(dendat,f)
{
# dendat is n*d
# f is n-vector of evaluations of the function at dendat

n<-dim(dendat)[1]
d<-dim(dendat)[2]
vol<-pi^(d/2)/gamma(d/2+1)
lkm<-n
parent<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
#center<-matrix(0,d,lkm)
#level<-matrix(0,lkm,1)
levset.radius<-matrix(0,lkm,1)

ford<-order(f) # indeksit pienimmasta suurimpaan
root<-ford[n]
leaf<-root
nindit<-nn.indit(dendat)
neig<-nindit[root,]
notvisited<-setdiff(seq(1,lkm),root)  # a[!a %in% b]  setdiff(a, b) 
radi<-sqrt(sum((dendat[neig[1],]-dendat[root,])^2))/2
volume[root]<-vol*radi^d
levset.radius[root]<-radi

cur<-1
for (i in 1:(n-1)){
    smaller<-ford[n-i]
    nearest<-neig[cur]
    if (smaller==nearest){
        parent[root]<-smaller
        dist.to.parent<-sqrt(sum((dendat[smaller,]-dendat[root,])^2))
        radi<-dist.to.parent+levset.radius[root]
        volume[smaller]<-max(vol*radi^d,volume[smaller])
        levset.radius[smaller]<-max(radi,levset.radius[smaller])
        #  pi*dist(rbind(dendat[root,],dendat[smaller,]))^2
        notvisited<-setdiff(notvisited,root) 
        root<-smaller
        cur<-cur+1
    }
    else{
        parent[root]<-nearest
        #volume[root]<-pi*sum((dendat[root,]-dendat[nearest,])^2)
        notvisited<-setdiff(notvisited,root)
        dist.to.parent<-sqrt(sum((dendat[nearest,]-dendat[root,])^2))
        radi<-dist.to.parent+levset.radius[root]
        volume[nearest]<-max(vol*radi^d,volume[nearest])
        levset.radius[nearest]<-max(radi,levset.radius[nearest])
# dist.nearest.to.root.bound<-2*sqrt(sum((dendat[root,]-dendat[nearest,])^2))
# newvolume<-pi*dist.nearest.to.root.bound^2
# volume[nearest]<-max(volume[nearest],newvolume)   
        root<-smaller
        leaf<-root
        visited<-setdiff(seq(1,lkm),notvisited)
        neig<-setdiff(nindit[root,],visited) #nindit[root,]
        if (volume[root]==0){
           radi<-sqrt(sum((dendat[neig[1],]-dendat[root,])^2))/2
           volume[root]<-vol*radi^d
           levset.radius[root]<-radi
        }
        cur<-1 
    }
}

#mt<-multitree(parent)
#for (i in 1:lkm){
#    if (mt$sibling[i]=0){
#       node<-parent[i]
#       volume[node]

lf<-list(
parent=parent,volume=volume,center=t(dendat),level=f,
#root=root,
#infopointer=infopointer,
#refe=refe,maxdis=maxdis,
dendat=dendat)

return(lf)
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
omamax<-function(v){
#v on vektori, palautetaan pienin arvo 
#
lkm<-length(v)
i<-1
while ((i<lkm) && (is.na(v[i]))) i<-i+1
if ((i==lkm) && (is.na(v[lkm]))) y<-NA
 else
 if ((i==lkm) && (!is.na(v[lkm]))) y<-v[lkm]
  else{
  apuu<-i
  valapu<-v[apuu]
  while (i<lkm){
    i<-i+1
    if ((!is.na(v[i])) && (v[i] > valapu)){
      apuu<-i
      valapu<-v[i]
    }
  }
y<-v[apuu]
  }
return(y)
}
omamin<-function(v){
#v on vektori, palautetaan pienin arvo 
#
lkm<-length(v)
i<-1
while ((i<lkm) && (is.na(v[i]))) i<-i+1
if ((i==lkm) && (is.na(v[lkm]))) y<-NA
 else
 if ((i==lkm) && (!is.na(v[lkm]))) y<-v[lkm]
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
y<-v[apuu]
  }
return(y)
}
omaord2<-function(a,b){
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





omaord<-function(values,recs,frekv=NULL){
#Jarjestaa paloittain vakion funktion palat funktion arvojen
#mukaan suuruusjarjestykseen
#
#palvak on lnum*(1+2*xlkm)-matriisi, missa lnum on laatikkojen lkm,
#matriisin ensimmainen sarake sisaltaa estimaatin values laatikoittain,
#naitten mukaan matriisin rivit jarjestetaan.
#Muut sarakkeet sis laatikoitten maaritykset, ts jokaista
#muuttujaa kohden vaihteluvali 
#ep on toleranssiparametri yhtasuuruuden testauksessa
#
#kutsuu: omaind
#
lnum<-length(values)        #length(recs[,1])     #laatikoitten lkm
ordrecs<-recs               #tahan oikea jarjestys
ordvalues<-values
if (is.null(frekv)){
 ordfrekv<-NULL
 i<-1
 while (i<=lnum){
   pienin<-omaind(values)
   ordrecs[i,]<-recs[pienin,]
   ordvalues[i]<-values[pienin]
   values[pienin]<-NA       #NA on plus aareton
   i<-i+1
 }
}
else{
ordfrekv<-frekv
i<-1
 while (i<=lnum){
   pienin<-omaind(values)
   ordrecs[i,]<-recs[pienin,]
   ordvalues[i]<-values[pienin]
   ordfrekv[i]<-frekv[pienin]
   values[pienin]<-NA       #NA on plus aareton
   i<-i+1
 }
}
return(list(values=ordvalues,recs=ordrecs,frekv=ordfrekv))
}





onko<-function(rivi,j){
#Checks whether j is in rivi
#
#rivi is vector where beginning is positive integers, rest NA
#j is positive inetger
#
#Returns TRUE is j is in rivi, FALSE otherwise
#
len<-length(rivi)
res<-FALSE
i<-1
while ((!is.na(rivi[i])) && (i<=len) && (rivi[i]<=j)){
  if (rivi[i]==j) res<-TRUE
  i<-i+1
}
return(res)
}

onsetissa<-function(kandi,h,delta,minim,
brnode,component,
index,
AtomlistAtom,AtomlistNext){
#
itis<-F
d<-length(minim)
# 
node<-brnode
compo<-component[node]
ato<-compo                          #ato is pointer to "value"
while ((ato>0) && !(itis)){
    inde<-index[AtomlistAtom[ato]]
    keski<-minim-h+delta*inde
    for (din in 1:d){
      if ((kandi[din]>=(keski[din]-delta[din]/2)) &&   
          (kandi[din]<=(keski[din]+delta[din]/2))){
               itis<-T
      }
    }
    ato<-AtomlistNext[ato]
}
#
return(itis)
}
paraclus<-function(dendat,algo="kmeans",k=2,method="complete",
scatter=FALSE,coordi1=1,coordi2=2,levelmethod="center",
startind=c(1:k),range="global",terminal=TRUE,coordi=1,
paletti=NULL,xaxt="s",yaxt="s",cex.axis=1,pch.paletti=NULL)
{
if (is.null(paletti)) paletti<-seq(1,2000)
if (is.null(pch.paletti)) pch.paletti<-rep(21,2000)
if (algo!="kmeans"){ 
      method<-algo
      algo<-"hclust"
}

n<-dim(dendat)[1]
d<-dim(dendat)[2]
colot<-c(colors()[2],colors()[3])

if (algo=="kmeans"){
    starters<-dendat[startind,]
    cl<-kmeans(dendat,k,centers=starters)
    ct<-cl$cluster
    centers<-cl$centers
}
else if (algo=="hclust"){
       dis<-dist(dendat)
       hc <- hclust(dis, method=method)
       ct<-cutree(hc,k=k)
       centers<-matrix(0,k,d)
       for (ij in 1:k) centers[ij,]<-mean(data.frame(dendat[(ct==ij),]))
}

# calculate innerlevel
innerlevel<-matrix(0,n,1)
maxlevel<-matrix(0,k,1)
for (i in 1:n){
  classlabel<-ct[i]
  cente<-centers[classlabel,]
  if (levelmethod=="random"){ 
        set.seed(i)
        luku<-runif(1)
  }
  else{ 
        luku<-sqrt(sum((dendat[i,]-cente)^2))
  }
  innerlevel[i]<-luku  
  maxlevel[classlabel]<-max(maxlevel[classlabel],luku)
}
# calculate classlevel
classlevel<-matrix(0,k,1)
i<-2
while (i<=k){
   if (levelmethod=="random"){ 
        classlevel[i]<-classlevel[i-1]+1
   }
   else{
        classlevel[i]<-classlevel[i-1]+maxlevel[i-1]
   }
   i<-i+1
}
# calculate level
level<-matrix(0,n,1)
for (i in 1:n){
    classlabel<-ct[i]
    level[i]<-innerlevel[i]+classlevel[classlabel]
}

if (d<=5){ 
   lkm<-d  
   times<-0
   reminder<-d
} 
else{
   lkm<-5
   times<-floor(d/lkm)
   reminder<-d-lkm*times
}
curcolo<-1
ymin<-0  #min(level)
ymax<-max(level)

if (!terminal){
       coordinate<-coordi
       x<-dendat[,coordinate]
       if (range=="global"){
          xmin<-min(dendat) 
          xmax<-max(dendat)
       }
       else{
          xmin<-min(x)
          xmax<-max(x)
       }
       plot(x="",y="",xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
            xaxt=xaxt,yaxt=yaxt,cex.axis=cex.axis)
       for (j in 1:k){
         if (curcolo==1) curcolo<-2 else curcolo<-1
         polygon(c(xmin,xmax,xmax,xmin),
                 c(classlevel[j],classlevel[j],
                   classlevel[j]+maxlevel[j],classlevel[j]+maxlevel[j]),
                 col=colot[curcolo]) 
       }
       points(x,level,col=paletti[ct],pch=pch.paletti[ct])
       if (scatter) plot(dendat[,coordi1],dendat[,coordi2], col = paletti[ct],
                         xaxt=xaxt,yaxt=yaxt,pch=pch.paletti[ct])

}
########################################################
else{

t<-1
while (t<=times){
   mat<-matrix(c(1:lkm),1,lkm)
   dev.new()
   layout(mat)
   for (i in 1:lkm){
       coordinate<-(times-1)*lkm+i
       x<-dendat[,coordinate]
       if (range=="global"){
          xmin<-min(dendat) 
          xmax<-max(dendat)
       }
       else{
          xmin<-min(x)
          xmax<-max(x)
       }
       plot(x="",y="",xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
            xaxt=xaxt,yaxt=yaxt,cex.axis=cex.axis)
       for (j in 1:k){
         if (curcolo==1) curcolo<-2 else curcolo<-1
         polygon(c(xmin,xmax,xmax,xmin),
                 c(classlevel[j],classlevel[j],
                   classlevel[j]+maxlevel[j],classlevel[j]+maxlevel[j]),
                 col=colot[curcolo]) 
       }
       points(x,level,col=ct)
   }
   t<-t+1
}
if (reminder>0){
   lkm<-reminder
   mat<-matrix(c(1:lkm),1,lkm)
   dev.new()
   layout(mat)
   for (i in 1:lkm){
       coordinate<-i
       x<-dendat[,coordinate]
       if (range=="global"){
          xmin<-min(dendat) 
          xmax<-max(dendat)
       }
       else{
          xmin<-min(x)
          xmax<-max(x)
       }
       plot(x="",y="",xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
            xaxt=xaxt,yaxt=yaxt)
       for (j in 1:k){
         if (curcolo==1) curcolo<-2 else curcolo<-1
         polygon(c(xmin,xmax,xmax,xmin),
                 c(classlevel[j],classlevel[j],
                   classlevel[j]+maxlevel[j],classlevel[j]+maxlevel[j]),
                 col=colot[curcolo]) 
       }
       points(x,level,col=ct)
   }
}

# scatter plot
if (scatter){
   dev.new()
   plot(dendat[,coordi1],dendat[,coordi2], col = ct, xaxt=xaxt, yaxt=yaxt)
}

} # if terminal

}


paracoor.dens<-function(dendat,type="classical",h=1,b=0.25,k=100,m=100,alpha=1)
{
# k<-1000  # grid lkm vaakatasossa
# m<-1000  # grid lkm pystytasossa

n<-dim(dendat)[1]

if (type=="new"){

 vals<-matrix(0,n,1)
 for (i in 1:n){
    arg<-dendat[i,]
    vals[i]<-kernesti.dens(arg,dendat,h=h)
 }
 w<-(vals-min(vals))/(max(vals)-min(vals))
 or<-order(w)
 w2<-(1-w)^b
 paletti<-grey(w2)[or]    
 x<-dendat[or,]
 paracoor(x,paletti=paletti) 

}

if (type=="classical"){

 d<-dim(dendat)[2]
 maks<-matrix(0,d,1)
 mini<-matrix(0,d,1)
 for (i in 1:d){
    maks[i]<-max(dendat[,i])
    mini[i]<-min(dendat[,i])
 }
 dendat2<-dendat
 for (i in 1:d) dendat2[,i]<-(dendat[,i]-mini[i])/(maks[i]-mini[i])
 pc<-matrix(0,m,k*(d-1))
 for (dd in 1:(d-1)){
   for (kk in 1:k){
      x1<-dendat2[,dd]
      x2<-dendat2[,dd+1]
      t<-kk/(k+1)
      datai<-(1-t)*x1+t*x2
      ind<-(dd-1)*k+kk
      for (mm in 1:m){
          arg<-mm/m
          pc[mm,ind]<-kernesti.dens(arg,datai,h=h)
      }
   }
 }
 pc2<-t(pc)^b
 colo<-grey(seq(0,1,0.1),alpha=alpha)
 image(pc2,col=colo)  #image(pc2,col=topo.colors(120))
 #image(pc2,col=terrain.colors(50))
 #heatmap(pc2)
 #contour(pc2)

}

}



paracoor<-function(X,Y=NULL,xmargin=0.1,
paletti=matrix("black",dim(X)[1],1),noadd=TRUE,verti=NULL,cex.axis=1,
points=TRUE,col.verti="black",col.verti.y="red",digits=3,
arg=NULL,colarg="red",lwd=1,cex=1,yaxt="s")
{
n<-dim(X)[1]
d<-dim(X)[2]
ylim<-c(min(X),max(X))
if (is.null(Y)) D<-d else D<-d+dim(Y)[2]

if (noadd)
plot(x="",y="",
xlim=c(1-xmargin,D+xmargin),ylim=ylim,
xlab="",ylab="",xaxt="n",cex.axis=cex.axis,yaxt=yaxt)

for (i in 1:n){
    if (points) points(X[i,],col=paletti[i],cex=cex)
    for (j in 1:(d-1)) segments(j,X[i,j],j+1,X[i,j+1],
                                col=paletti[i],lwd=lwd)
}
#if (points) for (i in 1:n) points(X[i,],col=paletti[i])

if (!is.null(Y)){
  miny<-min(Y)
  maxy<-max(Y)
  z<-matrix(0,n,dim(Y)[2])
  for (i in 1:n){
     for (j in 1:dim(Y)[2]){
         coeff<-(Y[i,j]-miny)/(maxy-miny)
         z[i,j]<-ylim[1]+coeff*(ylim[2]-ylim[1])
     }
  }
  for (i in 1:n){
     j<-2
     while (j<=dim(Y)[2]){ 
        if (points) points(d+j,z[i,j],col=paletti[i],cex=cex)  
        segments(d+j-1,z[i,j-1],d+j,z[i,j],col=paletti[i],lwd=lwd)
        j<-j+1
     }
     if (points){
         points(d+1,z[i,1],col=paletti[i],cex=cex)
         points(d,X[i,d],col=paletti[i],cex=cex)
     }
     segments(d,X[i,d],d+1,z[i,1],col=paletti[i],lwd=lwd)
  }
  #if (points) for (i in 1:n) points(d+1,z[i],col=paletti[i])
  segments(d+0.5,ylim[1],d+0.5,ylim[2],col=col.verti.y,lwd=lwd)
  text(d+dim(Y)[2]+xmargin/2,ylim[1],format(miny,digits=digits))
  text(d+dim(Y)[2]+xmargin/2,ylim[2],format(maxy,digits=digits))
  text(d+dim(Y)[2]+xmargin/2,ylim[1]+(ylim[2]-ylim[1])/2,
       format(miny+(maxy-miny)/2,digits=digits))
}

if (!is.null(verti)) segments(verti,ylim[1],verti,ylim[2],col=col.verti,lwd=lwd)

if (!is.null(arg)){
    if (points) points(arg,col=colarg,cex=cex)
    for (j in 1:(d-1)) segments(j,arg[j],j+1,arg[j+1],col=colarg,lwd=lwd)
}

}

pcf.boundary<-function(dendat,N=rep(10,dim(dendat)[2]-1),
rho=0,m=dim(dendat)[1],seed=1)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

set.seed(seed)
mc<-max(1,round(m/n))
M<-n*mc
data<-matrix(0,M,d-1)
distat<-matrix(0,M,1)
dendat.mc<-matrix(0,M,d)
for (i in 1:n){
    obs<-dendat[i,]
    for (j in 1:mc){
        diro<-2*pi*runif(1)
        riro<-rho*runif(1)
        newobs<-obs+riro*sphere.map(diro)
        len<-sqrt(sum(newobs^2))
        ii<-mc*(i-1)+j
        data[ii,]<-sphere.para(newobs/len)
        distat[ii]<-len
        dendat.mc[ii,]<-newobs
    }
}

support<-matrix(0,2*(d-1),1)
for (i in 1:(d-1)){
    support[2*i-1]<-min(data[,i])
    support[2*i]<-max(data[,i])
}

step<-matrix(0,d-1,1)
for (i in 1:(d-1)) step[i]<-(support[2*i]-support[2*i-1])/N[i]
recnum<-prod(N)
rowpointer<-matrix(0,recnum,1)

value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d-1)

inde<-matrix(0,d-1,1)
numpositive<-0
for (i in 1:M){
    # find the right rectangle
    point<-data[i,]
    for (k in 1:(d-1)) inde[k]<-min(floor((point[k]-support[2*k-1])/step[k]),N[k]-1)
    # inde[k] should be between 0 and N[k]-1

    # find the right row (if already there)
    recnum<-0
    for (kk in 1:(d-1)){
        if (kk==1) tulo<-1 else tulo<-prod(N[1:(kk-1)])
        recnum<-recnum+inde[kk]*tulo
    }
    recnum<-recnum+1
    row<-rowpointer[recnum]

    # update the value or create a new row
    if (row>0) value[row]<-max(value[row],distat[i])
    else{
         numpositive<-numpositive+1
         rowpointer[recnum]<-numpositive
         value[numpositive]<-distat[i]
         index[numpositive,]<-inde
    }
}
value<-value[1:numpositive]
index<-index[1:numpositive,]
if (d==2) index<-matrix(index,length(index),1)
down<-index
high<-index+1

pcf<-list(
value=value,index=NULL,
down=down,high=high,  #step=delta,
support=support,N=N,data=data,dendat.mc=dendat.mc)
return(pcf)
}






pcf.func<-function(func, N,
sig=rep(1,length(N)), support=NULL, theta=NULL, 
g=1, M=NULL, p=NULL, mul=3, t=NULL, 
marginal="normal", r=0,
mu=NULL, xi=NULL, Omega=NULL, alpha=NULL, df=NULL, 
a=0.5, b=0.5, distr=FALSE, std=1, lowest=0) # contrast="loglik")   
{
# t<-rep(1,length(N))

d<-length(N)

if (d>1){

  if (marginal=="unif") support<-c(0,sig[1],0,sig[2])

  recnum<-prod(N)
  value<-matrix(0,recnum,1)
  index<-matrix(0,recnum,d)

  # new ############################################

  if (func=="mixt"){ 

     if (is.null(support)){
       support<-matrix(0,2*d,1)
       for (i in 1:d){
           support[2*i-1]<-min(M[,i]-mul*sig[,i])
           support[2*i]<-max(M[,i]+mul*sig[,i])
       }
     }
     lowsuppo<-matrix(0,d,1)
     for (i in 1:d) lowsuppo[i]<-support[2*i-1]
     step<-matrix(0,d,1)
     for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]
     mixnum<-length(p)

     numpositive<-0
     for (i in 1:recnum){
        inde<-digit(i-1,N)+1
        point<-lowsuppo+step*inde-step/2
 
        if (!is.null(theta)){
           rotmat<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
           point<-rotmat%*%point
        }

        valli<-0
        for (mi in 1:mixnum){
            evapoint<-(point-M[mi,])/sig[mi,]
            valli<-valli+p[mi]*evanor(evapoint)/prod(sig[mi,])
        }
        if (valli>lowest){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive,]<-inde
        }
     }
     value<-value[1:numpositive]
     index<-index[1:numpositive,]
     down<-index-1
     high<-index
  }


  else if (func=="student"){ 
     lowsuppo<-matrix(0,d,1)
     for (i in 1:d) lowsuppo[i]<-support[2*i-1]
     step<-matrix(0,d,1)
     for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

     numpositive<-0
     for (i in 1:recnum){
        inde<-digit(i-1,N)+1
        x<-lowsuppo+step*inde-step/2

        #valli<-eva.student(x,t,marginal,sig,r,df)

        margx<-matrix(0,d,1)
        u<-matrix(0,d,1)

        if (marginal=="unif"){
           for (j in 1:d){
             u[j]<-x[j]/sig[j]  #+1/2
             margx[j]<-1/sig[j]
           }
        }
        if ((marginal=="normal")||(marginal=="gauss")){
           for (j in 1:d){
             u[j]<-pnorm(x[j]/sig[j])
             margx[j]<-evanor(x[j]/sig[j])/sig[j]
           }
        }
        if (marginal=="student"){
          for (j in 1:d){
             u[j]<-pt(x[j]/sig[j],df=t[j])
             margx[j]<-dt(x[j]/sig[j],df=t[j])/sig[j]
          }
        }
        
        x1<-qt(u[1],df=df)
        x2<-qt(u[2],df=df)

        d<-2
        vakio<-gamma((df+d)/2)*gamma(df/2)/gamma((df+1)/2)^2
        nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
        prod<-(1+x1^2/df)^((1+df)/2)*(1+x2^2/df)^((1+df)/2)
        copuval<-vakio*(1-r^2)^(-1/2)*prod*(1+nelio/df)^(-(df+d)/2)

        valli<-copuval*margx[1]*margx[2]

        ###############################################

        if (valli>0){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive,]<-inde
        }
     }
     value<-value[1:numpositive]
     index<-index[1:numpositive,]
     down<-index-1
     high<-index
  }

  else if (func=="gauss"){ 
     lowsuppo<-matrix(0,d,1)
     for (i in 1:d) lowsuppo[i]<-support[2*i-1]
     step<-matrix(0,d,1)
     for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

     numpositive<-0
     for (i in 1:recnum){
        inde<-digit(i-1,N)+1
        x<-lowsuppo+step*inde-step/2

        #valli<-eva.copula(x,type="gauss",marginal=marginal,sig=sig,r=r,t=t)

        margx<-matrix(0,d,1)
        u<-matrix(0,d,1)

        if (marginal=="unif"){
           for (j in 1:d){
             u[j]<-x[j]/sig[j]  #+1/2
             margx[j]<-1/sig[j]
           }
        }
        if ((marginal=="normal")||(marginal=="gauss")){
           for (j in 1:d){
             u[j]<-pnorm(x[j]/sig[j])
             margx[j]<-evanor(x[j]/sig[j])/sig[j]
           }
        }
        if (marginal=="student"){
          for (j in 1:d){
             u[j]<-pt(x[j]/sig[j],df=t[j])
             margx[j]<-dt(x[j]/sig[j],df=t[j])/sig[j]
          }
        }
        
        x1<-qnorm(u[1],sd=1)
        x2<-qnorm(u[2],sd=1)

        nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
        copuval<-(1-r^2)^(-1/2)*exp(-nelio/2)/exp(-(x1^2+x2^2)/2)

        valli<-copuval*margx[1]*margx[2]

        ########################################

        if (valli>0){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive,]<-inde
        }
     }
     value<-value[1:numpositive]
     index<-index[1:numpositive,]
     down<-index-1
     high<-index
  }


else{

# old #########################################################

if (is.null(support)){

   if (func=="epan"){
      if (is.null(sig)) sig<-c(1,1)
      support<-matrix(0,2*d,1)
      for (i in 1:d){
          support[2*i-1]<--sig[i]
          support[2*i]<-sig[i]
      }
   }

}

if ((marginal=="unif")) support<-c(0,sig[1],0,sig[2])
# && (is.null(support))) 
#support<-c(-sig[1]/2,sig[1]/2,-sig[2]/2,sig[2]/2)


lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

numpositive<-0
for (i in 1:recnum){
    inde<-digit(i-1,N)+1
    #if ((inde[1]==0) && (inde[2]==N[2])) inde<-c(0,0)
    point<-lowsuppo+step*inde-step/2

    if (!is.null(theta)){
         rotmat<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
         point<-rotmat%*%point
    }

    if (func=="prod") valli<-eva.prod(point,marginal,g)
    if (func=="skewgauss") valli<-eva.skewgauss(point,mu,sig,alpha)
    #if (func=="dmsn") valli<-dmsn(point,xi,Omega,alpha)
    if (func=="gumbel") valli<-eva.copula(point,
        type="gumbel",marginal=marginal,sig=sig,r=r,t=t,g=g)
    if (func=="frank") valli<-eva.copula(point,
        type="frank",marginal=marginal,sig=sig,t=t,g=g)
    if (func=="plackett") valli<-eva.plackett(point,t,marginal,sig)
    if (func=="clayton2") valli<-eva.clayton(point,t,marginal,sig,df)
    if (func=="clayton") valli<-eva.copula(point,
        type="clayton",marginal=marginal,sig=sig,r=r,t=t,g=g)
    if (func=="cop6") valli<-eva.cop6(point,t,marginal,sig)
    if (func=="epan") valli<-epan(point)
    if (func=="normal") 
        valli<-eva.gauss(point,t=t,marginal=marginal,sig=sig,r=r)   
    if (func=="hat") valli<-eva.hat(point,a=a,b=b)

    if (valli>0){
       numpositive<-numpositive+1
       value[numpositive]<-valli
       index[numpositive,]<-inde
    }
}

value<-value[1:numpositive]
index<-index[1:numpositive,]
down<-index-1
high<-index

}


pcf<-list(
value=value,index=index,
down=down,high=high,  #step=delta,
support=support,N=N)

  #pcf<-eval.func.dD(func,N,
  #sig=sig,support=support,theta=theta,g=g,
  #M=M,p=p,mul=mul,
  #t=t,marginal=marginal,r=r, 
  #mu=mu,xi=xi,Omega=Omega,alpha=alpha,df=df,a=a,b=b)

}

else{  # (d==1){ ######################################################

  pcf<-eval.func.1D(func,N,
  support=support,g=g,std=std,distr=distr,
  M=M,sig=sig,p=p,
  a=a,b=b,d=2)

}


return(pcf)

}

pcf.histo<-function(dendat,N,weights=rep(1,dim(dendat)[1]))
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]
support<-matrix(0,2*d,1)
for (i in 1:d){
       support[2*i-1]<-min(dendat[,i])
       support[2*i]<-max(dendat[,i])
}

step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]
recnum<-prod(N)
rowpointer<-matrix(0,recnum,1)

value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

inde<-matrix(0,d,1)
numpositive<-0
for (i in 1:n){
    # find the right rectangle
    point<-dendat[i,]
    weight<-weights[i]
   for (k in 1:d) inde[k]<-min(floor((point[k]-support[2*k-1])/step[k]),N[k]-1)
    # inde[k] should be between 0 and N[k]-1

    # find the right row (if already there)
    recnum<-0
    for (kk in 1:d){
        if (kk==1) tulo<-1 else tulo<-prod(N[1:(kk-1)])
        recnum<-recnum+inde[kk]*tulo
    }
    recnum<-recnum+1
    row<-rowpointer[recnum]

    # update the value or create a new row
    if (row>0) value[row]<-value[row]+weight
    else{
         numpositive<-numpositive+1
         rowpointer[recnum]<-numpositive
         value[numpositive]<-weight
         index[numpositive,]<-inde
    }
}
value<-value[1:numpositive]
index<-index[1:numpositive,]
down<-index
high<-index+1

pcf<-list(
value=value,index=NULL,
down=down,high=high,  #step=delta,
support=support,N=N)
return(pcf)
}




pcf.kernC<-function(dendat,h,N,kernel="epane",hw=NULL)
# creates piecewise constant function object
{
keva<-kereva(dendat,h,N,kernel=kernel,hw=hw)

d<-length(N)
recnum<-dim(keva$index)[1]
down<-matrix(0,recnum,d)
high<-matrix(0,recnum,d)
for (i in 1:recnum){
     down[i,]<-keva$index[i,]-1
     high[i,]<-keva$index[i,]
}

return(list(value=keva$value,down=down,high=high,N=N,support=keva$suppo,
index=keva$index))
}


pcf.kern<-function(dendat,h,N,kernel="gauss",weights=NULL,support=NULL,
lowest=0,radi=0)
{
d<-length(N)

if (d>1){

if (length(h)==1) h<-rep(h,d)

if (kernel=="bart") 
   ker<-function(xx,d){ 
         musd<-2*pi^(d/2)/gamma(d/2)
         c<-d*(d+2)/(2*musd)
         return( c*(1-rowSums(xx^2))*(rowSums(xx^2) <= 1) ) 
   }
if (kernel=="gauss") 
   ker<-function(xx,d){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx,d){ 
         c<-gamma(d/2+1)/pi^(d/2) 
         return( (rowSums(xx^2) <= 1) ) 
   } 
if (kernel=="epane") 
   ker<-function(xx,d){ 
      c<-(3/4)^d 
      xxx<-(1-xx^2)*(1-xx^2>=0)
      return( c*apply(xxx,1,prod) ) 
   } 

if (is.null(radi)) if (kernel=="gauss") radi<-2*h else radi<-h

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (is.null(support)){
  support<-matrix(0,2*d,1)
  for (i in 1:d){
     support[2*i-1]<-min(dendat[,i])-radi
     support[2*i]<-max(dendat[,i])+radi
  }
}
lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

numpositive<-0
for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     argu<-matrix(arg,dim(dendat)[1],d,byrow=TRUE)
#     neigh<-(rowSums((argu-x)^2) <= radi^2)
#     if (sum(neigh)>=2){     # if there are obs in the neigborhood
#
#       xred<-dendat[neigh,]
#       argu<-matrix(arg,dim(xred)[1],d,byrow=TRUE)

       xxx<-sweep(dendat-argu,2,h,"/")
       w<-ker(xxx,d)/prod(h)
       valli<-mean(w)
       if (!is.null(weights)) valli<-t(weights)%*%w
#     }
#     else valli<-mean(y)

      if (valli>lowest){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive,]<-inde
      }
      #value[i]<-valli
      #index[i,]<-inde

}

value<-value[1:numpositive]
index<-index[1:numpositive,]
down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}
else{  # d==1  #########################################

d<-1
x<-matrix(dendat,length(dendat),1)

if (kernel=="gauss") ker<-function(xx,d){ return( (2*pi)^(-1/2)*exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx,d){ return( (abs(xx) <= 1) ) }

index<-seq(1:N)
len<-length(index)

value<-matrix(0,N,1)
if (is.null(support)){
   support<-matrix(0,2,1)
   support[1]<-min(x)
   support[2]<-max(x)
}
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

numpositive<-0
for (i in 1:N){
     inde<-i
     argu<-lowsuppo+step*inde-step/2
     w<-ker((x-argu)/h,1)/h
     if (!is.null(weights)) valli<-t(weights)%*%w else valli<-mean(w)
     if (valli>lowest){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive]<-inde
     }
}

value<-value[1:numpositive]
index<-index[1:numpositive]

down<-matrix(0,numpositive,1)
high<-matrix(0,numpositive,1)
down[,1]<-index-1
high[,1]<-index

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

}

return(pcf)
}
pcf.kern.vech<-function(dendat,h,N,kernel="gauss",weights=NULL,support=NULL,
lowest=0,radi=0)
{
d<-length(N)

if (d>1){

if (kernel=="bart") 
   ker<-function(xx,d){ 
         musd<-2*pi^(d/2)/gamma(d/2)
         c<-d*(d+2)/(2*musd)
         return( c*(1-rowSums(xx^2))*(rowSums(xx^2) <= 1) ) 
   }
if (kernel=="gauss") 
   ker<-function(xx,d){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx,d){ 
         c<-gamma(d/2+1)/pi^(d/2) 
         return( (rowSums(xx^2) <= 1) ) 
   } 
if (kernel=="epane") 
   ker<-function(xx,d){ 
      c<-(3/4)^d 
      xxx<-(1-xx^2)*(1-xx^2>=0)
      return( c*apply(xxx,1,prod) ) 
   } 

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (is.null(support)){
  support<-matrix(0,2*d,1)
  for (i in 1:d){
     support[2*i-1]<-min(dendat[,i])-radi
     support[2*i]<-max(dendat[,i])+radi
  }
}
lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

numpositive<-0
for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     argu<-matrix(arg,dim(dendat)[1],d,byrow=TRUE)

     w<-ker((dendat-argu)/h,d)/prod(h)
     valli<-mean(w)
     if (!is.null(weights)) valli<-t(weights)%*%w

     if (valli>lowest){
         numpositive<-numpositive+1
         value[numpositive]<-valli
         index[numpositive,]<-inde
      }
}

value<-value[1:numpositive]
index<-index[1:numpositive,]
down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}
else{  # d==1  #########################################

d<-1
x<-matrix(dendat,length(dendat),1)

if (kernel=="gauss") ker<-function(xx,d){ return( (2*pi)^(-1/2)*exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx,d){ return( (abs(xx) <= 1) ) }

index<-seq(1:N)
len<-length(index)

value<-matrix(0,N,1)
if (is.null(support)){
   support<-matrix(0,2,1)
   support[1]<-min(x)
   support[2]<-max(x)
}
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

numpositive<-0
for (i in 1:N){
     inde<-i
     argu<-lowsuppo+step*inde-step/2
     w<-ker((x-argu)/h,1)/h
     if (!is.null(weights)) valli<-t(weights)%*%w else valli<-mean(w)
     if (valli>lowest){
           numpositive<-numpositive+1
           value[numpositive]<-valli
           index[numpositive]<-inde
     }
}

value<-value[1:numpositive]
index<-index[1:numpositive]

down<-matrix(0,numpositive,1)
high<-matrix(0,numpositive,1)
down[,1]<-index-1
high[,1]<-index

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

}

return(pcf)
}
pcf.matrix<-function(A)
{
d<-2
num<-dim(A)[1]
N<-c(num,num)
recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

for (i in 1:recnum){
    inde<-digit(i-1,N)+1
    value[i]<-A[inde[1],inde[2]]
    index[i,]<-inde
}
down<-index-1
high<-index
#support<-c(0,num+1,0,num+1)
support<-c(1,num,1,num)

pcf<-list(
value=value,index=index,
down=down,high=high,  #step=delta,
support=support,N=N)
return(pcf)
}

perspec.dyna<-function(x,y,z,col="black",phi=10,theta=0)
{
persp(x=x,y=y,z=z,col=col,
xlab="level",ylab="h",zlab="",ticktype="detailed",
phi=phi,theta=theta)

loc<-locator(1)
ycor<-loc$y 

alaraja<--0.4
while (loc$y>=alaraja){

     if (loc$x>=0) theta<-theta+10 else theta<-theta-10
     if (loc$y>=0) phi<-phi+10 else phi<-phi-10

     persp(x=x,y=y,z=z,col=col,
     xlab="level",ylab="h",zlab="",ticktype="detailed",
     phi=phi,theta=theta)

     loc<-locator(1)
}
dev.off()
}

pituus<-function(x){
#laskee euklid pituuden nelion matriisien x riveille
#
d<-length(x[1,])
lkm<-length(x[,1])
vast<-matrix(0,lkm,1)
i<-1
while (i<=lkm){
  j<-1
  while (j<=d){
    vast[i]<-vast[i]+(x[i,j])^2
    j<-j+1
  }
  i<-i+1
}
return(t(vast))
}
plotbary<-function(lst,coordi=1,
plot=TRUE,data=FALSE,crit=NULL,orderrule="distcenter",
modelabel=FALSE,ptext=0,leimat=NULL,symbo=NULL,
info=NULL,infolift=0,infopos=0,
xmarginleft=0,xmarginright=0,ymargin=0,
xlim=NULL,ylim=NULL,xaxt="s",yaxt="s",
nodesymbo=20,col=NULL,col.axis="black",collines=NULL,paletti=NULL,
shift=0,shiftindex=NULL,
modlabret=FALSE,modecolo=NULL,modepointer=NULL,colometh="lst",
colothre=min(lst$level),lines=TRUE,wedge=FALSE,lty.wedge=2,
title=TRUE,titletext="coordinate",
cex=NULL,nodemag=NULL,cex.sub=1,cex.axis=1,newtitle=FALSE,cex.lab=1,
lowest="dens",subtitle=NULL)
{

parent<-lst$parent
center<-lst$center
level<-lst$level

if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])
if (is.null(col)) 
   if (colometh=="lst")
            col<-colobary(parent,paletti,
                 modecolo=modecolo,modepointer=modepointer)
   else col<-colobary.roots(lst$parent,lst$level,paletti=paletti,
                            colothre=colothre)

if (is.null(collines)) collines<-col

nodenum<-length(parent)
xcoordinate<-center[coordi,]

if (is.null(xlim))
xlim<-c(min(xcoordinate)-xmarginleft,max(xcoordinate)+xmarginright)
if (lowest=="dens") lowesti<-0 else lowesti<-min(lst$level)
ylim<-c(lowesti,max(level)+ptext+ymargin)

if (newtitle) xlab<-paste(titletext,as.character(coordi))
else xlab<-""
plot(xcoordinate,level,xlab=xlab,ylab="",
xlim=xlim,ylim=ylim,xaxt=xaxt,yaxt=yaxt,
pch=nodesymbo,col=col,col.axis=col.axis,cex=nodemag,
cex.axis=cex.axis,cex.lab=cex.lab) 
if (!is.null(subtitle)){ 
   title<-FALSE
   title(sub=subtitle,cex.sub=cex.sub)
}
if (title) title(sub=paste(titletext,as.character(coordi)),cex.sub=cex.sub)


if (lines){
   for (i in 1:nodenum){
       if (parent[i]>0){
           xchild<-xcoordinate[i]
           ychild<-level[i]
           xparent<-xcoordinate[parent[i]]
           yparent<-level[parent[i]]
           if (length(collines)>1) colli<-collines[i] else colli<-collines
           segments(xparent,yparent,xchild,ychild,col=colli) 
        }
   }
}

if (wedge){
  maxx<-max(xcoordinate)
  minx<-min(xcoordinate)
  righthigh<-maxx-lst$refe[coordi]
  lefthigh<-lst$refe[coordi]-minx
  segments(lst$refe[coordi],0,maxx,righthigh,lty=lty.wedge)
  segments(lst$refe[coordi],0,minx,lefthigh,lty=lty.wedge)
}

#########
#########
if (modlabret) modelabel<-TRUE
if (modelabel){

data<-plotprof(lst,plot=F,data=T,cutlev=NULL,ptext=NULL,info=NULL,
infolift=0,infopos=0,crit=crit,orderrule=orderrule)
vecs<-data$vecs
mlkm<-moodilkm(parent)
modloc<-mlkm$modloc 

nodenum<-length(vecs[,1])
xcoor<-matrix(0,2*nodenum,1)
ycoor<-matrix(0,2*nodenum,1)

for (i in 1:nodenum){
 xcoor[2*i-1]<-vecs[i,1]
 xcoor[2*i]<-vecs[i,3]
 ycoor[2*i-1]<-vecs[i,2]
 ycoor[2*i]<-vecs[i,4]
}                     
moodinum<-length(modloc)
modelocx<-matrix(0,moodinum,1)
modelocy<-matrix(0,moodinum,1)
if (is.null(leimat)){
   if (is.null(symbo)){
       labels<-paste("M",1:moodinum,sep="")
   }
   else{
         if (symbo=="empty") labels<-paste("",1:moodinum,sep="")
         else  labels<-paste(symbo,1:moodinum,sep="")
   }
}
else{
   labels<-leimat
}
xcor<-matrix(0,moodinum,1)
for (i in 1:moodinum){
    loc<-modloc[i]
    xcor[i]<-xcoor[2*loc-1]
}
modloc<-omaord2(modloc,xcor)
for (i in 1:moodinum){
    loc<-modloc[i]
    modelocx[i]<-xcoordinate[loc]
    modelocy[i]<-level[loc]+ptext
}
if (!is.null(shiftindex)) modelocx[shiftindex]<-modelocx[shiftindex]+shift
text(modelocx,modelocy,labels,cex=cex)         

if (modlabret){ 
   d<-dim(lst$center)[1]
   modelocat<-matrix(0,moodinum,d)
   for (i in 1:moodinum){
       loc<-modloc[i]
       modelocat[i,]<-lst$center[,loc]
   }   
   return(list(modelocat=modelocat,labels=labels))
}

}
############################################
}











plotbary.slide<-function(tt)
{
d<-dim(tt$center)[1]
coordi<-1
plotbary(tt,paletti=seq(1:1000),coordi=coordi)
loc<-locator(1)
while (loc$y>=0){
    if (coordi==d) coordi<-1 else coordi<-coordi+1         
    plotbary(tt,paletti=seq(1:1000),coordi=coordi)          
    loc<-locator(1)
}

}

plotbranchmap<-function(bm,phi=55,theta=30)
{

persp(x=bm$level,y=bm$h,z=bm$z, 
xlab="level",ylab="h",zlab="excess mass",
ticktype="detailed", border=NA, shade=0.75,
col=bm$col,phi=phi,theta=theta) 

}


plot.complex<-function(complex,dendat,xlab="",ylab="",cex.lab=1,cex.axis=1,pch=19,
col=NULL,border="black")
{
plot(dendat,xlab=xlab,ylab=ylab,cex.lab=cex.lab,cex.axis=cex.axis,pch=pch)
lkm<-dim(complex)[1]
for (i in 1:lkm){
    cur<-complex[i,]
    x<-dendat[cur,1]
    y<-dendat[cur,2]
    polygon(x,y,col=col,border=border)
}

}

plotdata<-function(roots,child,sibling,sibord,levels,volumes,vecs)
{
#plots level-set profile

#parents<-c(0,1,1,0,4,2)
#levels<-c(1,2,2,1,2,3)
#volumes<-c(4,2,1,2,1,1)

itemnum<-length(volumes)

#vecs<-matrix(NA,itemnum,4)
#vecs<-alloroot(vecs,roots,sibord,levels,volumes)

rootnum<-length(roots)
left<-child
right<-sibling

for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    pino[1]<-roots[i]  
    pinin<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        if (left[cur]>0){     #if not leaf (root may be leaf)
           vecs<-allokoi(vecs,cur,child,sibling,sibord,levels,volumes)   
        }
        if (right[cur]>0){    #if right exists, put to stack
            pinin<-pinin+1
            pino[pinin]<-right[cur]
        }
        while (left[cur]>0){    #go to leaf and put right nodes to stack
             cur<-left[cur]
             if (left[cur]>0){  #if not leaf
                vecs<-allokoi(vecs,cur,child,sibling,sibord,levels,volumes)
             }
             if (right[cur]>0){ #if right exists, put to stack
                pinin<-pinin+1
                pino[pinin]<-right[cur]
             }
        }
    }
}       
#
return(vecs)
}







plotdelineator<-function(shtseq,coordi=1,ngrid=40,shift=0.05,
volumefunction=NULL,redu=TRUE,type="l")
{
if (is.null(volumefunction)){
   lnum<-length(shtseq$level)
   st<-shtseq$shtseq[[1]]
   td<-treedisc(st,shtseq$pcf,ngrid=ngrid)
   #td<-prunemodes(td,exmalim=0.5)$lst
   reduseq<-list(td)
   for (i in 2:lnum){
       st<-shtseq$shtseq[[i]]
       td<-treedisc(st,shtseq$pcf,ngrid=ngrid)
       #td<-prunemodes(td,exmalim=0.00001)$lst
       reduseq<-c(reduseq,list(td))
   }
   estiseq<-list(lstseq=reduseq,hseq=shtseq$level)
   mg<-modegraph(estiseq)
   plotmodet(mg,coordi=coordi,shift=shift)
}
else{
    vd<-volumefunction
    if (redu){
       x<-vd$delineator.redu[,coordi]
       y<-vd$delineatorlevel.redu
       or<-order(x)
       x1<-x[or]
       y1<-y[or]
       plot(x1,y1,type=type,
            ylab="level",xlab=paste("coordinate",as.character(coordi)))
    }
    else
       plot(vd$delineator[,coordi],vd$delineatorlevel,ylab="level")
}    

}


plotexmap<-function(sp,mt,
xaxt="n",
lift=0.1,leaflift=0.1,ylim=NULL,
leafcolors=NULL
)
{
if (is.null(leafcolors)) lc<-mt$colot
c2s<-colo2scem(sp,mt,lc)

plotvecs(sp$bigvecs,sp$bigdepths,
lift=lift,xaxt="n",
ylim=ylim,
#ylim=c(horilines[length(horilines)],horilines[1]),   #hseq[1]),
leafcolors=c2s,leaflift=leaflift)                        #log="y")

}

plot.histdata<-function(dendat,col,pcf,i1=1,i2=2,i3=3,
simple=FALSE,cut=dim(dendat)[1],
xlab="",xlim=NULL,ylim=NULL,cex.axis=1)
{
d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]=(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i]

ord<-order(dendat[,i3])#,decreasing=TRUE)
ord<-ord[1:cut]
deudat<-dendat[ord,c(i1,i2)]

if (is.null(xlim)){
   xmin<-min(deudat[,1])
   xmax<-max(deudat[,1])
   xlim<-c(xmin,xmax)
}
if (is.null(ylim)){
   ymin<-min(deudat[,2])
   ymax<-max(deudat[,2])
   ylim=c(ymin,ymax)
}

if (simple) 
plot(deudat,pch=19,col=col[ord],xlab=xlab,ylab="",cex.axis=cex.axis,
xlim=xlim,ylim=ylim)

if (!simple){
pointx<-(xlim[1]+xlim[2])/2
pointy<-(ylim[1]+ylim[2])/2
plot(pointx,pointy,type="n",ylab="",xlim=xlim,ylim=ylim,,cex.axis=cex.axis,
pch=20,xlab=xlab)
for (i in 1:cut){
     mu1<-deudat[i,1]
     mu2<-deudat[i,2]
     x1<-mu1-step[i1]/2
     x2<-mu1+step[i1]/2
     y1<-mu2-step[i2]/2
     y2<-mu2+step[i2]/2
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=col[i])
}
}

}


plot.histo<-function(pcf,col=NULL,cex.axis=1,cex.lab=1,ylab="",xlab="",
xaxt="s",yaxt="s")
{
if (is.null(col)){
   f0<-sqrt(pcf$value)  #f0<-pcf$value
   colo<-1-(f0-min(f0)+0.5)/(max(f0)-min(f0)+0.5)
   #colo<-1-(f0-min(f0)+0.02)/(max(f0)+0.05-min(f0)+0.02)
   #colo<-1-(f0-min(f0))/(max(f0)-min(f0))
   col<-gray(colo)
}

d<-length(pcf$N)
step<-matrix(0,d,1)
for (i in 1:d) step[i]=(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i];

xmin<-pcf$support[1]
xmax<-pcf$support[2]
ymin<-pcf$support[3]
ymax<-pcf$support[4]

plot(xmin,ymin,type="n",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
pch=20,cex.axis=cex.axis,cex.lab=cex.lab,ylab=ylab,xlab=xlab,xaxt=xaxt,yaxt=yaxt)

lenni<-length(pcf$value)
for (i in 1:lenni){
     x1<-pcf$support[1]+step[1]*pcf$down[i,1]
     x2<-pcf$support[1]+step[1]*pcf$high[i,1] 
     y1<-pcf$support[3]+step[2]*pcf$down[i,2]
     y2<-pcf$support[3]+step[2]*pcf$high[i,2] 

     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=col[i],lty="blank")
}

}



plotinfo<-function(vecs,info,pos=0,adj=NULL,lift=0,digits=3){
#
nodenum<-length(vecs[,1])
#
#remain<-data$remain
#if (!is.null(remain)){  #if we have cutted, cut also info
#   lenrem<-length(remain)
#   newinfo<-matrix(0,lenrem,1) 
#   for (i in 1:lenrem){
#      point<-remain[i]
#      newinfo[i]<-info[point]
#   }
#   info<-newinfo
##  orinodenum<-length(info)
##  newinfo<-matrix(0,orinodenum,1)
##  ind<-1
##  for (i in 1:orinodenum){
##     if (remain[i]==1){  
##        newinfo[ind]<-info[i]
##        ind<-ind+1
##     }
##  }
##  info<-newinfo[1:nodenum]
#}
##
infolocx<-matrix(nodenum,1)
infolocy<-matrix(nodenum,1)
#
for (i in 1:nodenum){
  infolocx[i]<-vecs[i,3]   #+(vecs[i,3]-vecs[i,1])/2  
  infolocy[i]<-vecs[i,2]+lift
}
info<-format(info,digits=digits)
text(infolocx,infolocy,info,pos,adj)
}


plot.kernscale<-function(scale,pnum=60,maxy0=0,dens=FALSE,cex.axis=1)
{
   hnum<-length(scale$hseq)
   for (i in 1:hnum){
     pk<-scale$pcfseq[[i]]
     dp<-draw.pcf(pk,dens=dens,pnum=pnum)
     if (i==1){ 
           minx<-min(dp$x) 
           miny<-min(dp$y)
           maxx<-max(dp$x) 
           maxy<-max(dp$y)  
     }
     else{ 
          minx<-min(minx,min(dp$x))
          miny<-min(miny,min(dp$y))
          maxx<-max(maxx,max(dp$x))
          maxy<-max(maxy,max(dp$y))
     }
   }
   maxy<-max(maxy,maxy0)
   plot(x="",y="",xlim=c(minx,maxx),ylim=c(miny,maxy),xlab="",ylab="",
   cex.axis=cex.axis)
   for (i in 1:hnum){
     pk<-scale$pcfseq[[i]]
     dp<-draw.pcf(pk,dens=dens,pnum=pnum)
     matpoints(dp$x,dp$y,type="l")
   }
}


plotmodet<-function(mt,coordi=1,colot=NULL,
shift=0,xlim=NULL,xlab="",ylab="",
horilines=NULL,
symbo=20,loga=NULL,lty="dashed",
cex.axis=1,title=TRUE,cex.sub=1,cex.lab=1,
xaxt="s",yaxt="s")
{
epsi<-0.0000001
if (!is.null(horilines)) horilines<-mt$hseq[horilines]

if (is.null(loga)) 
   if (!is.null(mt$type)){
       if (mt$type=="greedy") loga<-"not"
       if (mt$type=="bagghisto") loga<-"not"
       if (mt$type=="carthisto") loga<-"not"
       if (mt$type=="kernel") loga<-"y"
   }
   else loga<-"y"

d<-dim(mt$xcoor)[2]

if (is.null(colot)){
    as<-mt$colot
}
else if (colot=="black"){
   lenni<-length(mt$ycoor)
   as<-matrix("black",lenni,1)
}
else as<-colot

if (d==1) xvec<-mt$xcoor
else xvec<-mt$xcoor[,coordi]
yvec<-mt$ycoor
len<-length(xvec)
for (i in 1:len){
  j<-i+1
  while (j<=len){
    if ((xvec[i]<=xvec[j]+epsi)&&(xvec[i]>=xvec[j]-epsi)&& 
        (yvec[i]<=yvec[j]+epsi)&&(yvec[i]>=yvec[j]-epsi)){
        #&&(as[i]!=as[j])){
             #xvec[j]<-xvec[j]+shift
             xvec[i]<-xvec[i]+shift
    }    
    j<-j+1
  }
}
if (loga=="y")
plot(xvec,yvec,col=as,xlim=xlim,xlab=xlab,ylab=ylab,pch=symbo,log=loga,
     cex.axis=cex.axis,cex.lab=cex.lab,xaxt=xaxt,yaxt=yaxt)   
else
plot(xvec,yvec,col=as,xlim=xlim,xlab=xlab,ylab=ylab,pch=symbo,
     cex.axis=cex.axis,cex.lab=cex.lab,xaxt=xaxt,yaxt=yaxt) 

if (title) title(sub=paste("coordinate",as.character(coordi)),cex.sub=cex.sub)


if (!is.null(horilines)){
  xmin<-min(xvec)
  xmax<-max(xvec)
  horilen<-length(horilines)
  for (i in 1:horilen){
    segments(xmin,horilines[i],xmax,horilines[i],lty=lty)
  }
}

itemnum<-length(mt$parent)
for (i in 1:itemnum){
    if (mt$parent[i]>0){
        xchild<-mt$xcoor[i,coordi]
        #if (loga=="y") ychild<-mt$ycoor[i] else 
        ychild<-mt$ycoor[i]
        xparent<-mt$xcoor[mt$parent[i],coordi]
        #if (loga=="y") yparent<-mt$ycoor[mt$parent[i]] else 
        yparent<-mt$ycoor[mt$parent[i]]
        collo<-mt$colot[i]  #mt$parent[i]]
        segments(xparent,yparent,xchild,ychild,col=collo)
     }
}

}

plotprof<-function(profile,length=NULL,
plot=TRUE,data=FALSE,crit=NULL,orderrule="distcenter",
modelabel=TRUE,ptext=0,leimat=NULL,symbo=NULL,
info=NULL,infolift=0,infopos=0,
xmarginleft=0,xmarginright=0,ymargin=0,
xlim=NULL,ylim=NULL,axes=TRUE,
col="black",col.axis="black",
cutlev=NULL,xaxt="n",exmavisu=NULL,cex.axis=1,cex=1)
{

#xaxs="e"    (extended)  not implemented?  xaxt="n"

parents<-profile$parent
levels<-profile$level
if (is.null(length)) length<-profile$volume
center<-profile$center

mut<-multitree(parents)
if (is.null(profile$roots)) roots<-mut$roots else roots<-profile$roots
child<-mut$child
sibling<-mut$sibling

d<-dim(center)[1]
if (is.null(crit)){
   crit<-rep(0,d)          #order so that 1st closest to origo
   if (d==1) crit<-max(center)
   if (!is.null(profile$refe)) crit<-profile$refe
}

if (orderrule=="distcenter") sibord<-siborder(mut,crit,profile$distcenter)
else sibord<-siborder(mut,crit,center)

itemnum<-length(parents)
vecs<-matrix(NA,itemnum,4)
vecs<-alloroot(vecs,roots,sibord,levels,length)
vecs<-plotdata(roots,child,sibling,sibord,levels,length,vecs)
orivecs<-vecs

if (!(is.null(cutlev))){
  cm<-cutmut(mut,cutlev,levels)              # cut the tree
  roots<-cm$roots
  sibling<-cm$sibling
  mut$roots<-roots
  if (orderrule=="distcenter") sibord<-siborder(mut,crit,profile$distcenter)
  else sibord<-siborder(mut,crit,center)
  rootnum<-length(roots) 
  apuvecs<-matrix(NA,itemnum,4)
  for (i in 1:rootnum){
     inde<-roots[i]
     apuvecs[inde,]<-vecs[inde,]
  }
  vecs<-apuvecs          #we give for the roots the previous positions
  vecs<-plotdata(roots,child,sibling,sibord,levels,length,vecs)
}

if (plot==TRUE){
   if (!(is.null(cutlev))){
     xlim<-c(omamin(vecs[,1])-xmarginleft,omamax(vecs[,3])+xmarginright)
     ylim<-c(omamin(vecs[,2]),omamax(vecs[,2])+ptext+ymargin)
   }
   else{
     xlim<-c(omamin(vecs[,1])-xmarginleft,omamax(vecs[,3])+xmarginright)
     if (is.null(ylim)) ylim<-c(0,omamax(vecs[,2])+ptext+ymargin)
   }
   plotvecs(vecs,segme=T,xlim=xlim,ylim=ylim,xaxt=xaxt,
   col=col,col.axis=col.axis,cex.axis=cex.axis,axes=axes)
   # use original vectors (numbering will be correct)
   if (modelabel){
      plottext(parents,orivecs,ptext,leimat,symbo,cex=cex)  
   }
   if (!is.null(info)){
      plotinfo(vecs,info,pos=infopos,adj=NULL,lift=infolift,digits=3)
   }
}
#
#
if (data==TRUE){
 return(list(sibord=t(sibord),vecs=vecs,parents=parents,levels=levels,
 length=length,center=center,remain=NULL))
}

}












plottext<-function(parents,vecs,lift=0,leimat=NULL,symbo=NULL,cex=NULL){
#
mlkm<-moodilkm(parents)
modloc<-mlkm$modloc
#
nodenum<-length(vecs[,1])
xcoor<-matrix(0,2*nodenum,1)
ycoor<-matrix(0,2*nodenum,1)
#
for (i in 1:nodenum){
 xcoor[2*i-1]<-vecs[i,1]
 xcoor[2*i]<-vecs[i,3]
 ycoor[2*i-1]<-vecs[i,2]
 ycoor[2*i]<-vecs[i,4]
}                          
#
#mindiff<-vecs[nodenum,2]-vecs[1,2]
#for (i in 1:(nodenum-1)){
#  diff<-vecs[(i+1),2]-vecs[i,2]
#  if (diff>0) mindiff<-min(diff,mindiff)  
#}
#lift<-mindiff/5
#
moodinum<-length(modloc)
modelocx<-matrix(0,moodinum,1)
modelocy<-matrix(0,moodinum,1)
if (is.null(leimat)){
   if (is.null(symbo)){
       labels<-paste("M",1:moodinum,sep="")
   }
   else{
         if (symbo=="empty") labels<-paste("",1:moodinum,sep="")
         else  labels<-paste(symbo,1:moodinum,sep="")
   }
} 
else{
   labels<-leimat
}
xcor<-matrix(0,moodinum,1)
for (i in 1:moodinum){
    loc<-modloc[i] 
    xcor[i]<-xcoor[2*loc-1] 
}
modloc<-omaord2(modloc,xcor)
for (i in 1:moodinum){
    loc<-modloc[i] 
    modelocx[i]<-(xcoor[2*loc-1]+xcoor[2*loc])/2
    modelocy[i]<-ycoor[2*loc-1]+lift
}
text(modelocx,modelocy,labels=labels,cex=cex)
return(list(modelocx=modelocx,labels=labels))
}





plottree<-function(lst,
plot=T,data=F,crit=NULL,orderrule="distcenter",
modelabel=TRUE,ptext=0,leimat=NULL,symbo=NULL,
info=NULL,infolift=0,infopos=0,infochar=NULL,
xmarginleft=0,xmarginright=0,ymargin=0,
xlim=NULL,ylim=NULL,
col="black",col.axis="black",linecol=rep("black",length(lst$parent)),
pch=21,dimen=NULL,yaxt="s",axes=T,
cex=NULL,nodemag=NULL,linemag=1,cex.axis=1,ylab="",cex.lab=1,
colo=FALSE,paletti=NULL,lowest="dens")
{ 
# create vector verticalPos
# find modes, number of modes, attach vertical position to modes
# position of parent is the mid of positions of children:
# use multitree to find siblings of node and "parent" to fine parent
#
#pch=19: solid circle, pch=20: bullet (smaller circle), 
#pch=21: circle, pch=22: square, 
#pch=23: diamond, pch=24: triangle point-up, 
#pch=25: triangle point down. 

if (colo){
  if (is.null(paletti))
    paletti<-c("red","blue","green",
    "orange","navy","darkgreen",
    "orchid","aquamarine","turquoise",
    "pink","violet","magenta","chocolate","cyan",
    colors()[50:657],colors()[50:657])

  col<-colobary(lst$parent,paletti,modecolo=NULL,modepointer=NULL)
  linecol<-col
}
else col<-rep(col,length(lst$parent))

parent<-lst$parent
level<-lst$level
center<-lst$center
if (is.null(center)){
   nodenum<-length(parent)
   dimen<-length(lst$refe)
   nodenum<-length(lst$parent)
   center<-matrix(1,dimen,nodenum)
}
#      
mut<-multitree(parent)    #create multitree 
roots<-mut$roots
child<-mut$child
sibling<-mut$sibling 

if (is.null(dimen)){
  d<-dim(center)[1]
}
else{
  d<-dimen
}

if (is.null(crit)){
   crit<-rep(0,d)          #order so that 1st closest to origo
   if (d==1) crit<-max(center)
   if (!is.null(lst$refe)) crit<-lst$refe
}
if (orderrule=="distcenter") sibord<-siborder(mut,crit,lst$distcenter)
else sibord<-siborder(mut,crit,center)

mlkm<-moodilkm(parent)
modloc<-mlkm$modloc   
#mlkm$modnodes
modenum<-mlkm$lkm  

lst$center<-center
modelinks<-siborToModor(lst)        #make links in right order

itemnum<-length(parent)    
verticalPos<-matrix(0,itemnum,1)

step<-1/modenum
curloc<-0
for (i in 1:modenum){
   curmode<-modelinks[i]   
   verticalPos[curmode]<-curloc
   curloc<-curloc+step
} 


for (i in 1:modenum){
   curnode<-modloc[i]
   par<-parent[curnode]
   while (par>0){
      #calculate mid of children of par
      #go to the end of sibling list
        chi<-child[par]
        summa<-verticalPos[chi]
        childNum<-1
        while(sibling[chi]>0){
           chi<-sibling[chi]
           summa<-summa+verticalPos[chi]
           childNum<-childNum+1
        }                            
        verticalPos[par]<-summa/childNum
        par<-parent[par]
   }
}

if (lowest=="dens") lowesti<-0 else lowesti<-min(lst$level)
if (is.null(ylim)) ylim<-c(lowesti-ymargin,max(level)+ptext+ymargin)
xlim<-c(min(verticalPos)-xmarginleft,max(verticalPos)+xmarginright)
#axes<-
plot(verticalPos,level,xlab="",ylab=ylab,xlim=xlim,ylim=ylim,xaxt="n",
col=col,col.axis=col.axis,pch=pch,yaxt=yaxt,axes=axes,cex=nodemag,
cex.axis=cex.axis,cex.lab=cex.lab)  

for (i in 1:itemnum){
    if (parent[i]>0){
        xchild<-verticalPos[i]
        ychild<-level[i]
        xparent<-verticalPos[parent[i]]
        yparent<-level[parent[i]]
        segments(xparent,yparent,xchild,ychild,col=linecol[i],lwd=linemag)
     }
}                
#
# lets plot info
#
if (!is.null(info)){
   nodenum<-itemnum
   infolocx<-matrix(nodenum,1)
   infolocy<-matrix(nodenum,1)
   #
   for (i in 1:nodenum){
     infolocx[i]<-verticalPos[i] 
     infolocy[i]<-level[i]+infolift
   }
   digits<-3
   info<-format(info,digits=digits)
   adj<-NULL
   pos<-infopos
   text(infolocx,infolocy,info,pos,adj,cex=cex)       
}
#
# lets plot character info
#
if (!is.null(infochar)){
   nodenum<-itemnum
   infolocx<-matrix(nodenum,1)
   infolocy<-matrix(nodenum,1)
   #
   for (i in 1:nodenum){
     infolocx[i]<-verticalPos[i] 
     infolocy[i]<-level[i]+infolift
   }
   pos<-infopos
   text(infolocx,infolocy,infochar,pos,cex=cex)       
}
#
# lets plot labels for modes
#
if (modelabel){
#
xcoor<-verticalPos
ycoor<-level
#
mlkm<-moodilkm(parent)
modloc<-mlkm$modloc  
modenum<-length(modloc)
modelocx<-matrix(0,modenum,1)
modelocy<-matrix(0,modenum,1)
if (is.null(leimat)){
   if (is.null(symbo)){
       labels<-paste("M",1:modenum,sep="")
   }
   else{
      labels<-paste(symbo,1:modenum,sep="")
   }
}
else{
   labels<-leimat
}
xcor<-matrix(0,modenum,1)                       
for (i in 1:modenum){
    loc<-modloc[i]
    xcor[i]<-xcoor[loc]
}
modloc<-omaord2(modloc,xcor)
for (i in 1:modenum){
    loc<-modloc[i]
    modelocx[i]<-xcoor[loc]
    modelocy[i]<-ycoor[loc]+ptext
}
text(modelocx,modelocy,labels,cex=cex)      
##
}
###############
}









plottwin<-function(tt,et,lev,bary,orde="furthest",ordmet="etaisrec")
{

#if (is.null(et$low)){
   d<-length(et$N)
   step<-matrix(0,d,1)
   for (i in 1:d) step[i]=(et$support[2*i]-et$support[2*i-1])/et$N[i];
   et$step<-step
   et$low<-et$down
   et$upp<-et$high
#}

pp<-plotprof(tt,plot=FALSE,data=TRUE)
vecs<-pp$vecs

d<-length(et$step)

# order the atoms for the level set with level "lev"

lenni<-length(et$value)
distat<-matrix(0,lenni,1)
infopointer<-matrix(0,lenni,1)
lkm<-0
for (i in 1:lenni){
  if (et$value[i]>=lev){
     lkm<-lkm+1
     nod<-i  #nod<-et$nodefinder[i]
     if (ordmet=="etaisrec"){
         recci<-matrix(0,2*d,1)
         for (jj in 1:d){
            recci[2*jj-1]<-et$support[2*jj-1]+et$step[jj]*et$low[nod,jj]
            recci[2*jj]<-et$support[2*jj-1]+et$step[jj]*et$upp[nod,jj]
         }
         distat[lkm]<-etaisrec(bary,recci)
     }
     else{
         lowi<-matrix(0,d,1)
         uppi<-matrix(0,d,1)
         for (jj in 1:d){
            lowi[jj]<-et$support[2*jj-1]+et$step[jj]*et$low[nod,jj]
            uppi[jj]<-et$support[2*jj-1]+et$step[jj]*et$upp[nod,jj]
         }
         baryc<-lowi+(uppi-lowi)/2
         distat[lkm]<-etais(baryc,bary)  #etais(baryc[lk m,],baryind)
     }
     infopointer[lkm]<-i
  }
}
distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   #pointe->et$value,et$nodefinder

ord<-order(distat)
infopointer<-infopointer[ord]

xmin<-et$support[1]
xmax<-et$support[2]
ymin<-et$support[3]
ymax<-et$support[4]
plot(x=bary[1],y=bary[2],xlab="",ylab="",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
pch=20,col="red")

i<-1
while (i<=lkm){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]   #ip<-et$nodefinder[infopointer[node]]

     x1<-et$support[1]+et$step[1]*et$low[ip,1]
     x2<-et$support[1]+et$step[1]*et$upp[ip,1] 
     y1<-et$support[3]+et$step[2]*et$low[ip,2]
     y2<-et$support[3]+et$step[2]*et$upp[ip,2] 
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="lightblue")

     i<-i+1
}

xmin2<-min(vecs[,1])
xmax2<-max(vecs[,3])
ymin2<-0
ymax2<-omamax(vecs[,2])
dev.new()
plot("","",xlab="",ylab="",xlim=c(xmin2,xmax2),ylim=c(ymin2,ymax2))

ycor<-ymax
i<-1
while ((i<=lkm) && (ycor>ymin2)){

     if (orde=="furthest") node<-lkm-i+1 else node<-i
     ip<-infopointer[node]  #ip<-et$nodefinder[infopointer[node]]

     x1<-et$support[1]+et$step[1]*et$low[ip,1]
     x2<-et$support[1]+et$step[1]*et$upp[ip,1] 
     y1<-et$support[3]+et$step[2]*et$low[ip,2]
     y2<-et$support[3]+et$step[2]*et$upp[ip,2] 
     dev.set(which = dev.next())
     polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col="blue")
     points(x=bary[1],y=bary[2],pch=20,col="red")

     ttnode<-node     
     vecci<-vecs[ttnode,]
     x0<-vecci[1]
     y0<-vecci[2]
     x1<-vecci[3]
     y1<-vecci[4] 
     dev.set(which = dev.next())
     segments(x0, y0, x1, y1)

     loc<-locator(1)
     ycor<-loc$y 

     i<-i+1
}



}

plotvecs<-function(vecs,
depths=NULL,segme=T,lift=NULL,
modetest=NULL,alpha=NULL,
axes=TRUE,xlim=NULL,ylim=NULL,xaxt=xaxt,col="black",col.axis="black",
modecolors=NULL,modethickness=1,
leafcolors=NULL,leaflift=0,leafsymbo=20,
modelabels=NULL,ptext=0,
yaxt="s",log="",cex.axis=1)
{
#Plots vectors in vec
#
#vecs is nodenum*4-matrix
#vecs[i,1] x-coordi alulle
#vecs[i,2] y-coordi alulle = vecs[i,4] y-coordi lopulle
#vecs[i,3] x-coordi lopulle
#
#plot(c(1,2),c(3,3))  
#segments(1,3,2,3)
#
#plot(c(1,2,3,4),c(3,3,2,2))  
#segments(1,3,2,3)
#segments(3,2,4,2)
#
#vecs<-matrix(0,3,4)    
#vecs[1,]<-c(1,1,4,1)
#vecs[2,]<-c(5,1,6,1)
#vecs[3,]<-c(2,2,3,2)
#
#plot(c(1,4,5,6,2,3),c(1,1,1,1,2,2))
#segments(1,1,4,1)
#segments(5,1,6,1)
#segments(2,2,3,2)

nodenum<-length(vecs[,1])
xcoor<-matrix(0,2*nodenum,1)
ycoor<-matrix(0,2*nodenum,1)

for (i in 1:nodenum){
 xcoor[2*i-1]<-vecs[i,1]
 xcoor[2*i]<-vecs[i,3]
 ycoor[2*i-1]<-vecs[i,2]
 ycoor[2*i]<-vecs[i,4]
}

#ylim<-c(0,max(ycoor)+ptext)
plot(xcoor,ycoor,xlab="",ylab="",axes=axes,xlim=xlim,ylim=ylim,xaxt=xaxt,
col=col,col.axis=col.axis,yaxt=yaxt,log=log,cex.axis=cex.axis)

if (!is.null(leafcolors)){
   xpoint<-matrix(0,nodenum,1)
   ypoint<-matrix(0,nodenum,1)
   leafcol<-matrix("",nodenum,1)
   zahl<-0
   for (no in 1:nodenum){
      if (leafcolors[no]!="black"){
          zahl<-zahl+1
          xpoint[zahl]<-vecs[no,1]+(vecs[no,3]-vecs[no,1])/2
         
          lif<-(depths[no]-1)*lift
          yc<-ycoor[2*no-1]+lif
          ypoint[zahl]<-yc+leaflift
         
          leafcol[zahl]<-leafcolors[no]
      }
   }
   xpoint<-xpoint[1:zahl]
   ypoint<-ypoint[1:zahl]
   leafcol<-leafcol[1:zahl]
   points(xpoint,ypoint,pch=leafsymbo,col=leafcol)
}

if (!is.null(modelabels)){
   for (no in 1:nodenum){
      if (modelabels[no]!=""){
          xpoint<-vecs[no,1]+(vecs[no,3]-vecs[no,1])/2
         
          lif<-(depths[no]-1)*lift
          yc<-ycoor[2*no-1]+lif
          ypoint<-yc+ptext
         
          label<-modelabels[no]
          
          text(xpoint,ypoint,label)

      }
   }
}

if (segme==T){
 
   thick<-1
   lif<-0
   col<-"black"   

   for (i in 1:nodenum){

        if (!is.null(depths))  lif<-(depths[i]-1)*lift
        if (!is.null(modecolors)){
                if (modecolors[i]!="black") thick=modethickness 
                                            #thick<-2.2^(depths[i]-1)  
                col<-modecolors[i]   
        }
        if (!is.null(modetest)){
             col<-4
             if (modetest[i]>0){
                if (modetest[i]>alpha)  col<-2   
                     #red, hyvaksytaan 0-hypoteesi=ei moodia
                     #nodes with red are not a real feature
                else col<-4   #blue
             } 
       }
           #testcrit<-modetest[i]*qnorm(1-alpha/2)
           #if (excmassa>testcrit)

        yc<-ycoor[2*i-1]+lif
        segments(xcoor[2*i-1],yc,xcoor[2*i],yc,col=col,lwd=thick)
        
        #lines(c(xcoor[2*i-1],xcoor[2*i]),c(ycoor[2*i-1],ycoor[2*i]),col=2) 
  
   }
}

#return(t(tc),t(em))
}











plotvolu2d<-function(vd,theta=NULL,phi=NULL,typer="flat")
{
# typer "dens"/"flat"

if (is.null(phi)) phi<-30

if (vd$type2=="slice"){

if (vd$type=="radius"){

 if (typer=="flat"){
     if (is.null(theta)) theta<-50
     persp(vd$x,vd$y,vd$z,
     xlab="level",ylab="",zlab="radius",ticktype="detailed",
     phi=phi,theta=theta)
 }
 else{
     levnum<-length(vd$x)
     ynumold<-length(vd$y)
     maksi<-max(vd$z)
     gnum<-100
     step<-maksi/(gnum-1)
     xnew<-seq(0,maksi,step)
     znew<-matrix(0,gnum,ynumold)
     for (i in 1:levnum){
        for (j in 1:ynumold){
           highness<-round(gnum*vd$z[i,j]/maksi)
           znew[1:highness,j]<-vd$x[i]  #level[i]
        }
     }
     if (is.null(theta)) theta<-40
     persp(xnew,vd$y,znew,
     xlab="radius",ylab="",zlab="level",
     ticktype="detailed",
     phi=phi,theta=theta)
 }
}

if (vd$type=="proba"){
if (is.null(theta)) theta<--130
if (vd$norma) xlab<-"normalized volume" else xlab<-"volume"
persp(vd$x,vd$y,vd$z,
xlab=xlab,ylab="",zlab="radius",ticktype="detailed",
phi=phi,theta=theta)
}

}   #type2=="slice"

else{ #type2=="boundary"

if (is.null(theta)) theta<-50
persp(vd$x,vd$y,vd$z,
xlab="",ylab="",zlab="level",ticktype="detailed",
phi=phi,theta=theta)

}


}




plotvolu.new<-function(lst,dens=TRUE)
{
mt<-multitree(lst$parent)
itemnum<-length(lst$volume)
rootnum<-length(mt$roots)
left<-mt$child
right<-mt$sibling
vecs<-matrix(0,itemnum,3)

sibord<-mt$siborder  #siborder.new(mt)

# allocate space for roots

rootsvolume<-0
for (i in 1:rootnum){
  now<-mt$roots[i]
  rootsvolume<-rootsvolume+lst$volume[now]
}
basis<-rootsvolume+rootsvolume/4
gaplen<-(basis-rootsvolume)/(rootnum+1)
rootlinks<-matrix(0,rootnum,1)  #make links in right order
{
if (rootnum==1){ 
  rootlinks[1]<-mt$roots[1]  #1
}
else{ 
     for (i in 1:rootnum){
         now<-mt$roots[i]
         roor<-sibord[now]
         rootlinks[roor]<-now
     }
}
xbeg<-0
xend<-0
for (i in 1:rootnum){
  now<-rootlinks[i]
  ycoo<-lst$level[now]
  xend<-xbeg+lst$volume[now]
  vecs[now,]<-c(xbeg,xend,ycoo)
  xbeg<-gaplen+xend
}
}
# allocate space for nonroots

for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    pino[1]<-mt$roots[i]  
    pinin<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        if (left[cur]>0){     #if not leaf (root may be leaf)
           vecs<-allokoi.new(cur,vecs,lst,left,right,sibord)   
        }
        if (right[cur]>0){    #if right exists, put to stack
            pinin<-pinin+1
            pino[pinin]<-right[cur]
        }
        while (left[cur]>0){    #go to leaf and put right nodes to stack
             cur<-left[cur]
             if (left[cur]>0){  #if not leaf
                vecs<-allokoi.new(cur,vecs,lst,left,right,sibord)
             }
             if (right[cur]>0){ #if right exists, put to stack
                pinin<-pinin+1
                pino[pinin]<-right[cur]
             }
        }
    }
} 
      
if (dens) firstlevel<-0 else firstlevel<-min(lst$level)
xlim<-c(min(vecs[,1]),max(vecs[,2]))
ylim<-c(firstlevel,max(lst$level))
plot(x="",y="",xlab="",ylab="",xlim=xlim,ylim=ylim)

for (i in 1:itemnum){

    yc<-vecs[i,3]

    pare<-lst$parent[i]
    if (pare==0) lowlev<-firstlevel else lowlev<-lst$level[pare]

    segments(vecs[i,1],lowlev,vecs[i,1],yc)#,col=col,lwd=thick)
    segments(vecs[i,2],lowlev,vecs[i,2],yc)#,col=col,lwd=thick)

    if (left[i]==0){  #we are in leaf

       segments(vecs[i,1],yc,vecs[i,2],yc)#,col=col,lwd=thick)

    }
    else{

       childnum<-1
       curchi<-mt$child[i]
       while (mt$sibling[curchi]!=0){
           curchi<-mt$sibling[curchi]
           childnum<-childnum+1
       }

       sibpointer<-matrix(0,childnum,1)
       curchi<-mt$child[i]
       sibpointer[sibord[curchi]]<-curchi
       while (mt$sibling[curchi]!=0){
           curchi<-mt$sibling[curchi]
           sibpointer[sibord[curchi]]<-curchi
       }

       curchi<-sibpointer[1]
       x1<-vecs[curchi,1]      
       segments(vecs[i,1],yc,x1,yc)#,col=col,lwd=thick)
       x0<-vecs[curchi,2]

       cn<-2
       while (cn<=childnum){
             curchi<-sibpointer[cn]
             x1<-vecs[curchi,1] 
             segments(x0,yc,x1,yc)#,col=col,lwd=thick)
             x0<-vecs[curchi,2] 
             cn<-cn+1
       }

       segments(x0,yc,vecs[i,2],yc)#,col=col,lwd=thick)

    }
}


}





plotvolu<-function(lst,length=NULL,
toplot=TRUE,data=FALSE,crit=NULL,orderrule="distcenter",
modelabel=FALSE,ptext=0,leimat=NULL,symbo=NULL,
info=NULL,infolift=0,infopos=0,
xmarginleft=0,xmarginright=0,ymargin=0,
xlim=NULL,ylim=NULL,
col="black",col.axis="black",
cutlev=NULL,xaxt="s",yaxt="s",
exmavisu=NULL,bg="transparent",tyyppi="n",
lty="solid",colo=FALSE,lowest="dens",proba=FALSE,
paletti=NULL,cex=NULL,modecolo=NULL,modepointer=NULL,upper=TRUE,
cex.axis=1,xlab="",ylab="",cex.lab=1,colothre=NULL,nodes=NULL)
{
if (upper) firstlevel<-min(lst$level) else firstlevel<-max(lst$level)
if (lowest=="dens") firstlevel<-0

parents<-lst$parent
levels<-lst$level
length<-lst$volume
if (proba) length<-lst$proba
center<-lst$center

mut<-multitree(parents)
if (is.null(lst$roots)) roots<-mut$roots else roots<-lst$roots
child<-mut$child
sibling<-mut$sibling

d<-dim(center)[1]
if (is.null(crit)){
   crit<-rep(0,d)          #order so that 1st closest to origo
   if (d==1) crit<-max(center)
   if (!is.null(lst$refe)) crit<-lst$refe
}

if (orderrule=="distcenter") sibord<-siborder(mut,crit,lst$distcenter)
else sibord<-siborder(mut,crit,center)

itemnum<-length(parents)
vecs<-matrix(NA,itemnum,4)
vecs<-alloroot(vecs,roots,sibord,levels,length)
vecs<-plotdata(roots,child,sibling,sibord,levels,length,vecs)
orivecs<-vecs

if (!(is.null(cutlev))){
  cm<-cutmut(mut,cutlev,levels)              # cut the tree
  roots<-cm$roots
  sibling<-cm$sibling
  mut$roots<-roots
  if (orderrule=="distcenter") sibord<-siborder(mut,crit,lst$distcenter)
  else sibord<-siborder(mut,crit,center)
  rootnum<-length(roots) 
  apuvecs<-matrix(NA,itemnum,4)
  for (i in 1:rootnum){
     inde<-roots[i]
     apuvecs[inde,]<-vecs[inde,]
     if (i==1) miniroot<-apuvecs[inde,1]
     else if (apuvecs[inde,1]<=miniroot) miniroot<-apuvecs[inde,1]
  }
  vecs<-apuvecs          #we give for the roots the previous positions
  vecs<-plotdata(roots,child,sibling,sibord,levels,length,vecs)
}

#####################################

depths<-NULL
segme<-T
lift<-NULL
modetest<-NULL
alpha<-NULL
axes<-T
modecolors<-NULL
modethickness<-1
leafcolors<-NULL
leaflift<-0
leafsymbo<-20
modelabels<-NULL
log<-""

nodenum<-length(vecs[,1])
xcoor<-matrix(0,2*nodenum,1)
ycoor<-matrix(0,2*nodenum,1)

for (i in 1:nodenum){
 xcoor[2*i-1]<-vecs[i,1]
 xcoor[2*i]<-vecs[i,3]
 ycoor[2*i-1]<-vecs[i,2]
 ycoor[2*i]<-vecs[i,4]
}

oriminnu<-min(orivecs[,1],na.rm=T)
minnu<-min(xcoor,na.rm=T)
if (is.null(cutlev)) xcoor<-xcoor-minnu
else xcoor<-xcoor-oriminnu

if (lowest=="dens") lowesti<-0 else lowesti<-min(lst$level)
#xlim<-c(min(vecs[,1],na.rm=T)-xmarginleft,max(vecs[,3],na.rm=T)+xmarginright)
if (is.null(ylim)){
    ylim<-c(lowesti,max(ycoor,na.rm=T)+ptext+ymargin)
    if (!is.null(cutlev)) 
    ylim<-c(cutlev,max(ycoor,na.rm=T)+ptext+ymargin)
}

if (toplot){
par(bg=bg)
plot(xcoor[order(xcoor)],ycoor[order(xcoor)],  #xcoor,ycoor,
xlab=xlab,ylab=ylab,axes=axes,xlim=xlim,ylim=ylim,xaxt=xaxt,
col=col,col.axis=col.axis,yaxt=yaxt,log=log,
type=tyyppi,lty=lty,cex.axis=cex.axis,cex.lab=cex.lab)
}
###########################################################

if ((tyyppi=="n") && (toplot)){

thick<-1
col<-col #"black"

for (i in 1:nodenum){
if (!is.na(ycoor[2*i-1])){

    yc<-ycoor[2*i-1]

    pare<-parents[i]
    if (pare==0) lowlev<-firstlevel else lowlev<-levels[pare]

    segments(xcoor[2*i-1],lowlev,xcoor[2*i-1],yc,col=col,lwd=thick)
    segments(xcoor[2*i],lowlev,xcoor[2*i],yc,col=col,lwd=thick)

    if (child[i]==0){  #we are in leaf

       segments(xcoor[2*i-1],yc,xcoor[2*i],yc,col=col,lwd=thick)

    }
    else{

       yc<-ycoor[2*i-1]

       childnum<-1
       curchi<-child[i]
       while (sibling[curchi]!=0){
           curchi<-sibling[curchi]
           childnum<-childnum+1
       }

       sibpointer<-matrix(0,childnum,1)
       curchi<-child[i]
       sibpointer[sibord[curchi]]<-curchi
       while (sibling[curchi]!=0){
           curchi<-sibling[curchi]
           sibpointer[sibord[curchi]]<-curchi
       }

       curchi<-sibpointer[1]
       x1<-xcoor[2*curchi-1]      
       segments(xcoor[2*i-1],yc,x1,yc,col=col,lwd=thick)
       x0<-xcoor[2*curchi] 

       cn<-2
       while (cn<=childnum){
             curchi<-sibpointer[cn]
             x1<-xcoor[2*curchi-1] 
             segments(x0,yc,x1,yc,col=col,lwd=thick)
             x0<-xcoor[2*curchi] 
             cn<-cn+1
       }

       segments(x0,yc,xcoor[2*i],yc,col=col,lwd=thick)

    }
}
}

for (i in 1:nodenum){
   if (is.null(cutlev)){
     orivecs[i,1]<-orivecs[i,1]-minnu
     orivecs[i,3]<-orivecs[i,3]-minnu
   }
   else{
     orivecs[i,1]<-orivecs[i,1]-oriminnu
     orivecs[i,3]<-orivecs[i,3]-oriminnu
   }
}   
if (modelabel) 
modelab<-plottext(parents,orivecs,ptext,leimat,symbo=symbo,cex=cex)  


}  #tyyppi = "n"


if (!is.null(lst$predictor.node)) 
segments(
xcoor[2*lst$predictor.node-1],
ycoor[2*lst$predictor.node-1],
xcoor[2*lst$predictor.node],
ycoor[2*lst$predictor.node])



############################################# exmavisu start

if (colo) exmavisu<-roots #1

if (!is.null(exmavisu)){

if (colo){
  if (is.null(paletti))
    paletti<-c("red","blue","green",
    "orange","navy","darkgreen",
    "orchid","aquamarine","turquoise",
    "pink","violet","magenta","chocolate","cyan",
    colors()[50:657],colors()[50:657])

  col<-colobary(lst$parent,paletti,modecolo=modecolo,modepointer=modepointer)

  if (!is.null(colothre))
  col<-colobary.merge(lst$parent,lst$level,colothre,paletti)
  if (!is.null(nodes))
  col<-colobary.nodes(lst$parent,nodes,paletti)

}
else col<-rep("blue",length(lst$parent))

for (i in 1:length(exmavisu)){

node<-exmavisu[i]

x1<-xcoor[2*node-1] 
x2<-xcoor[2*node]
lev<-levels[node]
if (parents[node]>0) lev0<-levels[parents[node]] else lev0<-firstlevel
polygon(c(x1,x2,x2,x1),c(lev0,lev0,lev,lev),col=col[node],lty="blank")

pino<-matrix(0,nodenum,1)
pino[1]<-child[node]
if (child[node]>0) pinoin<-1 else pinoin<-0

while (pinoin>0){
   node<-pino[pinoin]
   pinoin<-pinoin-1   

   x1<-xcoor[2*node-1] 
   x2<-xcoor[2*node]
   lev<-levels[node]
   if (parents[node]>0) lev0<-levels[parents[node]] else lev0<-firstlevel
   polygon(c(x1,x2,x2,x1),c(lev0,lev0,lev,lev),col=col[node],lty="blank")

   if (sibling[node]>0){
         pinoin<-pinoin+1
         pino[pinoin]<-sibling[node] 
   }

   while (child[node]>0){    #go to left and put right nodes to stack
         node<-child[node]

         x1<-xcoor[2*node-1] 
         x2<-xcoor[2*node]
         lev<-levels[node]
         if (parents[node]>0) lev0<-levels[parents[node]] else lev0<-firstlevel
         polygon(c(x1,x2,x2,x1),c(lev0,lev0,lev,lev),col=col[node],lty="blank")

         if (sibling[node]>0){
            pinoin<-pinoin+1
            pino[pinoin]<-sibling[node] 
         }
   }
}
}
}
####################### exmavisu end

if (data) return(list(xcoor=xcoor,ycoor=ycoor))


}





point.eval<-function(tr,x)
{
# tr is an evaluation tree

d<-length(tr$support)/2
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(tr$support[2*i]-tr$support[2*i-1])/tr$N[i]

# if x is not in the support, then ans=0
insupport<-1
for (i in 1:d){
    if ((x[i]>=tr$support[2*i]) || (x[i]<=tr$support[2*i-1])){
       ans<-0
       insupport<-0
    }
}
if (insupport==1){
  node<-1
  while (tr$left[node]>0){
      dir<-tr$direc[node]
      spl<-tr$split[node]
      realspl<-tr$support[2*dir-1]+spl*step[dir]
      if (x[dir]>realspl) node<-tr$right[node]
      else node<-tr$left[node]
  }
  #loc<-tr$infopointer[node]
  #ans<-tr$value[loc]
  ans<-tr$mean[node]
}

return(ans)
}
posipart<-function(pcf)
{
pcf$value<-pmax(pcf$value,0)
return(pcf)
}

pp.plot<-function(dendat=NULL,compa="gauss",basis="gauss",mean=0,sig=1,df=1,
gnum=1000,d=1,R=3,pptype="1d",cex.lab=1,cex.axis=1,col="blue",lwd=1)
# basis is either data (dendat) or a theoretical distribution
{
if (pptype=="1d"){
   p<-dendat[order(dendat)]
   if (compa=="gauss") y<-pnorm(p,mean=mean,sd=sig)
   if (compa=="student") y<-pt((p-mean)/sig,df=df)
   if (compa=="unif") y<-punif((p-mean)/sig)
   if (compa=="exp") y<-pexp((p-mean)/sig)
   if (compa=="doubleexp") 
      y<-0.5*(1-pexp(-(p-mean)/sig))+0.5*pexp((p-mean)/sig)
   n<-length(dendat) #dim(dendat)[1]
   x<-seq(1:n)/n
   tyyppi<-"p"
   xlab<-"empirical distribution function"
   ylab<-"compared distribution function"
}
if (pptype=="v2p"){
      rp<-tailfunc(R,d,type=compa,gnum=gnum,sig=sig,nu=df)
      y<-rp$proba
      rp2<-tailfunc(R,d,type=basis,gnum=gnum,sig=sig,nu=df)
      x<-rp2$proba
      tyyppi="l"
      xlab<-"empirical"
      ylab<-"model"
}
if (pptype=="ddplot"){
}

plot(x,y,
type=tyyppi,
xlim=c(0,1),ylim=c(0,1),
xlab=xlab,ylab=ylab,cex.lab=cex.lab,cex.axis=cex.axis)
segments(0,0,1,1,col=col,lwd=lwd)
}

preprocess<-function(dendat, type="copula")
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]
prodendat<-matrix(0,n,d)

if (type=="sphering"){

   cova<-cov(dendat)
   eig<-eigen(cova,symmetric=TRUE)
   sigsqm<-eig$vectors%*%diag(eig$values^(-1/2)) 
   prodendat<-t(t(sigsqm)%*%t(dendat-mean(dendat)))   # dendat%*%sigsqm 

}
else if (type=="sd"){
   for (ii in 1:d){
        prodendat[,ii]<-(dendat[,ii]-mean(dendat[,ii]))/sd(dendat[,ii])
   }

}
else if (type=="standardcopula"){
   for (ii in 1:d){
        or<-order(dendat[,ii])
        mones<-matrix(0,n,1)
        for (i in 1:n) mones[or[i]]<-i
        prodendat[,ii]<-mones/n
   }

}
else{
   for (ii in 1:d){
        or<-order(dendat[,ii])
        mones<-matrix(0,n,1)
        for (i in 1:n) mones[or[i]]<-i
        prodendat[,ii]<-qnorm(mones/n)
   }
}

return(prodendat)
}
prof2vecs<-function(profile,level,n=NULL,crit,motes=NULL){

parents<-profile$parent
nodenum<-length(parents)
centers<-profile$center
 
nodenum<-length(parents)   
levels<-matrix(level,nodenum,1) #all will be plotted at same lev(=logh)
excma<-excmas(profile)       #instead of volumes, we use excesss mass
                             #to determine the lengths of the vectors
#motes<-mtest(profile,n)

mut<-multitree(parents)

# let us make a vector where modes are labelled with the order, others=0
# later we handle "mlabel" similarily as "motes"
mlabel<-matrix(0,nodenum,1)
mlkm<-moodilkm(parents)      #mlkm$lkm, mlkm$modloc 
for (run in 1:mlkm$lkm){
   alku<-mlkm$modloc[run]
   while ((parents[alku]>0) && 
          (mut$sibling[mut$child[parents[alku]]]==0)){
      alku<-parents[alku]
   }
   mlabel[alku]<-run  
}

mt<-pruneprof(mut)
depths<-depth(mt)
roots<-mt$roots
child<-mt$child
sibling<-mt$sibling

sibord<-siborder(mt,crit,centers)

itemnum<-length(parents)
vecs<-matrix(NA,itemnum,4)  
vecs<-alloroot(vecs,roots,sibord,levels,excma) 
vecs<-plotdata(roots,child,sibling,sibord,levels,excma,vecs)
vecnum<-length(vecs[,1])      #vecs has four columns

#  remove pruned

if (is.null(motes)) motes<-matrix(0,vecnum,1)

tempvecs<-matrix(0,vecnum,4)
tempdepths<-matrix(0,vecnum,1)
tempmotes<-matrix(0,vecnum,1)
tempmlabel<-matrix(0,vecnum,1)
ind<-0
for (i in 1:vecnum){
       if (!(is.na(vecs[i,1]))){
             ind<-ind+1
             tempvecs[ind,]<-vecs[i,]
             tempdepths[ind]<-depths[i]
             tempmotes[ind]<-motes[i]
             tempmlabel[ind]<-mlabel[i]
         }
}
vecs<-tempvecs[1:ind,]
depths<-tempdepths[1:ind]
motes<-tempmotes[1:ind]
mlabel<-tempmlabel[1:ind]

return(list(vecs=vecs,depths=depths,motes=motes,mlabel=mlabel))
}                        





profgene<-function(values,recs,frekv=NULL,cvol=TRUE,ccen=TRUE,cfre=FALSE,
outlsets=TRUE,invalue=TRUE)
{

cu<-cumu(values,recs,frekv)
levels<-cu$levels
lsets<-cu$lsets
atoms<-cu$atoms
binfrek<-cu$frekv  #kullekin suorakaiteelle frekvenssi

alkublokki<-200
blokki<-50
links<-toucrec(atoms,alkublokki,blokki)

alkublokki2<-200
blokki2<-50
dentree<-decom(lsets,levels,links,alkublokki2,blokki2)
seplsets<-dentree$lsets
sepval<-dentree$levels
parents<-dentree$parents

if (cfre) nodefrek<-cfrekv(seplsets,binfrek) else nodefrek<-NULL 

if (ccen==TRUE) cvol<-TRUE
if (cvol){
  volum<-cvolum(seplsets,atoms)
  kerroin<-cinte(sepval,volum,parents) 
  sepvalnor<-sepval/kerroin
} 
else{  
  volum<-NULL
  sepvalnor<-NULL
}

if (ccen && cvol) centers<-ccente(seplsets,atoms,volum) else centers<-NULL

if (!(outlsets)) seplsets<-NULL
if (!(invalue)) sepval<-NULL

return(list(parent=parents,level=sepvalnor,invalue=sepval,
volume=volum,center=centers,nodefrek=nodefrek,lsets=seplsets))
#values: normeeratut arvot
#invalues: alkuperaiset frekvenssit/arvot 
#nodefrek: kunkin solmun frekvenssi
}

















profhist<-function(dendat,binlkm,cvol=TRUE,ccen=TRUE,cfre=FALSE)
{
#esim. dendat<-matrix(rnorm(20),10) on 10*2 matriisi

epsi<-0
hi<-histo(dendat,binlkm,epsi)
recs<-hi$recs
hisfrekv<-hi$values

pr<-profgene(values=hisfrekv,recs=recs,frekv=hisfrekv,cvol=cvol,ccen=ccen,
cfre=cfre)

return(list(parent=pr$parent,level=pr$level,invalue=pr$invalue,
volume=pr$volum,center=pr$center,nodefrek=pr$nodefrek,recs=recs,
hisfrekv=t(hisfrekv),lsets=pr$lsets))
}










profkernC<-function(dendat,h,N,Q,cvol=TRUE,ccen=TRUE,#cfre=FALSE,
numofallR=10000){

#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(8,8)
#Q<-3

n<-dim(dendat)[1]
d<-length(N)
hnum<-length(h)
mnn<-maxnodenum(dendat,h,N,n,d)
extMaxnode<-mnn$maxnode
extMaxvals<-mnn$maxpositive

if (hnum>1){
 inh<-matrix(0,hnum+1,1)
 inh[2:(hnum+1)]<-h
}
else{
 inh<-h
}
inN<-matrix(0,d+1,1)
inN[2:(d+1)]<-N

dentree<-.C("kerprofC",as.integer(extMaxnode),
                  as.integer(extMaxvals),
                  as.double(dendat),
                  as.double(inh),
                  as.integer(inN),
                  as.integer(n),
                  as.integer(hnum),
                  as.integer(d),
                  as.integer(Q),
                  as.integer(numofallR),
                  level = double(numofallR+1),
                  parent = integer(numofallR+1),
                  component = integer(numofallR+1),
                  volume = double(numofallR+1),
                  center = double(d*numofallR+1),
                  efek = integer(1),
PACKAGE="denpro")

invalue<-dentree$level[2:(dentree$efek+1)]
parent<-dentree$parent[2:(dentree$efek+1)]
volume<-dentree$volume[2:(dentree$efek+1)]
kerroin<-cinte(invalue,volume,parent) 
sepvalnor<-invalue/kerroin
veccenter<-dentree$center[2:(d*dentree$efek+1)]
center<-matrix(0,dentree$efek,d)
for (i in 1:dentree$efek){
  for (j in 1:d){
     center[i,j]<-veccenter[(i-1)*d+j]
  }
}
center<-t(center)

#if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL 

#clus<-F
#if (clus){
#   clustervecs<-cluskern(compo,component,AtomlistAtom,AtomlistNext,kg,dendat,
#   h,N)
#}
#else{
#   clustervecs<-NULL
#}

return(list(parent=parent,level=sepvalnor,invalue=invalue,
volume=volume,center=center))
#,nodefrek=nodefrek,clustervec=clustervecs))
#
#values: normeeratut arvot
#invalues: normeeraamattomat arvot 
#nodefrek: kunkin solmun frekvenssi
}









profkernCRC<-function(dendat,h,N,Q,cvol=TRUE,ccen=TRUE,#cfre=FALSE,
kernel="epane",compoinfo=FALSE,trunc=3,threshold=0.0000001,katka=NULL,hw=NULL)
{
#dyn.load("/home/jsk/kerle/kerleCversio")
#pk2<-profkernCRC(dendat,h,N,Q)
#
#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(8,8)
#Q<-3
#
n<-dim(dendat)[1]
d<-length(N)

if (is.null(hw)) weig<-rep(1/n,n) 
else{
   weig<-weightsit(n,hw)

   dendatnew<-dendat
   weignew<-weig
   cumul<-0
   for (i in 1:n){
        if (weig[i]>0){
            cumul<-cumul+1
            dendatnew[cumul,]<-dendat[i,]
            weignew[cumul]<-weig[i] 
        }
   }
   dendat<-dendatnew[1:cumul,]
   weig<-weignew[1:cumul]
   n<-cumul
}
inweig<-matrix(0,n+1,1)
inweig[2:(n+1)]<-weig

hnum<-length(h)
mnn<-maxnodenum(dendat,h,N,n,d)
extMaxnode<-mnn$maxnode
extMaxvals<-mnn$maxpositive
#
if (hnum>1){
 inh<-matrix(0,hnum+1,1)
 inh[2:(hnum+1)]<-h
}
else{
 inh<-h
}
inN<-matrix(0,d+1,1)
inN[2:(d+1)]<-N

if (kernel=="epane") kertype<-1
else kertype<-2  # gaussian

kg<-.C("kergrid",
               as.integer(extMaxnode),
               as.integer(extMaxvals),
               as.double(dendat),
               as.double(inh),
               as.integer(inN),
               as.integer(n),
               as.integer(hnum),
               as.integer(d),
               as.integer(kertype),
               as.double(trunc),
               as.double(threshold),  
               as.double(inweig),        
               ioleft = integer(extMaxnode+1),
               ioright = integer(extMaxnode+1),
               ioparent = integer(extMaxnode+1),
               infopointer = integer(extMaxnode+1),
               iolow = integer(extMaxnode+1),
               ioupp = integer(extMaxnode+1),
               value = double(hnum*extMaxvals),
               index = integer(d*extMaxvals),
               nodefinder = integer(extMaxvals),
               numpositive = integer(1),
               numnode = integer(1),
PACKAGE="denpro")

left<-kg$ioleft[2:(kg$numnode+1)]
right<-kg$ioright[2:(kg$numnode+1)]
parent<-kg$ioparent[2:(kg$numnode+1)]
infopointer<-kg$infopointer[2:(kg$numnode+1)]
iolow<-kg$iolow[2:(kg$numnode+1)]
ioupp<-kg$ioupp[2:(kg$numnode+1)]

value<-kg$value[2:(kg$numpositive+1)]
nodefinder<-kg$nodefinder[2:(kg$numpositive+1)]
vecindex<-kg$index[2:(d*kg$numpositive+1)]
index<-matrix(0,kg$numpositive,d)
for (i in 1:kg$numpositive){
  for (j in 1:d){
     index[i,j]<-vecindex[(i-1)*d+j]
  }
}

nodenumOfDyaker<-length(left)

maxval<-max(value)
minval<-min(value)
step<-(maxval-minval)/Q
levseq<-seq(from=minval,to=maxval-step,by=step)

levfrekv<-matrix(0,Q,1)
atomnum<-length(value)
for (i in 1:atomnum){
   for (j in 1:Q){
       if (value[i]>=levseq[j]){
          levfrekv[j]<-levfrekv[j]+1
       }
   }
}
numofall<-sum(levfrekv)
levnum<-length(levseq)
    
inlevseq<-matrix(0,length(levseq)+1,1)
inlevseq[2:(length(levseq)+1)]<-levseq
inN<-matrix(0,d+1,1)
inN[2:(d+1)]<-N
inleft<-matrix(0,length(left)+1,1)
inleft[2:(length(left)+1)]<-left
inright<-matrix(0,length(left)+1,1)
inright[2:(length(left)+1)]<-right
inparent<-matrix(0,length(left)+1,1)
inparent[2:(length(left)+1)]<-parent
invalue<-matrix(0,length(value)+1,1)
invalue[2:(length(value)+1)]<-value
#inindex<-matrix(0,dim(kg$index)[1]+1,dim(kg$index)[2]+1)
#for (i in 1:dim(kg$index)[1]){
#  inindex[i+1,]<-c(0,kg$index[i,])
#}
innodefinder<-matrix(0,length(nodefinder)+1,1)
innodefinder[2:(length(nodefinder)+1)]<-nodefinder

# Tama koodi on jo kergrid:ssa, lasketaan volume of one atom
minim<-matrix(0,d,1)  #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-min(dendat[,i])
  maxim[i]<-max(dendat[,i])
}
delta<-(maxim-minim+2*h)/(N+1)  
volofatom<-prod(delta)

inminim<-matrix(0,d+1,1)
inminim[2:(d+1)]<-minim
indelta<-matrix(0,d+1,1)
indelta[2:(d+1)]<-delta

if (!is.null(katka)){
   invalue2<-invalue
   lenni<-length(invalue)
   for (i in 1:lenni){
      if (invalue[i]>=katka) invalue2[i]<-katka
   }
   invalue<-invalue2
}

dentree<-.C("decomdyaC",
               as.integer(numofall),
               as.integer(atomnum),
               as.double(inlevseq),
               as.integer(inN),
               as.integer(d),         
               as.integer(levnum),   
               as.double(volofatom),
               as.double(inminim),
               as.double(h),
               as.double(indelta),
               as.integer(nodenumOfDyaker),
               as.integer(inleft),
               as.integer(inright),
               as.integer(inparent), 
               as.double(invalue),
               as.integer(index),
               as.integer(innodefinder),
               level = double(numofall+1),
               parent = integer(numofall+1),
               component = integer(numofall+1),
               volume = double(numofall+1),
               center = double(d*numofall+1),
               efek = integer(1),
               AtomlistAtom = integer(numofall+1),
               AtomlistNext = integer(numofall+1),
               numOfAtoms = integer(1),
PACKAGE="denpro")

AtomlistAtom<-dentree$AtomlistAtom[2:(dentree$numOfAtoms+1)]
AtomlistNext<-dentree$AtomlistNext[2:(dentree$numOfAtoms+1)]

invalue<-dentree$level[2:(dentree$efek+1)]
parent<-dentree$parent[2:(dentree$efek+1)]
volume<-dentree$volume[2:(dentree$efek+1)]
component<-dentree$component[2:(dentree$efek+1)]
kerroin<-cinte(invalue,volume,parent) 
sepvalnor<-invalue/kerroin
veccenter<-dentree$center[2:(d*dentree$efek+1)]
center<-matrix(0,dentree$efek,d)
for (i in 1:dentree$efek){
  for (j in 1:d){
     center[i,j]<-veccenter[(i-1)*d+j]
  }
}
center<-t(center)

#if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL 

#clus<-F
#if (clus){
#   clustervecs<-cluskern(compo,component,AtomlistAtom,AtomlistNext,kg,dendat,
#   h,N)
#}
#else{
#   clustervecs<-NULL
#}

if (compoinfo)

  return(list(parent=parent,level=sepvalnor,invalue=invalue,
  volume=volume,center=center,
  component=component,
  AtomlistAtom=AtomlistAtom,AtomlistNext=AtomlistNext,index=index))

else

  return(list(parent=parent,level=sepvalnor,invalue=invalue,
  volume=volume,center=center,n=n))

#,nodefrek=nodefrek,clustervec=clustervecs))
#
#values: normeeratut arvot
#invalues: normeeraamattomat arvot 
#nodefrek: kunkin solmun frekvenssi

}












profkern<-function(dendat,h,N,Q,cvol=TRUE,ccen=TRUE,cfre=FALSE,kernel="epane",
compoinfo=FALSE,trunc=3,threshold=0.0000001,sorsa="crc",hw=NULL)
{

if (kernel=="gauss") h<-h*trunc   #trunc<-3

hnum<-length(h)
hrun<-1
while (hrun<=hnum){
   hcur<-h[hrun]

   if (sorsa=="crc")
   curtree<-profkernCRC(dendat,hcur,N,Q,kernel=kernel,compoinfo=compoinfo,
            trunc=trunc,threshold=threshold,hw=hw)
   else
   curtree<-profkernC(dendat,hcur,N,Q)

   if (hrun==1){
      if (hnum==1){
          treelist<-curtree
      }
      else{
          treelist=list(curtree)
      }
   }
   else{
      treelist=c(treelist,list(curtree))
   }
   hrun<-hrun+1
}
#
return(treelist)
}
profkernR<-function(kg,dendat,h,N,Q,frekv=NULL,cvol=TRUE,ccen=TRUE,cfre=FALSE){

#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(4,4)
#Q<-3

#kg<-kergrid(dendat,h,N)

nodenumOfDyaker<-length(kg$left)

value<-kg$value
maxval<-max(value)
minval<-min(value)
step<-(maxval-minval)/Q
levseq<-seq(from=minval,to=maxval-step,by=step)

levfrekv<-matrix(0,Q,1)
atomnum<-length(value)
for (i in 1:atomnum){
   for (j in 1:Q){
       if (value[i]>=levseq[j]){
          levfrekv[j]<-levfrekv[j]+1
       }
   }
}
numofall<-sum(levfrekv)
    
dentree<-decomdya(numofall,atomnum,levseq,kg,N,nodenumOfDyaker)
invalue<-dentree$level
parent<-dentree$parent
component<-dentree$component
AtomlistAtom<-dentree$AtomlistAtom
AtomlistNext<-dentree$AtomlistNext

# Tama koodi on jo kergrid:ssa, lasketaan volume of one atom
d<-length(N)
minim<-matrix(0,d,1)  #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-min(dendat[,i])
  maxim[i]<-max(dendat[,i])
}
delta<-(maxim-minim+2*h)/(N+1)  
volofatom<-prod(delta)

#if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL 

if (ccen==TRUE) cvol<-TRUE
if (cvol){
  volume<-cvolumdya(volofatom,component,AtomlistNext)
  kerroin<-cinte(invalue,volume,parent) 
  sepvalnor<-invalue/kerroin
} 
else{  
  volume<-NULL
  sepvalnor<-NULL
}

if (ccen && cvol){
  index<-kg$index
  d<-dim(dendat)[2]
  center<-ccentedya(volofatom,component,AtomlistNext,AtomlistAtom,
                    volume,minim,h,delta,index,d)
  }
  else{
      center<-NULL
  }

return(list(parent=parent,level=sepvalnor,invalue=invalue,
volume=volume,center=center))#,nodefrek=nodefrek))
#values: normeeratut arvot
#invalues: alkuperaiset frekvenssit/arvot 
#nodefrek: kunkin solmun frekvenssi
}








proftree<-function(tr,
Q=NULL,frekv=NULL,cvol=TRUE,ccen=TRUE,cfre=FALSE)
{

d<-dim(tr$upp)[2]

if (tr$left[1]==0){
  parent=c(0)
  sepvalnor=c(tr$mean[1])
  invalue=c(tr$mean[1])
  volume=c(tr$volume[1])
  rec<-matrix(0,2*d,1)
  for (j in 1:d){
     rec[2*j-1]<-tr$suppo[2*j-1]+tr$low[1,j]*tr$step[j]
     rec[2*j]<-  tr$suppo[2*j-1]+tr$upp[1,j]*tr$step[j]
  }
  center=t(cenone(rec))
}

else{

nodenumOfTree<-length(tr$left)

# make parent
parent<-makeparent(tr$left,tr$right)

mi<-makeinfo(tr$left,tr$right,tr$mean,tr$low,tr$upp)
#infopointer<-mi$infopointer
#terminalnum<-mi$terminalnum
#low<-mi$low
#upp<-mi$upp
#nodefinder<-mi$nodefinder
#value<-mi$value

{
if (!is.null(Q)){
   maxval<-max(mi$value)
   minval<-min(mi$value)
   step<-(maxval-minval)/Q
   levseq<-seq(from=minval,to=maxval-step,by=step)
}
else{
   eppsi<-0        #0.0000001
   levseq<-matrix(0,length(mi$value),1)
   ordu<-order(mi$value)
   ru<-1
   laskuri<-1
   car<-ordu[ru]
   levseq[laskuri]<-mi$value[car]-eppsi
   while (ru < length(mi$value)){
       carnew<-ordu[ru+1]
       if (mi$value[carnew]>mi$value[car]){
          laskuri<-laskuri+1
          levseq[laskuri]<-mi$value[carnew]-eppsi
      }
      ru<-ru+1
   }
   levseq<-levseq[1:laskuri]
   Q<-laskuri
}
}

levfrekv<-matrix(0,Q,1)
atomnum<-length(mi$value)   #=mi$terminalnum
for (i in 1:atomnum){
   for (j in 1:Q){
      if (mi$value[i]>=levseq[j]){
         levfrekv[j]<-levfrekv[j]+1
      }
   }
}
numofall<-sum(levfrekv)

inlevseq<-matrix(0,Q+1,1)
inlevseq[2:(Q+1)]<-levseq
insuppo<-matrix(0,2*d+1,1)
insuppo[2:(2*d+1)]<-tr$suppo
instep<-matrix(0,d+1,1)
sc<-matrix(0,d,1)
for (i in 1:d){
    step[i]<-(tr$support[2*i]-tr$support[2*i-1])/tr$N[i]
}
instep[2:(d+1)]<-sc    #stepcalc(tr$support,tr$N)    #tr$step
inleft<-matrix(0,nodenumOfTree+1,1)
inleft[2:(nodenumOfTree+1)]<-tr$left
inright<-matrix(0,nodenumOfTree+1,1)
inright[2:(nodenumOfTree+1)]<-tr$right
inparent<-matrix(0,nodenumOfTree+1,1)
inparent[2:(nodenumOfTree+1)]<-parent
inval<-matrix(0,nodenumOfTree+1,1)
inval[2:(nodenumOfTree+1)]<-tr$mean  #tr$val
invec<-matrix(0,nodenumOfTree+1,1)
invec[2:(nodenumOfTree+1)]<-tr$direc

for (i in 1:(nodenumOfTree+1)){
  if (is.na(inval[i])){
       inval[i]<-0
       invec[i]<-0
  }
}

ininfopointer<-matrix(0,nodenumOfTree+1,1)
ininfopointer[2:(nodenumOfTree+1)]<-mi$infopointer

invalue<-matrix(0,atomnum+1,1)
invalue[2:(atomnum+1)]<-mi$value
inlow<-matrix(0,atomnum*d+1,1)
inupp<-matrix(0,atomnum*d+1,1)
for (i in 1:atomnum){
   for (j in 1:d){
       inlow[1+(i-1)*d+j]=mi$low[i,j]
       inupp[1+(i-1)*d+j]=mi$upp[i,j]
   }
}
innodefinder<-matrix(0,atomnum+1,1)
innodefinder[2:(atomnum+1)]<-mi$nodefinder

inlowtr<-matrix(0,nodenumOfTree*d+1,1)
inupptr<-matrix(0,nodenumOfTree*d+1,1)
for (i in 1:nodenumOfTree){
   for (j in 1:d){
       inlowtr[1+(i-1)*d+j]=tr$low[i,j]
       inupptr[1+(i-1)*d+j]=tr$upp[i,j]
   }
}

# we have tree with "nodenumOfTree" nodes
# we hae assocoated info with "atomnum" elements => info for each leaf
#   that is, atomnum = number of leaves

dentree<-.C("proftreeC",
               as.integer(numofall),
               as.integer(atomnum),
               as.double(inlevseq),
               as.integer(d),
               as.integer(Q),
               as.double(instep),
               as.double(insuppo), 
               as.integer(nodenumOfTree),
               as.integer(inleft),
               as.integer(inright),
               as.integer(inparent),
               as.integer(inval),
               as.integer(invec),
               as.integer(ininfopointer),
               as.integer(inlowtr),
               as.integer(inupptr),
               as.double(invalue),
               as.integer(inlow),
               as.integer(inupp),
               as.integer(innodefinder),
               level = double(numofall+1),
               parent = integer(numofall+1),
               volume = double(numofall+1),
               center = double(d*numofall+1),
               efek = integer(1),
PACKAGE="denpro")
               #component = integer(numofall+1),
               #AtomlistAtomOut = integer(numofall+1),
               #AtomlistNextOut = integer(numofall+1),
               #numOfAtoms = integer(1))

#               lapu = double(numofall+1))

efek<-dentree$efek
numOfAtoms<-dentree$numOfAtoms
invalue<-dentree$level[2:(efek+1)]
parent<-dentree$parent[2:(efek+1)]
#component<-dentree$component[2:(efek+1)]
#AtomlistAtom<-dentree$AtomlistAtom[2:(numOfAtoms+1)]
#AtomlistNext<-dentree$AtomlistNext[2:(numOfAtoms+1)]

#if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL 
if (ccen==TRUE) cvol<-TRUE
if (cvol){
#  volume<-cvolumbag(component=component,AtomlistAtom=AtomlistAtom,AtomlistNext=AtomlistNext,low=tr$low,upp=tr$upp,steppi=tr$step)
   volume<-dentree$volume[2:(efek+1)]
#  kerroin<-cinte(invalue,volume,parent) 
#  sepvalnor<-invalue/kerroin
   sepvalnor<-invalue
} 
else{  
  volume<-NULL
  sepvalnor<-NULL
}

if (ccen && cvol){
  #center<-ccentebag(component,AtomlistAtom,AtomlistNext,tr$low,tr$upp,volume,
  #                  tr$step,tr$suppo)
  outcenter<-dentree$center[2:(d*efek+1)]
  center<-matrix(0,efek,d)
  for (i in 1:efek){
     for (j in 1:d){
        center[i,j]<-outcenter[(i-1)*d+j]
     }
  }
  }
  else{
      center<-NULL
  }


} #else (tr$left[1]>0)


return(list(parent=parent,level=sepvalnor,invalue=invalue,
volume=volume,center=t(center)))    #nodefrek=nodefrek))
#values: normeeratut arvot
#invalues: alkuperaiset frekvenssit/arvot 
#nodefrek: kunkin solmun frekvenssi
}










proftreeR<-function(tr,
Q=NULL,frekv=NULL,cvol=TRUE,ccen=TRUE,cfre=FALSE)
{
#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(4,4)
#Q<-3

d<-dim(tr$upp)[2]

nodenumOfTree<-length(tr$left)

low<-tr$low
upp<-tr$upp
val<-tr$val
#low<-matrix(0,nodenumOfTree,d)
#upp<-matrix(0,nodenumOfTree,d)
#val<-matrix(NA,nodenumOfTree,1)
#for (i in 1:nodenumOfTree){
#  dimu<-tr$vec[i]
#  if (!is.na(dimu) && (dimu>0)) 
#       val[i]<-tr$suppo[2*dimu-1]+tr$val[i]*tr$step[dimu]
#  for (j in 1:d){
#      low[i,j]<-tr$suppo[2*j-1]+tr$low[i,j]*tr$step[j]
#      upp[i,j]<-tr$suppo[2*j-1]+tr$upp[i,j]*tr$step[j]
#  }
#}

# make parent
parent<-matrix(0,length(tr$left),1)
node<-1
while (node<=length(tr$left)){
   if ((!is.na(tr$left[node])) && (tr$left[node]!=0)){
        parent[tr$left[node]]<-node
   }
   if ((!is.na(tr$right[node])) && (tr$left[node]!=0)){
        parent[tr$right[node]]<-node
   }
   node<-node+1
}

mi<-makeinfo(tr$left,tr$right,tr$mean,low,upp)
infopointer<-mi$infopointer
terminalnum<-mi$terminalnum
low<-mi$low
upp<-mi$upp
nodefinder<-mi$nodefinder
value<-mi$value

{
if (!is.null(Q)){
   maxval<-max(value)
   minval<-min(value)
   step<-(maxval-minval)/Q
   levseq<-seq(from=minval,to=maxval-step,by=step)
}
else{
   eppsi<-0        #0.0000001
   levseq<-matrix(0,length(value),1)
   ordu<-order(value)
   ru<-1
   #car<-ordu[ru]
   #while ((ru <= length(value)) && (value[car]==0)){
   #     ru<-ru+1
   #     car<-ordu[ru]
   #}  # we have found first non zero
   laskuri<-1
   car<-ordu[ru]
   levseq[laskuri]<-value[car]-eppsi
   while (ru < length(value)){
       carnew<-ordu[ru+1]
       if (value[carnew]>value[car]){
          laskuri<-laskuri+1
          levseq[laskuri]<-value[carnew]-eppsi
      }
      ru<-ru+1
   }
   levseq<-levseq[1:laskuri]
   Q<-laskuri
}
}

levfrekv<-matrix(0,Q,1)
atomnum<-length(value)
for (i in 1:atomnum){
   for (j in 1:Q){
      if (value[i]>=levseq[j]){
         levfrekv[j]<-levfrekv[j]+1
      }
   }
}
numofall<-sum(levfrekv)

dentree<-decombag(numofall,levseq,
tr$left,tr$right,val,tr$vec,infopointer,parent,
nodenumOfTree,terminalnum,
value,low,upp,nodefinder,
d)

invalue<-dentree$level
parent<-dentree$parent
component<-dentree$component
AtomlistAtom<-dentree$AtomlistAtom
AtomlistNext<-dentree$AtomlistNext

#if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL 

if (ccen==TRUE) cvol<-TRUE
if (cvol){
  volume<-cvolumbag(component,AtomlistAtom,AtomlistNext,tr$low,tr$upp,
                    steppi=tr$step)
  kerroin<-cinte(invalue,volume,parent) 
  sepvalnor<-invalue/kerroin
} 
else{  
  volume<-NULL
  sepvalnor<-NULL
}

if (ccen && cvol){
  center<-ccentebag(component,AtomlistAtom,AtomlistNext,tr$low,tr$upp,volume,
                    tr$step,tr$suppo)
  }
  else{
      center<-NULL
  }

return(list(parent=parent,level=sepvalnor,invalue=invalue,
volume=volume,center=center))    #nodefrek=nodefrek))
#values: normeeratut arvot
#invalues: alkuperaiset frekvenssit/arvot 
#nodefrek: kunkin solmun frekvenssi
}









prunemodes<-function(lst,modenum=1,num=NULL,exmalim=NULL,maxnum=NULL)
{
# prunes from a level set tree "lst" the modes with "num" 
# smallest excess masses 
# or the modes with smaller excess mass than "exmalim"

if (is.null(num)){
    curmodenum<-moodilkm(lst$parent)$lkm
    num<-curmodenum-modenum
}

go.on<-TRUE
nn<-1
while (go.on){

  len<-length(lst$parent)
  child.frekve<-matrix(0,len,1)
  for (i in 1:len){
     if (lst$parent[i]>0) 
     child.frekve[lst$parent[i]]<-child.frekve[lst$parent[i]]+1
  }

  ml<-moodilkm(lst$parent)
  mode.list<-ml$modloc
  roots.of.modes<-matrix(0,length(mode.list),1)
  for (aa in 1:length(mode.list)){
      node<-mode.list[aa]
      while ((lst$parent[node]>0) && (child.frekve[lst$parent[node]]==1)){ 
          node<-lst$parent[node]
      }
      roots.of.modes[aa]<-node
  }

  em<-excmas(lst)
  or<-order(em[roots.of.modes])
  smallest<-ml$modloc[or[1]]
  if (nn==1) exma.of.modes<-em[roots.of.modes]

  node<-smallest
  emsmallest<-em[node]

  if ((is.null(exmalim)) || ((!is.null(exmalim)) && (emsmallest<=exmalim))){

     rem.list<-c(node)
     while ((lst$parent[node]>0) && (child.frekve[lst$parent[node]]==1)){ 
           node<-lst$parent[node]
           rem.list<-c(rem.list,node)
     }

     for (kk in 1:length(rem.list)){
        remo<-rem.list[kk]
        for (ll in 1:length(lst$parent)){
            if (lst$parent[ll]>remo) lst$parent[ll]<-lst$parent[ll]-1
        }
        lst$parent<-lst$parent[-remo]
     }
     lst$level<-lst$level[-rem.list]
     lst$volume<-lst$volume[-rem.list]
     lst$center<-lst$center[,-rem.list]
     lst$distcenter<-lst$distcenter[,-rem.list]
     lst$proba<-lst$proba[-rem.list]
     lst$infopointer<-lst$infopointer[-rem.list]
  }
  else if ((!is.null(exmalim)) && (emsmallest>exmalim)) go.on<-FALSE

  nn<-nn+1
  if ((nn>num) && (is.null(exmalim))) go.on<-FALSE
  if ((!is.null(maxnum)) && (nn>maxnum)) go.on<-FALSE 
}

lst$exma.of.modes<-exma.of.modes

return(lst=lst)
}


pruneprof<-function(mt){
#prunes profile so that only root and nodes with siblings are left
#
#mt is a result from multitree
#
roots<-mt$roots
child<-mt$child
sibling<-mt$sibling
siborder<-mt$siborder
#
itemnum<-length(child)
newchild<-matrix(0,itemnum,1)
#
rootnum<-length(roots)
#
for (i in 1:rootnum){
    pino<-matrix(0,itemnum,1)
    pino[1]<-roots[i]  
    pinin<-1
    while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        if (sibling[cur]>0){
              pinin<-pinin+1
              pino[pinin]<-sibling[cur]
        }
        while (child[cur]>0){    #go to left and put right nodes to stack
             candi<-child[cur]
             while ((child[candi]>0) && (sibling[candi]==0)){
                 candi<-child[candi]
             }
             if (sibling[candi]>0){  #if candi has siblings
                newchild[cur]<-candi
                pinin<-pinin+1
                pino[pinin]<-sibling[candi]
             } 
             cur<-candi
        }
    }
}
return(list(roots=roots,child=newchild,sibling=sibling,siborder=siborder))
}







qq.plot<-function(dendat=NULL,compa="gauss",basis="gauss",
mean=0,sig=1,df=1,
gnum=1000,d=1,R=3,qqtype="1d",cex.lab=1,cex.axis=1,col="blue",lwd=1,flip=FALSE,
xlab="compared quantiles",ylab="empirical quantiles")
{
if (qqtype=="1d"){
   n<-length(dendat) #dim(dendat)[1]
   p<-(seq(1:n)-1/2)/n
   if (compa=="gauss") x<-qnorm(p,mean=mean,sd=sig)
   if (compa=="student") x<-sig*qt(p,df=df)+mean
   if (compa=="unif") x<-sig*qunif(p)+mean
   if (compa=="exp") x<-sig*qexp(p)+mean
   if (compa=="doubleexp"){
       x<-sig*qexp(p)+mean
       alku<-which(p<0.5)
       loppu<-which(p>=0.5)
       x[alku]<--sig*qexp(1-2*p[alku])+mean
       x[loppu]<-sig*qexp(2*p[loppu]-1)+mean
   }
   y<-dendat[order(dendat)]
   tyyppi<-"p"
}
if (qqtype=="lower"){
   n<-length(dendat) #dim(dendat)[1]
   p<-(seq(1:n)-1/2)/n
   if (compa=="gauss") x<-qnorm(p/2,mean=mean,sd=sig)
   if (compa=="student") x<-sig*qt(p/2,df=df)+mean
   if (compa=="unif") x<-sig*qunif(p/2)+mean
   if (compa=="exp") x<-sig*qexp(p/2)+mean
   y<-dendat[order(dendat)]
   tyyppi<-"p"
}
if (qqtype=="p2v"){
     rp<-tailfunc(R,d,type=compa,gnum=gnum,sig=sig,nu=df)
     x<-rp$volu
     rp2<-tailfunc(R,d,type=basis,gnum=gnum,sig=sig,nu=df)
     y<-rp2$volu
     tyyppi="l"
     ylab<-"empirical"
     xlab<-"model"
}

if (!flip){
plot(x,y,type=tyyppi,ylab=ylab,xlab=xlab,cex.lab=cex.lab,cex.axis=cex.axis)
maxxy<-max(max(x),max(y))
minxy<-min(min(x),min(y))
segments(minxy,minxy,maxxy,maxxy,col=col,lwd=lwd)
}
if (flip){
 plot(y,x,type=tyyppi,ylab=xlab,xlab=ylab,cex.lab=cex.lab,cex.axis=cex.axis)
 maxxy<-max(max(x),max(y))
 minxy<-min(min(x),min(y))
 segments(minxy,minxy,maxxy,maxxy,col=col,lwd=lwd)
}

}





quanti<-function(values,lkm,base){
#Quantises a vecor of values
#
#values is len-vector
#lkm is positive integer
#base>0
#
#returns len-vector
#
ma<-max(values)
askel<-ma/(lkm-1)
len<-length(values)
ans<-matrix(0,len,1)
for (i in 1:len){
  inv<-base^(values[i]*log(ma+1,base)/ma)-1
  ind<-round(inv/askel)+1
  diskr<-ma*seq(0,lkm-1,1)/(lkm-1)
  disinv<-diskr[ind]
  ans[i]<-ma*log(disinv+1,base)/log(ma+1,base)
}
return(ans)
}
rota.seq<-function(dendat,col,pcf,ste=3,cut=dim(dendat)[1],simple=TRUE)
{
i1<-1
i2<-2
i3<-3
aa<-seq(0,2*pi,ste)
bb<-seq(0,2*pi,ste)
cc<-seq(0,2*pi,ste)
ii<-1
while (ii<=length(aa)){
  alpha<-aa[ii]
  jj<-1
  while (jj<=length(bb)){
    beta<-bb[jj]
    kk<-1
    while (kk<=length(cc)){
      gamma<-cc[kk]
      dexdat<-rotation3d(dendat,alpha,beta,gamma)
      plot.histdata(dexdat,col,pcf,i1,i2,i3,simple=simple,cut=cut,
      xlab=paste(as.character(ii),as.character(jj),as.character(kk)))
      kk<-kk+1
    }
  jj<-jj+1
  }
ii<-ii+1
}


}
rotation2d<-function(dendat,alpha){
  Rx<-matrix(0,2,2)
  Rx[1,]<-c(cos(alpha),-sin(alpha))
  Rx[2,]<-c(sin(alpha),cos(alpha)) 
  detdat<-Rx%*%t(dendat)
  detdat<-t(detdat)
  return(detdat)
}

rotation3d<-function(dendat,alpha,beta,gamma){
  Rx<-matrix(0,3,3)
  Rx[1,]<-c(1,0,0)
  Rx[2,]<-c(0,cos(alpha),-sin(alpha))
  Rx[3,]<-c(0,sin(alpha),cos(alpha))
  Ry<-matrix(0,3,3)
  Ry[1,]<-c(cos(beta),0,sin(beta))
  Ry[2,]<-c(0,1,0)
  Ry[3,]<-c(-sin(beta),0,cos(beta))
  Rz<-matrix(0,3,3)
  Rz[1,]<-c(cos(gamma),-sin(gamma),0)
  Rz[2,]<-c(sin(gamma),cos(gamma),0)
  Rz[3,]<-c(0,0,1)
 
  detdat<-Rx%*%Ry%*%Rz%*%t(dendat)
  detdat<-t(detdat)
  return(detdat)
}






rotation<-function(t,d=2,basis=FALSE)
{

if (d==2){

  rota<-matrix(0,2,2)
  rota[1,1]<-cos(t)
  rota[1,2]<-sin(t)
  rota[2,1]<--sin(t)
  rota[2,2]<-cos(t)

}

if ((d==2) && (basis)){

  rota<-matrix(0,2,2)
  basis1<-c(1,1)
  basis2<-c(-1,1)
  basis1<-basis1/sqrt(sum(basis1^2))
  basis2<-basis2/sqrt(sum(basis2^2))
  rota[,1]<-basis1
  rota[,2]<-basis2

}

if (d==4){

  rotxy<-matrix(0,4,4)
  for (i in 1:4) rotxy[i,i]<-1
  rotxy[1,1]<-cos(t)
  rotxy[1,2]<-sin(t)
  rotxy[2,1]<--sin(t)
  rotxy[2,2]<-cos(t)

  rotyz<-matrix(0,4,4)
  for (i in 1:4) rotyz[i,i]<-1
  rotyz[2,2]<-cos(t)
  rotyz[2,3]<-sin(t) 
  rotyz[3,2]<--sin(t)
  rotyz[3,3]<-cos(t)

  rotzx<-matrix(0,4,4)
  for (i in 1:4) rotzx[i,i]<-1
  rotzx[1,1]<-cos(t)
  rotzx[1,3]<--sin(t)
  rotzx[3,1]<-sin(t)
  rotzx[3,3]<-cos(t)

  rotxw<-matrix(0,4,4)
  for (i in 1:4) rotxw[i,i]<-1
  rotxw[1,1]<-cos(t)
  rotxw[1,4]<-sin(t)
  rotxw[4,1]<--sin(t)
  rotxw[4,4]<-cos(t)

  rotyw<-matrix(0,4,4)
  for (i in 1:4) rotyw[i,i]<-1
  rotyw[2,2]<-cos(t)
  rotyw[2,4]<--sin(t)
  rotyw[4,2]<-sin(t)
  rotyw[4,4]<-cos(t)

  rotzw<-matrix(0,4,4)
  for (i in 1:4) rotzw[i,i]<-1
  rotzw[3,3]<-cos(t)
  rotzw[2,4]<--sin(t)
  rotzw[4,3]<-sin(t)
  rotzw[4,4]<-cos(t)

  rota<-rotxy%*%rotyz%*%rotzx%*%rotxw%*%rotyw%*%rotzw

}

return(rota)
}

scaletable<-function(estiseq,paletti=NULL,shift=0,ptext=0,ptextst=0,
bm=NULL,#mt=NULL,
levnum=60,levnumst=60,redu=TRUE,
volu.modelabel=TRUE,volu.colo=TRUE,st.modelabel=FALSE,st.colo=TRUE
)
{
# preparation
if ((length(estiseq$hseq)>1) && (estiseq$hseq[1]<estiseq$hseq[2])){  
    hnum<-length(estiseq$hseq)
    estiseq$hseq<-estiseq$hseq[seq(hnum,1)]
    apuseq<-list(estiseq$lstseq[[hnum]])
    i<-2
    while (i <= hnum){
         apuseq<-c(apuseq,list(estiseq$lstseq[[hnum-i+1]]))
         i<-i+1 
   }
   estiseq$lstseq<-apuseq
}

if (estiseq$type=="carthisto")  smootseq<--estiseq$leaf
else if (estiseq$type=="greedy") smootseq<--estiseq$hseq
else if (estiseq$type=="bagghisto") smootseq<--estiseq$hseq
else smootseq<-estiseq$hseq
hnum<-length(smootseq)
d<-dim(estiseq$lstseq[[hnum]]$center)[1]

if (estiseq$type=="carthisto") redu<-FALSE
if (is.null(estiseq$stseq)) levnumst<-NULL

if (is.null(paletti))
paletti<-c("red","blue","green","turquoise","orange","navy",
"darkgreen","orchid",colors()[50:100])

# prune the level set trees
if ((!is.null(levnum)) && (redu)){
   for (i in 1:hnum){
      lf<-treedisc(estiseq$lstseq[[i]],estiseq$pcfseq[[i]],levnum) 
      if (i==1){
           if (hnum==1){
               reduseq<-lf
           }
           else{
               reduseq<-list(lf)
           }
      }
      else{
          reduseq<-c(reduseq,list(lf))
      }
  }
  estiseq$lstseq<-reduseq
}

# prune the shape trees
if ((!is.null(levnumst)) && (redu)){
   for (i in 1:hnum){
      lf<-treedisc(estiseq$stseq[[i]],estiseq$pcfseq[[i]],levnumst) 
      if (i==1){
           if (hnum==1){
               reduseq<-lf
           }
           else{
               reduseq<-list(lf)
           }
      }
      else{
          reduseq<-c(reduseq,list(lf))
      }
  }
}
else reduseq<-estiseq$stseq

#if (is.null(mt)) 
mt<-modegraph(estiseq)
if (is.null(bm)) bm<-branchmap(estiseq)

####################################
devicontrol<-2
devibranch<-3
devibary<-4
devimodet<-5
devivolu<-6
deviradi<-7
deviloca<-8

xmin<--0.5
xmax<-0.5
ymin<--1
ymax<-1
lkm<-7
step<-(ymax-ymin)/(lkm-1)
heig<-seq(ymin,ymax,step)
xloc<-0

# control window
dev.new(width=2,height=6)
plot(x="",y="",xlab="",ylab="",xaxt="n",yaxt="n",
xlim=c(xmin,xmax),ylim=c(ymin,ymax))
text(xloc,heig[lkm],"I")      #"Mode graph")
text(xloc,heig[lkm-1],"II")   #"Map of branches")
text(xloc,heig[lkm-2],"III")  #"Volume plot")
text(xloc,heig[lkm-3],"IV")   #"Barycenter plot")
text(xloc,heig[lkm-4],"V")    #"Radius plot")
text(xloc,heig[lkm-5],"VI")   #"Location plot")
text(xloc,heig[lkm-6],"STOP")

devit<-matrix(0,lkm,1)
devit[1]<-devimodet
devit[2]<-devibranch
devit[3]<-devivolu
devit[4]<-devibary
devit[5]<-deviradi
devit[6]<-deviloca
devit[7]<-devicontrol

# choose estimate
indeksi<-1
pr<-estiseq$lstseq[[indeksi]]
pcf<-estiseq$pcfseq[[indeksi]]
if (!is.null(levnumst)){ 
      st<-estiseq$stseq[[indeksi]]
      stredu<-reduseq[[indeksi]]
}
hcur<-estiseq$hseq[indeksi]

# branch map
dev.new(width=4,height=5)
phi<-40
theta<-10
persp(x=bm$level,y=bm$h,z=bm$z, xlab="level",ylab="h",zlab="",
ticktype="detailed",col=bm$col,phi=phi,theta=theta)
title(main="II Map of branches")

# barycenter plot
dev.new(width=3.5,height=4)
coordi<-1
icolo<-mt$colot[mt$low[1]:mt$upp[1]]
inodes<-mt$nodepointer[mt$low[1]:mt$upp[1]]
modlab<-plotbary(pr,coordi=coordi,ptext=ptext,
        modlabret=TRUE,modecolo=icolo,modepointer=inodes)
title(main="IV Barycenter plot",
      sub=paste("coordinate",as.character(coordi)))

# mode tree
dev.new(width=4,height=5)
coordi<-1
plotmodet(mt,coordi=coordi)
modelocx<-modlab$modelocat[,coordi]+shift
modelocy<-smootseq[1]
labels<-modlab$labels
text(modelocx,modelocy,labels)
title(main="I Mode graph",sub=paste("coordinate",as.character(coordi)))

# volume plot
dev.new(width=3.5,height=4)
icolo<-mt$colot[mt$low[1]:mt$upp[1]]
inodes<-mt$nodepointer[mt$low[1]:mt$upp[1]]
plotvolu(pr,ptext=ptext,modelabel=volu.modelabel,colo=volu.colo,
         modecolo=icolo,modepointer=inodes)
title(main="III Volume plot",
      sub=paste("h=",as.character(round(hcur,digits=3))))

# radius plot
if (!is.null(levnumst)){
  refelab<-"moodi"
  st.moodi<-st

  lev<-0.1*max(pcf$value)  
  refe<-st$bary
  st.bary<-leafsfirst(pcf,lev=lev,refe=refe)

  dev.new(width=3,height=4)
  plotvolu(stredu,ptext=ptextst,symbo="T",
           modelabel=st.modelabel,colo=st.colo)
  title(main="V Radius plot",
        sub=paste("level=",as.character(round(lev,digits=3)),
        ", ref.point=mode"))
}

# location plot
if (!is.null(levnumst)){
  dev.new(width=3,height=4)
  lcoordi<-1
  plotbary(stredu,coordi=lcoordi,ptext=ptextst,symbo="T")
  title(main="VI Location plot",
        sub=paste("coordinate",as.character(lcoordi)))
}


##################################################################
ylow<-matrix(0,lkm,1)
yupp<-matrix(0,lkm,1)
for (i in 1:lkm){
    ylow[i]<-heig[i]-step/2
    yupp[i]<-heig[i]+step/2
}
ylow<-ylow[lkm:1]
yupp<-yupp[lkm:1]

dev.set(which = devicontrol)
loc<-locator(1)
while (loc$y>=yupp[lkm]){
  for (i in 1:lkm){
      if ((loc$y>ylow[i]) && (loc$y<=yupp[i])){
                 devi<-devit[i]
      }
  }
  dev.set(which = devi)
  loc<-locator(1)

  # interaction in modegraph
  if (devi==devimodet){
       alaraja<-smootseq[length(smootseq)]
       while (loc$y>=alaraja){
          coordi<-1
          ylamodet<-smootseq[1]
          while (loc$y>=ylamodet){
              if (coordi<=(d-1)) coordi<-coordi+1 else coordi<-1
              plotmodet(mt,coordi=coordi)
              modelocx<-modlab$modelocat[,coordi]+shift
              modelocy<-smootseq[indeksi]
              labels<-modlab$labels
              text(modelocx,modelocy,labels)
              title(main="I Mode graph",
                    sub=paste("coordinate",as.character(coordi)))

              loc<-locator(1)
          }
          
          if (loc$y>=alaraja){
             alamidi<-(smootseq[1]+smootseq[1+1])/2
             if (loc$y>=alamidi) indeksi<-1
             for (i in 2:(hnum-1)){
                alamidi<-(smootseq[i]+smootseq[i+1])/2
                ylamidi<-(smootseq[i-1]+smootseq[i])/2

                if ((loc$y>=alamidi) && (loc$y<ylamidi)) indeksi<-i
             }
             ylamidi<-(smootseq[hnum-1]+smootseq[hnum])/2
             if (loc$y<ylamidi) indeksi<-hnum

             pr<-estiseq$lstseq[[indeksi]]
             pcf<-estiseq$pcfseq[[indeksi]]
             hcur<-estiseq$hseq[[indeksi]]
             if (!is.null(levnumst)){ 
                    lev<-0.1*max(pcf$value)
                    st<-estiseq$stseq[[indeksi]]
                    st.moodi<-st
                    st.bary<-NULL
                    stredu<-reduseq[[indeksi]]
             }

             dev.set(which = devivolu)
             icolo<-mt$colot[mt$low[indeksi]:mt$upp[indeksi]]
             inodes<-mt$nodepointer[mt$low[indeksi]:mt$upp[indeksi]] 
             plotvolu(pr,ptext=ptext,modelabel=volu.modelabel,colo=volu.colo,
                      modecolo=icolo,modepointer=inodes)
             title(main="III Volume plot",
                   sub=paste("h=",as.character(round(hcur,digits=3))))
  
             dev.set(which = devibary) 
             coordi<-1
             modlab<-plotbary(pr,coordi=coordi,ptext=ptext,
                              modlabret=T,modecolo=icolo,modepointer=inodes)
             title(main="IV Barycenter plot",
                   sub=paste("coordinate",as.character(coordi)))

             dev.set(which = deviradi) 
             plotvolu(stredu,ptext=ptextst,symbo="T",
                      modelabel=st.modelabel,colo=st.colo)
             title(main="V Radius plot",
                   sub=paste("level=",as.character(round(lev,digits=3)),
                             ", ref.point=mode"))
 
             dev.set(which = deviloca) 
             lcoordi<-1
             plotbary(stredu,coordi=lcoordi,ptext=ptextst,symbo="T")
             title(main="VI Location plot",
                   sub=paste("coordinate",as.character(lcoordi)))
             
             dev.set(which = devimodet)

             modelocx<-modlab$modelocat[,coordi]+shift
             modelocy<-smootseq[indeksi]
             labels<-modlab$labels
             text(modelocx,modelocy,labels)

             loc<-locator(1)
          }
       }
       #dev.set(which = devicontrol)
  }

  # interaction in volume plot
  if (devi==devivolu){
     alaasso<-0
     ylaasso<-max(pr$level)
     alax<-0
     ylax<-pr$volume[1]
     while (loc$y>=alaasso){
        if (loc$x>=0){
           if (loc$y>ylaasso) plotvolu(pr)
           else if (loc$x>0){
              keskip<-alax+(ylax-alax)/2
              if (loc$x >= keskip) ylax<-loc$x
              else                 alax<-loc$x
              icolo<-mt$colot[mt$low[indeksi]:mt$upp[indeksi]]
              inodes<-mt$nodepointer[mt$low[indeksi]:mt$upp[indeksi]] 
              plotvolu(pr,xlim=c(alax,ylax),ptext=ptext,
                       modelabel=volu.modelabel,colo=volu.colo,
                       modecolo=icolo,modepointer=inodes)
           }
           title(main="III Volume plot",
                 sub=paste("h=",as.character(round(hcur,digits=3))))
        }     
        else if (!is.null(levnumst)){
             maksi<-max(pr$level)
             mode<-locofmax(pcf)
             lev<-min(max(0,loc$y),maksi)
             st<-leafsfirst(pcf,refe=mode,lev=lev)
             st.moodi<-st
             st.bary<-NULL
             if (redu) stredu<-treedisc(st,pcf,ngrid=levnumst) else stredu<-st
             refelab<-"moodi"

             dev.set(which = deviradi) 
             plotvolu(stredu,ptext=ptextst,symbo="T",
                      modelabel=st.modelabel,colo=st.colo)
             title(main="V Radius plot",
                   sub=paste("level=",as.character(round(lev,digits=3)),
                             ", ref.point=mode"))

             dev.set(which = deviloca) 
             lcoordi<-1
             plotbary(stredu,coordi=lcoordi,ptext=ptextst,symbo="T")
             title(main="VI Location plot",
                   sub=paste("coordinate",as.character(lcoordi)))

             dev.set(which = devivolu)       
        }
        loc<-locator(1)
     }
  }

  # interaction in barycenter plot
  if (devi==devibary){
      coordi<-1
      icolo<-mt$colot[mt$low[indeksi]:mt$upp[indeksi]]
      inodes<-mt$nodepointer[mt$low[indeksi]:mt$upp[indeksi]]
      modlab<-plotbary(pr,coordi=coordi,ptext=ptext,
                       modlabret=T,modecolo=icolo,modepointer=inodes)
      title(sub=paste("barycenter plot, coordinate",as.character(coordi)))
      alaasso<-0
      while (loc$y>=alaasso){
         if (coordi<=(d-1)) coordi<-coordi+1 else coordi<-1
         plotbary(pr,coordi=coordi,ptext=ptext,modecolo=icolo,
                  modepointer=inodes,modelabel=TRUE)
         title(main="IV Barycenter plot",
               sub=paste("coordinate",as.character(coordi)))

         loc<-locator(1)
      }
  }

  # interaction in radius plot
  if (devi==deviradi){
       alaraja<-0
       while (loc$y>=alaraja){

          ylaraja<-max(st$level)
          if (loc$y>=ylaraja){
              if (refelab=="moodi"){ 
                   refelab<-"bary"
                   if (is.null(st.bary)){  
                       refe<-st$bary
                       st.bary<-leafsfirst(pcf,lev=lev,refe=refe)
                   }
                   st<-st.bary
              }
              else{ 
                   refelab<-"moodi"
                   st<-st.moodi
              }
              if (redu) stredu<-treedisc(st,pcf,ngrid=levnumst) else stredu<-st
              plotvolu(stredu,ptext=ptextst,symbo="T",
                       modelabel=st.modelabel,colo=st.colo)
              if (refelab=="moodi") 
                 title(main="V Radius plot",
                       sub=paste("level=",as.character(round(lev,digits=3)),
                       ",  mode=refe'nce point"))
              else 
                 title(main="V Radius plot",
                       sub=paste("level=",as.character(round(lev,digits=3)),
                       ", ref.point= barycenter"))

              dev.set(which = deviloca) 
              lcoordi<-1
              plotbary(stredu,coordi=lcoordi,ptext=ptextst,symbo="T")
              title(main="VI Location plot",
                    sub=paste("coordinate",as.character(lcoordi)))

              dev.set(which = deviradi)      
              loc<-locator(1)
          }
          else{
              sarmilkm<-moodilkm(stredu$parent)$lkm
              streduredu<-stredu
              while ((loc$y<ylaraja) && (loc$y>alaraja)){
                   cursarmilkm<-moodilkm(streduredu$parent)$lkm
                   if (cursarmilkm>=2) newsarmilkm<-cursarmilkm-1 
                   else newsarmilkm<-sarmilkm 
                   streduredu<-prunemodes(stredu,modenum=newsarmilkm)
                   plotvolu(streduredu,ptext=ptextst,symbo="T",
                            modelabel=st.modelabel,colo=st.colo)
                   if (refelab=="moodi") 
                         title(main="V Radius plot",
                         sub=paste("level=",as.character(round(lev,digits=3)),
                         ",  mode=refe'nce point"))
                   else 
                       title(main="V Radius plot",
                       sub=paste("level=",as.character(round(lev,digits=3)),
                       ", ref.point=barycenter"))

                   dev.set(which = deviloca) 
                   lcoordi<-1
                   plotbary(streduredu,coordi=lcoordi,ptext=ptextst,symbo="T")
                   title(main="VI Location plot",
                         sub=paste("coordinate",as.character(lcoordi)))
 
                   dev.set(which = deviradi)      
                   loc<-locator(1)
              }
          }
          loc<-locator(1)
      }
  }

  # interaction in location plot
  if (devi==deviloca){
      coordi<-1
      plotbary(stredu,coordi=coordi,ptext=ptextst,symbo="T")
      title(main="VI Location plot",
            sub=paste("coordinate",as.character(coordi)))
      alaasso<-0
      while (loc$y>=alaasso){
         if (coordi<=(d-1)) coordi<-coordi+1 else coordi<-1
         plotbary(stredu,coordi=coordi,ptext=ptextst,symbo="T")
         title(main="VI Location plot",
               sub=paste("coordinate",as.character(coordi)))

         loc<-locator(1)
      }
  }

  # interaction in branching map
  if (devi==devibranch){
      alaraja<--0.4
      while (loc$y>=alaraja){

          if (loc$x>=0) theta<-theta+10 else theta<-theta-10
          if (loc$y>=0) phi<-phi+10 else phi<-phi-10

          persp(x=bm$level,y=bm$h,z=bm$z,col=bm$col,
          xlab="level",ylab="h",zlab="",ticktype="detailed",
          phi=phi,theta=theta)
          title(main="II Map of branches")

         loc<-locator(1)
      }
  }

  # end
  dev.set(which = devicontrol)
  loc<-locator(1)
}

if (!is.null(levnumst)) devlkm<-lkm
else devlkm<-lkm-2
for (i in 1:devlkm) dev.off()

}



shape2d<-function(shtseq, gnum=500, type="radius", type2="slice", 
gnum2=1000, ngrid=30, norma=FALSE, xmax=10, modelim=2, exmalim=NULL,
maxnum=NULL)
{
# type "proba"    type2 "boundary"
lkm<-length(shtseq$level)
d<-length(shtseq$pcf$N)

if (type2=="slice"){

  if (type=="radius") x<-shtseq$level else x<-matrix(0,lkm,1)

  td<-shtseq$shtseq[[1]]
  if (type=="proba") td$volume<-td$proba
  td<-treedisc(td,shtseq$pcf,ngrid=ngrid)
  xy<-lst2xy(td,gnum=gnum)
  ylen<-length(xy$x)
  ystep<-1/(ylen-1)
  y<-seq(0,1,ystep)   #matrix(0,xlen,1)
  z<-matrix(0,length(x),length(y))
  delineator<-matrix(0,10*length(x),d)
  delinrun<-1
  delineatorlevel<-matrix(0,10*length(x),1)
  
  delineator.redu<-matrix(0,4*length(x),d)
  dr.redu<-1
  delineatorlevel.redu<-matrix(0,4*length(x),1)

  for (i in 1:lkm){
     td<-shtseq$shtseq[[i]]

     if (type=="proba"){ 
         tdvolume<-td$volume
         td$volume<-td$proba
         indi<-lkm-i+1
         voluu<-max(tdvolume)  #[1]  #root=1
         if (norma) x[indi]<-(voluu/volball(1,d))^(1/d)
         else x[indi]<-voluu 
     }
     else indi<-i
     td<-treedisc(td,shtseq$pcf,ngrid=ngrid)
     if (length(td$parent)==1) ynew<-0
     else{
        xy<-lst2xy(td,gnum=gnum)   #ma<-matchxy(xy$x,xy$y,y)

        ## normalize
        volu<-xy$x[length(xy$x)]-xy$x[1]
        int<-0
        step<-xy$x[2]-xy$x[1]
        for (j in 1:length(xy$x)){
            int<-int+step*xy$y[j]
        }
        if (norma){
            normavolu<-(volu/volball(1,d))^(1/d)
            b<-volu*normavolu/int
        }
        else b<-volu^2/int
        ynew<-b*xy$y
        ## end normalize

        # location
        ml<-moodilkm(td$parent)
        mc<-t(td$center[,ml$modloc])  #modecent(td)
        modenum<-dim(mc)[1]
        delineator[delinrun:(delinrun+modenum-1),]<-mc     
        delineatorlevel[delinrun:(delinrun+modenum-1)]<-x[indi] 
        delinrun<-delinrun+modenum

        if (modenum>modelim){
            prunum<-modenum-modelim
            pru<-prunemodes(td,prunum,exmalim,num=maxnum)
        }
        else pru<-td 
        ml<-moodilkm(pru$parent)
        mc<-t(pru$center[,ml$modloc])  #modecent(td)
        modenum<-dim(mc)[1]
        delineator.redu[dr.redu:(dr.redu+modenum-1),]<-mc     
        delineatorlevel.redu[dr.redu:(dr.redu+modenum-1)]<-x[indi] 
        dr.redu<-dr.redu+modenum
     }
     z[indi,]<-ynew   
  }
  delineator<-delineator[1:(delinrun-1),]
  delineatorlevel<-delineatorlevel[1:(delinrun-1)]

  delineator.redu<-delineator.redu[1:(dr.redu-1),]
  delineatorlevel.redu<-delineatorlevel.redu[1:(dr.redu-1)]
}

else{ #type2=="boundary"

if (is.null(xmax)){
    td<-shtseq$shtseq[[1]]
    if (type=="proba") td$volume<-td$proba
    xmax<-max(td$volume)
}

ymax<-xmax
step<-2*xmax/(gnum-1)
x<-seq(-xmax,xmax,step)
y<-x
z<-matrix(0,length(x),length(y))

for (i in 1:lkm){
  td<-shtseq$shtseq[[i]]
  if (type=="proba") td$volume<-td$proba
  xy<-lst2xy(td,gnum=gnum2,type=type)  

  ## normalize
  volu<-xy$x[length(xy$x)]-xy$x[1]
  int<-0
  step<-xy$x[2]-xy$x[1]
  for (j in 1:length(xy$x)){
       int<-int+step*xy$y[j]
  }
  b<-volu^2/int
  ynew<-b*xy$y
  ## end normalize

  for (j in 1:length(x)){
      for (k in 1:length(y)){
          len<-sqrt(x[j]^2+y[k]^2)
          xn<-x[j]/len
          yn<-y[k]/len
          th2<-atan(xn/yn)
          if (yn<0) th2<-atan(xn/yn)+pi else if (xn<0) th2<-atan(xn/yn)+2*pi
          propo<-th2/(2*pi) 
          dirind<-max(1,round( propo*length(xy$x) ))
          rho<-ynew[dirind]
          if (len<=rho) z[j,k]<-shtseq$level[i]
      }
  }
}

}

return(list(x=x,y=y,z=z,type=type,type2=type2,norma=norma,
            delineator=delineator,delineatorlevel=delineatorlevel,
            delineator.redu=delineator.redu,
            delineatorlevel.redu=delineatorlevel.redu))
}


shapetree<-function(et,lev,bary,ordmet="etaisrec",levmet="proba")
{
# et is an evaluation tree

d<-length(et$step)

# order the atoms for the level set with level "lev"

lenni<-length(et$value)
distat<-matrix(0,lenni,1)
infopointer<-matrix(0,lenni,1)
lkm<-0
for (i in 1:lenni){
  if (et$value[i]>=lev){
     lkm<-lkm+1
     nod<-i  #nod<-et$nodefinder[i]
     if (ordmet=="etaisrec"){
         recci<-matrix(0,2*d,1)
         for (jj in 1:d){
            recci[2*jj-1]<-et$support[2*jj-1]+et$step[jj]*et$low[nod,jj]
            recci[2*jj]<-et$support[2*jj-1]+et$step[jj]*et$upp[nod,jj]
         }
         distat[lkm]<-etaisrec(bary,recci)
     }
     else{
        lowi<-matrix(0,d,1)
        uppi<-matrix(0,d,1)
        for (jj in 1:d){
            lowi[jj]<-et$support[2*jj-1]+et$step[jj]*et$low[nod,jj]
            uppi[jj]<-et$support[2*jj-1]+et$step[jj]*et$upp[nod,jj]
        }
        baryc<-lowi+(uppi-lowi)/2  #et$low[nod,]+(et$upp[nod,]-et$low[nod,])/2  
        distat[lkm]<-etais(baryc,bary)
     }
     infopointer[lkm]<-i
  }
}
distat<-distat[1:lkm]
infopointer<-infopointer[1:lkm]   #pointe->et$value,et$nodefinder

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
center<-matrix(0,lkm,d)
radius<-matrix(0,lkm,1)

proba<-matrix(0,lkm,1)
ekamome<-matrix(0,lkm,d)
highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0
# volume calculation
vol<-1
k<-1
ip<-infopointer[node]  #et$nodefinder[infopointer[node]]
while (k<=d){
    vol<-vol*(et$upp[ip,k]-et$low[ip,k])*et$step[k]
    k<-k+1
}
volume[node]<-vol
ip2<-infopointer[node]
proba[node]<-et$value[ip2]*vol
radius[node]<-distat[ord[node]]
# ekamome calculation
newcente<-matrix(0,d,1)
for (j in 1:d){
  volmin<-1
  k<-1
  while (k<=d){
      if (k!=j){
         volmin<-volmin*(et$upp[ip,k]-et$low[ip,k])*et$step[k]
      }
      k<-k+1
  }
  ala<-et$support[2*j-1]+et$step[j]*et$low[ip,j]
  yla<-et$support[2*j-1]+et$step[j]*et$upp[ip,j]
  newcente[j]<-volmin*(yla^2-ala^2)/2
}
ekamome[node,]<-newcente

beg<-node             #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-et$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-et$low[note,i]   
  boundrec[node,2*i]<-et$upp[note,i]  #et$index[infopointer[node],i]
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)  #luo sigleton
    note<-infopointer[node]  #note<-et$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-et$low[note,i]  
         rec1[2*i]<-et$upp[note,i] 
    }
    boundrec[node,]<-rec1

    # volume calculation
    vol<-1
    k<-1
    ip<-infopointer[node]    #et$nodefinder[infopointer[node]]
    while (k<=d){
          vol<-vol*(et$upp[ip,k]-et$low[ip,k])*et$step[k]
          k<-k+1
    }
    volume[node]<-vol
    ip2<-infopointer[node]
    proba[node]<-et$value[ip2]*vol
    radius[node]<-distat[ord[node]]
    # ekamome calculation
    newcente<-matrix(0,d,1)
    for (jj in 1:d){
         volmin<-1
         k<-1
         while (k<=d){
            if (k!=jj){
                volmin<-volmin*(et$upp[ip,k]-et$low[ip,k])*et$step[k]
            }
            k<-k+1
         }
         ala<-et$support[2*jj-1]+et$step[jj]*et$low[ip,jj]
         yla<-et$support[2*jj-1]+et$step[jj]*et$upp[ip,jj]
         newcente[jj]<-volmin*(yla^2-ala^2)/2
    }
    ekamome[node,]<-newcente

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
        istouch<-touchstep(node,curroot,boundrec,child,sibling,
                           infopointer,et$low,et$upp)
        if (istouch==1){
{
           # paivita parent, child, sibling, volume ekamome
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }
           volume[node]<-volume[node]+volume[curroot]
           proba[node]<-proba[node]+proba[curroot]
           ekamome[node,]<-ekamome[node,]+ekamome[curroot,]
           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]
}
        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter
for (i in 1:lkm){
   for (j in 1:d){
       center[i,j]<-ekamome[i,j]/volume[i]
   }
}

if (levmet=="proba")
level<-taillevel(root,#child,sibling,
parent,volume,proba)
else level<-radius

return(list(
parent=parent,volume=volume,center=t(center),level=level,
root=root,
#child=child,sibling=sibling,  #virhe??
infopointer=infopointer,
proba=proba,radius=radius,
bary=bary,maxdis=distat[ord[length(ord)]]))

}





siborder.new<-function(mt)
{
#mt is multitree

itemnum<-length(mt$child)
sibord<-matrix(0,itemnum,1)

#order first roots

rootnum<-length(mt$roots)
for (i in 1:rootnum) sibord[mt$roots[i]]<-i

# then order the other

for (i in 1:rootnum){
   curroot<-mt$roots[i]
   if (mt$child[curroot]>0){
      pino<-matrix(0,itemnum,1)
      pino[1]<-mt$child[curroot]
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          # if not yet ordered, order siblings
          if (sibord[cur]==0){
              indu<-1
              sibord[cur]<-indu
              runner<-cur
              while (mt$sibling[runner]>0){
                  sibord[mt$sibling[runner]]<-indu
                  indu<-indu+1  
                  runner<-mt$sibling[runner] 
              }
          }
          # put to the stack 
          if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
          }
          # go to left and put right nodes to the stack
          while (mt$child[cur]>0){   
             cur<-mt$child[cur]
             # if not yet ordered, order siblings
             if (sibord[cur]==0){
                 indu<-1
                 sibord[cur]<-indu
                 runner<-cur
                 while (mt$sibling[runner]>0){
                     sibord[mt$sibling[runner]]<-indu
                     indu<-indu+1  
                     runner<-mt$sibling[runner]
                 }
             }
             if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
             }
           }
       }
   }
}                                  
return(sibord)
}















siborder<-function(mt,crit,centers)
{
#mt is multitree

roots<-mt$roots
child<-mt$child
sibling<-mt$sibling

itemnum<-length(child)
sibord<-matrix(0,itemnum,1)

#order first roots

rootnum<-length(roots)
if (rootnum==1){
  sibord[roots[1]]<-1
}
else{
  rootlink<-matrix(0,itemnum,1)
  for (i in 1:(rootnum-1)){
     inde<-roots[i]
     rootlink[inde]<-roots[i+1]
  }
  sibord<-levord(roots[1],rootlink,sibord,centers,crit)
}

# then order the other

for (i in 1:rootnum){
   curroot<-roots[i]
   if (child[curroot]>0){
      pino<-matrix(0,itemnum,1)
      pino[1]<-child[curroot]
      pinin<-1
      while (pinin>0){
          cur<-pino[pinin]      #take from stack
          pinin<-pinin-1
          # if not yet ordered, order siblings
          if (sibord[cur]==0){
              sibord<-levord(cur,sibling,sibord,centers,crit)
          }
          # put to the stack 
          if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
          }
          # go to left and put right nodes to the stack
          while (child[cur]>0){   
             cur<-child[cur]
             # if not yet ordered, order siblings
             if (sibord[cur]==0){
                sibord<-levord(cur,sibling,sibord,centers,crit)
             }
             if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
             }
           }
       }
   }
}                                  
return(sibord)
}















siborToModor<-function(tree){
#
# From ordering in siblings to ordering of modes
# We have the right ordering in profile
#
data<-plotprof(tree,plot=F,data=T,cutlev=NULL,ptext=0,info=NULL,
infolift=0,infopos=0)
vecs<-data$vecs
#
parent<-tree$parent
mlkm<-moodilkm(parent)
modloc<-mlkm$modloc
moodinum<-mlkm$lkm    #length(modloc)
#
xcor<-matrix(0,moodinum,1)
for (i in 1:moodinum){
    loc<-modloc[i]
    xcor[i]<-vecs[loc,1]    
}
modloc<-omaord2(modloc,xcor)       
#
return(modloc)
}
sim.1d2modal<-function(n=NULL,seed=1,N=NULL,distr=FALSE)
{
d<-1
M<-c(0,2,4)
mixnum<-length(M)
sig<-matrix(1,mixnum,d)
sig[1]<-0.3
p<-matrix(1,mixnum,1)
p[2]<-2
p<-p/sum(p)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed,d=1)
   return(dendat)
}

if (!is.null(N)){
    xala<--2
    xyla<-7
    support<-c(xala,xyla)
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p,support=support,distr=distr)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


sim.claw<-function(n=NULL,seed=1,N=NULL)
{
d<-1

M<-c(0,-1,-0.5,0,0.5,1)
sig<-c(1,0.1,0.1,0.1,0.1,0.1)
p<-c(0.5,0.1,0.1,0.1,0.1,0.1)
mixnum<-length(M)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed,d=1)
   return(dendat)
}

if (!is.null(N)){
    support<-c(-3,3)
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p,support=support)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


sim.cross<-function(n=NULL,seed=1,N=NULL,sig1=0.5,sig2=1.5)
{
d<-2
mixnum<-2
M<-matrix(0,mixnum,d)
M[1,]<-c(0,0)      
M[2,]<-c(0,0)      
sig<-matrix(1,mixnum,d)
sig[1,1]<-sig1 
sig[1,2]<-sig2   
sig[2,1]<-sig2   
sig[2,2]<-sig1  
p<-matrix(1,mixnum,1)
p<-p/sum(p)

if (!is.null(n)){
   dendat<-simmix(n,M,sig,p,seed=seed)  
   theta<-pi/4
   rotmat<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
   dendat<-dendat%*%rotmat
   return(dendat)
}

if (!is.null(N)){
    theta<-pi/4
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p,theta=theta)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p,theta=theta))
}

}


sim.data<-function(n=NULL,seed=1,N=NULL,type="mulmod",
M=NULL,sig=NULL,p=NULL,d=NULL,
cova=NULL,marginal=NULL,t=NULL,df=NULL,distr=FALSE, noisedim=1,
sig1=0.5,sig2=1.5,diff=0.1,dist=4)
{
if (type=="mixt") return( simmix(n,M,sig,p,seed,d) )

if (type=="mulmod") return( sim.mulmod(n=n,seed=seed,N=N) )

if (type=="fox") return( sim.fox(n=n,seed=seed,N=N) )

if (type=="tetra3d") return( sim.tetra3d(n=n,seed=seed,N=N) )

if (type=="penta4d") return( sim.penta4d(n=n,seed=seed,N=N,dist=dist) )

if (type=="cross") return( sim.cross(n=n,seed=seed,N=N,sig1=sig1,sig2=sig2) )

if (type=="1d2modal") return( sim.1d2modal(n=n,seed=seed,N=N,distr=distr) )

if (type=="claw") return( sim.claw(n=n,seed=seed,N=N) )

if (type=="fssk") return( sim.fssk(n=n,noisedim=noisedim,seed=seed) )

if (type=="nested") return( sim.nested(n=n,seed=seed,N=N) )

if (type=="peaks") return( sim.peaks(n=n,seed=seed,N=N) )

if (type=="mulmodII") return( sim.mulmodII(n=n,seed=seed,N=N) )

if (type=="gauss"){
   eig<-eigen(cova,symmetric=TRUE)
   sigsqm<-eig$vectors%*%diag(eig$values^{1/2})  
   set.seed(seed)
   symmedata<-matrix(rnorm(2*n),n,2)
   dendat<-t(sigsqm%*%t(symmedata))
   if (!is.null(marginal)){
      dendat[,1]<-pnorm(dendat[,1],sd=sqrt(cova[1,1]))
      dendat[,2]<-pnorm(dendat[,2],sd=sqrt(cova[2,2]))
      if (marginal=="student") dendat<-qt(dendat, df=t)
      if (marginal=="gauss") dendat<-qnorm(dendat)
   }
   return(dendat)
}

if (type=="student"){
   eig<-eigen(cova,symmetric=TRUE)
   sigsqm<-eig$vectors%*%diag(eig$values^{1/2})  
   set.seed(seed)
   symmedata<-matrix(rt(2*n,df=df),n,2)
   dendat<-t(sigsqm%*%t(symmedata))
   if (!is.null(marginal)){
       dendat<-pt(dendat,df=df)
       if (marginal=="gauss") dendat<-qnorm(dendat)
   }
   return(dendat)
}

if (type=="gumbel"){
  link<-function(y,g){ return ( (-log(y))^g ) }
  linkinv<-function(y,g){ return ( exp(-y^(1/g)) ) }
  der1<-function(y,g){ return ( -g*(-log(y))^(g-1)/y ) }
  der1inv<-function(y,g){ return ( y ) }
}

if (type=="diff1d"){
   xala<--0
   xyla<-1
   support<-c(xala,xyla)
   d<-1
   M<-c(0.5-diff,0.5+diff)
   mixnum<-length(M)
   sig<-matrix(sig1,mixnum,d)
   p<-matrix(1,mixnum,1)
   p<-p/sum(p)
   pcf<-pcf.func("mixt",N=N,sig=sig,M=M,p=p,support=support,distr=distr)
   return(pcf)
}

}

sim.fox<-function(n=NULL,seed=1,N=NULL)
{
d<-2
mixnum<-14
D<-1.8
M<-matrix(0,mixnum,d)
M[1,]<-c(0,0)      #c(0,0)

M[2,]<-c(D,0)      #c(D1,0)
M[3,]<-c(2*D,0)

M[4,]<-c(0,D)
M[5,]<-c(0,2*D)
M[6,]<-c(0,3*D)

M[7,]<-c(0,-D)
M[8,]<-c(0,-2*D)
M[9,]<-c(0,-3*D)

M[10,]<-c(1.5,3.9*D)
M[11,]<-c(-1.5,3.7*D)
M[12,]<-c(-1.5,4.2*D)
M[13,]<-c(-1.5,4.5*D)
M[14,]<-c(-1.5,4.7*D)

sig<-matrix(1,mixnum,d)
sig[10,1]<-0.7
sig[11,1]<-0.7
sig[12,1]<-0.7
sig[13,1]<-0.7
sig[14,1]<-0.7
p<-matrix(1,mixnum,1)
p[6]<-0.6
p[10]<-0.3
p[11]<-0.25
p[12]<-0.1
p[13]<-0.05
p[14]<-0.05
p<-p/sum(p)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed)
   return(dendat)
}

if (!is.null(N)){
    #eg<-evalgrid(M,sig,p,N)
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


sim.fssk<-function(n,noisedim,seed)
{
# makes n*d data matrix, d=2+noisedim
# 3 moodia, (c,0), (-c,3), (-c,-3)

d<-2+noisedim
hajo<-1
noisehajo<-sqrt(7)
c<-3^(3/2)/2
set.seed(seed)
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






























simmix1d<-function(n,M,sig,p,seed){
#Simulates a mixture of l normal distributions in R^1,
#
#n is the sample size
#M is l-vector, rows are the means
#sig is l-vector, for l:th mixture variance
#p is l-vector, proportion for each mixture
#
#returns n*d-matrix
#
set.seed(seed) 
l<-length(M)
d<-1
data<-rnorm(n)        #n-vektori valkoista kohinaa 
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
data<-t(t(data))  #we make a n*1 matrix
return(data)
}
simmix<-function(n,M,sig,p,seed,d=NULL)
{
#Simulates a mixture of l normal distributions in R^d, l>1
#with diagonal cov matrices

#n is the sample size
#M is l*d-matrix, rows are the means
#sig is l*d-matrix, for l:th mixture d covariances
#p is l-vector, proportion for each mixture

#returns n*d-matrix

if (is.null(d)) d<-dim(M)[2] 

set.seed(seed) 
#if (dim(M)[2]==1) d<-1 else d<-length(M[1,]) 
if (d==1){
  data<-simmix1d(n,M,sig,p,seed)
  }
else{
l<-length(M[,1])
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
return(data)
}
sim.mulmodII<-function(n=NULL,seed=1,N=NULL)
{
d<-2
mixnum<-3
D<-4
M<-matrix(0,mixnum,d)
M[1,]<-c(0,0)   
M[2,]<-c(D,0)  
M[3,]<-c(D/2,D*sqrt(3)/2)   
sig<-matrix(1,mixnum,d)
p<-c(.2,.35,.45)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed)
   return(dendat)
}

if (!is.null(N)){
    #eg<-evalgrid(M,sig,p,N)
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


sim.mulmod<-function(n=NULL,seed=1,N=NULL)
{
d<-2
mnum<-3
D<-3.5  #3 
shift<-0.2  
M<-matrix(0,mnum,d)
M[1,]<-c(0,0)
M[2,]<-c(D,shift)   #c(D,0)
M[3,]<-c(D/2,D)    #   #c(D/2,D*sqrt(3)/2)
sig<-matrix(1,mnum,d)
p<-c(.25,.35,.45)
p<-p/sum(p)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed)
   return(dendat)
}

if (!is.null(N)){
    #eg<-evalgrid(M,sig,p,N)
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


sim.nested<-function(n=NULL,seed=1,N=NULL)
{
d<-2
con<-3^(3/2)/2
shift<-0.8
std<-0.3
mnum<-5
M<-matrix(0,mnum,d)
M[1,]<-c(0,con)
M[2,]<-c(3,-con)
M[3,]<-c(-3,-con) 
M[4,]<-c(shift,con)
M[5,]<-c(-shift,con)

sig<-matrix(1,mnum,d)
sig[4,]<-c(std,std)
sig[5,]<-c(std,std)
p<-c(5/7/3,5/7/3,5/7/3,1/7,1/7)
p<-p/sum(p)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed)
   return(dendat)
}

if (!is.null(N)){
    #eg<-evalgrid(M,sig,p,N)
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


sim.peaks<-function(n=NULL,seed=1,N=NULL)
{
d<-2
mnum<-5
D<-3.5  
shift<-0.6
M<-matrix(0,mnum,d)
M[1,]<-c(0,0)
M[2,]<-c(D,0)
M[3,]<-c(D/2,D)       #c(D/2,D*sqrt(3)/2)
M[4,]<-c(shift+D/2,D)
M[5,]<-c(-shift+D/2,D)

sig<-matrix(1,mnum,d)
std<-0.3
sig[4,]<-c(std,std)
sig[5,]<-c(std,std)

#p<-c(.25,.35,.45)
p<-c(6/8/3,6/8/3,6/8/3,1/8,1/8)
p<-p/sum(p)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed)
   return(dendat)
}

if (!is.null(N)){
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


sim.penta4d<-function(n=NULL,seed=1,N=NULL,dist=4)
{
d<-4
moodi<-5
M<-matrix(0,moodi,d)
#dist<-4     # determine the distance between vertices of the pentahedron
M[1,]<-dist*c(1/2, 0,0,0)
M[2,]<-dist*c(-1/2,0,0,0)
M[3,]<-dist*c(0,sqrt(3)/2,0,0)
M[4,]<-dist*c(0,1/(2*sqrt(3)),sqrt(2/3),0)
M[5,]<-dist*c(0,1/(2*sqrt(3)),1/(2*sqrt(6)),sqrt(15/24))
sig<-matrix(1,moodi,d)
p0<-1/moodi
p<-p0*rep(1,moodi)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed)
   return(dendat)
}

if (!is.null(N)){
    #eg<-evalgrid(M,sig,p,N)
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


sim.tetra3d<-function(n=NULL,seed=1,N=NULL)
{
dist<-3    # determine the distance between vertices of the tetrahedron

d<-3
moodi<-4
M<-matrix(0,moodi,d)
height<-sqrt(3)/2           # sqrt(3)/2 = 0.8660254
len<-1/(2*sqrt(3))          # 1/(2*sqrt(3)) = 0.2886751
kor<-sqrt(2/3)              # sqrt(2/3) = 0.8164966
M[1,]<-dist*c(1/2,0,0)      # ( 1.5, 0.0, 0.0)
M[2,]<-dist*c(-1/2,0,0)     # (-1.5, 0.0, 0.0)
M[3,]<-dist*c(0,height,0)   # ( 0.0, 2.6, 0.0)
M[4,]<-dist*c(0,len,kor)    # ( 0.0, 0.9, 2.4)
sig<-matrix(1,moodi,d)
p0<-1/moodi
p<-p0*rep(1,moodi)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed)
   return(dendat)
}

if (!is.null(N)){
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


simukern<-function(n,d,seed,newvalue,index,delta,minim,h){
#
set.seed(seed)
dendat<-matrix(0,n,d)
volofatom<-prod(delta)
newvalue<-newvalue/sum(volofatom*newvalue)
#
for (i in 1:n){
   uni<-runif(1)
   cumyla<-newvalue[1]*volofatom
   run<-1
   while (cumyla<uni){
      run<-run+1 
      cumyla<-cumyla+newvalue[run]*volofatom
   }
   #run leads to the right bin
   inde<-index[run,]
   uni2<-runif(1)                #add blur in bin
   obse<-minim-h+delta*inde+delta*(uni2-1/2)
   dendat[i,]<-obse
}
return(dendat)
}
slicing<-function(pcf,vecci,d1=1,d2=NULL)
# 1D slice: we calculate the slice parallel to the direction d1
# when d1=1, we calculate f(t,vecci) 
{
lenni<-length(pcf$value)
value<-matrix(0,lenni,1)
d<-length(pcf$N)
step<-matrix(0,d,1)
for (kk in 1:d) step[kk]<-(pcf$support[2*kk]-pcf$support[2*kk-1])/pcf$N[kk]

if  (is.null(d2)){  #1D slice

index<-matrix(0,lenni,1)

efek<-0
for (i in 1:lenni){

  currec<-matrix(0,2*d,1)
  for (kk in 1:d){
     currec[2*kk-1]<-pcf$support[2*kk-1]+pcf$down[i,kk]*step[kk]
     currec[2*kk]<-pcf$support[2*kk-1]+pcf$high[i,kk]*step[kk]
  }
  
  dimcal<-0
  onvalissa<-T
  j<-1
  while (j<=d){

     if (j!=d1){    
         ala<-currec[2*j-1]
         yla<-currec[2*j]
         if ((ala>vecci[j-dimcal]) || (yla<vecci[j-dimcal])) onvalissa<-F
     }
     else dimcal<-dimcal+1
     j<-j+1
  }
  if (onvalissa){
     efek<-efek+1
     value[efek]<-pcf$value[i]
     index[efek,1]<-pcf$high[i,d1]
  }
}

value<-value[1:efek]
index<-index[1:efek]
support<-pcf$support[(2*d1-1):(2*d1)]
N<-pcf$N[d1]
down<-matrix(0,efek,1)
high<-matrix(0,efek,1)
down[,1]<-index-1
high[,1]<-index
#down<-index-1
#high<-index

return(list(value=value,index=index,support=support,N=N,down=down,high=high))

}
else{ # 2D slice

if (is.null(d2)) d2<-2

down<-matrix(0,lenni,2)
high<-matrix(0,lenni,2)

efek<-0
for (i in 1:lenni){

  currec<-matrix(0,2*d,1)
  for (kk in 1:d){
     currec[2*kk-1]<-pcf$support[2*kk-1]+pcf$down[i,kk]*step[kk]
     currec[2*kk]<-pcf$support[2*kk-1]+pcf$high[i,kk]*step[kk]
  }
  
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
     value[efek]<-pcf$value[i]
     down[efek,1]<-pcf$down[i,d1]
     down[efek,2]<-pcf$down[i,d2] 
     high[efek,1]<-pcf$high[i,d1]
     high[efek,2]<-pcf$high[i,d2] 
  }
}

value<-value[1:efek]
down<-down[1:efek,]
high<-high[1:efek,]
support<-c(pcf$support[(2*d1-1):(2*d1)],pcf$support[(2*d2-1):(2*d2)])

return(list(value=value,down=down,high=high,
support=support,N=c(pcf$N[d1],pcf$N[d2])))
} #else 2D slice

}




sphere.map<-function(theta)
{
x<-matrix(0,2,1)
x[1]<-sin(theta)
x[2]<-cos(theta)

return(x)
}


sphere.para<-function(x)
{
d<-length(x)

if (d==2){
   if (x[1]>=0) theta<-acos(x[2])
   else theta<-acos(x[2])+pi
}

return(theta)
}

stseq<-function(N,lnum,
refe=NULL,func=NULL,dendat=NULL,
h=NULL,Q=NULL,kernel="epane",weights=NULL,
sig=rep(1,length(N)),support=NULL,theta=NULL,
M=NULL,p=NULL,mul=3,
t=rep(1,length(N)),marginal="normal",r=0,
mu=NULL,xi=NULL,Omega=NULL,alpha=NULL,df=NULL,g=1,
base=10
)
{
#lnum<-length(lseq)
level<-matrix(0,lnum,1)
volume<-matrix(0,lnum,1)
if (!is.null(dendat)) 
  pcf<-pcf.kern(dendat,h,N,kernel=kernel,weights=weights)
else
  pcf<-pcf.func(func,N,   #eval.func.dD
  sig=sig,support=support,theta=theta,
  M=M,p=p,mul=mul,
  t=t,marginal=marginal,r=r, 
  mu=mu,xi=xi,Omega=Omega,alpha=alpha,df=df,g=g)

maksi<-max(pcf$value)
l1<-maksi/(lnum+1) 
lmax<-maksi*lnum/(lnum+1)
level<-hgrid(l1,lmax,lnum,base=base)
level<-level[length(level):1]

for (i in 1:lnum){   
      #lev<-maksi*i/(lnum+1) 
      #level[i]<-lev 
      lev<-level[i]
      if (is.null(refe)) refe<-locofmax(pcf)
      st<-leafsfirst(pcf,lev=lev,refe=refe)
      volume[i]<-max(st$volume)
      if (i==1){
           if (lnum==1){ 
               istseq<-st
           }
           else{
               stseq<-list(st)
           }
      }
      else{
          stseq<-c(stseq,list(st))
      }
}
return(list(shtseq=stseq,level=level,volume=volume,pcf=pcf))
}



support<-function(dendat,epsi=0)
{
#estimoi kantajan tih datan perusteella
#dendat on n*xlkm matriisi
#epsi on tekn parametri
#kantajaksi estimoidaan [min-epsi,max+epsi]
#palauttaa xlkm*2-matriisin

xlkm<-length(dendat[1,])    #dendat matr sarakk lkm on muuttujien lkm
vast<-matrix(0,xlkm,2)  
i<-1
while (i<=xlkm){
    vast[i,1]<-min(dendat[,i])-epsi     #sis valien alkupisteet
    vast[i,2]<-max(dendat[,i])+epsi     #sis valien paatepisteet
    i<-i+1
}
return(vast)
}


tailfunc<-function(R,d,type,gnum=1000,sig=1,nu=1)
{
volball<-function(r,d){ return(r^d*pi^(d/2)/gamma(d/2+1)) }
volsphere<-function(d){ return(2*pi^(d/2)/gamma(d/2)) }

if (type=="bartlett"){
   norma<-d*(d+2)/(2*volsphere(d))
   funni<-function(t,d=d,nu=nu){ return( t^(d-1)*(1-t^2) ) }
   levfun<-function(t,d=d,sig=sig,nu=nu){ return( 1-(t/sig)^2 ) }
}
if (type=="gauss"){
   norma<-(2*pi)^(-d/2)
   funni<-function(t,d=d,nu=nu){ return( t^(d-1)*exp(-t^2/2) ) }
   levfun<-function(t,d=d,sig=sig,nu=nu){ return( exp(-(t/sig)^2/2) ) }
}
if (type=="student"){
   norma<-gamma((nu+d)/2)/((pi*nu)^(d/2)*gamma(nu/2))
   funni<-function(t,d=d,nu=nu){ return( t^(d-1)*(1+t^2/nu)^(-(d+nu)/2) ) }
   levfun<-function(t,d=d,sig=sig,nu=nu){ return( (1+(t/sig)^2/nu)^(-(d+nu)/2) ) }
}

# probability calc (numerical integral)
# y[r] = int_0^(r/sig) funni(t) dt
stepy<-R/sig/gnum
radiy<-seq(stepy,R/sig,stepy)
y<-matrix(0,length(radiy),1)
y[1]<-stepy*funni(radiy[1],d=d,nu=nu)
for (i in 2:length(y)){
    y[i]<-y[i-1]+stepy*funni(radiy[i],d=d,nu=nu)
}

# level calc
step<-R/gnum
radi<-seq(step,R,step)
level<-matrix(0,length(radi),1)
for (i in 1:length(level)){
    level[i]<-levfun(radi[i],d=d,sig=sig,nu=nu)
}

intrad2lev<-0
step<-radi[2]-radi[1]
for (i in 1:length(level)){
   intrad2lev<-intrad2lev+step*level[i]
}
levelnorma<-level/intrad2lev/2

proba<-norma*volsphere(d)*y
volu<-volball(radi,d)
level<-sig^(-d)*norma*level

return(list(radi=radi,proba=proba,volu=volu,level=level,levelnorma=levelnorma))
}









taillevel<-function(root,#child,sibling,
parent,volume,proba)
{
mt<-multitree(parent)
child<-mt$child
sibling<-mt$sibling

nodenum<-length(child)
level<-matrix(0,nodenum,1)
pino<-matrix(0,nodenum,1)

pino[1]<-root
pinin<-1
while (pinin>0){
      cur<-pino[pinin]      #take from stack
      pinin<-pinin-1

      chi<-child[cur]
      pare<-parent[cur]

      prochi<-0
      nexchi<-chi
      while (nexchi>0){
           prochi<-prochi+proba[nexchi]
           nexchi<-sibling[nexchi]
      }
      if (pare==0) levelpare<-0 else levelpare<-level[pare]

      level[cur]<-levelpare+(proba[cur]-prochi)/volume[cur]

      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
      }
      while (child[cur]>0){    #go to left and put right nodes to stack
            cur<-child[cur]

            chi<-child[cur]
            pare<-parent[cur]

            prochi<-0
            nexchi<-chi
            while (nexchi>0){
                 prochi<-prochi+proba[nexchi]
                 nexchi<-sibling[nexchi]
            }
            if (pare==0) levelpare<-0 else levelpare<-level[pare]

            level[cur]<-levelpare+(proba[cur]-prochi)/volume[cur]

            if (sibling[cur]>0){  #if candi has siblings
                pinin<-pinin+1
                pino[pinin]<-sibling[cur]
            } 
      }
}
return(level)

}





tail.plot.dens<-function(denmat,h=1,k=100,b=0.25,alpha=1,type="left.tail",
minx=-0.2,plot=TRUE)
{
# log="y",cex.axis=1,pch=20,pchs=rep(20,1000))
lkm<-dim(denmat)[2]
n<-dim(denmat)[1]

if (type=="left.tail"){

m<-floor(n/2)
detmat<-matrix(0,m,lkm)
for (i in 1:lkm){
       dencur<-denmat[,i]
       ordi<-order(dencur)
       dendat.ord<-dencur[ordi]
       detmat[,i]<-dendat.ord[1:m]
       #split<-median(dencur)
       #redu.ind<-(dencur<split) 
       #dendat.redu<-dencur[redu.ind]
       #ordi<-order(dendat.redu)
       #dendat.ord<-dendat.redu[ordi]  #nredu<-length(dendat.redu)
       #detmat[,i]<-dendat.ord[1:m]
}
minu<-min(detmat,na.rm=TRUE)
maki<-max(detmat,na.rm=TRUE)

x<-matrix(0,k,1)
pc<-matrix(0,k,m)
for (mm in 1:m){
     datai<-detmat[mm,]
     if (is.null(h)){
        expon<-1/(1+4)
        sdev<-sd(datai,na.rm=TRUE)
        h<-(4/(1+2))^expon*sdev*n^(-expon)
     }  
     ini<-!is.na(datai)
     dataj<-datai[ini]   
     for (kk in 1:k){
       arg<-minx+(maki-minx)*kk/(k+1) 
       x[kk]<-arg
       pc[kk,mm]<-kernesti.dens(arg,dataj,h=h)
   }
}
y<-log(seq(1,m))
pc2<-(pc)^b
colo<-grey(seq(1,0,-0.01),alpha=alpha)

if (plot) image(x,y,pc2,col=colo)  #image(pc2,col=topo.colors(120))
else return(list(x=x,y=y,pc=pc,colo=colo))
}


}


tail.plot<-function(dendat,type="both",split=median(dendat),
col="black",denmat=NULL,paletti=NULL,xlim=NULL,cex.axis=1,
pch=20,pchs=rep(20,1000),log="y")
{

if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])

if (type=="right.tail"){
   redu.ind<-(dendat>split) 
   dendat.redu<-dendat[redu.ind]
   ordi<-order(dendat.redu)
   dendat.ord<-dendat.redu[ordi]
   nredu<-length(dendat.redu)
   level<-seq(nredu,1)
   plot(dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,xlim=xlim,
   pch=pch)

   if (!is.null(denmat)){
      lkm<-dim(denmat)[2]
      for (i in 1:lkm){
          dencur<-denmat[,i]
          split=median(dencur)
          redu.ind<-(dencur>split) 
          dendat.redu<-dencur[redu.ind]
          ordi<-order(dendat.redu)
          dendat.ord<-dendat.redu[ordi]
          nredu<-length(dendat.redu)
          level<-seq(nredu,1)
          matplot(dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,
          add=TRUE,col=paletti[i],pch=pchs[i])
      }
   }
}

if (type=="left.tail"){
    redu.ind<-(dendat<split)
    dendat.redu<--dendat[redu.ind]
    ordi<-order(dendat.redu)
    dendat.ord<-dendat.redu[ordi]
    nredu<-length(dendat.redu)
    level<-seq(nredu,1)
    plot(-dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,xlim=xlim,
    pch=pch)

    if (!is.null(denmat)){
      lkm<-dim(denmat)[2]
      for (i in 1:lkm){
          dencur<-denmat[,i]
          split=median(dencur)
          redu.ind<-(dencur<split) 
          dendat.redu<--dencur[redu.ind]
          ordi<-order(dendat.redu)
          dendat.ord<-dendat.redu[ordi]
          nredu<-length(dendat.redu)
          level<-seq(nredu,1)
          matplot(-dendat.ord,level,log=log,xlab="",ylab="",cex.axis=cex.axis,
          add=TRUE,col=paletti[i],pch=pchs[i])
      }
    }
}

if (type=="both"){
    redu.ind<-(dendat<split)
    dendat.redu<--dendat[redu.ind]
    ordi<-order(dendat.redu)
    dendat.ord<-dendat.redu[ordi]
    nredu<-length(dendat.redu)
    level<-seq(nredu,1)

    redu.ind<-(dendat>split)
    dendat.redu<-dendat[redu.ind]
    ordi<-order(dendat.redu)
    dendat.ord.right<-dendat.redu[ordi]
    nredu<-length(dendat.redu)
    level.right<-seq(nredu,1)

    plot(c(-dendat.ord,dendat.ord.right),c(level,level.right),log=log,
    xlab="",ylab="",cex.axis=cex.axis,xlim=xlim,pch=pch)

    if (!is.null(denmat)){
      lkm<-dim(denmat)[2]
      for (i in 1:lkm){
          dencur<-denmat[,i]
          split=median(dencur)

          redu.ind<-(dencur<split) 
          dendat.redu<--dencur[redu.ind]
          ordi<-order(dendat.redu)
          dendat.ord<-dendat.redu[ordi]
          nredu<-length(dendat.redu)
          level<-seq(nredu,1)
          matplot(-dendat.ord,level,xlab="",ylab="",cex.axis=cex.axis,log=log,
          add=TRUE,col=paletti[i],pch=pchs[i])

          redu.ind<-(dencur>split) 
          dendat.redu<-dencur[redu.ind]
          ordi<-order(dendat.redu)
          dendat.ord<-dendat.redu[ordi]
          nredu<-length(dendat.redu)
          level<-seq(nredu,1)
          matplot(dendat.ord,level,xlab="",ylab="",cex.axis=cex.axis,log=log,
          add=TRUE,col=paletti[i],pch=pchs[i])
      }
    }
}

}











til1<-function(runi){
#
if (dim(t(runi))[1]==1) lkm<-1 else lkm<-length(runi[,1])
d<-length(runi[1,])/2
#
if (lkm>=2){
 parimat<-matrix(NA,lkm,lkm) #rivilla i lueteltu ne jotka leikk suorakaid i
                      # blokit tahan !!!!!!!!!!!!!!!!
 touchlkm<-matrix(0,lkm,1)   #kuinka monta kosketusta riv. i olevalle kaiteelle
# parimat2<-matrix(0,lkm,lkm) #rivilla i saralla j on 1 jos i ja j suork.leikk.
 l<-choose(lkm,2)
 curkosk<-matrix(0,l,2)      # blokit tahan !!!!!!!!!!!!!!!!
 currecs<-matrix(0,l,2*d)
 ind<-0
 for (i in 1:lkm){
   viite<-1 
   j<-i+1
   while (j<=lkm){
     ise<-leikkaa(runi[i,],runi[j,])
     if (!is.na(ise)){
       ind<-ind+1
       curkosk[ind,]<-c(i,j)
       currecs[ind,]<-ise 
       touchlkm[i]<-touchlkm[i]+1
       touchlkm[j]<-touchlkm[j]+1    
       parimat[i,touchlkm[i]]<-j
       parimat[j,touchlkm[j]]<-i
#       parimat2[i,j]<-1
#       parimat2[j,i]<-1
     }
     j<-j+1
   }
 }
}
if (ind==1){             #jos oli vain yksi leikkaus
  curkosk<-t(curkosk[1:ind,])
  currecs<-t(currecs[1:ind,])
  }
else if (ind>=2){      #jos oli useampi kuin yksi leikkaus
    curkosk<-curkosk[1:ind,]
    currecs<-currecs[1:ind,]
}
# supistetaan parimat
maxkosk<-max(touchlkm)
parimat<-parimat[,1:maxkosk]
if (maxkosk==1) parimat<-t(t(parimat))
return(list(ind=ind,curkosk=curkosk,currecs=currecs,parimat=parimat))
}






til2<-function(runi,curkosk,currecs,parimat,kosk){
#
blokki<-100
bloknum<-1
curpit<-blokki  # curpit<-lkm*curlkm #pahimillaan jokainen leikkaa jokaista
#
#if (dim(t(parimat))[1]==1) lkm<-1 else lkm<-length(parimat[1,]) 
lkm<-length(parimat[1,])    #kaiteiden maara
curlkm<-length(curkosk[,1]) #curlkm on edellisten kosketusten maara 
edkosk<-length(curkosk[1,]) #aikaisemmin haettiin kaikki leikkaukset
                            #edkosk:n kaiteen valilla  
d<-length(currecs[1,])/2
uuskosk<-matrix(0,curpit,kosk)     #matrix(0,lkm*curlkm,kosk) 
uusrecs<-matrix(0,curpit,2*d)      #matrix(0,lkm*curlkm,2*d)
#
ind<-0
for (i in 1:curlkm){  #kaydaan lapi kaikki nykyiset suorakaiteet
  vipu<-curkosk[i,edkosk] #uuden leikkaavan pitaa leikata esim viim curkosk:ssa
                  #to fasten the algorithm we do not consider every rec:
	          #only those who intersec with 1. in curkosk 
  ehdind<-1      
  ehdokas<-parimat[vipu,ehdind]  #ehdokkaat ovat ne jotka leikkaavat vipua
  while ((!is.na(ehdokas)) && (ehdind<=lkm)){ 
                  #kayd lapi ne jotka leikk vipua
                  #ehdokkaan pitaa leikata kaikkia muitakin 
                  #curkosk:n i:nnella rivilla olevia
    if (ehdokas>vipu){  #hetaan vain suuremmista kuin vipu   
      j<-1     
      touch<-TRUE
      olimuita<-FALSE
      while ((j<=(edkosk-1)) && (touch)){ #kayd lapi muut kuin vipu
        muu<-curkosk[i,j]
        if (!(ehdokas==muu)){
          olimuita<-TRUE
          curkoske<-parimat[ehdokas,]   #ne joihin ehdokas koskettaa
          touch<-onko(curkoske,muu)     #onko muu rivilla "curkoske"  
           #if (parimat2[ehdokas,muu]==0) touch<-FALSE
        }
              #jos ehdokas ja muu eivat kosketa ja ovat eri
	      #jos ehdokas=muu, niin parimat2[ehdokas,muu]=0  
        j<-j+1
      }
      if ((touch) && !(olimuita)) touch<-FALSE
      if (touch){  #jos ehdokas kosketti kaikkia muita
          ind<-ind+1   #lisataan uusien leikkausten laskuria
          if (ind<=curpit){  #jos ei tarvita uutta blokkia     
            uuskosk[ind,]<-c(curkosk[i,],ehdokas)  #???????
            uusrecs[ind,]<-leikkaa(currecs[i,],runi[ehdokas,])
          }
          else{
            bloknum<-bloknum+1
            uuspit<-bloknum*blokki
            apukosk<-matrix(0,uuspit,kosk)
            apurecs<-matrix(0,uuspit,2*d)
            apukosk[1:curpit,]<-uuskosk[1:curpit,]
            apurecs[1:curpit,]<-uusrecs[1:curpit,]
            apukosk[ind,]<-c(curkosk[i,],ehdokas)  #???????
            apurecs[ind,]<-leikkaa(currecs[i,],runi[ehdokas,])
            uuskosk<-apukosk
            uusrecs<-apurecs
            curpit<-uuspit
          }             
      }
    }  
    ehdind<-ehdind+1
    if (ehdind<=lkm) ehdokas<-parimat[vipu,ehdind] 
                                     #otetaan uusi vipua leikkaava
  }
}
if (ind>0){ 
       curkosk<-uuskosk[1:ind,]
       currecs<-uusrecs[1:ind,]
}     
return(list(ind=ind,currecs=currecs,curkosk=curkosk))
}








til<-function(runi){
#
if (dim(t(runi))[1]==1) lkm<-1 else lkm<-length(runi[,1])
masses<-matrix(0,lkm,1)
#
masses[1]<-sum(massat(runi))    #kaiteiden massojen summa
#
if (lkm>=2){
 apu<-til1(runi)
 ind<-apu$ind
 curkosk<-apu$curkosk
 currecs<-apu$currecs
 parimat<-apu$parimat
 if (ind>0){ #jos oli parittaisia leikkauksia
    masses[2]<-sum(massat(currecs)) #parittaisten leikkausten massojen summa
    kosk<-3
    while (ind>1){
      write(ind,file="apu",append=TRUE)
      apu2<-til2(runi,curkosk,currecs,parimat,kosk)
      ind<-apu2$ind
      if (ind>0){
        currecs<-apu2$currecs
        curkosk<-apu2$curkosk
        masses[kosk]<-sum(massat(currecs))
      }
      kosk<-kosk+1
    }
 }
}
res<-0                 # res<-til3(masses)
for (i in 1:lkm){
  res<-res+(-1)^(i-1)*masses[i]
}                                 
return(res)
}
touchi.boundary<-function(rec1,rec2,rho=0)
{
#Checks whether rectangles rec1, rec2 touch.
#rec1,rec2 are 2*d vectors, discrete rectangles (grid)

#Returns 0 if intersection is empty

d<-length(rec1)/2
if (length(rho)==1) rho<-rep(rho,d)

tulos<-1
i<-1
while ((i<=d) && (tulos==1)){  
    ala<-max(rec1[2*i-1],rec2[2*i-1])
    yla<-min(rec1[2*i],rec2[2*i])

    ala2<-min(rec1[2*i-1],rec2[2*i-1])
    yla2<-max(rec1[2*i],rec2[2*i])
    if ((ala2==0)&&(yla2==2*pi)) isboundary<-TRUE
    else isboundary<-FALSE

    if ((!isboundary)&&(yla+2*rho[i]<ala)) tulos<-0
    i<-i+1
}
return(tulos)
}




touchi.dela<-function(rec1,rec2,cate,dendat)
{
# returns 0 if intersection is empty.
# rec1 is (d+1)-vector
# if cate=simplex, then rec2 is (d+1)-vector (simplex) 
# if cate=rec,     then rec2 is 2d vector (rectangle)

d<-length(rec1)-1

if (cate=="rec"){  # rec2 is 2*d vector (rectangle)
 
    # make a bounding box of the simplex
    rec<-matrix(0,2*d,1)
    vertices<-matrix(0,d+1,d)
    for (dd in 1:(d+1)) vertices[dd,]<-dendat[rec1[dd]]    
    for (dd in 1:d){
        rec[2*dd-1]<-min(vertices[,dd])
        rec[2*dd]<-max(vertices[,dd])
    }
    # compare rec and rec2
    tulos<-1
    i<-1
    while ((i<=d) && (tulos==1)){
       ala<-max(rec[2*i-1],rec2[2*i-1])
       yla<-min(rec[2*i],rec2[2*i])
       if (yla<ala) tulos<-0
       i<-i+1
    }

}
else{    # comparison of simpleces: rec2 is d+1-vector

   tulos<-0
   i<-1
   while ( (i<=(d+1)) && (tulos==0) ){
      v1<-rec1[i]
      j<-1
      while ( (j<=(d+1)) && (tulos==0) ){
        v2<-rec2[j]
        if (v1==v2) tulos<-1
        j<-j+1
      }
      i<-i+1
   }

   #simp1<-dendat[rec1,]
   #simp2<-dendat[rec2,] 
   #if (tulos==0) tulos<-intersec.simpces(simp1,simp2)

}

return(tulos)
}

touchi<-function(rec1,rec2,rho=0)
{
#Checks whether rectangles rec1, rec2 touch.
#rec1,rec2 are 2*d vectors, discrete rectangles (grid)

#Returns 0 if intersection is empty

d<-length(rec1)/2
if (length(rho)==1) rho<-rep(rho,d)

tulos<-1
i<-1
while ((i<=d) && (tulos==1)){  
    ala<-max(rec1[2*i-1],rec2[2*i-1])
    yla<-min(rec1[2*i],rec2[2*i])
    if (yla+2*rho[i]<ala) tulos<-0
    i<-i+1
}
return(tulos)
}




touchi.simp<-function(rec1,rec2,cate,dendat)
{
# returns 0 if intersection is empty.
# rec1 is (d+1)-vector
# if cate=simplex, then rec2 is (d+1)-vector (simplex) 
# if cate=rec,     then rec2 is 2d vector (rectangle)

d<-length(rec1)-1

if (cate=="rec"){  # rec2 is 2*d vector (rectangle)
 
    # make a bounding box of the simplex
    rec<-matrix(0,2*d,1)
    vertices<-matrix(0,d+1,d)
    for (dd in 1:(d+1)) vertices[dd,]<-dendat[rec1[dd]]    
    for (dd in 1:d){
        rec[2*dd-1]<-min(vertices[,dd])
        rec[2*dd]<-max(vertices[,dd])
    }
    # compare rec and rec2
    tulos<-1
    i<-1
    while ((i<=d) && (tulos==1)){
       ala<-max(rec[2*i-1],rec2[2*i-1])
       yla<-min(rec[2*i],rec2[2*i])
       if (yla<ala) tulos<-0
       i<-i+1
    }

}
else{    # comparison of simpleces: rec2 is d+1-vector

   tulos<-0
   i<-1
   while ( (i<=(d+1)) && (tulos==0) ){
      v1<-rec1[i]
      j<-1
      while ( (j<=(d+1)) && (tulos==0) ){
        v2<-rec2[j]
        if (v1==v2) tulos<-1
        j<-j+1
      }
      i<-i+1
   }

   simp1<-dendat[rec1,]
   simp2<-dendat[rec2,] 
   if (tulos==0) tulos<-intersec.simpces(simp1,simp2)

}

return(tulos)
}

touchi.tail<-function(rec1,rec2,r1,r2=NULL,dist.type="euclid")
{
# Returns 0 if intersection is empty.
# rec1 is d-vector
# rec2 is d-vector or 2*d vector (rectangle)

d<-length(rec1)

if (dist.type=="euclid"){

if (length(rec2)==2*d){  # rec2 is 2*d vector (rectangle)
 
   point<-rec1
   rec<-rec2
   dist<-0
   for (i in 1:d){
      if (point[i]>rec[2*i]) 
          dist<-dist+(point[i]-rec[2*i])^2
      else if (point[i]<rec[2*i-1]) 
          dist<-dist+(point[i]-rec[2*i-1])^2
   }
   dist<-sqrt(dist)
   if (dist>r1) tulos<-0 else tulos<-1

}
else{                    # rec2 is d-vector

   dista<-sqrt(sum((rec1-rec2)^2))
   if (dista>r1+r2) tulos<-0 else tulos<-1

}

}
else{   # dist.type=="recta"

if (length(rec2)==2*d){  # rec2 is 2*d vector (rectangle)

    tulos<-1
    i<-1
    while ((i<=d) && (tulos==1)){
       ala<-max(rec1[i]-r1,rec2[2*i-1])
       yla<-min(rec1[i]+r1,rec2[2*i])
       if (yla<ala) tulos<-0
       i<-i+1
    }

}
else{                    # rec2 is d-vector

    tulos<-1
    i<-1
    while ((i<=d) && (tulos==1)){
       ala<-max(rec1[i]-r1,rec2[i]-r2)
       yla<-min(rec1[i]+r1,rec2[i]+r2)
       if (yla<ala) tulos<-0
       i<-i+1
    }

}

}

return(tulos)
}






touch<-function(rec1,rec2,epsi=0.000001)
{
# Checks whether rectangles rec1, rec2 touch.
# rec1,rec2 are 2*d vectors

# Returns FALSE if intersection is empty

d<-length(rec1)/2
tulos<-TRUE
i<-1
while ((i<=d) && (tulos)){  
    ala<-max(rec1[2*i-1],rec2[2*i-1])
    yla<-min(rec1[2*i],rec2[2*i])
    if (yla+epsi<ala) tulos<-FALSE
    i<-i+1
}
return(tulos)
}




touchstep.boundary<-function(node,curroot,boundrec,child,sibling,infopointer,
low,upp,rho=0)
{
# Checks whether "node" touches some of the leafs of the branch whose
# root is "curroot". Goes through the branch starting at "curroot".
# "comprec" is associated with the "node"
# "currec" is the bounding box of "cur"
# "pointrec" is associated with "cur"

d<-length(low[1,])
comprec<-matrix(0,2*d,1)
note<-infopointer[node]   #nodefinder[infopointer[node]]
for (i in 1:d){
    comprec[2*i-1]<-low[note,i]   #index[infopointer[node],i]
    comprec[2*i]<-upp[note,i]     #index[infopointer[node],i]
}

itemnum<-length(child)
pino<-matrix(0,itemnum,1)

potetouch<-1
istouch<-0
pino[1]<-curroot
pinin<-1
while ((pinin>0) && (istouch==0)){
      cur<-pino[pinin]      #take from stack
      pinin<-pinin-1

      # create currec and pointrec
      currec<-boundrec[cur,]
      pointrec<-matrix(0,2*d,1)
      note<-infopointer[cur]   #nodefinder[infopointer[cur]]
      for (i in 1:d){
         pointrec[2*i-1]<-low[note,i] #index[infopointer[cur],i]
         pointrec[2*i]<-upp[note,i]   #index[infopointer[cur],i]
      }
      # find touches                        
      potetouch<-touchi.boundary(comprec,currec,rho) 
      istouch<-touchi.boundary(comprec,pointrec,rho)

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
      }
      # go to left and put right nodes to the stack
      while ((child[cur]>0) && (istouch==0) && (potetouch==1)){
            cur<-child[cur]

            # create currec and pointrec
            currec<-boundrec[cur,]
            pointrec<-matrix(0,2*d,1)
            note<-infopointer[cur]   #nodefinder[infopointer[cur]]
            for (i in 1:d){
               pointrec[2*i-1]<-low[note,i] #index[infopointer[cur],i]
               pointrec[2*i]<-upp[note,i]   #index[infopointer[cur],i]
            }
            # find touches                        
            potetouch<-touchi.boundary(comprec,currec,rho) 
            istouch<-touchi.boundary(comprec,pointrec,rho)
 
            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
            }
      }
}

return(istouch)
}                                    
touchstep.complex<-function(node,curroot,boundrec,child,sibling,infopointer,
low,upp,dendat,complex)
{
# Checks whether "node" touches some of the leafs of the branch whose
# root is "curroot". Goes through the branch starting at "curroot".
# "comprec" is associated with the "node"
# "currec" is the bounding box of "cur"
# "pointrec" is associated with "cur"

d<-length(low[1,])

note<-infopointer[node]   #nodefinder[infopointer[node]]
comprec<-complex[note,]

itemnum<-length(child)
pino<-matrix(0,itemnum,1)
potetouch<-1
istouch<-0
pino[1]<-curroot
pinin<-1
while ((pinin>0) && (istouch==0)){
      cur<-pino[pinin]      #take from stack
      pinin<-pinin-1

      # create currec and pointrec
      currec<-boundrec[cur,]
      note<-infopointer[cur]   #nodefinder[infopointer[cur]]
      pointrec<-complex[note,]
      # find touches  
      potetouch<-touchi.simp(comprec,currec,cate="rec",dendat) 
      istouch<-touchi.simp(comprec,pointrec,cate="simplex",dendat)

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
      }
      # go to left and put right nodes to the stack
      while ((child[cur]>0) && (istouch==0) && (potetouch==1)){
            cur<-child[cur]

            # create currec and pointrec
            currec<-boundrec[cur,]
            note<-infopointer[cur]   #nodefinder[infopointer[cur]]
            pointrec<-complex[note,]
            # find touches                        
            potetouch<-touchi.simp(comprec,currec,cate="rec",dendat) 
            istouch<-touchi.simp(comprec,pointrec,cate="simplex",dendat)
 
            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
            }
      }
}

return(istouch)
}                                    

touchstep.delaunay<-function(node,curroot,boundrec,child,sibling,infopointer,
low,upp,dendat,complex)
{
# Checks whether "node" touches some of the leafs of the branch whose
# root is "curroot". Goes through the branch starting at "curroot".
# "comprec" is associated with the "node"
# "currec" is the bounding box of "cur"
# "pointrec" is associated with "cur"

d<-length(low[1,])

note<-infopointer[node]   #nodefinder[infopointer[node]]
comprec<-complex[note,]

itemnum<-length(child)
pino<-matrix(0,itemnum,1)
potetouch<-1
istouch<-0
pino[1]<-curroot
pinin<-1
while ((pinin>0) && (istouch==0)){
      cur<-pino[pinin]      #take from stack
      pinin<-pinin-1

      # create currec and pointrec
      currec<-boundrec[cur,]
      note<-infopointer[cur]   #nodefinder[infopointer[cur]]
      pointrec<-complex[note,]
      # find touches  
      potetouch<-touchi.dela(comprec,currec,cate="rec",dendat) 
      istouch<-touchi.dela(comprec,pointrec,cate="simplex",dendat)

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
      }
      # go to left and put right nodes to the stack
      while ((child[cur]>0) && (istouch==0) && (potetouch==1)){
            cur<-child[cur]

            # create currec and pointrec
            currec<-boundrec[cur,]
            note<-infopointer[cur]   #nodefinder[infopointer[cur]]
            pointrec<-complex[note,]
            # find touches                        
            potetouch<-touchi.dela(comprec,currec,cate="rec",dendat) 
            istouch<-touchi.dela(comprec,pointrec,cate="simplex",dendat)
 
            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
            }
      }
}

return(istouch)
}                                    

touchstep<-function(node,curroot,boundrec,child,sibling,infopointer,
low,upp,rho=0)
{
# Checks whether "node" touches some of the leafs of the branch whose
# root is "curroot". Goes through the branch starting at "curroot".
# "comprec" is associated with the "node"
# "currec" is the bounding box of "cur"
# "pointrec" is associated with "cur"

d<-length(low[1,])
comprec<-matrix(0,2*d,1)
note<-infopointer[node]   #nodefinder[infopointer[node]]
for (i in 1:d){
    comprec[2*i-1]<-low[note,i]   #index[infopointer[node],i]
    comprec[2*i]<-upp[note,i]     #index[infopointer[node],i]
}

itemnum<-length(child)
pino<-matrix(0,itemnum,1)

potetouch<-1
istouch<-0
pino[1]<-curroot
pinin<-1
while ((pinin>0) && (istouch==0)){
      cur<-pino[pinin]      #take from stack
      pinin<-pinin-1

      # create currec and pointrec
      currec<-boundrec[cur,]
      pointrec<-matrix(0,2*d,1)
      note<-infopointer[cur]   #nodefinder[infopointer[cur]]
      for (i in 1:d){
         pointrec[2*i-1]<-low[note,i] #index[infopointer[cur],i]
         pointrec[2*i]<-upp[note,i]   #index[infopointer[cur],i]
      }
      # find touches                        
      potetouch<-touchi(comprec,currec,rho) 
      istouch<-touchi(comprec,pointrec,rho)

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
      }
      # go to left and put right nodes to the stack
      while ((child[cur]>0) && (istouch==0) && (potetouch==1)){
            cur<-child[cur]

            # create currec and pointrec
            currec<-boundrec[cur,]
            pointrec<-matrix(0,2*d,1)
            note<-infopointer[cur]   #nodefinder[infopointer[cur]]
            for (i in 1:d){
               pointrec[2*i-1]<-low[note,i] #index[infopointer[cur],i]
               pointrec[2*i]<-upp[note,i]   #index[infopointer[cur],i]
            }
            # find touches                        
            potetouch<-touchi(comprec,currec,rho) 
            istouch<-touchi(comprec,pointrec,rho)
 
            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
            }
      }
}

return(istouch)
}                                    
touchstep.tail<-function(node,curroot,boundrec,child,sibling,infopointer,
low,upp,rho,dendat,dist.type="euclid")
{
# Checks whether "node" touches some of the leafs of the branch whose
# root is "curroot". Goes through the branch starting at "curroot".
# "comprec" is associated with the "node"
# "currec" is the bounding box of "cur"
# "pointrec" is associated with "cur"

d<-length(low[1,])

note<-infopointer[node]   #nodefinder[infopointer[node]]
comprec<-dendat[note,]
r1<-rho[note]

itemnum<-length(child)
pino<-matrix(0,itemnum,1)
potetouch<-1
istouch<-0
pino[1]<-curroot
pinin<-1
while ((pinin>0) && (istouch==0)){
      cur<-pino[pinin]      #take from stack
      pinin<-pinin-1

      # create currec and pointrec
      currec<-boundrec[cur,]
      note<-infopointer[cur]   #nodefinder[infopointer[cur]]
      pointrec<-dendat[note,]
      # find touches  
      r2<-rho[note]       
      potetouch<-touchi.tail(comprec,currec,r1,dist.type=dist.type) 
      istouch<-touchi.tail(comprec,pointrec,r1,r2,dist.type=dist.type)

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
      }
      # go to left and put right nodes to the stack
      while ((child[cur]>0) && (istouch==0) && (potetouch==1)){
            cur<-child[cur]

            # create currec and pointrec
            currec<-boundrec[cur,]
            note<-infopointer[cur]   #nodefinder[infopointer[cur]]
            pointrec<-dendat[note,]
            # find touches                        
            r2<-rho[note]
            potetouch<-touchi.tail(comprec,currec,r1,dist.type=dist.type) 
            istouch<-touchi.tail(comprec,pointrec,r1,r2,dist.type=dist.type)
 
            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
            }
      }
}

return(istouch)
}                                    
toucrec<-function(atoms,alkublokki,blokki){
#Finds which atoms touch each other
#
#items is atomnum*(2*d)-matrix
#alkublokki is an estimate to the maximum number of touches
#
#Returns links, atomnum*maxtouches-matrix
#
if (dim(t(atoms))[1]==1) m<-1 else m<-length(atoms[,1]) #m is number of atoms
len<-alkublokki
links<-matrix(NA,m,len)
maara<-matrix(0,m,1)
# merkitaan kosketukset linkit-matriisiin
i<-1
while (i<=m){
  j<-i+1
  while (j<=m){
    rec1<-atoms[i,]
    rec2<-atoms[j,]
    crit<-touch(rec1,rec2)
    if (crit){ #jos suorakaiteet koskettavat
        maari<-maara[i]+1
        maarj<-maara[j]+1
        if ((maari>len) || (maarj>len)){
            links<-blokitus2(links,blokki)
            len<-len+blokki
        }
        links[i,maari]<-j
        maara[i]<-maari
        links[j,maarj]<-i
        maara[j]<-maarj         
    }
    j<-j+1 
  }
  i<-i+1
} 
return(links)
}






travel.tree<-function(parent,node)
{
mt<-multitree(parent) #roots<-mt$roots child<-mt$child sibling<-mt$sibling
itemnum<-length(parent)
nodes<-matrix(0,itemnum,1)

curroot<-node
counter<-0
if (mt$child[curroot]>0){
   pino<-matrix(0,itemnum,1)
   pino[1]<-mt$child[curroot]
   pinin<-1
   while (pinin>0){
        cur<-pino[pinin]      #take from stack
        pinin<-pinin-1
        counter<-counter+1
        nodes[counter]<-cur   
        # put to the stack 
        if (mt$sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-mt$sibling[cur]
        }
        # go to left and put right nodes to the stack
        while (mt$child[cur]>0){   
            cur<-mt$child[cur]
            counter<-counter+1
            nodes[counter]<-cur   
            if (mt$sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-mt$sibling[cur]
            }
        }
   }#while (pinin>0)
}
       
nodes<-nodes[1:counter]
return(nodes)
}

treedisc.ada<-function(lst, pcf, ngrid=NULL, r=NULL, type=NULL, lowest="dens")
{
# r is vector of radiuses, we prune shapetree "lst" so that
# its radiuses are given by r

if (lowest=="dens") lowest<-0 else lowest<-min(lst$level)

if (is.null(type)){
   if (is.null(lst$refe)) type<-"lst"
   else type<-"shape"
}

if (is.null(r)){
  if (type=="shape"){
      stepsi<-lst$maxdis/ngrid
      r<-seq(0,lst$maxdis,stepsi)
  }
  else{  #type=="lst"
      stepsi<-lst$maxdis/(ngrid+1)    
      r<-seq(lowest+stepsi,lst$maxdis-stepsi,stepsi)
  }
}

mt<-multitree(lst$parent)
child<-mt$child
sibling<-mt$sibling

d<-dim(lst$center)[1]
itemnum<-length(lst$parent)

################################################

parent<-matrix(NA,itemnum,1)

pino<-matrix(0,itemnum,1)
pinoparent<-matrix(0,itemnum,1)
pinorad<-matrix(0,itemnum,1)

pino[1]<-1
pinoparent[1]<-0
pinorad[1]<-1
pinin<-1
curradind<-1

while (pinin>0){ # && (curradind<=length(r))){
      cur<-pino[pinin]      #take from stack
      curpar<-pinoparent[pinin]
      curradind<-pinorad[pinin]
      pinin<-pinin-1

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
            pinoparent[pinin]<-curpar
            pinorad[pinin]<-curradind
      }

      note<-lst$infopointer[cur] #cur
      if (type=="lst")
         etai<-pcf$value[note]
      else{
         recci<-matrix(0,2*d,1)
         downi<-pcf$down[lst$infopointer[note],]
         highi<-pcf$high[lst$infopointer[note],]
         for (jj in 1:d){
             recci[2*jj-1]<-pcf$grid[downi[jj],jj]
             recci[2*jj]<-pcf$grid[highi[jj],jj]
         }
         etai<-sqrt(etaisrec(lst$refe,recci))
      }

      if (curradind<=length(r)) currad<-r[curradind] else currad<-1000000
      if (etai>currad){
          parent[cur]<-curpar
          curpar<-cur
          curradind<-curradind+1
      }

      # go to left and put right nodes to the stack
      while (child[cur]>0){  # && (curradind<=length(r))){
            cur<-child[cur]

            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
                 pinoparent[pinin]<-curpar
                 pinorad[pinin]<-curradind
            }
 
            note<-lst$infopointer[cur] #cur
            if (type=="lst")
                etai<-pcf$value[note]
            else{
                recci<-matrix(0,2*d,1)
                downi<-pcf$down[lst$infopointer[note],]
                highi<-pcf$high[lst$infopointer[note],]
                for (jj in 1:d){
                    recci[2*jj-1]<-pcf$grid[downi[jj],jj]
                    recci[2*jj]<-pcf$grid[highi[jj],jj]
                }
                etai<-sqrt(etaisrec(lst$refe,recci))
            }

            if (curradind<=length(r)) currad<-r[curradind] else currad<-1000000
            if (etai>currad){
                parent[cur]<-curpar
                curpar<-cur
                curradind<-curradind+1 
            }
 
      }

}

# Prune ##################################

newparent<-matrix(0,itemnum,1)
newcenter<-matrix(0,d,itemnum)
newvolume<-matrix(0,itemnum,1)
newlevel<-matrix(0,itemnum,1)
newpointer<-matrix(0,itemnum,1)
#newdistcenter<-matrix(0,d,itemnum)
#newproba<-matrix(0,itemnum,1)

i<-1
newlkm<-0
while (i<=itemnum){
    if (!is.na(parent[i])){
         newlkm<-newlkm+1
         newpointer[i]<-newlkm
         if (parent[i]==0)  newparent[newlkm]<-0
         else newparent[newlkm]<-newpointer[parent[i]]
         newcenter[,newlkm]<-lst$center[,i]
         newlevel[newlkm]<-lst$level[i]
         newvolume[newlkm]<-lst$volume[i]
         #newdistcenter[,newlkm]<-lst$distcenter[,i]
         #newproba[newlkm]<-lst$proba[i]
    }
    i<-i+1
}

newparent<-newparent[1:newlkm]
if (newlkm<=1) newcenter<-matrix(newcenter[,1],d,1) 
else newcenter<-newcenter[,1:newlkm]
newvolume<-newvolume[1:newlkm]
newlevel<-newlevel[1:newlkm]
#if (newlkm<=1) newdistcenter<-matrix(newdistcenter[,1],d,1) 
#else newdistcenter<-newdistcenter[,1:newlkm]
#newproba<-newproba[1:newlkm]
newpointer<-newpointer[1:newlkm]

return(list(parent=newparent,level=newlevel,volume=newvolume,center=newcenter,
#distcenter=newdistcenter,  #branchradius=newbranchradius,
#proba=newproba,
refe=lst$refe,bary=lst$bary,root=1,infopointer=newpointer))

}   

treedisc.intpol<-function(lst, f, ngrid=NULL, r=NULL, lowest="dens")
{
# r is vector of radiuses, we prune shapetree "lst" so that
# its radiuses are given by r

if (lowest=="dens") lowest<-0 else lowest<-min(lst$level)

if (is.null(r)){
      stepsi<-(lst$maxdis-lowest)/(ngrid+1)    
      r<-seq(lowest+stepsi,lst$maxdis-stepsi,stepsi)
}

mt<-multitree(lst$parent)
child<-mt$child
sibling<-mt$sibling

d<-dim(lst$center)[1]
itemnum<-length(lst$parent)

################################################

parent<-matrix(NA,itemnum,1)

rootnum<-length(mt$root)
for (rr in 1:rootnum){

pino<-matrix(0,itemnum,1)
pinoparent<-matrix(0,itemnum,1)
pinorad<-matrix(0,itemnum,1)

pino[1]<-mt$root[rr]
pinoparent[1]<-0
pinorad[1]<-1
pinin<-1
curradind<-1

while (pinin>0){ # && (curradind<=length(r))){
      cur<-pino[pinin]      #take from stack
      curpar<-pinoparent[pinin]
      curradind<-pinorad[pinin]
      pinin<-pinin-1

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
            pinoparent[pinin]<-curpar
            pinorad[pinin]<-curradind
      }

      note<-lst$infopointer[cur] #cur
      
      etai<-f[note]

      if (curradind<=length(r)) currad<-r[curradind] else currad<-1000000
      if (etai>currad){
          parent[cur]<-curpar
          curpar<-cur
          curradind<-curradind+1
      }

      # go to left and put right nodes to the stack
      while (child[cur]>0){  # && (curradind<=length(r))){
            cur<-child[cur]

            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
                 pinoparent[pinin]<-curpar
                 pinorad[pinin]<-curradind
            }
 
            note<-lst$infopointer[cur] #cur

            etai<-f[note]

            if (curradind<=length(r)) currad<-r[curradind] else currad<-1000000
            if (etai>currad){
                parent[cur]<-curpar
                curpar<-cur
                curradind<-curradind+1 
            }
 
      }

}

}

# Prune ##################################

newparent<-matrix(0,itemnum,1)
newcenter<-matrix(0,d,itemnum)
newvolume<-matrix(0,itemnum,1)
newlevel<-matrix(0,itemnum,1)
newpointer<-matrix(0,itemnum,1)
newdistcenter<-matrix(0,d,itemnum)
newproba<-matrix(0,itemnum,1)

i<-1
newlkm<-0
while (i<=itemnum){
    if (!is.na(parent[i])){
         newlkm<-newlkm+1
         newpointer[i]<-newlkm
         if (parent[i]==0)  newparent[newlkm]<-0
         else newparent[newlkm]<-newpointer[parent[i]]
         newcenter[,newlkm]<-lst$center[,i]
         newlevel[newlkm]<-lst$level[i]
         newvolume[newlkm]<-lst$volume[i]
         #newdistcenter[,newlkm]<-lst$distcenter[,i]
         #newproba[newlkm]<-lst$proba[i]
    }
    i<-i+1
}

newparent<-newparent[1:newlkm]
if (newlkm<=1) newcenter<-matrix(newcenter[,1],d,1) 
else newcenter<-newcenter[,1:newlkm]
newvolume<-newvolume[1:newlkm]
newlevel<-newlevel[1:newlkm]
if (newlkm<=1) newdistcenter<-matrix(newdistcenter[,1],d,1) 
else newdistcenter<-newdistcenter[,1:newlkm]
newproba<-newproba[1:newlkm]
newpointer<-newpointer[1:newlkm]

return(list(parent=newparent,level=newlevel,volume=newvolume,center=newcenter,
distcenter=newdistcenter,  #branchradius=newbranchradius,
proba=newproba,
refe=lst$refe,bary=lst$bary,root=1,infopointer=newpointer))

}   

treedisc<-function(lst, pcf, ngrid=NULL, r=NULL, type=NULL, lowest="dens")
{
# r is vector of radiuses, we prune shapetree "lst" so that
# its radiuses are given by r

if (lowest=="dens") lowest<-0 else lowest<-min(lst$level)

if (is.null(type)){
   if (is.null(lst$refe)) type<-"lst"
   else type<-"shape"
}

if (is.null(r)){
  if (type=="shape"){
      stepsi<-lst$maxdis/ngrid
      r<-seq(0,lst$maxdis,stepsi)
  }
  else{  #type=="lst"
      stepsi<-lst$maxdis/(ngrid+1)    
      r<-seq(lowest+stepsi,lst$maxdis-stepsi,stepsi)
  }
}

mt<-multitree(lst$parent)
child<-mt$child
sibling<-mt$sibling

d<-dim(lst$center)[1]
itemnum<-length(lst$parent)

if (is.null(pcf$step)){
    step<-matrix(0,d,1)
    for (i in 1:d) step[i]=(pcf$support[2*i]-pcf$support[2*i-1])/pcf$N[i];
    pcf$step<-step
}

################################################

parent<-matrix(NA,itemnum,1)

pino<-matrix(0,itemnum,1)
pinoparent<-matrix(0,itemnum,1)
pinorad<-matrix(0,itemnum,1)

pino[1]<-1
pinoparent[1]<-0
pinorad[1]<-1
pinin<-1
curradind<-1

while (pinin>0){ # && (curradind<=length(r))){
      cur<-pino[pinin]      #take from stack
      curpar<-pinoparent[pinin]
      curradind<-pinorad[pinin]
      pinin<-pinin-1

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
            pinoparent[pinin]<-curpar
            pinorad[pinin]<-curradind
      }

      note<-lst$infopointer[cur] #cur
      if (type=="lst")
         etai<-pcf$value[note]
      else{
         recci<-matrix(0,2*d,1)
         for (jj in 1:d){
            recci[2*jj-1]<-pcf$support[2*jj-1]+pcf$step[jj]*pcf$down[note,jj]
            recci[2*jj]<-pcf$support[2*jj-1]+pcf$step[jj]*pcf$high[note,jj]  
         }
         etai<-sqrt(etaisrec(lst$refe,recci))
      }

      if (curradind<=length(r)) currad<-r[curradind] else currad<-1000000
      if (etai>currad){
          parent[cur]<-curpar
          curpar<-cur
          curradind<-curradind+1
      }

      # go to left and put right nodes to the stack
      while (child[cur]>0){  # && (curradind<=length(r))){
            cur<-child[cur]

            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
                 pinoparent[pinin]<-curpar
                 pinorad[pinin]<-curradind
            }
 
            note<-lst$infopointer[cur] #cur
            if (type=="lst")
                etai<-pcf$value[note]
            else{
                recci<-matrix(0,2*d,1)
                for (jj in 1:d){
                   recci[2*jj-1]<-pcf$support[2*jj-1]+pcf$step[jj]*pcf$down[note,jj]
                   recci[2*jj]<-pcf$support[2*jj-1]+pcf$step[jj]*pcf$high[note,jj]  
                }
                etai<-sqrt(etaisrec(lst$refe,recci))
            }

            if (curradind<=length(r)) currad<-r[curradind] else currad<-1000000
            if (etai>currad){
                parent[cur]<-curpar
                curpar<-cur
                curradind<-curradind+1 
            }
 
      }

}

#lst$roots<-c(1)
#lst$parent<-parent
#return(lst)

# Prune ##################################

newparent<-matrix(0,itemnum,1)
newcenter<-matrix(0,d,itemnum)
newvolume<-matrix(0,itemnum,1)
newlevel<-matrix(0,itemnum,1)
newpointer<-matrix(0,itemnum,1)
newdistcenter<-matrix(0,d,itemnum)
newproba<-matrix(0,itemnum,1)

#newparent[1]<-0
#newcenter[,1]<-lst$center[,1]
#newvolume[1]<-lst$volume[1]
#newlevel[1]<-lst$level[1]
#newpointer[1]<-1
#newdistcenter[,1]<-lst$distcenter[,1]

i<-1
newlkm<-0
while (i<=itemnum){
    if (!is.na(parent[i])){
         newlkm<-newlkm+1
         newpointer[i]<-newlkm
         if (parent[i]==0)  newparent[newlkm]<-0
         else newparent[newlkm]<-newpointer[parent[i]]
         newcenter[,newlkm]<-lst$center[,i]
         newlevel[newlkm]<-lst$level[i]
         newvolume[newlkm]<-lst$volume[i]
         newdistcenter[,newlkm]<-lst$distcenter[,i]
         newproba[newlkm]<-lst$proba[i]
    }
    i<-i+1
}

newparent<-newparent[1:newlkm]
if (newlkm<=1) newcenter<-matrix(newcenter[,1],d,1) 
else newcenter<-newcenter[,1:newlkm]
newvolume<-newvolume[1:newlkm]
newlevel<-newlevel[1:newlkm]
if (newlkm<=1) newdistcenter<-matrix(newdistcenter[,1],d,1) 
else newdistcenter<-newdistcenter[,1:newlkm]
newproba<-newproba[1:newlkm]
newpointer<-newpointer[1:newlkm]

return(list(parent=newparent,level=newlevel,volume=newvolume,center=newcenter,
distcenter=newdistcenter,  #branchradius=newbranchradius,
proba=newproba,
refe=lst$refe,bary=lst$bary,root=1,infopointer=newpointer))

}   

tree.segme<-function(tt,paletti=NULL,pcf=NULL)
{

if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])

colors<-colobary(tt$parent,paletti)
if (is.null(pcf)) segme<-colors
else{
  lenni<-length(pcf$value)
  segme<-matrix(0,lenni,1)
}
for (i in 1:length(colors)) segme[tt$infopointer[i]]<-colors[i]

return(segme)
}






vectomatch<-function(vec1,vec2)
{

d<-dim(vec1)[2]
prenum<-dim(vec1)[1]
curnum<-dim(vec2)[1]
parento<-matrix(0,curnum,1)

smallernum<-min(prenum,curnum)
greaternum<-max(prenum,curnum)

dista<-matrix(NA,smallernum,greaternum)
for (ap in 1:smallernum){
    for (be in 1:greaternum){
           if (d==1){
               if (prenum<=curnum){
               precenter<-vec1[ap]
               curcenter<-vec2[be]
               }
               else{
               precenter<-vec2[ap]
               curcenter<-vec1[be]
               }
           }
           else{
               if (prenum<=curnum){
               precenter<-vec1[ap,]
               curcenter<-vec2[be,]
               }
               else{
               precenter<-vec2[ap,]
               curcenter<-vec1[be,]
               }
           }
           dista[ap,be]<-etais(curcenter,precenter)
    }
}

match<-matrix(0,smallernum,1)  #for each mode the best match
findtie<-TRUE

# find the best match for all and check whether there are ties
match<-matrix(0,smallernum,1)
for (bm in 1:smallernum){
      minimi<-min(dista[bm,],na.rm=TRUE)
      match[bm]<-which(minimi==dista[bm,])[1]
}
findtie<-FALSE
bm<-1
while ((bm<=smallernum) && (findtie==FALSE)){
         koe<-match[bm]
         bm2<-bm+1
         while (bm2<=smallernum){
            if (koe==match[bm2]){
                  findtie<-TRUE
            }
            bm2<-bm2+1
         }
         bm<-bm+1
}
    
onkayty<-FALSE
tiematch<-matrix(0,smallernum,1)

while (findtie){

      onkayty<-TRUE
      
      # find the best match for all
      bestmatch<-matrix(0,smallernum,1)
      for (bm in 1:smallernum){
          allna<-TRUE
          am<-1
          while ((am<=greaternum) && (allna)){
             if (!is.na(dista[bm,am])) allna<-FALSE
             am<-am+1
          }
          if (!(allna)){
             minimi<-min(dista[bm,],na.rm=TRUE)
             bestmatch[bm]<-which(minimi==dista[bm,])[1]
          }
          else bestmatch[bm]<-tiematch[bm]
      }

      # find the first tie
      findtie<-FALSE

      tieset<-matrix(0,smallernum,1)
      bm<-1
      while ((bm<=smallernum) && (findtie==FALSE)){
         koe<-bestmatch[bm]
         bm2<-bm+1
         while (bm2<=smallernum){
            if (koe==bestmatch[bm2]){
                  findtie<-TRUE
                  tieset[bm]<-1
                  tieset[bm2]<-1
            }
            bm2<-bm2+1
         }
         bm<-bm+1
      }

      # solve the first tie
      if (findtie==TRUE){
         numofties<-sum(tieset)
         kavelija<-0
         tiepointer<-matrix(0,numofties,1) 
         # find the second best
         secondbest<-matrix(0,smallernum,1)
         for (bm in 1:smallernum){
            if (tieset[bm]==1){
               redudista<-dista[bm,]
               redudista[bestmatch[bm]]<-NA
               minimi<-min(redudista,na.rm=TRUE)
               secondbest[bm]<-which(minimi==redudista)[1]

               kavelija<-kavelija+1
               tiepointer[kavelija]<-bm
            }
         }
         # try different combinations       
         # try all subsets of size 2 from the set of ties
         numofsubsets<-choose(numofties,2)
            #gamma(numofties+1)/gamma(numofties-2+1)
         valuelist<-matrix(0,numofsubsets,1)
         vinnerlist<-matrix(0,numofsubsets,1)
         matchlist<-matrix(0,numofsubsets,1)
         runneri<-1
         eka<-1
         while (eka<=numofties){
            ekapo<-tiepointer[eka]
            toka<-eka+1
            while (toka<=numofties){
               tokapo<-tiepointer[toka]
               # try combinations for this subset (there are 2)
               # 1st combination
               fvinner<-ekapo
               fvinnermatch<-bestmatch[fvinner]
               floser<-tokapo
               flosermatch<-secondbest[floser]
               fvalue<-dista[fvinner,fvinnermatch]+dista[floser,flosermatch]
                # 2nd combination
               svinner<-tokapo
               svinnermatch<-bestmatch[svinner]
               sloser<-ekapo
               slosermatch<-secondbest[sloser]
               svalue<-dista[svinner,svinnermatch]+dista[sloser,slosermatch]
               # tournament
               if (fvalue<svalue){
                   valuelist[runneri]<-fvalue
                   vinnerlist[runneri]<-fvinner
                   matchlist[runneri]<-fvinnermatch
               }
               else{ 
                   valuelist[runneri]<-svalue
                   vinnerlist[runneri]<-svinner
                   matchlist[runneri]<-svinnermatch
               }
               runneri<-runneri+1 
               # 
               toka<-toka+1
            }
            eka<-eka+1
         }
         minimi<-min(valuelist,na.rm=TRUE)
         bestsub<-which(minimi==valuelist)[1]
         vinnerson<-vinnerlist[bestsub]
         matcherson<-matchlist[bestsub]

         tiematch[vinnerson]<-matcherson
         dista[vinnerson,]<-NA
         dista[,matcherson]<-NA

      }

}  #while (findtie)

if (onkayty){  #there was one tie
          
          for (sepo in 1:smallernum){
               if (tiematch[sepo]!=0) match[sepo]<-tiematch[sepo]
               else match[sepo]<-bestmatch[sepo]
          }
}

newnode<-matrix(0,curnum,1)
if (prenum>curnum) parento<-match
else{
    for (i in 1:prenum){   #kaannetaan linkit
        linko<-match[i]
        parento[linko]<-i
    }    
    for (j in 1:curnum){
        if (parento[j]==0){   #jos ei linkkia, haetaan lahin vanhemmaksi
             newnode[j]<-1    #we label the rest-labels
             distvec<-dista[,j]  #sarake antaa etaisyyden
             minimi<-min(distvec,na.rm=TRUE)
             parento[j]<-which(minimi==distvec)[1]
        }
    }
}

return(list(parent=parento,newnode=newnode))
}



















volball<-function(r,d){ return(r^d*pi^(d/2)/gamma(d/2+1)) }

vols.complex<-function(complex,dendat,meto="voltriangle")
{
# complex is lkm*(d+1) matrix
# dendat is n*d matrix

lkm<-dim(complex)[1]
vols<-matrix(0,lkm,1)

if (meto=="voltriangle"){
for (i in 1:lkm){
    ind<-complex[i,]
    simp<-dendat[ind,]
    vols[i]<-voltriangle(simp)
}}
else{
for (i in 1:lkm){
    ind<-complex[i,]
    simp<-dendat[ind,]
    vols[i]<-volsimplex(simp)
}}

return(vols)
}


volsimplex<-function(simp)
{
# simp is (d+1)*d matrix  / 3*2 matrix

d<-dim(simp)[2]
M<-matrix(0,d,d)
for (i in 1:d) M[i,]<-simp[i+1,]-simp[1,]
vol<-abs(det(M))/factorial(d)

#v1<-simp[1,]
#v2<-simp[2,]
#v3<-simp[3,]
#a<-sqrt( sum((v1-v2)^2) )
#b<-sqrt( sum((v1-v3)^2) )
#c<-sqrt( sum((v2-v3)^2) )
#s<-(a+b+c)/2
#vol<-sqrt( s*(s-a)*(s-b)*(s-c) )

return(vol)
}

voltriangle<-function(simp)
{
# simp is (d+1)*d matrix  / 3*2 matrix
# Heron's formula

v1<-simp[1,]
v2<-simp[2,]
v3<-simp[3,]

a<-sqrt( sum((v1-v2)^2) )
b<-sqrt( sum((v1-v3)^2) )
c<-sqrt( sum((v2-v3)^2) )

s<-(a+b+c)/2
vol<-sqrt( s*(s-a)*(s-b)*(s-c) )

return(vol)
}

weightsit<-function(n,h,katka=4)
{
#normvakio<-(sqrt(2*pi)*h)^{-1}
resu<-matrix(0,n,1)
zumma<-0
for (i in 1:n){
    eta<-(n-i)
    if (eta/h>katka) tulos<-0 else tulos<-exp(-eta^2/(2*h^2))#*normvakio
    resu[i]<-tulos
    zumma<-zumma+tulos
}

resu<-resu/zumma

return(resu)
}



