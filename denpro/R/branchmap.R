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





