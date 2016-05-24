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














