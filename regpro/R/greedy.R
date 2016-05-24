greedy<-function(x,y,arg,M,m,splitfreq=1)
{
d<-length(arg)
n<-length(y)
#s<-splitsearch(x,y,arg,splitfreq=splitfreq)

inx<-matrix(0,n*d+1,1)
for (i in 1:n){
    for (j in 1:d){
        inx[1+(i-1)*d+j]=x[i,j]
    }
}
iny<-matrix(0,n+1,1)
iny[2:(n+1)]<-y
inarg<-matrix(0,d+1,1)
inarg[2:(d+1)]<-arg
kg<-1
#kg<-.C("splitSearch",
#           as.double(inx),
#           as.double(iny),
#           as.double(inarg),
#           as.double(splitfreq),
#           as.integer(n),
#           as.integer(d),
#           indeces = integer(n+1),
#           lkm = integer(1)
#)
s<-kg$indeces[2:(kg$lkm+1)]

xnew<-matrix(x[s,],length(s),d)
ynew<-matrix(y[s],length(s),1)
num<-length(s)
lkm<-2
while ((lkm<=M) && (num>m)){
    #s<-splitsearch(xnew,ynew,arg,splitfreq=splitfreq)    

    n<-length(ynew)
    inx<-matrix(0,n*d+1,1)
    for (i in 1:n){
       for (j in 1:d){
          inx[1+(i-1)*d+j]=xnew[i,j]
       }
    }
    iny<-matrix(0,n+1,1)
    iny[2:(n+1)]<-ynew
    inarg<-matrix(0,d+1,1)
    inarg[2:(d+1)]<-arg
    #kg<-.C("splitSearch",
    #           as.double(inx),
    #           as.double(iny),
    #           as.double(inarg),
    #           as.double(splitfreq),
    #           as.integer(n),
    #           as.integer(d),
    #           indeces = integer(n+1),
    #           lkm = integer(1)
    #)
    s<-kg$indeces[2:(kg$lkm+1)]

    num<-length(s)
    if (num>=m){
       xnew<-matrix(xnew[s,],length(s),d)
       ynew<-matrix(ynew[s],length(s),1)
    }
    lkm<-lkm+1
}

val<-mean(ynew)
return(list(val=val,x=xnew,y=ynew))
}


