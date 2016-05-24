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
