no_bounded <-
function(Dx,C){
cp<-Dx%*%C
zERO<-1E-18
zeros<-0
negs<-0
big<-0
for(i in 1:length(cp))
if(abs(cp[i])>1e150||(abs(cp[i])<1e-150&&cp[i]!=0)) big <-1
for(i in 1:length(cp)){
if(abs(cp[i])<zERO) zeros<-zeros+1
if(cp[i]<0 && abs(cp[i])>=zERO) negs<-negs+1}
if((zeros!=length(cp) && (zeros+negs)==length(cp))||big==1) a<-1 else a<-0
return(a)}
