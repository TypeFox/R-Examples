
drgen<-function(fd,fr,e){

gpp<-1-prod((1-fd-fr)^2+2*fr*(1-fd-fr))

r<-matrix(NA,9,3) 
#1:mzt,2:parent-offspring,3:dzt,4:sibling,5:2direct,6:3rd,7:3d,8:4th,9:4d

r[1,1]<-1-gpp+gpp*(1-e)^2
r[1,2]<-gpp*2*e*(1-e)
r[1,3]<-gpp*e^2

gnn<-prod(2*fr*(1-fr-fd)*(1/2*(1-fr-fd)+1/2*(1-fd))+(1-fr-fd)^2*(1-fd))
gnd<-prod(fd^2*0+2*fd*(1-fr-fd)*1/2*(1-fd)+2*fd*fr*1/2*(1-fr-fd)+fr^2*(1-fr-fd)+2*fr*(1-fr-fd)*(1/2*(1-fr-fd)+1/2*(1-fd))+(1-fr-fd)^2*(1-fd))-prod((1-fr-fd)^2*(1-fd)+2*fr*(1-fr-fd)*(1/2*(1-fr-fd)+1/2*(1-fd)))+prod((1-fr-fd)^2+2*fr*(1-fr-fd))-gnn
t1<-c(fd[1]^2,2*fd[1]*fr[1],2*fd[1]*(1-fr[1]-fd[1]),fr[1]^2,2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(fd[i]^2,2*fd[i]*fr[i],2*fd[i]*(1-fr[i]-fd[i]),fr[i]^2,2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(0,1/2*(1-fr[1]-fd[1]),1/2*(1-fd[1]),1-fr[1]-fd[1],1/2*(1-fr[1]-fd[1])+1/2*(1-fd[1]),1-fd[1])
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(0,1/2*(1-fr[i]-fd[i]),1/2*(1-fd[i]),1-fr[i]-fd[i],1/2*(1-fr[i]-fd[i])+1/2*(1-fd[i]),1-fd[i]))
temp1<-t1*(1-t2)
t1<-c(2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/2*(1-fr[1]-fd[1])+1/2*(1-fd[1]),1-fd[1])
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/2*(1-fr[i]-fd[i])+1/2*(1-fd[i]),1-fd[i]))
temp2<-t1*(1-t2)
gdd<-sum(temp1)-sum(temp2)
r[2,1]<-gnn+gnd*(1-e)+gdd*(1-e)^2
r[2,2]<-gnd*e+gdd*2*e*(1-e)
r[2,3]<-gdd*e^2

gnn<-prod(2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))
gnd<-prod(fd^2*1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/4*(1-fr-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(2/4*(1-fr-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-gnn
t1<-c(fd[1]^2,2*fd[1]*fr[1],2*fd[1]*(1-fr[1]-fd[1]),fr[1]^2,2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(fd[i]^2,2*fd[i]*fr[i],2*fd[i]*(1-fr[i]-fd[i]),fr[i]^2,2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fr[1]-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),2/4*(1-fr[1]-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4+1/4*(1-fr[1]-fd[1])+1/4*(1-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4+2/4*(1-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fr[i]-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),2/4*(1-fr[i]-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4+1/4*(1-fr[i]-fd[i])+1/4*(1-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4+2/4*(1-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp1<-t1*(1-t2)
t1<-c(2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/4+1/4*(1-fr[1]-fd[1])+1/4*(1-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4+2/4*(1-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/4+1/4*(1-fr[i]-fd[i])+1/4*(1-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4+2/4*(1-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp2<-t1*(1-t2)
gdd<-sum(temp1)-sum(temp2)
r[3,1]<-gnn+gnd*(1-e)+gdd*(1-e)^2
r[3,2]<-gnd*e+gdd*2*e*(1-e)
r[3,3]<-gdd*e^2

gnn<-prod(2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))
gnd<-prod(fd^2*1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/4*(1-fr-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(2/4*(1-fr-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-gnn
t1<-c(fd[1]^2,2*fd[1]*fr[1],2*fd[1]*(1-fr[1]-fd[1]),fr[1]^2,2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(fd[i]^2,2*fd[i]*fr[i],2*fd[i]*(1-fr[i]-fd[i]),fr[i]^2,2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fr[1]-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),2/4*(1-fr[1]-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4+1/4*(1-fr[1]-fd[1])+1/4*(1-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4+2/4*(1-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fr[i]-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),2/4*(1-fr[i]-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4+1/4*(1-fr[i]-fd[i])+1/4*(1-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4+2/4*(1-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp1<-t1*(1-t2)
t1<-c(2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/4+1/4*(1-fr[1]-fd[1])+1/4*(1-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4+2/4*(1-fd[1])+1/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/4+1/4*(1-fr[i]-fd[i])+1/4*(1-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4+2/4*(1-fd[i])+1/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp2<-t1*(1-t2)
gdd<-sum(temp1)-sum(temp2)
r[4,1]<-gnn+gnd*(1-e)+gdd*(1-e)^2
r[4,2]<-gnd*e+gdd*2*e*(1-e)
r[4,3]<-gdd*e^2

gnn<-prod(2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))
gnd<-prod(fd^2*1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/4*(1-fr-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(1/2*(1-fr-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-gnn
t1<-c(fd[1]^2,2*fd[1]*fr[1],2*fd[1]*(1-fr[1]-fd[1]),fr[1]^2,2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(fd[i]^2,2*fd[i]*fr[i],2*fd[i]*(1-fr[i]-fd[i]),fr[i]^2,2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fr[1]-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/2*(1-fr[1]-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fr[1]-fd[1])+1/4*(1-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/2*(1-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fr[i]-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/2*(1-fr[i]-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fr[i]-fd[i])+1/4*(1-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/2*(1-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp1<-t1*(1-t2)
t1<-c(2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/4*(1-fr[1]-fd[1])+1/4*(1-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/2*(1-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/4*(1-fr[i]-fd[i])+1/4*(1-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/2*(1-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp2<-t1*(1-t2)
gdd<-sum(temp1)-sum(temp2)
r[5,1]<-gnn+gnd*(1-e)+gdd*(1-e)^2
r[5,2]<-gnd*e+gdd*2*e*(1-e)
r[5,3]<-gdd*e^2

gnn<-prod(2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))
gnd<-prod(fd^2*1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/4*(1-fr-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(1/2*(1-fr-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-gnn
t1<-c(fd[1]^2,2*fd[1]*fr[1],2*fd[1]*(1-fr[1]-fd[1]),fr[1]^2,2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(fd[i]^2,2*fd[i]*fr[i],2*fd[i]*(1-fr[i]-fd[i]),fr[i]^2,2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fr[1]-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/2*(1-fr[1]-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fr[1]-fd[1])+1/4*(1-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/2*(1-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fr[i]-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/2*(1-fr[i]-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fr[i]-fd[i])+1/4*(1-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/2*(1-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp1<-t1*(1-t2)
t1<-c(2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/4*(1-fr[1]-fd[1])+1/4*(1-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/2*(1-fd[1])+1/2*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/4*(1-fr[i]-fd[i])+1/4*(1-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/2*(1-fd[i])+1/2*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp2<-t1*(1-t2)
gdd<-sum(temp1)-sum(temp2)
r[6,1]<-gnn+gnd*(1-e)+gdd*(1-e)^2
r[6,2]<-gnd*e+gdd*2*e*(1-e)
r[6,3]<-gdd*e^2

gnn<-prod(2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))
gnd<-prod(fd^2*3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/8*(1-fr-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(1/4*(1-fr-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-gnn
t1<-c(fd[1]^2,2*fd[1]*fr[1],2*fd[1]*(1-fr[1]-fd[1]),fr[1]^2,2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(fd[i]^2,2*fd[i]*fr[i],2*fd[i]*(1-fr[i]-fd[i]),fr[i]^2,2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/8*(1-fr[1]-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/8*(1-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fr[1]-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/8*(1-fr[1]-fd[1])+1/8*(1-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/8*(1-fr[i]-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/8*(1-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fr[i]-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/8*(1-fr[i]-fd[i])+1/8*(1-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp1<-t1*(1-t2)
t1<-c(2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/8*(1-fr[1]-fd[1])+1/8*(1-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/8*(1-fr[i]-fd[i])+1/8*(1-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp2<-t1*(1-t2)
gdd<-sum(temp1)-sum(temp2)
r[7,1]<-gnn+gnd*(1-e)+gdd*(1-e)^2
r[7,2]<-gnd*e+gdd*2*e*(1-e)
r[7,3]<-gdd*e^2

gnn<-prod(2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))
gnd<-prod(fd^2*3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/8*(1-fr-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(1/4*(1-fr-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-gnn
t1<-c(fd[1]^2,2*fd[1]*fr[1],2*fd[1]*(1-fr[1]-fd[1]),fr[1]^2,2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(fd[i]^2,2*fd[i]*fr[i],2*fd[i]*(1-fr[i]-fd[i]),fr[i]^2,2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/8*(1-fr[1]-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/8*(1-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fr[1]-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/8*(1-fr[1]-fd[1])+1/8*(1-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/8*(1-fr[i]-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/8*(1-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fr[i]-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/8*(1-fr[i]-fd[i])+1/8*(1-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp1<-t1*(1-t2)
t1<-c(2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/8*(1-fr[1]-fd[1])+1/8*(1-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/4*(1-fd[1])+3/4*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/8*(1-fr[i]-fd[i])+1/8*(1-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/4*(1-fd[i])+3/4*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp2<-t1*(1-t2)
gdd<-sum(temp1)-sum(temp2)
r[8,1]<-gnn+gnd*(1-e)+gdd*(1-e)^2
r[8,2]<-gnd*e+gdd*2*e*(1-e)
r[8,3]<-gdd*e^2

gnn<-prod(2*fr*(1-fr-fd)*(1/16*(1-fr-fd)+1/16*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/8*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd))))
gnd<-prod(fd^2*7/8*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/16*(1-fr-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/16*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(1/8*(1-fr-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/16*(1-fr-fd)+1/16*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/8*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/16*(1-fr-fd)+1/16*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/8*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-gnn
t1<-c(fd[1]^2,2*fd[1]*fr[1],2*fd[1]*(1-fr[1]-fd[1]),fr[1]^2,2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(fd[i]^2,2*fd[i]*fr[i],2*fd[i]*(1-fr[i]-fd[i]),fr[i]^2,2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(7/8*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/16*(1-fr[1]-fd[1])+7/8*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/16*(1-fd[1])+7/8*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/8*(1-fr[1]-fd[1])+7/8*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/16*(1-fr[1]-fd[1])+1/16*(1-fd[1])+7/8*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/8*(1-fd[1])+7/8*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(7/8*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/16*(1-fr[i]-fd[i])+7/8*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/16*(1-fd[i])+7/8*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/8*(1-fr[i]-fd[i])+7/8*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/16*(1-fr[i]-fd[i])+1/16*(1-fd[i])+7/8*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/8*(1-fd[i])+7/8*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp1<-t1*(1-t2)
t1<-c(2*fr[1]*(1-fr[1]-fd[1]),(1-fr[1]-fd[1])^2)
if(length(fd)>=2) for(i in 2:length(fd)) t1<-outer(t1,c(2*fr[i]*(1-fr[i]-fd[i]),(1-fr[i]-fd[i])^2))
t2<-c(1/16*(1-fr[1]-fd[1])+1/16*(1-fd[1])+7/8*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])),1/8*(1-fd[1])+7/8*((1-fr[1]-fd[1])^2+2*fr[1]*(1-fr[1]-fd[1])))
if(length(fd)>=2) for(i in 2:length(fd)) t2<-outer(t2,c(1/16*(1-fr[i]-fd[i])+1/16*(1-fd[i])+7/8*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i])),1/8*(1-fd[i])+7/8*((1-fr[i]-fd[i])^2+2*fr[i]*(1-fr[i]-fd[i]))))
temp2<-t1*(1-t2)
gdd<-sum(temp1)-sum(temp2)
r[9,1]<-gnn+gnd*(1-e)+gdd*(1-e)^2
r[9,2]<-gnd*e+gdd*2*e*(1-e)
r[9,3]<-gdd*e^2

list(ge=gpp*e,r=r)

}

