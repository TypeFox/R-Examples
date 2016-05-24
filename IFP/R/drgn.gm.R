
drgn<-function(fd,fr){

pp<-1-prod((1-fd-fr)^2+2*fr*(1-fd-fr))

r<-matrix(NA,9,3) 
#1:mzt,2:parent-offspring,3:dzt,4:sibling,5:2direct,6:3rd,7:3d,8:4th,9:4d

r[1,1]<-1-pp
r[1,2]<-0
r[1,3]<-pp

r[2,1]<-prod(2*fr*(1-fr-fd)*(1/2*(1-fr-fd)+1/2*(1-fd))+(1-fr-fd)^2*(1-fd))
r[2,2]<-prod(fd^2*0+2*fd*(1-fr-fd)*1/2*(1-fd)+2*fd*fr*1/2*(1-fr-fd)+fr^2*(1-fr-fd)+2*fr*(1-fr-fd)*(1/2*(1-fr-fd)+1/2*(1-fd))+(1-fr-fd)^2*(1-fd))-prod((1-fr-fd)^2*(1-fd)+2*fr*(1-fr-fd)*(1/2*(1-fr-fd)+1/2*(1-fd)))+prod((1-fr-fd)^2+2*fr*(1-fr-fd))-r[2,1]
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
r[2,3]<-sum(temp1)-sum(temp2)

r[3,1]<-prod(2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))
r[3,2]<-prod(fd^2*1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/4*(1-fr-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(2/4*(1-fr-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-r[3,1]
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
r[3,3]<-sum(temp1)-sum(temp2)

r[4,1]<-prod(2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))
r[4,2]<-prod(fd^2*1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/4*(1-fr-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(2/4*(1-fr-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/4+1/4*(1-fr-fd)+1/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4+2/4*(1-fd)+1/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-r[4,1]
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
r[4,3]<-sum(temp1)-sum(temp2)

r[5,1]<-prod(2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))
r[5,2]<-prod(fd^2*1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/4*(1-fr-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(1/2*(1-fr-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-r[5,1]
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
r[5,3]<-sum(temp1)-sum(temp2)

r[6,1]<-prod(2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))
r[6,2]<-prod(fd^2*1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/4*(1-fr-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(1/2*(1-fr-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/4*(1-fr-fd)+1/4*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/2*(1-fd)+1/2*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-r[6,1]
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
r[6,3]<-sum(temp1)-sum(temp2)

r[7,1]<-prod(2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))
r[7,2]<-prod(fd^2*3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/8*(1-fr-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(1/4*(1-fr-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-r[7,1]
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
r[7,3]<-sum(temp1)-sum(temp2)

r[8,1]<-prod(2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))
r[8,2]<-prod(fd^2*3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/8*(1-fr-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(1/4*(1-fr-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/8*(1-fr-fd)+1/8*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/4*(1-fd)+3/4*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-r[8,1]
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
r[8,3]<-sum(temp1)-sum(temp2)

r[9,1]<-prod(2*fr*(1-fr-fd)*(1/16*(1-fr-fd)+1/16*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/8*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd))))
r[9,2]<-prod(fd^2*7/8*((1-fr-fd)^2+2*fr*(1-fr-fd))+2*fd*fr*(1/16*(1-fr-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fd*(1-fr-fd)*(1/16*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+fr^2*(1/8*(1-fr-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+2*fr*(1-fr-fd)*(1/16*(1-fr-fd)+1/16*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/8*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd))))-prod(2*fr*(1-fr-fd)*(1/16*(1-fr-fd)+1/16*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd)))+(1-fr-fd)^2*(1/8*(1-fd)+7/8*((1-fr-fd)^2+2*fr*(1-fr-fd))))+prod(2*fr*(1-fr-fd)+(1-fr-fd)^2)-r[9,1]
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
r[9,3]<-sum(temp1)-sum(temp2)

list(pp=pp,r=r)

}

