snapread <-
function(file){
data = file(file,'rb')
#first header block
block=readBin(data,'integer',n=1)
Npart=readBin(data,'integer',n=6)
Massarr=readBin(data,'numeric',n=6,size=8)
Time=readBin(data,'numeric',n=1,size=8)
z=readBin(data,'numeric',n=1,size=8)
FlagSfr=readBin(data,'integer',n=1)
FlagFeedback=readBin(data,'integer',n=1)
Nall=readBin(data,'integer',n=6)
FlagCooling=readBin(data,'integer',n=1)
NumFiles=readBin(data,'integer',n=1)
BoxSize=readBin(data,'numeric',n=1,size=8)
OmegaM=readBin(data,'numeric',n=1,size=8)
OmegaL=readBin(data,'numeric',n=1,size=8)
h=readBin(data,'numeric',n=1,size=8)
FlagAge=readBin(data,'integer',n=1)
FlagMetals=readBin(data,'integer',n=1)
NallHW=readBin(data,'integer',n=6)
flag_entr_ics=readBin(data,'integer',n=1)
readBin(data,'integer',n=256-241)
#last head block
block=readBin(data,'integer',n=1)

#1 data block = Positions
block=readBin(data,'integer',n=1)
posall=readBin(data,'numeric',n=block/4,size=4)
block=readBin(data,'integer',n=1)
#2 data block = Velocities
block=readBin(data,'integer',n=1)
velall=readBin(data,'numeric',n=block/4,size=4)
block=readBin(data,'integer',n=1)
#3 data block = IDs
block=readBin(data,'integer',n=1)
ID=readBin(data,'integer',n=block/4)
block=readBin(data,'integer',n=1)
#4 data block = Masses
block=readBin(data,'integer',n=1)
if(length(block)>0){
    Mass=readBin(data,'numeric',n=block/4,size=4)
}else{
    counter=1
    Mass=rep(NA,sum(Npart))
    whichmass=which(Npart>0)
        for(i in 1:length(whichmass)){
            N=Npart[whichmass[i]]
            Mass[ID>=counter & ID<=counter+N]=Massarr[whichmass[i]]
            counter=counter+N
        }
}
block=readBin(data,'integer',n=1)
#Extra blocks
extra=0
extramat={}
while(length(block)>0){
block=readBin(data,'integer',n=1)
	if(length(block)>0){
		extramat=cbind(extramat,readBin(data,'numeric',n=block/4,size=4))
		block=readBin(data,'integer',n=1)
		extra=extra+1
	}
}

close(data)

extract=((1:sum(Npart))*3)-2
part=data.frame(ID=ID,x=posall[extract],y=posall[extract+1],z=posall[extract+2],vx=velall[extract],vy=velall[extract+1],vz=velall[extract+2],Mass=Mass)

return(list(part=part,head=list(Npart = Npart, Massarr= Massarr, Time= Time, z= z, FlagSfr= FlagSfr, FlagFeedback= FlagFeedback, Nall= Nall, FlagCooling= FlagCooling, NumFiles= NumFiles, BoxSize= BoxSize, OmegaM= OmegaM, OmegaL= OmegaL,h=h, FlagAge= FlagAge, FlagMetals= FlagMetals, NallHW= NallHW,flag_entr_ics=flag_entr_ics),extra=extra,extramat=extramat))}
