snapwrite <-
function(part,head,file){
data = file(file,'wb')
#first header block
writeBin(as.integer(256),data)
writeBin(as.integer(head$Npart),data)
writeBin(as.numeric(head$Massarr),data,size=8)
writeBin(as.numeric(head$Time),data,size=8)
writeBin(as.numeric(head$z),data,size=8)
writeBin(as.integer(head$FlagSfr),data)
writeBin(as.integer(head$FlagFeedback),data)
writeBin(as.integer(head$Nall),data)
writeBin(as.integer(head$FlagCooling),data)
writeBin(as.integer(head$NumFiles),data)
writeBin(as.numeric(head$BoxSize),data,size=8)
writeBin(as.numeric(head$OmegaM),data,size=8)
writeBin(as.numeric(head$OmegaL),data,size=8)
writeBin(as.numeric(head$h),data,size=8)
writeBin(as.integer(head$FlagAge),data)
writeBin(as.integer(head$FlagMetals),data)
writeBin(as.integer(head$NallHW),data)
writeBin(as.integer(head$flag_entr_ics),data)
writeBin(as.integer(rep(0,length=as.integer(256-241))),data)
#last head block
writeBin(as.integer(256),data)

#1 data block = Positions
writeBin(as.integer(sum(head$Npart)*3*4),data)
posall=writeBin(as.numeric(t(part[,c('x','y','z')])),data,size=4)
writeBin(as.integer(sum(head$Npart)*3*4),data)
#2 data block = Velocities
writeBin(as.integer(sum(head$Npart)*3*4),data)
posall=writeBin(as.numeric(t(part[,c('vx','vy','vz')])),data,size=4)
writeBin(as.integer(sum(head$Npart)*3*4),data)
#3 data block = IDs
writeBin(as.integer(sum(head$Npart)*4),data)
writeBin(as.integer(part[,'ID']),data)
writeBin(as.integer(sum(head$Npart)*4),data)
#4 data block = Masses
if('Mass' %in% colnames(part)){
writeBin(as.integer(sum(head$Npart)*4),data)
writeBin(as.numeric(part[,'Mass']),data,size=4)
writeBin(as.integer(sum(head$Npart)*4),data)
}
close(data)
}
