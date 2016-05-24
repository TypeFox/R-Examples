addhead <-
function(part, Npart=2, Massarr=0, Time=0, z=0, FlagSfr=0, FlagFeedback=0, FlagCooling=0, BoxSize=0, OmegaM=0, OmegaL=0, h=1, FlagAge=0, FlagMetals=0, NallHW=0, flag_entr_ics=0){
N=length(part[,1])

NpartVec=rep(0,6)
NpartVec[Npart]=N

MassarrVec=rep(0,6)
MassarrVec[Npart]=Massarr

NallVec=rep(0,6)
NallVec[Npart]=N

NallHWVec=rep(0,6)
NallHWVec[Npart]=NallHW

return=list(part=part,head=list(Npart=NpartVec, Massarr=MassarrVec, Time=Time, z=z, FlagSfr=FlagSfr, FlagFeedback=FlagFeedback, Nall=NallVec, FlagCooling=FlagCooling, NumFiles=1, BoxSize=BoxSize, OmegaM=OmegaM, OmegaL=OmegaL, h=h, FlagAge=FlagAge, FlagMetals=FlagMetals, FlagMetals=FlagMetals, NallHW=NallHWVec, flag_entr_ics=flag_entr_ics))
}
