LVZSmp <-
function(oldmodel)
{
lvzextra=EmptyModel()

criticalz=NULL


a=AnalyzeLVZ(oldmodel$vp,oldmodel$vs,oldmodel$z,oldmodel$rp)
newvp=a[[1]]
newvpsec=a[[2]]
newzp=a[[3]]
critzp=a[[4]]

a=AnalyzeLVZ(oldmodel$vs,oldmodel$vp,oldmodel$z,oldmodel$rp)
newvs=a[[1]]
newvssec=a[[2]]
newzs=a[[3]]
critzs=a[[4]]

criticalz=c(critzp,critzs)

lvzextra=InterpModel(oldmodel,unique(c(newzp, newzs)),'simple');
newzplen=length(newzp)
for(newzpcnt in 1:newzplen){
    indy=which(lvzextra$z==newzp[newzpcnt])
    lvzextra$vp[indy]=newvp[newzpcnt]
    lvzextra$vs[indy]=newvpsec[newzpcnt]
}

newzslen=length(newzs)
for(newzscnt in 1:newzslen){
    indy=which(lvzextra$z==newzs[newzscnt])
    lvzextra$vp[indy]=newvssec[newzscnt]
    lvzextra$vs[indy]=newvs[newzscnt]
}

return(list(lvzextra=lvzextra,criticalz=criticalz))
}

