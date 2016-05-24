ConvVdepth2p <-
function(model,z)

{
radian=pi/180

newmodel=InterpModel(model,z,'simple')


zeroelements=which(newmodel$vp==0)
nonzeros=which(newmodel$vp!=0)
prayp=newmodel$z * 0
prayp[nonzeros]=radian*(newmodel$rp-newmodel$z[nonzeros])/newmodel$vp[nonzeros]
prayp[zeroelements]=prayp[zeroelements]+NaN


zeroelements=which(newmodel$vs==0)
nonzeros=which(newmodel$vs!=0)
srayp=0*newmodel$z
srayp[nonzeros]=radian*(newmodel$rp-newmodel$z[nonzeros])/newmodel$vs[nonzeros]
srayp[zeroelements]=srayp[zeroelements]+NaN


newz=newmodel$z

return(list(prayp=prayp,srayp=srayp,newz=newz))
}

