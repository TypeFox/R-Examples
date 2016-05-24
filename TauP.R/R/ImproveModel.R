ImproveModel <-
function(oldmodel){
### Define internal functions:
CheckExtreme=function(y){
eps=1.04

a = log(y[2:3]/y[1:2])
b = sign(a)*(abs(a)>log(eps)) 

yesno = b[1]!=b[2]

if(is.na(yesno)){return(FALSE)
}else{return(yesno)}
}

  
DerivSmp=function(oldmodel)
{
criticalz=NULL 
flatnewz=NULL 




vp=TransformS2Fz(oldmodel$vp,oldmodel$z,oldmodel$rp)[[1]]
a=TransformS2Fz(oldmodel$vs,oldmodel$z,oldmodel$rp) 
vs=a[[1]]
z=a[[2]]

             

dvp=diff(vp)/diff(z)
dz = diff(z)

dvs=diff(vs)/diff(z) 



dvp=abs(c(0,dvp,0))
dvs=abs(c(0,dvs,0))
dz=c(0,dz,0)
dvpdz=dvp/dz 
dvsdz=dvs/dz 
derivlen=length(dz) 
dvextrema=NULL
flatnewz=NULL

extrema=NULL
for(indy in 2:(derivlen-1) ){

    isextremum = (CheckExtreme(dvp[(-1:1)+indy])==1) |( CheckExtreme(dvs[(-1:1)+indy])==1)
    extrema=c(extrema,isextremum)
    if(isextremum){
       
       
      dvextrema=c(dvextrema,z[indy])
    } 
}

criticalz=dvextrema

flatnewz=TransformF2Sz(flatnewz,flatnewz,oldmodel$rp)[[2]]
criticalz=TransformF2Sz(criticalz,criticalz,oldmodel$rp)[[2]]

return(list(newdepths=flatnewz,criticalz=criticalz))

}

DiscontinuitySmp=function(model){

criticalz=NULL 
newdepths=NULL 


zepsilon=0.001 


z=model$z;
deltaz=diff(z)
disconindies=which(deltaz==0) 
disconcnt=length(disconindies) 
discondepths=z[disconindies] 
criticalz=c(criticalz, z[c(disconindies,disconindies+1)]) 

disconextension=c(discondepths-zepsilon, discondepths+zepsilon)

return(list(criticalz=criticalz,newdepths=disconextension))
}

### Done defining internal functions

zepsilon=0.001 

newmodel=oldmodel; 

a=DiscontinuitySmp(oldmodel)
disconnewdepths=a[[1]]
discriticalz=a[[2]]

a=LVZSmp(oldmodel)
lvzextra=a[[1]]
lvzcriticalz=a[[2]]


a=DerivSmp(oldmodel)
derivnewdepths=a[[1]]
derivcriticalz=a[[2]]

newdepths=c(disconnewdepths, derivnewdepths) 
criticaldepths=c(discriticalz, derivcriticalz, lvzcriticalz) 


alldepths=sort(c(oldmodel$z,newdepths)) 
newmodel=InterpModel(newmodel,alldepths,'preserve')


lvzanz=length(lvzextra$z) 

for(indy in 1:lvzanz){
	 samplesbefore=which(newmodel$z<lvzextra$z[indy])
	 samplesafter=(max(samplesbefore)+1):length(newmodel$z)

	 if(lvzextra$z[indy]!=newmodel$z[samplesafter[1]]){
		
		newmodel$z=c(newmodel$z[samplesbefore], lvzextra$z[indy], newmodel$z[samplesafter]) 
		newmodel$vp=c(newmodel$vp[samplesbefore],lvzextra$vp[indy],newmodel$vp[samplesafter]) 
		newmodel$vs=c(newmodel$vs[samplesbefore],lvzextra$vs[indy],newmodel$vs[samplesafter]) 
		newmodel$rho=c(newmodel$rho[samplesbefore],lvzextra$rho[indy],newmodel$rho[samplesafter]) 
		newmodel$qp=c(newmodel$qp[samplesbefore],lvzextra$qp[indy],newmodel$qp[samplesafter]) 
		newmodel$qs=c(newmodel$qs[samplesbefore],lvzextra$qs[indy],newmodel$qs[samplesafter])
		}
	}



a=ConvVdepth2p(newmodel,criticaldepths); 
prayp=a[[1]]
srayp=a[[2]]
raypz=a[[3]]

criticalrays=list()
criticalrays$z=raypz
criticalrays$p=prayp 
criticalrays$s=srayp 


if(FALSE){
crdf = data.frame(criticalrays)
keep = TRUE
for(i in 2:length(criticalrays$z)){
  keep[i] = criticalrays$z[i] != criticalrays$z[i-1]
}
criticalrays$z = round(criticalrays$z[keep], 3)
criticalrays$p = round(criticalrays$p[keep], 3)
criticalrays$s = round(criticalrays$s[keep], 3)

print(keep)
}
newmodel$criticalrays=criticalrays

return(list(newmodel=newmodel,newdepths=newdepths))
}

