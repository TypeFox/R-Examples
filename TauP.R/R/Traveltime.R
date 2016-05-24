Traveltime <-
function(phase,delta,h,model,pscan=NULL)
{
if(!('criticalrays' %in% names(model))){model=ImproveModel(model)[[1]]}



a = FindP4Dist(phase,delta,h,model,pscan)

p = a[[1]]
angles = a[[2]]
dists = a[[3]]


nonnans = which(!is.na(p))
if(length(nonnans) == 0){
  tt = NaN
  dists = NaN
  angles = NaN
}else{
	p = p[nonnans] 
	angles = angles[nonnans]
	dists = dists[nonnans] 
	
	lengthp = length(p) 
	tt = p*0+NaN; 
	for(indy in 1:lengthp){
          tt[indy] = FindTime4p(phase,h,p[indy],model,takeoff=angles[indy])[[1]]
        }

}
nonnans = which(!is.na(tt))
if(length(nonnans) == 0){
  p = NaN
  tt = NaN
  dists = NaN
  angles = NaN
} else {
  tt = tt[nonnans]
  p = p[nonnans]
  dists = dists[nonnans]
  angles = angles[nonnans]
}

return(list(tt=tt,p=p,angles=angles,dists=dists))
}

