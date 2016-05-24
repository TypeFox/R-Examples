q012=function(p=c(0.5,0.5)){
ped = nuclearPed(1, sex = 2)
lp=length(p)
m=marker(ped,alleles=1:lp,afreq=p)
q0=oneMarkerDistribution(ped,c(1,2),partialmarker=m,verbose=FALSE)
q1=oneMarkerDistribution(ped,c(1,3),partialmarker=m,verbose=FALSE)
q2=oneMarkerDistribution(ped,1,partialmarker=m,,verbose=FALSE)
list(q0=q0,q1=q1,q2=q2)
}
