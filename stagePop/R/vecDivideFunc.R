vecDivideFunc=function(definePop){

  ns=definePop$numSpecies
  nst=definePop$numStages
  nstr=definePop$numStrains
  dependD=definePop$timeDependLoss
  dependT=definePop$timeDependDuration
  
  numNs=0*seq(1,ns);  numPs=0*seq(1,ns);  numTs=0*seq(1,ns)
  for (species in seq(1,ns)){
    numNs[species]=nst[species]*nstr[species]
    if (dependD[species] | dependT[species]){
      numPs[species]=(nst[species]-1)*nstr[species]
    }else{numPs[species]=0}
    if (dependT[species]){numTs[species]=(nst[species]-1)*nstr[species]}else{numTs[species]=0}
  }
  
  Ndivs=matrix(0,nrow=ns,ncol=2)
  Ndivs[1,]=c(1,nst[1]*nstr[1])
  if (ns>1){for (species in seq(2,ns)){Ndivs[species,]=c(Ndivs[species-1,2]+1,Ndivs[species-1,2]+nst[species]*nstr[species])}}
  
  Pdivs=matrix(0,nrow=ns,ncol=2)
  if (numPs[1]>0){
    Pdivs[1,]=c(sum(numNs)+1,sum(numNs)+numPs[1])
    ct=1
  }else{Pdivs[1,]=c(0,0); ct=0}
  if (ns>1){
    for (species in seq(2,ns)){
      if (numPs[species]>0){
        if (ct==0){
          Pdivs[species,]=c(sum(numNs)+1,sum(numNs)+numPs[species])
          ct=species
        }else{
          ct=ct+1
          Pdivs[species,]=c(Pdivs[ct-1,2]+1,Pdivs[ct-1,2]+numPs[species])}}}}
  
  Tdivs=matrix(0,nrow=ns,ncol=2)
  if (numTs[1]>0){
    Tdivs[1,]=c(sum(numNs)+sum(numPs)+1,sum(numNs)+sum(numPs)+numTs[1])
    ct=1
  }else{Tdivs[1,]=c(0,0); ct=0}
  if (ns>1){
    for (species in seq(2,ns)){
      if (numTs[species]>0){
        if (ct==0){
          Tdivs[species,]=c(sum(numNs)+sum(numPs)+1,sum(numNs)+sum(numPs)+numTs[species])
          ct=species
        }else{
          ct=ct+1
          Tdivs[species,]=c(Tdivs[ct-1,2]+1,Tdivs[ct-1,2]+numTs[species])}}}
  }
  mat=cbind(Ndivs,Pdivs,Tdivs)
#  print(Ndivs)
  return(mat)
}
