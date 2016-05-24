ALINE.segments<-function(result,Syllabic=5,Place=40,Stop=50,Voice=10,Nasal=10,Retroflex=10,Lateral=10,Aspirated=5,Long=1,High=5,Back=5,Round=5,sk=10){
  s1<-result[[3]]
  s2<-result[[4]]
  m1<-strsplit(s1,split=" ")[[1]];m1=m1[-which(m1=="")]
  m2<-strsplit(s2,split=" ")[[1]];m2=m2[-which(m2=="")]
  ind<-which(m1=="|")
  sim.vec<-vector()
  p=1
  for(i in (ind[1]+1):(ind[2]-1)){
    if(m1[i]=='-' | m1[i]=='<' | m2[i]=='-' | m2[i]=='<') sim.vec[p]=0
    else{
      sim.vec[p]=raw.alignment(c(m1[i],m2[i]),Syllabic=5,Place=40,Stop=50,Voice=10,Nasal=10,Retroflex=10,Lateral=10,Aspirated=5,Long=1,High=5,Back=5,Round=5,sk=10)[[2]]
    }
    p=p+1
  }
  if(sum(sim.vec)==result[[2]])
    return(sim.vec)
  else return("Error")
}