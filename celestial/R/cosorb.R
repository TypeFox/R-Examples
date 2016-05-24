cosorbVisViva=function(M=1e12,Rad=162.635,SemiMajRad=162.635,Munit=1,Lunit=1e3,Vunit=1){
  G=6.67384e-11
  msol_to_kg=1.98892e30
  pc_to_m=3.08568e16
  kms_to_ms=1e3
  g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
  g = g*Munit/(Lunit*Vunit^2)
  return(sqrt(g*M*(2/Rad-1/SemiMajRad)))
}

cosorbFreeFall=function(M1=1e12,M2=1,Rad=162.635,Munit=1,Lunit=1e3,Vunit=1,Tunit=1e9){
  G=6.67384e-11
  msol_to_kg=1.98892e30
  pc_to_m=3.08568e16
  kms_to_ms=1e3
  g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
  g = g*Munit/(Lunit*Vunit^2)
  time=(pi/sqrt(g*(M1+M2)))*(Rad/2)^(3/2)*(Lunit/1e3)
  time=time*1e9/(Tunit*0.9778139)
  return(time)
}

cosorbRocheRad=function(M1=1e12,M2=1e10,Size=35.03865,Rfac=2.44){
  return(Rfac*Size*(M1/M2)^(1/3))
}

cosorbRocheSize=function(M1=1e12,M2=1e10,Rad=396.8294,Rfac=2.44){
  return((Rad/Rfac)*(M2/M1)^(1/3))
}
