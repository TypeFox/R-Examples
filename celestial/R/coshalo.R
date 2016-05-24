coshaloMvirToSigma=function(Mvir=1e12,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,H0=100){
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
rho_crit=(3*H0^2*(kms_to_ms/(1e6*pc_to_m))^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit)
return((Mvir*sqrt(32*pi*g^3*DeltaVir*rho_crit/3))^(1/3))
}

coshaloSigmaToMvir=function(Sigma=230,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,H0=100){
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
H0 = H0*kms_to_ms/(1e6*pc_to_m)
rho_crit=(3*H0^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit) #in user units
return(Sigma^3*sqrt(3/(32*pi*g^3*DeltaVir*rho_crit)))
}

coshaloMvirToRvir=function(Mvir=1e12,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,H0=100){
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
H0 = H0*kms_to_ms/(1e6*pc_to_m)
rho_crit=(3*H0^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit) #in user units
return(((3*Mvir)/(4*pi*DeltaVir*rho_crit))^(1/3))
}

coshaloRvirToMvir=function(Rvir=162.635,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,H0=100){
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
H0 = H0*kms_to_ms/(1e6*pc_to_m)
rho_crit=(3*H0^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit) #in user units
return((4*pi/3)*DeltaVir*rho_crit*Rvir^3)
}

coshaloSigmaToRvir=function(Sigma=230,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,H0=100){
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
H0 = H0*kms_to_ms/(1e6*pc_to_m)
rho_crit=(3*H0^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit) #in user units
return(Sigma*(27/(512*pi^3*g^3*DeltaVir^3*rho_crit^3))^(1/6))
}

coshaloRvirToSigma=function(Rvir=162.635,Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,H0=100){
G=6.67384e-11
msol_to_kg=1.98892e30
pc_to_m=3.08568e16
kms_to_ms=1e3
g = G*msol_to_kg/(pc_to_m*kms_to_ms^2)
g = g*Munit/(Lunit*Vunit^2)
H0 = H0*kms_to_ms/(1e6*pc_to_m)
rho_crit=(3*H0^2)/(8*pi*G) #in kg/m^3
rho_crit=rho_crit*(pc_to_m*Lunit)^3/(msol_to_kg*Munit) #in user units
return(Rvir*((512*pi^3*g^3*DeltaVir^3*rho_crit^3)/27)^(1/6))
}
