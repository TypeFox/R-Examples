# testing all the attitude representation conversion methods


Q1 <- c(-0.1677489, -0.7369231, -0.3682588, 0.5414703)
Q2 <- c(-0.8735598,  0.1145235, -0.2093062, 0.4242270)
Q3 <- c(0.426681700, -0.20287610,  0.43515810, -0.76643420)
Q4 <- matrix(c(-0.1677489, -0.7369231, -0.3682588, 0.5414703,
-0.8735598,  0.1145235, -0.2093062, 0.4242270,
0.426681700, -0.20287610,  0.43515810, -0.76643420),3,4,byrow=TRUE)

EV1 <- c(-0.1995301, -0.8765382, -0.4380279, 114.4324)
EV2 <- c(-9.646669e-001, 1.264676e-001, -2.311356e-001, 1.297965e+002)
EV3 <- c(6.642793e-001, -3.158476e-001, 6.774757e-001, 2.800695e+002)
EV4 <- matrix(c(-1.995301e-001, -8.765382e-001, -4.380280e-001, 1.144324e+002,
-9.646669e-001, 1.264676e-001, -2.311356e-001, 1.297965e+002,
6.642793e-001, -3.158476e-001, 6.774757e-001, 2.800695e+002),3,4,byrow=TRUE)

DCM1 <- matrix(c(-0.3573404, -0.1515663, 0.9215940, 0.6460385, 0.6724915, 0.3610947, -0.6744939, 0.7244189, -0.1423907),3,3,byrow=TRUE)
DCM2 <- matrix(c(0.88615060, -0.3776729,  0.2685150,-0.02249957, -0.6138316, -0.7891163,0.46285090,  0.6932344, -0.5524447),3,3,byrow=TRUE)
DCM3 <- matrix(c(0.5389574, -0.8401672,  0.06036564,0.4939131, 0.2571603, -0.83061330,0.6823304, 0.4774806, 0.55356800),3,3,byrow=TRUE)
DCM4 <- array(c(-0.35734040,  0.64603850, -0.67449390, -0.15156630,  0.67249150, 0.72441890,
0.92159400,  0.36109470, -0.14239070,  0.88615060, -0.02249957,  0.46285090,
-0.37767290, -0.61383160,  0.69323440,  0.26851500, -0.78911630, -0.55244470,
0.53895740,  0.49391310,  0.68233040, -0.84016720,  0.25716030,  0.47748060,
0.06036564, -0.83061330,  0.55356800),dim=c(3,3,3))

EAxyx <- matrix(c(0.3734309, 1.427920, 2.3205212,-1.2428130, 0.985502, 0.9821003,-1.4982479, 2.157439, 0.6105777),3,3,byrow=TRUE)
EAxyz <- matrix(c(  2.07603606, -0.7402790, -1.376712,-0.02538478,  0.4812086, -0.897943, 0.74181498,  0.7509456, -2.429857),3,3,byrow=TRUE)
EAxzy <- matrix(c( -2.9199163, -0.8101911, -1.3627437,-0.5515727, -0.7659671,  0.6973821,-1.8678238, -0.4977850,  2.2523838),3,3,byrow=TRUE)
EAyzx <- matrix(c( 1.4175036,  0.3694416,  2.3762543,0.4524244, -0.9093689, -0.0366379,3.0329736, -0.9802081,  2.0508342),3,3,byrow=TRUE)
EAyxz <- matrix(c( -2.0579912,  0.70238301,  2.6488233, 0.4813408, -0.02250147, -0.9096943, 0.9022582,  0.51658433, -1.8710396),3,3,byrow=TRUE)
EAzxy <- matrix(c( -2.3190385, -0.1521527, 1.9406908,-0.8460724, -0.3872818, 0.2942186,-2.0648275, -0.9975913, 0.1115396),3,3,byrow=TRUE)
EAzyx <- matrix(c(  1.1951869, 1.17216717, -2.7404413,-0.9600165, 0.27185113, -0.4028824,-2.1586537, 0.06040236, -1.0004280),3,3,byrow=TRUE)
EAxzx <- matrix(c( -1.197365, 1.427920, -2.391868,-2.813609, 0.985502,  2.552897,-3.069044, 2.157439,  2.181374),3,3,byrow=TRUE)
EAyxy <- matrix(c( 1.777046, 2.3083664,  0.5096786,2.069637, 0.9098912, -1.5993010,2.624796, 1.8308788, -1.0343298),3,3,byrow=TRUE)
EAyzy <- matrix(c( -2.935343, 2.3083664, -1.061118,-2.642752, 0.9098912,  3.113088,-2.087593, 1.8308788, -2.605126),3,3,byrow=TRUE)
EAzxz <- matrix(c( -0.8069433, 1.9362151, -1.733798, 1.6193689, 0.4818253, -2.523541, 0.9442344, 1.0015975, -3.069866),3,3,byrow=TRUE)
EAzyz <- matrix(c( -2.3777396, 1.9362151, -0.1630019, 0.0485726, 0.4818253, -0.9527442,-0.6265619, 1.0015975, -1.4990700),3,3,byrow=TRUE)
EAall <- rbind(EAxyx,EAxyz,EAxzy,EAyzx,EAyxz,EAzxy,EAzyx,EAxzx,EAyxy,EAyzy,EAzxz,EAzyz)
EAvct <- c('xyx','xyz','xzy','yzx','yxz','zxy','zyx','xzx','yxy','yzy','zxz','zyz')

# test Qnorm and Qnormalize
print('Qnorm and Qnormalize')
Qnormalize(Q1)
Qnormalize(Q4)
Qnorm(Q1)
Qnorm(Q4)

# test Q2EV and EV2Q
print('Q2EV')
Q2EV(Q1)
Q2EV(Q2)
Q2EV(Q3)
Q2EV(Q4)
print('EV2Q')
EV2Q(EV1,1e-7)
EV2Q(EV2,1e-7)
EV2Q(EV3,1e-7)
EV2Q(EV4,1e-7)

# test DCM2Q and Q2DCM
print('DCM2Q')
DCM2Q(DCM1)
DCM2Q(DCM2)
DCM2Q(DCM3)
DCM2Q(DCM4)
print('Q2DCM')
Q2DCM(Q1)
Q2DCM(Q2)
Q2DCM(Q3)
Q2DCM(Q4)

# test Q2EA and EA2Q
print('Q2EA')
for (EAv in EAvct) { print (Q2EA(Q1,EAv));print (Q2EA(Q2,EAv));print (Q2EA(Q3,EAv)) }
for (EAv in EAvct) print (Q2EA(Q4,EAv))
print('EA2Q')
n <- 1;for (EAv in rep(EAvct,each=3)) { print (EA2Q(EAall[n,],EAv));n <- n+1 }
n <- 1;for (EAv in EAvct) { print (EA2Q(EAall[n:(n+2),],EAv));n <- n+3 }

# test DCM2EV and EV2DCM
print('DCM2EV')
DCM2EV(DCM1)
DCM2EV(DCM2)
DCM2EV(DCM3)
DCM2EV(DCM4)
print('EV2DCM')
EV2DCM(EV1,1e-7)
EV2DCM(EV2,1e-7)
EV2DCM(EV3,1e-7)
EV2DCM(EV4,1e-7)

# test DCM2EA and EA2DCM
print('DCM2EA')
for (EAv in EAvct) { print (DCM2EA(DCM1,EAv));print (DCM2EA(DCM2,EAv));print (DCM2EA(DCM3,EAv)) }
for (EAv in EAvct) print (DCM2EA(DCM4,EAv))
print('EA2DCM')
print(EA2DCM(EAxyx,'xyx',1e-7))
print(EA2DCM(EAxyz,'xyz',1e-7))
print(EA2DCM(EAxzy,'xzy',1e-7))
print(EA2DCM(EAyzx,'yzx',1e-7))
print(EA2DCM(EAyxz,'yxz',1e-7))
print(EA2DCM(EAzxy,'zxy',1e-7))
print(EA2DCM(EAzyx,'zyx',1e-7))
print(EA2DCM(EAxzx,'xzx',1e-7))
print(EA2DCM(EAyxy,'yxy',1e-7))
print(EA2DCM(EAyzy,'yzy',1e-7))
print(EA2DCM(EAzxz,'zxz',1e-7))
print(EA2DCM(EAzyz,'zyz',1e-7))

# test EA2EV an EV2EA
print('EA2EV')
print(EA2EV(EAxyx,'xyx',1e-7))
print(EA2EV(EAxyz,'xyz',1e-7))
print(EA2EV(EAxzy,'xzy',1e-7))
print(EA2EV(EAyzx,'yzx',1e-7))
print(EA2EV(EAyxz,'yxz',1e-7))
print(EA2EV(EAzxy,'zxy',1e-7))
print(EA2EV(EAzyx,'zyx',1e-7))
print(EA2EV(EAxzx,'xzx',1e-7))
print(EA2EV(EAyxy,'yxy',1e-7))
print(EA2EV(EAyzy,'yzy',1e-7))
print(EA2EV(EAzxz,'zxz',1e-7))
print(EA2EV(EAzyz,'zyz',1e-7))
print('EV2EA')
for (EAv in EAvct) { print (EV2EA(EV1,EAv,1e-7));print (EV2EA(EV2,EAv,1e-7));print (EV2EA(EV3,EAv,1e-7)) }
for (EAv in EAvct) print (EV2EA(EV4,EAv,1e-7))

# test Q2GL
print('Q2GL')
Q2GL(Q1)
Q2GL(Q2)
Q2GL(Q3)
Q2GL(Q4)

# test EA2EA
# xyx to ...
print('EA2EA')
EA2EA(EAxyx,'xyx','xyz')
EAxyz
EA2EA(EAxyx,'xyx','xzy')
EAxzy
EA2EA(EAxyx,'xyx','yzx')
EAyzx
EA2EA(EAxyx,'xyx','yxz')
EAyxz
EA2EA(EAxyx,'xyx','zxy')
EAzxy
EA2EA(EAxyx,'xyx','zyx')
EAzyx
EA2EA(EAxyx,'xyx','xzx')
EAxzx
EA2EA(EAxyx,'xyx','yxy')
EAyxy
EA2EA(EAxyx,'xyx','yzy')
EAyzy
EA2EA(EAxyx,'xyx','zxz')
EAzxz
EA2EA(EAxyx,'xyx','zyz')
EAzyz
