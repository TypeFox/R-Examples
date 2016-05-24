source("../R/longData.r")
source("../R/longDataPlot.r")

cat("\n####################################################################
########################## Test  LongData ##########################
############################# Creation #############################
####################################################################\n")

cleanProg(.LongData.validity,,,0)
cleanProg(longData)
cleanProg(as.longData)
### Constructeurs
new("LongData")
#new("LongData",traj=array(c(1,2,3,1,4,6,1,8,10),dim=c(3,3)))
#new("LongData",traj=array(c(1,2,3,1,4,6,1,8,10),dim=c(3,3)),time=c(2,4,8),varName="T")

tr1 <- tr1n <- array(c(1,2,3,1,4, 3,6,1,8,10, 1,2,1,3,2, 4,2,5,6,3, 4,3,4,4,4, 7,6,5,5,4),
            dim=c(3,5,2),
            dimnames=list(c(101,102,104),c("T1","T2","T4","T8","T16"),c("P","A"))
            )
tr1n[1,2,1] <- NA; tr1n[2,4,2] <- NA; tr1n[3,1,2] <- NA; tr1n[3,3,2] <- NA;

new("LongData",
    traj=tr1,
    idAll=as.character(c(100,102)),
    idFewNA=as.character(c(101,102,104)),
    time=c(1,2,4,8,16),
    varNames=c("P","A"),
    maxNA=3
    )


tr2 <- array(c(1,2,3, 1,4,3, 6,1,8, 10,1,2,
              6,1,8, 10,1,2, 1,3,2, 4,2,5,
              1,3,2, 4,2,5, 6,3,4, 3,4,4,
              4,7,6, 5,5,4,  4,7,6, 5,5,4),
            dim=c(4,3,4),
            dimnames=list(c("I1","I2","I3","I4"),c("T1","T2","T4"),c("P","A","E","R"))
            )
new("LongData",
    traj=tr2,
    idFewNA=c("I1","I2","I3","I4"),
    idAll=c("I1","I2","I3","I4"),
    time=c(1,2,4),
    varNames=c("P","A","E","R"),
    maxNA=2
    )


tr3n <- array(c(1,NA,NA, 1,4,3,
              NA,1,8, 10,NA,2,
              4,NA,6, NA,5,4),
            dim=c(3,3,2),
            dimnames=list(c("I1","I2","I3"),c("T1","T2","T4"),c("P","A"))
            )
new("LongData",
    traj=tr3n,
    idAll=c("I1","I2","I3"),
    idFewNA=c("I1","I2","I3"),
    time=c(1,2,4),
    varNames=c("P","A"),
    maxNA=2
    )


tr4n <- array(c(NA,NA,2,
               NA,NA,1,
               NA,NA,NA,
               NA,3,2,
               2,NA,1,
               1,2,1,

               3,NA,NA,
               4,NA,6,
               5,2,4,
               2,NA,NA,
               NA,4,2,
               1,NA,2),
            dim=c(3,6,2),
            dimnames=list(c("T1","T2","T4"),c(101,102,103,105,106,107),c("P","A"))
            )
tr4n <- aperm(tr4n,c(2,1,3))

new("LongData",
    traj=tr4n,
    idAll=as.character(c(101,102,103,105,106,107)),
    idFewNA=as.character(c(101,102,103,105,106,107)),
    time=c(1,2,4),
    varNames=c("P","A"),
    maxNA=c(1,1)
    )


cat("\n###################################################################
########################## Test  LongData #########################
########################### Constructeur ##########################
###################################################################\n")

cleanProg(.LongData.constructor,,,1) #all

longData()
#longData(traj=array(c(1,2,3,1,4,6,1,8,10),dim=c(3,3)))

longData(traj=tr1,idAll=as.character(c(101,102,104)),time=c(1,2,4,8,16),varNames=c("P","A"),maxNA=3)
longData(traj=tr2,idAll=as.character(c(1,2,3,4)),time=c(1,2,4),varNames=c("P","A","E","R"),maxNA=2)
longData(traj=tr3n,idAll=as.character(c(1,2,3)),time=c(1,2,4),varNames=c("P","A"),maxNA=2)
longData(traj=tr4n,idAll=c(1,2,3,4,5,6)+100,time=c(1,2,4),varNames=c("P","A"),maxNA=2)

### Vérification de l'exclusion des manquantes
longData(traj=tr4n,idAll=c(1,2,3,4,5,6)+100,time=c(1,2,4),varNames=c("P","A"),maxNA=1)
longData(traj=tr4n,idAll=c(1,2,3,4,5,6)+100,time=c(1,2,4),varNames=c("P","A"),maxNA=2)
longData(traj=tr4n,idAll=c(1,2,3,4,5,6)+100,time=c(1,2,4),varNames=c("P","A"),maxNA=c(1,1))
longData(traj=tr4n,idAll=c(1,2,3,4,5,6)+100,time=c(1,2,4),varNames=c("P","A"),maxNA=c(2,2))
longData(traj=tr4n,idAll=c(1,2,3,4,5,6)+100,time=c(1,2,4),varNames=c("P","A"),maxNA=c(2,1))


### Base de données
#cleanProg(as.longData.data.frame,,,0)
#cleanProg(as.longData.array,,,0)

ld0 <- longData()
ld1 <- longData(traj=tr1,idAll=c(101,102,104),time=c(1,2,4,8,16),varNames=c("Pa","Av"),maxNA=3)
ld1n <- longData(traj=tr1n,idAll=c(101,102,104),time=c(1,2,4,8,16),varNames=c("Pa","Av"),maxNA=3)

data <- read.csv2("example.csv")
ld2 <- as.longData(data,timeDataFrame=list(A=c(2,4),P=c(5,7)),time=2:3,varNames=c("Av","Pe"))
ld2 <- as.longData(data,timeDataFrame=list(A=c(2,4),P=c(5,7)),time=2:3)
try(ld2 <- as.longData(data,timeDataFrame=list(c(2,4),c(5,7)),time=2:3))
ld2 <- as.longData(data,timeDataFrame=list(A=c(2,4),P=c(5,7)),time=2:3,varNames=c("Av","Pe"))
ld2n <- as.longData(data,timeDataFrame=list(A=c(2,NA,4),P=c(NA,5,7)),time=2:4)
ld2 <- as.longData(data,timeDataFrame=list(V21=c(2,3,4),V4=c(5,6,7)),time=c(11,13,14))

dn3 <- read.csv("DatasetKML.csv")[1:200,,]
ld3n <- as.longData(dn3,time=1:6,timeDataFrame=list(cred=3:8,creq=9:14,croq=c(24:28,NA)))
ld3 <- ld3n
imputation(ld3)
#for(i in 1:2400){dn3[floor(runif(1,1,244)),floor(runif(1,2,29))]<-NA}
#ld3n <- as.longData(dn3,timeCol=list(T=2:28),timeReal=list(0:26))

ld4 <- as.longData(data=array(rnorm(30*5*12),dim=c(30,5,12)))


time=c(1,2,3,4,8,12,16,20)
id1=1:12
id2=1:200
id3=1:1200
#varNames=c("A","B")
f <- function(id,t)((id-1)%%3-1) * t
g <- function(id,t)(id%%2+1)*t
tra1 <- array(cbind(outer(id1,time,f),outer(id1,time,g))+rnorm(12*8*2,0,2),dim=c(12,8,2))
ld6 <- as.longData(tra1)
tra2 <- array(cbind(outer(id2,time,f),outer(id2,time,g))+rnorm(200*8*2,0,3),dim=c(200,8,2))
ld7 <- as.longData(tra2)
tra3 <- array(cbind(outer(id3,time,f),outer(id3,time,g))+rnorm(1200*8*2,0,4),dim=c(1200,8,2))
ld8 <- as.longData(tra3)

for(i in 1:64){tra1[floor(runif(1,1,12)),floor(runif(1,1,9)),floor(runif(1,1,3))]<-NA}
ld6n <- as.longData(tra1,maxNA=c(2,6))
for(i in 1:640){tra2[floor(runif(1,1,200)),floor(runif(1,1,9)),floor(runif(1,1,3))]<-NA}
ld7n <- as.longData(tra2,maxNA=3)
for(i in 1:6400){tra3[floor(runif(1,1,1200)),floor(runif(1,1,9)),floor(runif(1,1,3))]<-NA}
ld8n <- as.longData(tra3,maxNA=4)
