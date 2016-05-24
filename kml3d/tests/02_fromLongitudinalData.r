cat("### Jeux de données de longitudinalData ###\n")


tr1 <- tr1n <- array(c(1,2,3,1,4, 3,6,1,8,10, 1,2,1,3,2, 4,2,5,6,3, 4,3,4,4,4, 7,6,5,5,4),
            dim=c(3,5,2),
            dimnames=list(c(101,102,104),c("t1","t2","t4","t8","t16"),c("P","A"))
            )
tr1n[1,2,1] <- NA; tr1n[2,4,2] <- NA; tr1n[3,1,2] <- NA; tr1n[3,3,2] <- NA;

tr2 <- array(c(1,2,3, 1,4,3, 6,1,8, 10,1,2,
              6,1,8, 10,1,2, 1,3,2, 4,2,5,
              1,3,2, 4,2,5, 6,3,4, 3,4,4,
              4,7,6, 5,5,4,  4,7,6, 5,5,4),
            dim=c(4,3,4),
            dimnames=list(c("i1","i2","i3","i4"),c("t1","t2","t4"),c("P","A","E","R"))
            )

tr3n <- array(c(1,NA,NA, 1,4,3,
                NA,1,8, 10,NA,2,
                4,NA,6, NA,5,4,
                3, 6,5, 5,5,6,
                NA,4,5, 5,NA,6,
                3,NA,5, 7,8,7,
                4,6,7,  4,NA,5,
                4,6,5,  NA,6,7
                ),
            dim=c(8,3,2),
            dimnames=list(c("i1","i2","i3","i4","i5","i6","i7","i8"),c("t1","t2","t4"),c("P","A"))
            )

tr4 <- tr4n <- array(rnorm(200*6*2),
            dim=c(200,6,2),
            dimnames=list(1:200,c("t1","t2","t4","t5","t8","t10"),c("P","A"))
            )
for(i in 1:640){tr4n[floor(runif(1,1,200)),floor(runif(1,1,7)),floor(runif(1,1,3))]<-NA}


time=c(1,2,3,4,8,12,16,20)
id1=1:8
id2=1:200
id3=1:2000
## #varNames=c("A","B")
f <- function(id,t)((id-1)%%3-1) * t
g <- function(id,t)(id%%2+1)*t

tra6 <- tra7n <- array(cbind(outer(id1,time,f),outer(id1,time,g))+rnorm(8*8*2,0,2),dim=c(8,8,2))
tra7 <- tra7n <- array(cbind(outer(id2,time,f),outer(id2,time,g))+rnorm(200*8*2,0,3),dim=c(200,8,2))

tra5n <- tra5 <- array(cbind(outer(id3,time,f),outer(id3,time,g))+rnorm(2000*8*2,0,4),dim=c(2000,8,2))
for(i in 1:6400){tra5n[floor(runif(1,1,1200)),floor(runif(1,1,9)),floor(runif(1,1,3))]<-NA}


dn6 <- read.csv("DatasetKML.csv")[1:200,,]

## LD6n <- longData3d(dn6,time=1:6,timeInData=list(cred=3:8,creq=9:14,croq=c(24:28,NA)),maxNA=5)
## traj6 <- LD6n["traj"]
## traj6[is.na(traj6)] <- 0
## LD6 <- longData3d(traj6)

## tra7 <- tra7n<- array(rnorm(8*5*12),dim=c(8,5,12))
## LD7 <- longData3d(tra7)
## for(i in 1:150){tra7n[floor(runif(1,1,9)),floor(runif(1,1,6)),floor(runif(1,1,13))]<-NA}
## LD7n <- longData3d(tra7n)





ld0 <- longData()

ld1 <- longData(traj=array(c(1,2,3,1,4,6,1,8,10),dim=c(3,3)),idAll=c(11,12,13),time=c(2,4,8),varNames="T")

dn2 <- data.frame(idAll=c(10,17,28,29,31),t1=c(5,4,2,1,0),t2=c(5,4,3,2,1),t4=c(4,2,3,1,1),t5=c(5,6,2,1,0))
ld2 <- longData(dn2[,-1],idAll=dn2[,1],time=c(1,2,5,6))

data3Imp <- rbind(c(1,2,3,4),
                  c(1,3,1,1),
                  c(2,3,4,5),
                  c(2,2,2,12),
                  c(3,4,4,6),
                  c(3,3,3,3),
                  c(2,4,8,5),
                  c(2,3,4,1))
dim(data3Imp) <- c(8,4)
ld3 <- longData(data3Imp)


dn4 <- dn4n <- read.csv2("../../01_longitudinalDataNewPlot/testsDev/divergingLines.csv")[c(1:180,101:120),]
ld4 <- longData(dn4[,-1],idAll=(1:200)*2,time=c(0,1,2,3,4,6,8,10,12,16,20))

dn5 <- read.csv2("../../../02_Recherche/Anorexie Tamara/trajectoires de soins.csv")[rep(1:200,10),-29]
ld5 <- longData(dn5[,-1],idAll=1:2000)




p0a <- p0b <- partition()

p1a <- partition(clusters=c("A","B","B"))
p1b <- partition(clusters=c("A","B","A"),ld1)
p1b <- partition(clusters=c("A","B","A"),ld1)
p1c <- partition(clusters=c("A","C","B"),ld1)
p1d <- partition(clusters=c("A","C","A"),ld1)
p1d <- partition(clusters=c("A","C","A"),ld1)

p2a <- partition(clusters=c("A","A","B","A","B"))
(p2b <- partition(clusters=c("A","B","A","A","B"),ld2,details=c(convergenceTime="3",algorithm="kml",aze=4,multiplicity="4")))
p2c <- partition(clusters=c("A","C","B","A","B"),ld2)

p2an <- partition(clusters=c("A",NA,"B",NA,"B"))
p2bn <- partition(clusters=c("A",NA,NA,"B",NA))
tryBug(p2cn <- partition(clusters=c(NA,NA,NA,NA,NA)))

p3a <- partition(clusters=rep(LETTERS[1:2],4))
p3b <- partition(clusters=rep(LETTERS[c(1:3,1:3,1:2)]),ld3,details=c(multiplicity="1"))
p3c <- partition(clusters=rep(LETTERS[1:4],each=2),ld3)
p3d <- partition(clusters=rep(c(1,2,1,3,2,3,3,2)),ld3)
p3e <- partition(clusters=rep(c(4,4,1,3,2,3,3,2)),ld3)
p3f <- partition(clusters=rep(c(1,1,1,3,3,3,3,2)),ld3)
p3g <- partition(clusters=rep(c(6,5,4,3,2,3,1,2)),ld3)
p3h <- partition(clusters=rep(c(2,2,1,2,2,3,3,2)),ld3)
p3i <- partition(clusters=rep(c(1,2,3,3,3,3,3,2)),ld3)
p3j <- partition(clusters=rep(c(1,2,1,3,4,3,4,2)),ld3)


p4a <- partition(clusters=c(rep(1:2,each=100)))
p4b <- partition(clusters=c(rep(1:3,c(50,30,120))),ld4)
tryBug(p4b["clusters"][1:6]<-"B")
tryBug(p4b["clusters"][7:9]<-"C")
p4c <- partition(clusters=c(rep(c(1:4,2:3,2:3),25)),ld4)
p4d <- partition(clusters=c(rep(c(3,2:4,1,3:5),25)),ld4)
p4e <- partition(clusters=c(rep(c(1:6,2,2),25)),ld4)
p4f <- partition(clusters=c(rep(c(1:7,4),25)),ld4)

p5a <- partition(clusters=LETTERS[rep(c(1,2),1000)])
p5b <- partition(clusters=LETTERS[c(rep(c(1,2,3),666),1:2)],ld5)
p5c <- partition(clusters=rep(1:4,500),ld5)
p7e <- partition(clusters=c(rep(1:6,333),1,2),ld5)
p7g <- partition(clusters=c(rep(1:8,250)),ld5)




cat("--- Fin jeux de données de longitudinalData ---\n")
