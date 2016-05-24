### R code from vignette source 'formatting.data.Rnw'

###################################################
### code chunk number 1: formatting.data.Rnw:20-24
###################################################
library (Frames2)
data(Dat)

head(Dat, 3)


###################################################
### code chunk number 2: formatting.data.Rnw:31-39
###################################################
attach(Dat)

DomainOnlyLandline <- Dat[Landline == 1 & Cell == 0,]
DomainBothLandline <- Dat[Drawnby == 1 & Landline == 1 & 
                            Cell == 1,]
DomainOnlyCell <- Dat[Landline == 0 & Cell == 1,]
DomainBothCell <- Dat[Drawnby == 2 & Landline == 1 & 
                        Cell == 1,]


###################################################
### code chunk number 3: formatting.data.Rnw:44-46
###################################################
FrameLandline <- rbind(DomainOnlyLandline, DomainBothLandline)
FrameCell <- rbind(DomainOnlyCell, DomainBothCell)


###################################################
### code chunk number 4: formatting.data.Rnw:51-58
###################################################
Domain <- c(rep("a", nrow(DomainOnlyLandline)), rep("ab", 
            nrow(DomainBothLandline)))
FrameLandline <- cbind(FrameLandline, Domain)

Domain <- c(rep("b", nrow(DomainOnlyCell)), rep("ba", 
            nrow(DomainBothCell)))
FrameCell <- cbind(FrameCell, Domain)


###################################################
### code chunk number 5: formatting.data.Rnw:63-66
###################################################
Hartley(FrameLandline$Opinion, FrameCell$Opinion, 
        FrameLandline$ProbLandline, FrameCell$ProbCell, 
        FrameLandline$Domain, FrameCell$Domain)


