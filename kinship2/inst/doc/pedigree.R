### R code from vignette source 'pedigree.snw'

###################################################
### code chunk number 1: desc
###################################################
options(width = 90)
desc <- packageDescription("kinship2")


###################################################
### code chunk number 2: pedigree.snw:45-46 (eval = FALSE)
###################################################
## library(kinship2)


###################################################
### code chunk number 3: pedigree.snw:49-50
###################################################
require(kinship2)


###################################################
### code chunk number 4: pedList
###################################################

data(sample.ped)
sample.ped[1:10,]

pedAll <- pedigree(id=sample.ped$id, 
                dadid=sample.ped$father, momid=sample.ped$mother, 
                sex=sample.ped$sex, famid=sample.ped$ped)

print(pedAll)



###################################################
### code chunk number 5: pedigree.snw:95-104
###################################################

ped1basic <- pedAll['1']
ped2basic <- pedAll['2']

print(ped1basic)
print(ped2basic)

plot(ped2basic)
# plot(ped1basic)


###################################################
### code chunk number 6: kinship
###################################################

kin2 <- kinship(ped2basic)

kin2



###################################################
### code chunk number 7: kinAll
###################################################

pedAll <- pedigree(id=sample.ped$id, 
                dadid=sample.ped$father, momid=sample.ped$mother, 
                sex=sample.ped$sex, famid=sample.ped$ped)

kinAll <- kinship(pedAll)

kinAll[1:14,1:14]

kinAll[40:43,40:43]

kinAll[42:46, 42:46]



###################################################
### code chunk number 8: censor
###################################################

df2 <- sample.ped[sample.ped$ped==2,]
names(df2)

df2$censor <- c(1,1, rep(0, 12))

ped2 <- pedigree(df2$id, df2$father, df2$mother, 
                 df2$sex, status=df2$censor)



###################################################
### code chunk number 9: affected
###################################################

ped2 <- pedigree(df2$id, df2$father, df2$mother, 
                 df2$sex, affected=df2$affected,
                 status=df2$censor)

aff2 <- data.frame(blue=df2$affected, 
                   bald=c(0,0,0,0,1,0,0,0,0,1,1,0,0,1))

ped2 <- pedigree(df2$id, df2$father, df2$mother, 
                 df2$sex, affected=as.matrix(aff2),
                 status=df2$censor)



###################################################
### code chunk number 10: twins
###################################################

## create twin relationships
relate2 <- matrix(c(210,211,1,
                   212,213,3), nrow=2, byrow=TRUE)

ped2 <- pedigree(df2$id, df2$father, df2$mother, 
                 df2$sex, affected=as.matrix(aff2),
                 status=df2$censor,
                 relation=relate2)


###################################################
### code chunk number 11: ped2update
###################################################
id2 <- paste(df2$id, c("John", "Linda", "Jack", "Rachel", "Joe", "Deb", 
                         "Lucy", "Ken", "Barb", "Mike", "Matt", 
                         "Mindy", "Mark", "George"), sep="\n")

plot(ped2, col=ifelse(df2$avail, 2, 1),
     id=id2)

pedigree.legend(ped2, location="topright", radius=.3) 


###################################################
### code chunk number 12: pedigree.snw:298-308
###################################################
df1<- sample.ped[sample.ped$ped==1,]
relate1 <- matrix(c(113, 114, 4), nrow=1)

ped1 <- pedigree(df1$id, df1$father, df1$mother, 
       df1$sex, affected=df1$affected, 
                 relation=relate1)

print(ped1)
plot(ped1, col=df1$avail+1)



###################################################
### code chunk number 13: ordering
###################################################

df1reord <- df1[c(35:41,1:34),]
ped1reord <- pedigree(df1reord$id, df1reord$father, df1reord$mother, 
       df1reord$sex, affected=df1reord$affected, relation=relate1)

plot(ped1reord, col=df1reord$avail+1)



###################################################
### code chunk number 14: ped2df
###################################################


dfped2 <- as.data.frame(ped2)
dfped2



###################################################
### code chunk number 15: subset
###################################################

ped2.rm210 <- ped2[-10]

data.frame(ped2.rm210)

ped2.rm210$relation

ped2$relation




###################################################
### code chunk number 16: trim
###################################################

ped2.trim210 <- pedigree.trim(210, ped2)

ped2.trim210$id
ped2.trim210$relation

ped2.trim.more <- pedigree.trim(c(212,214), ped2.trim210)
ped2.trim.more$id

ped2.trim.more$relation



###################################################
### code chunk number 17: pedigree.snw:426-435
###################################################

shrink1.B30 <- pedigree.shrink(ped=ped1,
                 avail=df1$avail, maxBits=30)

print(shrink1.B30)
ped1.B30 <- shrink1.B30$pedObj
plot(ped1.B30, col=shrink1.B30$avail + 1)




###################################################
### code chunk number 18: pedigree.snw:450-459
###################################################

set.seed(10)
shrink1.B25 <- pedigree.shrink(ped=ped1, avail=df1$avail, 
                               maxBits=25)
print(shrink1.B25)
ped1.B25 <- shrink1.B25$pedObj

plot(ped1.B25, col=shrink1.B25$avail + 1)



###################################################
### code chunk number 19: unrelated
###################################################

df2<- sample.ped[sample.ped$ped==2,]

ped2 <- pedigree(df2$id, df2$father, df2$mother, 
       df2$sex, affected=df2$affected)

set.seed(10)
set1 <- pedigree.unrelated(ped2, avail=df2$avail)
set1
set2 <- pedigree.unrelated(ped2, avail=df2$avail)
set2



###################################################
### code chunk number 20: unrelVerify
###################################################

df2
kin2[df2$avail==1,df2$avail==1]



###################################################
### code chunk number 21: pedigree.snw:514-515
###################################################
toLatex(sessionInfo())


