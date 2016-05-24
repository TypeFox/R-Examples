
##
## example data and test steps from pedigree.shrink
## Jason Sinnwell
##

require(kinship2)


data(minnbreast)
pedMN <- with(minnbreast, pedigree(id, fatherid, motherid, sex,famid=famid,
                         affected=cbind(cancer, bcpc, proband)))



## this pedigree as one person with cancer. The pedigree is not informative
## if they are the only available, so pedigree.shrink trims all.
## This caused an error in pedigree.shrink before kinship2. v1.2.8. Now fixed
mn2 <- pedMN[2]
#plot(mn2)

## breaks in pedigree.trim
shrink.mn2 <- pedigree.shrink(mn2,
               avail=ifelse(is.na(mn2$affected[,1]), 0, mn2$affected[,1]))
shrink.mn2

mnf8 <- pedMN['8']
#plot(mnf8)
shrink.mnf8 <- pedigree.shrink(mnf8,
              avail=ifelse(is.na(mnf8$affected[,1]), 0, mnf8$affected[,1]))

shrink.mnf8


## use sample.ped from the package
data(sample.ped)

pedAll <- pedigree(sample.ped$id, sample.ped$father, sample.ped$mother, 
       sample.ped$sex,
       affected=cbind(sample.ped$affected, sample.ped$avail), 
       famid=sample.ped$ped)


ped1 <- pedAll['1']

color1 <- sample.ped$avail[1:41] + 1

if(0) {
  pdf("sampleped1.pdf")
   plot(ped1, col=color1)
   plot(ped1, col=color1, pconnect=1.5)
  dev.off()
}


ped2 <- pedAll['2']

ped2$sex[c(13,12)] <- c("unknown", "terminated")

col2 <- ped2$affected[,2] + 1  ## red for availability, 

## set 2nd col of affected to NA
ped2$affected[c(7,9),2] <- NA

## show diamond and triangle for different sex codes.
## also show 2 shadings of affected, with ? as NA in affected matrix.
if(0) {
plot(ped2, col=col2)
}

set.seed(10)
shrink1.avail.B32 <- pedigree.shrink(ped=ped1, avail=ped1$affected[,2], maxBits=32)

set.seed(10)
shrink1.avail.B25 <- pedigree.shrink(ped=ped1, avail=ped1$affected[,2], maxBits=25)

shrink1.avail.B32$idTrimmed
## 101 102 107 108 111 121 122 123 131 132 134 139
shrink1.avail.B25$idTrimmed  
## 101 102 107 108 111 121 122 123 131 132 134 139 125 126

print(shrink1.avail.B32)
print(shrink1.avail.B25)


#Pedigree Size:
#                 N.subj Bits
#Original             41   49
#Only Informative     29   31
#Trimmed              26   25

# Unavailable subjects trimmed:
# 101 102 107 108 111 121 122 123 131 132 134 139 
#
# Informative subjects trimmed:
# 125 126

ped1df <- as.data.frame(ped1)

ped1df$idchar <- gsub("^1","A-", as.character(ped1df$id))
ped1df$dadidchar <- gsub("^1","A-", as.character(ped1df$dadid))
ped1df$momidchar <- gsub("^1","A-", as.character(ped1df$momid))
#ped1df$dadidchar <- ifelse(ped1df$dadidchar=="0", NA, ped1df$dadidchar)
#ped1df$momidchar <- ifelse(ped1df$momidchar=="0", NA, ped1df$momidchar)
ped1char <- with(ped1df, pedigree(idchar, dadidchar, momidchar, sex, affected,missid=c("0")))

options(stringsAsFactors=TRUE)
set.seed(10)
shrink1.p1char.B32 <- pedigree.shrink(ped=ped1char, avail=ped1char$affected[,2], maxBits=32)
shrink1.p1char.B32$idTrimmed
shrink1.avail.B32$idTrimmed

set.seed(10)
shrink1.p1char.B25 <- pedigree.shrink(ped=ped1char, avail=ped1char$affected[,2], maxBits=25)
shrink1.p1char.B25$idTrimmed
shrink1.avail.B25$idTrimmed
