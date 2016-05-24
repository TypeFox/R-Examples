
## test pedigree from bioinformatics manuscript
## try x-chrom kinship
## also has inbreeding and twins, for quick check
require(kinship2)
ped2mat <- matrix(c(1,1,0,0,1,
                    1,2,0,0,2,
                    1,3,1,2,1,
                    1,4,1,2,2,
                    1,5,0,0,2,
                    1,6,0,0,1,
                    1,7,3,5,2,
                    1,8,6,4,1,
                    1,9,6,4,1,
                    1,10,8,7,2),ncol=5,byrow=TRUE)

ped2df <- as.data.frame(ped2mat)
names(ped2df) <- c("fam", "id", "dad", "mom", "sex")
                  ## 1 2  3 4 5 6 7 8 9 10,11,12,13,14,15,16
ped2df$disease=   c(NA,NA,1,0,0,0,0,1,1,1)
ped2df$smoker=     c(0,NA,0,0,1,1,1,0,0,0)
ped2df$availstatus=c(0,0, 1,1,0,1,1,1,1,1)
ped2df$vitalstatus=c(1,1, 1,0,1,0,0,0,0,0)

ped2 <- with(ped2df, pedigree(id, dad, mom, sex, status=vitalstatus,
         affected=cbind(disease,smoker, availstatus), relation=matrix(c(8,9,1),ncol=3)))

## regular kinship matrix
kinship(ped2)

kinship(ped2, chr="X")

ped2$sex[9] <- "unknown"

## regular again, should be same as above
kinship(ped2)

## now with unknown sex, gets NAs
kinship(ped2, chrtype="X")

ped2$sex[9]="unknown"
kinship(ped2, chrtype="x")


# all descendants of sex=unknown to be NAs as well
ped2$sex[8]="unknown"
kinship(ped2, chr="X")


## testing kinship2 on pedigreeList when only one subject in a family
peddf <- rbind(ped2df, c(2,1,0,0,1,1,0,1,0)) 

peds <- with(peddf, pedigree(id, dad, mom, sex, status=vitalstatus,fam=fam,
         affected=cbind(disease,smoker, availstatus)))

kinfam <- kinship(peds)

kinfam

## now add two more for ped2, and check again
peddf <- rbind(peddf, c(2,2,0,0,2,1,0,1,0),c(2,3,1,2,1,1,0,1,0))

peds <- with(peddf, pedigree(id, dad, mom, sex, status=vitalstatus,fam=fam,
         affected=cbind(disease,smoker, availstatus)))

kin2fam <- kinship(peds)

kin2fam
