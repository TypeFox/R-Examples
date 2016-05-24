
require(kinship2)
#
# Here is a case where the levels fail to line up properly
#

#data(testped2)
data(sample.ped)



# rearrange the founders to get a nicer plot
df1<- sample.ped[sample.ped$ped==1,]

ped1 <- with(df1, pedigree(id, father, mother, sex, affected))

plot(ped1)

df1reord <- df1[c(35:41,1:34),]
ped1reord <- with(df1reord, pedigree(id, father, mother, 
       sex, affected=affected))


plot(ped1reord, col=df1reord$avail+1)


# Two brothers married two sisters, which is currently "too much" for
#  the kindepth routine.
