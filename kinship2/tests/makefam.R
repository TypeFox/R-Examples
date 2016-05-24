library(kinship2)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
#
# A pedigree that can cause problems, due to a double marriage
#
id <- 1:20
mom<- c(0,0,0,2,2,2,0,2,0, 0,2,2,0,2,0,2, 7,7, 11,14)
dad<- c(0,0,0,1,1,1,0,1,0, 0,3,3,0,3,0,3, 8,8, 10,13)

temp<- makefamid(id, mom, dad)
aeq(temp, pmin(temp,1))
