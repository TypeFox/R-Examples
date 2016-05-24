#
# The test case from Lange, chapter 5
#
library(kinship2)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

id <- 1:6
momid <- c(0,0,2,2,4,4)
dadid <- c(0,0,1,1,3,3)
sex <- c(1,2,1,2,1,2)
kmat1 <- kinship(id, dadid, momid)
aeq(as.matrix(kmat1), c(4,0,2,2,2,2, 0,4,2,2,2,2, 2,2,4,2,3,3, 
          2,2,2,4,3,3, 2,2,3,3,5,3, 2,2,3,3,3,5) /8)

kmat1x <- kinship(id, dadid, momid, sex, chrtype='X')
aeq(kmat1x, c(8,0,0,4,4,2, 0,4,4,2,2,3, 0,4,8,2,2,5, 4,2,2,4,4,3,
              4,2,2,4,8,3, 2,3,5,3,3,5)/8)


# And here is an an odd one with cross marriages, but no inbreeding
#
test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
#with(test1, plot(pedigree(id, dad, mom, sex)))

temp <- diag(14)
temp[1:2, 5:6] <- temp[3:4, 7:8] <- temp[c(3,6),9] <- 4  #children
temp[c(2,7),10] <- temp[11:12,13] <- temp[c(10,13),14] <- 4 #children
temp[5,6] <- temp[7,8] <- 4   #siblings
temp[1:2, 9] <- temp[3:4,10] <- temp[c(2,7),14]<- temp[11:12,14] <-2 #grandkids
temp[3:4, 14] <- 1            #great-grandchild
temp[5,9] <- temp[8,10] <- 2  #uncle or aunt
temp[8,14] <- 1               #great aunt
temp[5:6, 10] <- temp[7:8, 9] <- 2  # one shared parent
temp[5:6, 14] <- 1  # 2 is parent to 5/6, grandparent to 14
temp[9,14] <-1  # related through both 3 and 2
temp[9,10] <- 2
temp <- temp + t(temp)
diag(temp) <- 8
kmat2 <- kinship(test1$id, test1$mom, test1$dad)
aeq(as.matrix(kmat2), temp/16)


#And now a pedigree list
test2 <- rbind(test1, data.frame(id=1:6, mom= c(0,0,2,2,4,4),
                                 dad=c(0,0,1,1,3,3), 
                                 sex=c(0,1,0,1,1,1)))
test2$famid <- rep(1:2, c(nrow(test1), 6))

ped <- with(test2, pedigree(id=id, mom=mom, dad=dad, sex=sex, fam=famid))

kmat3 <- kinship(ped)
n1 <- nrow(test1)
aeq(as.matrix(kmat2), as.matrix(kmat3[1:n1, 1:n1]))
aeq(as.matrix(kmat1), as.matrix(kmat3[1:6+n1, 1:6+n1]))
all(as.matrix(kmat3[1:n1, 1:6+n1]) ==0)


kmat4 <- with(test2, makekinship(famid, paste(famid, id, sep='/'), 
                                 paste(famid, dad, sep='/'),
                                 paste(famid, mom, sep='/')))
all.equal(kmat3, kmat4)


# scramble the people, to make sure that the routine
#  keeps labels correctly
set.seed(1953)
temp <- order(runif(nrow(test2)))
ped2 <- with(test2[temp,], 
             pedigree(id=id, mom=mom, dad=dad, sex=sex, fam=famid)) 
kmat5 <- kinship(ped2)

indx <- match(dimnames(kmat3)[[2]], dimnames(kmat5)[[2]])
all.equal(kmat5[indx, indx], kmat3)
