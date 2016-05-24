
# Example cornell0: raw data, 7 inputs, 1 response, no missing
#
library("sivipm")

# -------------------------------------
# READ DATA (raw data)
# cornell0 <-read.table("../inst/extdata/cornell0.txt",  header=TRUE)# essai direct
nvar <- 7
XCornell0 <- cornell0[,1:nvar]
YCornell0 <-as.data.frame( cornell0[,8])
colnames(YCornell0)="Y"
# -------------------------------------
# CREATION OF THE  polyX OBJECT STRUCTURE 

# 1/   MONOME VECTOR: monomes coded by the inputs numbers
monomes <- c("1","2","3", "4","5", "6", "7",
             "1*3", "2*2", "2*4", "3*4", "5*5",   
             "6*6", "7*7*7")
#  Creation of an object of class 'polyX'
P1 <- vect2polyX(XCornell0, monomes)

# 2/ MONOMES OF STANDARD TYPE
# Example: complete polynome degree 2
# Creation of an object of class 'polyX'
PP3 <- crpolyX(XCornell0,2, type="full")
print(PP3)
# complete polynome degree 3
PP4 <- crpolyX(XCornell0,3, type="full")
print(PP4)
# -------------------------------------
# Illustration of  methods of class 'polyX'
 summary(P1)
P2 <- vect2polyX(XCornell0, c(as.character(1:7), "3*3*3", "3*3"))
z2 <- bind.polyX(P1, P2 )
print(z2)


# -------------------------------------
# CALCULATIONS
nc <- 2

print("TSIVIP without alea")
A <- sivipm(YCornell0, P1,  nc, options="tsivip")
getNames(A)
print(A, all=TRUE)
# compute tsivip Y by Y: here same results, because there is only one Y
print(apply(YCornell0, 2, sivipm, P1, nc))
print("TSIVIP without alea, full polynome degree 2")
A <- sivipm(YCornell0, PP3,  nc, options=c("tsivip", "fo.isivip", "simca"," lazraq"))
print(A, all=TRUE)
print("TSIVIP without alea, full polynome degree 3")
A <- sivipm(YCornell0, PP4,  nc, options=c("tsivip", "fo.isivip", "simca"," lazraq"))
print(A, all=TRUE)

print("TSIVIP with alea")
set.seed(15)
A <- sivipm(YCornell0, P1,  nc,  alea=TRUE, options="tsivip")
print(A, all=TRUE)


print("ISIVIP ")
A <- sivipm(YCornell0, P1,  nc, options="fo.isivip")
print(A, all=TRUE)


print("SIMCARULE")
A <- sivipm(YCornell0, P1,  nc, options=c("simca"," lazraq"))
print(A, all=TRUE)



print("ALL RESULTS IN A SINGLE INVOKATION")
print(sivipm(YCornell0, P1,  nc), all=TRUE)


print("BOOTSTRAP")
B=2
set.seed(15)
A <- sivipboot(YCornell0,P1, B , nc, alpha=0.05)
print(A)

set.seed(15)
B <- sivipboot(YCornell0,P1, B , nc, fast=T,  alpha=0.05)
print(all.equal(A,B))


