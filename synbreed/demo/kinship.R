#############################################
## Demo for pedigree-based kinship coefficients
##
##
## author : Valentin Wimmer
## date : 2012 - 03 - 07
##
##############################################



### Example from Legarra et al. (2009), J. Dairy Sci. 92: p. 4660
id <- 1:17
par1 <- c(0,0,0,0,0,0,0,0,1,3,5,7,9,11,4,13,13)
par2 <- c(0,0,0,0,0,0,0,0,2,4,6,8,10,12,11,15,14)

# create object of class pedigree
ped <- create.pedigree(id,par1,par2)
summary(ped)

# create object of class gpData
gp <- create.gpData(pedigree=ped)

# plot of the pedigree structure in a graph
plot(ped)

# gametic relationship
G <- kin(gp,ret="gam")
dim(G)
# 2n x 2n

# additive relationship
A <- kin(gp,ret="add")
dim(A)
# n x n

# extract inbreeding coefficients
diag(A) - 1


# connection to gametic relationship
A[16,17]
# sum over allele IBD
0.5*sum(G[c("16_1","16_2"),c("17_1","17_2")])


# dominance relationship
D <- kin(gp,ret="dom")
dim(D)
# n x n

# connection to gametic relationship
D[16,17]
# sum over allele IBD
G["16_1","17_1"] *  G["16_2","17_2"] + G["16_1","17_2"] *  G["16_1","17_2"]

plot(G)
plot(A)
plot(D)