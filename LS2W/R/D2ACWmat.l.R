`D2ACWmat.l` <-
function(J, filter.number = 10., family = "DaubLeAsymm", OPLENGTH = 100000.)
{
#
#
#This function is a latter version of D2ACWbmat
#in which the no. of rows = the no. of cols!
#
#However, in contrast to DACWmat, we now do things 
#within scale (cf. direction previously)
####################################################
#
# Initial setup of the DACW  matrix
#
####################################################
J <-  - J
P <- D2ACW( - J, filter.number = filter.number, family = family, 
OPLENGTH = OPLENGTH)
nc <- ncol(P[[3. * J]])
nr <- 3. * J * nc
nrlocal <- nc
##############
#Scale -J entries
##############
#
#Start off by placing the level -J coefficients in first.
#
###
#Diagonal
m <- #
matrix(0., nrow = nr, ncol = nc)
tmp <- matrix(0., nrow = nrlocal, ncol = nc)
tmp <- P[[3. * J]]
#
###
#Horizontal
m[((3. * J - 1.) * nc + 1.):(3. * J * nc),  ] <- tmp[1.:nc,  ]
tmp <- matrix(0., nrow = nrlocal, ncol = nc)
tmp <- P[[2. * J]]
#
###
#Vertical
m[((3. * J - 2.) * nc + 1.):((3. * J - 1.) * nc),  ] <- tmp[1.:nc,
]
tmp <- matrix(0., nrow = nrlocal, ncol = nc)
tmp <- P[[J]]
#
#####################
#All the remaining scales
#####################
#
#We do this direction by direction for several reasons.
#Firstly, it'll save on the memory. Ie three small loops 
#are better than one small one. Secondly, it's far more logical 
#to do it this way. This is because 
#
#vertical level j  > location 3j-2
#horizontal level j  > location 3j-1
#diagonallevel j  > location 3j
#
#Vertical
m[((3. * J - 3.) * nc + 1.):((3. * J - 2.) * nc),  ] <- tmp[1.:nc,
]
#####
# Horizontal
for(j in 1.:(J - 1.)) {
   nrj <- nrow(P[[j]])
   ncj <- ncol(P[[j]])
   ncz <- (nc - ncj)/2.
   nrz <- (nrlocal - nrj)/2.
   z1 <- matrix(0., nrow = nrj, ncol = ncz)
   z2 <- matrix(0., nrow = nrz, ncol = nc)
   tmp1 <- matrix(0., nrow = nrj, ncol = ncj)
   tmp1 <- P[[j]]
   m[((3. * j - 3.) * nc + 1.):((3. * j - 2.) * nc),  ] <- rbind(z2, cbind(z1, tmp1[1.:nrj,  ], z1), z2)
}
######
#Diag
for(j in (J + 1.):(2. * J - 1.)) {
   i <- j - J
   nrj <- nrow(P[[j]])
   ncj <- ncol(P[[j]])
   ncz <- (nc - ncj)/2.
   nrz <- (nrlocal - nrj)/2.
   z1 <- matrix(0., nrow = nrj, ncol = ncz)
   z2 <- matrix(0., nrow = nrz, ncol = nc)
   tmp1 <- matrix(0., nrow = nrj, ncol = ncj)
   tmp1 <- P[[j]]
   m[((3. * i - 2.) * nc + 1.):((3. * i - 1.) * nc),  ] <- rbind(
   z2, cbind(z1, tmp1[1.:nrj,  ], z1), z2)
}
for(j in (2. * J + 1.):(3. * J - 1.)) {
   i <- j - (2. * J)
   nrj <- nrow(P[[j]])
   ncj <- ncol(P[[j]])
   ncz <- (nc - ncj)/2.
   nrz <- (nrlocal - nrj)/2.
   z1 <- matrix(0., nrow = nrj, ncol = ncz)
   z2 <- matrix(0., nrow = nrz, ncol = nc)
   tmp1 <- matrix(0., nrow = nrj, ncol = ncj)
   tmp1 <- P[[j]]
   m[((3. * i - 1.) * nc + 1.):(3. * i * nc),  ] <- rbind(z2,
   cbind(z1, tmp1[1.:nrj,  ], z1), z2)
}
m
}

