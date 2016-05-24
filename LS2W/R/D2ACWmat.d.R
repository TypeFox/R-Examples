`D2ACWmat.d` <-
function(J, filter.number = 10., family = "DaubLeAsymm", OPLENGTH = 100000.)
{
#
#
#This function is a latter version of D2ACWamat
#in which the no. of rows = the no. of cols!
#Thus it can assign HUGE matricies
#
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
###############
# Diagonal entries
###############
#
# Start off by placing the level -J coefficients in first.
#
m <- #
matrix(0., nrow = nr, ncol = nc)
tmp <- matrix(0., nrow = nrlocal, ncol = nc)
tmp <- P[[3. * J]]
#
#
# Next, we place the level -1, ..., -(J-1) in the matrix.
# More accurately, we stick them in the centre of the matrix.
#
m[((3. * J - 1.) * nc + 1.):(3. * J * nc),  ] <- tmp[1.:nc,  ]
####################################################
#We now do the same for the horizontal and diagonal directions.
####################################################
#
# Horizontal
#
for(j in (2. * J + 1.):(3. * J - 1.)) {
   nrj <- nrow(P[[j]])
   ncj <- ncol(P[[j]])
   ncz <- (nc - ncj)/2.
   nrz <- (nrlocal - nrj)/2.
   z1 <- matrix(0., nrow = nrj, ncol = ncz)
   z2 <- matrix(0., nrow = nrz, ncol = nc)
   tmp1 <- matrix(0., nrow = nrj, ncol = ncj)
   tmp1 <- P[[j]]
   m[((j - 1.) * nc + 1.):(j * nc),  ] <- rbind(z2, cbind(z1,
   tmp1[1.:nrj,  ], z1), z2)
}
tmp <- matrix(0., nrow = nrlocal, ncol = nc)
tmp <- P[[2. * J]]
m[((2. * J - 1.) * nc + 1.):(2. * J * nc),  ] <- tmp[1.:nc,  ]
#
# Vertical
#
for(j in (J + 1.):(2. * J - 1.)) {
   nrj <- nrow(P[[j]])
   ncj <- ncol(P[[j]])
   ncz <- (nc - ncj)/2.
   nrz <- (nrlocal - nrj)/2.
   z1 <- matrix(0., nrow = nrj, ncol = ncz)
   z2 <- matrix(0., nrow = nrz, ncol = nc)
   tmp1 <- matrix(0., nrow = nrj, ncol = ncj)
   tmp1 <- P[[j]]
   m[((j - 1.) * nc + 1.):(j * nc),  ] <- rbind(z2, cbind(z1,
   tmp1[1.:nrj,  ], z1), z2)
}
tmp <- matrix(0., nrow = nrlocal, ncol = nc)
tmp <- P[[J]]
m[((J - 1.) * nc + 1.):(J * nc),  ] <- tmp[1.:nc,  ]
for(j in 1.:(J - 1.)) {
   nrj <- nrow(P[[j]])
   ncj <- ncol(P[[j]])
   ncz <- (nc - ncj)/2.
   nrz <- (nrlocal - nrj)/2.
   z1 <- matrix(0., nrow = nrj, ncol = ncz)
   z2 <- matrix(0., nrow = nrz, ncol = nc)
   tmp1 <- matrix(0., nrow = nrj, ncol = ncj)
   tmp1 <- P[[j]]
   m[((j - 1.) * nc + 1.):(j * nc),  ] <- rbind(z2, cbind(z1,
   tmp1[1.:nrj,  ], z1), z2)
}
m
}

