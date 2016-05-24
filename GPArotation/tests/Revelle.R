# This tests fix for an error caused by an exact initial setting.
#  (from William Revelle)

require("GPArotation")

f3 <- structure(c(0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0,0),
      .Dim = c(6L, 3L), .Dimnames = list(NULL, c("PC1", "PC2", "PC3")))

f3

#      PC1 PC2 PC3
#[1,]   0   0   1
#[2,]   0   1   0
#[3,]   1   0   0
#[4,]   0   0   1
#[5,]   0   1   0
#[6,]   1   0   0

# These previously gave object 'VgQt' not found
GPForth(f3)
Varimax(f3)
