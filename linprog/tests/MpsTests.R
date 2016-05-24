library( linprog )

## example of Steinhauser, Langbehn and Peters (1992)
## Production activities
cvec <- c(1800, 600, 600)  # gross margins
names(cvec) <- c("Cows","Bulls","Pigs")

## Constraints (quasi-fix factors)
bvec <- c(40, 90, 2500)  # endowment
names(bvec) <- c("Land","Stable","Labor")

## Needs of Production activities
Amat <- rbind( c(  0.7,   0.35,   0 ),
               c(  1.5,   1,      3 ),
               c( 50,    12.5,   20 ) )

## solve the model
result1a <- solveLP( cvec, bvec, Amat, TRUE )

## Write to a (virtual) MPS file
mpsFile <- file()
writeMps( mpsFile, cvec, bvec, Amat, "Steinhauser" )

## write the lines of this file to the output file
mpsLines <- readLines( mpsFile )
close( mpsFile )
print( mpsLines )

## Write to a (virtual) MPS file again (for readMps)
mpsFile <- file()
writeMps( mpsFile, cvec, bvec, Amat, "Steinhauser" )

## delete all LP objects
rm( cvec, bvec, Amat )

## Read LP data from MPS file and solve it.
lpModel <- readMps( mpsFile, TRUE, TRUE )
close( mpsFile )

## Print the model and its result
lpModel
all.equal( result1a, lpModel$res )


## example 1.1.3 of Witte, Deppe and Born (1975)
## Two types of Feed
cvec <- c(2.5, 2 )  # prices of feed
names(cvec) <- c("Feed1","Feed2")

## Constraints (minimum (<0) and maximum (>0) contents)
bvec <- c(-10, -1.5, 12)
names(bvec) <- c("Protein","Fat","Fibre")

## Matrix A
Amat <- rbind( c( -1.6,  -2.4 ),
               c( -0.5,  -0.2 ),
               c(  2.0,   2.0 ) )

## solve the model
result2a <- solveLP( cvec, bvec, Amat )

## Write to a (virtual) MPS file
mpsFile <- file()
writeMps( mpsFile, cvec, bvec, Amat, "Steinhauser" )

## write the lines of this file to the output file
mpsLines <- readLines( mpsFile )
close( mpsFile )
print( mpsLines )

## Write to a (virtual) MPS file again (for readMps)
mpsFile <- file()
writeMps( mpsFile, cvec, bvec, Amat, "Steinhauser" )

## delete all LP objects
rm( cvec, bvec, Amat )

## Read LP data from MPS file and solve it.
lpModel <- readMps( mpsFile, TRUE )
close( mpsFile )

## Print the model and its result
lpModel
all.equal( result2a, lpModel$res )
