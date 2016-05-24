###
### $Id: padarray.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.padarray <- function(input, expected) {
    output <- do.call(getFromNamespace("padarray", "matlab"), input)
    identical(output, expected)
}


## Add three elements of padding to beginning of vector
## The padding elements contain mirror copies of the array
padarray.expected.mat4x4 <- matrix(1:4, 4, 4, byrow = TRUE)
test.padarray(list(A = 1:4,
                   padsize = 3,
                   padval = "symmetric",
                   direction = "pre"),
              padarray.expected.mat4x4)

## Add three elements of padding to the end of the first dimension of array
## and two elements of padding to the end of the second dimension.
## Use value of the last array element as the padding value.
Amat.2x2 <- matrix(as.numeric(1:4), nrow = 2, ncol = 2, byrow = TRUE)
padarray.expected.mat5x4 <- matrix(c(1, 2, 2, 2,
                                     3, 4, 4, 4,
                                     3, 4, 4, 4,
                                     3, 4, 4, 4,
                                     3, 4, 4, 4),
                                   nrow = 5, ncol = 4, byrow = TRUE)
test.padarray(list(A = Amat.2x2,
                   padsize = c(3, 2),
                   padval = "replicate",
                   direction = "post"),
              padarray.expected.mat5x4)

## Add three elements of padding to each dimension of a three-dimensional array.
## Each pad element contains the value zero.
Bmat.2x2 <- matrix(5:8, nrow = 2, ncol = 2, byrow = TRUE)
Carr.2x2x2 <- array(c(Amat.2x2, Bmat.2x2), c(2, 2, 2))
padarray.expected.arr8x8x2 <- {
                                  A <- array(0, c(8, 8, 2))
                                  A[4:5, 4:5, 1] <- Amat.2x2
                                  A[4:5, 4:5, 2] <- Bmat.2x2
                                  A
                              }
test.padarray(list(A = Carr.2x2x2,
                   padsize = c(3, 3),
                   padval = 0,
                   direction = "both"),
              padarray.expected.arr8x8x2)

## Add three elements of padding to end of vector
## The padding elements contain a mirror copy of the vector
padarray.expected.vec <- c(letters[1:5], rev(letters[3:5]))
test.padarray(list(A = letters[1:5],
                   padsize = c(0, 3),
                   padval = "symmetric",
                   direction = "post"),
              padarray.expected.vec)

