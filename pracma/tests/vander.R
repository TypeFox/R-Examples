###
### VANDER.R  +++ Test suite +++
###

test.vander <- function(input, expected) {
    output <- do.call(getFromNamespace("vander", "pracma"), input)
    identical(output, expected)
}

vander.expected.empty <- matrix(as.numeric(NA), nrow=0, ncol=0)
vander.expected.scalar <- matrix(1.0, nrow=1, ncol=1)
vander.expected.seq3 <- matrix(c( 1,       1,     1,    1,   1,
                                  5.0625,  3.375, 2.25, 1.5, 1,
                                 16,       8,     4,    2,   1,
                                 39.0625, 15.625, 6.25, 2.5, 1,
                                 81,      27,     9,    3,   1),
                               nrow=5, ncol=5, byrow=TRUE)
vander.expected.complex <- matrix(c(-1, 0.0 + 1i, 1.0,
                                    -4, 0.0 + 2i, 1.0,
                                    -9, 0.0 + 3i, 1.0),
                                  nrow=3, ncol=3, byrow=TRUE)

test.vander(list(x=numeric()), vander.expected.empty)
test.vander(list(x=1), vander.expected.scalar)
test.vander(list(x=seq(from=1, to=3, by=0.5)), vander.expected.seq3)
test.vander(list(x=c(1,2,3)*1i), vander.expected.complex)
