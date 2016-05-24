###
### $Id: meshgrid.R 51 2014-02-05 21:22:28Z plroebuck $
###


##-----------------------------------------------------------------------------
test.meshgrid <- function(input, expected) {
    output <- do.call(getFromNamespace("meshgrid", "matlab"), input)
    identical(output, expected)
}

x1 <- 1:3
y1 <- 10:14
meshgrid.expected.xy <- list(x = matrix(rep(x1, length(y1)),
                                        nrow = length(y1),
                                        ncol = length(x1), byrow = TRUE),
                             y = matrix(rep(y1, length(x1)),
                                        nrow = length(y1),
                                        ncol = length(x1)))

test.meshgrid(list(x = 0), meshgrid.expected.xy)
test.meshgrid(list(x = x1, y = y1), meshgrid.expected.xy)

x2 <- 5:8
y2 <- 10:14
z2 <- 2:3
meshgrid.expected.xyz <- list(x = array(matrix(rep(x2, length(y2)),
                                               nrow = length(y2),
                                               ncol = length(x2), byrow = TRUE),
                                        c(length(y2), length(x2), length(z2))),
                              y = array(rep(y2, length(x2)),
                                        c(length(y2), length(x2), length(z2))),
                              z = array(sapply(z2,
                                               function(val, len) rep(val, len),
                                               length(y2) * length(x2)),
                                        c(length(y2), length(x2), length(z2))))

test.meshgrid(list(x = x2, y = y2, z = z2, nargout = 3), meshgrid.expected.xyz)

