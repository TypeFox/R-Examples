##
##  r e c t i n t . R  Tests
##

rectint <- pracma::rectint

x <- matrix(c(0, 0, 1, 1), ncol = 4)
y <- matrix(c( 0.75,-0.25, 0.5, 0.5,
               0.75, 0.25, 0.5, 0.5,
               0.75, 0.25, 0.2, 0.5,
               0.75, 0.75, 0.5, 0.5,
               0.75,-0.25, 0.5, 1.5), ncol = 4, byrow = TRUE)

all.equal(rectint(x, y),
          matrix(c(0.0625, 0.125, 0.1, 0.0625, 0.25), nrow = 1))
