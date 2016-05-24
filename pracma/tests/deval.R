##
##  d e v a l . R  Test suite
##

deval <- pracma::deval
deeve <- pracma::deeve

x <- seq(0, 10*pi, len=100)
y <- zapsmall(sin(x))

all.equal(deval(x, y, c(-1e-5, 0, 1, 5, 10, 15, 20, 25, 30, x[100], 40)),
          as.matrix(c(        NA,
                       0.0000000,
                       0.8358028,
                      -0.9499175,
                      -0.5372202,
                       0.6442378,
                       0.9117673,
                      -0.1307134,
                      -0.9756776,
                       0.0000000,
                              NA)),
          tolerance = 1e-5
)

all.equal(deeve(x, y),
          c(0.000000, 3.141206, 6.282671, 9.424329, 12.566114, 15.707963,
            18.849812, 21.991597, 25.133255, 28.274720, 31.415927),
          tolerance = 1e-5)
