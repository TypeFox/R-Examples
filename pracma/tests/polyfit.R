###
### polyfit.R  +++ Test suite +++
###


test.polyfit <- function(input, expected) {
   output <- do.call(getFromNamespace("polyfit", "pracma"), input)
   identical(all.equal(output,
                       expected,
                       tolerance=1e-7),
             TRUE
    )
}

polyfit.expected.n1  <- c(1, 0)
polyfit.expected.n23 <- c(0, 1, 1, 1)
polyfit.expected.n4  <- c(-1, 0, 7, 0, 0) / 6
polyfit.expected.mat <- c(0, 1, -14, 65, -112, 60) / 12

test.polyfit(list(x=c(1,2,3), y=c(1,2,3)), polyfit.expected.n1)
test.polyfit(list(x=c(-2,-1,0,1,2), y=c(3,1,1,3,7), n=3), polyfit.expected.n23)
test.polyfit(list(x=c(-2,-1,0,1,2), y=c(2,1,0,1,2), n=4), polyfit.expected.n4)
test.polyfit(list(x=matrix(1:6, nrow=2, ncol=3, byrow=TRUE),
                  y=matrix(c(0,0,1,1,0,0), nrow=2, ncol=3, byrow=TRUE),
                  n=5), polyfit.expected.mat)

