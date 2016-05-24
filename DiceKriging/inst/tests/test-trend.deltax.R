d <- 3
X <- expand.grid(x1=seq(0,1,length=4), x2=seq(0,1,length=4), x3=seq(0,1,length=4))
fun <- function(x){
  (x[1]+2*x[2]+3*x[3])^2
}
y <- apply(X, 1, fun) 

x <- c(0.2, 0.4, 0.6)
coef.cov=c(0.5, 0.9, 1.3); coef.var=3

precision <- 1e-10  # the three following tests should work with it, since the computations are analytical

m <- km(~1, design=X, response=y, coef.cov=coef.cov, coef.var=coef.var)
grad.trend <- trend.deltax(x, m)
grad.trueValue <- rep(0,3)
test_that(desc="test trend.delta.x, formula=~1", expect_true(max(abs(grad.trend-grad.trueValue)) < precision))

m <- km(~., design=X, response=y, coef.cov=coef.cov, coef.var=coef.var)
grad.trend <- trend.deltax(x, m)
grad.trueValue <- rbind(matrix(0, ncol=d), diag(d))
test_that(desc="test trend.delta.x, formula=~.", expect_true(max(abs(grad.trend-grad.trueValue)) < precision))

m <- km(~.^2, design=X, response=y, coef.cov=coef.cov, coef.var=coef.var)
grad.trend <- trend.deltax(x, m)
grad.trueValue <- rbind(matrix(0, ncol=d), diag(d), c(x[2], x[1], 0), c(x[3], 0, x[1]), c(0, x[3], x[2]))
test_that(desc="test trend.delta.x, formula=~.^2", expect_true(max(abs(grad.trend-grad.trueValue)) < precision))

precision <- 1e-6   

m <- km(~I(cos(x1)*x2)+I(x1*sin(x2))+I(x2*x3^3), design=X, response=y, coef.cov=coef.cov, coef.var=coef.var)
grad.trend <- trend.deltax(x, m, h=1e-8)
grad.trueValue <- rbind(matrix(0, ncol=d), c(-sin(x[1])*x[2], cos(x[1]), 0),  
                        c(sin(x[2]), x[1]*cos(x[2]), 0), c(0, x[3]^3, x[2]*3*x[3]^2))
test_that(desc="test trend.delta.x, formula=~I(cos(x1)*x2)+I(x1*sin(x2))+I(x2*x3^3)", expect_true(max(abs(grad.trend-grad.trueValue)) < precision))