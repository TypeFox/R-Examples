require(DoE.wrapper)
set.seed(1212)
## test utility functions from DoE.base for designs that are not in
## DoE.base

plan <- bbd.design(3)
y <- rnorm(nrow(plan))
lm(y~., plan)
lm(y~., qua.design(plan))
lm(y~., qua.design(plan, quantitative="none",
      contrasts=c(A="contr.treatment",B="contr.treatment",C="contr.treatment")))

plan <- ccd.design(4)
y <- rnorm(nrow(plan))
lm(y~., plan)
lm(y~., qua.design(plan))
lm(y~., qua.design(plan, quantitative="none",
      contrasts=c(X1="contr.treatment",X2="contr.treatment",
            X3="contr.treatment",X4="contr.helmert")))

plan <- FrF2(16,4,blocks=2)
y <- rnorm(nrow(plan))
lm(y~., plan)
lm(y~., qua.design(plan))
lm(y~., qua.design(plan, quantitative="none",
      contrasts=c(A="contr.treatment",B="contr.treatment",
            C="contr.treatment",D="contr.helmert")))

plan <- add.response(plan, y)
lm(y~., plan)
lm(y~., qua.design(plan))
lm(y~., qua.design(plan, quantitative="none",
      contrasts=c(A="contr.treatment",B="contr.treatment",
            C="contr.treatment",D="contr.helmert")))
