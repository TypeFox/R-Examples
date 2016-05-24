## packages
library("partykit")
library("rpart")

## data-generating process
dgp <- function(n)
    data.frame(y = gl(4, n), x1 = rnorm(4 * n), x2 = rnorm(4 * n))

## rpart check
learn <- dgp(100)
fit <- as.party(rpart(y ~ ., data = learn))
test <- dgp(100000)
system.time(id <- fitted_node(node_party(fit), test))
system.time(yhat <- predict_party(fit, id = id, newdata = test))

### predictions in info slots
tmp <- data.frame(x = rnorm(100))
pfit <- party(node = partynode(1L, split = partysplit(1L, breaks = 0), 
              kids = list(partynode(2L, info = -0.5), partynode(3L, info = 0.5))), data = tmp)
pfit
p <- predict(pfit, newdata = tmp)
p
table(p, sign(tmp$x))
