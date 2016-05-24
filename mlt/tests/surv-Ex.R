
library("mlt")
library("survival")
set.seed(29)

y <- rexp(10)
d <- rep(c(TRUE, FALSE), length = length(y))
### right
(s <- Surv(y, d, type = "right"))
attr(s, "type")
cbind(s, R(s))

### right censored, left truncated
(s <- Surv(y, y + 1, d))
attr(s, "type")
cbind(s, R(s))

### left censored
(s <- Surv(y, d, type = "left"))
attr(s, "type")
cbind(s, R(s))

### interval
dd <- rep(0:3, length = length(y))
(s <- Surv(y, y + 1, dd, type = "interval"))
attr(s, "type")
cbind(s, R(s))

### interval2
(s <- Surv(y, y + 1, type = "interval2"))
attr(s, "type")
cbind(s, R(s))

(s <- Surv(y, ifelse(d, Inf, y + 1), type = "interval2"))
attr(s, "type")
cbind(s, R(s))

(s <- Surv(y, ifelse(d, NA, y + 1), type = "interval2"))
attr(s, "type")
cbind(s, R(s))
