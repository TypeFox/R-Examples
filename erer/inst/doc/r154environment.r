# A. Function environment
library(erer)
test <- function(data, ...) {
  y <- data - 1
  dog <- function(z) {data + z}           # free variable "data"
  k1 <- dog(z = y)
  k2 <- round(x = data, digits = digits)  # "digits" from global
  k3 <- round(x = data, ...)              # "digits" from local
  out <- data.frame(y, k1, k2, k3, data, glob.digit = digits, ...)
  return(out)
}
digits <- 2
(res1 <- test(data = 6.1234, digits = 3))
(res2 <- test(data = 6.1234))

class(test); typeof(test); environment(test)
class(res1); typeof(res1); environment(res1)

en <- environment(fun = aiStaFit); en
ls(name = environment(fun = test)); ls(name = en); search()

# -------------------------------------------------------------------------
# B. deparse(); substitute(); get()
numb <- c(4, 10, 23)
lett <- c("aa", "bb", "cc")
numb2 <- get(x = "numb"); identical(numb, numb2)  # TRUE

deparse(numb); str(deparse(numb))  # content; character
substitute(numb); str(substitute(numb))  # name; symbol
deparse(substitute(numb)); str(deparse(substitute(numb))) # name; character

wood <- function(num, let) {
  n1 <- deparse(num)
  n2 <- substitute(num)
  n3 <- deparse(substitute(num))
  n4 <- match.call()
  n5 <- sapply(X = n4, FUN = deparse)
  return(listn(n1, n2, n3, n4, n5))
}
(ch <- wood(num = numb, let = lett))

# -------------------------------------------------------------------------
# C. match.arg(); match.call(); update()
deer <- function(c1, c2, c3, c4, c5 = c("Mon", "Wed", "Fri"), ...) {
  key <- c1 + c2 + c3 + c4
  c5 <- match.arg(arg = c5, choices = c("Mon", "Wed", "Fri"))
  c6 <- match.arg(c5)
  c7 <- match.call(definition = deer, expand.dots = TRUE)
  c8 <- match.call(definition = deer, expand.dots = FALSE)
  call <- sys.call()  # This has to be named as "call".
  result <- listn(key, c5, c6, c7, c8, call)
  return(result)
}
(sh <- deer(c1 = 1, c2 = 2, c3 = 3, c4 = 4, extra = "fish"))
str(sh); sh2 <- update(sh, c1 = 71); sh2[c("key", "call")]

# -------------------------------------------------------------------------
# D. model.frame(); terms()
data(daIns)
ff <- Y ~ Injury + HuntYrs + Edu + Inc
ra <- lm(formula = ff, data = daIns)
round(coef(ra), 5)

butter <- function(formula, data) {
  small <- model.frame(formula, data = data)
  y <- model.response(small)
  x <- model.matrix(object = formula, data = small)
  coeff <- solve(t(x) %*% x) %*% t(x) %*% y
  out <- listn(small, y, x, coeff)
  return(out)
}
ww <- butter(formula = ff, data = daIns)
round(t(ww[["coeff"]]), 5)
names(attributes(ww$small))
terms(ww$small)

# -------------------------------------------------------------------------
# D. expression(); eval()
math <- expression(x ^ 2 + y); class(math); str(math)
x <- 1; y <- 3; ma <- eval(math) + 30; ma
x <- 5; y <- 6; mb <- eval(math) + 50; mb