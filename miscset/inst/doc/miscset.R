## ----eval=FALSE----------------------------------------------------------
#  # from CRAN
#  install.packages("miscset")
#  # from github - latest version
#  devtools::install_github("setempler/miscset", build_vignettes = TRUE)

## ------------------------------------------------------------------------
library(miscset)
library(ggplot2)
d <- data.frame(a=c(2,1,3,NA,1), b=2:6, c=5:1)
m <- matrix(letters[1:9], 3, 3, dimnames = list(1:3,1:3))

## ----fig.height=3, fig.width=4, fig.align='center'-----------------------
ciplot(d)

## ----fig.height=2, fig.width=6, fig.align='center'-----------------------
ggplotGrid(list(
  ggplot(d, aes(x=b,y=-c,col=b)) + geom_line(),
  ggplot(d, aes(x=b,y=-c,shape=factor(b))) + geom_point()),
  ncol = 2)

## ----comment=NA, fig.height=3, fig.width=4, fig.align='center'-----------
n <- length(d)
gghcl(n)
ciplot(d, col = gghcl(n))

## ----fig.height=2, fig.width=2, fig.align='center'-----------------------
plotn()

## ---- comment=NA---------------------------------------------------------
d
sort(d, by = c("a", "c"))

## ----comment=NA----------------------------------------------------------
d[1:3,]
do.rbind(list(first=d[1:2,], second=d[1:3,]))

## ----comment=NA----------------------------------------------------------
m
enpaire(m)

## ----comment=NA----------------------------------------------------------
m[-1,]
squarematrix(m[-1,])

## ----comment=NA----------------------------------------------------------
textable(d, caption = 'miscset vignette example data.frame', as.document = TRUE)

## ----comment=NA----------------------------------------------------------
lsall()

## ----comment=NA----------------------------------------------------------
mgrepl(c("a","b"), c("ab","ac","bc"), any)
mgrepl(c("a","b"), c("ab","ac","bc"), all)
mgrepl(c("a","b"), c("ab","ac","bc"), all, use.which = TRUE)
mgrepl(c("a","b"), c("ab","ac","bc"), identity)

## ----comment=NA----------------------------------------------------------
gregexprind(c("a"), c("ababa","ab","xyz",NA), 1)
gregexprind(c("a"), c("ababa","ab","xyz",NA), 2)
gregexprind(c("a"), c("ababa","ab","xyz",NA), "last")

## ----comment=NA----------------------------------------------------------
leading0(c(9, 112, 5009))

## ----comment=NA----------------------------------------------------------
strextr("xa,xb,xn,ya,yb", "n$", ",")
strextr("xa,xb,xn,ya,yb", "^x", ",", mult=T)

## ----comment=NA----------------------------------------------------------
strpart("xa,xb,xn,ya,yb", ",", 3)

## ----comment=NA----------------------------------------------------------
strrev(c("olleH", "!dlroW"))

## ----comment=NA----------------------------------------------------------
data.frame(
  duplicate = d$a,
  ".d" = duplicated(d$a),
  ".s" = duplicates(d$a),
  ".i" = duplicatei(d$a))

## ----comment=NA----------------------------------------------------------
p2star(c(0.003, 0.049, 0.092, 0.431))

## ----comment=NA----------------------------------------------------------
d$a
confint(d$a, ret.attr = FALSE)

## ----comment=NA----------------------------------------------------------
ntri(12)

## ----comment=NA----------------------------------------------------------
d$c
scale0(d$c)
scaler(d$c, c(2, 6), b = c(1, 10))

## ----comment=NA----------------------------------------------------------
d$a
nunique(d$a)
nunique(d$a, FALSE)
uniquei(d$a)
uniquei(d$a, FALSE)

