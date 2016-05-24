## ----setup, echo=FALSE, message=FALSE------------------------------------
require(knitr)
require(kfigr)
opts_chunk$set(message=FALSE, warning=FALSE)

## ----first-chunk, anchor="figure"----------------------------------------
require(ggplot2)
qplot(rnorm(100), geom="histogram")

## ----second-chunk, anchor="figure"---------------------------------------
qplot(runif(100), geom="density")

## ----third-chunk, anchor="table", results='asis'-------------------------
kable(head(iris, 6))

## ----fourth-chunk, eval=FALSE--------------------------------------------
#  x = 1:20
#  y = x + rnorm(20)
#  lm(y~x)

## ----fourth-chunk, anchor="block"----------------------------------------
x = 1:20
y = x + rnorm(20)
lm(y~x)

## ----fifth-chunk, echo=FALSE, eval=FALSE---------------------------------
#  df <- data.frame(x=x, y=y)
#  ggplot(df, aes(x=x, y=y)) + geom_smooth(method="lm") +
#  geom_point(pch=21, color="black", fill="red")

## ----fifth-code, ref.label='fifth-chunk', anchor="block", eval=FALSE-----
#  df <- data.frame(x=x, y=y)
#  ggplot(df, aes(x=x, y=y)) + geom_smooth(method="lm") +
#  geom_point(pch=21, color="black", fill="red")

## ----fifth-plot, ref.label='fifth-chunk', anchor="figure", echo=FALSE----
df <- data.frame(x=x, y=y)
ggplot(df, aes(x=x, y=y)) + geom_smooth(method="lm") + 
geom_point(pch=21, color="black", fill="red")

## ----last-chunk, anchor="block"------------------------------------------
anchors("index")
anchors("history")

