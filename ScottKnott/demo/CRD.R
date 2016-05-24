##
## Examples: Completely Randomized Design (CRD)
##

## The parameters can be: vectors, design matrix and the response variable,
## data.frame or aov

## Example 1
library(ScottKnott)
data(CRD1)

## From: vectors x and y
sk1 <- with(CRD1,
            SK(x=x,
               y=y,
               model='y ~ x',
               which='x'))
summary(sk1)
plot(sk1)

## From: design matrix (dm) and response variable (y)
sk2 <- with(CRD1,
            SK(x=dm,
               y=y,
               model='y ~ x',
               which='x',
               dispersion='s'))
summary(sk2)
plot(sk2,
     pch=15,
     col=c('blue', 'red'),
     mm.lty=4,
     ylab='Response',
     title=NULL)

## From: data.frame (dfm)
sk3 <- with(CRD1,
            SK(x=dfm,
               model='y ~ x',
               which='x',
               dispersion='se'))
summary(sk3)
plot(sk3,
     mm.lty=3,
     id.col=FALSE,
     title=NULL)

## From: aov
av1 <- with(CRD1,
            aov(y ~ x,
                data=dfm))
summary(av1)

sk4 <- SK(x=av1,
          which='x')
summary(sk4)
plot(sk4, title=NULL)


## Example 2
library(ScottKnott)
data(CRD2)

## From: vectors x and y
sk5 <- with(CRD2,
            SK(x=x,
               y=y,
               model='y ~ x',
               which='x'))
summary(sk5)
plot(sk5,
     id.las=2,
     rl=FALSE)

## From: design matrix (dm) and response variable (y)
sk6 <- with(CRD2,
            SK(x=dm,
               y=y,
               model='y ~ x',
               which='x',
               sig.level=0.005))
summary(sk6)
plot(sk6,
     col=rainbow(max(sk6$groups)),
     mm.lty=3,
     id.las=2,
     rl=FALSE,
     title='sig.level=0.005', )

## From: data.frame (dfm)
sk7 <- with(CRD2,
            SK(x=dfm,
               model='y ~ x',
               which='x'))
summary(sk7)
plot(sk7,
     col=rainbow(max(sk7$groups)),
     id.las=2,
     id.col=FALSE,
     rl=FALSE)

## From: aov
av2 <- with(CRD2,
            aov(y ~ x,
                data=dfm))
summary(av2)

sk8 <- SK(x=av2,
          which='x')
summary(sk8)
plot(sk8,
     col=rainbow(max(sk8$groups)),
     rl=FALSE, 
     id.las=2,
     id.col=FALSE,
     title=NULL)
