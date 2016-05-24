##
## Examples: Completely Randomized Design (CRD)
##

## The parameters can be: vectors, design matrix and the response variable,
## data.frame or aov

## Example 1
library(TukeyC)
data(CRD1)

## From: vectors x and y - balanced
tk1 <- with(CRD1,
            TukeyC(x=x,
                   y=y,
                   model='y ~ x',
                   which='x'))
tk1
summary(tk1)
plot(tk1)

## From: vectors x and y - unbalanced
tk1u <- with(CRD1,
             TukeyC(x=x[-1],
                    y=y[-1],
                    model='y ~ x',
                    which='x',
                    dispersion='s'))
tk1u
summary(tk1u)
plot(tk1u)

## From: design matrix (dm) and response variable (y) - balanced
tk2 <- with(CRD1,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ x',
                   which='x',
                   dispersion='se'))
tk2
summary(tk2)

## From: design matrix (dm) and response variable (y) - unbalanced
tk2u <- with(CRD1,
             TukeyC(x=dm[-1,],
                    y=y[-1],
                    model='y ~ x',
                    which='x'))
tk2u
summary(tk2u)

## From: data.frame (dfm) - balanced
tk3 <- with(CRD1,
            TukeyC(x=dfm,
                   model='y ~ x',
                   which='x'))
tk3
summary(tk3)

## From: data.frame (dfm) - unbalanced
tk3u <- with(CRD1,
             TukeyC(x=dfm[-1,],
                    model='y ~ x',
                    which='x'))
tk3u
summary(tk3u)

## From: aov - balanced
av1 <- with(CRD1,
            aov(y ~ x,
                data=dfm))
summary(av1)

tk4 <- TukeyC(x=av1,
              which='x')
tk4
summary(tk4)

## From: aov - unbalanced
av1u <- with(CRD1,
             aov(y ~ x,
                 data=dfm[-1,]))
summary(av1u)

tk4u <- TukeyC(x=av1u,
               which='x')
tk4u
summary(tk4u)


## Example 2
library(TukeyC)
data(CRD2)

## From: vectors x and y - balanced
tk5 <- with(CRD2,
            TukeyC(x=x,
                   y=y,
                   model='y ~ x',
                   which='x', round=3))
tk5
summary(tk5)
plot(tk5,
     id.las=2,
     rl=FALSE)

## From: vectors x and y - unbalanced
tk5u <- with(CRD2,
             TukeyC(x=x[-1],
                    y=y[-1],
                    model='y ~ x',
                    which='x',
                    round=3))
tk5u
summary(tk5u)

## From: design matrix (dm) and response variable (y) - balanced
tk6 <- with(CRD2,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ x',
                   which='x'))
tk6
summary(tk6)

## From: design matrix (dm) and response variable (y) - unbalanced
tk6u <- with(CRD2,
             TukeyC(x=dm[-1,],
                    y=y[-1],
                    model='y ~ x',
                    which='x'))
tk6u
summary(tk6u)
plot(tk6u,
     id.las=2,
     rl=FALSE)

## From: data.frame (dfm) - balanced
tk7 <- with(CRD2,
            TukeyC(x=dfm,
                   model='y ~ x',
                   which='x'))
tk7
summary(tk7)

## From: data.frame (dfm) - unbalanced
tk7u <- with(CRD2,
             TukeyC(x=dfm[-1,],
                    model='y ~ x',
                    which='x'))
tk7u
summary(tk7u)
plot(tk7u,
     id.las=2,
     rl=FALSE)

## From: aov - balanced
av2 <- with(CRD2,
            aov(y ~ x ,
                data=dfm))
summary(av2)

tk8 <- TukeyC(x=av2,
              which='x')
tk8
summary(tk8)

## From: aov - unbalanced
av2u <- with(CRD2,
             aov(y ~ x,
                 data=dfm[-1,]))
summary(av2u)

tk8u <- TukeyC(x=av2u,
               which='x')
tk8u
summary(tk8u)
plot(tk8u,
     id.las=2,
     rl=FALSE)
