##
## Example: Randomized Complete Block Design (RCBD)
##

## The parameters can be: design matrix and the response variable,
## data.frame or aov

library(TukeyC)
data(RCBD)

## Design matrix (dm) and response variable (y)
tk1 <- with(RCBD,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ blk + tra',
                   which='tra'))
summary(tk1)
plot(tk1)

## From: data.frame (dfm), which='tra'
tk2 <- with(RCBD,
            TukeyC(x=dfm,
                   model='y ~ blk + tra',
                   which='tra'))
summary(tk2)

## From: data.frame (dfm), which='blk'
tk3 <- with(RCBD,
            TukeyC(x=dfm,
                   model='y ~ blk + tra',
                   which='blk'))
summary(tk3)
plot(tk3)

## From: aov
av1 <- with(RCBD,
            aov(y ~ blk + tra,
                data=dfm))
summary(av1)

## From: aov, which='blk' implicit
tk4 <- TukeyC(x=av1)
summary(tk4)

## From: aov, which='blk' explicit
tk5 <- TukeyC(x=av1,
              which='blk')
summary(tk5)

## From: aov, which='tra' explicit
tk6 <- TukeyC(x=av1,
              which='tra')
summary(tk6)
