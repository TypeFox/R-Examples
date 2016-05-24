##
## Example: Latin Squares Design (LSD)
##

## The parameters can be: design matrix and the response variable,
## data.frame or aov

library(TukeyC)
data(LSD)

## From: design matrix (dm) and response variable (y)
tk1 <- with(LSD,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ rows + cols + tra',
                   which='tra'))
summary(tk1)  # p-value E-A = .051
plot(tk1)

tk1 <- with(LSD,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ rows + cols + tra',
                   which='tra',
                   sig.level=.052))
summary(tk1)
plot(tk1)

## From: data.frame
tk2 <- with(LSD,
            TukeyC(x=dfm,
                   model='y ~ rows + cols + tra',
                   which='tra',
                   sig.level=.052))
summary(tk2)

## From: aov
av1 <- with(LSD,
            aov(y ~ rows + cols + tra,
                data=dfm))
summary(av1)

tk3 <- TukeyC(av1,
              which='tra',
              sig.level=.052)
summary(tk3)
