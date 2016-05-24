##
## Example: Split-plot Experiment (SPE)
##

## The parameters can be: design matrix and the response variable,
## data.frame or aov

## Note: Upper case for factors and lowercase for levels

library(TukeyC)
data(SPE)

## From: design matrix (dm) and response variable (y)
## Main factor: P
tk1 <- with(SPE,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ blk + P*SP + Error(blk/P)',
                   which='P',
                   error='blk:P'))
summary(tk1)
plot(tk1)

## Main factor: SP
tk2 <- with(SPE,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ blk + P*SP + Error(blk/P)',
                   which='SP',
                   error ='Within'))
summary(tk2)
plot(tk2)

## Nested: p1/SP
## Testing SP inside of level one of P
tkn1 <- with(SPE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + P*SP + Error(blk/P)',
                         which='P:SP',
                         error='Within',
                         fl1=1))
summary(tkn1)


## From: data.frame
## Main factor: SP
tk3 <- with(SPE,
            TukeyC(dfm,
                   model='y ~ blk + P*SP + Error(blk/P)',
                   which='SP',
                   error='Within',
                   sig.level=0.025))
summary(tk3)

## Nested: p1/SP
tkn2 <- with(SPE,
             TukeyC.nest(dfm,
                         model='y ~ blk + P*SP + Error(blk/P)',
                         which='P:SP',
                         error='Within', 
                         fl1=1 ))
summary(tkn2)
plot(tkn2)

## From: aovlist
av1 <- with(SPE,
            aov(y ~ blk + P*SP + Error(blk/P),
                data=dfm))
summary(av1)

## Main factor: SP
tk4 <- TukeyC(av1,
              which='SP',
              error='Within',
              sig.level=0.1)
summary(tk4)

## Main factor: P
tk5 <- TukeyC(av1,
              which='P',
              error='blk:P')

summary(tk5)

## Nested: p1/SP
## Testing SP inside of level one of P
tkn3 <- TukeyC.nest(av1,
                    which='P:SP',
                    error='Within',
                    fl1=1)
summary(tkn3)

## Nested: p2/SP
tkn4 <- TukeyC.nest(av1,
                    which='P:SP',
                    error='Within',
                    fl1=2)
summary(tkn4)

## Nested: p3/SP
tkn5 <- TukeyC.nest(av1,
                    which='P:SP',
                    error='Within',
                    fl1=3)
summary(tkn5)
