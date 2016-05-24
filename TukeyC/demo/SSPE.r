##
## Example: Split-split-plot Experiment (SSPE)
##

## Note: Upper case for factors and lowercase for levels

library(TukeyC)
data(SSPE)

## From: design matrix (dm) and response variable (y)
## Main factor: P
tk1 <- with(SSPE,
            TukeyC(dm,
                   y,
                   model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
                   which='P',
                   error='blk:P'))
summary(tk1)
plot(tk1)

# Main factor: SP
tk2 <- with(SSPE,
            TukeyC(dm,
                   y,
                   model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
                   which='SP',
                   error='blk:P:SP'))
summary(tk2)
plot(tk2)

# Main factor: SSP
tk3 <- with(SSPE, 
            TukeyC(dm,
                   y,
                   model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
                   which='SSP',
                   error='Within'))
summary(tk3)
plot(tk3)

## Nested: p1/SP
tkn1 <- with(SSPE,
             TukeyC.nest(dm,
                         y,
                         model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
                         which='P:SP',
                         error='blk:P:SP',
                         fl1=1))
summary(tkn1)


## From: data.frame
## Main factor: P
tk4 <- with(SSPE,
            TukeyC(dfm,
                   model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
                   which='P',
                   error='blk:P'))
summary(tk4)

## Nested: p2/SP
tkn2 <- with(SSPE,
             TukeyC.nest(dfm,
                         model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
                         which='P:SP',
                         error='blk:P:SP',
                         fl1=2))
summary(tkn2)

## Nested: p2/SSP
tkn3 <- with(SSPE,
             TukeyC.nest(dfm,
                         model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
                         which='P:SSP',
                         error='Within',
                         fl1=2))
summary(tkn3)
plot(tkn3)


## From: aovlist
av <- with(SSPE,
           aov(y ~  blk + P*SP*SSP + Error(blk/P/SP),
               data=dfm))
summary(av)

## Main factor: P 
tk5 <- TukeyC(av,
              which='P',
              error='blk:P')
summary(tk5)

## Main factor: SSP
tk6 <- TukeyC(av,
              which='SSP',
              error='Within')
summary(tk6)

## Nested: p1/SP
tkn4 <- TukeyC.nest(av,
                    which='P:SP',
                    error='blk:P:SP',
                    fl1=1)
summary(tkn4)

## Nested: p2/SP
tkn5 <- TukeyC.nest(av,
                    which='P:SP',
                    error='blk:P:SP', 
                    fl1=2)
summary(tkn5)

## Nested: Pi/SPi/SSP (at various levels of P and SP)
tkn6 <- TukeyC.nest(av,
                    which='P:SP:SSP',
                    error='Within',
                    fl1=1,
                    fl2=1)
summary(tkn6)
plot(tkn6)

tkn7 <- TukeyC.nest(av,
                    which='P:SP:SSP',
                    error='Within',
                    fl1=2,
                    fl2=1)
summary(tkn7)

tkn8 <- TukeyC.nest(av,
                    which='P:SP:SSP',
                    error='Within',
                    fl1=3,
                    fl2=3)
summary(tkn8)

tkn9 <- TukeyC.nest(av,
                    which='P:SP:SSP',
                    error='Within',
                    fl1=2,
                    fl2=3)
summary(tkn9)
