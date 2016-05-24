##
## Example: Randomized Complete Block Design (RCBD)
##

## The parameters can be: design matrix and the response variable,
## data.frame or aov

library(ScottKnott)
data(RCBD)

## Design matrix (dm) and response variable (y)
sk1 <- with(RCBD,
            SK(x=dm,
               y=y,
               model='y ~ blk + tra',
               which='tra'))
summary(sk1)
plot(sk1)

## From: data.frame (dfm), which='tra'
sk2 <- with(RCBD,
            SK(x=dfm,
               model='y ~ blk + tra',
               which='tra',
               dispersion='s'))
summary(sk2)
plot(sk2,
     mm.lty=3,
     title='Factor levels')

## From: data.frame (dfm), which='blk'
sk3 <- with(RCBD,
            SK(x=dfm,
               model='y ~ blk + tra',
               which='blk'))
summary(sk3)
plot(sk3,
     id.lab=paste('Block',
                  1:length(sk3$groups),
                  sep='_'), 
     title='Blocks')

## From: aov
av1 <- with(RCBD,
            aov(y ~ blk + tra,
                data=dfm))
summary(av1)

## From: aov, which='blk' implicit
sk4 <- SK(x=av1)
summary(sk4)
plot(sk4,
     id.lab=paste('Block',
                  1:length(sk4$groups),
                  sep='_'),
     title='Blocks')

## From: aov, which='blk' explicit
sk5 <- SK(x=av1,
          which='blk')
summary(sk5)
plot(sk5,
     id.lab=paste('Block',
                  1:length(sk5$groups),
                  sep='_'), 
     title='Blocks')

## From: aov, which='tra' explicit
sk6 <- SK(x=av1,
          which='tra')
summary(sk6)
plot(sk6,
     xlab='',
     title='Factor levels')

