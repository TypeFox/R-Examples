##
## Example: Split-split-plot Experiment (SSPE)
##

## Note: The factors are in uppercase and its levels in lowercase!

library(ScottKnott)
data(SSPE)

## From: design matrix (dm) and response variable (y)
## Main factor: P
sk1 <- with(SSPE,
            SK(dm,
               y,
               model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
               which='P',
               error='blk:P'))
summary(sk1)
plot(sk1)

# Main factor: SP
sk2 <- with(SSPE,
            SK(dm,
               y,
               model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
               which='SP',
               error='blk:P:SP',
               dispersion='s'))
summary(sk2)
plot(sk2,
     title='Main effect: SP')

# Main factor: SSP
sk3 <- with(SSPE,
            SK(dm,
               y,
               model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
               which='SSP', 
               error='Within',
               dispersion='se'))
summary(sk3)
plot(sk3,
     col=heat.colors(max(sk3$groups)),
     title='Main effect: SSP')

## Nested: p1/SP
skn1 <- with(SSPE,
             SK.nest(dm,
                     y,
                     model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
                     which='P:SP',
                     error='blk:P:SP',
                     fl1=1))
summary(skn1)
plot(skn1, 
     col='darkgray',
     title='Effect: p1/SP')


## From: data.frame
## Main factor: P
sk4 <- with(SSPE,
            SK(dfm,
               model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
               which='P',
               error='blk:P'))
summary(sk4)
plot(sk4,
     title='Main effect: P')

## Nested: p2/SP
skn2 <- with(SSPE,
             SK.nest(dfm,
                     model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
                     which='P:SP',
                     error='blk:P:SP',
                     fl1=2))
summary(skn2)
plot(skn2,
     title='Effect: p2/SP')

## Nested: p2/SP
skn3 <- with(SSPE,
             SK.nest(dfm,
                     model='y ~ blk + P*SP*SSP + Error(blk/P/SP)',
                     which='P:SP',
                     error='Within',
                     fl1=2))
summary(skn3)
plot(skn3,
     title='Effect: p2/SP')


## From: aovlist
av <- with(SSPE,
           aov(y ~  blk + P*SP*SSP + Error(blk/P/SP), 
               data=dfm))
summary(av)

## Main factor: P 
sk5 <- SK(av, 
          which='P',
          error='blk:P')
summary(sk5)
plot(sk5,
     title='Main effect: P')

## Main factor: SSP
sk6 <- SK(av,
          which='SSP',
          error='Within')
summary(sk6)
plot(sk6, 
     col=c('black', 'darkgray', 'gray'),
     title='Main effect = SSP')

## Nested: p1/SP
skn4 <- SK.nest(av,
                which='P:SP',
                error='blk:P:SP',
                fl1=1)
summary(skn4)
plot(skn4, 
     title='Effect: p1/SP')

## Nested: p2/SP
skn5 <- SK.nest(av,
                which='P:SP',
                error='blk:P:SP',
                fl1=2)
summary(skn5)
plot(skn5, title='Effect: p2/SP')

## Nested: P/SP/SSP (at various levels of SP and P)
skn6 <- SK.nest(av,
                which='P:SP:SSP',
                error='Within',
                fl1=1,
                fl2=1)
summary(skn6)
plot(skn6, 
     col=c('red', 'darkgreen'),
     title='Effect: p1/sp1/SSP')

skn7 <- SK.nest(av,
                which='P:SP:SSP', 
                error='Within',
                fl1=1,
                fl2=2)
summary(skn7)
plot(skn7, 
     col='darkgreen', 
     title='Effect: p1/sp2/SSP')

skn8 <- SK.nest(av,
                which='P:SP:SSP', 
                error='Within', 
                fl1=3, 
                fl2=3)
summary(skn8)
plot(skn8, 
     col='darkgreen',
     title='Effect: p3/sp3/SSP')

skn9 <- SK.nest(av,
                which='P:SP:SSP',
                error='Within',
                fl1=3, 
                fl2=2)
summary(skn9)
plot(skn9,
     id.lab=paste('test', 
                  1:length(skn9$groups),
                  sep='_'),
     title='Effect: p3/sp2/SSP')
