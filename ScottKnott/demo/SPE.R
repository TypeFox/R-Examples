##
## Example: Split-plot Experiment (SPE)
##

## The parameters can be: design matrix and the response variable,
## data.frame or aov

## Note: The factors are in uppercase and its levels in lowercase!

# library(ScottKnott)
data(SPE)

## From: design matrix (dm) and response variable (y)
## Main factor: P
sk1 <- with(SPE,
            SK(x=dm, 
               y=y, 
               model='y ~ blk + P*SP + Error(blk/P)',
               which='P',
               error='blk:P'))
summary(sk1)
plot(sk1)

## Main factor: SP
sk2 <- with(SPE,
            SK(x=dm,
               y=y,
               model='y ~ blk + P*SP + Error(blk/P)',
               which='SP',
               error='Within',
               dispersion='s'))
summary(sk2)
plot(sk2,
     title='Main effect: SP')

## Nested: p1/SP
## Studing SP inside the level one of P
skn1 <- with(SPE,
             SK.nest(x=dm,
                     y=y,
                     model='y ~ blk + P*SP + Error(blk/P)',
                     which='P:SP',
                     error='Within',
                     fl1=1,
                     dispersion='se'))
summary(skn1)
plot(skn1,
     title='Effect: p1/SP')


## From: data.frame
## Main factor: SP
sk3 <- with(SPE,
            SK(dfm,
               model='y ~ blk + P*SP + Error(blk/P)',
               which='SP',
               error='Within'))
summary(sk3)
plot(sk3, 
     title='Main effect: SP')

## Nested: p1/SP
skn2 <- with(SPE,
             SK.nest(dfm,
                     model='y ~ blk + P*SP + Error(blk/P)',
                     which='P:SP',
                     error ='Within',
                     fl1=1 ))
summary(skn2)
plot(skn2,
     title='Effect: p1/SP')


## From: aovlist
av1 <- with(SPE,
            aov(y ~ blk + P*SP + Error(blk/P),
                data=dfm))
summary(av1)

## Main factor: SP
sk4 <- SK(av1,
          which='SP',
          error='Within')
summary(sk4)
plot(sk4, title='Main effect: SP')

## Main factor: P
sk5 <- SK(av1,
          which='P',
          error='blk:P')
summary(sk5)
plot(sk5, 
     title='Main effect: P')

## Nested: p1/SP
## Studing SP inside of level one of P
skn3 <- SK.nest(av1,
                which='P:SP',
                error='Within',
                fl1=1)
summary(skn3)
plot(skn3,
     title='Effect: p1/SP')

## Nested: p2/SP
skn4 <- SK.nest(av1,
                which='P:SP',
                error='Within',
                fl1=2)
summary(skn4)
plot(skn4,
     title='Effect: p2/SP')

## Nested: p3/SP
skn5 <- SK.nest(av1,
                which='P:SP',
                error='Within',
                fl1=3)
summary(skn5)
plot(skn5,
     title='Effect: p3/SP')
