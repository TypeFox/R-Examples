##
## Example: Split-plot Experiment in time (SPET)
##

## The parameters can be: design matrix and the response variable,
## data.frame or aov

## Note: The factors are in uppercase and its levels in lowercase!

library(ScottKnott)
data(SPET)

## From: design matrix (dm) and response variable (y)
## Main factor: tra
sk1 <- with(SPET,
            SK(x=dm, 
               y=y, 
               model='y ~ blk + tra*year + Error(blk/tra)',
               which='tra',
               error='blk:tra'))
summary(sk1)

plot(sk1)

## Main factor: year
sk2 <- with(SPET,
            SK(x=dm,
               y=y,
               model='y ~ blk + tra*year + Error(blk/tra)',
               which='year',
               error='Within'))
summary(sk2)
plot(sk2,
     title='Main effect: year')

## Nested: crotjuncea/year
skn1 <- with(SPET,
             SK.nest(x=dm,
                     y=y,
                     model='y ~ blk + tra*year + Error(blk/tra)',
                     which='tra:year',
                     error='Within',
                     fl1=2))
summary(skn1)
plot(skn1,
     title='Effect: crotjuncea/year')


## From: data.frame
## Main factor: year
sk3 <- with(SPET,
            SK(dfm,
               model='y ~ blk + tra*year + Error(blk/tra)',
               which='year',
               error='Within'))
summary(sk3)
plot(sk3, 
     title='Main effect: year')

## Nested: crotjuncea/year
skn2 <- with(SPET,
             SK.nest(dfm,
                     model='y ~ blk + tra*year + Error(blk/tra)',
                     which='tra:year',
                     error ='Within',
                     fl1=2))
summary(skn2)
plot(skn2,
     title='Effect: crotjuncea/year')

## Nested: year_1/tra
skn3 <- with(SPET,
             SK.nest(dfm,
                     model='y ~ blk + tra*year + Error(blk/tra)',
                     which='year:tra',
                     error ='Within',
                     fl1=1))
summary(skn3)
plot(skn3,
     title='Effect: year_1/tra')



## From: aovlist
av1 <- with(SPET,
            aov(y ~ blk + tra*year + Error(blk/tra),
                data=dfm))
summary(av1)

## Main factor: year
sk4 <- SK(av1,
          which='year',
          error='Within')
summary(sk4)
plot(sk4, title='Main effect: year')

## Main factor: tra
sk5 <- SK(av1,
          which='tra',
          error='blk:tra')
summary(sk5)
plot(sk5, 
     title='Main effect: tra')

## Nested: crotjuncea/year
skn4 <- SK.nest(av1,
                which='tra:year',
                error='Within',
                fl1=2)
summary(skn4)
plot(skn4,
     title='Effect: crotjuncea/year')

## Nested: guandu/year
skn5 <- SK.nest(av1,
                which='tra:year',
                error='Within',
                fl1=4)
summary(skn5)
plot(skn5,
     title='Effect: guandu/year')


## Nested: year_1/tra
skn6 <- SK.nest(av1,
                which='year:tra',
                error='Within',
                fl1=1)
summary(skn6)
plot(skn6,
     title='Effect: year_1/tra')

## Nested: year_2/tra
skn7 <- SK.nest(av1,
                which='year:tra',
                error='Within',
                fl1=2)
summary(skn7)
plot(skn7,
     title='Effect: year_2/tra')

