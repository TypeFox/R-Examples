##
## Example: Split-plot Experiment in time (SPET)
##

## The parameters can be: design matrix and the response variable,
## data.frame or aov

## Note: The factors are in uppercase and its levels in lowercase!

library(TukeyC)
data(SPET)

## From: design matrix (dm) and response variable (y)
## Main factor: tra
tk1 <- with(SPET,
            TukeyC(x=dm, 
                   y=y, 
                   model='y ~ blk + tra*year + Error(blk/tra)',
                   which='tra',
                   error='blk:tra'))
summary(tk1)

plot(tk1)

## Main factor: year
tk2 <- with(SPET,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ blk + tra*year + Error(blk/tra)',
                   which='year',
                   error='Within'))
summary(tk2)
plot(tk2,
     title='Main effect: year')

## Nested: crotjuncea/year
tkn1 <- with(SPET,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + tra*year + Error(blk/tra)',
                         which='tra:year',
                         error='Within',
                         fl1=2))
summary(tkn1)
plot(tkn1,
     title='Effect: crotjuncea/year')


## From: data.frame
## Main factor: year
tk3 <- with(SPET,
            TukeyC(dfm,
                   model='y ~ blk + tra*year + Error(blk/tra)',
                   which='year',
                   error='Within'))
summary(tk3)
plot(tk3, 
     title='Main effect: year')

## Nested: crotjuncea/year
tkn2 <- with(SPET,
             TukeyC.nest(dfm,
                         model='y ~ blk + tra*year + Error(blk/tra)',
                         which='tra:year',
                         error ='Within',
                         fl1=2))
summary(tkn2)
plot(tkn2,
     title='Effect: crotjuncea/year')

## Nested: year_1/tra
tkn3 <- with(SPET,
             TukeyC.nest(dfm,
                         model='y ~ blk + tra*year + Error(blk/tra)',
                         which='year:tra',
                         error ='Within',
                         fl1=1))
summary(tkn3)
plot(tkn3,
     title='Effect: year_1/tra')



## From: aovlist
av1 <- with(SPET,
            aov(y ~ blk + tra*year + Error(blk/tra),
                data=dfm))
summary(av1)

## Main factor: year
tk4 <- TukeyC(av1,
              which='year',
              error='Within')
summary(tk4)
plot(tk4, title='Main effect: year')

## Main factor: tra
tk5 <- TukeyC(av1,
              which='tra',
              error='blk:tra')
summary(tk5)
plot(tk5, 
     title='Main effect: tra')

## Nested: crotjuncea/year
tkn4 <- TukeyC.nest(av1,
                    which='tra:year',
                    error='Within',
                    fl1=2)
summary(tkn4)
plot(tkn4,
     title='Effect: crotjuncea/year')

## Nested: guandu/year
tkn5 <- TukeyC.nest(av1,
                    which='tra:year',
                    error='Within',
                    fl1=4)
summary(tkn5)
plot(tkn5,
     title='Effect: guandu/year')


## Nested: year_1/tra
tkn6 <- TukeyC.nest(av1,
                    which='year:tra',
                    error='Within',
                    fl1=1)
summary(tkn6)
plot(tkn6,
     title='Effect: year_1/tra')

## Nested: year_2/tra
tkn7 <- TukeyC.nest(av1,
                    which='year:tra',
                    error='Within',
                    fl1=2)
summary(tkn7)
plot(tkn7,
     title='Effect: year_2/tra')

