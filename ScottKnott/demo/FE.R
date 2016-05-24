##
## Example: Factorial Experiment (FE)
##

## The parameters can be: design matrix and the response variable,
## data.frame or aov

## Note: The factors are in uppercase and its levels in lowercase!

library(ScottKnott)
data(FE)

## From: design matrix (dm) and response variable (y)
## Main factor: N
sk1 <- with(FE,
            SK(x=dm,
               y=y,
               model='y ~ blk + N*P*K', 
               which='N'))
summary(sk1)
plot(sk1,
     title='Main effect: N')

## Main factor: P
sk2 <- with(FE,
            SK(x=dm,
               y=y,
               model='y ~ blk + N*P*K',
               which='P',
               dispersion='s'))
summary(sk2)
plot(sk2,
     title='Main effect: P')

## Main factor: K
sk3 <- with(FE,
            SK(x=dm,
               y=y,
               model='y ~ blk + N*P*K',
               which='K',
               dispersion='se'))
summary(sk3)
plot(sk3, 
     title='Main effect: K')

## Nested: p1/N
## Studing N inside of the level one of P
nsk1 <- with(FE,
             SK.nest(x=dm,
                     y=y,
                     model='y ~ blk + N*P*K',
                     which='P:N', 
                     fl1=1,
                     dis='s'))
summary(nsk1)
plot(nsk1,
     title='Effect: p1/N')

## Nested: p2/N
nsk2 <- with(FE,
             SK.nest(x=dm,
                     y=y,
                     model='y ~ blk + N*P*K',
                     which='P:N',
                     fl1=2))
summary(nsk2)
plot(nsk2,
     title='Effect: p2/N')

## Nested: k1/N
nsk3 <- with(FE,
             SK.nest(x=dm,
                     y=y,
                     model='y ~ blk + N*P*K',
                     which='K:N',
                     fl1=1))
summary(nsk3)
plot(nsk3,
     title='Effect: k1/N')

## Nested: k2/N
nsk4 <- with(FE,
             SK.nest(x=dm,
                     y=y,
                     model='y ~ blk + N*P*K',
                     which='K:N',
                     fl1=2))
summary(nsk4)
plot(nsk4,
     title='Effect: k2/N')

## Nested: k1/P
nsk5 <- with(FE,
             SK.nest(x=dm,
                     y=y,
                     model='y ~ blk + N*P*K',
                     which='K:P',
                     fl1=1))
summary(nsk5)
plot(nsk5,
     title='Effect: k1/P')

## Nested: k2/P
nsk6 <- with(FE,
             SK.nest(x=dm, 
                     y=y,
                     model='y ~ blk + N*P*K',
                     which='K:P',
                     fl1=2))
summary(nsk6)
plot(nsk6,
     title='Effect: k2/P')

## Nested: k1/p1/N 
## Studing N inside of the level one of K and level one of P
nsk7 <- with(FE,
             SK.nest(x=dm,
                     y=y,
                     model='y ~ blk + N*P*K',
                     which='K:P:N',
                     fl1=1,
                     fl2=1))
summary(nsk7)
plot(nsk7,
     title='Effect: k1/p1/N')

## Nested: k2/p2/N
nsk8 <- with(FE,
             SK.nest(x=dm,
                     y=y,
                     model='y ~ blk + N*P*K',
                     which='K:P:N',
                     fl1=2,
                     fl2=2))
summary(nsk8)
plot(nsk8,
     title='Effect: k2/p2/N')

## Nested: k1/n1/P
nsk9 <- with(FE,
             SK.nest(x=dm,
                     y=y,
                     model='y ~ blk + N*P*K',
                     which='K:N:P',
                     fl1=1,
                     fl2=1))
summary(nsk9)
plot(nsk9,
     title='Effect: k1/n1/P')

## Nested: k2/n2/P
nsk10 <- with(FE,
              SK.nest(x=dm,
                      y=y,
                      model='y ~ blk + N*P*K',
                      which='K:N:P',
                      fl1=2,
                      fl2=2))
summary(nsk10)
plot(nsk10,
     title='Effect: k2/n2/P')

## Nested: p1/n1/K
nsk11 <- with(FE,
              SK.nest(x=dm,
                      y=y,
                      model='y ~ blk + N*P*K',
                      which='P:N:K',
                      fl1=1,
                      fl2=1))
summary(nsk11)
plot(nsk11, 
     title='Effect: p1/n1/K')

## Nested: p2/n2/K
nsk12 <- with(FE,
              SK.nest(x=dm,
                      y=y,
                      model='y ~ blk + N*P*K',
                      which='P:N:K',
                      fl1=2,
                      fl2=2))
summary(nsk12)
plot(nsk12,
     title='Effect: p2/n2/K')


## From: data.frame
## Nested: k2/p1/N
nsk13 <- with(FE,
              SK.nest(x=dm,
                      y=y, 
                      model='y ~ blk + N*P*K',
                      which='K:P:N',
                      fl1=2,
                      fl2=1))
summary(nsk13)
plot(nsk13,
     title='Effect: k2/p1/N')

## Nested: k1/p2/N
nsk14 <- with(FE,
              SK.nest(x=dm,
                      y=y,
                      model='y ~ blk + N*P*K',
                      which='K:P:N',
                      fl1=1,
                      fl2=2))
summary(nsk14)
plot(nsk14,
     title='Effect: k1/p2/N')


## From aov
nav1 <- with(FE, 
             aov(y ~ blk + N*P*K,
                 data=dfm))
summary(nav1)

## Main factor: N
nsk15 <- SK(nav1,
            which='N')
summary(nsk15)
plot(nsk15,
     title='Main effect: N')

## Nested: k1/P
nsk16 <- SK.nest(nav1,
                 which='K:P',
                 fl1=1)
summary(nsk15)
plot(nsk16, title='Effect: k1/P')

## Nested: k2/p1/N
nsk17 <- SK.nest(nav1,
                 which='K:P:N',
                 fl1=2,
                 fl2=1)
summary(nsk17)
plot(nsk17,
     title='Effect: k2/p1/N')

