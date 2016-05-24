##
## Example: Factorial Experiment (FE)
##

## The parameters can be: design matrix and the response variable,
## data.frame or aov

## Note: Upper case for factors and lowercase for levels

library(TukeyC)
data(FE)

## From: design matrix (dm) and response variable (y)
## Main factor: N
tk1 <- with(FE,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ blk + N*P*K', 
                   which='N'))
summary(tk1)
plot(tk1)

## Main factor: P
tk2 <- with(FE,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ blk + N*P*K',
                   which='P',
                   dispersion='s'))
summary(tk2)
plot(tk2)

## Main factor: K
tk3 <- with(FE,
            TukeyC(x=dm,
                   y=y,
                   model='y ~ blk + N*P*K',
                   which='K',
                   dispersion='se'))
summary(tk3)
plot(tk3)

## Nested: p1/N
## Testing N inside of level one of P
ntk1 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + N*P*K',
                         which='P:N',
                         fl1=1))
summary(ntk1)
plot(ntk1)

## Nested: p2/N
ntk2 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + N*P*K',
                         which='P:N',
                         fl1=2))
summary(ntk2)
plot(ntk2)

## Nested: k1/N
ntk3 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + N*P*K',
                         which='K:N',
                         fl1=1))
summary(ntk3)

## Nested: k2/N
ntk4 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + N*P*K',
                         which='K:N',
                         fl1=2))
summary(ntk4)

## Nested: k1/P
ntk5 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + N*P*K',
                         which='K:P',
                         fl1=1))
summary(ntk5)

## Nested: k2/P
ntk6 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + N*P*K',
                         which='K:P',
                         fl1=2))
summary(ntk6)

## Nested: k1/p1/N
## Testing N inside of level one of K and level one of P
ntk7 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + N*P*K',
                         which='K:P:N',
                         fl1=1,
                         fl2=1))
summary(ntk7)   # p-value=.053 to be different
plot(ntk7)

ntk7 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + N*P*K',
                         which='K:P:N',
                         fl1=1,
                         fl2=1,
                         sig.level=.054))
summary(ntk7)   # The p-value is very informative!
plot(ntk7)

## Nested: k2/p2/N
ntk8 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + N*P*K',
                         which='K:P:N',
                         fl1=2,
                         fl2=2))
summary(ntk8)
plot(ntk8)

## Nested: k1/n1/P
ntk9 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + P*N*K',
                         which='K:N:P',
                         fl1=1,
                         fl2=1))
summary(ntk9)

## Nested: k2/n2/P
ntk8 <- with(FE,
             TukeyC.nest(x=dm,
                         y=y,
                         model='y ~ blk + P*N*K',
                         which='K:N:P',
                         fl1=2,
                         fl2=2))
summary(ntk8)

## Nested: p1/n1/K
ntk10 <- with(FE,
              TukeyC.nest(x=dm,
                          y=y,
                          model='y ~ blk + K*N*P',
                          which='P:N:K',
                          fl1=1,
                          fl2=1))
summary(ntk10)

## Nested: p2/n2/K
ntk11 <- with(FE,
              TukeyC.nest(x=dm,
                          y=y,
                          model='y ~ blk + K*N*P',
                          which='P:N:K',
                          fl1=2,
                          fl2=2))
summary(ntk11)


## From: data.frame
## Nested: k2/p1/N
ntk12 <- with(FE,
              TukeyC.nest(x=dfm,
                          model='y ~ blk + N*P*K',
                          which='K:P:N',
                          fl1=1,
                          fl2=2))
summary(ntk12)

## Nested: k1/p2/N
ntk13 <- with(FE,
              TukeyC.nest(x=dfm,
                          model='y ~ blk + N*P*K',
                          which='K:P:N',
                          fl1=2,
                          fl2=1))
summary(ntk13)


## From aov
nav1 <- with(FE,
             aov(y ~ blk + N*P*K ,
                 data=dfm))
summary(nav1)

## Main factor: N
ntk14 <- TukeyC(nav1,
                which='N')
summary(ntk14)

## Nested: k1/P
ntk15 <- TukeyC.nest(nav1,
                     which='K:P',
                     fl1=1)
summary(ntk15)

## Nested: k2/p1/N
ntk16 <- TukeyC.nest(nav1,
                     which='K:P:N',
                     fl1=1,
                     fl2=2)
summary(ntk16)

# Changing the order of factors (test)
nav2 <- with(FE,
             aov(y ~ blk + K*N* P,
                 data=dfm))
summary(nav2)

## Nested: p1/n1/K
ntk17 <- TukeyC.nest(nav2,
                     which='P:N:K',
                     fl1=1,
                     fl2=1)
summary(ntk17)
