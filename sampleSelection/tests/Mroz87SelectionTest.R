library( "sampleSelection" )
library( "lmtest" )
data( "Mroz87" )
options( digits = 3 )

## Greene( 2003 ): example 22.8, page 786
## 2-step estimation
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
greene <- heckit( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 )
print( greene )
print(summary( greene ))
print(summary( greene$lm ) )
print( summary( greene$probit ) )
print( greene$sigma )
print( greene$rho )
print( coef( greene ) )
print( coef( greene, part = "outcome" ) )
print( coef( summary( greene ) ) )
print( coef( summary( greene ), part = "outcome" ) )
print( vcov( greene ) )
print( vcov( greene, part = "outcome" ) )
nobs( greene )
nObs( greene )
fitted( greene, part = "outcome" )
fitted( greene, part = "selection" )
residuals( greene, part = "outcome" )
all.equal( residuals( greene ), residuals( greene, part = "outcome" ) )
residuals( greene, part = "selection" )
all.equal( residuals( greene, part = "selection" ),
   residuals( greene, part = "selection", type = "deviance" ) )
all.equal( residuals( greene$probit ),
   residuals( greene, part = "selection" ) )
residuals( greene, part = "selection", type = "pearson" )
all.equal( residuals( greene$probit, type = "pearson" ),
   residuals( greene, part = "selection", type = "pearson" ) )
residuals( greene, part = "selection", type = "response" )
all.equal( residuals( greene$probit, type = "response" ),
   residuals( greene, part = "selection", type = "response" ) )
all.equal( residuals( greene, part = "selection", type = "response" ),
   Mroz87$lfp - fitted( greene, part = "selection" ) )
print( model.matrix( greene, part = "outcome" ) )
print( model.matrix( greene, part = "selection" ) )
print( model.frame( greene ) )
try( logLik( greene ) )

# ML estimation
greeneMl <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 )
print( greeneMl )
print(summary( greeneMl ))
print( coef( greeneMl ) )
print( coef( greeneMl, part = "outcome" ) )
print( coef( summary( greeneMl ) ) )
print( coef( summary( greeneMl ), part = "outcome" ) )
print( vcov( greeneMl ) )
print( vcov( greeneMl, part = "outcome" ) )
all.equal( nobs( greene ), nobs( greeneMl ) )
all.equal( nObs( greene ), nObs( greeneMl ) )
fitted( greeneMl, part = "outcome" )
fitted( greeneMl, part = "selection" )
residuals( greeneMl, part = "outcome" )
all.equal( residuals( greeneMl ), residuals( greeneMl, part = "outcome" ) )
residuals( greeneMl, part = "selection" )
all.equal( residuals( greeneMl, part = "selection" ),
   residuals( greeneMl, part = "selection", type = "deviance" ) )
!isTRUE( all.equal( residuals( greene, part = "selection" ),
   residuals( greeneMl, part = "selection" ) ) )
residuals( greeneMl, part = "selection", type = "pearson" )
!isTRUE( all.equal( residuals( greene, part = "selection", type = "pearson" ),
   residuals( greeneMl, part = "selection", type = "pearson" ) ) )
residuals( greeneMl, part = "selection", type = "response" )
!isTRUE( all.equal( residuals( greene, part = "selection", type = "response" ),
   residuals( greeneMl, part = "selection", type = "response" ) ) )
all.equal( residuals( greeneMl, part = "selection", type = "response" ),
   Mroz87$lfp - fitted( greeneMl, part = "selection" ) )
all.equal( model.matrix( greene, part = "outcome" )[,-6],
   model.matrix( greeneMl, part = "outcome" ),
   check.attributes = FALSE )
model.matrix( greeneMl, part = "outcome" )
all.equal( model.matrix( greene, part = "selection" ),
   model.matrix( greeneMl, part = "selection" ),
   check.attributes = FALSE )
model.matrix( greeneMl, part = "selection" )
all.equal( model.frame( greene )[,-c(2,12)], model.frame( greeneMl ),
   check.attributes = FALSE )
model.frame( greeneMl )
logLik( greeneMl )

# LR tests
greeneMl1 <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ 1, Mroz87 )
lrtest( greeneMl, greeneMl1 )
greeneMl2 <- selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ educ, Mroz87 )
lrtest( greeneMl, greeneMl2 )
greeneMl3 <- selection( lfp ~ faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 )
lrtest( greeneMl, greeneMl3 )

## Wooldridge( 2003 ): example 17.5, page 590
## 2-step estimation
data( "Mroz87" )
wooldridge <- heckit( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
   kids5 + kids618, log( wage ) ~ educ + exper + I( exper^2 ), Mroz87 )
print( wooldridge )
print( summary( wooldridge ) )
print( summary( wooldridge$lm ) )
print( summary( wooldridge$probit ) )
print( wooldridge$sigma )
print( wooldridge$rho )
print( coef( wooldridge ) )
print( coef( wooldridge, part = "outcome" ) )
print( coef( summary( wooldridge ) ) )
print( coef( summary( wooldridge ), part = "outcome" ) )
print( vcov( wooldridge ) )
print( vcov( wooldridge, part = "outcome" ) )
nobs( wooldridge )
nObs( wooldridge )
round( fitted( wooldridge, part = "outcome" ), digits = 3 )
round( fitted( wooldridge, part = "selection" ), digits = 3 )
round( residuals( wooldridge, part = "outcome" ), digits = 3 )
all.equal( residuals( wooldridge ), residuals( wooldridge, part = "outcome" ) )
round( residuals( wooldridge, part = "selection" ), digits = 3 )
all.equal( residuals( wooldridge, part = "selection" ),
   residuals( wooldridge, part = "selection", type = "deviance" ) )
all.equal( residuals( wooldridge$probit ),
   residuals( wooldridge, part = "selection" ) )
round( residuals( wooldridge, part = "selection", type = "pearson" ), digits = 3 )
all.equal( residuals( wooldridge$probit, type = "pearson" ),
   residuals( wooldridge, part = "selection", type = "pearson" ) )
round( residuals( wooldridge, part = "selection", type = "response" ), digits = 3 )
all.equal( residuals( wooldridge$probit, type = "response" ),
   residuals( wooldridge, part = "selection", type = "response" ) )
all.equal( residuals( wooldridge, part = "selection", type = "response" ),
   Mroz87$lfp - fitted( wooldridge, part = "selection" ) )
print( model.matrix( wooldridge, part = "outcome" ) )
print( model.matrix( wooldridge, part = "selection" ) )
print( model.frame( wooldridge ) )
try( logLik( wooldridge ) )

# ML estimation
wooldridgeMl <- selection( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
      kids5 + kids618, log( wage ) ~ educ + exper + I( exper^2 ), Mroz87 )
print( wooldridgeMl )
print( summary( wooldridgeMl ) )
print( coef( wooldridgeMl ) )
print( coef( wooldridgeMl, part = "outcome" ) )
print( coef( summary( wooldridgeMl ) ) )
print( coef( summary( wooldridgeMl ), part = "outcome" ) )
print( vcov( wooldridgeMl ) )
print( vcov( wooldridgeMl, part = "outcome" ) )
nobs( wooldridgeMl )
nObs( wooldridgeMl )
round( fitted( wooldridgeMl, part = "outcome" ), digits = 3 )
round( fitted( wooldridgeMl, part = "selection" ), digits = 3 )
round( residuals( wooldridgeMl, part = "outcome" ), digits = 3 )
all.equal( residuals( wooldridgeMl ),
   residuals( wooldridgeMl, part = "outcome" ) )
round( residuals( wooldridgeMl, part = "selection" ), digits = 3 )
all.equal( residuals( wooldridgeMl, part = "selection" ),
   residuals( wooldridgeMl, part = "selection", type = "deviance" ) )
!isTRUE( all.equal( residuals( wooldridge, part = "selection" ),
   residuals( wooldridgeMl, part = "selection" ) ) )
round( residuals( wooldridgeMl, part = "selection", type = "pearson" ), digits = 3 )
!isTRUE( all.equal( residuals( wooldridge, part = "selection", type = "pearson" ),
   residuals( wooldridgeMl, part = "selection", type = "pearson" ) ) )
round( residuals( wooldridgeMl, part = "selection", type = "response" ), digits = 3 )
!isTRUE( all.equal( residuals( wooldridge, part = "selection", type = "response" ),
   residuals( wooldridgeMl, part = "selection", type = "response" ) ) )
all.equal( residuals( wooldridgeMl, part = "selection", type = "response" ),
   Mroz87$lfp - fitted( wooldridgeMl, part = "selection" ) )
all.equal( model.matrix( wooldridge, part = "outcome" )[,-5],
   model.matrix( wooldridgeMl, part = "outcome" ),
   check.attributes = FALSE )
model.matrix( wooldridgeMl, part = "outcome" )
all.equal( model.matrix( wooldridge, part = "selection" ),
   model.matrix( wooldridgeMl, part = "selection" ),
   check.attributes = FALSE )
model.matrix( wooldridgeMl, part = "selection" )
all.equal( model.frame( wooldridge )[,-c(2,11)], model.frame( wooldridgeMl ),
   check.attributes = FALSE )
print( model.frame( wooldridgeMl ) )
logLik( wooldridgeMl )

# LR tests
wooldridgeMl1 <- selection( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
      kids5 + kids618, log( wage ) ~ 1, Mroz87 )
lrtest( wooldridgeMl1, wooldridgeMl )
wooldridgeMl2 <- selection( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
      kids5 + kids618, log( wage ) ~ educ, Mroz87 )
lrtest( wooldridgeMl2, wooldridgeMl )
wooldridgeMl3 <- selection( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
      kids5, log( wage ) ~ educ + exper, Mroz87 )
lrtest( wooldridgeMl3, wooldridgeMl )
