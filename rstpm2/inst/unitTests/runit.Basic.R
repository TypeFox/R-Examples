.setUp <- function() {
  suppressMessages( require(rstpm2) )
  data(brcancer)
  fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=5)
  fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                logH.formula=~nsx(log(rectime),df=3,stata=TRUE))
}

test.basic.1 <- function() {
  suppressMessages( require(rstpm2) )
  data(brcancer)
  fit1 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=5)
  checkEqualsNumeric(round(coef(fit1)[2],5), -0.36457)
}

test.basic.stata <- function() {
  suppressMessages( require(rstpm2) )
  data(brcancer)
  fit2 <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                logH.formula=~nsx(log(rectime),df=3,stata=TRUE))
  checkEqualsNumeric(round(coef(fit2)[2],5), -0.36144)
}



## .tearDown <- function() {
##   rm(fit1,fit2)
## }
