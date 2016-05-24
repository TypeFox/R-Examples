require(lmerTest)

## TODO uncomment the following whenever fixed deepcopy thing
## with the lme4
ifTest <- TRUE

if(ifTest){
  m <- lmer(Coloursaturation ~ TVset*Picture+
              (1|Assessor)+(1|Assessor:TVset), data=TVbo)
  
  f1 <- function(x){
    m1 <- refit(object=m, newresp=x, 
                rename.response = TRUE)
    m1
  }
  
  f2 <- function(x){
    m2 <- as(refit(object=m, newresp=x, 
                   rename.response = TRUE), "merModLmerTest")
    m2
  }
  
  f4 <- function(x){
    
    assign("x", x, env=environment(formula(m)))
    rf <- refit(object=m, newresp=x, 
                rename.response = TRUE)
    step(rf)
  }
  
  tools::assertError(update(f1(TVbo$Colourbalance)))
  tools::assertError(step(f2(TVbo$Colourbalance)))
  
  # the following does not work for lme4 1.1-8
  #lapply(TVbo[, 7, drop=FALSE], f4)
  
  ## after the assignment the error disappears
  ## Why?! seems like x becomes attached... 
  ## is it OK to do like that within a package?
  
  ## the following does not work for lme4 1.1-8
  ## update(f1(TVbo$Colourbalance))
  ## step(f2(TVbo$Colourbalance))
  
  
  
  ## example from R-sig mixed by Ben
  fit <- lm(sr ~ ., data = LifeCycleSavings)
  anova(fit)
  
  
  ## construct lmer model with near-zero variance
  LC2 <- transform(LifeCycleSavings,f=factor(1:2)) ## bogus
  
  form <- sr ~ pop15 + pop75 + dpi + ddpi + (1|f)  ## hack
  ## to avoid (Error in terms.formula(formula(x, fixed.only = TRUE)) :
  ##   '.' in formula and no 'data' argument)
  
  lmod <- lFormula(form, data=LC2)
  d2 <- lmer(form, data=LC2,devFunOnly=TRUE)
  llik <- d2(1e-5)
  fit2 <- mkMerMod(environment(d2),opt=list(par=1e-5,
                                            fval=llik,
                                            feval=1,
                                            conv=0,
                                            message=NULL),
                   lmod$reTrms, fr = lmod$fr)
  all.equal(coef(fit),fixef(fit2))
  anova(fit2)   ## practically equal to anova(fit) above
  
  ## anova does not work in this case
  an1 <- anova(as(fit2, "merModLmerTest"))
  
  tools::assertError(stopifnot(ncol(an1) == 6))
  
  ## anova works in the following case
  d2 <- lmer(form, data=LC2)
  an2 <- anova(d2)
  stopifnot(ncol(an2) == 6)  
}



