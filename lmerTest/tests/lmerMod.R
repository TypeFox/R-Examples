require(lmerTest)

## rand step lsmeans difflsmeans should work only with inheritance of lmerMod class

## TODO uncomment the following whenever fixed deepcopy thing
## with the lme4
ifTest <- TRUE

if(ifTest){
  gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               data = cbpp, family = binomial)
  
  ## should not work with glmer models
  tools::assertError(rand(gm1))
  tools::assertError(step(gm1))
  tools::assertError(lsmeans(gm1))
  tools::assertError(difflsmeans(gm1))
  
  ## should not work with nlmer models
  # startvec <- c(Asym = 200, xmid = 725, scal = 350)
  # nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
  #               Orange, start = startvec)
  # tools::assertError(rand(nm1))
  # tools::assertError(step(nm1))
  # tools::assertError(lsmeans(nm1))
  # tools::assertError(difflsmeans(nm1))
  
  
  ## should wotk with lmer from lme4 package (class lmerMod)
  m <- lme4::lmer(Coloursaturation ~ TVset*Picture+
                    (1|Assessor), data=TVbo)
  rand(m)
  step(m)
  lsmeans(m)
  difflsmeans(m)
}


