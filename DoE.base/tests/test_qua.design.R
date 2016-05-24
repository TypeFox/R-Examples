require(DoE.base)

set.seed(12445)

## NA
plan1 <- fac.design(4,nlevels=2)
 lm(rnorm(16)~., plan1)

plan2 <- oa.design(nlevels=c(2,6,2))
 lm(rnorm(12)~., plan2)

### NA to quantitative
quantplan1 <- qua.design(plan1, quantitative="all")
 lm(rnorm(16)~., quantplan1)

quantplan2 <- qua.design(plan2, quantitative="all")
 lm(rnorm(12)~., quantplan2)
 
### NA to qualitative, no contrasts given
### (does not change anything)
qualplan1 <- qua.design(plan1, quantitative="none")
 lm(rnorm(16)~., qualplan1)

qualplan2 <- qua.design(plan2, quantitative="none")
 lm(rnorm(12)~., qualplan2)

### quantitative to qualitative, no contrasts given
qualplan1 <- qua.design(quantplan1, quantitative="none")
 lm(rnorm(16)~., qualplan1)

qualplan2 <- qua.design(quantplan2, quantitative="none")
 lm(rnorm(12)~., qualplan2)

### quantitative to NA, no contrasts given
qualplan1 <- qua.design(quantplan1, quantitative=NA)
 lm(rnorm(16)~., qualplan1)

qualplan2 <- qua.design(quantplan2, quantitative=NA)
 lm(rnorm(12)~., qualplan2)

### quantitative to NA, contrasts given
### contrasts are ignored
qualplan1 <- qua.design(quantplan1, quantitative=NA, contrasts=c(B="contr.treatment"))
 lm(rnorm(16)~., qualplan1)

qualplan2 <- qua.design(quantplan2, quantitative=NA, contrasts=c(B="contr.treatment"))
 lm(rnorm(12)~., qualplan2)

### NA to qualitative, contrasts given
qualplan1 <- qua.design(plan1, quantitative="none", contrasts=c(B="contr.treatment"))
 lm(rnorm(16)~., qualplan1)

qualplan2 <- qua.design(plan2, quantitative="none", contrasts=c(B="contr.treatment"))
 lm(rnorm(12)~., qualplan2)

### quantitative to qualitative, contrasts given
qualplan1 <- qua.design(quantplan1, quantitative="none", contrasts=c(B="contr.treatment"))
 lm(rnorm(16)~., qualplan1)

qualplan2 <- qua.design(quantplan2, quantitative="none", contrasts=c(B="contr.treatment"))
 lm(rnorm(12)~., qualplan2)

plan3 <- oa.design(factor.names=list(X=c(1,2),Y=c(1,2),Z=c(1,2)), ID=L4.2.3)
desnum(qua.design(qua.design(cross.design(plan1,plan3),quantitative="all")))

desnum(change.contr(quantplan1, "contr.helmert"))
desnum(change.contr(quantplan2, "contr.helmert"))
