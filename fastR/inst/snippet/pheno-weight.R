pheno.lm <- lm(log(weight) ~ log(waist) + log(height), pheno)
###hop:3-9
summary(pheno.lm)
