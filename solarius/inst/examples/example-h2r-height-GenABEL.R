gkin <- ibs(ge03d2ex.clean)

gkin[upper.tri(gkin)] <- gkin[lower.tri(gkin)]
diag(gkin) <- 1
gkin <- gkin/2
#a <- ibs(ge03d2ex.clean,snps=sample(autosomal(ge03d2ex.clean),1000,replace=FALSE))
#a[upper.tri(a)] <- gkin[lower.tri(a)]
#diag(a) <- 1

plotKinship2(2*gkin, "hist")

df <- ge03d2ex.clean@phdata
df <- mutate(df,
  sex = sex + 1)
  
mod <- solarPolygenic(height ~ sex + age, df, kinship = gkin)
summary(mod)
