TukeyPvalue <-
function(hfobj)
{
additive.out <- lm(y~rows+cols,data=hfobj$tall)
psq <- fitted(additive.out)^2
singledf.out <- lm(y~rows + cols + psq,data=hfobj$tall)
singledf.anova <- anova(singledf.out)
#print(singledf.anova)
#pvalue <- prettyNum(singledf.anova$P[3],digits=4)
pvalue <- singledf.anova$P[3]
partial.r2 <- singledf.anova$S[3]/sum(singledf.anova$S)
#list(pvalue=pvalue,Fstat=anova(singledf.out)$F[3],partial.r2=partial.r2)
list(pvalue=pvalue,singledf.out=singledf.out)
}
