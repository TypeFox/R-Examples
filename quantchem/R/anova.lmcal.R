"anova.lmcal" <-
function (object,...) 
{

obj=object;

res=list()

res$mandel=anova(obj$models$p1,obj$models$p2)
res$logmandel=anova(obj$models$l1,obj$models$l2)
res$table=anova(obj$models$p1,obj$models$p2,obj$models$p3,obj$models$p4)
res$wtable=anova(obj$models$P1,obj$models$P2,obj$models$P3,obj$models$P4)

resp = res;
for (i in 1:4) attr(resp[[i]],"heading") = NULL;

cat("\nMandel's test:\n")
print(resp$mandel);
cat("\nMandel's test on log-log models:\n")
print(resp$logmandel);
cat("\nANOVA table for unweighted models:\n")
print(resp$table);
cat("\nANOVA table for weighted models:\n")
print(resp$wtable);

invisible(res);

}

