"anova.nlscal" <-
function (object,...) 
{

obj=object

res=list()

res$a=try(anova(obj$models$a1,obj$models$a2))
res$g=try(anova(obj$models$g1,obj$models$g2))

resp = res;
for (i in 1:2) attr(resp[[i]],"heading") = NULL;

if (inherits(res$a,"anova")) { cat("\nANOVA for asymptotic models:\n");
print(resp$a); }
if (inherits(res$g,"anova")) { cat("\nANOVA for logistic models:\n");
print(resp$g); }

invisible(res);

}

