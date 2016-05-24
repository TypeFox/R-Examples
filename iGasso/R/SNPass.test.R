SNPass.test <-
function(geno, pheno)
{
mydata = na.omit(cbind(geno, pheno))
nn = nrow(mydata)
P = rowsum(rep(1, nn), mydata[,1])/nn
DF = ncol(mydata) - 1
    names(DF) = "df"

    if (length(P) < 3)     # at least one genotype is absent
    {
    aaa = anova(glm(as.factor(geno) ~ pheno, family=binomial(link="logit")))
    S = as.vector(aaa["pheno", "Deviance"])
    p.value = 1-pchisq(S, DF)
    names(S) = "LRT statistic"
    method = "Logistic regression method (as one genotype is absent)"
    }
    else 
    {
    yy = rowsum(mydata[,-1], mydata[,1])
    w = as.vector(c(1-P[1], P[3]-P[1], -(1-P[3])) %*% yy)
    S = as.vector((w %*% solve(var(mydata[,-1]), w)))/nn/prod(1-P)
    p.value = 1-pchisq(S, DF)

    names(S) = "score statistic"
    method = "The proportional odds model method (Wang, 2012)"
    }
    
    
structure(list(statistic = S, p.value = p.value, method = method, parameter = DF,
              data.name = paste(deparse(substitute(geno)), "~", deparse(substitute(pheno)))),
              class = "htest") 
}
