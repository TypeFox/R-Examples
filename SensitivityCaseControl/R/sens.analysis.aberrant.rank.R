sens.analysis.aberrant.rank <-
function(cases.exposed,referents.exposed,case.rank,no.referents,Gamma){
mi=cases.exposed+referents.exposed;
Ji=no.referents+1;
pi.bar=mi/(mi+Gamma*(Ji-mi));
pi.double.bar=Gamma*mi/(Gamma*mi+(Ji-mi));
di=case.rank; # Score for Mantel-Haenszel test
teststat=sum(di*cases.exposed);
lower.bound.pval=1-pnorm((teststat-sum(di*pi.bar))/sqrt(sum(di^2*pi.bar*(1-pi.bar))));
upper.bound.pval=1-pnorm((teststat-sum(di*pi.double.bar))/(sqrt(sum(di^2*pi.double.bar*(1-pi.double.bar)))));
list(lower.bound.pval=lower.bound.pval,upper.bound.pval=upper.bound.pval);
}

