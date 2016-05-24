X.design <-
function(formula, data, Xmethod = "treatments")
{

if(Xmethod == "treat.first.0") {
options(contrasts=c("contr.treatment", "contr.treatment"))
X <- model.matrix(formula, data = data)
options(contrasts=c("contr.treatment", "contr.poly"))
}

if(Xmethod == "sum.to.0") {
options(contrasts=c("contr.sum", "contr.sum"))
X <- model.matrix(formula, data = data)
options(contrasts=c("contr.treatment", "contr.poly"))
}
if(Xmethod == "over.par") {
terms.list <- attr(terms(formula),"term.labels")
X <- NULL
if(attr(terms(formula),"intercept") == 1)
X <- rep(1, times = length(data[,1]))
for (i in 1:length(terms.list)) {
form1 <- formula(paste("~ ", terms.list[i], " - 1"))
X <- cbind(X,model.matrix(form1, data = data))
}
}
X
}

