OrdinalBoostControl<-function (nue=0.1, lin=NULL, katvar=NULL, start=NULL, q_start=NULL, 
                               OPT=TRUE, sel.method="aic", steps=100, method="EM", 
                               maxIter=500,print.iter.final=FALSE,eps.final=1e-5)
{

list(nue = nue, lin = lin, katvar = katvar, start = start, q_start = q_start, OPT = OPT,
        sel.method = sel.method, steps = steps, method = method, maxIter = maxIter,
        print.iter.final=print.iter.final, eps.final=eps.final)
}
