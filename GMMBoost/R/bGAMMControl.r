bGAMMControl<-function (nue=0.1,add.fix=NULL,start=NULL,q_start=NULL, OPT=TRUE,nbasis=20,spline.degree=3,diff.ord=2,sel.method="aic",steps=500,method="EM",overdispersion=FALSE)
{
list(nue = nue, add.fix = add.fix, start = start, q_start = q_start, OPT = OPT,
     nbasis = nbasis, spline.degree = spline.degree, diff.ord = diff.ord,
     sel.method = sel.method, steps = steps, method = method, overdispersion = overdispersion)
}
