stima.control <-
function(minbucket=NULL, crit = "f2", mincrit = 0.001, predtrunk = NULL, ref=1, sel = "none", ksel =2, predsel = NULL, cvvec=NULL,seed = 3){
	list(minbucket=minbucket,crit = crit, mincrit= mincrit, predtrunk = predtrunk, ref=ref, sel = sel,ksel=ksel, predsel = predsel, cvvec=cvvec,seed = seed)
}
