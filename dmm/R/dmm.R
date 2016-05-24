dmm <-
function(mdf,fixform = Ymat ~ 1,components=c("VarE(I)","VarG(Ia)"),cohortform=NULL,posdef=T,gls=F,glsopt=list(maxiter=200,bdamp=0.8,stoptol=0.01), dmeopt="qr",ncomp.pcr="rank",relmat="inline",dmekeep=F,dmekeepfit=F,traitspairwise=F,traitsblockwise=F,...) 
{
# dmm()
   UseMethod("dmm")
}
