working.cltsk.par <-
function(query,obs,th,bins,vth,vlen,llim,verbose,Large,future)
{
 	## working function for cltsk parallel
	apply(query,1,working.cltsk,obs=obs,th=th,bins=bins,vth=vth,vlen=vlen,llim=llim,verbose=verbose,Large=Large,future=future)
}
