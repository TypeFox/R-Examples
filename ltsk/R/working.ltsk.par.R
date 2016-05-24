working.ltsk.par <-
function(query,obs,th,vth,vlen,llim,verbose,Large,future)
{
 	## working function for ltsk parallel
	apply(query,1,working.ltsk,obs=obs,th=th,vth=vth,vlen=vlen,llim=llim,verbose=verbose,Large=Large,future=future)
}
