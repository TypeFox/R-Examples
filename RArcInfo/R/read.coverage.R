read.coverage<-function(datadir, coverage)
{
	infodir<-paste(c(datadir,"/info/"), collapse="")
	covdir<-paste(c(datadir, "/",coverage), collapse="")

	cov.arc<-get.arcdata(datadir, coverage)
	cov.cnt<-get.cntdata(datadir, coverage)
	cov.bnd<-get.bnddata(infodir, paste( c(casefold(coverage, upper=TRUE), ".BND"), collapse="") )
	cov.lab<-get.labdata(datadir, coverage)

	palfiles<-dir(covdir, pattern="pal.adf")
	palfiles<-c(palfiles, dir(covdir, pattern="*.pal") )

	if( length(palfiles)==1 )
		cov.pal<-get.paldata(datadir, coverage, palfiles[1])
	else
		cov.pal<-NULL

	
	if( length(dir(covdir, pattern="tol.adf"))>0 )
		cov.tol<-get.toldata(datadir, coverage)
	else
	{
		cov.tol<-get.toldata(datadir, coverage, "par.adf")
	}
	
	tblnames<-get.tablenames(infodir)
	
	pattern<-paste( c(casefold(coverage, upper=TRUE), ".*"), collapse="")

	cov.tables<-tblnames[grep(pattern, tblnames[[1]]),]

	cov<-list(datadir=datadir, coverage=coverage, arc=cov.arc, bnd=cov.bnd, cnt=cov.cnt,  lab=cov.lab, pal=cov.pal, palfiles=palfiles, tblnames=cov.tables, tol=cov.tol)

	return(cov)
}
