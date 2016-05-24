prepcombat.single.indep <-
function (ds1,ds2, i, j, batchID)
{
	batchID = batchID[c(grep(i,batchID),grep(j,batchID))]

	writeSamples(rbind(ds1, ds2),  batchID,"sampleFile")
	writeGeno(rbind(ds1,ds2), "genoFile")

	mat = compute.combat ("genoFile", "sampleFile")
	mat = t(mat)

	return(mat)
}

