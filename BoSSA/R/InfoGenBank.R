`InfoGenBank` <-
function(X,tsleep=3)
{	
	gbout <- paste("AccNb","Organism","Isolate","Taxonomy","DateSub","DateEch","Host","HostTaxoID","HostFamily","HostGenus","HostSpecies","MispellHost","Location","GPS","Authors","Title","Journal","PubmedURL",sep="\t")
	for(i in 1:length(X))
	{
		gbi <- InfoGB(X[i],tsleep)
		gbout[i+1] <- gbi
	}
gbout
}

