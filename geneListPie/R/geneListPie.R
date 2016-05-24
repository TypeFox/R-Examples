##read a single gene list from file in a column or a row format
load.genelist<-function(file="genelist.txt", format="column", sep="\t")
{
	if (format=="row")
	{	glist<-readLines(file) #gene list
		glist<-strsplit(glist[1], sep)[[1]]
	}
	else if(format=="column")
	{
		glist<-read.delim(file, as.is = TRUE)[[1]]
	}
	else{
		print ("Gene list has to be a single row or a single column")
	}
	glist<-unique(toupper(glist))
	return (glist)
}

##read gene sets in MSigDB format
load.genesets<-function(file="geneset.gmt.txt")
{
	gs<-readLines(file) #gene sets
	return (gs)
}

##compare the gene list against all the gene sets
##return the genes in each category
geneListProfile<-function(gs, glist, threshold=10)
{
	num.gs<-length(gs)
	labels<-NULL
	sizes<-NULL
	symbols<-NULL
	othercommon<-NULL
	for (i in 1: num.gs){
		line<-strsplit(gs[i],"\t")[[1]]
		sym<-toupper(line[3:length(line)])
		common<-intersect(glist, sym)
		size<-length(common)
		if(size<threshold){
			othercommon=union(othercommon, common)
		}
		else{
			labels<-c(labels, line[1])
			sizes<-c(sizes,length(common))
			symbols<-c(symbols, list(common))
		}
	}
	if(length(othercommon)>0)
	{
		labels<-c(labels, "Others")
		sizes<-c(sizes, length(othercommon))
		symbols<-c(symbols, list(othercommon))
	}
	return (list(labels=labels, sizes=sizes, symbols=symbols))
}

##Print the overlapped genes in each category in MSigDB 
##or two column (category, size) format for plotting in e.g. Excel
printGeneListProfile<-function(r, file="", format=NULL)
{
	if(is.character(format))
	{
		if(toupper(format)=="MSIGDB")
		{
			for (i in 1:length(r$labels))
			{
				cat(paste(r$labels[i],"\t"), file=file, append=TRUE)
				cat(paste("na","\t"), file=file, append=TRUE)
				for (j in 1: length(r$symbols[[i]]))
				{
					cat(paste(r$symbols[[i]][j],"\t"), file=file, append=TRUE)
				}
				cat("\n", file=file, append=TRUE)
			}
		}
	}
	else
		{
		for (i in 1:length(r$labels))
		{
			cat(paste(r$labels[i],"\t"), file=file, append=TRUE)
			cat(paste(r$sizes[i],"\n"), file=file, append=TRUE)
		}
	}
}


