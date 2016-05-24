autoevenn <-function(FolderPath="", pathRes="", annot=TRUE, ud=TRUE, prop=FALSE,
		overlaps=FALSE, display=FALSE, couleurs="", VennBar=FALSE, transp=0.5,
		Solid=TRUE, Profils=FALSE, ColorTxt="", colBlack=FALSE, Ptest=FALSE,
		tUD=NULL, tUDp=NULL, tnoUD=NULL, Gtype="png",lw=1, NutShell=TRUE, VennClust=TRUE)
{
	if(FolderPath=="")
	{
		write(paste("FolderPath is empty, it must contain the directory where are placed the Folders > Lists" , sep=""), file="")
		flush.console()
		FolderPath = choose.dir() 
		break
	}
	
	Venns = list.dirs(FolderPath, full.names=TRUE)
	
	if(pathRes!="")
	{
		FolderDest = pathRes
	}else{
		FolderDest = paste(getwd(), "/Venn.diagrams/", sep="")
	}	
	dir.create(file.path(FolderDest), showWarnings = FALSE)
	
	Venns = Venns[Venns!=FolderPath]  # filtre le repertoire principal de la liste
	write(paste(length(Venns), " Venn", if(length(Venns)>1){paste("s", sep="")}, sep=""), file="")
	for(V in 1:length(Venns))
	{
		write(paste(V, " / ", length(Venns), ": ", basename(Venns[V]), sep=""), file="")
		evenn(annot=annot, pathRes=FolderDest, pathLists=Venns[V], ud=ud, prop=prop, overlaps=overlaps,
				display=display, couleurs=couleurs,	VennBar=VennBar, transp=transp, Solid=Solid,
				Profils=Profils, ColorTxt=ColorTxt, colBlack=colBlack, Ptest=Ptest, tUD=tUD, 
				tUDp=tUDp, tnoUD=tnoUD, Gtype=Gtype, lw=lw, NutShell=NutShell, VennClust=VennClust)
	}
}
