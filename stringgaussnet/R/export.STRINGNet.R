export.STRINGNet <-
function (x, dirname, overwrite=F, ...)
{
	if (file.exists(dirname))
	{
		if(!overwrite){stop(paste(dirname,"already exists."))}
		unlink(dirname,recursive=T)
	}
	dir.create(dirname)
	write.table(x$Edges,paste(dirname,"Edges_attributes.txt",sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
	write.table(cbind(NodeName=rownames(x$DEGenes),x$DEGenes),paste(dirname,"DEGenes_attributes.txt",sep="/"),col.names=T,row.names=F,sep="\t",quote=F)
	if (!is.null(x$Annotations)){write.table(cbind(NodeName=rownames(x$Annotations),x$Annotations),paste(dirname,"Annotations_attributes.txt",sep="/"),col.names=T,row.names=F,sep="\t",quote=F)}
}
