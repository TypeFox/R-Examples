writeGeno <-
function (x, fileName)
{
	vec = NULL
	vec = c("geneinfo", rownames(x))
	 write.table(t(vec), file = fileName, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
	write.table (t(x), fileName,sep = "\t", quote = FALSE, col.names = FALSE, append = TRUE)
}

