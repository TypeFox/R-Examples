writeSamples <-
function(x, batchID, fileName)
{
	sampleInfo = cbind(rownames(x),batchID)
        colnames(sampleInfo) = c("Array.name", "Batch")
        write.table (sampleInfo, fileName, sep = "\t", quote = FALSE, row.names = FALSE)

}

