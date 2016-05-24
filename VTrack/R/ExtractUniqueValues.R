ExtractUniqueValues <-
function(sInputFile,iFieldToExtract)
{
ivals <- unique(sInputFile[,iFieldToExtract])
ivals[order(ivals)]
}
