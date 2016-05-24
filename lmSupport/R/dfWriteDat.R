dfWriteDat <-
function(Data, File, SubID='SubID')
#writes a data frame as a tab-delimited data file with standard Curtin lab format
#Will add rownames as a variable named SubID
#2011-02-03, JJC
{
  if(!is.element(SubID, names(Data)) && !is.null(SubID))
  {
    dNew = as.data.frame(rownames(Data))
    names(dNew)[1] = SubID
    dNew[2:(ncol(Data)+1)] = Data
    Data = dNew
  }

  write.table(Data, file = File, append = FALSE, quote = FALSE, sep = '\t',
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE)
}