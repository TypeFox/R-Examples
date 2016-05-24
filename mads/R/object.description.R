object.description <- function(){
  cat("\nDESCRIPTION OF MA OBJECT")
  cat("\n------------------------\n")
  cat("\nAn object of class ma comtains three components. The first component \nis an object of class ma.analysis named analysis.options, this gives \na summary of the analysis options selected. The second component \nis an object of class ma.allspecies named species, this provides \nsummaries relating to both the detection function models as well as \nestimated of density, abundance and where relevant expected cluster \nsize. The third component is an object of class ma.allunid named \nunidentified, this gives summaries of the models for the detection \nfunctions selected for the unidentified categories.\n")
}

