cpg.everything.character <-
function(x,...) {
  Phenotype<-x
  if(!grepl("\\]",Phenotype)) {
    Phenotype<-strsplit(Phenotype,"[:()$:]")[[1]]
    Phenotype<-Phenotype[length(Phenotype)]           
               }
    Phenotype
          }
