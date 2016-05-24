variants_compile <-
function(omim,clinvar,uniprot,localPDB = paste(getwd(),"localPDB",sep="/")){
    clinvarDB <- read.delim(gzfile(paste(localPDB,"variant_summary.txt.gz",sep="/"))) 
    var.omim <- setdiff(omim[,"clinvarAccessions"],clinvar$RCVaccession) 
    omim.add <- clinvar.add <- uniprot.add <- c()   
    
# compile omim variants to clinvar 
    omim.rcvacces <- unlist(lapply(omim[,"clinvarAccessions"],function(x) unique(unlist(strsplit(x,";")))))
    clinvar.rcvacces <- unlist(lapply(clinvar[,"RCVaccession"],function(x) unique(unlist(strsplit(as.character(x),";")))))
    var.omim <- setdiff(omim.rcvacces,clinvar.rcvacces)
    var.omim.inclinvar <- var.omim[var.omim != "" & !is.na(var.omim)]
    for(i in var.omim.inclinvar){
        clinvarDB.i <- clinvarDB[grep(i,clinvarDB$RCVaccession),]
        clinvar.add <- rbind(clinvar.add,clinvarDB.i)
    }
    clinvar.add <- unique(clinvar.add)
    omim.add <- omim[omim[,"clinvarAccessions"] == "" | is.na(omim[,"clinvarAccessions"]),]
    if(is.matrix(omim.add)) {
         omim.add <- omim.add[omim.add[,"status"] == "live",]
         }else{
            omim.add <- omim.add[omim.add["status"] == "live"]
    }     
      
# compile uniprot variants to clinvar
    var.uniprot <- paste(uniprot[,1],unlist(lapply(uniprot[,4],str_trim)),sep=":")
    var.clinvar <- paste(clinvar[,"GeneSymbol"], unlist(lapply(as.character(clinvar[,"HGVS.p.."]), function(x) unlist(strsplit(x,":"))[2])),sep=":")
    var.uniprot.1 <- setdiff(var.uniprot,var.clinvar)
    uniprot.add <- uniprot[is.element(var.uniprot,var.uniprot.1),]

#add the additional variants into clinvar, to get the final variant set.  
    var2pheno <- clinvar
    var2pheno$Mutation.add <- ""     
    if(is.matrix(omim.add)) {
         nrow.omim.add <- nrow(omim.add)
         }else{
            if(length(omim.add) == 0){
               nrow.omim.add <- 0
               }else{
                 nrow.omim.add <- 1
            }     
    }     

    if(is.data.frame(uniprot.add)) {
         nrow.uniprot.add <- nrow(uniprot.add)
         }else{
            nrow.uniprot.add <- 1
    }     
    
    var.add <- matrix(,nrow= sum(nrow.omim.add, nrow(clinvar.add), nrow.uniprot.add),ncol=ncol(var2pheno))
    colnames(var.add) <- colnames(var2pheno)
    if(!is.null(clinvar.add)) 
        var.add[1:nrow(clinvar.add),1:29] <- as.matrix(clinvar.add)
    if(nrow.omim.add > 1 ) {
        var.add[sum(nrow(clinvar.add),1):sum(nrow(clinvar.add),nrow.omim.add),c("GeneSymbol","Chromosome","Cytogenetic","omim.phenotype","OtherIDs","Mutation.add","RS...dbSNP.")] <- 
                 as.matrix(omim.add[,c("Approved.Symbol","Chromosome","cytoLocation","Phenotype","variants.ID","mutations","dbSNPs")])
        }else if(nrow.omim.add == 1){
             var.add[sum(nrow(clinvar.add),1):sum(nrow(clinvar.add),nrow.omim.add),c("GeneSymbol","Chromosome","Cytogenetic","omim.phenotype","OtherIDs","Mutation.add","RS...dbSNP.")] <- 
                   as.matrix(omim.add[c("Approved.Symbol","Chromosome","cytoLocation","Phenotype","variants.ID","mutations","dbSNPs")])        
    }    
                 
    if(nrow.uniprot.add > 1) { 
        var.add[sum(nrow(clinvar.add),nrow.omim.add,1):nrow(var.add),c("GeneSymbol","omim.phenotype","HGVS.p..","RS...dbSNP.","ClinicalSignificance")] <- 
                 as.matrix(uniprot.add[,c("GeneSymbol","DiseaseName","AA.change","dbSNP","type")])
         }else if(nrow.uniprot.add == 1) { 
        var.add[sum(nrow(clinvar.add),nrow.omim.add,1):nrow(var.add),c("GeneSymbol","omim.phenotype","HGVS.p..","RS...dbSNP.","ClinicalSignificance")] <- 
                 as.matrix(uniprot.add[c("GeneSymbol","DiseaseName","AA.change","dbSNP","type")]) 
    }
    var2pheno <- rbind(var2pheno,var.add)
    return(var2pheno)
}
