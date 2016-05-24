localPDB <-
function(localPDB = paste(getwd(),"localPDB",sep="/"),PDB="all", omim.url = NULL){
    download.path = localPDB
    if(!file.exists(download.path))
          dir.create(download.path )
    options(timeout = 1000)
    if(PDB == "all"){
      HPO <- "yes"; Orphanet <- "yes"; MedGen <- "yes"; HGNC <- "yes"; GeneReview <- "yes"; ClinVar <- "yes"; Uniprot <- "yes"
    }
    
       refFlat <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz"
       if( !file.exists(paste(download.path,"refFlat.txt.gz",sep="/")))
           download.file(refFlat,paste(download.path,"refFlat.txt.gz",sep="/"),method="auto")

   if(is.null(omim.url)){ 
      print("Warning: please make sure you have localized the OMIM file morbidmap! if NOT, you should apply for an OMIM account and get the URL from http://omim.org/downloads.")
      }else if(!is.null(omim.url)){
        morbidmap <- paste(omim.url,"morbidmap", sep = "/")
        if( !file.exists(paste(download.path,"morbidmap",sep="/")))
           download.file(morbidmap,paste(download.path,"morbidmap",sep="/"),method="auto")   
   }
   
   if(HPO == "yes" | toupper(PDB) == "HPO"){
       HPO <- "http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab"
       diseases_to_genes <- "http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/diseases_to_genes.txt"
       if( !file.exists(paste(download.path,"phenotype_annotation.tab",sep="/")))
           download.file(HPO,paste(download.path,"phenotype_annotation.tab",sep="/"),method="auto")
       if( !file.exists(paste(download.path,"diseases_to_genes.txt",sep="/")))
           download.file(diseases_to_genes,paste(download.path,"diseases_to_genes.txt",sep="/"),method="curl")     
   }
      
   if(Orphanet == "yes" | toupper(PDB) == "ORPHANET"){
       orphanet <- "http://www.orphadata.org/data/xml/en_product6.xml"
       if(!file.exists(paste(download.path,"en_product6.xml",sep="/")))
            download.file(orphanet,paste(download.path,"en_product6.xml",sep="/"),method="auto")          
   }

   if(MedGen == "yes" | toupper(PDB) == "MEDGEN"){
       medgene.names <- "ftp://ftp.ncbi.nlm.nih.gov/pub/medgen/csv/NAMES.csv.gz"
       if( !file.exists(paste(download.path,"NAMES.csv.gz",sep="/")))
           download.file(medgene.names,paste(download.path,"NAMES.csv.gz",sep="/"),method="auto")
   }
      
   if(HGNC == "yes" | toupper(PDB) == "HGNC"){
       hgnc <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc_complete_set.txt.gz"
       if( !file.exists(paste(download.path,"hgnc_complete_set.txt.gz",sep="/")))
           download.file(hgnc,paste(download.path,"hgnc_complete_set.txt.gz",sep="/"),method="auto")
   }
   
   if(GeneReview == "yes" | toupper(PDB) == "GENEREVIEW"){
       genereview <- "ftp://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/GRtitle_shortname_NBKid.txt"
       if( !file.exists(paste(download.path,"GRtitle_shortname_NBKid.txt",sep="/")))
           download.file(genereview,paste(download.path,"GRtitle_shortname_NBKid.txt",sep="/"),method="auto")
   }
   
   if(ClinVar == "yes" | toupper(PDB) == "CLINVAR"){
       clinvar <- "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
       if( !file.exists(paste(download.path,"variant_summary.txt.gz",sep="/")))
           download.file(clinvar,paste(download.path,"variant_summary.txt.gz",sep="/"),method="auto")
       gene2dis <- "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id"
       if(!file.exists(paste(download.path,"gene_condition_source_id",sep="/")) )
           download.file(gene2dis,paste(download.path,"gene_condition_source_id",sep="/"),method="auto")
   }
   
   if(Uniprot == "yes" | toupper(PDB) == "UNIPROT"){
       uniprot <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
       if( !file.exists(paste(download.path,"humsavar.txt",sep="/")))
           download.file(uniprot,paste(download.path,"humsavar.txt",sep="/"),method="auto")
   }
   
   print(paste("Congratulations! Public databases have been localized the in ",localPDB,".",sep=""))
}
