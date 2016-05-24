#'70 pathways from MSigDB c2CP
#'
#'@name cancer_pathways
#'@rdname cancer_pathways
#'@aliases cancer_pathways vijver
#'
#'@references MJ van de Vijver,YD He, LJ van't Veer, H Dai, AAM Hart, DW Voskuil,
#'A gene-expression signature as a predictor of survival in breast cancer,
#'\emph{The New England Journal of Medicine}, 347(25):1999-2009, 2002.
#'
#'@references T Cai, G Tonini, X Lin, Kernel Machine Approach to Testing the
#'Significance of Multiple Genetic Markers for Risk Prediction, \emph{Biometrics},
#'67(3):975-986, 2011.
#'
#'@usage data("cancer_pathways")
#'
#'@format a list of 70 relevant pathways from an old version of MSigDB c2CP contening the Entrez IDs.
#'
#'@examples
#'data("cancer_pathways")
#'
#'##get the data from Vijver publication
#'
#'#clinical data
#'import_xls_from_zip <- function(urlPath, filename, zipname, skip=0){
#'  zipFile <- paste0(zipname, ".zip")
#'  download.file(paste0(urlPath, zipFile), zipFile)
#'  unzip(zipFile, exdir="./temp_unzip")
#'  xlsFile <- paste0("./temp_unzip/", filename, ".xls")
#'  res <- gdata::read.xls(xlsFile, skip=skip)
#'  unlink(zipFile)
#'  unlink("./temp_unzip", recursive=TRUE)
#'  return(res)
#'}
#'
#'BC_dat_clin <- import_xls_from_zip(urlPath="http://ccb.nki.nl/data/",
#'                                   filename="Table1_ClinicalData_Table",
#'                                   zipname="nejm_table1",
#'                                   skip=2
#'                                   )
#'BC_dat_clin <- BC_dat_clin[order(BC_dat_clin$SampleID), ]
#'col2rmv <- 1:ncol(BC_dat_clin)
#'BC_dat_clin$ID <- paste0("S", BC_dat_clin$SampleID)
#'rownames(BC_dat_clin) <- BC_dat_clin$ID
#'BC_dat_clin$evdeath <- BC_dat_clin$EVENTdeath
#'BC_dat_clin$tsurv <- BC_dat_clin$TIMEsurvival
#'BC_dat_clin$evmeta <- BC_dat_clin$EVENTmeta
#'BC_dat_clin$tmeta<- pmin(BC_dat_clin$TIMEsurvival, BC_dat_clin$TIMEmeta, na.rm=TRUE)
#'samples2rmv <- c("S28", "S122", "S123", "S124", "S133", "S138", "S139", "S141", "S221", "S222",
#'                 "S224", "S226", "S227", "S228", "S229", "S230", "S231", "S237", "S238", "S240",
#'                 "S241", "S248", "S250", "S251", "S252", "S254", "S292", "S317", "S342", "S371",
#'                 "S379", "S380", "S397", "S398", "S401")
#'BC_dat_clin <- BC_dat_clin[-which(BC_dat_clin$ID %in% samples2rmv), -col2rmv]
#'head(BC_dat_clin)
#'
#'
#'\dontrun{
#'
#'#import genomics data
#'urlPath="http://ccb.nki.nl/data/"
#'zipFile <- paste0("ZipFiles295Samples", ".zip")
#'download.file(paste0(urlPath, zipFile), zipFile)
#'unzip(zipFile, exdir="./temp_unzip")
#'unlink(zipFile)
#'unlink("./temp_unzip/Readme.txt", recursive=FALSE)
#'txtfiles <- list.files("./temp_unzip/")
#'BC_dat_exp <- NULL
#'for(f in txtfiles){
#'  temp_exp <- read.delim(paste0("./temp_unzip/", f))
#'  if(f==txtfiles[1]){
#'    gene_id <- as.character(temp_exp[-1, 1])
#'    gene_symbol <- as.character(temp_exp[-1, 2])
#'  }
#'  temp_exp <- temp_exp[-1, grep("Sample.", colnames(temp_exp))]
#'  colnames(temp_exp) <- gsub("Sample.", "S", colnames(temp_exp))
#'  if(f==txtfiles[1]){
#'    BC_dat_exp <- temp_exp
#'  }else{
#'    BC_dat_exp <- cbind(BC_dat_exp, temp_exp)
#'  }
#'}
#'BC_dat_exp_all <- cbind.data.frame("SYMBOL"=gene_symbol, BC_dat_exp[,  BC_dat_clin$ID])
#'unlink("./temp_unzip", recursive=TRUE)
#'
#'# translating the pathways from Entrez ID to gene symbol
#'library(org.Hs.eg.db)
#'x <- org.Hs.egSYMBOL
#'mapped_genes <- mappedkeys(x)
#'xx <- as.list(x[mapped_genes])
#'cancer_pathways_Symbol <- lapply(cancer_pathways, function(v){unlist(xx[v])})
#'sapply(cancer_pathways, function(x){length(intersect(x, rownames(BC_dat_exp)))/length(x)})
#'}
#'
#' @keywords datasets
#' @docType data
NULL