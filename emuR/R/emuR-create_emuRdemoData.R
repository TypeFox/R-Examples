##' Create demo data for the emuR package
##' 
##' Create a folder within the folder specified
##' by the dir argument called emuR_demoData.
##' This folder contains the folders:
##' \itemize{
##' \item{ae_emuDB: }{Containing an emuDB that adheres to the new format specification 
##' (as expected by the \code{\link{load_emuDB}} function). See \code{vignette(emuDB)} 
##' for more information on this database format.}
##' \item{BPF_collection: }{Containing a BAS Patitur Format (BPF) file collection (as 
##' expected by the \code{\link{convert_BPFCollection}} function)}
##' \item{legacy_ae: }{Containing a legacyEmuDB (as expected by the 
##' \code{\link{convert_legacyEmuDB}} function)}
##' \item{TextGrid_collection: }{Containing a TextGrid file collection 
##' (as expected from the \code{\link{convert_TextGridCollection}} function)}
##' }
##' @param dir directory to create demo data in (default= tempdir())
##' @param precache creates an on-file-system cache for the ae emuDB to allow fast loading
##' (see \code{load_emuDB} for details about the emuDB file cache)
##' @export
##' @examples 
##' \dontrun{
##' 
##' # create demo data directory in directory
##' # provided by the tempdir function
##' create_emuRdemoData(dir = tempdir())
##' }
create_emuRdemoData <- function(dir = tempdir(), precache = FALSE){
  
  ddPath = file.path(dir,"emuR_demoData")
  
  path2data = system.file("extdata", package = "emuR")
  
  if(file.exists(ddPath)){
    stop("Path '", ddPath,"' already exists!")
  }
  
  dir.create(ddPath)
  #################################
  # create ae
  configPath = list.files(path2data, pattern = "DBconfig.json$", recursive = T, full.names = T)
  wavPaths = list.files(path2data, pattern = ".wav$", recursive = T, full.names = T)
  annotPaths = list.files(path2data, pattern = "_annot.json$", recursive = T, full.names = T)
  aePath = file.path(ddPath, paste0("ae", emuDB.suffix))
  
  dir.create(aePath)
  
  file.copy(configPath, aePath)
  
  sesPath = file.path(aePath, "0000_ses")
  dir.create(sesPath)
  
  for(p in wavPaths){
    bndlName = gsub(".wav$", "", basename(p))
    bndlPath = file.path(sesPath, paste0(bndlName, "_bndl"))
    dir.create(bndlPath)
    
    file.copy(p, bndlPath)
    idx = grep(paste0(bndlName, "_annot.json$"), annotPaths)
    file.copy(annotPaths[idx], bndlPath)
    
  }
  
  # calc dft and fms files
  wps = list.files(sesPath, pattern = ".wav$", recursive = T, full.names = T)
  wrassp::dftSpectrum(wps)
  wrassp::forest(wps)
  
  
  # generate cache of ae emuDB
  if(precache){
    dbHandle = load_emuDB(aePath, inMemoryCache = F, verbose = F)
  }
  
  ####################################
  # create TextGrid_collection and BPF_collection
  fpltgc = create_filePairList(path2data, path2data, "wav", "TextGrid")
  fplbpf_original = create_filePairList(path2data, path2data, "wav", "par")
  fplbpf_manipulated = create_filePairList(path2data, path2data, "wav", "parmanipulated")
  tgcPath = file.path(ddPath, "TextGrid_collection")
  bpfPath_original = file.path(ddPath, "BPF_collection")
  
  dir.create(tgcPath)
  dir.create(bpfPath_original)
  
  file.copy(fpltgc[,1], tgcPath)
  file.copy(fpltgc[,2], tgcPath)
  file.copy(fplbpf_original[,1], bpfPath_original)
  file.copy(fplbpf_original[,2], bpfPath_original)
  
  #################################
  # create legacyEmuDB
  tplPath = list.files(path2data, pattern = ".tpl$", recursive = T, full.names = T)
  wavPaths = list.files(path2data, pattern = ".wav$", recursive = T, full.names = T)
  hlbPaths = list.files(path2data, pattern = "hlb$", recursive = T, full.names = T)
  labPaths = list.files(path2data, pattern = "lab$", recursive = T, full.names = T)
  tonePaths = list.files(path2data, pattern = "tone$", recursive = T, full.names = T)
  
  legacyAePath = file.path(ddPath, "legacy_ae")
  dir.create(legacyAePath)
  labelsPath = file.path(legacyAePath, "labels")
  dir.create(labelsPath)
  signalsPath = file.path(legacyAePath, "signals")
  dir.create(signalsPath)
  
  # copy files
  file.copy(tplPath, legacyAePath)
  file.copy(wavPaths, signalsPath)
  file.copy(hlbPaths, labelsPath)
  file.copy(labPaths, labelsPath)
  file.copy(tonePaths, labelsPath)
  
  # calc dft and fms files
  wps = list.files(signalsPath, pattern = ".wav$", recursive = T, full.names = T)
  wrassp::dftSpectrum(wps)
  wrassp::forest(wps)
  
}

## create manipulated BPF_collection
##
## @param dir directory in that the BPF_collection is created
create_BPFcollectionManipulated = function(dir){
  
  path2data = system.file("extdata", package = "emuR")
  
  bpfPath_manipulated = file.path(dir, "BPF_collection_manipulated")
  
  if(file.exists(bpfPath_manipulated)){
    stop("Path '", bpfPath_manipulated,"' already exists!")
  }
  
  dir.create(bpfPath_manipulated)
  
  fplbpf_manipulated = create_filePairList(path2data, path2data, "wav", "parmanipulated")
  
  dir.create(file.path(bpfPath_manipulated, "0000"))
  file.copy(fplbpf_manipulated[,1], file.path(bpfPath_manipulated, "0000"))
  file.copy(fplbpf_manipulated[,2], file.path(bpfPath_manipulated, "0000"))
  
}

########################
# FOR DEVELOPMENT 
# unlink(file.path(tempdir(),"emuR_demoData"), recursive = T)
# create_emuRdemoData(precache = T)
