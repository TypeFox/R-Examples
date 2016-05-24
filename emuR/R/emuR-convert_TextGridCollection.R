requireNamespace("RSQLite", quietly = T)

##' Convert a TextGridCollection (e.g. .wav & .TextGrid files) to emuDB
##' 
##' Converts a TextGridCollection to an emuDB by searching a given directory for .wav & .TextGrid files (default
##' extensions) with the same base name. First, the function generates a file pair list 
##' containing paths to files with the same base
##' name. It then generates an emuDB DBconfig based on the first TextGrid in this list which specifies 
##' the allowed level names and types in the new emuDB. After this it converts all file pairs to the new format,
##' checking whether they comply to the newly generated database configuration. For 
##' more information on the emuDB format see \code{vignette{emuDB}}.
##' Note that since Praat TextGrids do not permit explicit hierarchical structures, the created emuDB does not contain
##' any links or link definitions. You can however use the \code{\link{autobuild_linkFromTimes}} function after the conversion process
##' to automatically build a hierarchal structure.
##' 
##' @param dir path to directory containing the TextGridCollection (nested directory structures are permitted as the 
##' function recursively searches through directories, generating the session names from dir. structure)
##' @param dbName name given to the new emuDB
##' @param targetDir directory where to save the new emuDB
##' @param tgExt extension of TextGrid files (default=TextGrid, meaning file names of the form baseName.TextGrid)
##' @param audioExt extension of audio files (default=wav, meaning file names of the form baseName.wav)
##' @param tierNames character vector containing names of tiers to extract and convert. If NULL (the default) all
##' tiers are converted.
##' @param verbose display infos & show progress bar
##' @export
##' @import tools
##' @return NULL
##' @examples 
##' \dontrun{
##' 
##' ##########################################################
##' # prerequisite: directory containing .wav & .TextGrid files
##' # (see \code{?create_emuRdemoData} how to create demo data)
##' 
##' # convert TextGridCollection and store 
##' # new emuDB in folder provided by tempdir()
##' convert_TextGridCollection(dir = "/path/to/directory/",
##'                            dbName = "myTGcolDB",
##'                            targetDir = tempdir())
##' 
##' 
##' # same as above but this time only convert 
##' # the information stored in the "Syllable" and "Phonetic" tiers
##' convert_TextGridCollection(dir = "/path/to/directory/", 
##'                            dbName = "myTGcolDB", 
##'                            targetDir = tempdir(),
##'                            tierNames = c("Syllable", "Phonetic"))
##'
##'} 
convert_TextGridCollection <- function(dir, dbName, 
                                       targetDir, tgExt = 'TextGrid', 
                                       audioExt = 'wav', tierNames = NULL, 
                                       verbose = TRUE){
  # normalize paths
  dir = suppressWarnings(normalizePath(dir))
  targetDir = suppressWarnings(normalizePath(targetDir))
  
  # check if dir exists
  if(!file.exists(dir)){
    stop("dir does not exist!")
  }
  # create
  if(!file.exists(targetDir)){
    res=dir.create(targetDir,recursive = TRUE)
    if(!res){
      stop("Could not create target directory: ",targetDir," !\n")
    }
  }
  
  basePath=file.path(targetDir, paste0(dbName, emuDB.suffix))
  # check if base path dir already exists
  if(file.exists(basePath)){
    stop('The directory ', basePath, ' already exists. Can not generate new emuDB if directory called ', dbName, ' already exists!')
  }else{
    res=dir.create(basePath)
    if(!res){
      stop("Could not create base directory: ",basePath," !\n")
    }
  }
  
  # generate file pair list
  fpl = create_filePairList(dir, dir, audioExt, tgExt)
  
  progress = 0
  
  if(verbose){
    cat("INFO: Loading TextGridCollection containing", length(fpl[,1]), "file pairs...\n")
    pb = utils::txtProgressBar(min = 0, max = length(fpl[,1]), initial = progress, style=3)
    utils::setTxtProgressBar(pb, progress)
  }
  
  # gereate DBconfig from first TextGrid in fpl
  DBconfig = create_DBconfigFromTextGrid(fpl[1,2], dbName, basePath,tierNames)
  
  # create tmp dbHandle
  dbHandle = emuDBhandle(dbName, basePath, DBconfig$UUID, connectionPath = ":memory:")
  
  # store to tmp DBI
  add_emuDbDBI(dbHandle)
  
  # store db DBconfig file
  store_DBconfig(dbHandle, DBconfig)
  
  # allBundles object to hold bundles without levels and links
  allBundles = list()
  
  
  # loop through fpl
  for(i in 1:dim(fpl)[1]){
    
    # create session name
    sesName = gsub('^_', '', gsub(.Platform$file.sep, '_', gsub(normalizePath(dir, winslash = .Platform$file.sep),'',dirname(normalizePath(fpl[i,1], winslash = .Platform$file.sep)))))
    
    # session file path
    if(sesName == ""){
      sesName = "0000"
      sfp = file.path(basePath, paste0(sesName, session.suffix))
    }else{
      sfp = file.path(basePath, paste0(sesName, session.suffix))
    }
    # create session dir
    if(!dir.exists(sfp)){
      res = dir.create(sfp)
      if(!res){
        stop("Could not create session directory: ", sfp, " !\n")
      }
    }

    # create session entry if it doesn't already exist
    sesDF = DBI::dbGetQuery(dbHandle$connection, paste0("SELECT * FROM session WHERE name = '", sesName, "'"))
    if(nrow(sesDF) == 0){    
      DBI::dbGetQuery(dbHandle$connection, paste0("INSERT INTO session VALUES('", dbHandle$UUID, "', '", sesName, "')"))
    }

    # media file
    mfPath = fpl[i,1]
    mfBn = basename(mfPath)
    
    # get sampleRate of audio file
    asspObj = wrassp::read.AsspDataObj(mfPath, begin = 0, end = 1, samples = T)
    sampleRate = attributes(asspObj)$sampleRate
    # create bundle name
    bndlName = tools::file_path_sans_ext(basename(fpl[i,1]))
    
    # parse TextGrid
    bundleAnnotDFs = TextGridToBundleAnnotDFs(fpl[i,2], sampleRate = sampleRate, name = bndlName, annotates = paste0(bndlName, ".wav"))
    
    # remove unwanted levels
    if(!is.null(tierNames)){
      # filter items
      bundleAnnotDFs$items = dplyr::filter_(bundleAnnotDFs$items, ~(level %in% tierNames))
      # filter labels
      bundleAnnotDFs$labels = dplyr::filter_(bundleAnnotDFs$labels, ~(name %in% tierNames))
    }
    
    # add to bundle table
    add_bundleDBI(dbHandle, sessionName = sesName, name = bndlName, annotates = bundleAnnotDFs$annotates, sampleRate = bundleAnnotDFs$sampleRate, MD5annotJSON = "")
    # add to items, links, labels tables
    store_bundleAnnotDFsDBI(dbHandle, bundleAnnotDFs, sessionName = sesName, bundleName = bndlName)
    
    
    # validate bundle
    valRes = validate_bundleDBI(dbHandle, session=sesName, bundle=bndlName)
    
    if(valRes$type != 'SUCCESS'){
      stop('Parsed TextGrid did not pass validator! The validator message is: ', valRes$message)
    }
    
    # create bndl folder
    bDir = paste0(bndlName, bundle.dir.suffix)
    bfp = file.path(sfp,bDir)
    res = dir.create(bfp)
    if(!res){
      stop("Could not create bundle directory ",bfp," !\n")
    }
    
    # store media file
    newMfPath = file.path(bfp, mfBn)
    if(file.exists(mfPath)){
      file.copy(from = mfPath, to = newMfPath)
    }else{
      stop("Media file :'", mfPath, "' does not exist!")
    }
    
    
    # update pb
    if(verbose){
      utils::setTxtProgressBar(pb, i)
    }
    
  }

  if(verbose){
    cat('\n') # hack to have newline after pb
  }

  # store all annotations
  rewrite_allAnnots(dbHandle, verbose = verbose)
  
  
  
}

# FOR DEVELOPMENT
# library('testthat')
# test_file('tests/testthat/test_aaa_initData.R')
# test_file('tests/testthat/test_emuR-convert_TextGridCollection.R')
