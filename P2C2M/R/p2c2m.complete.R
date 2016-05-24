p2c2m.complete <-
function (path="/home/user/Desktop/", xml.file="beast.xml", 
                           descr.stats="LCWT,NDC", beast.vers="1.8", 
                           single.allele=c("O"), num.reps=1, 
                           use.sorted=FALSE, use.mpi=FALSE, 
                           save.metadata=FALSE, verbose=FALSE,
                           dbg=FALSE) {
  # Descr:  initiates the entire package by starting the wrapper
  # Deps:   p2c2m.readstarb, p2c2m.analyze
  # I/p:    path = absolute path to Working directory
  #         xml.file = name of Beauti XML infile
  #         descr.stats = list of descriptive statistics
  #         beast.vers = version fo BEAST
  #         single.allele = name of outgroup
  #         num.reps = number of replicates
  #         use.sorted = decision if to use sorting
  #         use.mpi = decision if to use MPI
  #         save.metadata = decision if metadata to be saved
  #         verbose = flag if status output printed on screen
  #         dbg = debugging flag

#########################
# 1. Set global options #
#########################
  # TFL avoids line-wrapping in command "print()"
  options(width=1000)
  # TFL avoids numbers in scientific e-notation (i.e. "1.2345e+03" = 1234.5)
  options(scipen=1000)

#######################
# 2. Parse user input #
#######################
  # Check if selected descriptive stats are valid
  descrStats = toupper(sort(unlist(strsplit(descr.stats, split=","))))
  slctn = which(descrStats %in% c("LCWT","COAL","GSI","NDC"))
  # Error handling
  if (length(descrStats) != length(slctn)) {
      stop(cat("\nERROR: Incorrect spec. of descript. stat.\n"))
  }

################################
# 3. Set environment variables #
################################
  # Note: I am not using true global variables (i.e. envir = .GlobalEnv), 
  #       because such could have interfered with user-space.
  assign("P2C2M_flg_xmlFile", xml.file, envir = P2C2M_globalVars)
  # Note: not descr.stats, but descrStats
  assign("P2C2M_flg_dscrStat", descrStats, envir = P2C2M_globalVars)
  assign("P2C2M_flg_beastV", beast.vers, envir = P2C2M_globalVars) 
  assign("P2C2M_flg_snglAl", single.allele, envir = P2C2M_globalVars)
  assign("P2C2M_flg_nReps", num.reps, envir = P2C2M_globalVars)
  assign("P2C2M_flg_srtBool", use.sorted, envir = P2C2M_globalVars)
  assign("P2C2M_flg_mpiBool", use.mpi, envir = P2C2M_globalVars)
  assign("P2C2M_flg_vrbBool", verbose, envir = P2C2M_globalVars)
  assign("P2C2M_flg_dbgBool", dbg, envir = P2C2M_globalVars)

##########################
# 4. Conducting analyses #
##########################
  # Generate indata
  inData = p2c2m.readstarb(path, xml.file)
  
  # Generate results
  metaData = p2c2m.analyze(xml.file, inData)
        
  # Calculate descriptive stats
  resultData = stats.main(path, xml.file, inData$loci$dta, metaData)

##########################
# 5. Summarizing outdata #
##########################
  outData = list()
  outData$results = resultData

  if (save.metadata) {
    outData$inData = inData
    outData$metaData = metaData
  }
  
return(outData)
}
