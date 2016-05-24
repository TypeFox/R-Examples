p2c2m.readstarb <-
function (apWd, inFn, prmFHndl) {
  # Descr:  coordinates the reading of all input
  # Deps:   readtree.stree
  #         readtree.phybase
  # I/p:    apWd = absolute path to working directory
  #         inFn = infile name (name of Beauti XML outd file)
  #         inLogFn
  #         prmFHndl

  beastVers = get("P2C2M_flg_beastV", envir=P2C2M_globalVars)
  descrStats = get("P2C2M_flg_dscrStat", envir=P2C2M_globalVars)
  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  verboseBool = get("P2C2M_flg_vrbBool", envir=P2C2M_globalVars)

  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> readstarb.read", fg="red"),
        sep="")
  }

#####################################
# 1. Parsing data via Python script #
#####################################
  loghelpers.prntmngr("Parsing data from xml inputfile", uprFlg=T)
  # Specify name of parsing script
  pySc = system.file("exec", paste("parseBEAUTiXML_", beastVers, ".py",
                                   sep=""), package="P2C2M")
  # Load infile in Python
  rPython::python.assign("inFn", paste(apWd, inFn, sep=""))
  rPython::python.load(pySc)

######################################
# 2. Saving parsed info to variables #
######################################
  # Loading data of Species-Allele matrix into variable "SAmtrx"
  SAmtrx = do.call(rbind, rPython::python.get("SAmtrx"))
  colnames(SAmtrx) = c("speciesName", "alleleName")
  # Loading data of Locus-Ploidy-Matrix into variable "LPmtrx"
  LPmtrx = do.call(rbind, rPython::python.get("LPmtrx"))
  colnames(LPmtrx) = c("locus", "ploidy")
  loci = LPmtrx[,1]
  ploidy = as.numeric(LPmtrx[,2])
  # CURRENTLY PUT ON HOLD:
  #    # Loading data of Gene-Partition-matrix into variable "GPmtrx"        
  #    GPmtrx = do.call(rbind, rPython::python.get("GPmtrx"))
  #    colnames(GPmtrx) = c("genes","partitions")
  # CURRENTLY PUT ON HOLD:
  #    # Loading alignment-dictionary into variable "ALNdict"
  #    alignm = list()
  #    ALNdict = rPython::python.get("ALNdict")
  #    for (i in 1:length(loci))
  #        {
  #        # "order(names(ALNdict[[i]]))" assures that genes and loci are 
  #        # arranged in the same order as in LPmtrx
  #        listorder = order(names(ALNdict[[i]]))
  #        alignm[[loci[i]]] = ALNdict[[i]][listorder]
  #        }

##################################################
# 3. Loading gene trees via function "readgtree" #
##################################################
  loghelpers.prntmngr("Loading gene trees", uprFlg=T)
  gTrees = list()
  for (locus in loci) {
    gTreeName = paste(apWd, locus, ".trees", sep="")                    # Actual file names have to match info in xml file
    loghelpers.prntmngr(paste("Loading empirical gene trees of locus '", 
                        locus, "'", sep=""))
    gTrees[[locus]] = readtree.gtree(gTreeName, locus)                  # Loading gene trees via function "readgtree"
    
    # DEBUGLINES:
    #cat("\ngTrees\n"); print(gTrees)
    
    if (verboseBool) {
      cat("\t", "N of gene trees loaded: ", length(gTrees[[locus]]),
          "\n", sep="")
    }
  }

####################################
# 4. Loading regular species trees #
####################################
  loghelpers.prntmngr("Loading species trees", uprFlg=T)
  streeName = paste(apWd, "species.trees", sep="")
  loghelpers.prntmngr("Loading regular species trees")
  # Load sTrees into string (not list!)
  sTrees = readtree.stree(streeName)
  
  # DEBUGLINES:
  #cat("\nsTrees\n"); print(sTrees)

####################################
# 5. Loading phybase species trees #
####################################
  # See if flag descrStats contains "COAL", for which phybase 
  # species trees are needed
  if ("COAL" %in% descrStats) {
    loghelpers.prntmngr("Loading phybase species trees")
    # Load pTrees as list
    pTrees = readtree.phybase(streeName)
    if (verboseBool | debugBool) {  
      if (length(sTrees) != length(pTrees$tree)) {
        cat("\n", xtermStyle::style("ERROR: N(regular species trees) 
                                    != \ N(phybase species trees)", 
                                    fg="red"), "\n", sep="")
      }
      if (length(sTrees) != length(pTrees$tree)) {
        cat("\t","N of species trees loaded: ", length(sTrees), "\n", 
            sep="")
      }
    }
  }
  if (!"COAL" %in% descrStats) {pTrees = NULL}

# CURRENTLY PUT ON HOLD:
# #######################
# # 6. Loading log file #
# #######################
#   loghelpers.prntmngr("Loading logfile", prmFHndl, tofileFlg=F, color="green")
#   log = read.table(paste(apWd, inLogFn, sep=""), header = T)

#####################################################
# 7. Setting up ordered allele-species associations #
#####################################################
  loghelpers.prntmngr("Generating ordered allele-species associations",
                      uprFlg=T)
  assoc = list()
  for (i in 1:length(gTrees)) {
    tips = gTrees[[i]][[1]]$tip.label                                   # Get tips from first tree
    assocD = c()
    for (tip in tips) {
        assocD = rbind(assocD, SAmtrx[which(SAmtrx[,2]==tip),])         # Extract all those allele-species pairs that are in tips
    }
    assoc[[i]] = as.matrix(as.data.frame(assocD))
  }
  # DEBUGLINES:
  #cat("\nassoc\n"); print(assoc)

#####################################
# 8. Reducing outD in debug mode #
#####################################
  if (debugBool) {
    # With 20 input trees, 5% are 2 trees
    if (length(sTrees) > 39) {
      cat("\n", xtermStyle::style("DEBUG> *** REDUCING DATA SET TO 5% ***",
          fg="red"), "\n", sep="")
      reduced = length(sTrees)*0.05
      pTrees = head(pTrees, reduced)
      sTrees = head(sTrees, reduced)
      gTrees = lapply(gTrees, head, reduced)
    }
  }

######################################
# 9. Saving data to variable outD #
######################################
  outD = list()
  outD$asoc$nme = "ordered allele-species associations"
  outD$asoc$dta = assoc
  outD$loci$nme = "loci"
  outD$loci$dta = loci
  outD$gtre$nme = "gene trees"
  outD$gtre$dta = gTrees
  outD$pldy$nme = "ploidy level"
  outD$pldy$dta = ploidy
  outD$ptre$nme = "phybase species trees"
  outD$ptre$dta = pTrees
  outD$stre$nme = "regular species trees"
  outD$stre$dta = sTrees
  # CURRENTLY PUT ON HOLD:
  #     class(outD) = "starbeast.data"
  #     outD$alignm = alignm
  #     outD$log = log
  #     outD$parts = GPmtrx
  
  # DEBUGLINES:
  #cat("\noutD$asoc$dta\n"); print(outD$asoc$dta)
  
  return(outD)
}
