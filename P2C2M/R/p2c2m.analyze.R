p2c2m.analyze <-
function (inFn, inData, prmFHndl) {
  # Descr:  core of script, coordinates the metric calculation
  # Deps:   (various)
  # I/p:    inFn
  #         inData
  #         prmFHndl

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  descrStats = get("P2C2M_flg_dscrStat", envir=P2C2M_globalVars)
  mpiBool = get("P2C2M_flg_mpiBool", envir=P2C2M_globalVars)
  nReps = get("P2C2M_flg_nReps", envir=P2C2M_globalVars)
  singleAllele = get("P2C2M_flg_snglAl", envir=P2C2M_globalVars)
  srtBool = get("P2C2M_flg_srtBool", envir=P2C2M_globalVars)
  verboseBool = get("P2C2M_flg_vrbBool", envir=P2C2M_globalVars)

  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> core.analyze", fg="red"), sep="")
  }

##########################
# 1. Preparatory actions #
##########################
  n_loci = length(inData$loci$dta)
  n_sTrees = length(inData$stre$dta)
  # generate empty matrices of dimensions 'matrix'
  post_dist = post_pred_dist = matrix(nrow=n_sTrees, ncol=n_loci)
  # label the columns of the empty matrices by the locus names
  colnames(post_dist) = colnames(post_pred_dist) = inData$loci$dta
  # Initilizing outdata 
  outD = list()

##########################################################
# 2. Metrics for the trees of the posterior distribution #
##########################################################
  loghelpers.prntmngr(paste("Metrics for the trees of the", 
                            "posterior distribution"), uprFlg=T)

 ##########################################################
 # 2.a. Looping through genes, calc. descript. statistics #
 ##########################################################
  # Loop through all gene tree names
  for (j in 1:n_loci) {
    loghelpers.prntmngr(paste("Locus '", inData$loci$dta[j], "'",
                              sep=""))
    loghelpers.prntmngr("Analyzing gene trees:\n", nwlFlg=F)
    # Loop through all species trees
    for (i in 1:n_sTrees) {
      # Display tree number on screen
      loghelpers.prntmngr(paste(" ", i, sep=""), nwlFlg=F)
      # Save tree number as replicate ID
      assign("P2C2M_flg_repID", paste("post_dist", i, sep="."), 
             envir = P2C2M_globalVars)

    # DEBUGLINES:
    #cat("\ninData$gtre$dta[[j]][[i]]\n"); print(inData$gtre$dta[[j]][[i]])
    #cat("\ninData$ptre$dta$tree[[i]]\n"); print(inData$ptre$dta$tree[[i]])
    #cat("\ninData$ptre$dta$names\n"); print(inData$ptre$dta$names)
    #cat("\ninData$stre$dta[[i]]\n"); print(inData$stre$dta[[i]])
    #cat("\ninData$asoc$dta\n"); print(inData$asoc$dta)
    #cat("\ninData$pldy$dta[j]\n"); print(inData$pldy$dta[j])

      post_dist[i,j] = corehelpers.metrics(inData$gtre$dta[[j]][[i]],
                                     inData$ptre$dta$tree[[i]],
                                     inData$ptre$dta$names,
                                     inData$stre$dta[[i]],
                                     inData$asoc$dta[[j]],
                                     inData$pldy$dta[j],
                                     descrStats,
                                     singleAllele)
    }
    loghelpers.prntmngr("\n", nwlFlg=F)
  }

 ###########################
 # 2.b. Assembling results #
 ###########################
  for (s in descrStats) {
    # Populating list "outD" with the results of the trees of the posterior distribution
    outD[[s]]$post_dist$unsort = as.data.frame(rtnlc(post_dist, which(s==descrStats)))
    # Assigning colnames
    colnames(outD[[s]]$post_dist$unsort) = inData$loci$dta
    # Generating a sorted result clone
    outD[[s]]$post_dist$sorted = apply(outD[[s]]$post_dist$unsort, MARGIN=2, sort)
  }

#####################################################################
# 3. Metrics for the trees of the posterior predictive distribution #
#####################################################################
  loghelpers.prntmngr(paste("Metrics for the trees of the", 
                            "post. predictive distribution"), uprFlg=T)

 ##################################
 # 3.a. Initiate MPI (if applic.) #
 ##################################
  # Determining number of processors (i.e. CPUs) available
  if (Sys.info()['sysname'] == "Linux") {
    nproc = as.numeric(system("nproc", intern=T))
  }
  # MacOS option
  if (Sys.info()['sysname'] == "Darwin") {
    nproc = as.numeric(system("sysctl -n hw.ncpu", intern=T))
  }

  if (mpiBool && nproc >= 3) {
    loghelpers.prntmngr(paste("\n", "Initialize parallelization", "|", 
                              "N of processors:", toString(nproc)))
    # Specifying one master and nproc-1 slaves
    #Rmpi:: mpi.spawn.Rslaves(nslaves=nproc-1, needlog=F)
    Rmpi:: mpi.spawn.Rslaves(nslaves=nproc-1)
    # Passing separate functions to MPI environment
    mpiEnvir = c(calc.lcwt,
                 calc.gsi,
                 calc.ndc,
                 calc.parse,
                 calc.coal,
                 calchelpers.brprob,
                 calchelpers.descend,
                 calchelpers.dmvparse,
                 calchelpers.gtreeparse,
                 calchelpers.nodetips,
                 calchelpers.probcalc,
                 corehelpers.metrics,
                 # LEGACY: corehelpers.repl, 
                 frmtMntse,
                 # LEGACY: apTreeshape::as.treeshape.phylo,
                 ape::branching.times,
                 ape::drop.tip,
                 genealogicalSorting::gsi,
                 phybase::loglikeSP)
    lapply(mpiEnvir, Rmpi::mpi.bcast.Robj2slave)
  }

 ##############################
 # 3.b. Looping through genes #
 ##############################
  simTrees_across_loci = list()
  for (j in 1:n_loci) {
    loghelpers.prntmngr(paste("Locus '", inData$loci$dta[j], "'",
                              sep=""))
    loghelpers.prntmngr("Analyzing the posterior predictive trees:\n",
                        nwlFlg=F)
    n_tips = length(inData$gtre$dta[[j]][[1]]$tip.label)
    simTrees_across_sTrees = list()
    for (i in 1:n_sTrees) {
      # Display tree number on screen
      loghelpers.prntmngr(paste(" ", i, sep=""), nwlFlg=F)
      # Save tree number as replicate ID
      assign("P2C2M_flg_repID", paste("post_pred_dist", j, i, sep="."), 
             envir = P2C2M_globalVars)

 ############################
 # 3.c. Simulation of trees #
 ############################
      # print(inData$assoc)
      # speciesName   alleleName
      # [1,] "A"      "a1"
      # [2,] "A"      "a10"
      # [3,] "A"      "a11"
      # [4,] "A"      "a12"
      # [5,] "A"      "a2"

      # A constituent-frequency table is required for the simulations in ms,
      # because ms saves the alleles as numbers such as "1", "2", "3", ...
      CFTable = readhelpers.makeCFtable(inData$asoc$dta[[j]])

      # print(CFTable)
      #       aList  V2
      # Var1  "A"    "1"
      # Var1  "A"    "2"
      # Var1  "A"    "3"

      # Simulate gene trees
      simTrees = ms.exec(CFTable, inData$stre$dta[[i]], n_tips, 
                         inData$pldy$dta[j], prmFHndl)
      class(simTrees) = "multiPhylo"

      ## Possible improvement in TFL: why not "length(simTrees)" instead of 
      ## "as.numeric(nReps)"
      simReps = matrix(nrow=as.numeric(nReps), ncol=1)

 ##################################################################
 # 3.d. Calculation of descript. statistics depend. on MPI status #
 ##################################################################
      if (mpiBool && nproc >= 3) {
        # Applying function "corehelpers.metrics" in parallel
        #  The second element in mpi.parLapply has to be the function that
        #  the data is to be applied
        simReps = Rmpi::mpi.parLapply(simTrees,
                                      corehelpers.metrics,
                                      inData$ptre$dta$tree[[i]],
                                      inData$ptre$dta$names,
                                      inData$stre$dta[[i]],
                                      CFTable,
                                      inData$pldy$dta[j],
                                      descrStats,
                                      singleAllele)
        simReps = as.matrix(as.character(simReps))
      }

      else {
        for (k in 1:nReps) {
        simReps[k,1] = corehelpers.metrics(simTrees[[k]],
                                           inData$ptre$dta$tree[[i]],
                                           inData$ptre$dta$names,
                                           inData$stre$dta[[i]],
                                           CFTable,
                                           inData$pldy$dta[j],
                                           descrStats,
                                           singleAllele)
        }
      }

 #########################################
 # 3.e. Parsing of average metric values #
 #########################################
      # Extracting the listcolumns of the matrix simReps via GSO.rtnlc
      sumD = list()
      for (s in descrStats) {
        sumD[[s]] = Mode(rtnlc(simReps, which(s==descrStats)))
      }
      sumD = unlist(sumD)
      names(sumD) = NULL
      
      # Save the replicate values to variable "post_pred_dist"
      post_pred_dist[i,j] = toString(sumD)

 #############################################
 # 3.f. Preparing simulated trees for saving #
 #############################################

      # Loop over simulated replicates
      ## Possible improvement: Use apply function instead of loop
      # Generate set of print-ready simtrees
      simTrees_wLabels = list()
      for (k in 1:nReps) {
        # LEGA!CY: # Generating a simTree version where tips have been replaced with the correct allele names
        # LEGA!CY: assoc_ids = unlist(lapply(inData$asoc$dta[[j]][,2], substring, first=2))
        # LEGA!CY: # TFL removes leading zeros
        # LEGA!CY: assoc_ids = gsub("^[0]+", "", assoc_ids)
        # LEGA!CY: # TFL removes all non-numeric characters
        # LEGA!CY: oldLbls = gsub("[^0-9]", "", simTrees[[k]]$tip.label)

        # LEGA!CY: # DEBUGLINES: 
        # LEGA!CY: #cat("\nassoc_ids\n"); print(assoc_ids)
        # LEGA!CY: #cat("\nsimTrees[[k]]$tip.label\n"); print(simTrees[[k]]$tip.label)
        # LEGA!CY: #cat("\noldLbls\n"); print(oldLbls)
        
        # LEGA!CY: match_indcs = match(oldLbls, assoc_ids)
        # LEGA!CY: # Error handling
        # LEGA!CY: if (NA %in% match_indcs | "NA" %in% match_indcs) {
        # LEGA!CY:   stop(cat("\nERROR: Cannot label the terminals of the posterior predictive trees correctly.\n"))
        # LEGA!CY: }
        # LEGA!CY: nwLbls = inData$asoc$dta[[j]][,2][match_indcs]
        # LEGA!CY: simTrees_wLabels[[k]] = simTrees[[k]]
        # LEGA!CY: simTrees_wLabels[[k]]$tip.label = mapply(toString, nwLbls)
        simLbls = gsub("[^0-9]", "", simTrees[[k]]$tip.label)
        nwLbls = inData$asoc$dta[[j]][as.integer(simLbls),2]
        simTrees_wLabels[[k]] = simTrees[[k]]
        simTrees_wLabels[[k]]$tip.label = mapply(toString, nwLbls)       
      class(simTrees_wLabels) = "multiPhylo"
      }
      simTrees_across_sTrees[[paste("sTree_", i, sep="")]] = simTrees_wLabels


    # End of looping through sTrees
    }
    simTrees_across_loci[[inData$loci$dta[[j]]]] = simTrees_across_sTrees
    loghelpers.prntmngr("\n", nwlFlg=F)

  # End of looping through loci
  }

 ##################
 # 3.g. Close MPI #
 ##################
  if (mpiBool && nproc >= 3) {
    # Closing slaves
    loghelpers.prntmngr(paste("\n", "Close parallelization", "\n",
                              sep=""))
    Rmpi::mpi.close.Rslaves(dellog=T)
  }
  
 #############################
 # 3.g. Save simulated trees #
 #############################
  outD$simTrees = simTrees_across_loci 

 ###########################
 # 3.h. Assembling results #
 ###########################
  for (s in descrStats) {
    # Populating list "outD" with the simulated results
    # Both "as.matrix" and "as.data.frame" are critical below
    outD[[s]]$post_pred_dist$unsort = as.matrix(as.data.frame(rtnlc(post_pred_dist, which(s==descrStats))))
    # Assigning colnames
    colnames(outD[[s]]$post_pred_dist$unsort) = inData$loci$dta
    # Generating a sorted result clone
    outD[[s]]$post_pred_dist$sorted = apply(outD[[s]]$post_pred_dist$unsort, MARGIN=2, sort)
  }

#####################################
# 4. Calculating metric differences #
#####################################
  loghelpers.prntmngr(paste("Calculating differences btw. the",
                            "posterior and the posterior predictive", 
                            "distribution"))
  
  for (s in descrStats) {
    # Populating list "outD" with the differences
    if (srtBool) {
      outD[[s]]$dif = statshelpers.diffrnce(outD[[s]]$post_dist$sorted, 
                                            outD[[s]]$post_pred_dist$sorted)
      # Remove the outD set not chosen
      outD[[s]]$post_dist$unsort = NULL
      outD[[s]]$post_pred_dist$unsort = NULL
    }
    else {
      outD[[s]]$dif = statshelpers.diffrnce(outD[[s]]$post_dist$unsort, 
                                            outD[[s]]$post_pred_dist$unsort)
      # Remove the outD set not chosen
      outD[[s]]$post_dist$sorted = NULL
      outD[[s]]$post_pred_dist$sorted = NULL
    }
  }

  return(outD)
}
