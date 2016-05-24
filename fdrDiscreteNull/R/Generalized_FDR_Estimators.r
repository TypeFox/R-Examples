#########################################################################
# Generalized Estimator of Proportion of True Nulls
# This file contains the main functions for Generalized FDR Estimators
# Definitions of modular functions are in companion file "Generalized_FDR_Estimators_ModularFuncs.r"
# Copyrith: Xiongzhi Chen, 2013 -
# Contact: xiongzhi.chen@gmail.com
####################################################################

### Start of main function  
GeneralizedFDREstimators = function(data=NULL, Test=c("Binomial Test", "Fisher's Exact Test","Exact Negative Binomial Test"),
                                    FET_via = c("PulledMarginals","IndividualMarginals"), ss=NULL, grp.ids = NULL,
                                    CommonDisp = NULL, nOptimize = 5, FDRlevel=0.05, lambda=0.5, epsilon=1) 
  { # start of main function
  
  #########################################################################
  # Section 1: check if arguments are correctly provided
  ####################################################################
  
  if (is.null(data))      stop("^^Please provide data.")
  # check test type
   chktest = identical(Test,"Binomial Test") + identical(Test,"Fisher's Exact Test") + identical(Test,"Exact Negative Binomial Test")
  if (chktest ==0 )      stop("^^Unsupported type of test statistics.")    
  
  # For Fisher's exact test, two column data is ok
  if (Test == "Fisher's Exact Test")   {
   if (is.null(FET_via)) stop("^^Please specify type of marginals for Fisher's exact test","\n")
   else {
        if (FET_via == "PulledMarginals" & ncol(data) != 2)          stop("Please reformat data into an m-by-2 array.")
        if (FET_via == "IndividualMarginals" & ncol(data) != 4)      stop("Please reformat data into an m-by-4 array.")
        } 
   }
  
  # For Binomial Test, two column data is ok      
  if (Test == "Binomial Test" & ncol(data) != 2)      stop("^^Please reformat data into an m-by-2 array.")
  
  if (is.null(FDRlevel))       stop("^^Please specify FDR level.")
  if (is.null(lambda) | is.null(epsilon))   stop("^^Please input tuning parameters.")  
   
  # get m, number of rows of data
  m = nrow(data);   n = ncol(data)
  
  ######################################################################### ##########
  # Section 2: Get p-values and their supports from Binomial Test 
  #################################################################### ########## ##########
  if (Test == "Binomial Test") {      
       cat("^^Using binomial test.","\n")
       # check for zero rows 
       if (any(rowSums(data) == 0)) stop("^^Please remove zero rows.","\n")
        
       sizePoi = data[,1] + data[,2]
       marginalsAndsizes = cbind(data[,1], sizePoi)
       # same structure as under
       cat("^^Computing supports of discrete p-value distributions ...","\n")
       simpvaluesupports  <- apply(marginalsAndsizes,1,pvalueByBinoSupport)
       simpvalues <- unlist(sapply(simpvaluesupports,'[',2))
             
        cat("^^Computing deviations of p-value distributions from the uniform...","\n") 
        mnew = length(simpvalues);         simdeltas <- double(mnew) 
        for (idxA in 1:mnew)   
            simdeltas[idxA] <- deviationsExtPvalSupp(lambda,simpvaluesupports[[idxA]])          
     } # end of case 1, if (Test == "Binomial Test")  
  
  ######################################################################### ##########
  # Section 3: Get p-values and their supports from Fisher's Exact Test 
  #################################################################### ########## ##########
  # different test will be used as per method chosen
  if (Test == "Fisher's Exact Test") {
  
        cat("^^Using Fisher's Exact Test.","\n")
        # check for zero rows
        if (any(rowSums(data) == 0)) stop("^^Please remove zero rows in data.","\n")
        
        if (FET_via == "PulledMarginals" & ncol(data) == 2) {  
           cat("^^Fisher's Exact Test is based on pulled marginal","\n")
          ## fulltable retruns a 3D array
          full_tables_data <- fulltable(data)
          # get cell counts and marginals for Fisher's exact test
          countveccontrol = data[,1];    countvectreat = data[,2]
          cellcountsmarginals <- getcellcountsandmarginals(countveccontrol,countvectreat)
          }
         
        if (FET_via == "IndividualMarginals" & ncol(data) == 4) {
            cat("^^Fisher's Exact Test is based on individual marginals","\n")
             cellcountsmarginals = getcellcountsandmarginals_DE(data)
           }
           
        # assign cell counts and marginals
        simallcellcounts <- cellcountsmarginals[[1]]
        simallmarginals <- cellcountsmarginals[[2]]
        
        # compute two-sides pvalues and their supports
        cat("^^Computing supports of discrete p-value distributions ...","\n")    
        simpvalues <- unlist(lapply(simallcellcounts, FUN = function(x) fisher.test(x)$p.value))
        simpvaluesupports  <- lapply(simallmarginals,pvalueSupport)
         
        cat("^^Computing deviations of p-value distributions from the uniform...","\n")
        mnew = length(simpvalues);         simdeltas <- double(mnew)  
        for (idxA in 1:mnew)    
           simdeltas[idxA] <- deviations(lambda,simpvaluesupports[[idxA]])
     }       # end of case 2, if (Test == "Fisher's Exact Test")

  ######################################################################### ####################
  # Section 4: Get p-values and their supports from Exact Negative Binomial Test 
  ####################################################################  ####################
  # for negative binomial, the size have to be estimated
  if (Test == "Exact Negative Binomial Test")  {  
        cat("^^Exact negative binomial test is used; normalization is needed.","\n")
       # check for zero rows 
       if (any(rowSums(data[,1:ss]) == 0) | any(rowSums(data[,(ss+1):(2*ss)]) == 0))
             stop("^^Please remove zero rows in data.","\n")         
            
        # group identifications
        ids_grps = sort(unique(grp.ids))
        # force it only allows two groups
        if (length(ids_grps) > 2) stop("Caution: more than 2 groups found.")
            
        # column indices for each group
        colIdx_grp1 = which(ids_grps[1] == grp.ids)
        colIdx_grp2 = which(ids_grps[2] == grp.ids)
        
        # balanced design required
        if (length(colIdx_grp1) != length(colIdx_grp2))
          stop("^^^Please provide data with equal sample size across two groups","\n")
                           
        # check if all counts for a group is zero
        if (any(rowSums(data[,colIdx_grp1]) ==0) | any(rowSums(data[,colIdx_grp2]) ==0)) 
         stop("^^Please remove rows for which total counts for a group is zero.","\n")
       
        #### normalizing by edgeR
        cat("^^Start normalizing by edgeR","\n") 
        cat("^^edgeR: preperating the DGElist object","\n")
        datatmp = DGEList_adjusted(counts = data, group=rep(ids_grps,each=length(grp.ids)/length(ids_grps)))
       # print(datatmp)
      
           if (CommonDisp == "Yes") {
              cat("^^edgeR: estimating common dispersion","\n")
              datatmp = estimateCommonDisp_adjusted(datatmp,nopsinestdisp = nOptimize,verbose=TRUE)  
                # args default nopsinestdisp = NULL, tol = 1e-06, rowsum.filter = 5
                  
              estdisp =  datatmp$common.dispersion    # save estimated dispersion
              if (is.null(estdisp)) stop("^^edgeR: failed to estimate common dispersion.","\n")
              
              marginalsAndsizes = CountsNormalizerInEdgeR(datatmp, pair = 1:2, 
                                                 dispersion = estdisp, prior.count.total = 0.5)
            } # end if (CommonDisp == "Yes")
           
            if (CommonDisp == "No") {
                datatmp = estimateCommonDisp_adjusted(datatmp, nopsinestdisp = nOptimize,verbose=TRUE) 
                  # estimateCommonDisp_adjusted replaces estimateCommonDisp
                cat("^^edgeR: start estimating tagwise dispersions","\n")
                datatmp = estimateTagwiseDisp(datatmp)  
                
                estdisp_tagwise =  datatmp$tagwise.dispersion 
                estdisp = median(estdisp_tagwise, na.rm = TRUE)   # get and save the median
        
                cat("^^edgeR: finished estimating tagwise dispersions","\n")
                if (any(is.null(estdisp_tagwise)))
                   stop("^^edgeR: failed to estimate tagwise dispersions.","\n")
                
                 marginalsAndsizes = CountsNormalizerInEdgeR(datatmp, pair = 1:2, 
                                                   dispersion = estdisp_tagwise, prior.count.total = 0.5)
              }  # if (CommonDisp == "No")
        
        cat("^^edgeR: normalization completed.","\n")  
        ### end of normalization by edgeR

        ################## compute p values, pvalue distributions, pvalue supports   ##################
        #   pvalueByNegativeBinoSupportApp uses parametrization in Di et al, whereas pvalueByNegativeBinoSupport 
        # uses another parametrization
        cat("^^Computing supports of discrete p-value distributions ...","\n")

        nrMS = nrow(marginalsAndsizes);         simpvaluesupports = vector("list",nrMS)
        for (i3 in 1:nrMS) { simpvaluesupports[[i3]] = pvalueByNegativeBinoSupportApp(marginalsAndsizes[i3,])  }
        
        # extract p values from the 2nd entry
        simpvalues <- unlist(sapply(simpvaluesupports,'[',2))
        if (any(is.na(simpvalues))) cat("^^NA in pvalues found ...","\n")
           if (any(is.nan(simpvalues))) cat("^^NAN in pvalues found ...","\n")
           
        cat("^^Computing deviations of discrete p-value distributions from the uniform...","\n")
        mnew = length(simpvalues);   simdeltas <- double(mnew)      
        for (idxA in 1:mnew)    
           simdeltas[idxA] <- deviationsExtPvalSupp(lambda,simpvaluesupports[[idxA]]) 
                     
     }     # end of case 3, if (Test == "Exact Negative Binomial Test")

   
   #########################################################################
  # Section 5: Estimate pi0 and FDP from p-values and their supports 
  #################################################################### 
      
   cat("^^Estimating pi0...","\n")
   twopisAtDefPair <- twoestimators(lambda,epsilon,simpvalues,simdeltas)
   
  ### FDP estimators
  cat("^^Implementing Generalized Estimators and Storey's estimators ...","\n")
  discoveriesAtDefPair <- StAndChenFDREstimatorApp(simpvalues,FDRlevel,twopisAtDefPair)
    
  ############ apply vanilla BH procedure    ###########
   cat("^^Implementing BH procedure ...","\n")
  BHrejAndTresh <- BHFDRApp(simpvalues,FDRlevel) 
   
  ##################################################################
  # Section 7: put estimates in a dataframe and return it 
  ####################################################################
  cat("^^Gathering analysis results","\n")
  # names for entries: thresh, FDP, number of discoveries, indices of discoveries
  GenEst = list("pi0Est" = twopisAtDefPair[2], 
                    "Threshold" = discoveriesAtDefPair[[3]][2], "FDPEst" = discoveriesAtDefPair[[3]][4], 
                    "NumberOfDiscoveries" = length(discoveriesAtDefPair[[2]]),"IndicesOfDiscoveries" = discoveriesAtDefPair[[2]])
  StoreyEst = list("pi0Est" = twopisAtDefPair[1],
                     "Threshold" = discoveriesAtDefPair[[3]][1], "FDPEst" =discoveriesAtDefPair[[3]][3],
                      "NumberOfDiscoveries" = length(discoveriesAtDefPair[[1]]),"IndicesOfDiscoveries" = discoveriesAtDefPair[[1]])
  BHEst = list("Threshold" = BHrejAndTresh[[2]],
                  # "pi0Est" = 1, "NominalFDR" = FDRlevel,
                   "NumberOfDiscoveries" = length(BHrejAndTresh[[1]][,2]),"IndicesOfDiscoveries" = BHrejAndTresh[[1]][,2])
  
  AnalysisResults = list("Generalized_Estimator" = GenEst,"Storey_Estimator" = StoreyEst,"BH_Procedure" = BHEst,"pvalues" = simpvalues,
                         "pvalSupp" = simpvaluesupports,"Deltas" = simdeltas) 
 
 cat("\n","^^---Finished implementing Generalized Estimators, Storey's and BH procedures ---^^^","\n")
 
 # rerturn results 
 if (Test == "Binomial Test" | Test == "Fisher's Exact Test")
       return(AnalysisResults)
    
 if (Test == "Exact Negative Binomial Test")  {
      anaRes = list("GeneralizedEstimator" = GenEst,"StoreyEstimator" = StoreyEst,"BHProcedure" = BHEst,
                     "pvalues" = simpvalues,"pvalSupp" = simpvaluesupports,"Deltas" = simdeltas,"DispEst" = estdisp)  
      return(anaRes)  }
      
} # end of main function