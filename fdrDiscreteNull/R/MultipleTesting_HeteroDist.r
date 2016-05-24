################################################################################
# This file contains the main functions for  Grouped and Weighted Generalized FDR Estimators
# Definitions of modular functions are in companion file "MultipleTesting_HeteroDist_ModFuncs.r"
# Copyrith: Xiongzhi Chen, 2013 -
# Contact: xiongzhi.chen@gmail.com
################################################################################


GeneralizedEstimatorsGrouped = function(data_in = NULL, test_in = NULL, FET_via_in = NULL,
                                          ss=NULL, gpids_in = NULL, cdisp_in = NULL, nopsinestdisp_in = 5,
                                          grpby = c("quantileOfRowTotal","kmeans","divergence"), ngrp_in = NULL,
                                          RefDivergence = NULL,eNetSize = NULL, unif_tol= 10^-3,
                                            FDRlevel_in = 0.05,lambda_in = 0.5,epsilon_in = 1) 
{
  if (test_in == "Exact Negative Binomial Test") {
     if (cdisp_in == "Yes")
       message("^^^Common dispersion will be used for Negative Binomial data")
     else
      message("^^^Heterogenous dispersions will be used for Negative Binomial data")
      }
  m <- nrow(data_in)
  #################  stat of step 1: UNGrouped   ################
  cat("\n","^^^Implement UNGROUPED generalized estimators ...^^^","\n","\n")
  UngroupedResults = GeneralizedFDREstimators(data_in, grp.ids = gpids_in,Test = test_in,
                                              FET_via = FET_via_in, nOptimize = nopsinestdisp_in,
                                            FDRlevel = FDRlevel_in,lambda = lambda_in,epsilon = epsilon_in,CommonDisp = cdisp_in,ss=ss)

  pval_vec_ungrouped = UngroupedResults[[4]]   # pvalue vector
  pval_Supp_ungrped =  UngroupedResults[[5]]   # pvalue supports list 
  deltas_ungrped =  UngroupedResults[[6]]         # 
 
  if (test_in == "Exact Negative Binomial Test")      
          dispest = UngroupedResults[[7]]  # disp or median of dispersions, scalor

  if (test_in == "Exact Negative Binomial Test" & grpby == "divergence" & is.null(eNetSize))
    cat("^^Grouping by divergence but size of epsilon-net NOT specified. Dynamic e-net size will be used","\n")

  ### start of gathering results 
   
       Dis_ugp = UngroupedResults[[1]][[5]]
       Dis_ST_ugp = UngroupedResults[[2]][[5]]
       Dis_BH_ugp = UngroupedResults[[3]][[3]]
  
  
  ### end of gathering results
   
  # pi0 estimated ungrouped by 2 methods
  pi0Est_ugp = UngroupedResults[[1]]$pi0Est;   pi0Est_st_ugp = UngroupedResults[[2]]$pi0Est

  ###### end of step 1: UNGrouped   ##################


  #################  start Step 2: GROUPED ONLY  #################

  #################### start: preparations ##################
  cat("\n","^^ ----- Implement GROUPWISE the procedures ------^^^","\n")

   ngrp = ngrp_in

   grpidx_list = vector("list",ngrp);    grpdata_list = vector("list",ngrp)
   grpres_list = vector("list",ngrp);    wgtPval_list = vector("list",ngrp)
   wgtSTPval_list = vector("list",ngrp)

    pi0Est_gped = double(0);   pi0Est_st_gped = double(0)
    namepiEstgp = double(0);   namepiEst_st_gp = double(0)
     pi0_ov_gp = 0;   pi0st_ov_gp = 0

    grpId_dis = double(0);  grpId_STdis = double(0);  totalCounts = rowSums(data_in)

      # caution ids_grps = sort(unique(grp.ids))
      diffCounts = abs(data_in[,1] - data_in[,2])  # if Poisson

     # data in reality has a differtn structure
     if (test_in == "Exact Negative Binomial Test") {
        nc = ncol(data_in)
       diffCounts = abs(rowSums(data_in[,1:(0.5*nc)]) - rowSums(data_in[,(0.5*nc+1):nc]))     }
     #################### end: preparations ##################

     ######################################################################################
     ############ Section 2: start: actual grouping and groupwise analysis   ################
     ######################################################################################

     # layout goruping strategy first

      if (grpby == "quantileOfRowTotal") {
          cat("^^Grouping by quantiles of row total counts","\n")
          brks = quantile(totalCounts, probs = seq(from=0,to=1,by=1/ngrp), na.rm = TRUE)
          }

     if (grpby == "kmeans") {
         cat("^^Grouping by kmeans on group differences and total counts","\n")
        kmgrps = kmeans(cbind(diffCounts,totalCounts),centers = ngrp)
        }

      ##  reference divergences
     if (grpby == "divergence") {
        
        ## allow reference divergence
        if (RefDivergence == "Yes") { 
             
             cat("^^Grouping by divergence: compute supremum norm w.r.t. Unif ...","\n")
             d_supp = Div_Ref_Unif(pval_Supp_ungrped,test_in)
             du_vec = d_supp[[1]]; ss_vec = d_supp[[2]]
             
             du_max = max(du_vec); du_min = min(du_vec)  #minimal is always zero
             ss_max = max(ss_vec); ss_min = min(ss_vec)
      
             cat("Max and Min RefNorm:", du_max,",", du_min, "; Max and Min supp size:",ss_max,",",ss_min,"\n")
             
            # div_ref = du_vec/du_max + ss_vec/ss_max  # option A # relaxed triangular inequality
             
             div_ref = du_vec/du_max   # option B #relaxed triangular inequality
             
             div_ref_max = max(div_ref) ;  div_ref_min = min(div_ref); diff_div = div_ref_max - div_ref_min
             
             div_brks = seq(from=div_ref_min, to=div_ref_max, by=diff_div/ngrp)
           } # end of reference divergence
         
         ## if not reference divergence  
         if (RefDivergence == "No") {
             
             cat("^^Grouping by pairwises divergences. Be patient ...","\n")

             # insert:  when  |F_i - id| <= m^-3, take it as uniform   
              pv_supp_And_id_appr_unif = Div_Appr_Unif(pval_Supp_ungrped,test_in, unif_tol)  #(pvSpList=NULL,test_used=NULL,appr_unif_tol = NULL) 
              pv_supp_formatted = pv_supp_And_id_appr_unif[[1]]
              id_appr_unif = pv_supp_And_id_appr_unif[[2]]
              
              lgX = length(id_appr_unif)
              if (lgX > 0)  { 
                  cat("^^",lgX,"of p-value cdf's are very close to and will be identified as Unif","\n")
                  
                  if (lgX < m) {
                    id_pvDist_far_unif = (1:m)[-id_appr_unif]
                    ng_div =  ngrp-1
                   # grpidx_list[[1]] = id_appr_unif   # first group via divergence
                    
                    lgXa =  length(id_pvDist_far_unif)
                    pv_supp_far_unif = vector("list",lgXa)
                     for (iX in 1:lgXa)  
                      pv_supp_far_unif[[iX]] = pv_supp_formatted[[iX]]
                    }  
               } else { 
                     cat("^^All p-value cdf's are NOT close enough to Unif","\n") 
                     ng_div =  ngrp 
                     pv_supp_far_unif = pval_Supp_ungrped
                     }
             # end insertion  
             
             ## group those that are far from Unif
             scalfac_in = 1
             Diff_mat = GetDivergenceMatrix(scalfac_in,pv_supp_far_unif)  # change scalfac when needed
             Divs_mat = Diff_mat[[1]]; chi_mat = Diff_mat[[2]]; infNorm_mat = Diff_mat[[3]]; chi_max = Diff_mat[[4]]
    
             chi_vec = as.vector(chi_mat[lower.tri(chi_mat,diag=FALSE)])
             infNorm_vec = as.vector(infNorm_mat[lower.tri(infNorm_mat,diag=FALSE)])
             Divs_vec = as.vector(Divs_mat[lower.tri(Divs_mat,diag=FALSE)])
    
             div_max = max(Divs_vec); div_min = 0 #min(Divs_vec)  #minimal is always zero
    
             cat("Max infNorm is", max(infNorm_vec), "; max ChiSD is",chi_max,"; max Div is",div_max,"\n")
    
             # boxplot may cause memory surge
    
             # two cases for eNet size
             if (is.null(eNetSize))  {
               rad = (div_max - div_min)/(2*ng_div) 
               # make ngrp-1 groups out of those cdf's far from Unif  
               grpidx_list_div = eNetFull(Divs_mat, ng_div, 30, rad, 30) #(data, ngrp, merge, rad, grpsize) 
               }  else {
                 grpidx_list_div = eNetFull(Divs_mat,ng_div, 30, eNetSize, 30)  # last arg, minimal gp size, 2nd von last merge size
               }  # end of eNet
    
            # ngp_div = length(grpidx_list_div)
             cat("^^Number of groups for cdf's far from Unif via divergence is:", ng_div,"\n")
           }  # end if not reference divergence  
         
         } #   end if divergence

     #################### start of grouping and groupwise analysis ####################  
     for (j in 1:ngrp) {
       # start of groupwise analysis by bincounts
        if (grpby == "quantileOfRowTotal") {
           if (j < ngrp)
              grpidx_list[[j]] = which(totalCounts >= brks[j] & totalCounts < brks[j+1])
           if (j == ngrp)
                grpidx_list[[j]] = which(totalCounts >= brks[j] & totalCounts <= brks[j+1])
         }

       # start of groupwise analysis by kmeans
       if (grpby == "kmeans")  grpidx_list[[j]] = which(kmgrps$cluster == j)

        # start of groupwise analysis by e-net divergence
        if ( test_in == "Exact Negative Binomial Test" & grpby == "divergence") {
         
         # if reference divergence 
          if (RefDivergence == "Yes")  {
             if (j < ngrp)
              grpidx_list[[j]] = which(div_ref >= div_brks[j] & div_ref < div_brks[j+1])
             if (j == ngrp)
                grpidx_list[[j]] = which(div_ref >= div_brks[j] & div_ref <= div_brks[j+1])
              }  # end if RefDivergences
          
          # if not reference divergence
             if (RefDivergence == "No") {
               # adjust index of group_lidx    
                if (lgX > 0 & lgX <m) { 
                    if (j == 1)  {
                        grpidx_list[[j]] = id_appr_unif
                      } else {
                       grpidx_list[[j]] = grpidx_list_div[[j-1]]  } 
                     } else {
                    grpidx_list[[j]] = grpidx_list_div[[j]]
                    } #
               }  # end if not reference divergence  
              
         } # end of grouping by divergence
         
       ############ start of  grouped data, groupwise anlysis  #######
       grpdata_list[[j]] = cbind(data_in[grpidx_list[[j]],],grpidx_list[[j]])

       # in order to group, original data is added with a col of total counts
       nc_kep = ncol(data_in)
       cat("\n","^^-----Estimating pi0 for group ",j,"----","\n")

       grpdata_in = grpdata_list[[j]][,1:nc_kep]
      # cat("^^","groupwise data is a matrix",is.matrix(grpdata_in),";", nrow(grpdata_in),"rows in groups",j,"\n")
       CheckIt(grpdata_in)   # check line

       ### when a group of data is non-empty, analyze them
       if (nrow(grpdata_in) == 0) {
          cat("^^^Group",j,"has no data ....","\n")
        } else {
         twopi0est_gp = twoestimators(lambda_in,epsilon_in,pval_vec_ungrouped[grpidx_list[[j]]],deltas_ungrped[grpidx_list[[j]]])
         
          pi0eg = twopi0est_gp[2]; pi0eSTg = twopi0est_gp[1]
          pi0Est_gped = c(pi0Est_gped, pi0eg)
          pi0Est_st_gped = c(pi0Est_st_gped,pi0eSTg)

          namepiEstgp = c(namepiEstgp,paste('pi0Est_gp',j,sep=""))
          namepiEst_st_gp = c(namepiEst_st_gp,paste('pi0Est_st_gp',j,sep=""))
       }   ##  only analyze non-empty grouped data
   } ######## end of loop for (j in 1:ngrp) ######
      
       ####################################################################################
       ######### step 3: start of weighting the p-values from the ungrouped procedure  ############
       ####################################################################################
       # grouping by divergence can create intersecting groups, for which HU 2010 JASA weighting does not work

       ### step 3.1: check if parition via divergence is reached
        non_partition = 0
        for (j1 in 2:ngrp) {
             for (j2 in 1:(j1-1)) {
                 non_partition = non_partition + any(grpidx_list[[j1]] %in% grpidx_list[[j2]])  }  }

       union_chk = length(unlist(grpidx_list))== m

        is_partition = 1- non_partition - !union_chk

        if (is_partition) {
               cat("\n","^^^----- Partition by grouping obtained. Ready for weighting ...","\n")
               }  else {
                cat("^^ ------Partition by grouping NOT obtained. Weighting will NOT be applied...","\n") }

    # if partition reached, then compute weights and following steps
    if (grpby != "divergence" | is_partition) {

         ### step 3.2: get overall pi0est and weight ungrouped pvalues
         for (j in 1:ngrp) {
            # pi0eg = grpres_list[[j]][[1]]$pi0Est; pi0eSTg = grpres_list[[j]][[2]]$pi0Est
              pi0eg = pi0Est_gped[j];    pi0eSTg =  pi0Est_st_gped[j]
             wgt = pi0eg/(1-pi0eg);   wgt_st = pi0eSTg/(1-pi0eSTg)

             if (is.na(wgt) | is.nan(wgt) | is.infinite(wgt))
                   cat("^^Pathological weight--", wgt,"--via gen est for group",j,"\n")
             if (is.na(wgt_st) | is.nan(wgt_st) | is.infinite(wgt_st))
                   cat("^^Pathological weight--", wgt_st,"--via Storey's est for group",j,"\n")

             ###### start of overall pi0est  #############
              pi0_ov_gp = pi0_ov_gp + pi0eg*length(grpidx_list[[j]])/m
              pi0st_ov_gp = pi0st_ov_gp + pi0eSTg*length(grpidx_list[[j]])/m

             ### weigthing pvalues from ungrouped analysis
             pval_vec_ungrouped_gp = pval_vec_ungrouped[grpidx_list[[j]]]
             wgtPval_list[[j]] =  pval_vec_ungrouped_gp*wgt
             wgtSTPval_list[[j]] =   pval_vec_ungrouped_gp*wgt_st

             }  # end: overall pi0est and weight ungrouped pvalues

       ### weighted/grouped storey's and gen est applicable regardless if  overall pi0est = 1
              ov_pi0_est = c(pi0st_ov_gp,pi0_ov_gp)
              if (any(ov_pi0_est >= 1))  cat("^^One of overall estimates of pi0 is one: ",ov_pi0_est,"\n")

              wgtPval_vec = unlist(wgtPval_list);  wgtSTPval_vec = unlist(wgtSTPval_list)

              if (any(is.na(wgtSTPval_vec)) | any(is.nan(wgtSTPval_vec)))  {
                cat("^^Undefined weighted pvalue via Storey est found. Filtering them out ..","\n")
                cat("^^ --- UNCOMPENSATED Filtering of NA or NAN p-values effects GROUPED AND WEIGHTED procedure ..","\n")

                wgtPval_vec_filtered = wgtPval_vec[!is.na(wgtPval_vec) & !is.nan(wgtPval_vec)]
              }   else {  wgtPval_vec_filtered = wgtPval_vec}


             if (any(is.na(wgtSTPval_vec)) | any(is.nan(wgtSTPval_vec))) {
                 cat("^^Undefined weighted pvalue via Storey est found. Filtering them out ..","\n")
                 cat("^^ --- UNCOMPENSATED Filtering of NA or NAN p-values effects GROUPED AND WEIGHTED procedure ..","\n")

                 wgtSTPval_vec_filtered = wgtSTPval_vec[!is.na(wgtSTPval_vec) & !is.nan(wgtSTPval_vec)]
             }   else { wgtSTPval_vec_filtered = wgtSTPval_vec }


       ######### step 3.4: start: apply grouped and weighted scheme to BH if needed #####
         ### Now on grouped and weighted BH
         if ( pi0_ov_gp >= 1)   {
           cat("^^Overall pi0 est by Gen Est is", pi0_ov_gp,". BH weighted and grouped is zero..", "\n")
            TDP_BHwg = 0;   FDP_BHwg = 0;   Dis_BH_gw = double(0) 
            }  else {
                 cat("^^Implement grouped BH by MEANINGFUL weighting via Gen Est of pi0 ...","\n")
                 FDRlevel_ada = FDRlevel_in/(1- pi0_ov_gp)
                 BHadaWG <- BHFDRApp(wgtPval_vec_filtered,FDRlevel_ada)
              
              # if not simulation
               Dis_BH_gw = BHadaWG[[1]][,2]   

            }  ## end weighted grouped BH by Gen est

          # weighted grouping BH via storey's
         if ( pi0st_ov_gp >= 1)  {
            cat("^^Overall pi0 est by Storey is", pi0st_ov_gp ,". BH weighted and grouped is zero..", "\n")
            TDP_BHwgST = 0;   FDP_BHwgST = 0;    Dis_BH_STgw = double(0) 
            } else {
               cat("^^Implement grouped BH by MEANINGFUL weighting via Storey's Est of pi0 ...","\n")
               FDRlevel_STada = FDRlevel_in/(1- pi0st_ov_gp)
               BHadaWG_ST <- BHFDRApp(wgtSTPval_vec_filtered,FDRlevel_STada)

                 #### start: results from grouping and weighting of BH via storey's #####
               # if not simulation
                Dis_BH_STgw = BHadaWG_ST[[1]][,2] 

            } # end else    weighted grouping via storey's

    } # end if partition, then compute weights

  #################### step 4: start: collection all results #####################
  # if partition not reached, overall pi0 est is undefined
  if (grpby == "divergence" & non_partition != 0) {
      pi0_ov_gp = NA;  pi0st_ov_gp = NA }

  pi0Es = c(pi0Est_ugp,pi0Est_st_ugp,pi0_ov_gp,pi0st_ov_gp,pi0Est_gped,pi0Est_st_gped)
  nmpi0es = c("pi0E_Gen","pi0E_Storey","pi0E_gwGen","pi0E_gwStorey",namepiEstgp,namepiEst_st_gp)

 # if grouped by diverence, all overall quanities via weighting in HU 2010 are undefined

## if it is not a simulation study, then no way to know TDP and FDP
 
    if (grpby == "divergence" & non_partition != 0) {
       Dis_BH_STgw = Dis_BH_gw = NA   }

  results_grp = list(Dis_ugp,Dis_ST_ugp,Dis_BH_ugp,Dis_BH_gw,Dis_BH_STgw,pval_vec_ungrouped,pval_Supp_ungrped,deltas_ungrped,pi0Es)
  nmDis = c("Gen","Storey","BH","GWGen","GWStorey","pval_ungrped","pSupp_ungrped","deltas_ungrped","pi0ests")
  names(results_grp) = nmDis 
  
  cat("\n","^^----Finished gathering grouped/weighted analysis results---^^","\n")
  
  if (test_in == "Binomial Test" | test_in == "Fisher's Exact Test")
    return(results_grp)
    
  if (test_in == "Exact Negative Binomial Test")  {
   results_grp_disp = list(Dis_ugp,Dis_ST_ugp,Dis_BH_ugp,Dis_BH_gw,Dis_BH_STgw,pval_vec_ungrouped,pval_Supp_ungrped,deltas_ungrped,
                                  pi0Es,dispest)
   names(results_grp_disp) = c(nmDis,"DispEst")
   return(results_grp_disp)  }
      

} # end of main function
