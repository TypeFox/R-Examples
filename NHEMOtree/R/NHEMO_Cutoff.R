NHEMO_Cutoff <-
function(data, CostMatrix, 
         gens=50, popsize=50, max_nodes=10, 
         ngens=14, bound=10^-10, 
         init_prob=0.8, 
         ps=c("tournament", "roulette", "winkler"), tournament_size=4, 
         crossover=c("standard", "brood", "poli"), brood_size=4, 
         crossover_prob=0.5, mutation_prob=0.5, 
         CV=5, vim=0,
         ncutoffs=10, opt=c("gini", "mcr")){
  
  # Checking
  ps       <- match.arg(ps)
  crossover<- match.arg(crossover)
  opt      <- match.arg(opt)
  
  # Miscellaneous
  N_Varis      <- ncol(data)-1
  temp         <- data[,2:ncol(data)]
  uS           <- min(temp)             # Lower limit for cutoffs
  oS           <- max(temp)             # Upper limit for cutoffs
  n_child_total<- popsize
  temp_mc_count<- 0
  
  # Costs
  MaxCosts        <- sum(CostMatrix[,2])
  MaxCosts_rounded<- ceiling(MaxCosts/10)*10
  
  # Variable importance measures
  vim1<- vim2<- vim3<- vim4<- vim5<- vim6<- rep(0, N_Varis)  
  vim1_Gen<- vim2_Gen<- vim3_Gen<- vim4_Gen<- vim5_Gen<- vim6_Gen<- matrix(nrow=gens, ncol=N_Varis)
  Weigths_lin<- Gewichte_lin(Max_Baumtiefe=max_nodes)
  Weigths_exp<- Gewichte_exp(Max_Baumtiefe=max_nodes)
  
  # Permutated data
  Daten_Perm<- Datenpermutation_Sim(Daten=data, Vari_N=N_Varis)
  
  # Initialization
  Initialisierung<- list()  
  if (opt=="mcr"){
    Initialisierung<- Init_RHaH_opt_Sim_FKR(Daten=data, N_Varis=N_Varis, Proteinkosten=CostMatrix, 
                                            Popsize=popsize, Max_Knoten=max_nodes, CV.Lauf=CV,
                                            Init_Wsk=init_prob, uS=uS, oS=oS,
                                            AnzCO=ncutoffs)
  }
  if (opt=="gini"){
    Initialisierung<- Init_RHaH_opt_Sim_Gini(Daten=data, N_Varis=N_Varis, Proteinkosten=CostMatrix, 
                                             Popsize=popsize, Max_Knoten=max_nodes, CV.Lauf=CV,
                                             Init_Wsk=init_prob, uS=uS, oS=oS,
                                             AnzCO=ncutoffs)
  }
  Baum           <- Initialisierung$Baum
  Fitness        <- cbind(Initialisierung$Rate, Initialisierung$Finanzen, Initialisierung$Anzahl)
  S_Metrik       <- c()
  S_Metrik[1]    <- dominated_hypervolume(t(Fitness[,1:2]), c(100, MaxCosts))
  Quantile_Knoten<- quantile(Initialisierung$Knoten, probs = c(0, 0.25, 0.5, 0.75, 1))
  Quantile_FKR   <- quantile(Initialisierung$Rate, probs = c(0, 0.25, 0.5, 0.75, 1))
  
  
  ##########################
  # Main loop of NHEMOtree #
  ##########################
  pWert1<- 1  # p value of the tests for the preceding n generations
  pWert2<- 1  # p value of the tests for the latest n generations, i.e. pWert2[g]=pWert1[g]
  g     <- 1
  
  while (g <= gens && pWert1 > 0.05 && pWert2 > 0.05){
    
    Nachkomme<- list()
    Rate<- c(); Finanzen<- c(); Anzahl<- c()
    j<- 1
    
    while(j <= n_child_total){
      #%#%#%#%#%#%#%#%#%#%
      # Parent selection #
      #%#%#%#%#%#%#%#%#%#%
      ELTERN<- c()
      ELTERN<- Parentselection(Type=ps, Tree=Baum, K=tournament_size)
      
      #%#%#%#%#%#%#
      # Crossover # 
      #%#%#%#%#%#%#
      # Selection of correct vim
      if (vim==0) VW_in<- rep(1, N_Varis) 
      if (vim==1) VW_in<- vim1_Gen[g-1,]     
      if (vim==2) VW_in<- vim2_Gen[g-1,]  
      if (vim==3) VW_in<- vim3_Gen[g-1,]  
      if (vim==4) VW_in<- vim4_Gen[g-1,]  
      if (vim==5) VW_in<- vim5_Gen[g-1,]
      if (vim==6) VW_in<- vim6_Gen[g-1,]
      if (g==1)   VW_in<- 0
      
      Nachkomme[[j]]<- Rekombination_VIM3(Daten=data, Proteinkosten=CostMatrix, 
                                          Eltern_Baum1=Baum[[ELTERN[1]]], Eltern_Baum2=Baum[[ELTERN[2]]],
                                          CV_Laeufe=CV, X_Wsk=crossover_prob, 
                                          X_Art=crossover, Brutgroesse=brood_size, 
                                          VIM=VW_in,
                                          Max_Knoten=max_nodes, N_Varis=N_Varis, uS=uS, oS=oS)
      
      #%#%#%#%#%#%#
      # Mutations #
      #%#%#%#%#%#%#
      # Tree structure
      Nachkomme[[j]]<- Mutation_overall(Tree=Nachkomme[[j]], N_Varis=N_Varis, 
                                        uS=uS, oS=oS, Max_Knoten=max_nodes, 
                                        Mutationswsk=mutation_prob, 
                                        Init_Wsk=init_prob)
      
      
      #%#%#%#%#%#%#%#%#%#
      # Optimal Cutoffs #
      #%#%#%#%#%#%#%#%#%#
      if (opt=="mcr"){
        Nachkomme[[j]]<- Cutoffwahl_fkr_Sim(Tree=Nachkomme[[j]], Daten=data, CV.Lauf=CV, AnzCO=ncutoffs)
      } else{
        Nachkomme[[j]]<- Cutoffwahl_gini_Sim(Tree=Nachkomme[[j]], Daten=data, CV.Lauf=CV, AnzCO=ncutoffs)
      }
      
      
      #%#%#%#%#%#%#%#%#%#%#%#%#
      # Deleting empty leaves #
      #%#%#%#%#%#%#%#%#%#%#%#%#
      Nachkomme[[j]]<- Delete_EL(Tree=Nachkomme[[j]], Daten=data, N_Varis=N_Varis, uS=uS, oS=oS)
      
      
      #%#%#%#%#%#%#%#%#%#%#%#%
      # Fitness calculations #
      #%#%#%#%#%#%#%#%#%#%#%#%
      P_N        <- Nachkomme[[j]]$Varis[,2]-1
      temp_N     <- NoDupKey(P_N)
      Rate[j]    <- 100*FKR_CV(Tree=Nachkomme[[j]], Daten=data, CV_Laeufe=CV)
      Finanzen[j]<- getCosts_Tree_sum_Sim(Kostenmatrix=CostMatrix, Protein_IDs=temp_N)
      Anzahl[j]  <- length(temp_N)
      
      Nachkomme[[j]]$FKR     <- Rate[j]
      Nachkomme[[j]]$Finanzen<- Finanzen[j]
      Nachkomme[[j]]$Anzahl  <- Anzahl[j]
      
      j<- j+1
      temp_rep<- NULL; temp_crossover<- NULL; temp_mb<- NULL; mb_sample<- NULL; temp_mkc<- NULL
    }
    
    #%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
    # Environmental selection - Fast-non-dominated-sort #
    #%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
    Baum_neu<- c(Baum, Nachkomme)
    FKR     <- c()
    Finanzen<- c()
    for (fv in 1:length(Baum_neu)){
      FKR[fv]     <- Baum_neu[[fv]]$FKR
      Finanzen[fv]<- Baum_neu[[fv]]$Finanzen
    }
    
    # Sorting by ranks
    Temp          <- rbind(1:length(Baum_neu), FKR, Finanzen)
    temp_O        <- rbind(nds_rank(Temp[2:3,]), Temp)
    ii            <- order(temp_O[1,], temp_O[3,], temp_O[4,], temp_O[2,])
    Fitness_sorted<- temp_O[,ii]      
    cd_rank       <- Fitness_sorted[1,popsize]
    which_cd_rank <- which(Fitness_sorted[1,]<=cd_rank)
    n_toomany     <- length(which_cd_rank)-popsize
    
    # Crowding distance NOT required
    if (n_toomany == 0){
      Fitness_sorted<- Fitness_sorted[,which(Fitness_sorted[1,]<=cd_rank)]
    }
    
    # Crowding distance required
    if (n_toomany > 0){
      FS_after_NDS<- Fitness_sorted[,which_cd_rank]
      
      # Calculation of the crowding distance
      Crowding_Set<- FS_after_NDS[,which(FS_after_NDS[1,]==cd_rank)]
      n_keep      <- ncol(Crowding_Set)-n_toomany
      
      if (ncol(Crowding_Set)==2){
        temp<- runif(1,0,1)
        if (temp <  0.5) Crowding_Set_kept<- Crowding_Set[,1]
        if (temp >= 0.5) Crowding_Set_kept<- Crowding_Set[,2]
      }
      
      if (ncol(Crowding_Set) > 2){
        Distance_overall<- matrix(rep(0, ncol(Crowding_Set)), ncol=1)
        for (m in 1:2){
          # Sorting
          ii              <- order(Crowding_Set[(m+2),])
          Crowding_Set    <- Crowding_Set[,ii]
          Distance_overall<- t(as.matrix(Distance_overall))[,ii]
          Distance        <- rep(0, ncol(Crowding_Set))
          Crowding_Set    <- rbind(Crowding_Set, Distance)
          
          # Distance
          Spannweite<- Crowding_Set[(m+2), ncol(Crowding_Set)]-Crowding_Set[(m+2),1]
          for (n in 2:(ncol(Crowding_Set)-1)){
            Crowding_Set[nrow(Crowding_Set),n]<- Crowding_Set[nrow(Crowding_Set),n]+(Crowding_Set[(m+2),(n+1)]-Crowding_Set[(m+2),(n-1)])/Spannweite
          }
          Crowding_Set[nrow(Crowding_Set), which.min(Crowding_Set[(m+2),])]<- Inf
          Crowding_Set[nrow(Crowding_Set), which.max(Crowding_Set[(m+2),])]<- Inf
          
          Distance_overall<- Distance_overall + Crowding_Set[nrow(Crowding_Set),]
        } 
        
        Crowding_Set     <- rbind(Crowding_Set, Distance_overall)
        iii              <- order(Crowding_Set[nrow(Crowding_Set),], decreasing = TRUE)
        Crowding_Set     <- Crowding_Set[,iii]
        Crowding_Set_kept<- Crowding_Set[1:(2+2),(1:n_keep)]
      }
      
      # New population
      Fitness_sorted<- cbind(Fitness_sorted[,which(Fitness_sorted[1,]<cd_rank)], Crowding_Set_kept)
    }
    
    # New trees in 'Baum'                    
    Baeume<- Fitness_sorted[2,]
    Baum  <- list()
    for (b in 1:length(Baeume)) Baum[[b]]<- Werte_Reduktion_ZHT2(Tree=Baum_neu[[Baeume[b]]])
    
    ######################################
    # S metric for convergence detection #
    ######################################
    Temp         <- t(Fitness_sorted[3:4,])
    Temp2        <- cbind(Temp, rep(1,nrow(Temp)))
    S_Metrik[g+1]<- dominated_hypervolume(t(Temp2[,1:2]), c(100, MaxCosts))
    
    # Variance test
    if (g >= ngens){
      pWert1<- pWert2
      pWert2<- Chi_Test(PI=S_Metrik[(g-ngens+2):(g+1)], Limit=bound)
    } 
    
    ########################
    # VIM for generation g #
    ########################
    for (b2 in 1:popsize){
      # Simple absolute frequency
      vim1<- VIM1(Tree=Baum[[b2]], VIM_alt=vim1)
      
      # Simple relative frequency
      vim2<- VIM2(Tree=Baum[[b2]], VIM_alt=vim2) 
      
      # Relative frequency
      vim3<- VIM3(Tree=Baum[[b2]], VIM_alt=vim3)
      
      # Linear weigthed relative frequency
      vim4<- VIM4_complete(Tree=Baum[[b2]], Weights=Weigths_lin, VIM_alt=vim4)
      
      # Exponentially weigthed relative frequency
      vim5<- VIM4_complete(Tree=Baum[[b2]], Weights=Weigths_exp, VIM_alt=vim5)
      
      # Permutation accuracy
      vim6<- PF_1perm(Tree=Baum[[b2]], VIM_alt=vim6, Daten_Perm=Daten_Perm, Daten=data, CV_Laeufe=CV)
    }
    
    # VIMs
    vim1_Gen[g,]<- round((vim1/(popsize*g)*100),2)
    vim2_Gen[g,]<- round((vim2/(popsize*g)*100),2)
    vim3_Gen[g,]<- round((vim3/(popsize*g)*100),2)
    vim4_Gen[g,]<- round((vim4/(popsize*g)*100),2)
    vim5_Gen[g,]<- round((vim5/(popsize*g)*100),2)
    vim6_Gen[g,]<- round((vim6/(popsize*g)*100),2)
    
    g<- g+1
  } 
  
  ##########
  # Output #
  ##########
  # Final population
    Misclassification             <- sapply(Baum, "[[", "FKR")
    Costs                         <- sapply(Baum, "[[", "Finanzen")
    FP_ges                        <- cbind(1:length(Baum), Misclassification, Costs) 
    ii                            <- order(FP_ges[,2], FP_ges[,3], decreasing = F)
    FP_ges                        <- FP_ges[ii,]
    S_ref                         <- dominated_hypervolume(t(FP_ges[,2:3]), c(100, MaxCosts))
    Output                        <- list()
    Output$S_Metric               <- S_ref
    Output$Misclassification_total<- range(FP_ges[,2])
    Output$Costs_total            <- range(FP_ges[,3])
  
  # Variable importance measures
    VIMS<- rbind(vim1_Gen[g-1,], vim2_Gen[g-1,], vim3_Gen[g-1,], vim4_Gen[g-1,], vim5_Gen[g-1,], vim6_Gen[g-1,])
  
  # Rename trees
    Misclassification<- sapply(Baum, "[[", "FKR")
    Costs            <- sapply(Baum, "[[", "Finanzen")
    for (i in 1:length(Baum)){
      Baum[[i]]$Misclassification<- Misclassification[i]
      Baum[[i]]$Costs<-             Costs[i]
    }
  
  # Output
    Out       <- c(Output, list(VIMS=VIMS, Trees=Baum, S_Metrik_temp = S_Metrik, 
                   MaxCosts_rounded=MaxCosts_rounded, method="NHEMO_Cutoff"))
    class(Out)<- "NHEMOtree"
    return(Out)
}
