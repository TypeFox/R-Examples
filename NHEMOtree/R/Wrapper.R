Wrapper <-
function(data, CostMatrix, 
         gens=50, popsize=50, 
         ngens=14, bound=10^-10, 
         init_prob=0.80,
         crossover_prob=0.5, mutation_prob=0.5,
         CV=5, ...){

  # Miscellaneous
    Ziel          <- names(data)[1]   # Grouping variable
    Variablennamen<- names(data)[-1]  # Explanatory variables
    Chromo_length <- ncol(data)-1     # Amount of explanatory variables
    n_child_total <- popsize
    
  # Costs
    MaxCosts        <- sum(CostMatrix[,2])
    MaxCosts_rounded<- ceiling(MaxCosts/10)*10

  # Initialization
    Ind<- Init_Ctree(Daten=data, KOSTEN=CostMatrix, popsize=popsize, 
                     in_percent=init_prob, Chromo_length=Chromo_length, 
                     CV.Lauf=CV, proteins=Variablennamen)

  # Convergence detection
    F1<- F2<- F3<- c()
    S_Metrik    <- c()
    for (i in 1:popsize){
         F1[i]<- Ind[[i]]$Misclassification
         F2[i]<- Ind[[i]]$Costs
         F3[i]<- Ind[[i]]$Anz
    }
    Fitness    <- cbind(F1, F2, F3)
    S_Metrik[1]<- dominated_hypervolume(t(Fitness[,1:2]), c(100, MaxCosts))


#############
# Main loop #
#############
pWert1<- 1  # p value of the tests for the preceding n generations
pWert2<- 1  # p value of the tests for the latest n generations, i.e. pWert2[g]=pWert1[g]
g     <- 1

###############
# Pruned tree #
###############
  while (g <= gens && pWert1 > 0.05&& pWert2 > 0.05){

  # Population
    Vectors<- NULL
    Misclassification<- c(); Costs<- c(); Anz<- c()
    for (i in 1:popsize){
         Vectors             <- rbind(Vectors, Ind[[i]]$Vector)
         Misclassification[i]<- Ind[[i]]$Misclassification
         Costs[i]            <- Ind[[i]]$Costs 
         Anz[i]              <- Ind[[i]]$Anz
    }
    Fitness<- rbind(Misclassification, Costs, Anz)

    j<- 1
    while(j <= n_child_total){

    #%#%#%#%#%#%#%#%#%#%
    # Parent selection #
    #%#%#%#%#%#%#%#%#%#%
      ELTERN<- BinaryTS_FINANZEN_oD(N_Parents=popsize, Fitness=Fitness, Vectors=Vectors)

    #%#%#%#%#%#%#
    # Crossover # 
    #%#%#%#%#%#%#
      Temp<- OnePoint_X_1child(Eltern_Chromos=ELTERN[[1]], crossover_prob=crossover_prob)

    #%#%#%#%#%#%#
    # Mutations #
    #%#%#%#%#%#%#
      temp            <- Mutation3(Chromosom=Temp, m=mutation_prob)
      Chromo_child_mut<- temp$Chromosom_mutiert
      Gene_mut        <- temp$Genes_flipped

    #%#%#%#%#%#%#%#%#%#%#%#%
    # Fitness calculations #
    #%#%#%#%#%#%#%#%#%#%#%#%
      Tupel_Init_Child<- Fitness_Calc_Ctree_Sim(Daten=data, CV.Lauf=CV, 
                                                Chromosom=Chromo_child_mut, proteins=Variablennamen)
      VectorOut_P     <- Chromo_bit_long(Varis=Tupel_Init_Child$P_Tree, Gesamtvariablen=Variablennamen)

    # Pruned trees without any variables
      if (sum(VectorOut_P[[1]])>0){

        # Misclassification rate of pruned variables
          FKR_P<- Fitness_Calc_Ctree_Sim(Daten=data, CV.Lauf=CV, 
                                         Chromosom=VectorOut_P[[1]], proteins=Variablennamen)$Misclassification
  
        # Costs of the pruned tree
          CostsOut_P    <- CostMatrix[VectorOut_P[[1]], 2]
          CostsOut_P_sum<- sum(CostsOut_P)
          CostsOut_P_Anz<- length(CostsOut_P)

        # New individuals saved in 'Ind'
          Individuum                  <- list()
          Individuum$Vector           <- VectorOut_P[[1]]
          Individuum$Misclassification<- FKR_P
          Individuum$Costs            <- CostsOut_P_sum
          Individuum$Anz              <- CostsOut_P_Anz
          Ind[[popsize+j]]            <- Individuum

        j<- j+1
        } 
    }

  #%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
  # Environmental selection - Fast-non-dominated-sort #
  #%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
    Vectors<- NULL
    Misclassification<- c(); Costs<- c(); Anz<- c()
    for (l in 1:length(Ind)){
         Vectors             <- rbind(Vectors, Ind[[l]]$Vector)
         Misclassification[l]<- Ind[[l]]$Misclassification
         Costs[l]            <- Ind[[l]]$Costs 
    }
    Fitness_Pruned<- rbind(Misclassification, Costs)

  # Sorting by ranks
    temp_P               <- rbind(nds_rank(Fitness_Pruned), 1:length(Ind), Fitness_Pruned)
    ii                   <- order(temp_P[1,], temp_P[3,], temp_P[4,], temp_P[2,])
    Fitness_sorted_Pruned<- temp_P[,ii]           
    cd_rank              <- Fitness_sorted_Pruned[1,popsize]
    which_cd_rank        <- which(Fitness_sorted_Pruned[1,]<=cd_rank)
    n_toomany            <- length(which_cd_rank)-popsize

  # Crowding distance NOT required
    if (n_toomany == 0){
        Fitness_sorted_Pruned<- Fitness_sorted_Pruned[,which(Fitness_sorted_Pruned[1,]<=cd_rank)]
    }

  # Crowding distance required
    if (n_toomany > 0){
        FS_after_NDS<- Fitness_sorted_Pruned[,which_cd_rank]

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
            Spannweite<- Crowding_Set[(m+2),ncol(Crowding_Set)]-Crowding_Set[(m+2),1]
            for (i in 2:(ncol(Crowding_Set)-1)){
                 Crowding_Set[nrow(Crowding_Set),i]<- Crowding_Set[nrow(Crowding_Set),i]+(Crowding_Set[(m+2),(i+1)]-Crowding_Set[(m+2),(i-1)])/Spannweite
            }
            Crowding_Set[nrow(Crowding_Set), which.min(Crowding_Set[(m+2),])]<- Inf
            Crowding_Set[nrow(Crowding_Set), which.max(Crowding_Set[(m+2),])]<- Inf

          Distance_overall<- Distance_overall + Crowding_Set[nrow(Crowding_Set),]
          }
 
          Crowding_Set     <- rbind(Crowding_Set, Distance_overall)
          iii              <- order(Crowding_Set[nrow(Crowding_Set),], decreasing = TRUE)
          Crowding_Set     <- Crowding_Set[,iii]
          Crowding_Set_kept<- Crowding_Set[1:4,(1:n_keep)]
      }

  # New population
    Fitness_sorted_Pruned<- cbind(Fitness_sorted_Pruned[,which(Fitness_sorted_Pruned[1,]<cd_rank)], Crowding_Set_kept)
  }

  # New individuals in 'Ind'
    Parents_ID<- Fitness_sorted_Pruned[2,]
    Temp      <- list()
    for (o in 1:length(Parents_ID)) Temp[[o]]<- Ind[[Parents_ID[o]]]
    Ind<- list(); Ind<- Temp

  ######################################
  # S metric for convergence detection #
  ######################################
    Temp         <- t(Fitness_sorted_Pruned[3:4,])
    Temp2        <- cbind(Temp, rep(1,nrow(Temp)))
    S_Metrik[g+1]<- dominated_hypervolume(t(Temp2[,1:2]), c(100, MaxCosts))

    # Variance test
      if (g >= ngens){   
          pWert1<- pWert2                                                 
          pWert2<- Chi_Test(PI=S_Metrik[(g-ngens+2):(g+1)], Limit=bound)
      } 
g<- g+1
}

  ##########
  # Output #
  ##########
  # Final population
    FP_ges                        <- cbind(1:length(Ind), Temp) 
    ii                            <- order(FP_ges[,2], FP_ges[,3], decreasing = F)
    FP_ges                        <- FP_ges[ii,]
    Output                        <- list()
    Output$S_Metric               <- S_Metrik[length(S_Metrik)]
    Output$Misclassification_total<- range(Temp[,1])
    Output$Costs_total            <- range(Temp[,2])

  # Tree with lowest misclassification rate
    Best_Tree<- Best_Tree_new(data=data, KOSTEN=CostMatrix, Ind=Ind, CV=CV)

  # Output
    Out       <- c(Output, list(Trees=Ind, Best_Tree=Best_Tree, 
                   S_Metrik_temp=S_Metrik, MaxCosts_rounded=MaxCosts_rounded, method="Wrapper"))
    class(Out)<- "NHEMOtree"
    return(Out)
}
