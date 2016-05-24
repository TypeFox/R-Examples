PSOSearch <- 
  function(pts,class,MaxOrder=3,Population=1000,Iteration=100, 
           c1=2,c2=2,TopSNP=10,measure=1,alpha=0){
    
    ################################################################
    #
    ## PSO based method for computing main effects of selected SNPs and interaction
    ## effects of selected SNP-combinations within the maximum order.
    #
    # input
    #    pts: Row -> Sample, Column -> SNP
    #         1 -> AA
    #         2 -> Aa
    #         3 -> aa
    #    class: Row -> 1, Column -> class label
    #         1 -> case
    #         2 -> control
    #    MaxOrder: the specified maximum order, must be setted as 1,2,3,4 or 5.
    #             By default, 3.
    #    Population: numeric. The number of particles.
    #              For example, Population=1000
    #    Iteration: numeric. The number of iterations.
    #             For example, Iteration=100
    #    c1: numeric. The acceleration factor of individual experience.
    #              For example, c1=2
    #    c2: numeric. The acceleration factor of global experience.
    #            For example, c2=2
    #    TopSNP: numeric. The selected SNPs with top indexes
    #            For example, TopSNP=10
    #    measure: the label of current used evaluation measure
    #            1 -> The classic co-information measure
    #            2 -> The Normalized co-information measure
    #            3 -> TingHu's Co-Information Measure
    #    alpha: the lower threshold of effects, either main effects or interaction
    #           effects, must be higher or equal to 0, and by default, 0.
    #
    # output
    #    SingleEffect: There are 2 columns. The first column saves all SNPs, and the 
    #                  second column saves their corresponding effects. Descending 
    #                  save according to their effects.                   
    #    TwoEffect: There are 3 columns. The first 2 columns save all 2-SNP combinations,
    #               and the last column saves their corresponding effects. Descending 
    #               save according to their interaction effects.
    #    ThreeEffect: There are 4 columns. The first 3 columns save all 3-SNP combinations,
    #                 and the last column saves their corresponding effects. Descending 
    #                 save according to their interaction effects.
    #    FourEffect: There are 5 columns. The first 4 columns save all 4-SNP combinations,
    #                and the last column saves their corresponding effects. Descending 
    #                save according to their interaction effects.
    #    FiveEffect: There are 6 columns. The first 5 columns save all 5-SNP combinations,
    #                and the last column saves their corresponding effects. Descending
    #                save according to their interaction effects.
    #
    ## Junliang Shang
    ## shangjunliang110@163.com
    ## 10.10/2014
    #
    ################################################################
    
    # Input Data
    
    cat("  PSO Initialization ...\n")
    # maximum particle position
    x_max <- ncol(pts)
    # minimum particle position
    x_min <- 1
    
    # maximum particle velocity
    v_max <- x_max-x_min
    # minimum particle velocity
    v_min <- -v_max
    
    SingleEffect <- array(0,dim=c(1,1,1))
    TwoEffect <- array(0,dim=c(1,1))
    ThreeEffect <- array(0,dim=c(1,1))
    FourEffect <- array(0,dim=c(1,1))
    FiveEffect <- array(0,dim=c(1,1))
    
    # exhaustive search 1-order SNP
    A <- combn(1:x_max,1)
    Effect <- sapply(A,EvaluationMeasure,pts,class,measure)
    Value <- sapply(A,function(x,Effect) Effect[[x]],Effect)
    SingleEffect <- t(rbind(A,Value))
    SingleEffect[,2] <- abs(SingleEffect[,2])
    dim(SingleEffect) <- c(length(SingleEffect)/2,2)
    colnames(SingleEffect) <- c("SNP","Value")
    
    SingleEffect <- SingleEffect[order(SingleEffect[,2],decreasing=T),]
    
    # PSO search high-order SNP combinations
    if (MaxOrder>1){
      # generate order
      order <- round(runif(Population, 2, MaxOrder))
      
      # initialize particle position
      position <- list()
      for (i in 1:Population){
        position[i] <- list(sample(x_min:x_max,order[i],replace=FALSE))
      }
      
      # initialize particle velocity
      velocity <- list()
      for (i in 1:Population){
        velocity[i] <- list(runif(order[i],v_min,v_max))
      }
      
      # initialize individual optimal solution
      pbest <- list()
      pbest_value <- numeric()
      for (i in 1:Population){
        pbest[i] <- list(sample(x_min:x_max,order[i],replace=FALSE))
        Effect <- EvaluationMeasure(pbest[[i]],pts,class,measure)
        pbest_value[i] <- Effect$Value
      }
      
      # initialize global optimal solution
      gbest <- list()
      gbest_value <- numeric()
      for (i in 1:(MaxOrder-1)){
        gbest[i] <- list(sample(x_min:x_max,i+1,replace=FALSE))
        Effect <- EvaluationMeasure(gbest[[i]],pts,class,measure)
        gbest_value[i] <- Effect$Value
      }
      
      # iteration
      cat("  PSO Iteration...\n")
      
      index <- numeric(x_max)
      
      for (i in 1:Iteration){
        
        # display current iteration
        cat("  Iteration: ", i, "\n") 
        
        #
        for (j in 1:Population){
          for (k in 1:order[j]){
            index[pbest[[j]][k]] <- index[pbest[[j]][k]]+1
          }
        }
        
        # update particle velocity
        for (j in 1:Population){
          for (k in 1:order[j]){
            velocity[[j]][k] <- (max(index)-index[position[[j]][k]])/(max(index)-min(index))*velocity[[j]][k]+
              c1*runif(1,0,1)*(position[[j]][k]-pbest[[j]][k])+
              c2*runif(1,0,1)*(position[[j]][k]-gbest[[order[j]-1]][k])
            if (velocity[[j]][k]>v_max){
              velocity[[j]][k]=runif(1,v_min,v_max)
            }
            if (velocity[[j]][k]<v_min){
              velocity[[j]][k]=runif(1,v_min,v_max)
            }
          }
        }
        
        # update particle position
        for (j in 1:Population){
          for (k in 1:order[j]){
            position[[j]][k] <- position[[j]][k]+round(velocity[[j]][k])
            if (position[[j]][k]>x_max){
              position[[j]][k]=sample(x_min:x_max,1,replace=FALSE)
            }
            if (position[[j]][k]<x_min){
              position[[j]][k]=sample(x_min:x_max,1,replace=FALSE)
            }
            while (length(position[[j]])!=length(unique(position[[j]]))){
              position[[j]][k]=sample(x_min:x_max,1,replace=FALSE)
            }
          }  
        }
        
        # update individual optimal solution
        for (j in 1:Population){
          Effect <- EvaluationMeasure(position[[j]],pts,class,measure)
          tmp <- Effect$Value
          if (tmp>pbest_value[[j]]){
            pbest[j] <- position[j]
            pbest_value[j] <- tmp
          }
          
          newposition <- position[j]
          for (k in 1:order[j]){
            newposition[[1]][k] <- x_max+x_min-position[[j]][k]
          }
          
          Effect <- EvaluationMeasure(newposition[[1]],pts,class,measure)
          tmp <- Effect$Value
          if (tmp>pbest_value[[j]]){
            pbest[j] <- newposition[1]
            pbest_value[j] <- tmp
          }
        }
        
        # update global optimal solution
        for (j in 1:Population){
          tmp <- length(pbest[[j]])
          if (pbest_value[[j]]>gbest_value[[tmp-1]]){
            gbest[tmp-1] <- pbest[j]
            gbest_value[tmp-1] <- pbest_value[j]
          }
        } 
      }
      
      
      if (TopSNP>ncol(pts)){
        TopSNP <- ncol(pts)
      }
      
      reTest <- order(index,decreasing=TRUE)[1:TopSNP]
      NewFacterNum <- 0;
      for (i in 2:MaxOrder){
        NewFacterNum <- NewFacterNum+choose(TopSNP,i)
        order <- c(order,rep(i,choose(TopSNP,i)))
      }
      
      postag <- Population+1
      
      for (i in 2:MaxOrder){
        A <- combn(reTest,i)
        for (j in 1:ncol(A)){
          pbest[postag] <- list(A[,j])
          Effect <- EvaluationMeasure(A[,j],pts,class,measure)
          tmp <- Effect$Value
          pbest_value[postag] <- tmp
          postag <- postag+1
        }
      }
      
      for (i in 1:length(pbest)){
        pbest[[i]] <- sort(pbest[[i]])
      }  
      
      if (alpha<0) alpha <- 0  
      
      if (MaxOrder>=2){
        Two <- which(order==2)
        TwoEffect <- array(0,dim=c(length(Two),3))
        for (i in 1:length(Two)){
          TwoEffect[i,1:2] <- pbest[[Two[i]]]
          TwoEffect[i,3] <- pbest_value[Two[i]]
        }
        TwoEffect <- TwoEffect[TwoEffect[,3]>alpha]
        dim(TwoEffect) <- c(length(TwoEffect)/3,3)
        TwoEffect <- TwoEffect[!duplicated(TwoEffect[,1:2])]
        dim(TwoEffect) <- c(length(TwoEffect)/3,3)
        colnames(TwoEffect) <- c("SNP1","SNP2","Value")
        
        TwoEffect <- TwoEffect[order(TwoEffect[,3],decreasing=T),]
      }
      
      if (MaxOrder>=3){
        Three <- which(order==3)
        ThreeEffect <- array(0,dim=c(length(Three),4))
        for (i in 1:length(Three)){
          ThreeEffect[i,1:3] <- pbest[[Three[i]]]
          ThreeEffect[i,4] <- pbest_value[Three[i]]
        }
        ThreeEffect <- ThreeEffect[ThreeEffect[,4]>alpha]
        dim(ThreeEffect) <- c(length(ThreeEffect)/4,4)
        ThreeEffect <- ThreeEffect[!duplicated(ThreeEffect[,1:3])]
        dim(ThreeEffect) <- c(length(ThreeEffect)/4,4)
        colnames(ThreeEffect) <- c("SNP1","SNP2","SNP3","Value")
        
        ThreeEffect <- ThreeEffect[order(ThreeEffect[,4],decreasing=T),]
      }
      
      if (MaxOrder>=4){
        Four <- which(order==4)
        FourEffect <- array(0,dim=c(length(Four),5))
        for (i in 1:length(Four)){
          FourEffect[i,1:4] <- pbest[[Four[i]]]
          FourEffect[i,5] <- pbest_value[Four[i]]
        }
        FourEffect <- FourEffect[FourEffect[,5]>alpha]
        dim(FourEffect) <- c(length(FourEffect)/5,5)
        FourEffect <- FourEffect[!duplicated(FourEffect[,1:4])]
        dim(FourEffect) <- c(length(FourEffect)/5,5)
        colnames(FourEffect) <- c("SNP1","SNP2","SNP3","SNP4","Value")
        
        FourEffect <- FourEffect[order(FourEffect[,5],decreasing=T),]
      }
      
      if (MaxOrder>=5){
        Five <- which(order==5)
        FiveEffect <- array(0,dim=c(length(Five),6))
        for (i in 1:length(Five)){
          FiveEffect[i,1:5] <- pbest[[Five[i]]]
          FiveEffect[i,6] <- pbest_value[Five[i]]
        }
        FiveEffect <- FiveEffect[FiveEffect[,6]>alpha]
        dim(FiveEffect) <- c(length(FiveEffect)/6,6)
        FiveEffect <- FiveEffect[!duplicated(FiveEffect[,1:5])]
        dim(FiveEffect) <- c(length(FiveEffect)/6,6)
        colnames(FiveEffect) <- c("SNP1","SNP2","SNP3","SNP4","SNP5","Value")
        
        FiveEffect <- FiveEffect[order(FiveEffect[,6],decreasing=T),]
      }   
    }
    
    #############################
    # Return results
    #############################
    list(SingleEffect=SingleEffect,TwoEffect=TwoEffect,ThreeEffect=ThreeEffect,
         FourEffect=FourEffect,FiveEffect=FiveEffect)   
  }