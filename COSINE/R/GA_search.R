GA_search <-
function(lambda,diff_expr,diff_coex, num_iter=1000, muCh=0.05, zToR=10){


## define the objective scoring function for condition specific subnetwork 

    my.fun<-function(vector){
       return(diff_coex[vector[1],vector[2]])
       }


    subset_score<-function(vector){    
      selected_subset<-which(vector==1)
      n<-length(selected_subset)    
      if(n<2){return(0)}
      else{
          node_score<-sum(diff_expr[selected_subset])/sqrt(n)
          edges<-combn(selected_subset,2)
          edge_score<-sum(apply(edges,2,my.fun))/sqrt(choose(n,2))
          total_score<- lam*edge_score + (1-lam)*node_score
       return (-total_score)
      }
     }


     monitor <- function(obj) {
     minEval = min(obj$evaluations);
     filter = obj$evaluations == minEval;
     bestObjectCount = sum(rep(1, obj$popSize)[filter]);
     if (bestObjectCount > 1) {
     bestSolution = obj$population[filter,][1,];
     } else {
     bestSolution = obj$population[filter,];
     }

     outputBest = paste(obj$iter, " #selected=", sum(bestSolution),
     " Best (Score=", -minEval, "):\n", sep="");
     print(outputBest)
     print(which(bestSolution==1))
     }



    ##Start search

    gene_num <- length(diff_expr)
    num_lam <- length(lambda)
    num_gene_selected <- rep(0,num_lam)
    best_score <- rep(0,num_lam)
    optimal_subnet <- vector(length=num_lam,mode="list")
    GA_result <- vector(length=num_lam,mode="list")
    for(i in 1:num_lam){
      lam <- lambda[i]
      print(paste("Working on lambda=",lam))
      GA_result[[i]] <- rbga.bin(size=gene_num,evalFunc=subset_score,iters=num_iter,mutationChance=muCh,monitorFunc=monitor,zeroToOneRatio=zToR)
      #plot(GA_result[[i]])
      a <- which.min(GA_result[[i]]$evaluations)
      final <- GA_result[[i]]$population[a,]
      b <- which(final==1)
      num_gene_selected[i] <- length(b)
      optimal_subnet[[i]] <- b
      best_score[i] <- (-1)*min(GA_result[[i]]$evaluations)
      print(paste("Finished lambda=",lam))
    }

    return(list( Subnet_size = num_gene_selected, Best_Scores = best_score, Subnet = optimal_subnet, GA_obj = GA_result))

}

