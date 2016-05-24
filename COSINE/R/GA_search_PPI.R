GA_search_PPI <-
function(lambda,scaled_node_score,scaled_edge_score,PPI,
num_iter=1000, muCh=0.05, zToR=10, minsize=10){


    ## define the objective scoring function for condition specific subnetwork 
    all_genes<-names(scaled_node_score)

    subset_score<-function(sub){    
      genes<-all_genes[sub==1]
      n<-length(genes)    
      if(n<minsize){return(10000)}
      else{
          node_score<-sum(scaled_node_score[genes])/sqrt(n)
          edges<- PPI[,1] %in% genes & PPI[,2] %in% genes
          m<-sum(edges)
          if(m==0)total_score<- -10000
          if(m>0){
          edge_score<-sum(scaled_edge_score[edges])/sqrt(m)
          total_score<- lambda*edge_score + (1-lambda)*node_score
          }
          return (-total_score)
      }
     }

     monitor <- function(obj) {
       minEval = min(obj$evaluations);
       filter = obj$evaluations == minEval;
    #   print(table(filter))
       bestObjectCount = sum(rep(1, obj$popSize)[filter]);
       if (bestObjectCount > 1) {
          bestSolution = obj$population[filter,][1,];
       } else {
       bestSolution = obj$population[filter,];
       }

     outputBest = paste(obj$iter, " #selected=", sum(bestSolution),
     " Best (Score=", -minEval, "):\n", sep="");
     print(outputBest)
     }

 
    ##Start search

     gene_num <- length(scaled_node_score)
 #   num_lam <- length(lambda)
 #   num_gene_selected <- rep(0,num_lam)
 #   best_score <- rep(0,num_lam)
 #   optimal_subnet <- vector(length=num_lam,mode="list")
 #   GA_result <- vector(length=num_lam,mode="list")
 #   for(i in 1:num_lam){
     print(paste("Working on lambda=",lambda)) 
     GA_result <- rbga.bin(size=gene_num,evalFunc=subset_score,iters=num_iter,mutationChance=muCh,monitorFunc=monitor,zeroToOneRatio=zToR)      
     a <- which.min(GA_result$evaluations)
     final <- GA_result$population[a,]
     b <- which(final==1)
     num_gene_selected <- length(b)
     optimal_subnet <- b
     best_score <- (-1)*min(GA_result$evaluations)
     print(paste("Finished lambda=",lambda))
#   }
        
    return(list( Subnet_size = num_gene_selected, Best_Scores = best_score, Subnet = optimal_subnet, GA_obj = GA_result))

}

