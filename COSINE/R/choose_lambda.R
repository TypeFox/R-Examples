choose_lambda <-
function(diff_expr,diff_coex,lambda,subnet_size,num_random_sampling,best_score){

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



    #Random sampling

    gene_num<-length(diff_expr)
    random_score<-matrix(0,nrow=5,ncol=num_random_sampling)

 
    for(i in 1:5){
    
      for(j in 1:num_random_sampling){
        lam<-lambda[i]
        random_net<-rep(0,gene_num)
        a<-sample(1:gene_num,subnet_size[i])
        random_net[a]<-1
        random_score[i,j]<- -subset_score(random_net)
        print(j)
      } 
    }

    mean<-apply(random_score,1,mean)
    sd<-apply(random_score,1,sd)
    adjusted_score<-(best_score-mean)/sd
    best_lambda <- lambda[which.max(adjusted_score)]
  
    return (list(Adj_score=adjusted_score, Best_lam=best_lambda, Random_Score=random_score))
}

