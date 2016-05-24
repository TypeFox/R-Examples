DataSimu <-
function(){


cat("Generating covariance matrix......")

set1_matrix<-matrix(0,500,500)
for(i in 1:500){
  for(j in 1:500){
    if(i==j)set1_matrix[i,j]=1
  }
}

set2_matrix<-matrix(0,500,500)
for(i in 1:500){
   for(j in 1:500){
      if(i==j)set2_matrix[i,j]=1
      if(i<=50 & j<=50 & i!=j)set2_matrix[i,j]=0.6
   }
}


set3_matrix<-set1_matrix

set4_matrix<-set2_matrix


set5_matrix<-matrix(0,500,500)
for(i in 1:500){
   for(j in 1:500){
      if(i==j)set5_matrix[i,j]=1
      if(i<=50 & j<=50 & i!=j){
         if((i<=25 & j<=25)|(i>25 & j>25) )set5_matrix[i,j]=0.6
         else set5_matrix[i,j]=-0.6
      }
   }
}

cat("Simulating data......")

set1_data<-mvrnorm(n = 20, rep(0,500), Sigma=set1_matrix)
set2_data<-mvrnorm(n = 20, c(rep(0.75,50),rep(0,450)), Sigma=set2_matrix)
set3_data<-mvrnorm(n = 20, c(rep(0.75,50),rep(0,450)), Sigma=set3_matrix)
set4_data<-mvrnorm(n = 20, rep(0,500),Sigma=set4_matrix)
set5_data<-mvrnorm(n = 20, c(rep(0.75,25),rep(-0.75,25),rep(0,450)),Sigma=set5_matrix)
set6_data<-cbind(set2_data[,1:10],set3_data[,1:10],set4_data[,1:10],set5_data[,1:10],set1_data[,1:460])
control_data<-mvrnorm(n = 20, rep(0,500), Sigma=set1_matrix)

cat("Done!")

return(list(set1_data,set2_data,set3_data,set4_data,set5_data,set6_data,control_data))

}

