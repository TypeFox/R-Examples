#-------------------- ITERATIVE K MEANS --------------------
#Iterative K means personal implementation
#NOTE: Good k estimation but time consuming
# x-vector to cluster
# minimum- Minimum number of clusters (default 2)
# maximum- Maximum number of clusters (default 10)
# returns -Optimal clustering
iterativeKmeans=function(x,minimum=2,maximum=10,choice=0.5)
  {
  n=length(x)
  
  #Input checking
  if(minimum<2)
    {
    print("Error: clustering must divide data in at least 2 clusters")
    print("       (minimum must value at least 2)")
    break
    }
  if(maximum>=n)
    {
    print("Error: clustering must divide data in as much as n-1 clusters")
    print("       (maximum value at max is length(x)-1)")
    break    
    }
  if(maximum<minimum)
    {
    print("Error: maximum value must be greater than minimum")
    break
    }  
    
  numClusterings=maximum-minimum+1;
  k=minimum
  ss=c()
  while(k<=maximum)
    {
    ss=c(ss,mean(kmeans(x,k,iter.max=20)$withinss))
    k=k+1
    }
  diferencias=diff(ss)
  md=mean(diferencias)
  optimo=1
  for(i in 1:length(diferencias))
    {
    if(diferencias[i]>md) 
      {
      optimo=i
      break
      }
    }
  optimo=optimo+minimum-1
  kmeans(x,optimo,iter.max=20)
  }