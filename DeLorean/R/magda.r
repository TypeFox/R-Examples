#
# This code copyright Magdalena Strauss: magdalena.strauss@mrc-bsu.cam.ac.uk
#
#Use function CombfuncPaths to 
#create path with a transition matrix of the form
#exp(-n*d), starting from around line 110, 
#functions before are only the components
#of this function

# library(seriation)
# library(fields) #"rdist"
#Constructing stochastic paths

stochasticPath1 <- function(P11_stoch, starting_point, number_paths)
{
  #generates stochastic paths with a given transition matrix P11_stoch
  #and one given starting point 

  t_max <- dim(P11_stoch)[1]
  visitedAll <- matrix(data=0,nrow=t_max,ncol=number_paths)

  apply_func <-  function(n)    
  {
    P1_stoch <- P11_stoch
    current_vertex <- starting_point
    P1_stoch[,current_vertex]=0
    P1_stoch[current_vertex,] <- P1_stoch[current_vertex,]/sum(P1_stoch[current_vertex,])
    visited <- matrix(data=0,nrow=t_max,ncol=1)
    visited[1] <- current_vertex
    not_visited <- (1:t_max)[-current_vertex]
    alpha<-runif(t_max,0,1)
    kk=2
    for (jj in not_visited)
    {
      P1_stoch[current_vertex,] <- P1_stoch[current_vertex,]/sum(P1_stoch[current_vertex,])
      if (alpha[jj]<P1_stoch[current_vertex,1])
      {
        index<-1
      }
      else
      {
        cumSum <- cumsum(P1_stoch[current_vertex,])
        if (alpha[jj]>cumSum[t_max])
        {
          index <- t_max
        }
        else
        {
          index <- which(cumSum> alpha[jj])[1]
        }
      }
      visited[kk] <- index
      kk=kk+1
      current_vertex <- index
      P1_stoch[,current_vertex]=0
      if (kk<t_max)
      {
        not_visited <- not_visited[not_visited!=current_vertex]}

    }
    return(visited)}
  output <- sapply(1:number_paths,apply_func)
  return(output)

}

stochasticPaths <- function(P11_stoch,starting_points,number_paths)
{
  #for several starting points
  #P11_stoch is the transition matrix
  #starting_points is a vector
  #number_paths is the number of paths for each starting point
  aaa<-sapply(starting_points,stochasticPath1,P11_stoch=P11_stoch, number_paths=number_paths,simplify="array")
  output <- matrix(data=aaa,nrow=dim(aaa)[1],ncol=dim(aaa)[2]*dim(aaa)[3])
  output<- output[,!(duplicated(t(output)))]
  return(output)
}
#######################
#Creating the matrix

rowStd <- function(M)
{
  #row-standardises a matrix 
  std<-function(v)
    {return (v/sum(v))}
  output<-t(apply(M,1,std))
  return(output)
}

expo <- function(M,n=1)
{
  row_Means<-rowMeans(M)
  output<-M
  for (i in 1:dim(M)[1])
  {
    output[i,]<-exp(-M[i,]*n+row_Means[i])
  }
  return(output)
}

P1 <- function(data,n)
{
  #data: rows=genes, cols=cells
  distMat <- fields::rdist(t(data),t(data))
  distMat <- expo(distMat,n)
  #computes exp(-n*distance)
  diag(distMat) <- 0
  distMat <- rowStd(distMat)
  return(distMat)
}


#######################################################################
#####Combined function
CombfuncPaths <- function(data,starting_points, number_paths=1000,n=10)
{
  #data: rows=genes, cols=cells
  #starting_points: vector
  #number_paths: number of paths per starting point
  #n as in exp(-n*distance)
  P11 <- P1(data,n)
  return(stochasticPaths(P11,starting_points,number_paths))
}


###Simply using TSP
TSP_order <- function(data)
{
  a<-seriation::seriate(as.dist(fields::rdist(t(data),t(data))),method="TSP")
  return(seriation::get_order(a))
}
