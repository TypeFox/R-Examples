factors.to.graph <- function(factors, sep = '(')
{
  # compute adjacency matrix from factor chain
  # eg: (1)(2)(3|1,2)
  # accepts '(' or '[' as separator character
  # nodes are numbers from 1 to N
  # DOES NOT CHECK FOR CORRECTNESS OF INPUT
  l <- list()
  if (sep == '(') {
    sep1 = '('
    sep2 = ')'
  } else if (sep == '[') {
    sep1 = '['
    sep2 = ']'    
  } else {
    # error
  }
  
  for (i in unlist(strsplit(factors, sep1, TRUE)))
    for (j in unlist(strsplit(i, sep2, TRUE)))
      l[[length(l)+1]] <- list(j)
  
  num_nodes = length(l)
  am = matrix(rep(0, num_nodes*num_nodes), c(num_nodes,num_nodes))
  for (i in l)
  {
    item <- unlist(strsplit(unlist(i), "|",TRUE))
    if (length(item) > 1)
    {
      to <- as.integer(item[1])
      for (j in unlist(strsplit(unlist(item[2]), ",", TRUE)))
      {
        from <- as.integer(j)
        am[from,to] <- 1
      }    
    }
  }
  return(am)
}

graph.to.factors <- function(am, sep = '(')
{
  # compute factor chain from adjacency matrix
  # accepts '(' or '[' as separator character
  # nodes are numbers from 1 to N
  # DOES NOT CHECK FOR CORRECTNESS OF INPUT
  l <- list()
  if (sep == '(') {
    sep1 = '('
    sep2 = ')'
  } else if (sep == '[') {
    sep1 = '['
    sep2 = ']'    
  } else {
    # error
  }
  
  factors <- c()
  
  # build up string node after node
  for (i in 1:nrow(am))
  {
    factor = c(sep1,i)
    parents <- which(am[,i] > 0)
    if (length(parents) > 0)
    {
      factor <- c(factor,'|')
      # build up parents
      while(length(parents) > 1)
      {
        factor <- c(factor, parents[1],',')
        parents <- parents[-1]
      }
      factor <- c(factor, parents[1])
    }
    factor <- c(factor,sep2)
    factors <- c(factors, factor)
  }
  factors <- paste(unlist(factors), collapse='')
  return(factors)
}

# create a file with scores in the format used by Jaakkola et al, Cussens.
# 1st row: number of nodes
# then, for each node:
# node, number of cpcs
# followed by score, length of cpc, cpc as list of nodes in ascending order, one per row
#
# start.cpcs.vars should contain also the start index of edge vars
save.scores <- function(instance.name, num.nodes, scores, cpcs.vars, start.cpcs.vars)
{
  file.name <- paste(instance.name, ".scores", sep="")
  print(file.name)
  row <- NULL
  k <- 1
  row[[k]] <- paste(num.nodes)
  k <- k + 1
  for (i in 1:num.nodes)
  {
    row[[k]] <- (paste(i, start.cpcs.vars[i+1] - start.cpcs.vars[i]))
    k <- k + 1
    for (j in 0:(start.cpcs.vars[i+1]-start.cpcs.vars[i]-1))
    {
      cpc <- c(unlist(cpcs.vars[[start.cpcs.vars[i]+j]]))
      scpc <- paste(cpc, collapse=" ")
      row[[k]] <- paste(scores[start.cpcs.vars[i]+j], length(cpc), scpc)
      k <- k + 1 #print(s)
    }
  }
  write(row, file=file.name)
}

# read a file of scores in the format used by Jaakkola et al, Cussens.
# return the informations the orientation.init method would have returned,
# plus the whole cplex problem filled
# TODO FINISH
read.scores <- function(file.name, cplex.problem)
{
#   l <- readLines(file.name)
#   num.nodes <- strtoi(l[1])
#   k <- 1
#   scores    <- c()
#   cpcs      <- NULL
#   cpcs.vars <- NULL
#   start.cpcs.vars <- c()
#   start.edge.vars <- -1
#   num.cpcs.vars   <- length(l) - 1 - num.nodes
#   
#   edge.vars <- matrix(rep(0, num.nodes*num.nodes), nrow = num.nodes)
#   
#   curr.line <- 2
#   for (i in 1:num.nodes)
#   {
#     start.cpcs.vars <- c(start.cpcs.vars, k)
#     line <- c(unlist(strsplit(l[curr.line], split = " ")))
#     num.i.vars <- strtoi(line[2])
#     curr.line <- curr.line + 1
#     local.cpcs.vars <- c()
#     for (j in 1:num.i.vars)
#     {
#       line <- c(unlist(strsplit(l[curr.line], split = " ")))
#       scores <- c(scores, as.numeric(line[1]))
#       len.cpc <- strtoi(line[2])
#       cpc.vs <- c()
#       if (len.cpc > 0)
#       {
#         for (m in 1:len.cpc)
#         {
#           cpc.vs <- c(cpc.vs, as.integer(line[2+m]))
#         }
#         local.cpcs.vars <- unique(c(local.cpcs.vars, cpc.vs))
#       }
#       cpcs.vars[[k]] <- as.list(cpc.vs)
#       curr.line <- curr.line + 1
#       k <- k + 1
#     }
#     for (u in local.cpcs.vars)
#       edge.vars[u, i] <- 1
#     cpcs[[i]] <- as.list(sort(local.cpcs.vars))
#   }
#   # print(edge.vars)
#   # print(cpcs)
#   # print(start.cpcs.vars)
#   # print(cpcs.vars)
#   
#   start.edge.vars <- k
#   
#   # create combinations
#   
#   # create constraints and fill into model
}