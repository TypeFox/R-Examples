#####################################################################################
#################################FCMapper###########################################
#####################################################################################
#Draft 1.1
#Feb 10, 2016
#Shaun Turney


#####################################################################################
#The first function, check.matrix,checks to see if all the values are between 
#-1 and 1 as well as making sure that the matrix is square. It will return warnings 
#if either of these checks are violated. The function checks whether there are self-
#loops and gives a warning if there are.
#####################################################################################

check.matrix = function (matrix) {
  
  ifelse (length(matrix[1,]) != length(matrix[,1]), 
          print("Warning: Matrix is not square. Matrix must be square for fuzzy cognitive mapping.")
          ,print("Matrix is square."))
  
  if (length(which(matrix>1)) > 0) { 
    problem_entries = paste(which(matrix>1), collapse = ", ")
    print(paste("Warning: The following entries are greater than 1:",problem_entries))
  }
  
  if (length(which(matrix< -1)) > 0) { 
    problem_entries = paste(which(matrix< -1), collapse = ", ")
    print(paste("Warning: The following entries are less than 1:",problem_entries))
  }
  
  if (length(which(matrix>1)) == 0 & length(which(matrix< -1)) == 0) { 
    print ("All values of the matrix are within -1 and 1.")
  }
  a=0
  for (b in 1:length(matrix[1,])) {
    if (matrix[b,b] != 0 & a==0) { 
      print("The diagonal is not equal to 0 (ie, there is a self-loop). Consider whether this is appropriate.")
      a=1
    }
  } 
}

#####################################################################################
#The second function, matrix.indices, gives the values of matrix-level indices.
#####################################################################################

matrix.indices = function(matrix) {
  
  Connections  = length(which(matrix!=0))
  Density = Connections/length(matrix)
  Concepts = length(matrix[1,])
  Transmitters = length(which(colSums(abs(matrix))==0 & rowSums(abs(matrix))!=0))  
  Receivers  = length(which(rowSums(abs(matrix))==0 & colSums(abs(matrix))!=0))
  NoConnection  = length(which(rowSums(abs(matrix))==0 & colSums(abs(matrix))==0))
  Ordinary  = Concepts - Transmitters - Receivers - NoConnection
  SelfLoops	= 0
  for (i in 1:Concepts) {
    if (matrix[i,i] != 0) { SelfLoops = SelfLoops + 1 }
  }
  Connectionspervariable = Connections/length(matrix[1,])
  Complexity = Receivers/Transmitters
  Outdegree = rowSums(abs(matrix))
  Hierarchy = (stats::var(Outdegree)) * 12/(length(matrix[1,])*(length(matrix[1,])+1)*(length(matrix[1,])+1))
  
  Value = c(Connections,Density,Concepts,Transmitters,Receivers,NoConnection,Ordinary,SelfLoops,
            Connectionspervariable,Complexity,Hierarchy)
  Index = c("Number of connections","Connection density","Number of concepts","Number of transmitters",
            "Number of receivers","Number of no connections","Number of ordinary","Number of self loops",
            "Connections/variable","Complexity (R/T)","Hierarchy")
  return (data.frame(Index,Value))
}

#####################################################################################
#The third function, concept.indices, gives the values of concept-level indices. The
#function must be fed both the matrix and the concept names.
#####################################################################################

concept.indices = function(matrix,concept.names) {
  
  Outdegree = rowSums(abs(matrix))
  Indegree = colSums(abs(matrix))
  Centrality = Outdegree + Indegree
  Concept = length(matrix[1,])
  Transmitter = numeric(length=Concept)
  Transmitter[which(colSums(abs(matrix))==0 & rowSums(abs(matrix))!=0)] = 1
  Receiver = numeric(length=Concept)
  Receiver[which(rowSums(abs(matrix))==0 & colSums(abs(matrix))!=0)] = 1
  Ordinary = numeric(length=Concept)
  Ordinary[which(rowSums(abs(matrix))!=0 & colSums(abs(matrix))!=0)] = 1
  NoConnection = numeric(length=Concept)
  NoConnection[which(rowSums(abs(matrix))==0 & colSums(abs(matrix))==0)] = 1
  
  Index = c("Concept","Outdegree","Indegree","Centrality","Transmitter","Receiver","Ordinary","No connection")
  Values = data.frame(concept.names,Outdegree,Indegree,Centrality,Transmitter,Receiver,Ordinary,NoConnection)
  colnames(Values) = Index
  return (Values)
  
}

#####################################################################################
#This function, nochanges.scenario, will iterate the matrix to find the equilibrium
#values of the concepts. It activates the matrix with a vector of 1s and squeezes the
#resulting vector with a logic function. The number of iterations must be fed to the
#function and can be increased if convergence is not reached. The function checks 
#for convergence and gives a warning if convergence isn't reached.
#####################################################################################

#no changes scenario

nochanges.scenario = function (matrix,concept.names,iter) {
  
  act_vector = matrix(0,nrow=iter,ncol=length(matrix[1,])) #Create activation vector of 1s
  act_vector[1,] = rep(1,length(matrix[1,]))
  
  
  for (i in 2:iter) { #Apply function iteratively
    act_vector[i,] = 1/(1+exp(-1*(act_vector[i-1,] + (act_vector[i-1,] %*% matrix))))
  }
  
  #check that convergence is reached
  
  if (all.equal (act_vector[iter,],act_vector[iter-1,]) != TRUE) {
    print("WARNING: Convergence not reached. Try increasing the number of iterations.")
  }
  
  results = data.frame(concept.names,act_vector[iter,])
  colnames(results) = c("Concept","Equilibrium_value")
  
  
  #plot change in concept values over time
  graphics::plot(act_vector[,1]~seq(1,iter,1),type="n",ylim=c(0,1),xlab="Iteration",ylab="Value")
  for (n in 1:length(matrix[1,])) {
    graphics::points(act_vector[,n]~seq(1,iter,1),type="l",col=n)
  }
  graphics::legend("topright",legend=concept.names,col=seq(1,n,1),lty=1)
  
  return (results)
}

#####################################################################################
#This similar function, changes.scenario, finds the equilibrium values of the concepts
#while fixing one or more concepts to a user-defined value between 0 and 1. In the 
#example below, "B" is set to 0.5. 
#####################################################################################

changes.scenario = function (matrix,concept.names,iter,set.concepts,set.values) {
  
  act_vector = matrix(0,nrow=iter,ncol=length(matrix[1,]))
  act_vector[1,] = rep(1,length(matrix[1,]))
  act_vector[1,which(concept.names %in% set.concepts == TRUE)] = set.values
  
  
  for (i in 2:iter) {
    act_vector[i,] = 1/(1+exp(-1*(act_vector[i-1,] + (act_vector[i-1,] %*% matrix))))
    act_vector[i,which(concept.names %in% set.concepts == TRUE)] = set.values
  }
  
  
  #check that convergence is reached
  
  if (all.equal (act_vector[iter,],act_vector[iter-1,]) != TRUE) {
    print("WARNING: Convergence not reached. Try increasing the number of iterations.")
  }
  
  results = data.frame(concept.names,act_vector[iter,])
  colnames(results) = c("Concept","Equilibrium_value")
  
  
  #plot
  graphics::plot(act_vector[,1]~seq(1,iter,1),type="n",ylim=c(0,1),xlab="Iteration",ylab="Value")
  for (n in 1:length(matrix[1,])) {
    graphics::points(act_vector[,n]~seq(1,iter,1),type="l",col=n)
  }
  graphics::legend("topright",legend=concept.names,col=seq(1,n,1),lty=1)
  
  return (results)
}

#####################################################################################
#The function comp.scenarios compares the equilibrium concept values of two scenarios.
#The function is fed the output of nochanges.scenario or changes.scenario. 
#####################################################################################

comp.scenarios = function (scenario1,scenario2) {
  
  #check that concepts are the same 
  if (paste(scenario1$Concept,collapse="")!=paste(scenario2$Concept,collapse="")) {
    print ("WARNING: The two scenarios do not share the same concepts and so are non-comparable.")
  }
  
  #output
  
  diff = scenario2[,2]-scenario1[,2]
  percent.change = (diff/scenario1[,2]) * 100  
  
  
  results = data.frame(scenario1,scenario2[,2],diff,percent.change)
  colnames(results) = c("Concept","Scenario_1","Scenario_2","Difference","Percent_change")
  return(results)
}

#####################################################################################
#The function comp.maps compares two concept maps by determining their similarity
#coefficients: S2 and Jaccard.
#####################################################################################

comp.maps = function (concept.names1,concept.names2) {
  allnames = c(concept.names1,concept.names2)
  S2= (length(allnames) - length(unique(allnames)))/length(allnames)
  
  a=0
  b=0
  c=0
  
  for (x in 1:length(allnames)) {
    if (length(which(concept.names1==unique(allnames)[x])) != 0 & length(which(concept.names2==unique(allnames)[x])) == 0) {
      c = c+1
    }
    if (length(which(concept.names1==unique(allnames)[x])) == 0 & length(which(concept.names2==unique(allnames)[x])) != 0) {
      b = b+1
    }
    if (length(which(concept.names1==unique(allnames)[x])) != 0 & length(which(concept.names2==unique(allnames)[x])) != 0) {
      a = a+1
    }
  }
  Jaccard = a/(a+b+c)
  
  results=data.frame(S2,Jaccard)
  colnames(results) = c("S2","Jaccard")
  return(results)
}


#####################################################################################
#The cognitive map is fed to the package igraph in order to plot the map.
#The size of the concepts and the edges are determined by their weights. The size of
#the concepts is fed to the function as their equilibrium values. Negative edges are
#red and positive edges are black. The labels are shown for the concepts.
#####################################################################################

graph.fcm = function (matrix,concept.sizes,concept.names) {
  
  matrix.plot = igraph::graph.adjacency(matrix,mode="directed",weighted=T) #put into format igraph can read
  
  igraph::V(matrix.plot)$size = concept.sizes * 40
  
  
  igraph::E(matrix.plot)$color = ifelse(igraph::E(matrix.plot)$weight<0,"red","black")
  edge.labels = ifelse(igraph::E(matrix.plot)$weight<0,"-","+")
  
  edge.curved = rep(0,length(igraph::E(matrix.plot))) #should the arrows be curved (yes, if arrows are in both directions)
  i=1
  for (x in 1:length(matrix[1,])) {
    for (y in 1:length(matrix[1,])) {
      if(matrix[x,y]!=0 & matrix[y,x]!=0) {
        edge.curved[i] = 0.5
      }  
      if(matrix[x,y]!=0) {
        i = i+1
      }
    }
  }
  
  igraph::tkplot(matrix.plot,edge.width=abs(igraph::E(matrix.plot)$weight*5), vertex.color="grey",
         vertex.label.color="blue",vertex.label=concept.names,vertex.label.dist=0,
         edge.curved=edge.curved)
  
}

#####################################################################################
#Multiple cognitive maps are combined together by the function combine.maps. The
#function is fed two matrices and the concept names of the two matrices. For edge
#values shared between the two maps the average value is taken. More than two maps
#can be aggregated by applying the function more than once.
#####################################################################################

combine.maps = function(matrix1,matrix2,concept.names1,concept.names2) {
  
  #What are the unique concepts?
  concept.names.agg = c(concept.names1,concept.names2)
  concept.names.agg = unique(concept.names.agg)
  
  #Matrix of correct size
  matrix.agg1 = matrix(0,ncol=length(concept.names.agg),nrow=length(concept.names.agg))
  
  #name rows and columns of all matrices
  colnames(matrix1) = concept.names1
  rownames(matrix1) = concept.names1
  colnames(matrix2) = concept.names2
  rownames(matrix2) = concept.names2
  colnames(matrix.agg1) = concept.names.agg
  rownames(matrix.agg1) = concept.names.agg
  
  matrix.agg2 = matrix.agg1
  
  #Apply matrix 1 and matrix 2 values to aggregated matrix
  for (x in 1:length(concept.names1)) {
    for (y in 1:length(concept.names1)) {
      matrix.agg1[concept.names1[x],concept.names1[y]] = matrix1[concept.names1[x],concept.names1[y]]
    }
  }
  
  for (x in 1:length(concept.names2)) {
    for (y in 1:length(concept.names2)) {
      matrix.agg2[concept.names2[x],concept.names2[y]] = matrix2[concept.names2[x],concept.names2[y]]
    }
  }
  
  #add together matrices
  matrix.agg = matrix.agg1 + matrix.agg2
  
  #take average of shared values
  for (n in 1:length(concept.names.agg)) {
    for (m in 1:length(concept.names.agg)) {
      if (matrix.agg1[n,m] != 0 & matrix.agg2[n,m] != 0) {
        matrix.agg[n,m] = mean(c(matrix.agg1[n,m],matrix.agg2[n,m]))
      }
    }
  }
  
  return(matrix.agg)
}


