# Function to compute the survival signature of an arbitrary system
# Define a graph of the system using the igraph package, with additional nodes named ``s'' and ``t'' for either end of the system diagram (see example below)
# Additionally, create a vertex attribute named compType which names the different component types
# Pass this to the graph object to the function to compute the signature
computeSystemSurvivalSignature <- function(graph, cutsets = NULL, frac = FALSE) {
  # Check that we have s,t nodes as a minimum
  if(sum(is.na(match(c("s","t"), V(graph)$name)))) {
    stop("The graph object must contain vertices named 's' and 't' to indicate the start and terminal points of the system being specified.")
  }
  
  # Check component type attribute has been set
  if(is.na(match("compType", list.vertex.attributes(graph)))) {
    stop("For survival signatures, the type of each component must be stored in a vertex attribute named 'compType'.")
  }
  # Define the component type of s,t nodes as unknown
  V(graph)$compType[match(c("s","t"), V(graph)$name)] <- NA
  # And check that otherwise the list is complete
  if(length(V(graph)$name[-match(c("s","t"), V(graph)$name)]) != length(na.omit(V(graph)$compType[-match(c("s","t"), V(graph)$name)]))) {
    stop("Some component types are missing in 'compType' vertex attribute.")
  }
  
  # Build component type to vertex mapping list
  CbyType <- list()
  for(type in unique(na.omit(V(graph)$compType))) {
    CbyType[[type]] <- as.integer(V(graph)[which(V(graph)$compType==type)])
  }
  if(length(names(CbyType)) > 1)
    CbyType <- CbyType[order(names(CbyType))]
  
  # Handy stuff to know
  numC <- length(na.omit(unique(V(graph)$compType))) # number of types of component
  Cnums <- table(na.omit(V(graph)$compType)) # number of each type of component
  totalC <- sum(Cnums) # total number of components
  
  # Construct the survival signature table
  survSig <- matrix(0, nrow=prod(Cnums+1), ncol=numC+1)
  tmpT <- 1
  tmpE <- prod(Cnums+1)
  for(i in 1:numC) {
    tmpE <- tmpE/(Cnums[i]+1)
    survSig[,i] <- rep(0:Cnums[i], times=tmpT, each=tmpE)
    tmpT <- tmpT*(Cnums[i]+1)
  }
  #survSig <- as.data.frame(survSig)
  
  # Now start finding the survival signature
  if(is.null(cutsets)) {
    cutsets <- minimalVertexPairCutSets(graph, "s", "t")
  }
  
  # Check all states
  comp <- as.integer(V(graph))[-match(c("s","t"), V(graph)$name)]
  for(x in 0:{2^totalC-1}) {
    state <- digitsBase(x, base=2, ndigits=totalC)
    failed <- comp[state==0]
    working <- comp[state==1]
    if( sum(vapply(cutsets, function(x) { prod(x %in% failed) == 1 }, TRUE)) == 0 ) { # TRUE only if the system is working
      #numWorkingByType <- sapply(CbyType, function(x) { sum(working%in%x) })
      #print(which(apply(survSig[,1:numC], 1, identical, y=as.double(vapply(CbyType, function(x) { sum(working%in%x) }, 1)))))      
      survSig[which(apply(survSig[,1:numC,drop=FALSE], 1, identical, y=as.double(vapply(CbyType, function(x) { sum(working%in%x) }, 1)))),numC+1] <- survSig[which(apply(survSig[,1:numC,drop=FALSE], 1, identical, y=as.double(vapply(CbyType, function(x) { sum(working%in%x) }, 1)))),numC+1]+1
    }
  }
  
  # Normalise
  survSig <- as.data.frame(survSig)
  if(is.null(names(CbyType))) {
    names(survSig) <- c(1:numC, "Probability")
  } else {
    names(survSig) <- c(names(CbyType), "Probability")
  }
  if(!frac) {
    survSig[,numC+1] <- survSig[,numC+1]/apply(survSig[,-(numC+1),drop=FALSE], 1, function(x) { prod(choose(Cnums, x)) })
  } else {
    survSig[,numC+1] <- apply(survSig, 1, function(x) {
      cd <- gcd(x[numC+1], prod(choose(Cnums, x[-(numC+1)])))
      num <- x[numC+1]/cd
      denom <- prod(choose(Cnums, x[-(numC+1)]))/cd
      if(num == 0) { return("0") }
      if(num == denom) { return("1") }
      paste(x[numC+1]/cd, "/", prod(choose(Cnums, x[-(numC+1)]))/cd, sep="") })
  }
  survSig
}
