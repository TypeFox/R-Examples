plot.MNS <-
function(x, view="pop", subID=NULL, ...){
  # Plotting function for MNS object
  # Can be used to plot either estimated networks given data or to directly plot
  # simulated networks
  #
  #
  # INPUT:
  #      - x: MNS object (obtained by running MNS(...))
  #      - view: plotting view. Either "pop": plot population network, "var": plot variance network or "sub": plot subject networks
  #      - subID: if view="sub", select which subjects to plot
  #      -
  #
  # NOTE THAT THROUGHOUT THE AND RULE IS USED!!
  
  if (names(x)[1] == "PresPop"){
    # plotting estimated networks!!
    
    if (view=="pop"){
      # plot POPULATION network:
      # define igraph objects:
      
      Adj = x$PresPop
      Adj = Adj * ((Adj!=0) * t(Adj!=0))
      Adj = Adj + t(Adj)
      
      G = graph.adjacency(Adj, weighted = TRUE, diag = FALSE, mode="undirected")
      E(G)$width = abs(E(G)$weight) * 10
      E(G)$color = sapply(E(G)$weight, FUN=function(x){ifelse(x>0, "blue", "red")})
      
      plot(G, layout=layout.circle)
    } else if (view=="var"){
      # plot VARIANCE network:
      
      Adj = x$PresRE
      Adj = Adj * ((Adj!=0) * t(Adj!=0))
      Adj = Adj + t(Adj)
      
      G = graph.adjacency(Adj, weighted = TRUE, diag = FALSE, mode="undirected")
      E(G)$width = abs(E(G)$weight) * 10
      E(G)$color = sapply(E(G)$weight, FUN=function(x){ifelse(x>0, "blue", "red")})
      
      plot(G, layout=layout.circle)
    } else if (view=="sub"){
      # plot SUBJECT networks
      
      if (is.null(subID)){
        minLen = min(3, dim(x$PresBLUP)[3])
        subID = seq(1, minLen)
      }
      
      # define igraph objects:
      G = vector("list", length(subID))
      for (i in 1:length(subID)){
        Adj_i = x$PresPop + x$PresBLUP[,,subID[i]]
        Adj_i = Adj_i * ((Adj_i!=0) * t(Adj_i!=0))
        Adj_i = Adj_i + t(Adj_i)
        G[[i]] = graph.adjacency(Adj_i, weighted = TRUE, diag = FALSE, mode="undirected")
        E(G[[i]])$width = abs(E(G[[i]])$weight) * 5
        E(G[[i]])$color = sapply(E(G[[i]])$weight, FUN=function(x){ifelse(x>0, "blue", "red")})
        
        ii = (apply(get.edgelist(G[[i]]),1, FUN=function(ee){
          ifelse(x$PresRE[ee[1], ee[2]]!=0, return(TRUE), return(FALSE))
        }))
        E(G[[i]])$lty = ifelse(ii, 2,1)#"red", "blue")
      }
      
      par(mfrow=c(1, length(subID)))
      for (i in 1:length(G)) plot(G[[i]], layout=layout.circle)
    } else {
      stop("Unrecognized plotting view \n  Must be one of pop, var or sub.")
    }
  } else {
    # plotting SIMULATED networks!
    
    if (length(x)<=2){
      # simulated networks using Danaher method! Can only plot subject networks!!
      if (view!="sub") warning('plotting networks simulated using Danaher method - view parameter changed to view="sub"\nSee documentation for further details')    
      view = "sub"
    }
    
    if (view=="pop"){
      # plot POPULATION network:
      # define igraph objects:
      
      Adj = x$PopNet
      
      G = graph.adjacency(Adj, weighted = TRUE, diag = FALSE, mode="undirected")
      E(G)$width = E(G)$weight * 10
      E(G)$color = sapply(E(G)$weight, FUN=function(x){ifelse(x>0, "blue", "red")})
      
      plot(G, layout=layout.circle)
    } else if (view=="var"){
      # plot VARIANCE network:
      
      Adj = x$RanNet
      
      G = graph.adjacency(Adj, weighted = TRUE, diag = FALSE, mode="undirected")
      E(G)$width = abs(E(G)$weight) * 2
      E(G)$color = sapply(E(G)$weight, FUN=function(x){ifelse(x>0, "blue", "red")})
      
      plot(G, layout=layout.circle)
    } else if (view=="sub"){
      # plot SUBJECT networks
      
      if (is.null(subID)){
        minLen = min(3, length(x$Networks))
        subID = seq(1, minLen)
      }
      
      # define igraph objects:
      G = vector("list", length(subID))
      for (i in 1:length(subID)){
        Adj_i = x$Networks[[i]]
        G[[i]] = graph.adjacency(Adj_i, weighted = TRUE, diag = FALSE, mode="undirected")
        E(G[[i]])$width = abs(E(G[[i]])$weight) * 5
        E(G[[i]])$color = sapply(E(G[[i]])$weight, FUN=function(x){ifelse(x>0, "blue", "red")})
        
        ii = (apply(get.edgelist(G[[i]]),1, FUN=function(ee){
          ifelse(x$RanNet[ee[1], ee[2]]!=0, return(TRUE), return(FALSE))
        }))
        E(G[[i]])$lty = ifelse(ii, 2,1)#"red", "blue")
      }
      
      par(mfrow=c(1, length(subID)))
      for (i in 1:length(G)) plot(G[[i]], layout=layout.circle)
    } else {
      stop("Unrecognized plotting view \n  Must be one of pop, var or sub.")
    }
  }
}