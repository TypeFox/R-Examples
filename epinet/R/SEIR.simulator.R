SEIR.simulator <- function(M, N, beta, ki, thetai, ke = ki, thetae = thetai, latencydist = "fixed", 
                           latencyperiod = 0)
{
  
  # First do some input format checking
  
  if (!is.null(M))
  {
    # M should be an edgelist matrix
    if(!is.matrix(M))
      stop("Input error: Network M must be an edgelist matrix.")
    
    # Make sure M has the correct dimensions
    if ( (length(dim(M)) != 2) || (dim(M)[2] != 2) )
      stop("Input error: Network M must an a 2-dimensional edgelist matrix.")
  }
  
  # convert to network data structure
  net = network::network.initialize(N, directed = FALSE)
  network::add.edges(net, head = M[,1], tail = M[,2])
  
  # t_net is a logical array indicating whether edge in net is between an SI pair  
  t_net = array(FALSE, network::network.edgecount(net))
  all_edge = unlist(net$mel)
  # v1 and v2 are the vertices associated with edges in net
  v1 = all_edge[seq(1,length(all_edge),by = 3)]
  v2 = all_edge[seq(2,length(all_edge),by = 3)]
  
  # lists of current states of state
  is_infectious = is_exposed = is_removed = array(FALSE, N)
  is_susceptible = array(TRUE, N)
  
  # get location of first infection
  init <- sample(1:N,1)  # Inital infected individual is chosen at random
  
  # update state vectors
  is_infectious[init] = TRUE
  is_susceptible[init] = FALSE
  
  # update live transmission graph
  # get all neighbours
  neighbour_id = network::get.edgeIDs(net,init)
  # find susceptible ones
  susc_neighbour_id = neighbour_id[is_susceptible[v1[neighbour_id]] | is_susceptible[v2[neighbour_id]]]
  # add them
  t_net[susc_neighbour_id] = TRUE
  
  # Keep a list of all upcoming transition and recovery times (t.time[i] = r.time = NA if i is susceptible)
  t.time = array(dim = N)
  r.time = array(dim = N)
  
  # Generate a transition time and recovery time for initial infected
  t.time[init] <- ifelse( latencydist=="fixed", latencyperiod, rgamma(1, ke, scale = thetae) )  
  r.time[init] = rgamma(1, ki, scale = thetai) + t.time[init]
  
  nextrec = init  	# Keep track of who, of current infecteds, is next to recover
  
  inf.list <- matrix(c(init, NA, 0, t.time[init], NA), nrow = 1)   # Keep track of initial infection
  
  time <- cm.time <- t.time[init]
  
  nexttrans = init # Temporary value
  t.time[init] <- Inf # Temporary value
  
  # keep track of number exposed/infectious
  count_i = 1
  count_e = 0
  
  # Maximum number of iterations is 3*N (1 infection ,1 transition from  exposed to infective, and 1 recovery for every node)
  for( i in 2:(N*3) )	{
    # count number of SI pairs
    n.si = sum(t_net)
    
    # Draw waiting times for the next removal, infection, and transition
    dwt <- ifelse( count_i > 0, r.time[nextrec] - cm.time, Inf ) 
    bwt <- ifelse( n.si != 0, rexp(1,beta*n.si), Inf)
    twt <- t.time[nexttrans] - cm.time
    
    ewt <- min(bwt, dwt, twt, na.rm = TRUE)
    time <- c(time, ewt)			# Increment time
    cm.time <- cm.time + ewt		# Increment cumulative time
    
    if (ewt == bwt) { 
      test <- "Infect"
    } else if (ewt == dwt) {
      test <- "removal" 
    } else {
      test <- "transition"
    }
    
    if (test == "Infect"){ 
      # Event is an infection	
      
      # choose the edge where transmission occured
      trans_edge_id = ifelse(n.si > 1, sample(which(t_net), 1), which(t_net))
      
      # see which end is new infection and which is the parent
      if ( is_susceptible[v1[trans_edge_id]] ){
        new.inf = v1[trans_edge_id]
        parent = v2[trans_edge_id]
      } else {
        new.inf = v2[trans_edge_id]
        parent = v1[trans_edge_id]
      }
      
      # Generate a transition time for new infected
      lat <- ifelse( latencydist == "fixed", latencyperiod, rgamma(1, ke, scale = thetae) )
      
      t.time[new.inf] <- cm.time + lat
      
      # Update state vectors
      is_exposed[new.inf] = TRUE
      is_susceptible[new.inf] = FALSE
      # update number of exposed 
      count_e = count_e + 1
      
      # update live transmission graph --- remove edges adjacent to newly exposed 
      t_net[network::get.edgeIDs(net, new.inf)] = FALSE
      
      inf.list <- rbind(inf.list,c(new.inf, parent, cm.time, NA, NA))
      
      nexttrans <- which(t.time == min(t.time, na.rm = TRUE))
      
    } else if (test == "removal")	# Event is a removal
    {								
      if(i==2)	# Only infected dies out
      {
        inf.list[1,5] <- cm.time
        break
      } 	
      
      new.rec <- nextrec
      
      # Update lists of infectious, removed
      is_infectious[new.rec] = FALSE
      is_removed[new.rec] = TRUE
      # update number of infectious
      count_i = count_i - 1
      
      # update live transmission graph --- remove edges adjacent to newly recovered 
      t_net[network::get.edgeIDs(net, new.rec)] = FALSE
      
      # record removal time
      inf.list[which(inf.list[,1]==new.rec), 5] <- cm.time
      
      # Update recovery times and next infected to recover
      r.time[nextrec] <- NA
      if(count_i > 0) 
        nextrec <- which( r.time == min(r.time, na.rm = TRUE) ) 
      else if (count_e > 0)
      {
        nextrec <- which( t.time == min(t.time, na.rm = TRUE) )
        r.time[nextrec] <- Inf
      }
      
    } else 		# Event is a transition
    { 	 
      new.trans <- nexttrans
      
      # Update lists of susceptible, exposed and infecteds	 
      is_exposed[new.trans] = FALSE
      is_infectious[new.trans] = TRUE
      
      # update count of infecteds/exposeds
      count_i = count_i + 1
      count_e = count_e - 1
      
      # update live transmission graph 
      # get all neighbours of newly infectious
      neighbour_id = network::get.edgeIDs(net, new.trans)
      # find susceptoble ones
      susc_neighbour_id = neighbour_id[is_susceptible[v1[neighbour_id]] | is_susceptible[v2[neighbour_id]]]
      # add them
      t_net[susc_neighbour_id] = TRUE
      
      # record time of transition
      inf.list[which(inf.list[,1]==new.trans),4] <- cm.time
      
      # Update transition times, recovery times, and assign a recovery time to the new transition
      t.time[nexttrans] <- NA      
      nexttrans <-  which(t.time == min(t.time,na.rm = TRUE))
      r.time[new.trans] <- cm.time + rgamma(1,ki, scale = thetai)
      if (r.time[new.trans] < r.time[nextrec]) nextrec <- new.trans
    }
    if (count_e + count_i == 0) {break}	# No more infectious or exposed members, so epidemic ends
  }		
  
  # reset time 0 to be time of first R
  inf.list[,3:5] <- inf.list[,3:5] - min(inf.list[,5])
  
  # add any never infecteds to list
  if (any(is_susceptible)) inf.list <- rbind(inf.list, cbind(which(is_susceptible), NA, NA, NA, NA) )
  
  # fix up column names
  colnames(inf.list) <- c("Node ID", "Parent", "Etime", "Itime", "Rtime")
  
  class(inf.list) <- c("epidemic", "matrix")
  
  return(inf.list)
}

# FUNCTION print.epidemic
# Print method for class epidemic

print.epidemic <- function(x, ...)
{
    class(x) <- "matrix"
    print(x)
    class(x) <- c("epidemic", "matrix")
}

# FUNCTION summary.epidemic
# Summary method for class epidemic

summary.epidemic <- function(object, ...)
{
    N = dim(object)[1]
    ninf <- min(which(is.na(object[,5])),dim(object)[1] + 1) - 1
    cat("Epidemic data: ", ninf, "infecteds,", N - ninf, "susceptibles, ", N, "total individuals. \n")
    cat("Length of epidemic: ", max(object[,5], na.rm = TRUE) - min(object[,3], na.rm = TRUE), "\n \n")
}


# FUNCTION plot.epidemic
# Plot method for class epidemic

plot.epidemic <- function(x, lwd = 1, leaf.labs = TRUE, leaf.cex = 0.75,
zero.at.start = FALSE, main = "Transmission Tree", xlab = "Time", ylab= "",
e.col = "black", i.col = "red", lty.transmission = 3, marktransitions = TRUE,
label.trans = "|", cex.trans = 0.5, ...){

    plotepitree(epi = x, lwd = lwd, leaf.labs = leaf.labs, leaf.cex = leaf.cex,
        zero.at.start = zero.at.start, main = main, xlab = xlab, ylab = ylab,
        e.col = e.col, i.col = i.col, lty.transmission = lty.transmission,
        marktransitions = marktransitions, label.trans = label.trans,
        cex.trans = cex.trans, ...)
}
