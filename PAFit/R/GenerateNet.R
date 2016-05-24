# function to generate simulated network  2015-3-11 Thong Pham
GenerateNet<-
function(N=1000,m = 1, mode = c(1,2,3), alpha = 1, beta = 1, sat_at = 100,
         offset = 1, rate = 0, shape = 0,num_seed = 2, prob_m = FALSE,increase = FALSE,log = FALSE){
   if (num_seed >= N)
       stop("num_seed too large")
   if (num_seed < 2)
      stop("num_seed too small")
   if ((alpha < 0) || (beta < 0) || (sat_at < 0) || (rate < 0) || (shape < 0) || (m <= 0))
       stop("The parameters must be non-negative")  
   if ((mode[1] != 1) && (mode[1] != 3) && (mode[1] != 2))
       stop("Mode must be 1, 2 or 3 ")  
    graph <- vector("list", N)
    # Node weights in the BA model. 
    if (shape*rate > 0)
        fitness <- rgamma(N,rate = rate,shape = shape)
    else fitness <- rep(1,N)  
  
    for (n in 2:num_seed)
         graph[[n]] <- n - 1
    
    degree        <- rep(0,num_seed)
    names(degree) <- 1:num_seed 
    for (i in 1:num_seed) {
        count_degree <-table(graph[[i]][graph[[i]] <= i])
        degree[labels(count_degree)[[1]]] <- degree[labels(count_degree)[[1]]] + 
                                             count_degree
    }

    P         <- degree
    P[P == 0] <- offset
    seed.graph.size <- length(P)
    count <- 0
    for(n in seed.graph.size:(N-1))    # n is the size of the sub-network.
    {	  if (mode[1] == 1) {
            P.sum   <- sum(P^alpha*fitness[1:n])
            node.weights <- P^alpha*fitness[1:n]/P.sum
        } else if (mode[1] == 2) {
            temp    <- pmin(P,sat_at)^alpha
            P.sum   <- sum(temp*fitness[1:n])
            node.weights <- temp*fitness[1:n]/P.sum
        } else {
            temp    <- alpha*(log(P))^beta + 1 
            P.sum   <- sum(temp*fitness[1:n])
            node.weights <- temp*fitness[1:n]/P.sum
        }

        #Allow duplication:  
        if (TRUE == increase) {
            count <- count + 1  
            nodes <- sort(sample(1:n,size = ifelse(log,max(round(log(count)),1),count),prob = node.weights, replace = TRUE))      
        }else {
        if (prob_m == TRUE) {
            num_edge_temp <- rpois(1,lambda = m)
            if (num_edge_temp > 0)
                nodes <- sort(sample(1:n,num_edge_temp,prob = node.weights, replace = TRUE))
            else nodes <- NULL
        }
        else
            nodes <- sort(sample(1:n,size = m,prob = node.weights, replace = TRUE)) 
        }
        ##########################################
        degree  <- c(degree,0)
        if (0 != length(nodes)) {
            temp    <- table(nodes)
            graph[[n+1]]  <- c(graph[[n+1]], nodes)
            for(i in 1:length(temp)) { 
                num_edge  <- as.numeric(temp[i]) 
                node_name <- as.numeric(labels(temp[i]))
                degree[node_name]   <- degree[node_name] + num_edge # Update degrees.
            }
        }
        P <- degree
        P[degree == 0] <- offset     
    }
    num_of_edge <- sum(unlist(lapply(graph,function(x) length(as.vector(x)))))
    edge_list   <- matrix(nrow = num_of_edge,ncol = 3,0)
    sum_m       <- 0
    for (i in 1:N) {
        m_t <- length(graph[[i]][graph[[i]] <= i])
        if (m_t > 0) { 
            temp  <- as.vector(graph[[i]])
            edge_list[(sum_m + 1):(sum_m + m_t),3]   <-  ifelse(i > num_seed,i-num_seed,0)
            edge_list[(sum_m + 1):(sum_m + m_t),1]   <-  i
            edge_list[(sum_m + 1):(sum_m  +  m_t),2] <- temp 
            sum_m <- sum_m + m_t
        }
    }
    return(list(graph = edge_list, fitness = fitness))
}
