# function to summarize statistics from a growing network  2015-3-11 Thong Pham
GetStatistics <-
function(data,net_type = c("directed","undirected"), only_PA = FALSE, Binning = FALSE,G = 1000, start_deg = 0,  
         deg_threshold = 1, CompressMode = c(0,1,2,3), CompressRatio = 0.5 , CustomTime = NULL){

    
    time_stamp        <- data[,3]
    in_node           <- data[,2]
    out_node          <- data[,1]
    node_id           <- sort(union(in_node,out_node))
   
    if (net_type[1] == "directed") {
        deg           <- table(in_node)
    } else
    if (net_type[1] == "undirected")
        deg           <- table(as.vector(as.matrix(data[,1:2])))        
    
    deg_new           <- rep(0,length(node_id))
    names(deg_new)    <- node_id
    deg_new[labels(deg)[[1]]] <- deg
    deg               <- deg_new
    final_deg         <- deg
    deg.max           <- as.numeric(max(deg))
    if (start_deg >= deg.max)
        stop("Starting degree too large!") else 
    if (start_deg < 0)
        stop("Negative starting degree!")
    unique_time       <- sort(unique(time_stamp))
    T                 <- length(unique_time)
    N                 <- length(node_id)
    ##############  Binning #########################
    #We have to cover from start_deg degree to deg.max degree, that is an interval with length deg.max - start_deg + 1
    if ((TRUE == Binning) && (G > 0) && (G <= deg.max - start_deg + 1)) {
        if (1 == G) {
            base            <- deg.max - start_deg + 1
            interval_length <- deg.max - start_deg + 1
        } else {
            #find the base that gives exactly G bin
            is.warn <- options()$warn 
            options(warn = -1) #temporily supress warning
            ff <- function(x){deg.max - start_deg + 1.0 - sum(floor(x^(0:(G - 1))))}
            base      <- uniroot(ff,interval = c(1 + 1e-15,deg.max - start_deg + G + 1.1),tol = .Machine$double.eps)$root
            options(warn = is.warn)
            interval_length <- floor(base^(0:(G-1)))
        }
    } else if ((FALSE == Binning) || (0 == G) || (G > deg.max - start_deg + 1)) {
        G               <- deg.max - start_deg + 1
        interval_length <- rep(1,G)
        base            <- 1
    }
    
    bin_vector   <-rep(G + start_deg - 1, deg.max + 1)  # degree 0 to deg.max
    # The right-end degree of the bins
    begin_deg   <- c(start_deg,start_deg + cumsum(interval_length)[-G])
    end_deg     <- begin_deg + interval_length - 1
    if (start_deg > 0)
        bin_vector[1:start_deg]  <- 0:(start_deg - 1)
    for (i in 1:G) 
        bin_vector[(begin_deg[i]:end_deg[i]) + 1]  <- i + start_deg - 1
    ########### Compress Time stamp #########################
    
    if (1 == CompressMode[1]) {
        T_compressed           <- round(CompressRatio*(T - 1))
        compressed_unique_time <- floor(seq(1,T - 1,length.out = T_compressed))
    } else if (2 == CompressMode[1]){
        edge_cumsum             <- cumsum(as.vector(table(time_stamp))) 
        edge_ratio              <- edge_cumsum/edge_cumsum[T]
        compressed_unique_time  <- unique_time[which(edge_ratio >= 1 - CompressRatio)]
        temp                    <- edge_ratio[which(edge_ratio < 1 - CompressRatio)]
        CompressRatio           <- temp[length(temp)]
        T_compressed            <- length(compressed_unique_time)
    }  else if (3 == CompressMode[1]) {
        compressed_unique_time <- sort(unique(CustomTime))
        CompressRatio          <- length(compressed_unique_time)/(T - 1)
    } else {
        #No Time compression
        CompressRatio   <- 1
        T_compressed    <- T
        compressed_unique_time <- unique_time
    }
    if (net_type[1] == "directed")
        first_edge       <- table(in_node[time_stamp == unique_time[1]]) else 
    if (net_type[1] == "undirected")
        first_edge       <- table(data[time_stamp == unique_time[1],1:2])
       
    first_deg        <- rep(0,N)
    increase         <- rep(0,N)
    names(increase)  <- node_id
    names(first_deg) <- node_id
    first_deg[labels(first_edge)[[1]]] <- first_edge
    # this is not the true number of new edges, due to the first appearance of a node together with some edges
    # but this can be used as an crude first step selection
    # the accurate selection can be done after
    inc              <- deg - first_deg  
    increase         <- inc
    initial_nodes    <- length(first_edge)
    if (FALSE == only_PA) {
        pos_temp          <- inc >= deg_threshold
        if (sum(pos_temp) == 0)
            stop("Degree threshold is too high. Please decrease degree threshold.")  
        f_position        <- node_id[pos_temp]
    } else f_position        <- NULL
    
    node_id_old <- node_id
    if (FALSE == only_PA)
        node_id     <- node_id[inc >= deg_threshold]
    
    N_new            <- length(node_id) 
    degree_appear    <- rep(0,deg.max + 1)
    Sum_m_k          <- rep(0,start_deg + G)
    n_tk             <- matrix(0,nrow = T_compressed - 1, ncol = start_deg + G)
    m_tk             <- matrix(0,nrow = T_compressed - 1, ncol = start_deg + G)
    m_t              <- rep(0,T_compressed - 1)
  
    
    if (FALSE == only_PA) {
        offset_tk               <- matrix(0,nrow = T_compressed - 1, ncol = start_deg + G) 
        z_j                     <- rep(0,N_new)
    
        node_degree <- matrix(-1,nrow = T_compressed - 1, ncol = N_new) 
       
    } else {
        offset_tk               <- matrix(0,0,0)
        z_j                     <- vector()  
        node_degree             <- matrix(0,0,0)
    }
   
    undirected  = 0;
    max_node_id = max(node_id_old);
    only_PA_num = ifelse(only_PA,1,0);
    .get_stats(time_stamp,unique_time,in_node,out_node,node_id_old,node_id,bin_vector, max_node_id, undirected, 
              only_PA_num,              
              compressed_unique_time,
              Sum_m_k,n_tk,m_tk,m_t,offset_tk,z_j,node_degree)

    if (FALSE == only_PA) {
      names(z_j)            <- node_id
      colnames(node_degree) <- node_id
    }
   
    names(node_id)    <- node_id
    #now perform the final selection
    true                           <- which(z_j >= deg_threshold)
    if (FALSE == only_PA) {
        if (length(true) == 0)
            stop("Degree threshold is too high. Please decrease degree threshold.")  
        increase[inc >= deg_threshold] <- z_j
        z_j                            <- z_j[true]
        f_position                     <- f_position[true]
    }
    node_degree                    <- node_degree[,true,drop=FALSE]
    

    result  <- list(offset_tk = offset_tk,net_type = net_type[1], n_tk = n_tk,m_tk = m_tk, bin_vector = bin_vector, Sum_m_k = Sum_m_k,
                    node_degree = node_degree,m_t = m_t,z_j = z_j, initial_nodes = initial_nodes,
                deg_thresh = deg_threshold, final_deg = final_deg, only_PA = only_PA, 
                increase = increase, start_deg = start_deg, Binning = Binning, G = G, 
                CompressMode = CompressMode[1], f_position = f_position, compressed_unique_time = compressed_unique_time, begin_deg = begin_deg, end_deg = end_deg,
                interval_length = interval_length,node_id = node_id_old, N = N, T = T, T_compressed = T_compressed,deg.max = deg.max, CompressRatio = CompressRatio , CustomTime = CustomTime)
    class(result) <- "PAFitData"
   
    return(result)
}

.onUnload <- function (libpath) {
  library.dynam.unload("PAFit", libpath)
}
