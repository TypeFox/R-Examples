dag.sim  <- function (dag, b = rep(0, nrow(dag$arc)), bxy = 0, n, 
                          mu = rep(0, length(dag$x)),
                      binary = rep(0, length(dag$x)),
                       stdev = rep(0, length(dag$x)), naming = 2, seed = NA, verbose = FALSE) 
{

  # arguments:
  #
  # b  - the regression coefficient (on linear scale)
  # mu - for continuous variables, the overall mean simulated
  #      for binary variables, proportion of successes simulated
  # stdev - for continuous variables without ancestors, the sd of the simulated distribution
  #         for other variables, the sd of the noise added to the calculated value
  #         (might not be very reasonable for binary variables, where it just biases towards 0,
  #          i.e. reducing the b should do the same; for continuous, it broadens the conf.int.)
  #         Note. the noise is added after shifting to the mean, so the mean will not be exact.
  #         Note. the noise is added before calculating descendant nodes, i.e. it is sort of
  #               true inter-individual variation, not measurement error.
  # binary - 1 for binary nodes, 0 for continuous nodes
  #         Note. the nodes are made binary before calculating descendant nodes, i.e. there
  #               will be a discrete difference--not a continuous change--in the descendants.


  logit<-function(p) log(p/(1-p));            # helper function;
  inv.logit<-function(l) exp(l)/(1+exp(l));   # helper function;

  if (length(b) != nrow(dag$arc)) {
      stop("number of coefficients does not fit number of arcs");
  }
      
  writeLines("Note: undirected arcs are discarded, an xy arc is appended.");
  dag$arc <- dag$arc[dag$arc.type == 0, ]; # discard undirected arcs;
  nodes <- length(dag$x);
  dag$arc <- rbind(dag$arc, c(1, nodes))   # append xy arc (node 1 -> last node);
    
  b <- b[dag$arc.type == 0];               # discard b's belonging to undirected arcs;
  b <- c(b, bxy);                          # append bxy to b's;
  
  # following two lines just so that the internal dag remains consistent in itself...
  dag$arc.type <- dag$arc.type[dag$arc.type == 0];  # discard arc.type of undirected arcs;
  dag$arc.type <- c(dag$arc.type, 0);               # append arc.type for xy arc;
    
   
  df <- data.frame(matrix(ncol = nodes, nrow = n)); # prepare data.frame to simulate;
  if(!is.na(seed)) set.seed(seed);
  all.simulated <- FALSE;
  counter <- 1;
  
  while (all.simulated == FALSE) {
    for (i in 1:nodes) {
      if (is.na(df[1, i])) {
        if (length(dag.ancestors(dag, A = i)) == 1) {   # does node have no ancestor?;
          if (binary[i]==0) {
            writeLines(paste(c("simulating node", i, "from Normal, with mean SD", mu[i], stdev[i]), collapse = " "));
            df[, i] <- rnorm(n, mean=mu[i], sd=stdev[i]);  # simulate continuous node with no ancestors;
          } else
          if (binary[i]==1) {
            writeLines(paste(c("simulating node", i, "from Binomial, with prob", mu[i]), collapse = " "));
            df[, i] <- rbinom(n=n, size=1, prob=mu[i]); # simulate binary node with no ancestors;
          }
        }
        else {
          if (all(!is.na(df[1, dag.ancestors(dag, i)[-1]]))) { # are all ancestors already simulated?;
            relevant.b <- b[which(dag$arc[, 2] == i)]; # reduce b's to those of arcs pointing to node i;
            relevant.x <- as.matrix(cbind(df[, dag$arc[dag$arc[,2] == i, 1]], rep(0, n)));
                                                       # reduce simulated nodes to ancestors of i;
            
            pred.value <- rowSums(relevant.x %*% c(relevant.b, 0)); # calculate new node "predicted value";
            noise2add  <- rnorm(mean=0, sd=stdev[i], n=n);             # noise to add;
            
            if(binary[i]==0) {
              raw.value <- pred.value + (mu[i]-mean(pred.value)) + noise2add; # shift and add noise!!!;
              df[, i]   <- raw.value;                                         # assign to data.frame;
            
              writeLines(paste(c("node", i,
                                 ": calculated from nodes", dag$arc[dag$arc[,2] == i, 1],
                                 "via arcs", which(dag$arc[, 2] == i),
                                 "with b", relevant.b, ", shifted to mean", mu[i], "(added noise SD:", stdev[i],");"),
                               collapse = " "));
              if (verbose) {
                writeLines("calculation values:");
                print(data.frame(relevant.x[,1:ncol(relevant.x)-1], pred.value, noise2add, raw.value, df[, i]));
              }              
            } else
            if(binary[i]==1) {
              raw.value <- pred.value + (logit(mu[i])-mean(pred.value)) + noise2add;             # shift and add noise!!!;
              df[, i]   <- sapply(raw.value, FUN = function(x)
                                                   {rbinom(n = 1, size = 1, prob = inv.logit(x))}); # assign binary to data.frame;
            
              writeLines(paste(c("node", i,
                                 ": calculated from nodes", dag$arc[dag$arc[,2] == i, 1],
                                 "via arcs", which(dag$arc[, 2] == i),
                                 "with exp(b)", exp(relevant.b),
                                 ", shifted to probability", mu[i], "(added noise SD:", stdev[i],");"), collapse = " "));
              if (verbose) {
                writeLines("calculation values:");
                print(data.frame(relevant.x[,1:ncol(relevant.x)-1], pred.value, noise2add, raw.value, df[, i]));
              }
            }
          }
        }
      }
    }
    counter <- counter + 1; # not really used for anything...;
    all.simulated <- all(!is.na(df[1, ])); # are all nodes already simulated?;
  }
  
  if(naming==2) names(df)<-dag$symbols; # for now, it's either X1,X2,... or the "symbols";
    
  rv <- df
  return(rv)
}
