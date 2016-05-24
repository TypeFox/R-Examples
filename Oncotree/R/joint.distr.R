"joint.distr" <-
  function(otree, with.errors=FALSE,
           edge.weights=if (with.errors) "estimated" else "observed"){
    edge.weights <- match.arg(edge.weights, c("observed","estimated"))
    if (with.errors & is.null(otree$eps)) 
      stop("Need false positive and negative rates")

    mut.names <- otree$parent$child
    N <- otree$nmut #includes root 
    
    # tree-based joint distribution
    # 1. least common ancestor table
    lca.tab <- matrix(NA, nrow=N, ncol=N)
    diag(lca.tab) <- 1:N #lca(v_i, v_i) = v_i
    lca.tab[1, ] <- "Root"    #lca(root, v_i) = root
    lca.tab[ , 1] <- "Root"   #lca(v_i, root) = root
    for (i in 2:N){
      for (j in 2:i){
        lca.tab[i, j] <- least.common.ancestor(otree, mut.names[i], mut.names[j])
        lca.tab[j, i] <- lca.tab[i, j]  #symmetric relationship
      }
    }
  
    # 2. tree-based marginal distribution
    pvec <- marginal.distr(otree, with.errors=FALSE, edge.weights=edge.weights)
    
    # 3. tree-based p_{ij} = p_i * p_j / p_k, where v_k = lca(v_i, v_j)
    pmat <- outer(pvec, pvec)/ matrix(pvec[lca.tab], nrow=N)
    
    # 4. error-adjustment p'_{ij} = p_{ij} (1-epos-eneg)^2 + (p_i+p_j)epos(1-epos-eneg) + epos^2
    if (with.errors){
      epos <- otree$eps["epos"]
      eneg <- otree$eps["eneg"]
      noerr <- 1 - epos - eneg
      pmat <- pmat * noerr^2 + (pvec[row(pmat)] + pvec[col(pmat)]) * epos * noerr + epos^2
      diag(pmat) <- pvec  * noerr + epos
    }
    
    pmat[-1, -1]  # remove the Root
   }
