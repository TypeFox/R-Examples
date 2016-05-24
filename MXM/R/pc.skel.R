
pc.skel <- function(dataset, method = "pearson", alpha = 0.05, rob = FALSE, R = 1, graph = FALSE) {
  ## dataset contains the data, it must be a matrix 
  ## type can be either "pearson" or "spearman" for continuous variables OR
  ## "cat" for categorical variables
  ## alpha is the level of significance, set to 0.05 by default
  ## rob is TRUE or FALSE and is supported by type = "pearson" only, i.e. it is to be used
  ## for robust estimation of Pearson correlation coefficient only
  ## if graph is true, the graph will appear

  alpha <- log(alpha)
  title <- deparse( substitute(dataset) )

  ### if you want to use Spearman, simply use Spearman on the ranks
  if (method == "spearman")  {
    dat = apply(dataset, 2, rank)
    rob = FALSE
  }
  
  if (method == "spearman" || method == "pearson") {
    if (R == 1) {    
      ci.test = condi 
      type = method
      rob = rob
      R = 1
    } else if (R > 1) {
      ci.test = condi.perm
    }
    
  } else {
    ci.test = cat.ci
    type = NULL
    rob = FALSE
  }
  
  n = ncol(dataset)
  m = nrow(dataset)
  k <- 0  ## initial size of the conditioning set
  G = matrix(2, n, n)  # 3 sep-set indicates a subset of variables which eliminate given edge  
  ## If an element has the number 2 it means there is connection, otherwiser it will have 0
  diag(G) = -100
  durat = proc.time()
  if ( method == "pearson" || method == "spearman") {
    
    if (R == 1) {
      if ( rob == FALSE ) {
        r = cor(dataset)
        if (type == "pearson") {
          stata = abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m - 3) )  ## absolute of the test statistic
        } else if (type == "spearman") {
          stata = abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m - 3) ) / 1.029563  ## absolute of the test statistic
        }
        pv = pvalue = log(2) + pt(stata, m - 3, lower.tail = FALSE, log.p = TRUE)  ## logged p-values 
        dof = matrix(m - 3, n, n)
        diag(stata) = diag(dof) = 0
        stadf = stata / dof

      } else {

        stat = pv = matrix(0, n, n)
        for ( i in 1:c(n - 1) ) {
          for ( j in c(i + 1):n ) {
            ro <- condi(i, j, 0, dataset, type = "pearson", rob = TRUE) 
            stat[i, j] = ro[1]
            pv[i, j] = ro[2]
          }
        }
        pvalue = pv + t(pv)  ## p-values
        stata = stat + t(stat)
        dof = matrix(m - 3, n, n)
        stadf = stata / dof
      }
      
    } else if (R > 1) {
      stat = pv = matrix(0, n, n)
      for ( i in 1:c(n - 1) ) {
        for ( j in c(i + 1):n ) {
          ro <- permcor(dataset[, c(i, j)], R = R) 
          stat[i, j] = ro$result[1]
          pv[i, j] = log( ro$result[2] )
        }
      }
      pvalue = pv + t(pv)  ## p-values
      stata = stat + t(stat)
      dof = matrix(m - 3, n, n)
      stadf = stata / dof
    }   

  } else { ## type = cat
    stat = pv = dof = matrix(0, n, n)
    for ( i in 1:c(n - 1) ) {
      for ( j in c(i + 1):n ) {
        ro <- cat.ci(i, j, 0, dataset) 
        stat[i, j] = ro[1]
        pv[i, j] = ro[2]
        dof[i, j] = ro[3] 
      }
    }
    pvalue = pv + t(pv)  ## p-values
    stata = stat + t(stat)
    diag(stata) <- 0
    dof = dof + t(dof)  ## p-values
    stadf = stata / dof
  }
    
  pv = pvalue 

  #stat[ lower.tri(stat) ] = 2
  pv[ lower.tri(pv) ] = 2 
  G[pvalue > alpha] <- 0   ## removes edges from non significantly related pairs
  if ( is.null( colnames(dataset) ) ) {
    colnames(G) = rownames(G) = paste("X", 1:n, sep = "")
  } else colnames(G) = rownames(G) = colnames(dataset)
  diag(pvalue) = diag(pv) = 0
  ina = 1:n 
  sep = list()
  n.tests = NULL
  pval = list()
  #### some more initial stuff 
   dial = which( pv <= alpha, arr.ind = T )
   zeu = cbind( dial, stadf[ dial ], pv[ dial ] )  ## all significant pairs of variables
   zeu = zeu[ order( - zeu[, 4], zeu[, 3] ), ] ## order of the pairs based on their strength
   if ( !is.matrix(zeu) ) zeu = matrix(zeu, nrow = 1)
   duo = nrow(zeu)  ## number of pairs to be checkd for conditional independence
   n.tests[1] = n * (n - 1) / 2

  #### main search

  if (duo == 0) {
    diag(G) = 0
    final = list(kappa = k, G = G) 
  } else {

    ell = 2

    ## Execute PC algorithm: main loop

    while ( k < ell & nrow(zeu) > 0 )  {
      k = k + 1   ## size of the seperating set will change now
      tes = 0
      met = matrix(nrow = nrow(zeu), ncol = k + 2)

      for ( i in 1:nrow(zeu) ) {

        adjx = ina[ G[ zeu[i, 1], ] == 2 ]  ;  lx = length(adjx)  ## adjacents to x
        adjy = ina[ G[ zeu[i, 2], ] == 2 ]  ;  ly = length(adjy)  ## adjacents to y
        if ( lx >= k )  {
          pvalx = pvalue[ zeu[i, 1], adjx ]
          infox = cbind( adjx, pvalx)
          infox = infox[ order( - pvalx ), ]
          if ( !is.matrix(infox) ) {
            samx = cbind( infox[1], infox[2] )
          } else  samx = cbind( t( combn(infox[, 1], k) ), t( combn(infox[, 2], k) ) )  ## factorial, all possible unordered pairs
        } 
        if ( ly >= k ) {
          pvaly = pvalue[ zeu[i, 2], adjy ]
          infoy = cbind(adjy, pvaly)
          infoy = infoy[ order( - pvaly ), ]
          if ( !is.matrix(infoy) ) {
            samy = cbind( infoy[1], infoy[2] )
          } else  samy = cbind( t( combn(infoy[, 1], k) ), t( combn(infoy[, 2], k) ) )  ## factorial, all possible unordered pairs
        }
        if ( !is.null(samx) ) sx = 1  else sx = 0
        if ( !is.null(samy) ) sy = 1  else sy = 0 
        sam = rbind( samx * sx, samy * sy ) 
        sam = as.matrix(sam)
        sam = unique(sam) 
        ## sam contains either the sets of k neighbours of X, or of Y or of both
        ## if the X and Y have common k neighbours, they are removed below
        rem = intersect( zeu[i, 1:2], sam )
        if ( length(rem) > 0 ) {
          pam = list()
          for ( j in 1:length(rem) ) {
            pam[[ j ]] = as.vector( which(sam == rem[j], arr.ind = TRUE)[, 1] ) 
          }
        }

        pam = unlist(pam)
        sam = sam[ - pam, ] 

        if ( !is.matrix(sam) ) {
          sam = matrix( sam, nrow = 1 ) 
        } else if ( nrow(sam) == 0 ) {
          G = G 
        } else { 
          if (k == 1) {
            sam = sam[ order( sam[, 2 ] ), ]
          } else {
            an <- t( apply(sam[, -c(1:2)], 1, sort, decreasing = TRUE) )
            sam <- cbind(sam[, 1:2], an)
            nc <- ncol(sam)
            sam2 <- as.data.frame( sam[, nc:1] )     
            sam2 <- sam2[ do.call( order, as.list( sam2 ) ), ] 
            sam <- as.matrix( sam2[, nc:1] )
          }
        }

        if ( nrow(sam) == 0 ) {
          G = G  
        } else {
          a = ci.test( zeu[i, 1], zeu[i, 2], sam[1, 1:k], dataset, type = type, rob = rob, R = R )
          if ( a[2] > alpha ) {
            G[ zeu[i, 1], zeu[i, 2] ] = 0  ## remove the edge between two variables
            G[ zeu[i, 2], zeu[i, 1] ] = 0  ## remove the edge between two variables 
            met[i, ] = c( sam[1, 1:k], a[1:2] )
            tes = tes + 1 
          } else {
            m = 1
            while ( a[2] < alpha  &  m < nrow(sam) ) {
              m = m + 1
              a = ci.test( zeu[i, 1], zeu[i, 2], sam[m, 1:k], dataset, type = type, rob = rob, R = R )
              tes = tes + 1
            }
            if (a[2] > alpha) {
              G[ zeu[i, 1], zeu[i, 2] ] = 0  ## remove the edge between two variables
              G[ zeu[i, 2], zeu[i, 1] ] = 0  ## remove the edge between two variables
              met[i, ] = c( sam[m, 1:k], a[1:2] ) 
            }
          }
        }

        sam = samx = samy = NULL

      }  
        ax = ay = list()
        lx = ly = numeric( nrow(zeu) )
        for ( i in 1:nrow(zeu) ) {
          ax[[ i ]] = ina[ G[ zeu[i, 1], ] == 2 ]  ;  lx[i] = length( ax[[ i ]] )
          ay[[ i ]] = ina[ G[ zeu[i, 2], ] == 2 ]  ;  ly[i] = length( ay[[ i ]] ) 
        }

      ell = max(lx, ly)
      id = which( rowSums(met) > 0 )
      if (length(id) == 1) {
         sep[[ k ]] = c( zeu[id, 1:2], met[id, ] )
      } else {
        sep[[ k ]] = cbind( zeu[id, 1:2], met[id, ] )
      }
      zeu = zeu[-id, ]  
      if ( class(zeu) != "matrix" ) {
        zeu <- matrix(zeu, ncol = 4)
      }
      n.tests[ k + 1 ] = tes

    }

    G <- G / 2
    diag(G) = 0
    durat = proc.time() - durat

    ###### end of the algorithm

    for ( l in 1:k ) { 
      if ( is.matrix(sep[[ l ]]) ) {
        if ( nrow(sep[[ l ]]) > 0) {  
          sepa = sep[[ l ]]
          colnames( sepa )[1:2] = c("X", "Y")
          colnames( sepa )[ 2 + 1:l ] = paste("SepVar", 1:l)
          colnames( sepa )[ c(l + 3):c(l + 4) ] = c("stat", "logged.p-value")
          sepa =  sepa[ order(sepa[, 1], sepa[, 2] ), ]
          sep[[ l ]] = sepa
        }
      } else {
        if ( length(sep[[ l ]]) > 0) { 
          names( sep[[ l ]] )[1:2] = c("X", "Y")
          names( sep[[ l ]] )[ 2 + 1:l ] = paste("SepVar", 1:l)
          names( sep[[ l ]] )[ c(l + 3):c(l + 4) ] = c("stat", "logged.p-value")
        } 
      }
    }
  
  }  

  n.tests = n.tests[ n.tests>0 ]
  k = length(n.tests) - 1
  sepset = list()
  if (k == 0) {
    sepset = NULL
  } else {
    for ( l in 1:k ) {
      if ( is.matrix( sep[[ l ]] ) )  {
        nu <- nrow( sep[[ l ]] )
        if ( nu > 0 ) sepset[[ l ]] = sep[[ l ]][1:nu, ]
      } else sepset[[ l ]] = sep[[ l ]]    
    }
  }  
  names(n.tests) = paste("k=", 0:k, sep ="")
  
  info <- summary( rowSums(G) )
  density <- sum(G) / ( n * ( n - 1 ) )
  
 if(graph == TRUE)
  {
    if(requireNamespace("Rgraphviz", quietly = TRUE, warn.conflicts = FALSE) == TRUE)
    {
      am.graph <- new("graphAM", adjMat = G, edgemode = "undirected")
      plot( am.graph, main = paste("Skeleton of the PC algorithm for", title ) ) 
    }else{
      warning('In order to plot the generated network, package Rgraphviz is required.')
    }
  }

  list(stat = stata, pvalue = pvalue, runtime = durat, kappa = k, n.tests = n.tests, density = density, info = info, G = G, sepset = sepset, title = title )
}