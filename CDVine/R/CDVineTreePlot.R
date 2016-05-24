#################################################
# CDvineTreePlot  			#
#						#
# Description: Plots all trees of a C-vine	#
# Input:					#
# data		data matrix			#
# family	copula families			#
# par		copula parameters		#
# par2		second copula parameters	#
# type		vine type			#
# tree		tree number or all trees	#
# Output:					#
# Plots; C-vine or D-vine trees			#
#################################################

CDVineTreePlot <- function(data = NULL, family, par = rep(0, length(family)), par2 = rep(0, length(family)), 
  names = NULL, type, method = "mle", max.df = 30, max.BB = list(BB1 = c(5, 6), BB6 = c(6, 6), BB7 = c(5, 
    6), BB8 = c(6, 1)), tree = "ALL", edge.labels = c("family"), P = NULL, ...) {
  if (type == "CVine") 
    type <- 1 else if (type == "DVine") 
    type <- 2
  if (type != 1 & type != 2) 
    stop("Vine model not implemented.")
  
  if (edge.labels[1] != FALSE & !all(edge.labels %in% c("family", "par", "par2", "theotau", "emptau"))) 
    stop("Edge label not implemented.")
  if (is.null(data) & any(edge.labels == "emptau")) 
    stop("Empirical Kendall's tau values cannot be obtained if no data is provided.")
  if (is.null(data) != TRUE && any(data > 1) || any(data < 0)) 
    stop("Data has be in the interval [0,1].")
  
  dd <- length(family)
  d <- (1 + sqrt(1 + 8 * dd))/2
  anz_trees <- d - 1
  
  if (max.df <= 2) 
    stop("The upper bound for the degrees of freedom parameter has to be larger than 2.")
  if (!is.list(max.BB)) 
    stop("'max.BB' has to be a list.")
  if (max.BB$BB1[1] < 0.001) 
    stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
  if (max.BB$BB1[2] < 1.001) 
    stop("The upper bound for the second parameter of the BB1 copula should be greater than 1.001 (lower bound for estimation).")
  if (max.BB$BB6[1] < 1.001) 
    stop("The upper bound for the first parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
  if (max.BB$BB6[2] < 1.001) 
    stop("The upper bound for the second parameter of the BB6 copula should be greater than 1.001 (lower bound for estimation).")
  if (max.BB$BB7[1] < 1.001) 
    stop("The upper bound for the first parameter of the BB7 copula should be greater than 1.001 (lower bound for estimation).")
  if (max.BB$BB7[2] < 0.001) 
    stop("The upper bound for the second parameter of the BB7 copula should be greater than 0.001 (lower bound for estimation).")
  if (max.BB$BB8[1] < 1.001) 
    stop("The upper bound for the first parameter of the BB1 copula should be greater than 0.001 (lower bound for estimation).")
  if (max.BB$BB8[2] < 0.001 || max.BB$BB8[2] > 1) 
    stop("The upper bound for the second parameter of the BB1 copula should be in the interval [0,1].")
  
  if (!is.null(names) & length(names) != d) 
    stop("Number of variable names incorrect.")
  
  if (is.null(P)) {
    P <- list()
    for (i in 1:anz_trees) P[[i]] <- 0
  }
  
  if (d < 3) 
    stop("Dimension has to be at least 2.")
  if (length(family) != d * (d - 1)/2) 
    stop("Number of copula families incorrect.")
  if (length(par) != d * (d - 1)/2) 
    stop("Number of copula parameters incorrect.")
  if (length(par2) != d * (d - 1)/2) 
    stop("Number of second copula parameters incorrect.")
  
  if (tree != "ALL" && tree > anz_trees) 
    stop("Selected tree does not exist.")
  
  for (i in 1:(d * (d - 1)/2)) {
    # fam richtig?
    if (!(family[i] %in% c(0, 1:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40))) 
      stop("Copula family not implemented.")
    # Parameterbereiche abfragen
    if (any(par != 0)) {
      if ((family[i] == 1 || family[i] == 2) && abs(par[i]) >= 1) 
        stop("The parameter of the Gaussian and t-copula has to be in the interval (-1,1).")
      if (family[i] == 2 && par2[i] <= 2) 
        stop("The degrees of freedom parameter of the t-copula has to be larger than 2.")
      if ((family[i] == 3 || family[i] == 13) && par[i] <= 0) 
        stop("The parameter of the Clayton copula has to be positive.")
      if ((family[i] == 4 || family[i] == 14) && par[i] < 1) 
        stop("The parameter of the Gumbel copula has to be in the interval [1,oo).")
      if ((family[i] == 6 || family[i] == 16) && par[i] <= 1) 
        stop("The copula parameter of the Joe copula has to be in the interval (1,oo).")
      if (family[i] == 5 && par[i] == 0) 
        stop("The parameter of the Frank copula has to be unequal to 0.")
      if ((family[i] == 7 || family[i] == 17) && par[i] <= 0) 
        stop("The first parameter of the BB1 copula has to be positive.")
      if ((family[i] == 7 || family[i] == 17) && par2[i] < 1) 
        stop("The second parameter of the BB1 copula has to be in the interval [1,oo).")
      if ((family[i] == 8 || family[i] == 18) && par[i] <= 0) 
        stop("The first parameter of the BB6 copula has to be in the interval [1,oo).")
      if ((family[i] == 8 || family[i] == 18) && par2[i] < 1) 
        stop("The second parameter of the BB6 copula has to be in the interval [1,oo).")
      if ((family[i] == 9 || family[i] == 19) && par[i] < 1) 
        stop("The first parameter of the BB7 copula has to be in the interval [1,oo).")
      if ((family[i] == 9 || family[i] == 19) && par2[i] <= 0) 
        stop("The second parameter of the BB7 copula has to be positive.")
      if ((family[i] == 10 || family[i] == 20) && par[i] < 1) 
        stop("The first parameter of the BB8 copula has to be in the interval [1,oo).")
      if ((family[i] == 10 || family[i] == 20) && (par2[i] <= 0 || par2[i] > 1)) 
        stop("The second parameter of the BB8 copula has to be in the interval (0,1].")
      if ((family[i] == 23 || family[i] == 33) && par[i] >= 0) 
        stop("The parameter of the rotated Clayton copula has to be negative.")
      if ((family[i] == 24 || family[i] == 34) && par[i] > -1) 
        stop("The parameter of the rotated Gumbel copula has to be in the interval (-oo,-1].")
      if ((family[i] == 26 || family[i] == 36) && par[i] >= -1) 
        stop("The parameter of the rotated Joe copula has to be in the interval (-oo,-1).")
      if ((family[i] == 27 || family[i] == 37) && par[i] >= 0) 
        stop("The first parameter of the rotated BB1 copula has to be negative.")
      if ((family[i] == 27 || family[i] == 37) && par2[i] > -1) 
        stop("The second parameter of the rotated BB1 copula has to be in the interval (-oo,-1].")
      if ((family[i] == 28 || family[i] == 38) && par[i] >= 0) 
        stop("The first parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
      if ((family[i] == 28 || family[i] == 38) && par2[i] > -1) 
        stop("The second parameter of the rotated BB6 copula has to be in the interval (-oo,-1].")
      if ((family[i] == 29 || family[i] == 39) && par[i] > -1) 
        stop("The first parameter of the rotated BB7 copula has to be in the interval (-oo,-1].")
      if ((family[i] == 29 || family[i] == 39) && par2[i] >= 0) 
        stop("The second parameter of the rotated BB7 copula has to be negative.")
      if ((family[i] == 30 || family[i] == 40) && par[i] > -1) 
        stop("The first parameter of the rotated BB8 copula has to be in the interval (-oo,-1].")
      if ((family[i] == 30 || family[i] == 40) && (par2[i] >= 0 || par2[i] < (-1))) 
        stop("The second parameter of the rotated BB8 copula has to be in the interval [-1,0).")
    }
  }
  
  theoTauVec <- numeric()
  empTauVec <- numeric()
  parVec <- numeric()
  par2Vec <- numeric()
  seqpar <- list()
  if (!is.null(data)) {
    seqpar <- CDVineSeqEst(data, family, method = method, max.df = max.df, max.BB = max.BB, type = type)
    
    for (i in 1:dd) {
      empTauVec[i] <- round(BiCopPar2Tau(family[i], seqpar$par[i], seqpar$par2[i]), 2)
    }
  }
  
  if (!all(par == 0)) {
    for (i in 1:dd) {
      theoTauVec[i] <- round(BiCopPar2Tau(family[i], par[i], par2[i]), 2)
    }
  }
  
  if (any(edge.labels == "par")) {
    if (!is.null(data) && all(par == 0) == TRUE) 
      parVec <- round(seqpar$par, 2) else parVec <- round(par, 2)
  }
  if (any(edge.labels == "par2")) {
    if (!is.null(data) && all(par == 0) == TRUE) 
      par2Vec <- round(seqpar$par, 2) else par2Vec <- round(par2, 2)
  }
  
  
  
  if (tree != "ALL" && tree > d - 1) 
    stop("Selected tree does not exist.")
  
  if (tree == "ALL") 
    tree <- 1:(d - 1)
  
  
  weight <- numeric()
  for (j in 1:dd) weight[j] <- ifelse(is.null(data), theoTauVec[j], empTauVec[j])
  
  if (edge.labels[1] != FALSE) {
    numlabels <- length(edge.labels)
    elabels <- list()
    for (j in 1:dd) elabels[[j]] <- rep(NA, numlabels)
  }
  
  # initial edge
  if (edge.labels[1] != FALSE) {
    for (jj in 1:numlabels) {
      for (j in 1:dd) {
        if (edge.labels[jj] == "family") 
          elabels[[j]][jj] <- BiCopName(family[j])
        if (edge.labels[jj] == "par") 
          elabels[[j]][jj] <- parVec[j]
        if (edge.labels[jj] == "par2") 
          elabels[[j]][jj] <- par2Vec[j]
        if (edge.labels[jj] == "theotau") 
          elabels[[j]][jj] <- theoTauVec[j]
        if (edge.labels[jj] == "emptau") 
          elabels[[j]][jj] <- empTauVec[j]
      }
    }
  }
  
  
  elabels2 <- rep("", dd)
  if (edge.labels[1] != FALSE) {
    for (j in 1:dd) {
      if (numlabels > 1) {
        elabels2[j] <- paste(elabels[[j]], collapse = ",")
      } else elabels2[j] <- elabels[[j]]
    }
  }
  
  
  nodes <- numeric()
  
  if (!is.null(names)) 
    no <- names else if (!is.null(data) && !is.null(colnames(data))) 
    no <- colnames(data) else {
    varnames <- rep(0, d)
    for (i in 1:d) varnames[i] <- paste("V", i, sep = "")
    no <- varnames
  }
  
  nodes <- node_names2(no, type)
  
  
  # Plotten
  
  for (t in tree) {
    anz <- length(nodes[[t]])
    adjmatrix <- matrix(0, ncol = anz, nrow = anz)
    
    kk <- 1
    if (t > 1) {
      for (o in 1:(t - 1)) kk <- kk + d - o
    }
    
    if (type == 1) {
      for (i in 2:anz) adjmatrix[1, i] <- weight[(kk + i - (t + 1))]
    } else if (type == 2) {
      for (i in 1:anz - 1) adjmatrix[i, i + 1] <- weight[(kk + i - t)]
    } else {
      stop("This vine type is not implemented.")
    }
    
    colnames(adjmatrix) <- nodes[[t]]
    
    g <- graph.adjacency(adjmatrix, weighted = TRUE, diag = FALSE)
    
    if (edge.labels[1] != FALSE) {
      elabel <- elabels2[[t]]
    } else {
      elabel <- NULL
    }
    
    if (all(P[[t]] == 0)) {
      P[[t]] <- layout.fruchterman.reingold(g)
    }
    
    if (!exists("main")) 
      main <- paste("Tree ", t, sep = "") else {
      if (main != paste("Tree ", (t - 1), sep = "")) 
        main <- paste("Tree ", t, sep = "")
    }
    
    plot(g, layout = P[[t]], vertex.label = V(g)$name, vertex.shape = "rectangle", vertex.size = max(strwidth(V(g)$name, 
      units = "figure")) * 500, vertex.label.family = "sans", edge.label.family = "sans", edge.label = elabels2[kk:(kk + 
      d - t - 1)], edge.width = (10 * abs(weight[kk:(kk + d - t - 1)]) + 0.5), edge.arrow.size = 0, 
      main = main)
    
    if (t != max(tree)) {
      par(ask = TRUE)
    } else {
      par(ask = FALSE)
    }
    
  }
  
  return(P)
}


#################################################################
# node_names  Function to name the nodes in a C- or D-vine	#
#								#
# Input:	col_names,type					#
# Output	names						#
#################################################################


node_names2 <- function(col_names, type) {
  d <- length(col_names)
  names <- list()
  
  # Fuer den C-vine
  if (type == 1) {
    # Der erste Baum mit seinen Konten ist bereits festgelegt und man bekommt in col_names einen Vektor mit
    # den Knoten des ersten Baumes
    names[[1]] <- col_names
    # Zweiter Baum hat keinen Bedingungsstrich '|'
    node_tree2 <- numeric()
    for (i in 1:(d - 1)) {
      node_tree2[i] <- paste(col_names[1], col_names[i + 1], sep = ",")
    }
    names[[2]] <- node_tree2
    
    # Die weiteren Baeume
    for (t in 2:(d - 1)) {
      nodes <- numeric()
      if (t > 2) {
        nachStrich <- col_names[1]
        for (j in 2:(t - 1)) nachStrich <- paste(nachStrich, col_names[j], sep = ",")
      } else nachStrich <- col_names[1]
      for (i in 1:(d - t)) {
        vorStrich <- paste(col_names[t], col_names[(t + i)], sep = ",")
        nodes[i] <- paste(vorStrich, nachStrich, sep = "|")
      }
      
      names[[(t + 1)]] <- nodes
    }
  } else if (type == 2) {
    # D-vine
    names[[1]] <- col_names
    # Zweiter Baum hat keinen Bedingungsstrich '|'
    node_tree2 <- numeric()
    for (i in 1:(d - 1)) {
      node_tree2[i] <- paste(col_names[i], col_names[i + 1], sep = ",")
    }
    names[[2]] <- node_tree2
    
    # Die weiteren Baeume
    for (t in 2:(d - 1)) {
      nodes <- numeric()
      for (i in 1:(d - t)) {
        vorStrich <- paste(col_names[i], col_names[(t + i)], sep = ",")
        nachStrich <- col_names[(i + 1)]
        if (t > 2) {
          for (j in (i + 2):(i + t - 1)) nachStrich <- paste(nachStrich, col_names[j], sep = ",")
        }
        nodes[i] <- paste(vorStrich, nachStrich, sep = "|")
      }
      
      names[[(t + 1)]] <- nodes
    }
  }
  
  
  return(names)
} 
