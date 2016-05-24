####### INTERNAL FUNCTIONS
networkmatrix = function(network, nodes, nettype) {
  if ((nettype %in% c('directed','undirected')) == FALSE) stop("nettype must be either 'directed' or 'undirected' ")
  requireNamespace("igraph", quietly = TRUE)
  # adjacency matrix: first converted to igraph
  if (class(network) == 'matrix' & nettype == 'directed') {
    colnames(network) = nodes
    network = igraph::graph.adjacency(network,"directed", diag=F)
  }
  if (class(network) == 'matrix' & nettype == 'undirected') {
    colnames(network) = nodes
    network = igraph::graph.adjacency(network,"undirected", diag=F)
  }
  # then, get edgelist from igraph
  if (class(network)==c('igraph')) {
    if (is.null(igraph::V(network)$name)) igraph::V(network)$name = nodes
    output = igraph::get.edgelist(network, names = T)
  }
  else if (class(network) == "dgCMatrix") {
    requireNamespace("Matrix", quietly = TRUE)
    indexes = as.matrix(Matrix::summary(network))
    output = matrix(nrow = dim(indexes)[1], ncol = 2)
    v1 = sort(unique(indexes[,1]))
    v2 = sort(unique(indexes[,2]))
    output[,1] = as.character(factor(indexes[,1], labels=nodes[v1]))
    output[,2] = as.character(factor(indexes[,2], labels=nodes[v2]))
  }
  return(output)
}

# compute pvalue
pvalue = function(x, ss, K, N) { #2*min(p>,p<)+p=
  p1 = phyper(x, m = K, n = N-K, k = ss, lower.tail = FALSE) # p > nab
  p2 = phyper(x-1, m = K, n = N-K, k = ss, lower.tail = TRUE) # p < nab
  p3 = 1-p1-p2
  pval = 2*min(p1,p2)+p3
  overenr = (p1<p2)
  return(list('p'=pval,'overenr'=overenr))
}


########## neat CLASS
neatc = function(x,...) UseMethod('neat', x)

neat = function(alist, blist = NULL, network, nettype, nodes, alpha = NULL, anames = NULL, bnames = NULL) {
  if (is.null(alpha) == FALSE) {if ( alpha<=0 | alpha>=1 ) stop('alpha must be in (0,1)')}
  if ((nettype %in% c('directed','undirected')) == FALSE) stop("nettype must be either 'directed' or 'undirected' ")
  if (nettype == 'directed' & is.null(blist)) stop("blist cannot be null when nettype == 'directed' (see manual)")
  # first case: igraph object
  if (class(network) == 'igraph') { 
    requireNamespace("igraph", quietly = TRUE)
    net = networkmatrix(network, nodes, nettype)
  }
  # second case A: adjacency matrix
  else if (class(network) == "matrix" & ncol(network)>2) {
    if (isSymmetric(network)==TRUE & nettype == 'directed') {
      warning('The adjacency matrix is symmetric. Should you set nettype = "undirected"?')
    }
    if (isSymmetric(network)==FALSE & nettype == 'undirected') {
      warning('The adjacency matrix is not symmetric. Should you set nettype = "directed"?')
    }
    net = networkmatrix(network, nodes, nettype)
  }
  # second case B: sparse adjacency matrix (class "dgCMatrix")
  else if (class(network) == "dgCMatrix") {
    requireNamespace("Matrix", quietly = TRUE)
    net = networkmatrix(network, nodes, nettype)
  }
  # third case: two-column matrix with labels
  else if (class(network) == "matrix" & ncol(network)==2) net=network 
  # NB: from this point on, 'net' is the network matrix to be used!!!
  if (is.factor(nodes) == T) {nodes = as.character(nodes)}
  oa = numeric(); ib = numeric(); nab = numeric() 
  alogic = vector("list", length(alist))
  netred = vector("list", length(alist))
  p = numeric(); expect=numeric()
  o = vector(); concl = vector(); from = vector(); to = vector()
  # UNDIRECTED NETWORKS:
  if (nettype == 'undirected') {
    eina = vector("list", length(alist))
    alogic2 = alogic
    if ( is.null( names(alist) ) ) names(alist) = anames
    for (i in 1:length(alist)) {
      if (is.factor(alist[[i]]) == T) {alist[[i]] = as.character(alist[[i]])}
      eina[[i]] = (net[,1] %in% alist[[i]])|(net[,2] %in% alist[[i]])
      if (sum(eina[[i]]) == 1) netred[[i]] = as.matrix(t(net[eina[[i]],]))
      else if (sum(eina[[i]]) == 0) stop(paste('There are no edges connected to genes in',
                            names(alist)[i],'. The test cannot be computed.'))
      else netred[[i]] = net[eina[[i]],]
      alogic[[i]] = (netred[[i]][,1] %in% alist[[i]])
      alogic2[[i]] = (netred[[i]][,2] %in% alist[[i]])
      oa[i] = sum(alogic[[i]])+sum(alogic2[[i]]) # total degree of A (undirected net!)
    }
    # first case: each A vs each B (blist provided!)
    if (is.null(blist) == FALSE) {
      blogic = vector("list", length(blist))
      blogic2 = blogic
      einb = vector("list", length(blist))
      if ( is.null( names(blist) ) ) names(blist) = bnames
      for (i in 1:length(blist)) {
        if (is.factor(blist[[i]]) == T) {blist[[i]] = as.character(blist[[i]])}
        einb[[i]] = (net[,1] %in% blist[[i]])|(net[,2] %in% blist[[i]])
        blogic[[i]] = (net[einb[[i]],][,1] %in% blist[[i]])
        blogic2[[i]] = (net[einb[[i]],][,2] %in% blist[[i]])
        ib[i] = sum(blogic[[i]]) + sum(blogic2[[i]]) # total degree of B (undirected net!)
      }
      k=1
      for (i in 1:length(alist)) {
        for (j in 1:length(blist)) {
          from[k] = names(alist)[i]
          to[k] = names(blist)[j]
          link = (alogic[[i]] & netred[[i]][,2] %in% blist[[j]])|(alogic2[[i]] & netred[[i]][,1] %in% blist[[j]])
          nab[k] = sum(link)
          temp = pvalue(nab[k], ss = oa[i], K = ib[j], N = 2*dim(net)[1])
          p[k] = temp$p
          expect[k] = round(ib[j] * oa[i] / (2*dim(net)[1]), 4)
          o[k] = temp$overenr
          if (is.null(alpha) == FALSE) {
            if (p[k] > alpha) concl[k] = 'No enrichment'
            else if (p[k] <= alpha & o[k]==T) concl[k] = 'Overenrichment'
            else if (p[k] <= alpha & o[k]==F) concl[k] = 'Underenrichment'
          }
          k = k+1
        }
      }
    }
    # second case: blist NOT provided
    if (is.null(blist) == TRUE) {
      ib = oa # undirected: A = B (and no distinction between indegree and outdegree)
      k=1
      for (i in 1:(length(alist)-1)) {
        for (j in (i+1):length(alist)) {
          from[k] = names(alist)[i]
          to[k] = names(alist)[j]
          link = (alogic[[i]] & netred[[i]][,2] %in% alist[[j]])|(alogic2[[i]] & netred[[i]][,1] %in% alist[[j]])
          nab[k] = sum(link)
          temp = pvalue(nab[k], ss = oa[i], K = ib[j], N = 2*dim(net)[1])
          p[k] = temp$p
          expect[k] = round(ib[j] * oa[i] / (2*dim(net)[1]), 4)
          o[k] = temp$overenr
          if (is.null(alpha) == FALSE) {
            if (p[k] > alpha) concl[k] = 'No enrichment'
            else if (p[k] <= alpha & o[k]==T) concl[k] = 'Overenrichment'
            else if (p[k] <= alpha & o[k]==F) concl[k] = 'Underenrichment'
          }
          k = k+1
        }
      }
    }
  } # end for part for undirected networks
  # DIRECTED NETWORKS:
  if (nettype == 'directed') {
    if ( is.null( names(alist) ) ) names(alist) = anames
    for (i in 1:length(alist)) {
      if (is.factor(alist[[i]]) == T) {alist[[i]] = as.character(alist[[i]])}
      alogic[[i]] = (net[,1] %in% alist[[i]])
      oa[i] = sum(alogic[[i]]) # outdegree of A
      if ( sum(alogic[[i]]) ==1 ) netred[[i]] = as.matrix(t(net[alogic[[i]],]))
      else if (sum(alogic[[i]]) == 0) stop(paste('There are no arrows going out from genes in',
                           names(alist)[i]),'. The test cannot be computed.')
      else netred[[i]] = net[alogic[[i]],]
    }
    # directed: only relevant case is each A vs B (blist provided!)
    if (is.null(blist) == FALSE) {
      blogic = vector("list", length(blist))
      if ( is.null( names(blist) ) ) names(blist) = bnames
      for (i in 1:length(blist)) {
        if (is.factor(blist[[i]]) == T) {blist[[i]] = as.character(blist[[i]])}
        blogic[[i]] = (net[,2] %in% blist[[i]])
        ib[i] = sum(blogic[[i]]) # indegree of B
      }
      k=1
      for (i in 1:length(alist)) {
        for (j in 1:length(blist)) {
          from[k] = names(alist)[i]
          to[k] = names(blist)[j]
          nab[k] = sum( netred[[i]][,2] %in%  blist[[j]])
          temp = pvalue(nab[k], ss = oa[i], K = ib[j], N = dim(net)[1])
          p[k] = temp$p
          expect[k] = round(ib[j] * oa[i] / dim(net)[1], 4)
          o[k] = temp$overenr
          if (is.null(alpha) == FALSE) {
            if (p[k] > alpha) concl[k] = 'No enrichment'
            else if (p[k] <= alpha & o[k]==T) concl[k] = 'Overenrichment'
            else if (p[k] <= alpha & o[k]==F) concl[k] = 'Underenrichment'
          }
          k = k+1
        }
      }
    }
  } # end of part for directed networks
  # final common code
  if (is.null(alpha) == FALSE) {
    results = data.frame(from, to, nab, expect, p, concl)
    names(results) = c('A', 'B', 'nab', 'expected_nab', 'pvalue', 'conclusion')
  }
  else {
    results = data.frame(from, to, nab, expect, p)
    names(results) = c('A', 'B', 'nab', 'expected_nab', 'pvalue')
  }
  class(results)=c('neat','data.frame')
  results
}

# summary, plot and print:
summary.neat = function(object, ...) {
  cat("Number of comparisons:",dim(object)[1],"\n")
  cat("Enrichments at 1% level:",sum(object$pvalue<0.01),"\n")
  cat("Enrichments at 5% level:",sum(object$pvalue<0.05),"\n")
  cat("Kolmogorov-Smirnov test for uniformity of p-values:",round(ks.test(x=object$pvalue, y='punif')$p.value,4),"\n")
}

plot.neat = function(x, nbreaks=10, ...) {
  par(mfrow=c(1,2))
  hist(x$pvalue,breaks=nbreaks, xlim=c(0,1), xlab='p-value', main='Histogram of p-values')
  plot(ecdf(x$pvalue),main='P-p plot of p-values',xlab='x',ylab ='F(x)',xlim=c(0.037,0.963),ylim=c(0.037,0.963))
  abline(v=1);abline(a=0,b=1,col='red')
  legend(0.6,0.2,c('Distribution of','p-values','Uniform','distribution'),col=c('black','white','red','white'),lty=1,cex=0.6)
  }

print.neat = function(x,nrows=10,...) {
  class(x) = 'data.frame'
  if (nrows == 'ALL') {
    nrows = dim(x)[1]
  }
  head(x,nrows)
}


