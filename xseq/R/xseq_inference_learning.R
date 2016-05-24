# xseq inference and learning algorithms
# 
# Date: 
#   Revised: February 25, 2015
#
# Author: Jiarui Ding <jiaruid@cs.ubc.ca>
#   Department of Computer Science, UBC
#   Department of Molecular Oncology, BC Cancer Agency 
#

#==============================================================================
LogPlus = function(x, y) {
  # Compute the log of pairwise sum of two log transformed vectors
  # 
  # Args:
  #   x: One of the two log-transformed vectors
  #   y: The other log-transformed vector. x and y must have the same length 
  #       and without missing values
  #       
  # Returns:
  #   The log of the pair-wise sum of x and y
  
  pmax(x, y) + log1p(exp(-abs(x - y)))
}

#' @useDynLib xseq sexp_log_sum_exp
#' @importFrom stats na.omit
LogSumExp = function(x) { 
  # The log-sum-exp trick, x is a vector. 
  x = as.numeric(x)
  x = na.omit(x)
  if (length(x) < 1) {
    stop("Error: x should be a vector of length greater than zero!")
  }
  
  res = .Call("sexp_log_sum_exp", x)
  
  return(res)
}

# Sample data from a dirichlet distribution
#' @importFrom stats rgamma
rdirichlet = function(n, alpha) {
  len = length(alpha)  
  x = matrix(rgamma(len*n, alpha), ncol=len, byrow=TRUE)
  
  return(x / rowSums(x))
}

#' @importFrom stats setNames
ReplicateData = function(data, time, name, temporal=FALSE) {
  if (time <= 0) {
    if(is.matrix(data)) {
      ept.list = setNames(vector("list", nrow(data)), rownames(data))
    } else {
      ept.list = NULL
    }
    return(ept.list)
  }
  
  nrow.data = nrow(data)
  if(temporal == TRUE) {
    data.back = data
    if(length(nrow.data)) { ## CPTs
      data = setNames(vector("list", nrow.data), rownames(data.back))
      for (i in 1:nrow.data) {
        tmp = replicate(time, data.back[i, ], simplify=FALSE)
        tmp = do.call(rbind, tmp)
        
        rownames(tmp) = name
        data[[i]] = tmp
      }
    } else { ## Priors
      data = replicate(time, data, simplify=FALSE)
      names(data) = name
      data = do.call(rbind, data)
    }
  }
  else { ## Priors
    data = replicate(time, data, simplify=FALSE)
    names(data) = name
  }
  
  return(data)
}


#==============================================================================
#' The datastructure to store the xseq models
#' 
#' @export
#' @param mut A data.frame of mutations. The data.frame should have 
#' three columns of characters: sample, hgnc_symbol, and variant_type. 
#' The variant_type column cat be either "HOMD", "HLAMP", 
#'   "MISSENSE", "NONSENSE", "FRAMESHIFT", "INFRAME", "SPLICE", "NONSTOP", 
#'   "STARTGAINED", "SYNONYMOUS", "OTHER", "FUSION", "COMPLEX".
#' @param expr A matrix of gene expression values where 
#' each row corresponds to a patient and each column is a gene
#' @param net A list of gene interaction networks
#' @param expr.dis The fitted gene expression distributions, 
#' output from \code{GetExpressionDistribution}
#' @param prior The prior for xseq, output from \code{SetXseqPrior}
#' @param cpd A list of conditional probability tables for xseq, 
#' output from \code{SetXseqPrior}
#' @param gene A character vector of gene names, 
#' default to all the genes with mutations
#' @param p.h The down-regulation probability list of each gene 
#' connected to a mutated gene, 
#' typically from running \code{LearnXseqParameter} 
#' on a discovery dataset
#' @param weight The weight list of each gene 
#' connected to a mutated gene, 
#' typically from running \code{LearnXseqParameter} 
#' on a discovery dataset
#' @param cis Logical, cis or trans analysis
#' @param debug Logical, whether to output debug information
#' @return A xseq model
#' 
#' @importFrom stats setNames median
#' 
InitXseqModel = function(mut, expr, net, expr.dis, prior, cpd, gene, 
                         p.h, weight, cis=FALSE, debug=FALSE) {
  
  ## Remove the cases without expression data (may add back in future)
  if (!missing(expr.dis)) {
    sample.expr = rownames(expr.dis[[1]])
    gene.expr   = colnames(expr.dis[[1]]) 
    
    id = mut[, "sample"] %in% sample.expr
    mut = mut[id, , drop=FALSE]
  }
  
  sample.mut = unique(mut[, "sample"])
  gene.mut   = unique(mut[, "hgnc_symbol"])
  
  if (missing(gene)) {
    gene = gene.mut
  } else {
    gene = intersect(gene, gene.mut)
  }
  
  if (length(gene) < 1 | length(sample.mut) < 1) {
    stop("Error: the number of samples or genes is smaller than 1!") 
  }
  
  ## 
  p.d = ReplicateData(cpd$p.d.share, length(gene), gene)
  
  p.df  = vector("list", length(gene))
  names(p.df)  = gene
  
  p.fg  = vector("list", length(gene))
  names(p.fg)  = gene
  
  if (!missing(expr.dis)) {
    p.gy  = vector("list", length(gene))
    names(p.gy)  = gene
    
    ##
    lambda = expr.dis$lambda
    weight.sample = t(expr.dis[[2]]) * lambda[, 2] / (
                    t(expr.dis[[2]]) * lambda[, 2] + 
                    t(expr.dis[[1]]) * lambda[, 1] + 
                    t(expr.dis[[3]]) * lambda[, 3])
    
    weight.sample = colMeans(weight.sample)
  }
  
  p.h.status = TRUE
  if (missing(p.h)) {
    p.h = setNames(vector("list", length(gene)), gene)
    p.h.status = FALSE
  }
  
  weight.status = TRUE
  if (missing(weight)) {
    weight = p.h
    weight.status = FALSE
  }
  
  ## =============
  for (gene.i in gene) {
    if (debug == TRUE) {
      cat(gene.i, "\n")
    }
    
    id = mut[, "hgnc_symbol"] %in% gene.i
    sample.mut.gene.i= unique(mut[id, "sample"])
    if (length(sample.mut.gene.i) == 0) { # no mutations
      next
    }
    
    p.df[[gene.i]] = ReplicateData(cpd$p.df.share, 
                                   temporal=TRUE,
                                   length(sample.mut.gene.i), 
                                   sample.mut.gene.i)
    ## For p.fg and p.gy 
    if (cis == TRUE) {
      inf.gene.name = gene.i
    } else {
      x = net[[gene.i]]
      inf.gene.name = names(x)
      inf.gene.name = setdiff(inf.gene.name, gene.i)
    }
    
    if (!missing(expr.dis)) {
      inf.gene.name = intersect(inf.gene.name, gene.expr) 
    }
    
    if (length(inf.gene.name) > 0) {
      num.neg.reg = t(expr[sample.mut.gene.i, inf.gene.name, drop=FALSE]) < 
        expr.dis$mu[inf.gene.name, 2]
      neg.reg.ratio = rowSums(num.neg.reg) / length(sample.mut.gene.i) 
    }
    
    ## Add directions
    if (p.h.status == TRUE) {
      inf.gene.name = intersect(names(p.h[[gene.i]]), inf.gene.name)
      p.h[[gene.i]] = p.h[[gene.i]][inf.gene.name]
    }
    
    if (weight.status == TRUE) {
      weight[[gene.i]] = weight[[gene.i]][inf.gene.name]
    } else {
      weight[[gene.i]] = setNames(seq(length(inf.gene.name)), inf.gene.name)
    }
    
    if (cis == TRUE) {
      weight[[gene.i]] = setNames(1, gene.i)
      p.h[[gene.i]] = setNames(1-.Machine$double.eps, gene.i) 
    } else {
      if (p.h.status==TRUE & weight.status==FALSE) {
        weight = sapply(p.h, function(z) max(z, 1-z))
      } else if (p.h.status==FALSE & !missing(expr.dis)) {
        if (length(sample.mut.gene.i) <= 3) {
          p.h[[gene.i]] = setNames(rep(0.5, length(inf.gene.name)), 
                                   inf.gene.name)
          
          if (weight.status == FALSE) {
            weight[[gene.i]] = p.h[[gene.i]]
          }
        }
        
        if (length(sample.mut.gene.i)>3 & length(inf.gene.name)>0) {
          down.share = expr.dis[[1]][sample.mut.gene.i, 
                                     inf.gene.name, drop=FALSE]     
          up.share   = expr.dis[[3]][sample.mut.gene.i, 
                                     inf.gene.name, drop=FALSE]
          
          down.weight   = down.share / (down.share + up.share)
          down.weight.median = apply(down.weight, 2, median)
          
          id = (down.weight.median > 0.5 & neg.reg.ratio < 0.5) | 
               (down.weight.median < 0.5 & neg.reg.ratio > 0.5)
          
          if (weight.status == FALSE) {
            weight[[gene.i]] = 
              pmax(down.weight.median, 1-down.weight.median)
            weight[[gene.i]] = weight[[gene.i]] * net[[gene.i]][inf.gene.name]
          }
          p.h[[gene.i]] = down.weight.median
          p.h[[gene.i]][id] = 0.5
        }
      }
    }
    
    p.fg[[gene.i]] = vector("list", length(sample.mut.gene.i)) 
    names(p.fg[[gene.i]]) = sample.mut.gene.i
    
    for (sample.i in sample.mut.gene.i) {      
      p.fg[[gene.i]][[sample.i]] = 
        ReplicateData(cpd$p.fg.share, temporal=TRUE, 
                      length(inf.gene.name), inf.gene.name)
      
      if (!missing(expr.dis))
        if (sample.i %in% sample.expr) {
          weight.gene = weight[[gene.i]]
          
          p.gy.down = expr.dis[[1]][sample.i, inf.gene.name, drop=FALSE]
          p.gy.down = t(t(p.gy.down) ^ weight.gene)
          p.gy.down = p.gy.down ^ weight.sample[sample.i]
          
          p.gy.neutral = expr.dis[[2]][sample.i, inf.gene.name, drop=FALSE]
          p.gy.neutral = t(t(p.gy.neutral) ^ weight.gene)
          p.gy.neutral = p.gy.neutral ^ weight.sample[sample.i]
          
          p.gy.up = expr.dis[[3]][sample.i, inf.gene.name, drop=FALSE]
          p.gy.up = t(t(p.gy.up) ^ weight.gene)
          p.gy.up = p.gy.up ^ weight.sample[sample.i]
          
          tmp = t(rbind(p.gy.down, p.gy.neutral, p.gy.up))
          colnames(tmp) = c(-1, 0, 1)          
          p.gy[[gene.i]][[sample.i]] = tmp
        }
    }
  }
  
  model = list(p.d=p.d, p.df=p.df, p.fg=p.fg, p.gy=p.gy, 
               p.h=p.h, weight=weight,
               
               beta.pd.share=prior$beta.pd.share, 
               beta.df.share=prior$beta.df.share, 
               dir.fg.share =prior$dir.fg.share, 
               
               p.d.share =cpd$p.d.share, 
               p.df.share=cpd$p.df.share, 
               p.fg.share=cpd$p.fg.share)
  return(model)
} 


#==============================================================================
## Potential for f
SummarizeVar = function(model, gene, direction=TRUE) {
  if(missing(gene)) {
    gene = names(model$p.d)
  } else {
    gene = intersect(gene, names(model$p.d))
  }
  
  ## Summerize G, note NA
  f.potential = vector("list", length=length(gene))
  names(f.potential) = gene
  h.potential = f.potential
  for (gene.i in gene) {
    sample.mut.gene.i = names(model$p.gy[[gene.i]])
    f.potential.i = matrix(0, nrow=length(sample.mut.gene.i), ncol=2, 
                           dimnames=list(sample.mut.gene.i, c("0", "1")))
    if (direction == TRUE) {
      gene.inf = rownames(model$p.gy[[gene.i]][[1]])
      h.potential.i = matrix(0, nrow=length(gene.inf), ncol=2, 
                             dimnames=list(gene.inf, c("0", "1")))
    }
    
    for (sample.i in sample.mut.gene.i) {
      if(nrow(model$p.gy[[gene.i]][[sample.i]]) == 0) {
        next
      }
      
      ## Summarize Y
      p.fg = model$p.fg[[gene.i]][[sample.i]]
      if (direction == TRUE) {
        p.h = model$p.h[[gene.i]]
        p.h = list(1, p.h*1.0, (1-p.h)*1.0)
        
        p.fg = lapply(seq(3), function(z) p.fg[[z]] * p.h[[z]])
      }
      potential.f.sample.up = sapply(p.fg,  function(x) 
        x * model$p.gy[[gene.i]][[sample.i]], simplify=FALSE)
      
      ## Summarize G
      tmp.gene.i = sapply(potential.f.sample.up, function(x) rowSums(x))
      if (is.matrix(tmp.gene.i)) {
        if (direction == TRUE) {
          h.potential.i = h.potential.i + tmp.gene.i[, 2:3]
          
          tmp.gene.i[, 2] = tmp.gene.i[, 2] + tmp.gene.i[, 3]
          tmp.gene.i = tmp.gene.i[, 1:2, drop=FALSE]
        } 
        ## Mulitiply for different connected genes
        tmp.gene.i = log(tmp.gene.i)
        tmp.gene.i = colSums(tmp.gene.i)
      } else {
        if (direction == TRUE) {
          h.potential.i = h.potential.i + tmp.gene.i[2:3]
          
          tmp.gene.i[2] = tmp.gene.i[2] + tmp.gene.i[3]
          tmp.gene.i = tmp.gene.i[1:2]
        } 
        tmp.gene.i = log(tmp.gene.i)
      }
      
      f.potential.i[sample.i, ] = tmp.gene.i
    }
    
    f.potential[[gene.i]] = f.potential.i 
    
    if (direction == TRUE) {
      h.potential[[gene.i]] = h.potential.i 
    }
  }
  
  return(list(f.potential=f.potential, h.potential=h.potential))
}

#==============================================================================
#' @importFrom stats setNames
CompPotentialH = function(model, gene, potential.f) {
  
  gene.inf = rownames(model$p.gy[[gene]][[1]])
  sample.mut.gene.i = names(model$p.gy[[gene]])
  h.potential.sample = setNames(vector("list", length(sample.mut.gene.i)), 
                                sample.mut.gene.i)
  
  for (sample.i in sample.mut.gene.i) {
    if(nrow(model$p.gy[[gene]][[sample.i]]) == 0) {
      next
    }
    
    p.fg = model$p.fg[[gene]][[sample.i]]
    
    p.h = model$p.h[[gene]]
    p.h = list(p.h*1.0, (1-p.h)*1.0)
    p.fg = lapply(seq(2), function(z) p.fg[[z+1]] * p.h[[z]])
    
    potential.f.sample.up = 
      sapply(p.fg, function(x) 
        x * model$p.gy[[gene]][[sample.i]], simplify=FALSE)
    
    tmp.gene.i = sapply(potential.f.sample.up, function(x) rowSums(x))
    tmp = log(tmp.gene.i) + potential.f[sample.i, 2]
    h.potential.sample[[sample.i]] = tmp
    
  }
  
  return(h.potential.sample)
}

#==============================================================================
## Potential for D
IntegratePrior = function(model, gene=NULL) {
  if(missing(gene)) {
    gene = names(model$beta.pd)
  } else {
    gene = intersect(gene, names(model$beta.pd))
  }
  
  pi.integral = vector("list", length(gene))
  names(pi.integral) = gene
  for (gene.i in gene) {
    gamma.value = gamma(c(model$beta.pd[[gene.i]], 
                          sum(model$beta.pd[[gene.i]])))
    pi.integral.gene.i = gamma.value[3] / gamma.value[1] / 
      gamma.value[2] / gamma(1+sum(model$beta.pd[[gene.i]]))
    
    pi.integral.1 = pi.integral.gene.i * gamma(1+model$beta.pd[[gene.i]][2]) *
      gamma(model$beta.pd[[gene.i]][1])
    pi.integral.0 = pi.integral.gene.i * gamma(1+model$beta.pd[[gene.i]][1]) * 
      gamma(model$beta.pd[[gene.i]][2])
    
    pi.integral.gene.i = t(c(pi.integral.0, pi.integral.1))
    colnames(pi.integral.gene.i) = names(model$beta.pd[[gene.i]])
    pi.integral[[gene.i]] = log(pi.integral.gene.i)
  }
  return(pi.integral)
}

#==============================================================================
#' Learn xseq parameters given an initialized model
#' 
#' @export
#' @param model An xseq model.
#' @param constraint A list of constraints on \eqn{\theta_{G|F}}. 
#' @param debug Logical, specifying whether debug information should be printed
#' @return The posterior probabilities of latent variables in xseq 
#' 
#' @importFrom stats setNames
#' 
InferXseqPosterior = function(model, constraint, debug=FALSE) {
  direction = TRUE
  if (nrow(model$p.fg.share) == 2) {
    direction = FALSE
  }
  
  f.summarize.g = SummarizeVar(model, direction=direction)
  f.summarize.g = f.summarize.g$f.potential
  
  gene = names(model$p.d)
  if(length(gene) == 0) {
    return
  }
  
  nstates.d = length(model$p.df[[1]])
  nstates.f = max(sapply(model$p.df, function(x) ncol(x[[1]])))
  nstates.g = sort(do.call(rbind, lapply(model$p.fg, function(x) 
    ncol(x[[1]][[1]])))[,1])[1]
  
  posterior.d = setNames(vector("list", length(gene)), gene)
  posterior.h = setNames(vector("list", length(gene)), gene)
  posterior.f = setNames(vector("list", length(gene)), gene)
  posterior.g = setNames(vector("list", length(gene)), gene)
  
  sufficient.df = setNames(vector("list", length(gene)), gene)
  sufficient.fg = setNames(vector("list", length(gene)), gene)
  
  loglik = sum(log(model$p.d.share) * (model$beta.pd.share - 1))   +
    sum(log(model$p.df.share[1,]) * (model$beta.df.share[1,] - 1)) +
    sum(log(model$p.df.share[2,]) * (model$beta.df.share[2,] - 1))
  
  loglik = loglik + 
    sum(log(model$p.fg.share[1,]) * (model$dir.fg.share[1,] - 1)) + 
    sum(log(model$p.fg.share[2,]) * (model$dir.fg.share[2,] - 1))
  
  if (direction == TRUE) {
    loglik = loglik + 
      sum(log(model$p.fg.share[3,]) * (model$dir.fg.share[3,] - 1))
  }
  
  loglik.gene = setNames(rep(0, length(gene)), gene)
  for (gene.i in gene) {
    no.connection = nrow(model$p.gy[[gene.i]][[1]]) == 0
    d.f.summarize.g = 
      sapply(model$p.df[[gene.i]], function(x) 
        log(x) + f.summarize.g[[gene.i]], simplify=FALSE)
    
    d.summarize.f = sapply(d.f.summarize.g, function(x) 
      apply(x, 1, LogSumExp), simplify=TRUE)
    
    ## Mulifply potential D for different sample
    if(is.matrix(d.summarize.f)) {
      d.summarize = colSums(d.summarize.f) 
    } else {
      d.summarize = d.summarize.f
      d.summarize.f = t(d.summarize.f)
    }
    
    posterior.d[[gene.i]]  = d.summarize + log(model$p.d[[gene.i]])    
    loglik.gene[gene.i] = LogSumExp(posterior.d[[gene.i]])
    loglik = loglik + loglik.gene[gene.i]
    posterior.d[[gene.i]]  = exp(posterior.d[[gene.i]] - 
                                   LogSumExp(posterior.d[[gene.i]]))
    
    sample.gene.i = rownames(f.summarize.g[[gene.i]])
    potential.d = sapply(seq(length(sample.gene.i)), 
                         function(id) d.summarize - d.summarize.f[id, ])
    potential.d = t(potential.d)
    rownames(potential.d) = sample.gene.i
    
    potential.d.prior = t(apply(potential.d, 1, 
                                function(x) x + log(model$p.d[[gene.i]])))
    
    d.f.summarize.g.suf = sapply(seq(2), function(id) 
      apply(d.f.summarize.g[[id]], 2, function(z) 
        z + potential.d.prior[, id]), simplify=FALSE)
    if (!is.matrix(d.f.summarize.g.suf[[1]])) {
      d.f.summarize.g.suf = sapply(d.f.summarize.g.suf, t, simplify=FALSE)
    }
    sufficient.df[[gene.i]] = d.f.summarize.g.suf
    
    
    potential.d.f.down = sapply(seq(nstates.d), function(id) 
      log(model$p.df[[gene.i]][[id]]) + potential.d.prior[, id], 
      simplify=FALSE)
    
    potential.f.down = sapply(seq(nstates.f), function(id) 
      LogPlus(potential.d.f.down[[1]][, id], potential.d.f.down[[2]][, id]))
    
    if (!is.matrix(potential.f.down)) {
      potential.f.down = t(potential.f.down)
      rownames(potential.f.down) = rownames(potential.d.f.down[[1]])
      colnames(potential.f.down) = colnames(potential.d.f.down[[1]])
    }
    
    potential.f = potential.f.down + f.summarize.g[[gene.i]]
    posterior.f[[gene.i]] = exp(potential.f - 
                                  apply(potential.f, 1, LogSumExp))    
    
    if (direction == TRUE) {
      potential.h.sample = CompPotentialH(model, gene=gene.i, potential.f.down)
      potential.h = Reduce("+", potential.h.sample)
      
      if (is.matrix(potential.h)) {
        sum.potential = LogPlus(potential.h[,1], potential.h[,2])
      } else {
        sum.potential = LogPlus(potential.h[1], potential.h[2])
      }
      
      posterior.h[[gene.i]] = exp(potential.h - sum.potential)
    }
    
    
    if (no.connection == TRUE) {
      next
    }
    gene.connection = rownames(model$p.gy[[gene.i]][[1]])
    
    sufficient.sample.i = setNames(vector("list", 
                                          length(sample.gene.i)), sample.gene.i)
    posterior.sample.i = setNames(vector("list", 
                                         length(sample.gene.i)), sample.gene.i)
    for (sample.i in sample.gene.i) {      
      p.fg = model$p.fg[[gene.i]][[sample.i]]
      if (direction == TRUE) {
        p.h = model$p.h[[gene.i]] 
        p.h = list(1, p.h*1.0, (1-p.h)*1.0)
        
        p.fg = lapply(seq(3), function(z) p.fg[[z]] * p.h[[z]])
      }
      potential.f.sample.up = 
        sapply(p.fg, function(x) 
          x * model$p.gy[[gene.i]][[sample.i]], simplify=FALSE)
      
      potential.f.sample.up.simplify = potential.f.sample.up
      if (direction == TRUE) {
        potential.f.sample.up.simplify[[2]] = 
          potential.f.sample.up.simplify[[2]] +
          potential.f.sample.up.simplify[[3]]
        potential.f.sample.up.simplify[[3]] = NULL
      }
      potential.f.sample.up.sum = log(sapply(potential.f.sample.up.simplify, 
                                             rowSums, simplify=TRUE))
      if(!is.matrix(potential.f.sample.up.sum)) {
        potential.f.sample.up.sum = t(potential.f.sample.up.sum)
      }
      potential.f.sample.up.sum = 
        t(sapply(seq(length(gene.connection)), function(id) 
          f.summarize.g[[gene.i]][sample.i, ] - potential.f.sample.up.sum[id, ]))
      
      potential.f.sample = t(apply(potential.f.sample.up.sum, 1, function(x) 
        x + potential.f.down[sample.i, ]))
      
      potential.f.g.sample = lapply(seq(2), function(id) 
        log(potential.f.sample.up.simplify[[id]]) + potential.f.sample[, id])
      if (direction == TRUE) {
        potential.f.g.sample.a = lapply(seq(2), function(id) 
          log(potential.f.sample.up[[id]]) + potential.f.sample[, id])
        
        potential.f.g.sample.a[[3]] = log(potential.f.sample.up[[3]]) + 
          potential.f.sample[, 2]
        
        sufficient.sample.i[[sample.i]] = potential.f.g.sample.a
      } else {
        sufficient.sample.i[[sample.i]] = potential.f.g.sample
      }
      
      potential.g.sample = sapply(seq(3), function(id) 
        LogPlus(potential.f.g.sample[[1]][, id], 
                potential.f.g.sample[[2]][, id]))
      
      ## Matrix
      if (!is.matrix(potential.g.sample)) {
        potential.g.sample = t(potential.g.sample)
      }
      posterior.sample.i[[sample.i]] = 
        t(apply(potential.g.sample, 1, function(x) exp(x - LogSumExp(x))))
    }
    sufficient.fg[[gene.i]] = sufficient.sample.i
    posterior.g[[gene.i]] = posterior.sample.i
  }    
  
  return(list(posterior.d = posterior.d, 
              posterior.h = posterior.h,
              posterior.f = posterior.f, 
              posterior.g = posterior.g, 
              sufficient.fg = sufficient.fg,
              sufficient.df = sufficient.df, 
              loglik = loglik, 
              loglik.gene = loglik.gene))
}


#==============================================================================
#' @importFrom stats setNames
Maximize = function(model, sufficient.statistics, constraint, 
                    cis=FALSE, debug=FALSE) {
  direction = TRUE
  if (nrow(model$p.fg.share) == 2) {
    direction = FALSE
  }
  
  gene.all = names(model$p.d)
  
  gene.connection = sapply(model$p.fg, function(x) length(x[[1]][[1]]))
  gene = names(which(gene.connection > 0))
  
  states.d = names(model$p.df[[1]])
  states.f = names(model$p.fg[[1]][[1]])
  
  num.mut  = sapply(sufficient.statistics$sufficient.df, 
                    function(x) nrow(x[[1]]))
  num.conn = sapply(sufficient.statistics$sufficient.fg, 
                    function(x) nrow(x[[1]][[1]]))
  if (is.list(num.conn)) {
    num.conn = do.call(rbind, num.conn)[,1]  
  }
  num.mut  = num.mut[gene]
  num.conn = num.conn[gene]
  
  ## Based on sufficient statistics
  reduce.d  = Reduce("+", sufficient.statistics$posterior.d)
  reduce.d  = reduce.d + (model$beta.pd.share - 1)
  p.d.share = reduce.d / sum(reduce.d)
  
  p.df.all = setNames(vector("list", length(gene)), gene)
  p.fg.all = setNames(vector("list", length(gene)), gene)
  for (gene.i in gene) {
    p.df = sapply(sufficient.statistics$sufficient.df[[gene.i]], function(x) 
      apply(x, 2, LogSumExp), simplify=TRUE)
    p.df.all[[gene.i]] = t(p.df) - sufficient.statistics$loglik.gene[gene.i]
    
    sample = names(sufficient.statistics$sufficient.fg[[gene.i]])
    p.fg.gene = setNames(vector("list", length(sample)), sample)
    for (sample.i in sample) {
      p.fg.gene[[sample.i]] = 
        sapply(sufficient.statistics$sufficient.fg[[gene.i]][[sample.i]], 
               function(x) apply(x, 2, LogSumExp), simplify=TRUE) 
    }
    p.fg.gene = Reduce("LogPlus", p.fg.gene)
    p.fg.gene = t(p.fg.gene) - sufficient.statistics$loglik.gene[gene.i]
    p.fg.all[[gene.i]] = p.fg.gene
  }
  
  if(length(gene) > 0) {  
    p.df.share = Reduce("LogPlus", p.df.all)  
    p.df.share = log(exp(p.df.share) + (model$beta.df.share - 1))
    
    p.df.share = t(apply(p.df.share, 1, function(x) exp(x - LogSumExp(x))))
    if (p.df.share[1, 1] <= p.df.share[1, 2]+0.1) {
      p.df.share[1, ] = model$p.df.share[1, ]
    } 
    if (p.df.share[2, 2] <= p.df.share[2, 1]+0.1) {
      p.df.share[2, ] = model$p.df.share[2, ]
    }
    
    if (debug == TRUE) {
      cat("p.df.share: \n")
      print(p.df.share)
    }
    
    p.fg.share = Reduce("LogPlus", p.fg.all)
    p.fg.share = log(exp(p.fg.share) + (model$dir.fg.share - 1)) 
    
    if (constraint$equal.fg == TRUE) {
      if (direction == TRUE) {
        p.fg.share[2, 1] = LogPlus(p.fg.share[2, 1], p.fg.share[3, 3])
        p.fg.share[2, 2] = LogPlus(p.fg.share[2, 2], p.fg.share[3, 2])
        p.fg.share[2, 3] = LogPlus(p.fg.share[2, 3], p.fg.share[3, 1])
        
        p.fg.share[3, ]  = rev(p.fg.share[2, ])
      } else {
        p.fg.share[2, 1] = LogPlus(p.fg.share[2, 1], 
                                   p.fg.share[2, 3]) - log(2)
        p.fg.share[2, 3] = p.fg.share[2, 1]
      }
    }
    
    p.fg.share = t(apply(p.fg.share, 1, function(x) exp(x - LogSumExp(x))))
    if (!is.null(constraint$baseline)) {
      if (p.fg.share[1, 2] <= constraint$baseline) {
        p.fg.share[1, 2] = model$p.fg.share[1, 2]
      }
      if (p.fg.share[2, 2] >= constraint$baseline) {
        p.fg.share[2, 2] = model$p.fg.share[2, 2]
      }
      
      if (direction == TRUE) {
        if (p.fg.share[3, 2] >= constraint$baseline) {
          p.fg.share[3, 2] = model$p.fg.share[3, 2]
        }
      }
      
      low.threshold = 0.4
      if (TRUE == cis) {
        low.threshold = 0.2
      }
      
      if (p.fg.share[2, 2] <= low.threshold) {
        p.fg.share[2, 2] = model$p.fg.share[2, 2]
      }
      if (direction == TRUE) {
        if (p.fg.share[3, 2] <= low.threshold) {
          p.fg.share[3, 2] = model$p.fg.share[3, 2]
        }
      }
      
      p.fg.share[1, c(1,3)] = p.fg.share[1, c(1,3)] / 
        sum(p.fg.share[1, c(1,3)]) * (1 - p.fg.share[1,2])
      
      p.fg.share[2, c(1,3)] = p.fg.share[2, c(1,3)] / 
        sum(p.fg.share[2, c(1,3)]) * (1 - p.fg.share[2,2])
      
      if (direction == TRUE) {
        p.fg.share[3, c(1,3)] = p.fg.share[3, c(1,3)] / 
          sum(p.fg.share[3, c(1,3)]) * (1 - p.fg.share[3,2])
      }
    }
    
    if (debug == TRUE) {
      cat("p.fg.share: \n")
      print(p.fg.share) 
    }
  }
  
  model.new = model
  model.new$p.d = sapply(model.new$p.d, function(x) 
    x=p.d.share, simplify=FALSE)
  for (gene.i in gene) {
    sample.mut.gene.i = rownames(model$p.df[[gene.i]][[2]])
    model.new$p.df[[gene.i]] = 
      ReplicateData(p.df.share, temporal=TRUE, 
                    length(sample.mut.gene.i), sample.mut.gene.i)
    
    sample = names(sufficient.statistics$sufficient.fg[[gene.i]])
    inf.gene.name = rownames(model$p.fg[[gene.i]][[1]][[1]])
    for (sample.i in sample) {
      model.new$p.fg[[gene.i]][[sample.i]] = 
        ReplicateData(p.fg.share, temporal=TRUE, 
                      length(inf.gene.name), inf.gene.name)
    }
  }
  
  model.new$p.d.share  = p.d.share
  model.new$p.df.share = p.df.share
  model.new$p.fg.share = p.fg.share
  
  return(model.new)
}

#==============================================================================
## Based on log-likelihood 
CheckConvergence = function(loglik, loglik.new, threshold=1e-5) {  
  loglik.diff = (loglik - loglik.new) / loglik
  
  ## Note numeric issues
  done = FALSE
  if (abs(loglik.diff) < threshold) {
    done = TRUE
  } else if (loglik.diff < 0) {
    warning("WARNING! The likelihood decreases", "\n")
    done = TRUE
  }
  
  return(done)
}

#==============================================================================
#' Learn xseq parameters given an initialized model
#' 
#' @export
#' @param model An xseq model
#' @param constraint A list of constraints on \eqn{\theta_{G|F}}. 
#' @param iter.max Maximum number of iterations in learning xseq parameters
#' @param threshold The threshold to stop learning paramters
#' @param cis Logical, cis-analysis or trans-analysis
#' @param debug Logical, specifying whether debug information should be printed
#' @param show.plot Logical, specifying whether to plot the Log-Likelihoods
#' @return A list including the learned xseq model
LearnXseqParameter = function(model, constraint, iter.max=20, threshold=1e-5, 
                              cis=FALSE, debug=FALSE,  show.plot=TRUE) {
  
  direction = TRUE
  if (nrow(model$p.fg.share) == 2) {
    direction = FALSE
  }
  
  iter = 1
  loglik = rep(0, length(iter.max))
  
  if (constraint$equal.fg == TRUE) {
    p.fg.share = model$p.fg.share
    p.fg.share[1, 1] = p.fg.share[1, 3] = (1 - p.fg.share[1, 2]) / 2
    
    if (direction == FALSE) {
      p.fg.share[2, 1] = p.fg.share[2, 3] = (1 - p.fg.share[2, 2]) / 2
    }
    
    model$p.fg.share = p.fg.share
    
    model$p.fg = lapply(model$p.fg, function(x) {
      lapply(x, function(z) {
        if(length(z[[1]]) == 0) {
          return(z)
        }
        
        z[[1]][1, 1] = z[[1]][1, 3] = p.fg.share[1, 1]
        if (direction == FALSE) {
          z[[2]][, 1] = z[[2]][, 3] = p.fg.share[2, 1]
        }
        
        return(z)
      })
    })
  }
  
  sufficient.statistics = 
    InferXseqPosterior(model, constraint=constraint, debug=debug)
  model.new = Maximize(model, sufficient.statistics, 
                       constraint=constraint, debug=debug)
  loglik[iter] = sufficient.statistics$loglik
  
  p.d.trace  = vector("list", iter.max)
  p.df.trace = vector("list", iter.max)
  p.fg.trace = vector("list", iter.max)
  
  p.d.trace[[iter]]  = model$p.d.share
  p.df.trace[[iter]] = model$p.df.share
  p.fg.trace[[iter]] = model$p.fg.share
  
  done = FALSE
  while(iter < iter.max & !done) {
    iter = iter + 1
    
    p.d.trace[[iter]]  = model.new$p.d.share
    p.df.trace[[iter]] = model.new$p.df.share
    p.fg.trace[[iter]] = model.new$p.fg.share
    
    model = model.new
    sufficient.statistics = 
      InferXseqPosterior(model, constraint=constraint, debug=debug)
    model.new = Maximize(model, sufficient.statistics, cis=cis, 
                         constraint=constraint, debug=debug)
    
    loglik[iter] = sufficient.statistics$loglik
    done = CheckConvergence(loglik[iter-1], loglik[iter], threshold)
  }
  p.d.trace[[iter+1]]  = model.new$p.d.share
  p.df.trace[[iter+1]] = model.new$p.df.share
  p.fg.trace[[iter+1]] = model.new$p.fg.share
  
  model = model.new
  
  if (show.plot == TRUE) {
    NeatPlot(loglik, type="o", col="dodgerblue", lwd=2.5, 
         ylab="Log-likelihood", xlab="Iteration")
  }
  
  para.trace = list(p.d.trace=p.d.trace, p.df.trace=p.df.trace, 
                    p.fg.trace=p.fg.trace)
  return(list(model=model, posterior=sufficient.statistics, 
              loglik=loglik, para.trace=para.trace))
}


# ==============================================================================
#' Convert xseq output to a data.frame
#' 
#' @export
#' @param posterior The posterior probabilities of mutations and mutated genes, 
#' output from \code{InferXseqPosterior}
#' @return A data.frame with sample, gene, 
#' probability of individual mutations and the probabilities of individual mutated genes
#' 
ConvertXseqOutput = function(posterior) {
  gene.all = names(posterior$posterior.d)
  count    = sapply(posterior$posterior.f, nrow)
  gene.all = rep(gene.all, times=count)
  
  prob.f   = do.call(rbind, posterior$posterior.f)
  mut.prob = data.frame(sample=rownames(prob.f), hgnc_symbol=gene.all,
                        driver_prob=prob.f[,2], stringsAsFactors=FALSE)
  
  driver.prob.gene = do.call(rbind, posterior$posterior.d)[, 2]
  driver_prob_gene = rep(driver.prob.gene, times=count)
  mut.prob = cbind(mut.prob, driver_prob_gene)
  mut.prob = mut.prob[order(mut.prob[, "driver_prob_gene"], decreasing=TRUE), ]
  colnames(mut.prob) = c("sample", "hgnc_symbol", "P(F)", "P(D)")
  
  return(mut.prob)
}


#==============================================================================
SpecifyPriorModel = function(p.d, p.df, p.fg) {
  if (missing(p.d)) {
    p.d = structure(c(0.95, 0.05), names=c(0, 1))
  }
  
  if (missing(p.df)) {
    p.df = matrix(c(0.90, 0.10, 0.10, 0.90), 2, 2, 
                  byrow=TRUE, dimnames=list(seq(2)-1, seq(2)-1))
  }
  
  if (missing(p.fg)) {
     p.fg = matrix(c(0.025, 0.95, 0.025, 0.20, 0.60, 0.20), 
                  2, 3, byrow=TRUE, dimnames=list(seq(2)-1, seq(3)-1))
  }
  
  cpd = list(p.d.share=p.d, p.df.share=p.df, p.fg.share=p.fg)
  return (cpd)
}


#==============================================================================
#' Set model paramerter priors
#' 
#' @export
#' @param mut A data.frame of mutations. The data.frame should have 
#' three columns of characters: sample, hgnc_symbol, and variant_type. 
#' The variant_type column cat be either "HOMD", "HLAMP", 
#'   "MISSENSE", "NONSENSE", "FRAMESHIFT", "INFRAME", "SPLICE", "NONSTOP", 
#'   "STARTGAINED", "SYNONYMOUS", "OTHER", "FUSION", "COMPLEX".
#' @param expr.dis A list, the outputs from calling \code{GetExpressionDistribution}
#' @param regulation.direction Logical, whether considering the directionality
#' , i.e., up-regulation or down-regulatin of genes, only used
#' when \code{cis=FALSE}.
#' @param net A list of gene interactions
#' @param cis Logical, cis analysis or trans analysis
#' @param mut.type Character, only used when \code{cis = TRUE}, and can be either
#' loss, gain or both
#' @param ... Reserved for extension
#' 
#' @importFrom stats setNames median sd
#' 
SetXseqPrior = function(mut, expr.dis, net, regulation.direction=TRUE, 
                    cis=TRUE, mut.type="loss", ...) {
  
  cpd = SpecifyPriorModel(...)
  virtual.data.fraction = 0.1  # effective virtual data size is 10%
  
  ## genes with both mutation, expression data and network connections
  gene  = intersect(unique(mut[, "hgnc_symbol"]), colnames(expr.dis[[1]]))
  
  if (cis == FALSE) {
    count.conn = sapply(net, length)
    gene = intersect(names(count.conn), gene)
    id  = mut[, "hgnc_symbol"] %in% gene
    mut = mut[id, ]
  } else {
    count.conn = setNames(seq(length(gene)), gene)
  }
  
  id  = paste(mut[, "sample"], mut[, "hgnc_symbol"], sep="_")
  id  = !duplicated(id)
  mut = mut[id, ]
  
  count.mut  = table(mut[, "hgnc_symbol"])
  
  ## Typically we expect less than 200 high probability genes, 
  beta.pd.share = ceiling(cpd$p.d.share * length(gene) * 
                          virtual.data.fraction)
  if (beta.pd.share[2] == 1) {
    beta.pd.share = beta.pd.share + c(1, 1)
  } else if (beta.pd.share[2] > 200) {
    tmp = beta.pd.share
    tmp[1] = sum(beta.pd.share) - 200
    tmp[2] = 200
    
    cpd$p.d.share = tmp / sum(tmp)
    beta.pd.share = tmp
  }
  
  median.mut = median(count.mut)
  beta.df.share = cpd$p.df.share * cpd$p.d.share * 
    length(gene) * median.mut * virtual.data.fraction
  beta.df.share = ceiling(beta.df.share)
  
  if(sum(beta.df.share==1) > 0) {
    beta.df.share = beta.df.share + 1
  }
  
  expr.dis.lambda = expr.dis[[4]]
  lambda = colMeans(expr.dis.lambda)
  lambda.sd = apply(expr.dis.lambda, 2, sd)
  
  p.fg.plus  = lambda + lambda.sd * c(0, 1, 0)
  p.fg.plus  = p.fg.plus / sum(p.fg.plus)
  
  if (regulation.direction == TRUE | 
      (mut.type %in% c("gain", "loss") & cis == TRUE)) {
    p.fg.h.plus  = p.fg.plus
    p.fg.h.plus[3] = p.fg.h.plus[2] = (1 - p.fg.h.plus[1]) / 2
    
    p.fg.h.minus = p.fg.plus
    p.fg.h.minus[1] = p.fg.h.minus[2] = (1 - p.fg.h.minus[3]) / 2
  } else {
    p.fg.minus = lambda + 5 * lambda.sd * c(1, 0, 1)
    p.fg.minus = p.fg.minus / sum(p.fg.minus)
    p.fg.share = rbind(p.fg.plus, p.fg.minus)
  }

  if (cis == TRUE) {
    if (mut.type == "loss") {
      p.fg.share = rbind(p.fg.plus, p.fg.h.minus)
    } else if (mut.type == "gain") {
      p.fg.share = rbind(p.fg.plus, p.fg.h.plus)
    }
  }
  
  if (regulation.direction == FALSE | cis == TRUE) {
    dir.fg.share = ceiling(cpd$p.df.share %*% 
                           p.fg.share * cpd$p.d.share * 
                           length(gene) * median.mut * median(count.conn) * 
                           virtual.data.fraction)
  } else {
    p.fg.share = rbind(p.fg.plus, p.fg.h.minus, p.fg.h.plus)
    dir.fg.minus = cpd$p.df.share %*% 
                   p.fg.share[1:2, ] * cpd$p.d.share * 
                   length(gene) * median.mut * median(count.conn) * 
                   virtual.data.fraction
    dir.fg.plus  = cpd$p.df.share %*% 
                   p.fg.share[c(1,3), ] * cpd$p.d.share * 
                   length(gene) * median.mut * median(count.conn) * 
                   virtual.data.fraction
    
    dir.fg.share = rbind((dir.fg.minus[1,] + dir.fg.plus[1,])/2, 
                         dir.fg.minus[2,], 
                         dir.fg.plus[2,])
    dir.fg.share = ceiling(dir.fg.share)
  }
  cpd$p.fg.share = p.fg.share
  
  if(sum(dir.fg.share==1) > 0) {
    dir.fg.share = dir.fg.share + 1
  }
  
  prior = setNames(vector("list", 3), 
                   c("beta.pd.share", "beta.df.share", "dir.fg.share"))
  prior$beta.pd.share = beta.pd.share
  prior$beta.df.share = beta.df.share
  prior$dir.fg.share  = dir.fg.share
  
  return(list(prior=prior, cpd=cpd, baseline=lambda[2]))
}

