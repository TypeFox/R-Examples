#==============================================================================
#' Get the conditional distributions for a set of genes
#' 
#' @export
#' @param expr A matrix of gene expression values where 
#' each row corresponds to a patient and each column is a gene.
#' @param mut A data.frame of mutations. The data.frame should have 
#' three columns of characters: sample, hgnc_symbol, and variant_type. 
#' The variant_type column cat be either "HOMD", "HLAMP", 
#'   "MISSENSE", "NONSENSE", "FRAMESHIFT", "INFRAME", "SPLICE", "NONSTOP", 
#'   "STARTGAINED", "SYNONYMOUS", "OTHER", "FUSION", "COMPLEX".
#' @param cna.call A matrix containing the copy number calls, where each element is coded:
#' \itemize{
#'   \item -2, homozygous deletions
#'   \item -1, hemizygous deletions
#'   \item 0, neutral
#'   \item 1, gain
#'   \item 2, amplifications
#'   }
#' @param gene A character vector of official HGNC gene names
#' @param type Character, either Gaussian  ("gauss") or 
#' Student ("student")
#' @param show.plot Logical, specifying whether to plot the fitted model
#' @return  A list containg the fitted expression distributions
#' 
#' @importFrom stats setNames dnorm
#' @importFrom graphics split.screen close.screen screen par
#' 
GetExpressionDistribution = function(expr, mut=NULL, cna.call=NULL, 
                     gene=NULL, type="student", show.plot=FALSE) {
  # Get the conditional distributions for a set of genes
  #
  # Args: 
  #   mut: a data.frame, the sample-mutation list
  #   expr: a matrix, the sample-expression matrix
  #
  # Note currently we use name matching for simplicity, i.e., 
  # we mach genes and patients in different matrices by the row/column names 
  #   
  # Returns:
  #   A list of the parameters and the probabilities
  #
  # Date: 
  #   updated date: 2013-05-23
  #   Revised: February 15, 2015
  #
  # Author: Jiarui Ding <jiaruid@cs.ubc.ca>
  #   Department of Computer Science, UBC
  #   Department of Molecular Oncology, BC Cancer Agency 
  #
  
  gene.with.expr = colnames(expr)
  
  sample.with.mut = NULL
  if (is.null(mut)) {
    gene.with.mut = NULL
  } else {
    id = mut[, "variant_type"] %in% c("HOMD", "HLAMP", "NONSENSE", 
         "FRAMESHIFT", "NOSTOP", "SPLICE", "STARTGAINED", "CNV", 
         "FUSION")
    mut.rm = mut[id, ]
    
    id = mut.rm[, "sample"] %in% rownames(expr)
    mut.rm = mut.rm[id, ]
    
    gene.with.mut = unique(mut.rm[, "hgnc_symbol"])
    sample.with.mut = unique(mut.rm[, "sample"])
    
    ## For plot purpose
    id = mut.rm[, "variant_type"] %in% c("HOMD", "HLAMP", 'CNV')
    mut.rm.loss = mut.rm[!id, ]
  }
    
  if (length(gene) == 0) {
    gene.all = gene.with.expr
  } else {
    gene.all = intersect(gene, gene.with.expr)
  }
  gene.all = setdiff(gene.all, "") # Remove the columns without gene name
  
  # No expression data, return
  if (length(gene.all) == 0) {
    print(paste("Gene", gene, "doesn't have expression data!"))
    return
  }
  
  sample.prob.down.regulate = matrix(0, nrow=nrow(expr), 
       ncol=length(gene.all), dimnames=list(rownames(expr), gene.all))
  sample.prob.up.regulate = sample.prob.neutral = sample.prob.down.regulate
  
  name.expr.dis = c("down.regulate", "neutral", "up.regulate")
  num.comp = length(name.expr.dis)
  lambda = matrix(0, nrow=length(gene.all), ncol=num.comp, 
                  dimnames=list(gene.all, name.expr.dis))
  nu = mu = sigma = lambda 
  AIC = BIC = setNames(rep(0, length(gene.all)), gene.all)
  
  len.sample = length(sample.with.mut)
  if (len.sample > 0) {
    mut.vec = rep(0, times=len.sample)
    names(mut.vec) = sample.with.mut
  } else {
    mut.vec = NULL
  }
  
  mut.col = GenerateColourSymbol("mutation")
  mut.pch = GenerateColourSymbol(type="symbol")
  sample.common = intersect(rownames(expr), rownames(cna.call))
  
  loss.mut.vec = mut.vec
  ## For each gene, compute the conditional distribution
  for (gene in gene.all) {    
    x = expr[, gene]  
    
    # With mutation data
    if(!is.null(mut.vec)) {
      mut.vec[mut.vec > 0] = 0
      loss.mut.vec = mut.vec
      if (gene %in% gene.with.mut) { 
        id  = mut.rm[, "hgnc_symbol"] %in% gene
        sample = unique(mut.rm[id, "sample"])
              
        if(length(sample) <= 0.7*len.sample) {
          mut.vec[sample] = 1
        }
        
        id  = mut.rm.loss[, "hgnc_symbol"] %in% gene
        sample = unique(mut.rm.loss[id, "sample"])
        loss.mut.vec[sample] = 1
      } 
    }
    
    if (show.plot == TRUE) {
      close.screen(all.screens=T)
      fig = rbind(c(0.01,  0.95, 0.70, 0.95),
                  c(0.01,  0.95, 0.01, 0.69))
      split.screen(fig)
      
      screen(1)
      par(mar=c(0, 3, 0, 0))
      mix = GetMixModel(x, mut.vec, loss.mut.vec, gene, label=FALSE,
                        show.plot=show.plot, type=type, xaxt="n")
      
      screen(2)
      par(mar=c(3, 3, 0, 0))
      id  = mut[, "hgnc_symbol"] %in% gene & 
        !(mut[, "variant_type"] %in% c("HOMD", "HLAMP"))
      if(gene %in% colnames(cna.call) & length(sample.common)>0) {
        boxjitter(expr[sample.common, gene], 
                  cna.call[sample.common, gene],  
                  mut=mut[id, , drop=FALSE], 
                  mut.col=mut.col, 
                  mut.pch=mut.pch,
                  xlab=paste(gene, "mRNA expression"), 
                  ylab=expression("Copy number call"), 
                  ylim=range(x, na.rm=TRUE)) 
      }
      
      close.screen(all.screens=T)
    } else {
      mix = GetMixModel(x, mut.vec, loss.mut.vec, gene, label=FALSE,
                        show.plot=show.plot, type=type, xaxt="n")
    }
    
    ## Summarize results
    lambda.norm = mix$lambda
    mu.norm     = mix$mu
    sigma.norm  = mix$sigma
    
    if(type == "gauss") {
      sample.prob.down.regulate[, gene] = 
        dnorm(x, mean=mu.norm[1], sd=sigma.norm[1])
      sample.prob.neutral[, gene]       = 
        dnorm(x, mean=mu.norm[2], sd=sigma.norm[2])
      sample.prob.up.regulate[, gene]   = 
        dnorm(x, mean=mu.norm[3], sd=sigma.norm[3]) 
    } else if (type == "student") {
      sample.prob.down.regulate[, gene] = exp(StudentLogProb(x, mu=mu.norm[1], 
        sigma=sigma.norm[1], nu=mix$nu[1]))
      sample.prob.neutral[, gene]       = exp(StudentLogProb(x, mu=mu.norm[2], 
        sigma=sigma.norm[2], nu=mix$nu[2]))
      sample.prob.up.regulate[, gene]   = exp(StudentLogProb(x, mu=mu.norm[3], 
        sigma=sigma.norm[3], nu=mix$nu[3]))
    } else {
      stop("Error: un-recognized distributions!")
    }
        
    lambda[gene, ] = lambda.norm
    mu[gene, ]     = mu.norm
    sigma[gene, ]  = sigma.norm

    AIC[gene] = mix$AIC
    BIC[gene] = mix$BIC

    if(type == "student") {
      nu[gene, ] = mix$nu 
    }
  }
  
  return(list(sample.prob.down.regulate=sample.prob.down.regulate, 
              sample.prob.neutral=sample.prob.neutral, 
              sample.prob.up.regulate=sample.prob.up.regulate, lambda=lambda, 
              mu=mu, sigma=sigma, nu=nu, AIC=AIC, BIC=BIC))
} 

#==============================================================================
#' @importFrom stats dnorm
GetMixModel = function(x, mut, loss.mut, gene, type=type, 
                       debug=FALSE, show.plot=FALSE, ...)
{
  # MAP estimation of the mixture model parameters for a gene
  #
  # Args: 
  #   x: The mRNA expression value vector of gene (the names of x are samples)
  #   mut: The mutation vector for gene with samples as the names
  #   gene: The name of the gene 
  #   show.plot: If TRUE, plot the estimated Gaussian models
  #
  # Return:
  #   mix: The estimated Gaussian mixture model
  # 
  # ToDo:
  #   Add mutation type
  # 
  
  if(debug == TRUE) {
    cat(gene, "\t")
  }

  init.model = InitExprDis(x=x, mut=mut)
  prior = list()
  prior$lambda = init.model$lambda
  prior$mu     = init.model$mu
  prior$sigma  = init.model$sigma
  prior$sigma  = rep(max(prior$sigma), 3)
  
  # To prevent a component with zero weight and label-switch
  mult = length(x) / 5
  prior$alpha = ceiling(mult * prior$lambda)
  if(any(prior$alpha <= 1)) {
    prior$alpha = prior$alpha * 2
  } 
       
  prior$kapp = 0.1
  prior$dof  = max(length(x) * 0.05, 3)

  if(type == "gauss") {
    mix = MixGaussFitEM(x, K=3, mu=init.model$mu, sigma=prior$sigma, 
                        lambda=init.model$lambda, prior=prior, 
                        sigma.equal=TRUE) 
  } else if(type == "student") {
    mix = MixStudentFitEM(x, K=3, mu=init.model$mu, sigma=prior$sigma, 
                        lambda=init.model$lambda, prior=prior, 
                        sigma.equal=TRUE, nu.max=20, nu.equal=TRUE) 
  }
  
  if(show.plot) {
    MixModelPlot(mix, main2="", xlab2="mRNA expression", 
                 xlim=range(x, na.rm=TRUE), 
                 loglik=FALSE, ...)
    
    x = x[names(which(loss.mut > 0))]
    x = x[!is.na(x)]
    if(length(x) > 0) {
      up.regulated   = mix$lambda[3]*dnorm(x, mean=mix$mu[3], sd=mix$sigma[3])
      down.regulated = mix$lambda[1]*dnorm(x, mean=mix$mu[1], sd=mix$sigma[1]) 
      id = up.regulated >= down.regulated
      if(sum(id) > 0) {
        sfsmisc::p.arrows(x1=x[id], y1=rep(-0.05, sum(id)), x2=x[id], 
                 y2=up.regulated[id], col="red4", fill="red4", size=0.5)
      }
      if(sum(id) < length(x)) {
        sfsmisc::p.arrows(x1=x[!id], y1=rep(0, sum(!id)), x2=x[!id], size=0.5,
          y2=down.regulated[!id], col="dodgerblue", fill="dodgerblue")        
      }
      
    }
  }
  
  return(mix)
}
