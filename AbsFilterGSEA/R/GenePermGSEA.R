if(is.loaded("OnetailedGSEA",type = "External")){print("Function 'OnetailedGSEA' is properly loaded")}
if(is.loaded("TwotailedGSEA",type = "External")){print("Function 'TwotailedGSEA' is properly loaded")}

Onetailed = function(tvalue, genesetfile, min, max, nPerm, cutoff ,q){
  .Call("OnetailedGSEA", tvalue, genesetfile, min, max, nPerm, cutoff ,q)
}

Twotailed = function(tvalue, genesetfile, min, max, nPerm, cutoff ,q){
  .Call("TwotailedGSEA", tvalue, genesetfile, min, max, nPerm, cutoff ,q)
}

snr = function(value, g1, g2){
  meandiff = mean(value[g1])-mean(value[g2])
  sdsum = sd(value[g1])+sd(value[g2])
  result = meandiff/sdsum
  return(result)
}

moderated_t = function(mat, g1, g2){
  arr = rep(0, length(g1)+length(g2))
  arr[g1] = 1
  design = cbind(Intersect = 1, x = arr)
  fit = lmFit(mat, design = design)
  eb = eBayes(fit)
  moderated_tvalue = eb$t[,2]
  return(moderated_tvalue)
}

foldchange = function(value, g1, g2)
{
  val1 = mean(value[g1])
  val2 = mean(value[g2])
  if(val1 == 0){val1 = val1+(0.0001/length(g1))}
  if(val2 == 0){val2 = val2+(0.0001/length(g2))}
  result = log2(val1/val2)
  return(result)
}

ranksum = function(value, g1, g2)
{
  ORDER = rank(value, ties.method='average')
  len1 = length(g1)
  len2 = (length(g1)+length(g2))+1
  Tvalue = sum(ORDER[g1])-(len1*len2/2)
  return(Tvalue)
}

#' Gene permuting GSEA with or without filtering by absolute GSEA.
#'
#' Gene-set enrichment analysis (GSEA) is popularly used to assess the enrichment of differential signal in a pre-defined gene-set without using a cutoff threshold for differential expression. The significance of enrichment is evaluated through sample- or gene-permutation method. Although the sample-permutation approach is highly recommended due to its good false positive control, we must use gene-permuting method if the number of samples is small. However, such gene-permuting GSEA (or preranked GSEA) generates a lot of false positive gene-sets as the inter-gene correlation in each gene set increases. These false positives can be successfully reduced by filtering with the one-tailed absolute GSEA results. This package provides a function that performs gene-permuting GSEA calculation with or without the absolute filtering. Without filtering, users can perform (original) two-tailed or one-tailed absolute GSEA.
#'
#' @param countMatrix Normalized RNA-seq read count matrix.
#'
#' @param GeneScoreType Type of gene score. Possible gene score is "moderated_t","SNR", "FC" (log fold change score) or "RANKSUM" (zero centered).
#'
#' @param idxCase Indices for case samples in the count matrix. e.g., 1:3
#'
#' @param idxControl Indices for control samples in the count matrix. e.g., 4:6
#'
#' @param GenesetFile File path for gene set file. Typical GMT file or its similar 'tab-delimited' file is available. e.g., "C:/geneset.gmt"
#'
#' @param normalization Type 'DESeq' if the input matrix is composed of raw read counts. It will normalize the raw count data using DESeq method. Or type 'AlreadyNormalized' if the input matrix is already normalized.
#'
#' @param minGenesetSize Minimum size of gene set allowed. Gene-sets of which sizes are below this value are filtered out from the analysis. Default = 10
#'
#' @param maxGenesetSize Maximum size of gene set allowed. Gene-sets of which sizes are larger this value are filtered out from the analysis. Default = 300
#'
#' @param q Weight exponent for gene score. For example, if q=0, only rank of gene score is reflected in calculating gene set score (preranked GSEA). If q=1, the gene score itself is used. If q=2, square of the gene score is used.
#'
#' @param nPerm The number of gene permutation. Default = 1000.
#'
#' @param GSEAtype Type of GSEA. Possible value is "absolute", "original" or "absFilter". "absolute" for one-tailed absolute GSEA. "original" for the original two-tailed GSEA. "absFilter" for the original GSEA filtered by the results from the one-tailed absolute GSEA.
#'
#' @param FDR FDR cutoff for the original or absolute GSEA. Default = 0.05
#'
#' @param FDRfilter FDR cutoff for the one-tailed absolute GSEA for absolute filtering (only working when GSEAtype is "absFilter"). Default = 0.05
#'
#' @param minCount Minimum median count of a gene to be included in the analysis. It is used for gene-filtering to avoid genes having small read counts. Default = 0
#'
#' @import Rcpp
#'
#' @import DESeq
#'
#' @importFrom limma lmFit
#'
#' @importFrom limma eBayes
#'
#' @importFrom stats sd
#'
#' @importFrom stats median
#'
#' @importFrom stats quantile
#'
#' @return GSEA result table sorted by FDR Q-value.
#'
#' @examples
#'
#' data(example)
#'
#' # Create a gene set file and save it to your local directory.
#' # Note that you can use your local gene set file (tab-delimited like *.gmt file from mSigDB).
#' # But here, we will generate a toy gene set file to show the structure of this gene set file.
#' # It consists of 50 gene sets and each contains 100 genes.
#'
#' for(Geneset in 1:50)
#' {
#'   GenesetName = paste("Geneset", Geneset, sep = "_")
#'   Genes = paste("Gene", (Geneset*100-99):(Geneset*100), sep="", collapse = '\t')
#'   Geneset = paste(GenesetName, Genes, sep = '\t')
#'   write(Geneset, file = "geneset.txt", append = TRUE, ncolumns = 1)
#' }
#'
#' # Run Gene-permuting GSEA
#' RES = GenePermGSEA(countMatrix = example, GeneScoreType = "moderated_t", idxCase = 1:3,
#'                     idxControl = 4:6, GenesetFile = 'geneset.txt', normalization = 'DESeq',
#'                     GSEAtype = "absFilter")
#' RES
#'
#' @export
#'
#' @details Typical usages are
#' GenePermGSEA(countMatrix = countMatrix, GeneScoreType = "moderated_t", idxCase = 1:3,
#'                    idxControl = 4:6, GenesetFile = 'geneset.txt', GSEAtype = "absFilter")
#'
#' @source Nam, D. Effect of the absolute statistic on gene-sampling gene-set analysis methods. Stat Methods Med Res 2015.
#' Subramanian, A., et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. P Natl Acad Sci USA 2005;102(43):15545-15550.
#' Li, J. and Tibshirani, R. Finding consistent patterns: A nonparametric approach for identifying differential expression in RNA-Seq data. Statistical Methods in Medical Research 2013;22(5):519-536.
#'
#' @references Nam, D. Effect of the absolute statistic on gene-sampling gene-set analysis methods. Stat Methods Med Res 2015.
#' Subramanian, A., et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. P Natl Acad Sci USA 2005;102(43):15545-15550.
#' Li, J. and Tibshirani, R. Finding consistent patterns: A nonparametric approach for identifying differential expression in RNA-Seq data. Statistical Methods in Medical Research 2013;22(5):519-536.
#' Simon Anders and Wolfgang Huber (2010): Differential expression analysis for sequence count data. Genome Biology 11:R106
#'
#' @useDynLib AbsFilterGSEA
GenePermGSEA = function(countMatrix, GeneScoreType, idxCase, idxControl, GenesetFile, normalization, minGenesetSize=10, maxGenesetSize=300, q=1, nPerm=1000, GSEAtype="absFilter", FDR=0.05, FDRfilter=0.05, minCount=3)
{
  dimMat = try(dim(countMatrix))
  if(is.null(dimMat)){stop("The dimension of input count matrix is NULL.")}
  if(dimMat[1]*dimMat[2] == 0){stop("Count matrix must have positive dimension.")}
  if(GeneScoreType!="moderated_t" & GeneScoreType!="SNR" & GeneScoreType!="FC" & GeneScoreType!="RANKSUM"){stop("Gene score type must be 'SNR', 'FC' or 'RANKSUM'.")}
  if(length(idxCase)<1 | length(idxControl)<1){stop("idxCase and idxControl must be positive integer.")}
  if(!file.exists(GenesetFile)){stop(paste(GenesetFile, " : Such gene set file does not exist.", sep=""))}
  if(GSEAtype!="absolute" & GSEAtype !="original" & GSEAtype!="absFilter"){stop("GSEAtype must be 'absolute' (for absolute one-tailed GSEA), 'original' (both up and down direction(=two-tailed GSEA)) or 'absFilter' (Result for two-tailed GSEA filtered by one-tailed GSEA result).")}
  if(FDR<0){stop('FDR must not be negative value.')}
  if(GSEAtype == "absFilter" & is.null(FDRfilter)){stop("FDRfilter must be set if GSEAtype is 'absFilter'.")}
  if(!is.null(FDRfilter)){if(FDRfilter < 0){stop("FDRfilter must not be negative value.")}}
  if(is.null(q)){stop('q (weight exponent of gene score) is empty')
    }else
      {
        if(q<0){stop("Negative value for q (Weight exponent of gene score) is not allowed.")
        }else if(q-as.integer(q) == 0){qType = "integer"
        }else if(q-as.integer(q) != 0){qType = "decimal"
        }else {stop("Input proper q (weight exponent of enrichment score.)")}
      }


  countMatrix = data.matrix(countMatrix)

  # DESeq normalization
  if(normalization == 'DESeq')
  {
    condition = array(rep(1, length(c(idxCase, idxControl))))
    condition[idxControl] = 2
    cds = newCountDataSet(countMatrix, condition)
    cds = estimateSizeFactors(cds)
    sf = sizeFactors(cds)
    countMatrix = sweep(countMatrix, 2, sf, "/")
  }

  # Filtering genes having small read number and add pseudocount (5% quantile of positive counts).
  medianReads = apply(countMatrix, 1, median)
  index.deletion = which(medianReads < minCount)
  if(length(index.deletion)>0){countMatrix = countMatrix[-index.deletion,]}
  pseudoCount = quantile(countMatrix, 0.05)
  countMatrix = countMatrix + pseudoCount

  # Gene score
  if(GeneScoreType == 'SNR')
  {
    FUNC = snr
  }
  if(GeneScoreType == 'FC')
  {
    FUNC = foldchange
  }
  if(GeneScoreType == 'RANKSUM')
  {
    FUNC = ranksum
  }
  genescore = if(sum(c("SNR","FC","RANKSUM")%in%GeneScoreType) == 1){
    try(apply(countMatrix, 1, FUN = FUNC, g1 = idxCase, g2 = idxControl), silent = T)
    }else if(GeneScoreType == 'moderated_t'){
      try(moderated_t(mat = countMatrix, g1 = idxCase, g2 = idxControl), silent = T)}

  if(class(genescore)=='try-error'){stop("Invalid gene scores. Please check idxCase and idxControl (indices for case and control samples, respectively).")}

  if(q!=0 & qType == "integer"){genescore = genescore^q}
  if(q!=0 & qType == "decimal"){genescore = abs(genescore)^q}
  if(q==0){genescore = genescore}

  genescore = sort(genescore, decreasing = TRUE)

  if(GSEAtype == "absolute" | GSEAtype == "absFilter"){genescore_abs = abs(genescore); genescore_abs = sort(genescore_abs, decreasing = TRUE)}

  # GSEA
  if(GSEAtype == "absolute")
  {
    Result_table = Onetailed(genescore_abs, GenesetFile, minGenesetSize, maxGenesetSize, nPerm, FDR, q)
    Result_table = Result_table[order(Result_table[[5]]),]
    #Result_table$NES = format(Result_table$NES, scientific = FALSE)
    Result_table$Nominal.P.value = as.numeric(Result_table$Nominal.P.value)
    Result_table$FDR.Q.value = signif(Result_table$FDR.Q.value,)
    return(Result_table)
  }

  if(GSEAtype == "original")
  {
    Result_table = Twotailed(genescore, GenesetFile, minGenesetSize, maxGenesetSize, nPerm, FDR, q)
    Result_table = Result_table[order(Result_table[[5]]),]
   # Result_table$NES = format(Result_table$NES, scientific = FALSE)
    Result_table$Nominal.P.value = as.numeric(Result_table$Nominal.P.value)
    Result_table$FDR.Q.value = signif(Result_table$FDR.Q.value, 4)
    return(Result_table)
  }

  if(GSEAtype == "absFilter")
  {
    Result_table_abs = Onetailed(genescore_abs, GenesetFile, minGenesetSize, maxGenesetSize, nPerm, FDRfilter, q)
    Result_table_ord = Twotailed(genescore, GenesetFile, minGenesetSize, maxGenesetSize, nPerm, FDR, q)
    Filtered = which(Result_table_ord$GenesetName%in%Result_table_abs$GenesetName)
    Result_table = Result_table_ord[Filtered,]
    Result_table = Result_table[order(Result_table[[5]]),]
    #Result_table$NES = format(Result_table$NES, scientific = FALSE)
    Result_table$Nominal.P.value = as.numeric(Result_table$Nominal.P.value)
    Result_table$FDR.Q.value = signif(Result_table$FDR.Q.value,4)
    return(Result_table)
  }
}
