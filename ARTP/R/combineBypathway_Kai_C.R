# History Nov 05 2010 Initial coding
#         Nov 17 2010 For R^2 option calculate 1 matrix using all snps in the pathway
#                     and then subset.
#         Jan 12 2011 Add function for the R package
#                     Remove pathway.combine and pathway.results
#         Jan 18 2011 Add check for the most significant snp/gene
#         Mar 04 2011 
#         Mar 22 2011 Let gene.list be NULL to allow all snps in same gene
#         Feb 23 2012 Add option to work with new method (perm.list = NULL)
#         Feb 04 2013 Fix bug in runPermutations with miss.rate option

#########################################################################################
#########################################################################################
#########################################################################################
# NOTE: The observed and permutation files must contain only numeric data (no NA, NaN, Inf, etc

# NOTE: To get the top set of genes for the pathway
# $p.for.minp: is the final pathway p-values
# $p.obs: is the vector of p-values for each check at the given truncation point (you 
#   need to know the configuration, i.e, percentage, and number of gene).
# $t.obs: is the p-value for each gene is the order given by $gene
# To find the best set of genes, check $p.obs to see which element is the smallest
# Let k be the index with the smallest value, based on the percentage, you can figure out
#   how many genes are used at the kth checkup, let this number be g.
# The top g ranked gene (according to $t.obs, the smaller the better) is the final set.
#########################################################################################
#########################################################################################
#########################################################################################


# Function to compute gene-level permutation p-values. 
combineByPathway_Kai <- function(perm.list, obs.list, gene.list,  
                          pathway.list, op=NULL) {

  ###################################################################
  # perm.list         A list with the following names
  #   file            If all permutations are in 1 file, then specify 
  #                   the complete path to the file.
  #                   The default is NULL.
  #   file.type       1 or 3
  #                   1: R object file created with the save() function.
  #                   The object should be a matrix with the column names
  #                   as the SNPs.
  #                   3: Flat file with the SNP names as the header.
  #                   The default is 3
  #   delimiter       The delimiter used in the files.
  #                   The default is " "
  #   pvalues         0 or 1 to denote whether the permutation values
  #                   are p-values.
  #                   The default is 1.
  #   nperm           Number of permutations (if known)
  #                   The default is NULL
  ###################################################################
  # obs.list         A list with the folowing names
  #   file            File containing the observed data
  #                   No default
  #   file.type       1 or 3
  #                   1: R object file created with the save() function.
  #                   The object must be a data frame.
  #                   3: Flat file with a header of variable names.
  #                   The default is 1
  #   delimiter       The delimiter used in the file. (Only for file.type=3)
  #                   The default is " "
  #   pvalues         0 or 1 to denote whether the observed values 
  #                   are p-values.
  #                   The default is 1.
  ###################################################################   
  # gene.list        A list with the following names
  #                  If NULL, then all snps belong to the same gene
  #   file            File containing the SNP names and gene names.
  #                   No default
  #   header          0 or 1 if the file contains a header
  #                   The default is 1
  #   gene.var        Variable name or column number of the gene variable
  #                   The default is "Gene"
  #   snp.var         Variable name or column number of the group variable.
  #                   The default is "SNP"
  #   file.type       1 or 3
  #                   The default is 3
  #   delimiter       The default is "\t"
  ###################################################################   
  # pathway.list      A list with the following names
  #   file            File containing the SNP names and gene names.
  #                   No default
  #   header          0 or 1 if the file contains a header
  #                   The default is 1
  #   gene.var        Variable name or column number of the gene variable
  #                   The default is "Gene"
  #   file.type       1 or 3
  #                   The default is 3
  #   delimiter       The default is "\t"
  #   name            Name of the pathway
  ###################################################################
  # op                     A list with the following names 
  #   out.file             Output file for pathway results.
  #                        The default is NULL.
  #   out.gene             Output file to store the summary info for each gene.
  #                        The default is NULL.
  #   temp.list         
  #   inspect.fun          The default is fun.ck
  #   inspect.snp.percent  The default is 1e-10
  #   inspect.snp.n        1 or 3-element vector for max(a[1], min(a[2], floor(N.SNPs/a[3])))
  #                        The default is 1
  #   inspect.gene.percent The default is 0.05
  #   inspect.gene.n       The default is 10  
  #   missing              Missing values in obs and permutation files
  #                        The default is NA
  #   snps.list            List of type file.list of SNPs to keep/remove
  #                        Has options toupper, tolower for uppercase/lowercase
  #                        The default is NULL
  #   genes                Subset of genes for parallel computing
  #                        The default is NULL  
  #   r2.genoFile          The name of the genotype file for computing R^2 matrices
  #                        This must be a file that GLU can read
  #                        The default is NULL
  #   r2.threshold         Only used if r2.genoFile is specified
  #                        The default value is 0.8
  #   r2.inspect           0 or 1 to change the values of inspect.snp.n and
  #                        inspect.snp.percent based on the number of blocks.
  #                        The default is 1.
  #   r2.min.MAF   
  #   ties.method          0 or 1  for breaking ties (0=random, 1=first)
  #                        The default is 0.
  #   pathway.method       0 or 1 for pathway p-value computation
  #                        
  #                        The default is 0
  #   parms.list           List of info for each gene parameter for Jianxin's method
  #     theta              Value of the theta parm. The default is 0.2.
  #     file               File containing phi for each gene. parm.var and gene.var must be specified.
  #    

  t00 <- proc.time()

  glFlag   <- !is.null(gene.list)
  permFlag <- !is.null(perm.list)
  if (!permFlag) stop("ERROR: perm.list is NULL")

  # Check the lists
  obs.list  <- check.obs.list(obs.list)
  if (glFlag) gene.list <- check.gene.list(gene.list)
  if (permFlag) {
    perm.list <- default.list(perm.list, 
      c("file.type", "delimiter", "pvalues"),
      list(3, ",", 1),
      error=c(0, 0, 0))
    dir <- perm.list[["dir", exact=TRUE]]
    if (!is.null(dir)) {
      pattern   <- perm.list[["pattern", exact=TRUE]]
      permFiles <- list.files(path=dir, pattern=pattern)
      dir       <- checkForSep(dir)
      permFiles <- paste(dir, permFiles, sep="")
      rm(dir, pattern)
    } else {
      permFiles <- perm.list[["file", exact=TRUE]]
      if (is.null(permFiles)) {
        stop("ERROR: perm.list$dir or perm.list$file must be specified")
      }
    }
  }
  pathway.list <- default.list(pathway.list, 
      c("gene.var", "file.type", "delimiter", "header", "name"),
      list(1, 3, "\t", 1, "pathway"),
      error=c(0, 0, 0, 0, 0))

  op <- default.list(op,
   c("inspect.snp.percent", "inspect.snp.n", "inspect.gene.percent", "inspect.gene.n", 
     "missing", "TEST",  "method", "r2.threshold", "TEST.gene", "r2.inspect", "r2.min.MAF", 
     "temp.list", "ties.method", "pathway.method", "mvn.maxn"), 
   list(1e-10, 1, 0.05, 10, NA, 0, 1, 0.8, 1, 1, 0.05, list(), 0, 0, 1))
  TEST <- op$TEST
  temp.list <- op$temp.list
  temp.list <- check.temp.list(temp.list)

  r2.genoFile <- op[["r2.genoFile", exact=TRUE]]
  r2Flag      <- !is.null(r2.genoFile)

  # Get the genes and snps
  if (glFlag) {
    data  <- getColumns(gene.list, c(gene.list$gene.var, gene.list$snp.var))
    genes <- as.character(data[[gene.list$gene.var]])
    snps  <- as.character(data[[gene.list$snp.var]])
    genes <- removeWhiteSpace(genes)
    snps  <- removeWhiteSpace(snps)

    rm(data, gene.list)
    temp <- gc()
  } else {
    snps  <- scan(obs.list$file, what="character", sep=obs.list$delimiter, nlines=1, quiet=TRUE)
    genes <- rep("gene", times=length(snps))  
  }
  
  subset <- op[["genes", exact=TRUE]]
  if (!is.null(subset)) {
    subset <- removeWhiteSpace(subset)
    temp   <- genes %in% subset
    genes  <- genes[temp]
    snps   <- snps[temp] 
  }

  # Read in the pathway file
  temp   <- pathway.list[["file", exact=TRUE]]
  if (!is.null(temp)) {
    temp   <- loadData.table(pathway.list)
    subset <- makeVector(temp[, pathway.list$gene.var])
    subset <- removeWhiteSpace(subset)
    subset <- unique(subset)
    ###################################################
    # For testing
    if ((op$TEST) && (r2Flag)) subset <- subset[op$TEST.gene]
    ###################################################
    temp   <- !(subset %in% c(NA, "", " "))
    subset <- subset[temp]
    temp   <- genes %in% subset
    genes  <- genes[temp]
    snps   <- snps[temp]
  }
  
  # Keep/remove snps
  temp <- op[["snps.list", exact=TRUE]]
  if (!is.null(temp)) {
    temp <- default.list(temp, c("file", "tolower", "toupper"), 
                 list("ERROR", 0, 0), error=c(1, 0, 0))
    subset <- removeWhiteSpace(getIdsFromFile(temp))
    if (temp$tolower) subset <- tolower(subset)
    if (temp$toupper) subset <- toupper(subset)
    temp   <- snps %in% subset
    genes  <- genes[temp]
    snps   <- snps[temp]
  }

  rm(subset)

  # Get the unique genes
  ugenes <- unique(genes)

  # Get the number of genes
  ngenes <- length(ugenes) 

  # Get the observed vectors
  fid <- file(obs.list$file, "r")
  temp <- scan(fid, what="character", sep=obs.list$delimiter, nlines=1, quiet=TRUE)
  obs  <- scan(fid, what=double(0), sep=obs.list$delimiter, nlines=1, quiet=TRUE)
  names(obs) <- temp
  close(fid)
  nSNPsInFile <- length(obs)
  snp0 <- temp

  # Get the size of row 1 and row 2
  rowsize <- integer(2)
  temp <- scan(obs.list$file, what="character", sep="\n", nlines=2, quiet=TRUE)
  rowsize[1] <- nchar(temp[1]) + 100
  rowsize[2] <- nchar(temp[2]) + 1000

  # Change pvalues to NA if > 1 or < 0
  if (obs.list$pvalues) {
    temp <- (obs > 1) | (obs < 0)
    temp[is.na(temp)] <- TRUE
    obs[temp] <- NA 
  }

  # Missing values
  missing <- op$missing
  obs     <- removeMiss(obs, miss=missing)

  # Pvalue will be log transformed, so watch out for zeros
  temp <- obs == 0
  obs[temp] <- 1e-300 

  # Get the snps for each gene
  temp <- geneSNPs.study(list(names(obs)), snps, genes, method=2)
  csnps.list <- temp$csnps.list

  rm(snps, genes)
  temp <- gc()

  # Change values if p-values
  if (obs.list$pvalues) obs <- -1.0*log(obs)
  
  # Get the column names in the perm file
  if (permFlag) {
    temp <- list(delimiter=perm.list$delimiter, return.cols=1)
    snp0 <- getNcols(permFiles[1], temp) 

    # Get the number of permuations
    nperm <- perm.list[["nperm", exact=TRUE]]
    if (is.null(nperm)) {
      print("Obtaining the number of permutations")
      nperm <- getNperm(perm.list)
      print(paste("nperm = ", nperm, sep=""))
    }
  } else {
    nperm <- op[["nperm", exact=TRUE]]
  }

  geneSnpMat   <- rep(0,2)
  obs.all      <- c()
  snp.all      <- c()
  mygene.start <- 1
  gene.used    <- character()

  for (gg in ugenes)
  {
    snps  <- csnps.list[[gg]]
    nsnps <- length(snps)
    if (nsnps < 1) next

    gene.used    <- c(gene.used, gg)
    mygene.end   <- mygene.start + nsnps -1
    geneSnpMat   <- rbind(geneSnpMat, c(mygene.start, mygene.end))   
    mygene.start <- mygene.end + 1
    obs.all      <- c(obs.all, obs[snps])
    snp.all      <- c(snp.all, snps)
  }

  geneSnpMat <- removeOrKeepRows(geneSnpMat, 1, which=-1)
  geneStart  <- geneSnpMat[, 1]
  geneStop   <- geneSnpMat[, 2] 
  num.gene   <- length(geneStart)

  # Get the positions of the snps in the observed and permutation files
  snpcols <- match(snp.all, snp0)
  if (any(is.na(snpcols))) stop("ERROR: with SNP column names")
  nsnpcols <- length(snpcols)

  rm(geneSnpMat, mygene.start, mygene.end, snp0, obs)
  gc()

  # Get the gene parms
  pathMethod <- op$pathway.method
  if (pathMethod) {
    gene.parms <- getGene.parms(op[["parms.list", exact=TRUE]], gene.used)
  } else {
    gene.parms <- 0
  }

  inspect0 <- op$inspect.snp.n
  inspect2Flag   <- length(inspect0) == 3
  if (inspect2Flag) op$r2.inspect <- 0
  nblocks <- NULL

  # Compute R^2 matrices. They will be ordered by the list of snps by default.
  r2.inspect <- NULL
  r2Files    <- ""
  if (r2Flag) {
    r2Files <- character(num.gene)
    if (op$r2.inspect) {
      r2.inspect <- integer(num.gene)
      op$min.MAF <- op$r2.min.MAF
    }
    temp <- unique(snp.all)
    mat  <- glu.ldMatrix(temp, r2.genoFile, op=op)
    i    <- nrow(mat)
    mat  <- as.integer(mat < op$r2.threshold)
    dim(mat) <- c(i, i)
    rownames(mat) <- temp
    colnames(mat) <- temp

    i <- 1
    tid <- temp.list$id
    for (gg in ugenes) {
      snps  <- csnps.list[[gg]]
      nsnps <- length(snps)
      if (nsnps < 1) next
    
      # Get the R^2 matrix
      #temp <- glu.ldMatrix(snps, r2.genoFile, op=op)
      #temp <- as.integer(temp < op$r2.threshold)
      #dim(temp) <- c(nsnps, nsnps)
      temp <- mat[snps, snps] 

      # Store in a temporary file
      r2Files[i] <- getTempfile(temp.list$dir, prefix=paste("r2mat_", tid, "_", i, "_", sep=""), ext=".txt")
      write.table(temp, file=r2Files[i], sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)

      if (op$r2.inspect) {
        nblocks <- glu.nBins(snps, r2.genoFile, op=op)
        if (nblocks < 1) nblocks <- 1
        r2.inspect[i] <- nblocks
      }
      i <- i + 1 
    }
    rm(mat, temp)
    temp <- gc()
  } 
  rm(gg, ugenes, csnps.list, nblocks, nsnps, temp, snps)
  gc()

  if (TEST) nperm  <- max(10, TEST)

  nper <- op$inspect.snp.percent
  if (nper > 1) nper <- nper/100
  nper  <- 1/nper
  
  inspect.vec.gene.list <- NULL
  iStart <- integer(num.gene)
  iStop  <- integer(num.gene)
  a      <- 1
  inspectVec <- integer(num.gene)
  for (i in 1:num.gene)
  {
     my.gene.size <- geneStop[i] - geneStart[i] + 1

     if (inspect2Flag) {
       snp.n <- max(inspect0[1], min(inspect0[2], floor(my.gene.size/inspect0[3])))
     } else if ((r2Flag) && (op$r2.inspect)) {
       snp.n <- r2.inspect[i]
     } else {
       snp.n <- op$inspect.snp.n
     }
     inspectVec[i] <- snp.n 
     my.step  <- floor(my.gene.size/nper)
     my.step  <- max(1, my.step) 
     up.limit <- min(my.gene.size, snp.n*my.step)
     up.limit <- max(up.limit, my.step)
     ## if want to use minP for each gene, then use the following
     #    inspect.vec.gene.list[[i]] <-c(1)
     temp <- seq(from=my.step, to=up.limit, by=my.step)
     inspect.vec.gene.list <- c(inspect.vec.gene.list, temp)
     iStart[i] <- a
     iStop[i]  <- a + length(temp) - 1
     a         <- iStop[i] + 1  
  }   

  nper.gene <- op$inspect.gene.percent
  if (nper.gene >= 1) nper.gene <- nper.gene/100
  nper.gene <- 1/nper.gene

  my.step  <- floor(num.gene/nper.gene)
  my.step  <- max(1, my.step)
  up.limit <- min(num.gene, op$inspect.gene.n*my.step)
  inspect.vec.pathway <- seq(from=my.step, to=up.limit, by=my.step)

  # Get a temporary file for the gene results
  out.gene <- op[["out.gene", exact=TRUE]]
  outgFlag <- !is.null(out.gene)
  tmpfile  <- getTempfile(temp.list$dir, prefix=paste("pathgene_", temp.list$id, "_", sep=""), ext=".txt")

  #if (outgFlag) {
  #  tmpfile <- getTempfile(temp.list$dir, prefix=paste("pathgene_", temp.list$id, "_", sep=""), ext=".txt")
  #} else {
  #  tmpfile <- ""
  #} 

  # Set up argument for the C function
  nperm       <- as.integer(nperm)
  nSNPsInFile <- as.integer(nSNPsInFile)
  snpcols     <- as.integer(snpcols - 1)
  nsnpcols    <- as.integer(nsnpcols)
  geneStart   <- as.integer(geneStart - 1)
  geneStop    <- as.integer(geneStop - 1)
  num.gene    <- as.integer(num.gene)
  rowsize     <- as.integer(rowsize)
  inspect.vec.gene.list <- as.integer(inspect.vec.gene.list - 1)
  llen        <- as.integer(length(inspect.vec.gene.list))
  iStart      <- as.integer(iStart - 1)
  iStop       <- as.integer(iStop - 1)
  inspect.vec.pathway <- as.integer(inspect.vec.pathway - 1)
  ilen        <- as.integer(length(inspect.vec.pathway))
  ret_p_obs   <- double(ilen)
  ret_minp    <- double(1)
  ret_nperm   <- integer(1)
  r2Flag      <- as.integer(r2Flag)
  ties.method <- as.integer(op$ties.method)
  pathMethod  <- as.integer(pathMethod)
  gene.parms  <- as.double(gene.parms)
  ret_code    <- as.integer(0)

  if (permFlag) {
    temp <- try(.C("ARTP_pathway", obs.list$file, perm.list$file, nperm, nSNPsInFile, snpcols, nsnpcols, geneStart, geneStop, num.gene, rowsize,
                  inspect.vec.gene.list, llen, iStart, iStop, inspect.vec.pathway, ilen, r2Files, r2Flag, tmpfile, ties.method,
                  pathMethod, gene.parms, ret_p_obs=ret_p_obs, ret_minp=ret_minp, ret_nperm=ret_nperm, PACKAGE="ARTP"))
  } else {
    maxn <- as.integer(op$mvn.maxn)
    func <- list(op$mvn.function)
    temp <- try(.C("pathway2", maxn, func, obs.list$file, nperm, nSNPsInFile, snpcols, nsnpcols, geneStart, geneStop, num.gene, rowsize,
                  inspect.vec.gene.list, llen, iStart, iStop, inspect.vec.pathway, ilen, r2Files, r2Flag, tmpfile, ties.method,
                  pathMethod, gene.parms, ret_p_obs=ret_p_obs, ret_minp=ret_minp, ret_nperm=ret_nperm, ret_code=ret_code))
    ret_code <- temp$ret_code
  }
  if ((checkTryError(temp, conv=0)) || (ret_code != 0)) {
    print(temp)
    stop("ERROR: calling C function")
  }
  ret_nperm <- temp$ret_nperm
  ret_minp  <- temp$ret_minp
  ret_p_obs <- temp$ret_p_obs
  rm(temp)
  gc()

  # Delete temporary files
  if ((r2Flag) && (temp.list$delete)) {
    for (i in 1:num.gene) file.remove(r2Files[i])
  }

  geneStart <- geneStart + 1
  geneStop  <- geneStop + 1 
  iStart    <- iStart + 1
  iStop     <- iStop + 1
  inspect.vec.pathway <- inspect.vec.pathway + 1
  inspect.vec.gene.list <- inspect.vec.gene.list + 1

  temp <- c("Gene", "N.SNP", "Pvalue", "Most.Significant.SNPs", "Inspect.SNP.N")
  final <- initDataFrame(num.gene, c("C", "N", "N", "C", "N"), colnames=temp)
  final[, "Gene"] <- gene.used
  final[, "Inspect.SNP.N"] <- inspectVec

  # Open the gene summary file
  fid <- getFID(tmpfile, list(open="r")) 
  path.t.obs <- double(num.gene)   

  # Obs.all are -log(pvalues)
  obs.all <- exp(-obs.all)
 
  # The order of the file is first row, final p-value, second row p_obs, third row p_perm
  for (i in 1:num.gene) {
    
    x <- scan(fid, what="character", sep="\n", nlines=3, quiet=TRUE)
    path.t.obs[i] <- as.numeric(getVecFromStr(x[3], delimiter=",")[1])

   
      final[i, "Pvalue"] <- as.numeric(x[1]) 
      t.obs <- obs.all[geneStart[i]:geneStop[i]]
      final[i, "N.SNP"] <- length(t.obs)     
   
      p.obs <- as.numeric(getVecFromStr(x[2], delimiter=","))
      kk    <- which.min(p.obs)
      ivec  <- inspect.vec.gene.list[iStart[i]:iStop[i]]
      ns    <- ivec[kk]
      rnk   <- rank(t.obs, ties.method="random")
      temp  <- rnk %in% 1:ns
      most.sig.snps <- names(t.obs)[temp]
      final[i, "Most.Significant.SNPs"] <- paste(most.sig.snps, collapse=", ", sep="")

      # Check for most sig snp
      k   <- which.min(t.obs)
      snp <- names(t.obs[k])
      #if (!(snp %in% most.sig.snps)) stop("ERROR: with most significant SNPs") 
    
  }
  close(fid)
  if (temp.list$delete) file.remove(tmpfile)
  if (outgFlag) writeTable(final, out.gene)
  path.t.obs <- exp(-path.t.obs)
  
  # Create return object
  ret                      <- list(final=final, p.obs=ret_p_obs, p.for.minp=ret_minp, t.obs=path.t.obs)
  ret$nperm                <- ret_nperm
  ret$gene                 <- gene.used
  ret$geneSnpMat           <- cbind(geneStart, geneStop)
  ret$snp.list             <- snp.all
  ret$inspect.vec.pathway  <- inspect.vec.pathway
  ret$inspect.snp.percent  <- op$inspect.snp.percent
  ret$inspect.snp.n        <- op$inspect.snp.n
  ret$inspect.gene.percent <- op$inspect.gene.percent
  ret$inspect.gene.n       <- op$inspect.gene.n

  # Get the best set of genes
  k     <- which.min(ret_p_obs)
  ng    <- inspect.vec.pathway[k]
  rnk   <- rank(ret$t.obs, ties.method="random")
  temp  <- rnk %in% 1:ng
  ret$most.sig.genes <- ret$gene[temp]

  # Check the most significant gene
  k  <- which.min(ret$t.obs)
  gg <- gene.used[k]
  #if (!(gg %in% ret$most.sig.genes)) stop("ERROR: in most significant genes") 

  outfile <- op[["out.file", exact=TRUE]]
  if (!is.null(outfile)) save(ret, file=outfile)

  if (op$TEST) print(ret)
  
  ret$nperm <- ret_nperm
  #print(paste("nperm = ", ret_nperm, sep=""))
  #print((proc.time()-t00)/60)
  
  ret

} # END: combineByPathway_Kai


# Function for the R package
ARTP_pathway <- function(obs.file, perm.file, nperm, temp.dir, gene.list=NULL, op=NULL) {

  if (!is.null(gene.list)) {
    gene.list$file.type <- 3
    gene.list <- default.list(gene.list, c("snp.var", "gene.var"), list("SNP", "Gene"))
    vars <- c(gene.list[["gene.var", exact=TRUE]], gene.list[["snp.var", exact=TRUE]])
    gene.list <- check.file.list(gene.list, op=list(exist=1, vars=vars))
  }
  if (is.null(op)) op <- list()
  op <- default.list(op, 
         c("inspect.snp.n", "inspect.snp.percent", "inspect.gene.n", "inspect.gene.percent"), 
         list(1, 0, 10, 0.05))
  temp.list <- list(dir=temp.dir, delete=1)
  op$temp.list <- temp.list
  perm.list <- list(file=perm.file, file.type=3, delimiter=",", nperm=nperm, pvalues=1)
  obs.list  <- list(file=obs.file, file.type=3, delimiter=",", pvalues=1)
  temp <- combineByPathway_Kai(perm.list, obs.list, gene.list,  
                          list(), op=op)
  
  ret <- list(pathway.pvalue=temp$p.for.minp)
  temp2 <- removeOrKeepCols(temp$final, 1:3)
  ret$gene.table <- temp2
  ret$nperm <- temp$nperm

  ret

} # END: ARTP_pathway

# Function to get the snp prior parms
getGene.parms <- function(parms.list, geneNames) {

  parms.list <- default.list(parms.list, c("theta", "phi"), list(0.2, 1.0))
  theta <- parms.list$theta

  f <- parms.list[["file", exact=TRUE]]
  if (!is.null(f)) {
    x     <- loadData.table(parms.list)
    genes <- removeWhiteSpace(x[, parms.list$gene.var])
    parms <- as.numeric(x[, parms.list$parm.var])
    temp  <- genes %in% geneNames
    genes <- genes[temp]
    parms <- parms[temp]
    parms <- parms*theta
    names(parms) <- genes 
  } else {
    parms <- rep(1.0, length(geneNames))
    parms <- parms*theta
    names(parms) <- geneNames
  }

  parms

} # END: getGene.parms

# Function to call test4
runPermutations <- function(snp.list, pheno.list, family, op=NULL) {

  # family 1 = binomial, 2 = gaussian

  if (!(family %in% 1:2)) stop("ERROR: family must be 1 (logistic regression) or 2 (linear regresion)")

  snp.list <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)

  if (is.null(op)) op <- list()
  op <- default.list(op, 
      c("nperm", "obs.outfile", "perm.outfile", "perm.method",
        "min.count", "miss.rate"),
        list(100, "obs.txt", "perm.txt", 2, 5, 0.20))
  if (family == 1) {
    op$family <- binomial()
    op$y.cont <- FALSE
  } else {
    op$family <- gaussian()
    op$y.cont <- TRUE
  }
  op$test <- 4
  op$test4.min.count <- op$min.count
  op$test4.missing.rate <- op$miss.rate

  ret <- permReg(snp.list, pheno.list, op=op)
  
  NULL

} # END: runPermutations

