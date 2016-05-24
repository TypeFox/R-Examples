##########
#map2class
##########


map2class <- function(coords, SNPs, extend.boundary = 0) {
  
  # Make sure coords object has correct columsn and format
  if(sum(!c("chr", "start", "stop", "class") %in% colnames(coords))) stop("coords object must include column names: 'chr', 'start', 'stop', and 'class'.")
  coords$class <- as.character(coords$class)
  
  # Make sure SNPs object has correct columsn and format
  if(sum(!c("SNP", "chr", "position") %in% colnames(SNPs))) stop("SNPs object must include column names: 'SNP', 'chr', 'position'.")
  SNPs$SNP <- as.character(SNPs$SNP)
  
  map2chr <- function(chrTemp){
    
    chr = position = NULL
    
    print(paste("Mapping SNPs to classes on chromosome ", chrTemp, ".", sep = ""))
    
    # Filter SNP and coordinate data for chromosome
    SNPchr <- filter(SNPs, chr == chrTemp)
    coordsChr <- filter(coords, chr == chrTemp)
    
    map2classTemp <- function(classTemp){
      
      # Filter coordinate data for class
      coordsSub <- filter(coordsChr, class == classTemp)
      
      # Extend coordinate boundaries
      coordsSub$start <- coordsSub$start - extend.boundary
      coordsSub$stop <- coordsSub$stop + extend.boundary
      
      # Filter SNPs data for class position
      SNPsub <- filter(SNPchr, position >= coordsSub$start)
      SNPsub <- filter(SNPsub, position <= coordsSub$stop)
      
      if(nrow(SNPsub)){
        return(data.frame(SNPsub, class = classTemp, stringsAsFactors = FALSE))
      } else {NULL}}
    
    # Run on one class at a time
    classes <- coordsChr$class
    
    classOut <- lapply(as.list(classes), map2classTemp)
    classOut1 <- do.call(rbind, classOut)
    return(classOut1)}
  
  # Run on one chromosome at a time
  chrTemp <- unique(coords$chr)
  chrTemp <- chrTemp[chrTemp%in%1:22 & chrTemp%in%SNPs$chr]
  
  chrOut <- lapply(as.list(chrTemp), map2chr)
  chrOut1 <- do.call(rbind, chrOut)
  return(chrOut1)}

#################
#GenCAT_manhattan
#################

GenCAT_manhattan <- function(GenCATout, sigThresh = NULL, highlightPosi = FALSE, labelPosi = FALSE, sepChr = 800000, plotTitle = "Manhattan Plot of GenCAT Results"){
  
  chr = position = pos = plotpos = negLogp = NULL
  
  # Check that input comes from a GenCAToutput
  if(class(GenCATout) != "GenCATtest" ) stop("GenCATout must be of the class, 'GenCATtest'.")
  
  # Create positions for each class from output based on median position of SNPs
  loc <- as.data.frame(GenCATout$Used %>% group_by(chr, class) %>% summarise(position = median(position)), stringsAsFactors = FALSE)
  
  # Subset GenCAT output for class and p-values and merge with positions
  GenCAT_data <- select_(GenCATout$GenCAT, ~one_of(c("class", "CsumP")))
  manhattanInput <- merge(loc, GenCAT_data)
  colnames(manhattanInput) <- c("class", "chr", "pos", "p")
  
  # Find -log10 of p-value for plotting
  manhattanInput$negLogp <- -log10(manhattanInput$p)
  
  # Order by chromsome and position
  manhattanInput <- arrange(manhattanInput, chr, pos)
  manhattanInput$plotpos <- manhattanInput$pos
  
  # Find and order chromosoms used in output
  chrs <- unique(manhattanInput$chr)
  chrs <- sort(chrs)
  # Create plot positions based on position and chromsome number
  for(i in 1:(length(chrs)-1)){
    maxChr<-max(manhattanInput$plotpos[manhattanInput$chr == chrs[i]]) + sepChr
    manhattanInput$plotpos[manhattanInput$chr == (chrs[i+1])]<-manhattanInput$plotpos[manhattanInput$chr == (chrs[i+1])] + maxChr}
  
  # Find the midpoint of the range of class positions on chromsome
  mid<-function(x) min(x) + (max(x)-min(x))/2
  
  # Find chromsome label position
  lab <- as.data.frame(manhattanInput %>% group_by(chr) %>% summarise(labpos1 = mid(plotpos), stringsAsFactors = FALSE))
  
  # Create manhattan plot
  p <- ggplot(manhattanInput) +
    
    # Add points for each class
    geom_point(aes(x = plotpos, y = negLogp, colour = as.factor(chr))) +
    
    # Add chromosome labels
    scale_x_continuous(labels = as.character(chrs), breaks = lab$labpos) +
    
    # Set chromosome colors to be black and grey
    scale_color_manual(values = rep(c('black', 'grey'), 12)) + 
    
    # Set plot background to be white
    theme_bw() +
    
    # Edit grid lines to be nice and clean
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = .4, color = "grey" ), panel.grid.minor = element_blank()) + 
    
    # No legend
    theme(legend.position = 'none') + 
    
    # Add axis labels and title
    xlab('Chromosome') + ylab('-log10(P-value)') + ggtitle(plotTitle)
  
  # If significance threshold is supplied add vertical line at this value
  if(!is.null(sigThresh)) {
    p <- p + geom_hline(yintercept = -log10(sigThresh), linetype = 1, col = 'black', lwd = 0.75)
    posiSub1 <- filter(manhattanInput, p < sigThresh)
    
    # Highlight significant classes
    if(highlightPosi){
      p <- p + geom_point(data = posiSub1, aes(x = plotpos, y = negLogp), colour = "blue")}
    # Label significant classes
    if(labelPosi){
      p <- p + geom_text(data = posiSub1, aes(x = plotpos, y = negLogp, label = class), size = 3, hjust = -0.08, vjust = 0)
    }
  }
  p}

################
#GenCAT Function
################

GenCAT<-function(SNPdata, genoData, snpInfo, pcCutoff = 0.95, workers = getOption("mc.cores",2L)){
  
  # Check that the required columns are included from GWAS data
  if(sum(!c("SNP", "effect_allele", "other_allele", "testStat", "chr", "class") %in% colnames(SNPdata))) stop("SNPdata must include column names: 'SNP', 'effect_allele', 'other_allele', 'testStat', 'chr', 'class'.")
  
  # Subset columns and format
  SNPdata<-select_(SNPdata, ~one_of(c("SNP", "effect_allele", "other_allele", "testStat", "chr", "class")))
  for(i in c("SNP", "effect_allele", "other_allele", 'class')) SNPdata[,i] <- as.character(SNPdata[,i])
  
  # If genotype data is supplied chec that snp information is also supplied
  # and in the correct format
  if(class(genoData) != "SnpMatrix") stop ("genoData must of of class SnpMatrix.")
  if(is.null(snpInfo)) stop("snpInfo must be supplied with genoData")
  if(sum(!c("chr", "SNP", "position", "A1", "A2") %in% colnames(snpInfo))) stop("snpInfo must include column names: 'chr', 'SNP', 'position', 'A1', 'A2'.")
  snpInfo <- select_(snpInfo, ~one_of(c("chr", "SNP", "position", "A1", "A2")))
  for(i in c("SNP", "A1", "A2")) snpInfo[,i] <- as.character(snpInfo[,i])
  snpInfo <- filter(snpInfo, snpInfo$SNP %in% colnames(genoData))
  
  # Run analysis for each chromosome
  chrs <- sort(unique(SNPdata$chr))
  
  GenCATchr <- function(chrTemp){
    
    chr = SNP = classTemp = effect_allele = other_allele = A1 = A2 = NULL
    
    # Create and report type of workers for parallel implementation
    cl <- makeCluster(workers)
    registerDoParallel(cl)
    
    # Subset dataframe for chromosome
    framChr <- filter(SNPdata, chr == chrTemp)
    
    
    # Subset for SNPs in chromosome from data
    chrSupport <- filter(snpInfo, chr == chrTemp)
    chrSupport <- filter(chrSupport, SNP %in% framChr$SNP)
    chrSub <- genoData[, colnames(genoData) %in% chrSupport$SNP]
    chrSub <- chrSub[, colnames(chrSub) %in% framChr$SNP]
    chrSub <- as(chrSub, "numeric")
    
    # Run GenCAT on classes of a particular chromosome
    classes <- unique(framChr$class)
    
    # Report which chromosome is being processed
    print(paste("Running GenCAT on ", length(classes), " classes on chromosome ", chrTemp, ".", sep = ""))
    
    pcCutoff <- 1-pcCutoff
    cl<-cl
    
    chrOut <- tryCatch(foreach(classTemp = classes, .combine = 'rbind', .packages = c("dplyr")) %dopar% {
      
      # Subset chr data.frame for gene
      framClass <- filter(framChr, class == classTemp)
      framClass <- select_(framClass, ~one_of(c("SNP", "effect_allele", "other_allele", "testStat", "class")))
      
      # Subset gene data.frame for SNPs in reference set
      geneSub_pre <- merge(framClass, chrSupport, by.x = "SNP", by.y = "SNP")
      
      # Check for duplicate SNPS
      if(sum(duplicated(geneSub_pre$SNP)>0)) stop(paste("There is a duplicate SNP in ", classTemp,".", sep=""))
      
      # Subset gene data.frame for SNPs that have matching alleles to reference set
      geneSub <- filter(geneSub_pre, effect_allele == A1 | other_allele == A1)
      geneSub <- filter(geneSub, effect_allele == A2 | other_allele == A2)
      
      #Find SNPs used in GenCAT
      Used <- as.matrix(geneSub)
      
      # Make sure that reference alleles match (change sign if they don't)
      geneSub$testStat[geneSub$effect_allele!= geneSub$A1] <- geneSub$testStat[geneSub$effect_allele!= geneSub$A1]*-1
      
      # Find removed SNP info
      notFound <- as.matrix(filter(framClass, !SNP %in% geneSub_pre$SNP))
      
      # Find unMatched SNP info
      unMatched <- filter(geneSub_pre, !SNP %in% geneSub$SNP)
      unMatched <- as.matrix(unMatched)
      
      
      if(nrow(geneSub)>0){
        
        # Read in reference SNP genotypes
        genoSub <- chrSub[,colnames(chrSub)%in%geneSub$SNP]
        
        # Find correlation and eigen matrix for genotypes
        cors <- cor(as.matrix(genoSub), use = "pairwise.complete.obs")
        
        # Check that all pairwise correlations are present
        if(sum(is.na(cors))>0) stop(paste("Pairwise correlation on class,", classTemp, ", found to be NA. Try filtering genotype data on more rigid minor allele frequency and/or call rate threshold."))
        
        e <- eigen(cors)
        e$values[e$values<0]<-0
        
        # Find m (The number of ordered eigenvalues that account for the specific cumulative proportion)
        eval <- e$values[rev(cumsum(rev(e$values)))/sum(e$values)>pcCutoff]
        m <- length(eval)
        
        # Match order of test statistics to correlation matrix
        geneSub <- geneSub[order(match(geneSub$SNP, colnames(genoSub))),]
        
        # Perform eigen decomposition
        H <- diag(1/sqrt(e$values[1:m]),m)%*%t(e$vectors[, (1:m)])
        
        # Obtain GenCAT test statistic and p-value
        test.new <- H%*%geneSub$testStat
        CsumT <- t(test.new)%*%test.new
        CsumP <- pchisq(CsumT, m, lower.tail = FALSE)
        
        # Create data frame of transformed test statistics
        transStats<-as.matrix(data.frame(class = classTemp, transStat = test.new))
        
        return(list(as.matrix(data.frame(class = classTemp, chr = chrTemp, n_SNPs = nrow(geneSub), n_Obs = m, CsumT = CsumT, CsumP = CsumP, stringsAsFactors = FALSE)), Used, notFound, unMatched, transStats))
        
      } else {return(list(NULL, Used, notFound, unMatched, NULL))}}
      , finally = stopCluster(cl))
  }
  
  
  GenCATout <- lapply(as.list(chrs), GenCATchr)
  GenCATout <- do.call(rbind, GenCATout)
  
  # Write out and format desired objects
  GenCATmain <- as.data.frame(do.call(rbind, GenCATout[,1]), stringsAsFactors = FALSE)
  for(i in c("chr", "n_SNPs", "n_Obs", "CsumT", "CsumP")) GenCATmain[,i] <- as.numeric(GenCATmain[,i])
  
  Used <- as.data.frame(do.call(rbind, GenCATout[,2]), stringsAsFactors = FALSE)
  for(i in c("testStat", "chr", "position")) Used[,i] <- as.numeric(Used[,i])
  
  notFound <- as.data.frame(do.call(rbind, GenCATout[,3]), stringsAsFactors = FALSE)
  notFound$testStat <- as.numeric(notFound$testStat)
  
  unMatched <- as.data.frame(do.call(rbind, GenCATout[,4]), stringsAsFactors = FALSE)
  for(i in c("testStat", "chr", "position")) unMatched[,i] <- as.numeric(unMatched[,i])
  
  TransStats <- as.data.frame(do.call(rbind, GenCATout[,5]), stringsAsFactors = FALSE)
  TransStats$transStat <- as.numeric(TransStats$transStat)
  
  # Specify format of output
  GenCATfinal <- list(GenCATmain, Used, notFound, unMatched, TransStats)
  names(GenCATfinal)<-c("GenCAT", "Used", "notFound", "unMatched", "TransStats")
  class(GenCATfinal) <- "GenCATtest"
  return(GenCATfinal)}

