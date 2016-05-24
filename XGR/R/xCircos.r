#' Function to visualise a network as a circos plot
#'
#' \code{xCircos} is used to visualise a network as a circos plot. The network must be a 'igraph' object. 
#'
#' @param g an object of class "igraph". For example, it stores semantic similarity results with nodes for genes/SNPs and edges for pair-wise semantic similarity between them 
#' @param entity the entity of similarity analysis for which results are being plotted. It can be either "SNP" or "Gene"
#' @param top_num the top number of similarity edges to be plotted
#' @param ideogram logical to indicate whether chromosome banding is plotted
#' @param chr.exclude a character vector of chromosomes to exclude from the plot, e.g. c("chrX", "chrY"). By defautl, it is 'auto' meaning those chromosomes without data will be excluded. If NULL, no chromosome is excluded
#' @param entity.label.cex the font size of genes/SNPs labels. Default is 0.8
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' a circos plot with edge weights between input snps/genes represented by the colour of the links
#' @note none
#' @export
#' @import RCircos
#' @seealso \code{\link{xSocialiserGenes}}, \code{\link{xSocialiserSNPs}}
#' @include xCircos.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' library(RCircos)
#' library(GenomicRanges)
#' 
#' # provide genes and SNPs reported in AS GWAS studies
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#' 
#' # 1) SNP-based similarity analysis using GWAS Catalog traits (mapped to EF)
#' ## Get lead SNPs reported in AS GWAS
#' example.snps <- names(ImmunoBase$AS$variants)
#' SNP.g <- xSocialiserSNPs(example.snps, include.LD=NA)
#' # Circos plot of the EF-based SNP similarity network
#' #out.file <- "SNP_Circos.pdf"
#' #pdf(file=out.file, height=12, width=12, compress=TRUE)
#' xCircos(g=SNP.g, entity="SNP")
#' #dev.off()
#'
#' # 2) Gene-based similarity analysis using Disease Ontology (DO)
#' ## Get genes within 10kb away from AS GWAS lead SNPs
#' example.genes <- names(which(ImmunoBase$AS$genes_variants<=10000))
#' gene.g <- xSocialiserGenes(example.genes, ontology=c("DO")
#' # Circos plot of the DO-based gene similarity network
#' #out.file <- "Gene_Circos.pdf"
#' #pdf(file=out.file, height=12, width=12, compress=TRUE)
#' xCircos(g=gene.g, entity="Gene", chr.exclude="chrY")
#' #dev.off()
#' } 

xCircos <- function(g, entity=c("SNP","Gene"), top_num=50, ideogram=T, chr.exclude="auto", entity.label.cex=0.8, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0")
{
  
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
	entity <- match.arg(entity)
	
  	## Check input g
  	if (class(g) != "igraph") {
    	stop("The function must apply to a 'igraph' object.\n")
  	}
  
	## Convert from igraph into data.frame
  	df <- igraph::get.data.frame(g, what="edges")

  	if(is.null(top_num)){
    	top_num <- nrow(df)
  	}
  	if(top_num > nrow(df)){
    	top_num <- nrow(df)
  	}
  	top_num <- as.integer(top_num)
  	df <- df[1:top_num, ]
  
  	## load positional information
	if(verbose){
		now <- Sys.time()
		message(sprintf("Loading positional information for %s (%s) ...", entity, as.character(now)), appendLF=T)
	}
  	if(entity=="SNP") {
    	pos <- xRDataLoader(RData.customised="RegulomeDB_SNPs", verbose=verbose, RData.location=RData.location)
  	}else if(entity == "Gene") {
    	pos <- xRDataLoader(RData.customised="UCSC_genes", verbose=verbose, RData.location=RData.location)
  	}else{
    	stop("Please indicate whether your analysis entity was SNPs or genes.\n")
  	}
  
  	## Convert into format required for Circos plot
  	if(entity == "SNP"){
		allnames <- names(pos)
    }else if(entity == "Gene"){
    	allnames <- mcols(pos)$Symbol
    }
    A <- match(df$from, allnames)
	B <- match(df$to, allnames)
	#flag <- complete.cases(cbind(A, B))
	flag <- !is.na(A) & !is.na(B)
	AA <- A[flag]
	BB <- B[flag]
	input.data.A <- GenomicRanges::as.data.frame(pos[AA], row.names=NULL)
	input.data.B <- GenomicRanges::as.data.frame(pos[BB], row.names=NULL)
	input.data <- cbind.data.frame(input.data.A[, 1:3], input.data.B[, 1:3])
	if(is.null(df$weight)){
		input.data$similarity <- rep(1, sum(flag))
	}else{
		input.data$similarity <- as.numeric(as.character(df$weight[flag]))
	}
	label.data <- rbind(input.data.A[, 1:3], input.data.B[, 1:3])
	label.data$Name <- c(df$from[flag], df$to[flag])
  	
  	## decide on which chromosomes will be excluded
  	if(!is.null(chr.exclude)){
  		chr.exclude <- chr.exclude[!is.na(chr.exclude)]
		if(length(chr.exclude)==0){
			chr.exclude <- NULL
		}else if(sum(chr.exclude=='auto')>0){
			flag <- levels(label.data$seqnames) %in% as.character(unique(label.data$seqnames))
			chr.exclude <- levels(label.data$seqnames)[!flag]
		}
  	}
  	
  	## Load human chromosome ideogram
	if(verbose){
		now <- Sys.time()
		message(sprintf("Loading human chromosome banding information (hg19) (%s) ...", as.character(now)), appendLF=T)
	}
	
	#data(UCSC.HG19.Human.CytoBandIdeogram, package="RCircos")
	eval(parse(text="data(UCSC.HG19.Human.CytoBandIdeogram)"))
	
  	cyto.info <- ""
  	eval(parse(text=paste("cyto.info <- UCSC.HG19.Human.CytoBandIdeogram", sep="")))
  	if(ideogram==F) {
    	cyto.info$Stain <- rep("gpos100", nrow(cyto.info))
  	}
  	
  	## Set RCircos core components
	if(verbose){
		now <- Sys.time()
		message(sprintf("Initialising RCircos Core Components (%s) ...", as.character(now)), appendLF=T)
	}
  	num.inside <- 1
  	num.outside <- 1
  	RCircos.Set.Core.Components(cyto.info, chr.exclude, num.inside, num.outside)
  
  	## Reset parameters
  	params <- RCircos.Get.Plot.Parameters()
  	params$track.padding <- 0 # 0.02
  	params$track.height <- 0.05 # 0.1
  	
  	params$chr.ideog.pos <- 1
  	params$highlight.pos <- 1.1
  	params$chr.name.pos <- 1.1
  	params$plot.radius <- 0.9
  	params$track.out.start <- 1.2
  	params$highlight.width <- 0
  	
  	params$text.size <- entity.label.cex
  	RCircos.Reset.Plot.Parameters(params)
  
  	## Initialise graphic device, plot chromosome ideogram
	if(verbose){
		now <- Sys.time()
		message(sprintf("Plotting chromosome ideogram (%s) ...", as.character(now)), appendLF=T)
	}
  	RCircos.Set.Plot.Area()
  	RCircos.Chromosome.Ideogram.Plot()
  
  	## Plot link data coloured according to the similarity output
	if(verbose){
		now <- Sys.time()
		message(sprintf("Plotting link data (%s) ...", as.character(now)), appendLF=T)
	}
  	input.data$PlotColor <- colorRampPalette(c("gray", "red"))(10)[as.numeric(cut(input.data$similarity, breaks=seq(0, 1, 0.1)))]
  	input.data <- input.data[order(input.data$similarity, decreasing=F), ]
  	RCircos.Link.Plot(input.data, track.num=1, FALSE)

  	## Label SNPs/genes in outside track
	if(verbose){
		now <- Sys.time()
		message(sprintf("Adding SNP or gene names (%s) ...", as.character(now)), appendLF=T)
	}
  	name.col <- "Name"
  	side <- "out"
  	track.num <- 1
  	label.data <- label.data[!duplicated(label.data$Name), ]
  	if(verbose){
  		RCircos.Gene.Name.Plot(label.data, name.col, track.num, side)
  	}else{
  		suppressMessages(RCircos.Gene.Name.Plot(label.data, name.col, track.num, side))
  	}
  	
  	invisible()
}
