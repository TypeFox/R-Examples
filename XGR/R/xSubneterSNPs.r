#' Function to identify a gene network from an input network given a list of seed SNPs together with the significance level (e.g. GWAS reported p-values)
#'
#' \code{xSubneterSNPs} is supposed to identify maximum-scoring gene subnetwork from an input graph with the node information on the significance (measured as p-values or fdr). To do so, it defines seed genes and their weights that take into account the distance to and the significance of input SNPs. It returns an object of class "igraph". 
#'
#' @param data a named input vector containing the sinificance level for nodes (dbSNP). For this named vector, the element names are dbSNP, the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for dbSNP, 2nd column for the significance level
#' @param include.LD additional SNPs in LD with Lead SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, LD SNPs will be included based on one or more of 26 populations and 5 super populations from 1000 Genomics Project data (phase 3). The population can be one of 5 super populations ("AFR", "AMR", "EAS", "EUR", "SAS"), or one of 26 populations ("ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"). Explanations for population code can be found at \url{http://www.1000genomes.org/faq/which-populations-are-part-your-study}
#' @param LD.r2 the LD r2 value. By default, it is 0.8, meaning that SNPs in LD (r2>=0.8) with input SNPs will be considered as LD SNPs. It can be any value from 0.8 to 1
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathways Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), and "STRING_medium" for interactions with medium confidence (confidence scores>=400). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addtion to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network
#' @param distance.max the maximum distance between genes and SNPs. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby SNPs per gene
#' @param seed.genes logical to indicate whether the identified network is restricted to seed genes (ie nearby genes that are located within defined distance window centred on lead or LD SNPs). By default, it sets to true
#' @param subnet.significance the given significance threshold. By default, it is set to NULL, meaning there is no constraint on nodes/genes. If given, those nodes/genes with p-values below this are considered significant and thus scored positively. Instead, those p-values above this given significance threshold are considered insigificant and thus scored negatively
#' @param subnet.size the desired number of nodes constrained to the resulting subnet. It is not nulll, a wide range of significance thresholds will be scanned to find the optimal significance threshold leading to the desired number of nodes in the resulting subnet. Notably, the given significance threshold will be overwritten by this option
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a subgraph with a maximum score, an object of class "igraph". It has ndoe attributes: signficance, score
#' @note The algorithm identifying a gene subnetwork that is likely modulated by input SNPs and/or their LD SNPs includes two major steps. The first step is to define and score nearby genes that are located within distance window of input and/or LD SNPs. The second step is to use \code{\link{xSubneterGenes}} for identifying a maximum-scoring gene subnetwork that contains as many highly scored genes as possible but a few lowly scored genes as linkers.
#' @export
#' @seealso \code{\link{xSubneterGenes}}
#' @include xSubneterSNPs.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' library(dnet)
#' library(GenomicRanges)
#'
#' # a) provide the seed SNPs with the weight info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' seeds.snps <- as.matrix(mcols(gr)[,c(1,3)])
#' 
#' # b) perform network analysis
#' # b1) find maximum-scoring subnet based on the given significance threshold
#' subnet <- xSubneterSNPs(data=seeds.snps, network="STRING_high", seed.genes=F, subnet.significance=0.01)
#' # b2) find maximum-scoring subnet with the desired node number=50
#' subnet <- xSubneterSNPs(data=data, network="STRING_high", subnet.size=50)
#'
#' # c) save subnet results to the files called 'subnet_edges.txt' and 'subnet_nodes.txt'
#' output <- igraph::get.data.frame(subnet, what="edges")
#' utils::write.table(output, file="subnet_edges.txt", sep="\t", row.names=FALSE)
#' output <- igraph::get.data.frame(subnet, what="vertices")
#' utils::write.table(output, file="subnet_nodes.txt", sep="\t", row.names=FALSE)
#'
#' # d) visualise the identified subnet
#' ## do visualisation with nodes colored according to the significance
#' xVisNet(g=subnet, pattern=-log10(as.numeric(V(subnet)$significance)), vertex.shape="sphere", colormap="wyr")
#' ## do visualisation with nodes colored according to transformed scores
#' xVisNet(g=subnet, pattern=V(subnet)$score, vertex.shape="sphere")
#' 
#' # e) visualise the identified subnet as a circos plot
#' library(RCircos)
#' xCircos(g=subnet, entity="Gene")
#' }

xSubneterSNPs <- function(data, include.LD=NA, LD.r2=0.8, network=c("STRING_highest","STRING_high","STRING_medium","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD"), network.customised=NULL, distance.max=200000, seed.genes=T, subnet.significance=5e-5, subnet.size=NULL, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0")
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    network <- match.arg(network)
    
    if(is.null(data)){
        stop("The input data must be not NULL.\n")
    }else{
		if (is.vector(data)){
			if(length(data)>1){
				# assume a vector
				if(is.null(names(data))){
					stop("The input data must have names with attached gene symbols.\n")
				}
			}else{
				# assume a file
				data <- utils::read.delim(file=data, header=F, row.names=NULL, stringsAsFactors=F)
			}
		}
		
		if (is.vector(data)){
			pval <- data
		}else if(is.matrix(data) | is.data.frame(data)){
			data <- as.matrix(data)
			data_list <- split(x=data[,2], f=as.character(data[,1]))
			res_list <- lapply(data_list, function(x){
				x <- as.numeric(x)
				x <- x[!is.na(x)]
				if(length(x)>0){
					min(x)
				}else{
					NULL
				}
			})
			pval <- unlist(res_list)
		}
		
		# force those zeros to be miminum of non-zeros
		#tmp <- as.numeric(format(.Machine)['double.xmin'])
		tmp <- min(pval[pval!=0])
		pval[pval < tmp] <- tmp
	}
	
	Lead_Sig <- data.frame(SNP=names(pval), Sig=pval, row.names=NULL, stringsAsFactors=F)
	leads <- Lead_Sig[,1]
	sigs <- Lead_Sig[,2]

	###########################
	## include additional SNPs that are in LD with input SNPs
	if(LD.r2>=0.8 & LD.r2<=1){
		default.include.LD <- c("ACB","AFR","AMR","ASW","BEB","CDX","CEU","CHB","CHS","CLM","EAS","ESN","EUR","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","SAS","STU","TSI","YRI")
		ind <- match(default.include.LD, include.LD)
		include.LD <- default.include.LD[!is.na(ind)]
	}
	
	if(length(include.LD) > 0){
		GWAS_LD <- xRDataLoader(RData.customised='GWAS_LD', RData.location=RData.location, verbose=verbose)
		res_list <- lapply(include.LD, function(x){
			data_ld <- ''
			eval(parse(text=paste("data_ld <- GWAS_LD$", x, sep="")))
			ind <- match(rownames(data_ld), leads)
			ind_lead <- which(!is.na(ind))
			ind_ld <- which(Matrix::colSums(data_ld[ind_lead,]>=LD.r2)>0)
		
			sLL <- data_ld[ind_lead, ind_ld]
			summ <- summary(sLL)
			res <- data.frame(Lead=rownames(sLL)[summ$i], LD=colnames(sLL)[summ$j], R2=summ$x, stringsAsFactors=F)
		})
		## get data frame (Lead LD R2)
		LLR <- do.call(rbind, res_list)
		
		###########################
		## also based on ImmunoBase
		if(1){
			ImmunoBase_LD <- xRDataLoader(RData.customised='ImmunoBase_LD', RData.location=RData.location, verbose=verbose)
			res_list <- lapply(include.LD, function(x){
				data_ld <- ''
				eval(parse(text=paste("data_ld <- ImmunoBase_LD$", x, sep="")))
				ind <- match(rownames(data_ld), leads)
				ind_lead <- which(!is.na(ind))
				ind_ld <- which(Matrix::colSums(data_ld[ind_lead,]>=LD.r2)>0)
		
				sLL <- data_ld[ind_lead, ind_ld]
				summ <- summary(sLL)
				res <- data.frame(Lead=rownames(sLL)[summ$i], LD=colnames(sLL)[summ$j], R2=summ$x, stringsAsFactors=F)
			})
			## get data frame (Lead LD R2)
			LLR_tmp <- do.call(rbind, res_list)
			LLR <- rbind(LLR, LLR_tmp)
		}
		###########################
				
		## get data frame (LD Sig)
		ld_list <- split(x=LLR[,-2], f=LLR[,2])
		res_list <- lapply(ld_list, function(x){
			ind <- match(x$Lead, leads)
			## power transformation of p-values X by R2, then keep the min
			min(sigs[ind] ^ x$R2)
		})
		vec <- unlist(res_list)
		LD_Sig <- data.frame(SNP=names(vec), Sig=vec, row.names=NULL, stringsAsFactors=F)
	
		## merge Lead and LD
		df <- rbind(Lead_Sig, as.matrix(LD_Sig))
		res_list <- split(x=df$Sig, f=df$SNP)
		res <- lapply(res_list, function(x){
			min(x)
		})
		vec <- unlist(res)
		SNP_Sig <- data.frame(SNP=names(vec), FDR=vec, row.names=NULL, stringsAsFactors=F)
	
	}else{
		SNP_Sig <- Lead_Sig
	}
	###########################
	
	pval <- as.numeric(SNP_Sig[,2])
	names(pval) <- SNP_Sig[,1]
	
	# transformed into scores according to log-likelihood ratio between the true positives and the false positivies
	scores <- log10(2) * dFDRscore(pval, fdr.threshold=NULL, scatter=F)
	scores[scores<0] <- 0
	seeds.snps <- scores
    
    ######################################################
    # Link to targets based on genomic distance
    ######################################################
    
  	## load positional information
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for SNPs (%s) ...", as.character(now)), appendLF=T)
	}
  	pos_SNP <- xRDataLoader(RData.customised="RegulomeDB_SNPs", RData.location=RData.location, verbose=verbose)
  	ind <- match(names(seeds.snps), names(pos_SNP))
  	ind <- ind[!is.na(ind)]
  	if(length(ind)){
  		gr_SNP <- pos_SNP[ind,]
  		
  		## append p-value weight
  		wS <- seeds.snps[names(gr_SNP)]
  		GenomicRanges::mcols(gr_SNP) <- data.frame(mcols(gr_SNP), wS=wS)
  		
		if(verbose){
			now <- Sys.time()
			message(sprintf("\t %d SNPs are used as seeds", length(gr_SNP)), appendLF=T)
		}
  	}
  	
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=T)
	}
  	gr_Gene <- xRDataLoader(RData.customised="UCSC_genes", RData.location=RData.location, verbose=verbose)
    
	# genes: get all UCSC genes within 500k away from variants
	maxgap <- distance.max
	minoverlap <- 1L # 1b overlaps
	subject <- gr_Gene
	query <- gr_SNP
	q2r <- as.matrix(suppressWarnings(GenomicRanges::findOverlaps(query=query, subject=subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
	
	list_gene <- split(x=q2r[,1], f=q2r[,2])
	ind_gene <- as.numeric(names(list_gene))
	res_list <- lapply(1:length(ind_gene), function(i){
		x <- subject[ind_gene[i],]
		y <- query[list_gene[[i]],]
		dists <- GenomicRanges::distance(x, y, select="all", ignore.strand=T)
		
		## weights according to distance away from lead SNPs
		#wD <- 1- dists/maxgap
		wD <- 10^(-1*dists/maxgap)
		## weights according to P-values
		wS <- mcols(y)$wS
		
		## seeds weights according to wD and WP
		res <- max(wD * wS)
		names(res) <- mcols(x)$Symbol
		res
	})
	dist.seeds.genes <- unlist(res_list)
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\t %d Genes are defined according to genomic distance", length(dist.seeds.genes)), appendLF=T)
	}
    
    ################################
	seeds.genes <- dist.seeds.genes	
	## take back to the p-value format
    pval <- 10^(-1*seeds.genes)
    ################################
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("\t\t minimum p-value: %1.2e; maximum p-value: %1.2e", min(pval), max(pval)), appendLF=T)
	}
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("xSubneterGenes is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
    subg <- xSubneterGenes(data=pval, network=network, network.customised=network.customised, seed.genes=seed.genes, subnet.significance=subnet.significance, subnet.size=subnet.size, verbose=verbose, RData.location=RData.location)

	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("xSubneterGenes has finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    return(subg)
}
