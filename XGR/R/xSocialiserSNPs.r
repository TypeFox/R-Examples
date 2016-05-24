#' Function to calculate pair-wise semantic similarity given a list of SNPs and the ontology in query
#'
#' \code{xSocialiserSNPs} is supposed to calculate pair-wise semantic similarity between a list of input SNPs and the ontology in query. It returns an object of class "igraph", a network representation of socialized SNPs. Now it supports analysis for SNPs using GWAS Catalog traits mapped to Experimental Factor Ontology. If required, additional SNPs that are in linkage disequilibrium (LD) with input SNPs are also be used for calculation. It first calculates semantic similarity between terms and then derives semantic similarity from term-term semantic similarity. Parallel computing is also supported for Linux or Mac operating systems.
#'
#' @param data an input vector. It contains a list of SNPs of interest
#' @param ontology the ontology supported currently. Now it is only "EF" for Experimental Factor Ontology (used to annotate GWAS Catalog SNPs). However, there are several subparts of this ontology to choose: 'EF_disease' for the subpart under the term 'disease' (EFO:0000408), 'EF_phenotype' for the subpart under the term 'phenotype' (EFO:0000651), 'EF_bp' for the subpart under the term 'biological process' (GO:0008150)
#' @param include.LD additional SNPs in LD with input SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, LD SNPs will be included based on one or more of 26 populations and 5 super populations from 1000 Genomics Project data (phase 3). The population can be one of 5 super populations ("AFR", "AMR", "EAS", "EUR", "SAS"), or one of 26 populations ("ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"). Explanations for population code can be found at \url{http://www.1000genomes.org/faq/which-populations-are-part-your-study}
#' @param LD.r2 the LD r2 value. By default, it is 0.8, meaning that SNPs in LD (r2>=0.8) with input SNPs will be considered as LD SNPs. It can be any value from 0.8 to 1
#' @param measure the measure used to derive semantic similarity between genes/SNPs from semantic similarity between terms. Take the semantic similartity between SNPs as an example. It can be "average" for average similarity between any two terms (one from SNP 1, the other from SNP 2), "max" for the maximum similarity between any two terms, "BM.average" for best-matching (BM) based average similarity (i.e. for each term of either SNP, first calculate maximum similarity to any term in the other SNP, then take average of maximum similarity; the final BM-based average similiary is the pre-calculated average between two SNPs in pair), "BM.max" for BM based maximum similarity (i.e. the same as "BM.average", but the final BM-based maximum similiary is the maximum of the pre-calculated average between two SNPs in pair), "BM.complete" for BM-based complete-linkage similarity (inspired by complete-linkage concept: the least of any maximum similarity between a term of one SNP and a term of the other SNP). When comparing BM-based similarity between SNPs, "BM.average" and "BM.max" are sensitive to the number of terms invovled; instead, "BM.complete" is much robust in this aspect. By default, it uses "BM.average"
#' @param method.term the method used to measure semantic similarity between terms. It can be "Resnik" for information content (IC) of most informative common ancestor (MICA) (see \url{http://dl.acm.org/citation.cfm?id=1625914}), "Lin" for 2*IC at MICA divided by the sum of IC at pairs of terms (see \url{https://www.cse.iitb.ac.in/~cs626-449/Papers/WordSimilarity/3.pdf}), "Schlicker" for weighted version of 'Lin' by the 1-prob(MICA) (see \url{http://www.ncbi.nlm.nih.gov/pubmed/16776819}), "Jiang" for 1 - difference between the sum of IC at pairs of terms and 2*IC at MICA (see \url{http://arxiv.org/pdf/cmp-lg/9709008.pdf}), "Pesquita" for graph information content similarity related to Tanimoto-Jacard index (ie. summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 (see \url{http://www.ncbi.nlm.nih.gov/pubmed/18460186}))
#' @param rescale logical to indicate whether the resulting values are rescaled to the range [0,1]. By default, it sets to true
#' @param force logical to indicate whether the only most specific terms (for each SNP) will be used. By default, it sets to true. It is always advisable to use this since it is computationally fast but without compromising accuracy (considering the fact that true-path-rule has been applied when running \code{\link{xDAGanno}})
#' @param fast logical to indicate whether a vectorised fast computation is used. By default, it sets to true. It is always advisable to use this vectorised fast computation; since the conventional computation is just used for understanding scripts
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to true
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' It returns an object of class "igraph", with nodes for input SNPs and edges for pair-wise semantic similarity between them. If no similarity is calculuated, it returns NULL.
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @importFrom Matrix colSums
#' @seealso \code{\link{xSocialiser}}
#' @include xSocialiserSNPs.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' 
#' # SNP-based similarity analysis using GWAS Catalog traits (mapped to EF)
#' # a) provide the input SNPs of interest (eg 8 randomly chosen SNPs)
#' anno <- xRDataLoader(RData='GWAS2EF')
#' allSNPs <- rownames(anno)
#' data <- sample(allSNPs,8)
#' data
#' 
#' # b) perform similarity analysis
#' sim <- xSocialiserSNPs(data=data)
#'
#' # b') optionally, enrichment analysis for input SNPs plus their LD SNPs
#' ## LD based on European population (EUR) with r2>=0.8
#' #sim <- xSocialiserSNPs(data=data, include.LD="EUR", LD.r2=0.8)
#' 
#' # c) save similarity results to the file called 'EF_similarity.txt'
#' output <- igraph::get.data.frame(sim, what="edges")
#' utils::write.table(output, file="EF_similarity.txt", sep="\t", row.names=FALSE)
#'
#' # d) visualise the SNP network
#' ## extract edge weight (with 2-digit precision)
#' x <- signif(as.numeric(E(sim)$weight), digits=2)
#' ## rescale into an interval [1,4] as edge width
#' edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
#' ## do visualisation
#' xVisNet(g=sim, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
#' }

xSocialiserSNPs <- function(data, ontology=c("EF","EF_disease","EF_phenotype", "EF_bp"), include.LD=NA, LD.r2=0.8, measure=c("BM.average","BM.max","BM.complete","average","max"), method.term=c("Resnik","Lin","Schlicker","Jiang","Pesquita"), rescale=TRUE, force=TRUE, fast=TRUE, parallel=TRUE, multicores=NULL, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=T, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    ontology <- match.arg(ontology)
    measure <- match.arg(measure)
    method.term <- match.arg(method.term)
    path.mode <- match.arg(path.mode)
    
    if (is.vector(data)){
        data <- unique(data)
    }else{
        stop("The input data must be a vector.\n")
    }
    
    if(!is.na(ontology)){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Load the ontology %s and its SNP annotations (%s) ...", ontology, as.character(now)), appendLF=T)
		}
		
		
		#########
		## load ontology information
		ig <- xRDataLoader(RData=paste('ig.EF', sep=''), RData.location=RData.location, verbose=verbose)
		
		if(ontology != 'EF'){
			if(ontology=='EF_disease'){
				node <- 'EFO:0000408'
			}else if(ontology=='EF_phenotype'){
				node <- 'EFO:0000651'
			}else if(ontology=='EF_bp'){
				node <- 'GO:0008150'
			}
            neighs.out <- igraph::neighborhood(ig, order=vcount(ig), nodes=node, mode="out")
            nodeInduced <- V(ig)[unique(unlist(neighs.out))]$name
            g <- igraph::induced.subgraph(ig, vids=nodeInduced)
        }else{
        	g <- ig
        }
        
		#########
		## load annotation information
		anno <- xRDataLoader(RData=paste('GWAS2EF', sep=''), RData.location=RData.location, verbose=verbose)
			
		#########
		## include additional SNPs that are in LD with input SNPs
		if(LD.r2>=0.8 & LD.r2<=1){
			default.include.LD <- c("ACB","AFR","AMR","ASW","BEB","CDX","CEU","CHB","CHS","CLM","EAS","ESN","EUR","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","SAS","STU","TSI","YRI")
			
			ind <- match(default.include.LD, include.LD)
			include.LD <- default.include.LD[!is.na(ind)]
			
			if(length(include.LD) > 0){
				GWAS_LD <- xRDataLoader(RData='GWAS_LD', RData.location=RData.location, verbose=verbose)
				
				ld.list <- lapply(include.LD, function(x){
					data_ld <- ''
					eval(parse(text=paste("data_ld <- GWAS_LD$", x, sep="")))
					## bugs
					#ind <- match(rownames(data_ld), data)
					#ld <- colnames(data_ld)[which(Matrix::colSums(data_ld[!is.na(ind),]>=LD.r2)>0)]
					
					ind <- match(colnames(data_ld), data)
					ld <- rownames(data_ld)[which(Matrix::rowSums(data_ld[,!is.na(ind)]>=LD.r2)>0)]
					
					ld
				})
				
				data <- base::union(unlist(ld.list), data)
			}
			
		}
	
	}else{
		stop("There is no input for the ontology.\n")
	}
	
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xSocialiser' is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
    res <- xSocialiser(data=data, annotation=anno, g=g, measure=measure, method.term=method.term, rescale=rescale, force=force, fast=fast, parallel=parallel, multicores=multicores, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose)
    
    ## the resulting graph has gene symbols (instead of Entrez GeneIDs) as nodes
    if(!is.null(res)){
    	sim_ig <- res
    	## sort (by weight) and order (from and to)
		tEdges <- igraph::get.data.frame(sim_ig, what="edges")
		sEdges <- tEdges[with(tEdges,order(-weight,from,to)), ]
		relations <- apply(sEdges, 1, function(x){
			if(x[1] > x[2]){
				return(cbind(x[2],x[1],x[3]))
			}else{
				x
			}
		})
		relations <- t(relations)
		colnames(relations) <- c("from","to","weight")
		res <- igraph::graph.data.frame(d=relations, directed=F)
    }
    

	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xSocialiser' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(res)
}
