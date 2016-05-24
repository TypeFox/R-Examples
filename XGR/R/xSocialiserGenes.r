#' Function to calculate pair-wise semantic similarity given a list of genes and the ontology in query
#'
#' \code{xSocialiserGenes} is supposed to calculate pair-wise semantic similarity between a list of input genes and the ontology in query. It returns an object of class "igraph", a network representation of socialized genes. Now it supports enrichment analysis using a wide variety of ontologies such as Gene Ontology and Phenotype Ontologies. It first calculates semantic similarity between terms and then derives semantic similarity from term-term semantic similarity. Parallel computing is also supported for Linux or Mac operating systems.
#'
#' @param data an input vector. It contains a list of Gene Symbols of interest
#' @param ontology the ontology supported currently. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPCM" for Human Phenotype Clinical Modifier, "HPMA" for Human Phenotype Mortality Aging, "MP" for Mammalian Phenotype
#' @param measure the measure used to derive semantic similarity between genes/SNPs from semantic similarity between terms. Take the semantic similartity between SNPs as an example. It can be "average" for average similarity between any two terms (one from SNP 1, the other from SNP 2), "max" for the maximum similarity between any two terms, "BM.average" for best-matching (BM) based average similarity (i.e. for each term of either SNP, first calculate maximum similarity to any term in the other SNP, then take average of maximum similarity; the final BM-based average similiary is the pre-calculated average between two SNPs in pair), "BM.max" for BM based maximum similarity (i.e. the same as "BM.average", but the final BM-based maximum similiary is the maximum of the pre-calculated average between two SNPs in pair), "BM.complete" for BM-based complete-linkage similarity (inspired by complete-linkage concept: the least of any maximum similarity between a term of one SNP and a term of the other SNP). When comparing BM-based similarity between SNPs, "BM.average" and "BM.max" are sensitive to the number of terms invovled; instead, "BM.complete" is much robust in this aspect. By default, it uses "BM.average"
#' @param method.term the method used to measure semantic similarity between terms. It can be "Resnik" for information content (IC) of most informative common ancestor (MICA) (see \url{http://dl.acm.org/citation.cfm?id=1625914}), "Lin" for 2*IC at MICA divided by the sum of IC at pairs of terms (see \url{https://www.cse.iitb.ac.in/~cs626-449/Papers/WordSimilarity/3.pdf}), "Schlicker" for weighted version of 'Lin' by the 1-prob(MICA) (see \url{http://www.ncbi.nlm.nih.gov/pubmed/16776819}), "Jiang" for 1 - difference between the sum of IC at pairs of terms and 2*IC at MICA (see \url{http://arxiv.org/pdf/cmp-lg/9709008.pdf}), "Pesquita" for graph information content similarity related to Tanimoto-Jacard index (ie. summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 (see \url{http://www.ncbi.nlm.nih.gov/pubmed/18460186}))
#' @param rescale logical to indicate whether the resulting values are rescaled to the range [0,1]. By default, it sets to true
#' @param force logical to indicate whether the only most specific terms (for each SNP) will be used. By default, it sets to true. It is always advisable to use this since it is computationally fast but without compromising accuracy (considering the fact that true-path-rule has been applied when running \code{\link{xDAGanno}})
#' @param fast logical to indicate whether a vectorised fast computation is used. By default, it sets to true. It is always advisable to use this vectorised fast computation; since the conventional computation is just used for understanding scripts
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to true
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' It returns an object of class "igraph", with nodes for input genes and edges for pair-wise semantic similarity between them. If no similarity is calculuated, it returns NULL.
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @seealso \code{\link{xSocialiser}}
#' @include xSocialiserGenes.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' 
#' # Gene-based similarity analysis using Mammalian Phenotype Ontology (MP)
#' # a) provide the input Genes of interest (eg 100 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg')
#' data <- as.character(sample(org.Hs.eg$gene_info$Symbol, 100))
#' data
#' 
#' # b) perform similarity analysis
#' sim <- xSocialiserGenes(data=data, ontology="MP")
#'
#' # c) save similarity results to the file called 'MP_similarity.txt'
#' output <- igraph::get.data.frame(sim, what="edges")
#' utils::write.table(output, file="MP_similarity.txt", sep="\t", row.names=FALSE)
#'
#' # d) visualise the gene network
#' ## extract edge weight (with 2-digit precision)
#' x <- signif(as.numeric(E(sim)$weight), digits=2)
#' ## rescale into an interval [1,4] as edge width
#' edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
#' ## do visualisation
#' xVisNet(g=sim, vertex.shape="sphere", edge.width=edge.width, edge.label=x, edge.label.cex=0.7)
#' }

xSocialiserGenes <- function(data, ontology=c("GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPCM","HPMA","MP"), measure=c("BM.average","BM.max","BM.complete","average","max"), method.term=c("Resnik","Lin","Schlicker","Jiang","Pesquita"), rescale=TRUE, force=TRUE, fast=TRUE, parallel=TRUE, multicores=NULL, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=T, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0")
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
			message(sprintf("Load the ontology %s and its gene annotations (%s) ...", ontology, as.character(now)), appendLF=T)
		}
                                  
		#########
		## load GS information
		GS <- xRDataLoader(RData=paste('org.Hs.eg', ontology, sep=''), RData.location=RData.location, verbose=verbose)
		
		#########
		## get annotation information
		anno <- GS$gs

		#########
		## get ontology information
		## check the eligibility for the ontology
		all.ontologies <- c("GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPCM","HPMA","MP")
		flag_ontology <- ontology %in% all.ontologies
    	
    	if(flag_ontology){
			g <- xRDataLoader(RData=paste('ig.', ontology, sep=''), RData.location=RData.location, verbose=verbose)
		}
	
	}else{
		stop("There is no input for the ontology.\n")
	}   
    
    #############################################################################################
    # A function converting from symbol to entrezgene
    symbol2entrezgene <- function(Symbol, check.symbol.identity, allGeneID, allSymbol, allSynonyms, verbose){
    
        ## correct for those symbols being shown as DATE format
        if(1){
            ## for those starting with 'Mar' in a excel-input date format
            a <- Symbol
            flag <- grep("-Mar$", a, ignore.case=T, perl=T, value=F)
            if(length(flag)>=1){
                b <- a[flag]
                c <- sub("-Mar$", "", b, ignore.case=T, perl=T)
                d <- sub("^0", "", c, ignore.case=T, perl=T)
                e <- sapply(d, function(x) paste(c("March",x), collapse=""))
                a[flag] <- e
                Symbol <- a
            }

            ## for those starting with 'Sep' in a excel-input date format
            a <- Symbol
            flag <- grep("-Sep$", a, ignore.case=T, perl=T, value=F)
            if(length(flag)>=1){
                b <- a[flag]
                c <- sub("-Sep$", "", b, ignore.case=T, perl=T)
                d <- sub("^0", "", c, ignore.case=T, perl=T)
                e <- sapply(d, function(x) paste(c("Sept",x), collapse=""))
                a[flag] <- e
                Symbol <- a
            }
        }
        
        ## case-insensitive
        match_flag <- match(tolower(Symbol),tolower(allSymbol))
        
        ## match vis Synonyms for those unmatchable by official gene symbols
        if(check.symbol.identity){
            ## match Synonyms (if not found via Symbol)
            na_flag <- is.na(match_flag)
            a <- Symbol[na_flag]

            ###
            tmp_flag <- is.na(match(tolower(allSymbol), tolower(Symbol)))
            tmp_Synonyms <- allSynonyms[tmp_flag]
            Orig.index <- seq(1,length(allSynonyms))
            Orig.index <- Orig.index[tmp_flag]
            ###

            b <- sapply(1:length(a), function(x){
                tmp_pattern1 <- paste("^",a[x],"\\|", sep="")
                tmp_pattern2 <- paste("\\|",a[x],"\\|", sep="")
                tmp_pattern3 <- paste("\\|",a[x],"$", sep="")
                tmp_pattern <- paste(tmp_pattern1,"|",tmp_pattern2,"|",tmp_pattern3, sep="")
                tmp_result <- grep(tmp_pattern, tmp_Synonyms, ignore.case=T, perl=T, value=F)
                ifelse(length(tmp_result)==1, Orig.index[tmp_result[1]], NA)
            })
            match_flag[na_flag] <- b
            
            if(verbose){
                now <- Sys.time()
                message(sprintf("\tAmong %d symbols of input data, there are %d mappable via official gene symbols, %d mappable via gene alias but %d left unmappable", length(Symbol), (length(Symbol)-length(a)), sum(!is.na(b)), sum(is.na(b))), appendLF=T)
            }
        }else{
            if(verbose){
                now <- Sys.time()
                message(sprintf("\tAmong %d symbols of input data, there are %d mappable via official gene symbols but %d left unmappable", length(Symbol), (sum(!is.na(match_flag))), (sum(is.na(match_flag)))), appendLF=T)
            }
        
        }
        
        ## convert into GeneID
        GeneID <- allGeneID[match_flag]
        
        return(GeneID)
    }
    #############################################################################################
    
    ## convert gene symbol to entrz gene for both input data of interest and the input background (if given)
    if(verbose){
		now <- Sys.time()
		message(sprintf("Do gene mapping from Symbols to EntrezIDs for (%s) ...", as.character(now)), appendLF=T)
	}
    
    ## load Enterz Gene information
	EG <- xRDataLoader(RData=paste('org.Hs.eg', sep=''), RData.location=RData.location, verbose=verbose)	
    
	allGeneID <- EG$gene_info$GeneID
	allSymbol <- as.vector(EG$gene_info$Symbol)
	allSynonyms <- as.vector(EG$gene_info$Synonyms)
    data <- symbol2entrezgene(Symbol=data, check.symbol.identity=F, allGeneID=allGeneID, allSymbol=allSymbol, allSynonyms=allSynonyms, verbose=verbose)
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xSocialiser' is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    
    res <- xSocialiser(data=data, annotation=anno, g=g, measure=measure, method.term=method.term, rescale=rescale, force=force, fast=fast, parallel=parallel, multicores=multicores, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose)
    
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xSocialiser' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ## the resulting graph has gene symbols (instead of Entrez GeneIDs) as nodes
    if(!is.null(res)){
    	sim_ig <- res
    	ind <- match(V(sim_ig)$name, allGeneID)
    	V(sim_ig)$name <- allSymbol[ind]
    	#res <- sim_ig
    	
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
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(res)
}
