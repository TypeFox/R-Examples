#' Function to conduct enrichment analysis given a list of genes and the ontology in query
#'
#' \code{xEnricherGenes} is supposed to conduct enrichment analysis given the input data and the ontology in query. It returns an object of class "eTerm". Enrichment analysis is based on either Fisher's exact test or Hypergeometric test. The test can respect the hierarchy of the ontology. Now it supports enrichment analysis using a wide variety of ontologies such as Gene Ontology and Phenotype Ontologies.
#'
#' @param data an input vector. It contains a list of Gene Symbols of interest
#' @param background a background vector. It contains a list of Gene Symbols as the test background. If NULL, by default all annotatable are used as background
#' @param ontology the ontology supported currently. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "PS" for phylostratific age information, "PS2" for the collapsed PS version (inferred ancestors being collapsed into one with the known taxonomy information), "SF" for domain superfamily assignments, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPCM" for Human Phenotype Clinical Modifier, "HPMA" for Human Phenotype Mortality Aging, "MP" for Mammalian Phenotype, and Drug-Gene Interaction database (DGIdb) for drugable categories, and the molecular signatures database (Msigdb, including "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7")
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 2000
#' @param min.overlap the minimum number of overlaps. Only those terms with members that overlap with input data at least min.overlap (3 by default) will be processed
#' @param which.distance which terms with the distance away from the ontology root (if any) is used to restrict terms in consideration. By default, it sets to 'NULL' to consider all distances
#' @param test the statistic test used. It can be "fisher" for using fisher's exact test, "hypergeo" for using hypergeometric test, or "binomial" for using binomial test. Fisher's exact test is to test the independence between gene group (genes belonging to a group or not) and gene annotation (genes annotated by a term or not), and thus compare sampling to the left part of background (after sampling without replacement). Hypergeometric test is to sample at random (without replacement) from the background containing annotated and non-annotated genes, and thus compare sampling to background. Unlike hypergeometric test, binomial test is to sample at random (with replacement) from the background with the constant probability. In terms of the ease of finding the significance, they are in order: hypergeometric test > binomial test > fisher's exact test. In other words, in terms of the calculated p-value, hypergeometric test < binomial test < fisher's exact test
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param ontology.algorithm the algorithm used to account for the hierarchy of the ontology. It can be one of "none", "pc", "elim" and "lea". For details, please see 'Note' below
#' @param elim.pvalue the parameter only used when "ontology.algorithm" is "elim". It is used to control how to declare a signficantly enriched term (and subsequently all genes in this term are eliminated from all its ancestors)
#' @param lea.depth the parameter only used when "ontology.algorithm" is "lea". It is used to control how many maximum depth is used to consider the children of a term (and subsequently all genes in these children term are eliminated from the use for the recalculation of the signifance at this term)
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to false
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return 
#' an object of class "eTerm", a list with following components:
#' \itemize{
#'  \item{\code{term_info}: a matrix of nTerm X 4 containing snp/gene set information, where nTerm is the number of terms, and the 4 columns are "id" (i.e. "Term ID"), "name" (i.e. "Term Name"), "namespace" and "distance"}
#'  \item{\code{annotation}: a list of terms containing annotations, each term storing its annotations. Always, terms are identified by "id"}
#'  \item{\code{data}: a vector containing input data in consideration. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{background}: a vector containing the background data. It is not always the same as the input data as only those mappable are retained}
#'  \item{\code{overlap}: a list of overlapped snp/gene sets, each storing snps overlapped between a snp/gene set and the given input data (i.e. the snps of interest). Always, gene sets are identified by "id"}
#'  \item{\code{zscore}: a vector containing z-scores}
#'  \item{\code{pvalue}: a vector containing p-values}
#'  \item{\code{adjp}: a vector containing adjusted p-values. It is the p value but after being adjusted for multiple comparisons}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note The interpretation of the algorithms used to account for the hierarchy of the ontology is:
#' \itemize{
#' \item{"none": does not consider the ontology hierarchy at all.}
#' \item{"lea": computers the significance of a term in terms of the significance of its children at the maximum depth (e.g. 2). Precisely, once snps are already annotated to any children terms with a more signficance than itself, then all these snps are eliminated from the use for the recalculation of the signifance at that term. The final p-values takes the maximum of the original p-value and the recalculated p-value.}
#' \item{"elim": computers the significance of a term in terms of the significance of its all children. Precisely, once snps are already annotated to a signficantly enriched term under the cutoff of e.g. pvalue<1e-2, all these snps are eliminated from the ancestors of that term).}
#' \item{"pc": requires the significance of a term not only using the whole snps as background but also using snps annotated to all its direct parents/ancestors as background. The final p-value takes the maximum of both p-values in these two calculations.}
#' \item{"Notes": the order of the number of significant terms is: "none" > "lea" > "elim" > "pc".}
#' }
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xEnricher}}
#' @include xEnricherGenes.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(XGR)
#' library(igraph)
#' 
#' # Gene-based enrichment analysis using Mammalian Phenotype Ontology (MP)
#' # a) provide the input Genes of interest (eg 100 randomly chosen human genes)
#' ## load human genes
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg')
#' data <- as.character(sample(org.Hs.eg$gene_info$Symbol, 100))
#' data
#' 
#' # optionally, provide the test background (if not provided, all human genes)
#' #background <- as.character(org.Hs.eg$gene_info$Symbol)
#' 
#' # b) perform enrichment analysis
#' eTerm <- xEnricherGenes(data=data, ontology="MP")
#'
#' # c) view enrichment results for the top significant terms
#' xEnrichViewer(eTerm)
#'
#' # d) save enrichment results to the file called 'MP_enrichments.txt'
#' res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
#' output <- data.frame(term=rownames(res), res)
#' utils::write.table(output, file="MP_enrichments.txt", sep="\t", row.names=FALSE)
#'
#' # e) visualise the top 10 significant terms in the ontology hierarchy
#' ## load ig.MP (an object of class "igraph" storing as a directed graph)
#' g <- xRDataLoader(RData='ig.MP')
#' g
#' nodes_query <- names(sort(eTerm$adjp)[1:10])
#' nodes.highlight <- rep("red", length(nodes_query))
#' names(nodes.highlight) <- nodes_query
#' subg <- dnet::dDAGinduce(g, nodes_query)
#' # color-code terms according to the adjust p-values (taking the form of 10-based negative logarithm)
#' dnet::visDAG(g=subg, data=-1*log10(eTerm$adjp[V(subg)$name]), node.info="both", zlim=c(0,2), node.attrs=list(color=nodes.highlight))
#' # color-code terms according to the z-scores
#' dnet::visDAG(g=subg, data=eTerm$zscore[V(subg)$name], node.info="both", colormap="darkblue-white-darkorange", node.attrs=list(color=nodes.highlight))
#' }

xEnricherGenes <- function(data, background=NULL, ontology=c("GOBP","GOMF","GOCC","PS","PS2","SF","DO","HPPA","HPMI","HPCM","HPMA","MP", "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7", "DGIdb"), size.range=c(10,2000), min.overlap=3, which.distance=NULL, test=c("hypergeo","fisher","binomial"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm=c("none","pc","elim","lea"), elim.pvalue=1e-2, lea.depth=2, path.mode=c("all_paths","shortest_paths","all_shortest_paths"), true.path.rule=F, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    ontology <- match.arg(ontology)
    test <- match.arg(test)
    p.adjust.method <- match.arg(p.adjust.method)
    ontology.algorithm <- match.arg(ontology.algorithm)
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
		## flag the simplified version of PS
		flag_PS2 <- FALSE
		if(ontology=="PS2"){
			flag_PS2 <- TRUE
			ontology <- "PS"
		}
		GS <- xRDataLoader(RData=paste('org.Hs.eg', ontology, sep=''), RData.location=RData.location, verbose=verbose)
		
		################
		if(flag_PS2){
			tmp <- as.character(unique(GS$set_info$name))
			inds <- sapply(tmp,function(x) which(GS$set_info$name==x))
		
			## new set_info
			set_info <- data.frame()
			for(i in 1:length(inds)){
				set_info<- rbind(set_info,as.matrix(GS$set_info[max(inds[[i]]),]))
			}
			## new gs
			gs <- list()
			for(i in 1:length(inds)){
				gs[[i]] <- unlist(GS$gs[inds[[i]]], use.names=F)
			}
			names(gs) <- rownames(set_info)
		
			## new GS
			GS$set_info <- set_info
			GS$gs <- gs
		}
		################
		
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
		}else{
			# force ontology.algorithm to be 'none'
			ontology.algorithm <- 'none'
		
			nodes <- data.frame(name=as.character(GS$set_info$setID), term_id=as.character(GS$set_info$setID), term_name=as.character(GS$set_info$name), term_distance=as.character(GS$set_info$distance), stringsAsFactors=F)
			nodes <- rbind(nodes, c('root','root','root','root'))
			relations <- data.frame(from='root', to=nodes$name)
			g <- igraph::graph.data.frame(d=relations, directed=T, vertices=nodes)
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
		message(sprintf("Do gene mapping from Symbols to EntrezIDs (%s) ...", as.character(now)), appendLF=T)
	}
    
    ## load Enterz Gene information
	EG <- xRDataLoader(RData=paste('org.Hs.eg', sep=''), RData.location=RData.location, verbose=verbose)	
    
	allGeneID <- EG$gene_info$GeneID
	allSymbol <- as.vector(EG$gene_info$Symbol)
	allSynonyms <- as.vector(EG$gene_info$Synonyms)
    data <- symbol2entrezgene(Symbol=data, check.symbol.identity=F, allGeneID=allGeneID, allSymbol=allSymbol, allSynonyms=allSynonyms, verbose=verbose)
    if(length(background)>0){
        background <- symbol2entrezgene(Symbol=background, check.symbol.identity=F, allGeneID=allGeneID, allSymbol=allSymbol, allSynonyms=allSynonyms, verbose=verbose)
    }
    
    #############################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=T))
        message(sprintf("'xEnricher' is being called (%s):", as.character(now)), appendLF=T)
        message(sprintf("#######################################################", appendLF=T))
    }
    eTerm <- xEnricher(data=data, annotation=anno, g=g, background=background, size.range=size.range, min.overlap=min.overlap, which.distance=which.distance, test=test, p.adjust.method=p.adjust.method, ontology.algorithm=ontology.algorithm, elim.pvalue=elim.pvalue, lea.depth=lea.depth, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose)
	
	# replace EntrezGenes with gene symbols	
	if(1 & class(eTerm)=="eTerm"){
		overlap <- eTerm$overlap
		overlap_symbols <- lapply(overlap,function(x){
			ind <- match(x, allGeneID)
			allSymbol[ind]
		})
		eTerm$overlap <- overlap_symbols
	}
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=T))
        message(sprintf("'xEnricher' has been finished (%s)!", as.character(now)), appendLF=T)
        message(sprintf("#######################################################\n", appendLF=T))
    }
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(eTerm)
}
