#' Function to conduct ontology enrichment analysis given a group of domains
#'
#' \code{dcEnrichment} is supposed to conduct enrichment analysis for an input group of domains using a specified ontology. It returns an object of S4 class "Eoutput". Enrichment analysis is based on either Fisher's exact test or Hypergeometric test. The test can respect the hierarchy of the ontology. The user can customise the background domains; otherwise, the function will use all annotatable domains as the test background
#'
#' @param data an input vector. It contains id for a list of domains, for example, sunids for SCOP domains
#' @param background a background vector. It contains id for a list of background domains, for example, sunids for SCOP domains. If NULL, by default all annotatable domains are used as background
#' @param domain the domain identity. It can be one of 'SCOP.sf' for SCOP superfamilies, 'SCOP.fa' for SCOP families, 'Pfam' for Pfam domains, 'InterPro' for InterPro domains, 'Rfam' for Rfam RNA families
#' @param ontology the ontology identity. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPON" for Human Phenotype ONset and clinical course, "MP" for Mammalian Phenotype, "EC" for Enzyme Commission, "KW" for UniProtKB KeyWords, "UP" for UniProtKB UniPathway. For details on the eligibility for pairs of input domain and ontology, please refer to the online Documentations at \url{http://supfam.org/dcGOR/docs.html}
#' @param sizeRange the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 1000
#' @param min.overlap the minimum number of overlaps. Only those terms that overlap with input data at least min.overlap (3 domains by default) will be processed
#' @param which_distance which distance of terms in the ontology is used to restrict terms in consideration. By default, it sets to 'NULL' to consider all distances
#' @param test the statistic test used. It can be "FisherTest" for using fisher's exact test, "HypergeoTest" for using hypergeometric test, or "BinomialTest" for using binomial test. Fisher's exact test is to test the independence between domain group (domains belonging to a group or not) and domain annotation (domains annotated by a term or not), and thus compare sampling to the left part of background (after sampling without replacement). Hypergeometric test is to sample at random (without replacement)  from the background containing annotated and non-annotated domains, and thus compare sampling to background. Unlike hypergeometric test, binomial test is to sample at random (with replacement) from the background with the constant probability. In terms of the ease of finding the significance, they are in order: hypergeometric test > binomial test > fisher's exact test. In other words, in terms of the calculated p-value, hypergeometric test < binomial test < fisher's exact test
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param ontology.algorithm the algorithm used to account for the hierarchy of the ontology. It can be one of "none", "pc", "elim" and "lea". For details, please see 'Note'
#' @param elim.pvalue the parameter only used when "ontology.algorithm" is "elim". It is used to control how to declare a signficantly enriched term (and subsequently all domains in this term are eliminated from all its ancestors)
#' @param lea.depth the parameter only used when "ontology.algorithm" is "lea". It is used to control how many maximum depth is uded to consider the children of a term (and subsequently all domains in these children term are eliminated from the use for the recalculation of the signifance at this term)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param domain.RData a file name for RData-formatted file containing an object of S4 class 'InfoDataFrame' (i.g. domain). By default, it is NULL. It is only needed when the user wants to customise enrichment analysis using their own data. See \code{\link{dcBuildInfoDataFrame}} for how to creat this object
#' @param ontology.RData a file name for RData-formatted file containing an object of S4 class 'Onto' (i.g. ontology). By default, it is NULL. It is only needed when the user wants to customise enrichment analysis using their own data. See \code{\link{dcBuildOnto}} for how to creat this object
#' @param annotations.RData a file name for RData-formatted file containing an object of S4 class 'Anno' (i.g. annotations). By default, it is NULL. It is only needed when the user wants to customise enrichment analysis using their own data. See \code{\link{dcBuildAnno}} for how to creat this object
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{dcRDataLoader}} for details
#' @return 
#' an object of S4 class \code{\link{Eoutput}}, with following slots:
#' \itemize{
#'  \item{\code{domain}: a character specifying the domain identity}
#'  \item{\code{ontology}: a character specifying the ontology used}
#'  \item{\code{term_info}: a matrix of nTerm X 5 containing term information, where nTerm is the number of terms in consideration, and the 5 columns are "term_id" (i.e. "Term ID"), "term_name" (i.e. "Term Name"), "namespace" (i.e. "Term Namespace"), "distance" (i.e. "Term Distance") and "IC" (i.e. "Information Content for the term based on annotation frequency by it")}
#'  \item{\code{anno}: a list of terms, each storing annotated domain members (also within the background domains). Always, terms are identified by "term_id" and domain members identified by their ids (e.g. sunids for SCOP domains)}
#'  \item{\code{data}: a vector containing input data in consideration. It is not always the same as the input data as only those mappable and annotatable are retained}
#'  \item{\code{background}: a vector containing background in consideration. It is not always the same as the input background as only those mappable/annotatable are retained}
#'  \item{\code{overlap}: a list of terms, each storing domains overlapped between domains annotated by a term and domains in the input data (i.e. the domains of interest). Always, terms are identified by "term_id" and domain members identified by their IDs (e.g. sunids for SCOP domains)}
#'  \item{\code{zscore}: a vector containing z-scores}
#'  \item{\code{pvalue}: a vector containing p-values}
#'  \item{\code{adjp}: a vector containing adjusted p-values. It is the p value but after being adjusted for multiple comparisons}
#' }
#' @note The interpretation of the algorithms used to account for the hierarchy of the ontology is:
#' \itemize{
#' \item{"none": does not consider the ontology hierarchy at all.}
#' \item{"lea": computers the significance of a term in terms of the significance of its children at the maximum depth (e.g. 2). Precisely, once domains are already annotated to any children terms with a more signficance than itself, then all these domains are eliminated from the use for the recalculation of the signifance at that term. The final p-values takes the maximum of the original p-value and the recalculated p-value.}
#' \item{"elim": computers the significance of a term in terms of the significance of its all children. Precisely, once domains are already annotated to a signficantly enriched term under the cutoff of e.g. pvalue<1e-2, all these domains are eliminated from the ancestors of that term).}
#' \item{"pc": requires the significance of a term not only using the whole domains as background but also using domains annotated to all its direct parents/ancestors as background. The final p-value takes the maximum of both p-values in these two calculations.}
#' \item{"Notes": the order of the number of significant terms is: "none" > "lea" > "elim" > "pc".}
#' }
#' @export
#' @importFrom dnet dDAGinduce visDAG dDAGlevel dDAGroot
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcDAGannotate}}, \code{\link{Eoutput-class}}, \code{\link{visEnrichment}}, \code{\link{dcConverter}}
#' @include dcEnrichment.r
#' @examples
#' \dontrun{
#' # 1) Enrichment analysis for SCOP domain superfamilies (sf)
#' ## 1a) load SCOP.sf (as 'InfoDataFrame' object)
#' SCOP.sf <- dcRDataLoader('SCOP.sf')
#' ### randomly select 50 domains as a list of domains of interest
#' data <- sample(rowNames(SCOP.sf), 50)
#' ## 1b) perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, domain="SCOP.sf", ontology="GOMF")
#' eoutput
#' ## 1c) view the top 10 significance terms 
#' view(eoutput, top_num=10, sortBy="pvalue", details=TRUE)
#' ## 1d) visualise the top 10 significant terms in the ontology hierarchy
#' ### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
#' visEnrichment(eoutput)
#' ## 1e) the same as above but using a customised background
#' ### randomly select 500 domains as background
#' background <- sample(rowNames(SCOP.sf), 500) 
#' ### perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, background=background, domain="SCOP.sf", ontology="GOMF")
#' eoutput
#' ### view the top 10 significance terms 
#' view(eoutput, top_num=10, sortBy="pvalue", details=TRUE)
#' ### visualise the top 10 significant terms in the ontology hierarchy
#' ### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
#' visEnrichment(eoutput)
#' 
#' ###########################################################
#' # 2) Enrichment analysis for Pfam domains (Pfam)
#' ## 2a) load Pfam (as 'InfoDataFrame' object)
#' Pfam <- dcRDataLoader('Pfam')
#' ### randomly select 100 domains as a list of domains of interest
#' data <- sample(rowNames(Pfam), 100)
#' ## 2b) perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, domain="Pfam", ontology="GOMF")
#' eoutput
#' ## 2c) view the top 10 significance terms 
#' view(eoutput, top_num=10, sortBy="pvalue", details=TRUE)
#' ## 2d) visualise the top 10 significant terms in the ontology hierarchy
#' ### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
#' visEnrichment(eoutput)
#' ## 2e) the same as above but using a customised background
#' ### randomly select 1000 domains as background
#' background <- sample(rowNames(Pfam), 1000)
#' ### perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, background=background, domain="Pfam", ontology="GOMF")
#' eoutput
#' ### view the top 10 significance terms 
#' view(eoutput, top_num=10, sortBy="pvalue", details=TRUE)
#' ### visualise the top 10 significant terms in the ontology hierarchy
#' ### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
#' visEnrichment(eoutput)
#' 
#' ###########################################################
#' # 3) Enrichment analysis for InterPro domains (InterPro)
#' ## 3a) load InterPro (as 'InfoDataFrame' object)
#' InterPro <- dcRDataLoader('InterPro')
#' ### randomly select 100 domains as a list of domains of interest
#' data <- sample(rowNames(InterPro), 100)
#' ## 3b) perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, domain="InterPro", ontology="GOMF")
#' eoutput
#' ## 3c) view the top 10 significance terms 
#' view(eoutput, top_num=10, sortBy="pvalue", details=TRUE)
#' ## 3d) visualise the top 10 significant terms in the ontology hierarchy
#' ### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
#' visEnrichment(eoutput)
#' ## 3e) the same as above but using a customised background
#' ### randomly select 1000 domains as background
#' background <- sample(rowNames(InterPro), 1000)
#' ### perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, background=background, domain="InterPro", ontology="GOMF")
#' eoutput
#' ### view the top 10 significance terms 
#' view(eoutput, top_num=10, sortBy="pvalue", details=TRUE)
#' ### visualise the top 10 significant terms in the ontology hierarchy
#' ### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
#' visEnrichment(eoutput)
#' 
#' ###########################################################
#' # 4) Enrichment analysis for Rfam RNA families (Rfam)
#' ## 4a) load Rfam (as 'InfoDataFrame' object)
#' Rfam <- dcRDataLoader('Rfam')
#' ### randomly select 100 RNAs as a list of RNAs of interest
#' data <- sample(rowNames(Rfam), 100)
#' ## 4b) perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, domain="Rfam", ontology="GOBP")
#' eoutput
#' ## 4c) view the top 10 significance terms 
#' view(eoutput, top_num=10, sortBy="pvalue", details=FALSE)
#' ## 4d) visualise the top 10 significant terms in the ontology hierarchy
#' ### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
#' visEnrichment(eoutput)
#' ## 4e) the same as above but using a customised background
#' ### randomly select 1000 RNAs as background
#' background <- sample(rowNames(Rfam), 1000)
#' ### perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, background=background, domain="Rfam", ontology="GOBP")
#' eoutput
#' ### view the top 10 significance terms 
#' view(eoutput, top_num=10, sortBy="pvalue", details=FALSE)
#' ### visualise the top 10 significant terms in the ontology hierarchy
#' ### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
#' visEnrichment(eoutput)
#'
#' ###########################################################
#' # 5) Advanced usage: customised data for domain, ontology and annotations
#' # 5a) create domain, ontology and annotations
#' ## for domain
#' domain <- dcBuildInfoDataFrame(input.file="http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt", output.file="domain.RData")
#' ## for ontology
#' dcBuildOnto(relations.file="http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_edges.txt", nodes.file="http://dcgor.r-forge.r-project.org/data/onto/igraph_GOMF_nodes.txt", output.file="ontology.RData")
#' ## for annotations
#' dcBuildAnno(domain_info.file="http://dcgor.r-forge.r-project.org/data/InterPro/InterPro.txt", term_info.file="http://dcgor.r-forge.r-project.org/data/InterPro/GO.txt", association.file="http://dcgor.r-forge.r-project.org/data/InterPro/Domain2GOMF.txt", output.file="annotations.RData")
#' ## 5b) prepare data and background
#' ### randomly select 100 domains as a list of domains of interest
#' data <- sample(rowNames(domain), 100)
#' ### randomly select 1000 domains as background
#' background <- sample(rowNames(domain), 1000)
#' ## 5c) perform enrichment analysis, producing an object of S4 class 'Eoutput'
#' eoutput <- dcEnrichment(data, background=background, domain.RData='domain.RData', ontology.RData='ontology.RData', annotations.RData='annotations.RData')
#' eoutput
#' ## 5d) view the top 10 significance terms 
#' view(eoutput, top_num=10, sortBy="pvalue", details=TRUE)
#' ### visualise the top 10 significant terms in the ontology hierarchy
#' ### color-coded according to 10-based negative logarithm of adjusted p-values (adjp)
#' visEnrichment(eoutput)
#' }

dcEnrichment <- function(data, background=NULL, domain=c(NA,"SCOP.sf","SCOP.fa","Pfam","InterPro","Rfam"), ontology=c(NA,"GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPON","MP","EC","KW","UP"), sizeRange=c(10,1000), min.overlap=3, which_distance=NULL, test=c("HypergeoTest","FisherTest","BinomialTest"), p.adjust.method=c("BH", "BY", "bonferroni", "holm", "hochberg", "hommel"), ontology.algorithm=c("none","pc","elim","lea"), elim.pvalue=1e-2, lea.depth=2, verbose=T, domain.RData=NULL, ontology.RData=NULL, annotations.RData=NULL, RData.location="http://dcgor.r-forge.r-project.org/data")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    domain <- match.arg(domain)
    ontology <- match.arg(ontology)
    test <- match.arg(test)
    p.adjust.method <- match.arg(p.adjust.method)
    ontology.algorithm <- match.arg(ontology.algorithm)
    
    if (is.vector(data)){
        data <- unique(data)
        data <- data[!is.null(data)]
        data <- data[!is.na(data)]
    }else{
        stop("The input data must be a vector.\n")
    }

    if(!is.na(domain) & !is.na(ontology)){
    
        if(verbose){
            now <- Sys.time()
            message(sprintf("First, load the ontology '%s', the domain '%s', and their associations (%s) ...", ontology, domain, as.character(now)), appendLF=T)
        }
        
        #########
        ## load ontology information
        g <- dcRDataLoader(paste('onto.', ontology, sep=''), RData.location=RData.location)
        if(class(g)=="Onto"){
            g <- dcConverter(g, from='Onto', to='igraph', verbose=F)
        }
    
        #########
        ## load domain information
        Domain <- dcRDataLoader(domain, RData.location=RData.location)
    
        #########
        ## load annotation information
        Anno <- dcRDataLoader(domain=domain, ontology=ontology, RData.location=RData.location)
        
    }else if(file.exists(domain.RData) & file.exists(ontology.RData) & file.exists(annotations.RData)){
    
        if(verbose){
            now <- Sys.time()
            message(sprintf("First, load customised ontology '%s', the domain '%s', and their associations '%s' (%s)...", ontology.RData, domain.RData, annotations.RData, as.character(now)), appendLF=T)
        }
    
        ## load ontology informatio
        g <- ''
        eval(parse(text=paste("g <- get(load('", ontology.RData,"'))", sep="")))
        if(class(g)=="Onto"){
            g <- dcConverter(g, from='Onto', to='igraph', verbose=F)
        }
        ontology <- ontology.RData
        
        ## load domain information
        Domain <- ''
        eval(parse(text=paste("Domain <- get(load('", domain.RData,"'))", sep="")))
        domain <- domain.RData
        
        ## load annotation information
        Anno <- ''
        eval(parse(text=paste("Anno <- get(load('", annotations.RData,"'))", sep="")))
    }else{
        stop("There is no input for domains and/or ontology and/or annotation.\n")
    }
    
    if(1){
        ########################
        if (is.vector(background)){
            background <- unique(background)
            background <- background[!is.null(background)]
            background <- background[!is.na(background)]
        }
        if(length(background)>0){
            ## check input background (only those domains in existance)
            ind <- match(background, rowNames(Domain))
            background <- background[!is.na(ind)]
            ## if input background cannot be found in annotable domains, then use the all annotatable domains as background
            if(length(background)>0){
                
                ## background should be: customised background plus input domains of interest
                background <- union(background, data)
                ###########################################
                
                ind <- match(domainNames(Anno), background)
                Anno <- Anno[!is.na(ind), ]
                #data <- intersect(data, background)
            }
        }
        ########################
    }
    
    #########
    ## obtain the induced subgraph according to the input annotation data
    ## based on all possible paths (i.e. the complete subgraph induced)
    dag <- dcDAGannotate(g, annotations=Anno, path.mode="all_paths", verbose=F)
    
    #####################################################
    
    ## check input data (only those domains in existance)
    ind <- match(data, rowNames(Domain))
    domains.group <- data[!is.na(ind)]
    
    ## filter based on "which_distance"
    distance <- V(dag)$term_distance
    if(!is.null(which_distance) & sum(is.na(distance))==0){
        set_filtered <- sapply(which_distance, function(x) {
            V(dag)$term_id[(distance==as.integer(x))]
        })
        set_filtered <- unlist(set_filtered)
    }else{
        set_filtered <- V(dag)$term_id
    }
    dag <- dnet::dDAGinduce(dag, nodes_query=set_filtered, path.mode="all_paths")
    
    ## filter based on "sizeRange"
    gs.length <- sapply(V(dag)$annotations, length)
    ind.length <- which(gs.length >= sizeRange[1] & gs.length <= sizeRange[2])
    dag <- dnet::dDAGinduce(dag, nodes_query=V(dag)$term_id[ind.length], path.mode="all_paths")
    
    if(length(V(dag))==0){
        stop("There is no term being used.\n")
    }

    ##############################################################################################
    ## Fisher's exact test: testing the independence between domain group (domains belonging to a group or not) and domain annotation (domains annotated by a term or not); thus compare sampling to the left part of background (after sampling without replacement)
    doFisherTest <- function(domains.group, domains.term, domains.universe){
        domains.hit <- intersect(domains.group, domains.term)
        # num of success in sampling
        X <- length(domains.hit)
        # num of sampling
        K <- length(domains.group)
        # num of success in background
        M <- length(domains.term)
        # num in background
        N <- length(domains.universe)
        ## Prepare a two-dimensional contingency table: #success in sampling, #success in background, #failure in sampling, and #failure in left part
        cTab <- matrix(c(X, K-X, M-X, N-M-K+X), nrow=2, dimnames=list(c("anno", "notAnno"), c("group", "notGroup")))
        p.value <- ifelse(all(cTab==0), 1, stats::fisher.test(cTab, alternative="greater")$p.value)
        return(p.value)
    }

    ## Hypergeometric test: sampling at random from the background containing annotated and non-annotated domains (without replacement); thus compare sampling to background
    doHypergeoTest <- function(domains.group, domains.term, domains.universe){
        domains.hit <- intersect(domains.group, domains.term)
        # num of success in sampling
        X <- length(domains.hit)
        # num of sampling
        K <- length(domains.group)
        # num of success in background
        M <- length(domains.term)
        # num in background
        N <- length(domains.universe)
    
        x <- X
        m <- M
        n <- N-M # num of failure in background
        k <- K
        p.value <- ifelse(m==0 || k==0, 1, stats::phyper(x,m,n,k, lower.tail=F, log.p=F))
        return(p.value)
    }
    
    
    ## Binomial test: sampling at random from the background with the constant probability of having annotated domains (with replacement)
    doBinomialTest <- function(domains.group, domains.term, domains.universe){
        domains.hit <- intersect(domains.group, domains.term)
        # num of success in sampling
        X <- length(domains.hit)
        # num of sampling
        K <- length(domains.group)
        # num of success in background
        M <- length(domains.term)
        # num in background
        N <- length(domains.universe)
    
        p.value <- ifelse(K==0 || M==0 || N==0, 1, stats::pbinom(X,K,M/N, lower.tail=F, log.p=F))
        return(p.value)
    }
    
    
    ##  Z-score from hypergeometric distribution
    zscoreHyper <- function(domains.group, domains.term, domains.universe){
        domains.hit <- intersect(domains.group, domains.term)
        # num of success in sampling
        X <- length(domains.hit)
        # num of sampling
        K <- length(domains.group)
        # num of success in background
        M <- length(domains.term)
        # num in background
        N <- length(domains.universe)
        
        ## calculate z-score
        if(1){
            ## Z-score based on theoretical calculation
            x.exp <- K*M/N
            var.exp <- K*M/N*(N-M)/N*(N-K)/(N-1)
            if(var.exp==0){
                z <- NA
            }else{
                suppressWarnings(z <- (X-x.exp)/sqrt(var.exp))
            }
            
        }else{
            ## Z-score equivalents for deviates from hypergeometric distribution
            x <- X
            m <- M
            n <- N-M # num of failure in background
            k <- K
            
            suppressWarnings(d <- stats::dhyper(x,m,n,k,log=TRUE)-log(2))
            suppressWarnings(pupper <- stats::phyper(x,m,n,k,lower.tail=FALSE,log.p=TRUE))
            suppressWarnings(plower <- stats::phyper(x-1,m,n,k,lower.tail=TRUE,log.p=TRUE))
            d[is.na(d)] <- -Inf
            pupper[is.na(pupper)] <- -Inf
            plower[is.na(plower)] <- -Inf

            # Add half point probability to upper tail probability preserving log-accuracy
            a <- pupper
            b <- d-pupper
            a[b>0] <- d[b>0]
            b <- -abs(b)
            pmidupper <- a+log1p(exp(b))
            pmidupper[is.infinite(a)] <- a[is.infinite(a)]

            # Similarly for lower tail probability preserving log-accuracy
            a <- plower
            b <- d-plower
            a[b>0] <- d[b>0]
            b <- -abs(b)
            pmidlower <- a+log1p(exp(b))
            pmidlower[is.infinite(a)] <- a[is.infinite(a)]

            up <- pmidupper<pmidlower
            if(any(up)) z <- stats::qnorm(pmidupper,lower.tail=FALSE,log.p=TRUE)
            if(any(!up)) z <- stats::qnorm(pmidlower,lower.tail=TRUE,log.p=TRUE)
        }
        
        return(z)
    }
    ##############################################################################################
    
    terms <- V(dag)$term_id
    gs <- V(dag)$annotations
    names(gs) <- terms
    domains.universe <- unique(unlist(V(dag)$annotations))
    domains.group <- intersect(domains.universe, domains.group)
    
    if(length(domains.group)==0){
        warnings("There is no domain being used.\n")
        return(F)
    }
    
    if(ontology.algorithm=="none"){
    
        if(verbose){
            now <- Sys.time()
            message(sprintf("Second, perform enrichment analysis using %s (%s) ...", test, as.character(now)), appendLF=T)
            if(is.null(which_distance)){
                message(sprintf("\tThere are %d terms being used, each restricted within [%s] annotations", length(terms), paste(sizeRange,collapse=",")), appendLF=T)
            }else{
                message(sprintf("\tThere are %d terms being used, each restricted within [%s] annotations and [%s] distance", length(terms), paste(sizeRange,collapse=","), paste(which_distance,collapse=",")), appendLF=T)
            }
        }
    
        pvals <- sapply(terms, function(term){
            domains.term <- unique(unlist(gs[term]))
            p.value <- switch(test,
                FisherTest =  doFisherTest(domains.group, domains.term, domains.universe),
                HypergeoTest = doHypergeoTest(domains.group, domains.term, domains.universe),
                BinomialTest = doBinomialTest(domains.group, domains.term, domains.universe)
            )
        })
        
        zscores <- sapply(terms, function(term){
            domains.term <- unique(unlist(gs[term]))
            zscoreHyper(domains.group, domains.term, domains.universe)
        })

    }else if(ontology.algorithm=="pc" || ontology.algorithm=="elim" || ontology.algorithm=="lea"){

        if(verbose){
            now <- Sys.time()
            message(sprintf("Third, perform enrichment analysis using %s based on %s algorithm to respect ontology structure (%s) ...", test, ontology.algorithm, as.character(now)), appendLF=T)
        }
        
        ###############################
        subg <- dag
        
        if(verbose){
            message(sprintf("\tThere are %d terms being used", length(V(subg))), appendLF=T)
        }
        
        level2node <- dnet::dDAGlevel(subg, level.mode="longest_path", return.mode="level2node")
        
        ## build a hash environment from the named list "level2node"
        ## level2node.Hash: key (level), value (a list of nodes/terms)
        level2node.Hash <- list2env(level2node)
        ## ls(level2node.Hash)
        nLevels <- length(level2node)
        
        ## create a new (empty) hash environment
        ## node2pval.Hash: key (node), value (pvalue)
        node2pval.Hash <- new.env(hash=T, parent=emptyenv())        
        ## node2zscore.Hash: key (node), value (zscore)
        node2zscore.Hash <- new.env(hash=T, parent=emptyenv())
        
        if(ontology.algorithm=="pc"){
        
            for(i in nLevels:2) {
                currNodes <- get(as.character(i), envir=level2node.Hash, mode="character")
    
                for(currNode in currNodes){
                    domains.term <- unique(unlist(gs[currNode]))
                
                    ## do test based on the whole domains as background
                    pvalue_whole <- switch(test,
                        FisherTest =  doFisherTest(domains.group, domains.term, domains.universe),
                        HypergeoTest = doHypergeoTest(domains.group, domains.term, domains.universe),
                        BinomialTest = doBinomialTest(domains.group, domains.term, domains.universe)
                    )
                    zscore_whole <- zscoreHyper(domains.group, domains.term, domains.universe)
            
                    ## get the incoming neighbors/parents (including self) that are reachable
                    neighs.in <- igraph::neighborhood(subg, order=1, nodes=currNode, mode="in")
                    adjNodes <- setdiff(V(subg)[unlist(neighs.in)]$name, currNode)
                
                    ## domains annotated in parents are as background
                    domains.parent <- unique(unlist(gs[adjNodes]))
        
                    ## make sure domains in group (domains in term) are also in parents
                    domains.group.parent <- intersect(domains.group, domains.parent)
                    domains.term.parent <- intersect(domains.term, domains.parent)

                    ## do test based on the domains in parents as background
                    pvalue_relative <- switch(test,
                        FisherTest =  doFisherTest(domains.group.parent, domains.term.parent, domains.parent),
                        HypergeoTest = doHypergeoTest(domains.group.parent, domains.term.parent, domains.parent),
                        BinomialTest = doBinomialTest(domains.group.parent, domains.term.parent, domains.parent)
                    )
                    zscore_relative <- zscoreHyper(domains.group.parent, domains.term.parent, domains.parent)
                
                    ## take the maximum value of pvalue_whole and pvalue_relative
                    pvalue <- max(pvalue_whole, pvalue_relative)
                    ## store the result (the p-value)
                    assign(currNode, pvalue, envir=node2pval.Hash)
                    
                    ## take the miminum value of zscore_whole and zscore_relative
                    zscore <- ifelse(pvalue_whole>pvalue_relative, zscore_whole, zscore_relative)
                    ## store the result (the z-score)
                    assign(currNode, zscore, envir=node2zscore.Hash)
                }
                
                if(verbose){
                    message(sprintf("\tAt level %d, there are %d nodes/terms", i, length(currNodes), appendLF=T))
                }
            }
            
            ## the root always has p-value=1 and z-score=0
            root <- dnet::dDAGroot(subg)
            assign(root, 1, envir=node2pval.Hash)
            assign(root, 0, envir=node2zscore.Hash)
        
        }else if(ontology.algorithm=="elim"){
        
            ## sigNode2pval.Hash: key (node called significant), value (pvalue)
            sigNode2pval.Hash <- new.env(hash=T, parent=emptyenv())
            ## ancNode2domain.Hash: key (node at ancestor), value (domains to be eliminated)
            ancNode2domain.Hash <- new.env(hash=T, parent=emptyenv())
            
            if(is.null(elim.pvalue) || is.na(elim.pvalue) || elim.pvalue>1 || elim.pvalue<0){
                elim.pvalue <- 1e-2
            }
            pval.cutoff <- elim.pvalue

            #pval.cutoff <- 1e-2 / length(V(subg))
            
            for(i in nLevels:1) {
                currNodes <- get(as.character(i), envir=level2node.Hash, mode="character")
                currAnno <- gs[currNodes]
    
                ## update "ancNode2domain.Hash" for each node/term
                for(currNode in currNodes){
                    domains.term <- unique(unlist(gs[currNode]))
        
                    ## remove the domains (if any already marked) from annotations by the current node/term
                    if(exists(currNode, envir=ancNode2domain.Hash, mode="numeric")){
                        domains.elim <- get(currNode, envir=ancNode2domain.Hash, mode="numeric")
                        domains.term <- setdiff(domains.term, domains.elim)
                        #message(sprintf("\t\t%d %d", length(domains.elim), length(domains.term)), appendLF=T)
                    }
        
                    ## do test
                    pvalue <- switch(test,
                        FisherTest =  doFisherTest(domains.group, domains.term, domains.universe),
                        HypergeoTest = doHypergeoTest(domains.group, domains.term, domains.universe),
                        BinomialTest = doBinomialTest(domains.group, domains.term, domains.universe)
                    )
                    zscore <- zscoreHyper(domains.group, domains.term, domains.universe)
                    
                    ## store the result (the p-value)
                    assign(currNode, pvalue, envir=node2pval.Hash)
                    ## store the result (the z-score)
                    assign(currNode, zscore, envir=node2zscore.Hash)
                    
                    ## condition to update "ancNode2domain.Hash"
                    if(pvalue < pval.cutoff) {
                        ## mark the significant node
                        assign(currNode, pvalue, envir=sigNode2pval.Hash)

                        ## retrieve domains annotated by the significant node for the subsequent eliminating
                        elimGenesID <- currAnno[[currNode]]

                        ## find all the ancestors of the significant node
                        dag.ancestors <- dnet::dDAGinduce(subg, currNode, path.mode="all_paths")
                        ancestors <- setdiff(V(dag.ancestors)$name, currNode)
            
                        ## get only those ancestors that are already in "ancNode2domain.Hash"
                        oldAncestors2GenesID <- sapply(ancestors, function(ancestor){
                            if(exists(ancestor, envir=ancNode2domain.Hash, mode="numeric")){
                                get(ancestor, envir=ancNode2domain.Hash, mode='numeric')
                            }
                        })

                        ## add the new GenesID to the ancestors
                        newAncestors2GenesID <- lapply(oldAncestors2GenesID, function(oldGenes){
                            union(oldGenes, elimGenesID)
                        })

                        ## update the "ancNode2domain.Hash" table
                        if(length(newAncestors2GenesID) > 0){
                            sapply(names(newAncestors2GenesID), function(ancestor){
                                assign(ancestor, newAncestors2GenesID[[ancestor]], envir=ancNode2domain.Hash)
                            })
                        }
                    }
                }
                
                if(verbose){
                    num.signodes <- length(ls(sigNode2pval.Hash))
                    num.ancnodes <- length(ls(ancNode2domain.Hash))
                    num.elimdomains <- length(unique(unlist(as.list(ancNode2domain.Hash))))
                    message(sprintf("\tAt level %d, there are %d nodes/terms: up to %d significant nodes, %d ancestral nodes changed (%d domains eliminated)", i, length(currNodes), num.signodes, num.ancnodes, num.elimdomains), appendLF=T)
                }
            }
            
        }else if(ontology.algorithm=="lea"){
        
            ## node2pvalo.Hash: key (node called significant), value (original pvalue)
            node2pvalo.Hash <- new.env(hash=T, parent=emptyenv())
        
            if(is.null(lea.depth) || is.na(lea.depth) || lea.depth<0){
                lea.depth <- 2
            }
            depth.cutoff <- as.integer(lea.depth)
            
            for(i in nLevels:1) {
                currNodes <- get(as.character(i), envir=level2node.Hash, mode="character")
                currAnno <- gs[currNodes]
                
                num.recalculate <- 0
                
                ## update "node2pval.Hash" for each node/term
                for(currNode in currNodes){
                    domains.term <- unique(unlist(gs[currNode]))
                    
                    ## do test
                    pvalue.old <- switch(test,
                        FisherTest =  doFisherTest(domains.group, domains.term, domains.universe),
                        HypergeoTest = doHypergeoTest(domains.group, domains.term, domains.universe),
                        BinomialTest = doBinomialTest(domains.group, domains.term, domains.universe)
                    )
                    zscore.old <- zscoreHyper(domains.group, domains.term, domains.universe)
                    
                    ## store the result (old pvalue)
                    assign(currNode, pvalue.old, envir=node2pvalo.Hash)
                    
                    ## get the outgoing neighbors/children (including self) that are reachable at most of given depth
                    neighs.out <- igraph::neighborhood(subg, order=depth.cutoff, nodes=currNode, mode="out")
                    adjNodes <- setdiff(V(subg)[unlist(neighs.out)]$name, currNode)
                        
                    if(length(adjNodes)!=0){
                        ## get children with the lower p-value
                        if(1){
                            pvalue.children <- sapply(adjNodes, function(child){
                                if(exists(child, envir=node2pvalo.Hash, mode="numeric")){
                                    get(child, envir=node2pvalo.Hash, mode="numeric")
                                }
                            })
                        }else{
                            pvalue.children <- sapply(adjNodes, function(child){
                                if(exists(child, envir=node2pval.Hash, mode="numeric")){
                                    get(child, envir=node2pval.Hash, mode="numeric")
                                }
                            })
                        }
                        
                        chNodes <- names(pvalue.children[pvalue.children < pvalue.old])
                        
                        ## whether there exist any children with the lower p-value
                        if(length(chNodes)>0){
                            num.recalculate <- num.recalculate + 1
                        
                            ## if yes, get domains that are annotated by children with the lower p-value
                            ## they will be removed
                            domains.elim <- unique(unlist(gs[chNodes]))
                            domains.term.new <- setdiff(domains.term, domains.elim)
                            
                            ## recalculate the significance
                            pvalue.new <- switch(test,
                                FisherTest =  doFisherTest(domains.group, domains.term.new, domains.universe),
                                HypergeoTest = doHypergeoTest(domains.group, domains.term.new, domains.universe),
                                BinomialTest = doBinomialTest(domains.group, domains.term.new, domains.universe)
                            )
                            zscore.new <- zscoreHyper(domains.group, domains.term.new, domains.universe)
                            
                            ## take the maximum value of pvalue_new and the original pvalue
                            pvalue <- max(pvalue.new, pvalue.old)
                            
                            ## take the minimum value of zscore_new and the original zscore
                            zscore <- ifelse(pvalue.new>pvalue.old, zscore.new, zscore.old)
                            
                        }else{
                            pvalue <- pvalue.old
                            zscore <- zscore.old
                        }
                        
                    }else{
                        pvalue <- pvalue.old
                        zscore <- zscore.old
                    }
                    
                    ## store the result (recalculated pvalue if have to)
                    assign(currNode, pvalue, envir=node2pval.Hash)
                    
                    ## store the result (recalculated zscore if have to)
                    assign(currNode, zscore, envir=node2zscore.Hash)
                }
    
                if(verbose){
                    message(sprintf("\tAt level %d, there are %d nodes/terms and %d being recalculated", i, length(currNodes), num.recalculate), appendLF=T)
                }
            
            }
        }
        
        pvals <- unlist(as.list(node2pval.Hash))
        zscores <- unlist(as.list(node2zscore.Hash))
    
    }

    if(verbose){
        now <- Sys.time()
        message(sprintf("Last, adjust the p-values using the %s method (%s) ...", p.adjust.method, as.character(now)), appendLF=T)
    }

    overlaps <- sapply(names(gs), function(term){
        domains.term <- unique(unlist(gs[term]))
        intersect(domains.group, domains.term)

    })
    ## for those with "min.overlap" overlaps will be processed and reported
    flag_filter <- sapply(overlaps, function(x) ifelse(length(x)>=min.overlap,T,F))
    
    if(sum(flag_filter)==0){
        warnings("It seems there are no terms meeting the specified 'sizeRange' and 'min.overlap'.\n")
        return(F)
    }
    
    gs <- gs[flag_filter]
    overlaps <- overlaps[flag_filter]
    
    ## common terms
    common <- intersect(names(gs), names(zscores))
    ind_gs <- match(common,names(gs))
    ind_zscores <- match(common, names(zscores))
    
    ## restrict to the common terms (and sorted too)
    gs <- gs[ind_gs[!is.na(ind_gs)]]
    overlaps <- overlaps[ind_gs[!is.na(ind_gs)]]
    zscores <- zscores[ind_zscores[!is.na(ind_zscores)]]
    pvals <- pvals[ind_zscores[!is.na(ind_zscores)]]
    
    ## remove those with zscores=NA
    flag <- !is.na(zscores)
    gs <- gs[flag]
    overlaps <- overlaps[flag]
    zscores <- zscores[flag]
    pvals <- pvals[flag]
    
    zscores <- signif(zscores, digits=3)
    pvals <- sapply(pvals, function(x) min(x,1))
    ## Adjust P-values for multiple comparisons
    adjpvals <- stats::p.adjust(pvals, method=p.adjust.method)
    
    pvals <- signif(pvals, digits=2)
    adjpvals <- sapply(adjpvals, function(x) min(x,1))
    adjpvals <- signif(adjpvals, digits=2)
    
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    set_info <- get.data.frame(dag, what="vertices")[names(gs), c(2:5,7)]
    annotations <- V(dag)$annotations
    names(annotations) <- V(dag)$term_id
    annotations <- annotations[names(gs)]
    
    ########################################
    if(1){
        ## overlaps
        overlaps <- lapply(overlaps, function(x){
            ind <- match(x, rowNames(Domain))
            names(x) <- as.character(Domain@data$description[ind])
            x
        })
    }
    ########################################
    
    eoutput <- new("Eoutput",
                    domain    = domain,
                    ontology  = ontology,
                    term_info = set_info,
                    anno      = annotations,
                    data      = domains.group,
                    background=domains.universe,
                    overlap   = overlaps,
                    zscore    = zscores,
                    pvalue    = pvals,
                    adjp      = adjpvals
                 )

    invisible(eoutput)
}