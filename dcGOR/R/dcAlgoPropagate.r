#' Function to propagate ontology annotations according to an input file
#'
#' \code{dcAlgoPropagate} is supposed to propagate ontology annotations, given an input file. This input file contains original annotations between domains/features and ontology terms, along with the hypergeometric scores (hscore) in support for their annotations. The annotations are propagated to the ontology root (either retaining the maximum hscore or additively accumulating the hscore). After the propogation, the ontology terms of increasing levels are determined based on the concept of Information Content (IC) to product a slim version of ontology. It returns an object of S3 class "HIS" with three components: "hscore", "ic" and "slim".
#'
#' @param input.file an input file used to build the object. This input file contains original annotations between domains/features and ontology terms, along with the hypergeometric scores (hscore) in support for their annotations. For example, a file containing original annotations between SCOP domain architectures and GO terms can be found in \url{http://dcgor.r-forge.r-project.org/data/Feature/Feature2GO.sf.txt}. As seen in this example, the input file must contain the header (in the first row) and three columns: 1st column for 'Feature_id' (here SCOP domain architectures), 2nd column for 'Term_id' (GO terms), and 3rd column for 'Score' (hscore). Alternatively, the input.file can be a matrix or data frame, assuming that input file has been read. Note: the file should use the tab delimiter as the field separator between columns
#' @param ontology the ontology identity. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPON" for Human Phenotype ONset and clinical course, "MP" for Mammalian Phenotype, "EC" for Enzyme Commission, "KW" for UniProtKB KeyWords, "UP" for UniProtKB UniPathway. For details on the eligibility for pairs of input domain and ontology, please refer to the online Documentations at \url{http://supfam.org/dcGOR/docs.html}. If NA, then the user has to input a customised RData-formatted file (see \code{RData.ontology.customised} below)
#' @param propagation how to propagate the score. It can be "max" for retaining the maximum hscore (by default), "sum" for additively accumulating the hscore
#' @param output.file an output file used to save the \code{HIS} object as an RData-formatted file (see 'Value' for details). If NULL, this file will be saved into "HIS.RData" in the current working local directory. If NA, there will be no output file
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.ontology.customised a file name for RData-formatted file containing an object of S4 class 'Onto' (i.g. ontology). By default, it is NULL. It is only needed when the user wants to perform customised analysis using their own ontology. See \code{\link{dcBuildOnto}} for how to creat this object
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{dcRDataLoader}} for details
#' @return 
#' an object of S3 class \code{HIS}, with following components:
#' \itemize{
#'  \item{\code{hscore}: a list of features, each with a term-named vector containing hscore}
#'  \item{\code{ic}: a term-named vector containing information content (IC). Terms are ordered first by IC and then by longest-path level, making sure that for terms with the same IC, parental terms always come first}
#'  \item{\code{slim}: a list of four slims, each with a term-named vector containing information content (IC). Slim '1' for very general terms, '2' for general terms, '3' for specific terms, '4' for very specific terms}
#' }
#' @note None
#' @export
#' @importFrom dnet dDAGinduce visDAG dDAGlevel dDAGroot
#' @seealso \code{\link{dcRDataLoader}}, \code{\link{dcConverter}}, \code{\link{dcAlgo}}, \code{\link{dcList2Matrix}}
#' @include dcAlgoPropagate.r
#' @examples
#' \dontrun{
#' # build an "HIS" object for GO Molecular Function
#' input.file <- "http://dcgor.r-forge.r-project.org/data/Feature/Feature2GO.sf.txt"
#' Feature2GOMF.sf <- dcAlgoPropagate(input.file=input.file, ontology="GOMF", output.file="Feature2GOMF.sf.RData")
#' names(Feature2GOMF.sf)
#' Feature2GOMF.sf$hscore[1]
#' Feature2GOMF.sf$ic[1:10]
#' Feature2GOMF.sf$slim[1]
#'
#' # extract hscore as a matrix with 3 columns (Feature_id, Term_id, Score)
#' hscore <- Feature2GOMF.sf$hscore
#' hscore_mat <- dcList2Matrix(hscore)
#' colnames(hscore_mat) <- c("Feature_id", "Term_id", "Score")
#' dim(hscore_mat)
#' hscore_mat[1:10,]
#' }

dcAlgoPropagate <- function(input.file, ontology=c(NA,"GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPON","MP","EC","KW","UP"), propagation=c("max","sum"), output.file="HIS.RData", verbose=T, RData.ontology.customised=NULL, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/dcGOR")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
    message("", appendLF=T)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    ontology <- match.arg(ontology)
    propagation <- match.arg(propagation)
    
    if(is.matrix(input.file) | is.data.frame(input.file)){
        if(verbose){
            now <- Sys.time()
            message(sprintf("Load the input file (%s) ...", as.character(now)), appendLF=T)
        }
        if(is.data.frame(input.file)){
            input <- cbind(as.character(input.file[,1]), as.character(input.file[,2]), as.character(input.file[,3]))
        }else{
            input <- input.file
        }
    }else if(is.character(input.file) & input.file!='' & !is.null(input.file) & !is.na(input.file)){
        if(verbose){
            now <- Sys.time()
            message(sprintf("Read the input file '%s' (%s) ...", input.file, as.character(now)), appendLF=T)
        }
        #tab <- read.delim(input.file, header=F, sep="\t", nrows=50, skip=1)
        #input <- read.table(input.file, header=F, sep="\t", skip=1, colClasses=sapply(tab,class))
        input <- utils::read.delim(input.file, header=T, sep="\t", colClasses="character")
    }else{
        stop("The file 'input.file' must be provided!\n")
    }

    if(is.null(output.file)){
        warnings("Since the output file is not provided, the function will use the default output file 'HIS.RData'!\n")   
        output.file <- "HIS.RData"
    }
    
    # how to deal with negative scores
    if(1){
        ## only retain all entries with positive score
        input <- input[as.numeric(input[,3])>0, ]
    }
    
    ## original annotations: Feature_id, Term_id, Score
    tmp_feature <- base::split(x=input[,1], f=input[,2], drop=T)
    tmp_score <- base::split(x=as.numeric(input[,3]), f=input[,2], drop=T)
    oAnnos <- lapply(1:length(tmp_score), function(i){
        x <- tmp_score[[i]]
        names(x) <- tmp_feature[[i]]
        return(x)
    })
    names(oAnnos) <- names(tmp_score)
    
    ## load ontology information
    if(!is.na(ontology)){
    
        if(verbose){
            now <- Sys.time()
            message(sprintf("Load the ontology '%s' (%s) ...", ontology, as.character(now)), appendLF=T)
        }
        
        #########
        ## load ontology information
        g <- dcRDataLoader(paste('onto.', ontology, sep=''), RData.location=RData.location)
        if(class(g)=="Onto"){
            g <- dcConverter(g, from='Onto', to='igraph', verbose=F)
        }
        
    }else if(file.exists(RData.ontology.customised)){
    
        if(verbose){
            now <- Sys.time()
            message(sprintf("Load customised ontology '%s' (%s)...", RData.ontology.customised, as.character(now)), appendLF=T)
        }
    
        ## load ontology informatio
        g <- ''
        eval(parse(text=paste("g <- get(load('", RData.ontology.customised,"'))", sep="")))
        if(class(g)=="Onto"){
            g <- dcConverter(g, from='Onto', to='igraph', verbose=F)
        }
        ontology <- RData.ontology.customised

    }else{
        stop("There is no input for ontology! Please input one of two parameters ('ontology' and 'RData.ontology.customised').\n")
    }
    
    ## generate a subgraph of a direct acyclic graph (DAG) induced by terms in input annotations
    dag <- dnet::dDAGinduce(g, names(oAnnos), path.mode="all_paths")
    allNodes <- V(dag)$name
    
    ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Do propagation via '%s' operation (%s) ...", propagation, as.character(now)), appendLF=T)
    }
    
    #######################
    ## Propagate using list
    if(1){
        ## initialise annotations    
        pAnnos <- oAnnos[allNodes]
        names(pAnnos) <- allNodes

        ## get the levels list
        level2node <- dnet::dDAGlevel(dag, level.mode="longest_path", return.mode="level2node")
        nLevels <- length(level2node)
        for(i in nLevels:2) {
            currNodes <- level2node[[i]]

            ## get the incoming neighbors (excluding self) that are reachable (i.e. nodes from i-1 level)
            adjNodesList <- lapply(currNodes, function(node){
                neighs.in <- igraph::neighborhood(dag, order=1, nodes=node, mode="in")
                setdiff(V(dag)[unlist(neighs.in)]$name, node)
            })
            names(adjNodesList) <- currNodes

            ## inherit the annotations from level i to level i - 1
            for(k in 1:length(currNodes)){
                node <- currNodes[k]
                ## get the domain annotations from this current node
                nowDomain <- pAnnos[[node]]
                nowDomain_mat <- cbind(names(nowDomain), as.numeric(nowDomain))
            
                ## assigin inherit annotations to all its adjacent nodes
                adjNodes <- adjNodesList[[node]]
                res <- lapply(adjNodes, function(adjNode){
                    ## get the domain annotations from this adjacent node
                    adjDomain <- pAnnos[[adjNode]]
                    adjDomain_mat <- cbind(names(adjDomain), as.numeric(adjDomain))
            
                    ### update: keep the largest score if overlap
                    all_mat <- rbind(nowDomain_mat, adjDomain_mat)
                    all_list <- base::split(x=as.numeric(all_mat[,2]), f=all_mat[,1])
                    output_list <- lapply(all_list, function(x){
                        if(propagation=='max'){
                            max(x)
                        }else if(propagation=='sum'){
                            sum(x)
                        }
                    })
                    x_mat <- base::do.call(base::rbind, output_list)
                    output <- as.vector(x_mat)
                    names(output) <- rownames(x_mat)
                    return(output)
                })
                pAnnos[adjNodes] <- res
            }
        
            if(verbose){
                message(sprintf("\tAt level %d, there are %d nodes, and %d incoming neighbors (%s).", i, length(currNodes), length(unique(unlist(adjNodesList))), as.character(Sys.time())), appendLF=T)
            }
        
        }
        

    #######################
    ## Propagate using environment
    ## no longer support!!
    #######################
    }else{
   
        if(verbose){
            now <- Sys.time()
            message(sprintf("\tassign original annotations (%s) ...", as.character(now)), appendLF=T)
        }
    
        ## node2domain.HoH: 1st key (node/term), 2nd key (domain), value (score)
        ### create a new (empty) hash environment
        node2domain.HoH <- new.env(hash=T, parent=emptyenv())
        ### assigin original annotations to "node2domain.HoH"
        tmp_trash <- lapply(allNodes, function(node){
            #message(sprintf("\t%s", node), appendLF=T)
            e <- new.env(hash=T, parent=emptyenv())
            if(node %in% names(oAnnos)){
                tmp <- oAnnos[[node]]
                if(length(tmp)>0){
                    for(i in 1:length(tmp)){
                        domain <- tmp[i]
                        assign(names(domain), as.numeric(domain), envir=e)
                    }
                }
            }
            assign(node, e, envir=node2domain.HoH)
        })
    
        if(verbose){
            now <- Sys.time()
            message(sprintf("\tpropagate annotations (%s) ...", as.character(now)), appendLF=T)
        }
    
        ## get the levels list
        level2node <- dnet::dDAGlevel(dag, level.mode="longest_path", return.mode="level2node")
        ## build a hash environment from the named list "level2node"
        ## level2node.Hash: key (level), value (a list of nodes/terms)
        level2node.Hash <- list2env(level2node)
        nLevels <- length(level2node)
        for(i in nLevels:2) {
            currNodes <- get(as.character(i), envir=level2node.Hash, mode='character')

            ## get the incoming neighbors (excluding self) that are reachable (i.e. nodes from i-1 level)
            adjNodesList <- lapply(currNodes, function(node){
                neighs.in <- igraph::neighborhood(dag, order=1, nodes=node, mode="in")
                setdiff(V(dag)[unlist(neighs.in)]$name, node)
            })
            names(adjNodesList) <- currNodes

            ## push the domains from level i to level i - 1
            lapply(currNodes, function(node){
                #message(sprintf("\t%s", node), appendLF=T)
            
                ## get the domains from this node
                nowDomain <- unlist(as.list(get(node, envir=node2domain.HoH, mode='environment')))
                domainsID <- names(nowDomain)
            
                ## assigin inherit annotations to "node2domain.HoH"
                lapply(adjNodesList[[node]], function(adjNode){
                    #message(sprintf("\t%s", adjNode), appendLF=T)
                    adjEnv <- get(adjNode, envir=node2domain.HoH, mode='environment')
                    ### domains from its adjacent nodes
                    adjDomain <- unlist(as.list(adjEnv))
                
                    ### no domains from its adjacent nodes
                    if(is.null(adjDomain)){
                        sapply(domainsID, function(domainID){
                            assign(domainID, as.numeric(nowDomain[domainID]), envir=adjEnv)
                        })
                    }else{
                        sapply(domainsID, function(domainID){
                            #message(sprintf("\t%s", domainID), appendLF=T)
                            if(is.na(adjDomain[domainID])){
                                ### missing
                                assign(domainID, as.numeric(nowDomain[domainID]), envir=adjEnv)
                            }else{
                                ### max
                                if(adjDomain[domainID] < nowDomain[domainID]){
                                    #message(sprintf("\t%s\t%s", node, adjNode), appendLF=T)
                                    assign(domainID, as.numeric(nowDomain[domainID]), envir=adjEnv)
                                }
                            }
                        })
                    }
                
                })
            })       
        
            if(verbose){
                message(sprintf("\tAt level %d, there are %d nodes, and %d incoming neighbors (%s).", i, length(currNodes), length(unique(unlist(adjNodesList))), as.character(Sys.time())), appendLF=T)
            }
        
        }
    
        ## get annotations after propagation to the root
        node2domains <- as.list(node2domain.HoH)[allNodes]
        pAnnos <- sapply(node2domains, function(node){
            unlist(as.list(node))
        })
    }
    
    ## reverse for a list of features, each containing Terms  
    res <- suppressMessages(dcList2Matrix(pAnnos))
    ### split into a list of features
    ### term, feature, score
    tmp_term <- split(x=res[,1], f=res[,2])
    tmp_score <- split(x=as.numeric(res[,3]), f=res[,2])
    fAnnos <- lapply(1:length(tmp_score), function(i){
        x <- tmp_score[[i]]
        names(x) <- tmp_term[[i]]
        return(x)
    })
    names(fAnnos) <- names(tmp_score)
    
    if(verbose){
        message(sprintf("\tafter propagation, there are %d features annotated by %d terms.", length(fAnnos), length(pAnnos)), appendLF=T)
    }
    
    #################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("Determining IC-based slim levels (%s) ...", as.character(now)), appendLF=T)
    }
    
    ## define IC
    num_f <- length(fAnnos)
    go_ic <- sapply(pAnnos, function(x){
        -1*log10(length(x)/num_f)
    })
    
    ## define the cutoff
    if(0){
        IC_cutoff <- seq(1,4)*0.5
        IC_min <- IC_cutoff - 0.25
        IC_max <- IC_cutoff + 0.25
    }else if(0){
        IC_cutoff <- max(go_ic) * seq(1,7,2)/8
        IC_min <- max(go_ic) * seq(0,6,2)/8
        IC_max <- max(go_ic) * seq(2,8,2)/8
    }else{
        nlev <- 4
        ninterval <- 1/(nlev*4) # interval units
        npoint <- seq(0,(nlev-1))/nlev + 1/(nlev*2) # central points for the levels
        naway <- c(ninterval*2, rep(ninterval, nlev-1)) # away from npoint: the wider 1st, the rest same
        IC_cutoff <- max(go_ic) * npoint
        IC_min <- max(go_ic) * (npoint - naway)
        IC_max <- max(go_ic) * (npoint + naway)
    }
    
    ## derive the IC-based levels
    level_ic <- list()
    ig <- dag
    for(i in 1:length(IC_cutoff)){
        ic_cf <- IC_cutoff[i]
        ic_min <- IC_min[i]
        ic_max <- IC_max[i]
        
        flag <- T
        marks <- vector()
        go_levels <- vector()
        while(flag){
            
            ## obtain terms with IC closest to the cutoff
            tmp <- abs(go_ic - ic_cf)
            ind <- which(is.na(match(names(tmp), marks))) # only those unmarked
            tmp_diff <- tmp[ind]
            
            ## obtain those terms closest
            go_closest <- names(tmp_diff)[which(tmp_diff==min(tmp_diff))]

            ## delete parental terms in the same path among those closest terms
            tmp_dist <- sapply(go_closest, function(from){
                vpaths <- suppressWarnings(igraph::get.shortest.paths(ig, from=from, to=go_closest, output="vpath"))
                length(unlist(vpaths))
            })
            go_rest <- names(tmp_dist[tmp_dist==1])
            
            ## mark all ancestors and descendants for @go_rest
            for(go_need in go_rest){

                to <- setdiff(V(ig)$name, marks)
                
                ### for all descendants
                neighs.out <- igraph::neighborhood(ig, order=vcount(ig), nodes=go_need, mode="out")
                nodeInduced <- V(ig)[unique(unlist(neighs.out))]$name
                if(length(nodeInduced) > 1){
                    marks <- union(marks,nodeInduced)
                }
                
                ### for all ancestors
                neighs.in <- igraph::neighborhood(ig, order=vcount(ig), nodes=go_need, mode="in")
                nodeInduced <- V(ig)[unique(unlist(neighs.in))]$name
                if(length(nodeInduced) > 1){
                    marks <- union(marks,nodeInduced)
                }
                
                if(go_ic[go_need]<ic_max & go_ic[go_need]>ic_min){
                    go_levels <- c(go_levels, go_ic[go_need])
                }
            }
            
            ## judge whether all are marked
            if(length(go_ic) == length(marks)){
                flag <- F
            }
            
        }
        
        level_ic[[i]] <- go_levels
        
        if(verbose){
            message(sprintf("\t%d level with %d terms with IC falling around %.2f (between %.2f and %.2f).", i, length(go_levels), ic_cf, ic_min, ic_max), appendLF=T)
        }
        
    }
    #sapply(level_ic,function(x) length(x))
    names(level_ic) <- paste("", 1:length(IC_cutoff), sep='')
    
    ########################
   
    # all terms in an order: first by ic and then by level (ie longest path)
    terms <- unlist(level2node, use.names=F)
    times <- sapply(level2node, function(x){
        length(x)
    })
    lvs <- rep(as.numeric(names(times)), times)
    df <- data.frame(ind=1:length(terms), terms=terms, levels=lvs, ic=as.numeric(go_ic[terms]))
    ordering <- df[base::order(df$ic,df$levels),]$ind
    ic <- df[ordering,]$ic
    names(ic) <- df[ordering,]$terms
    
    x <- list(hscore = fAnnos,
              ic     = ic,
              slim  = level_ic
              )
    class(x) <- "HIS"
    
    if(!is.na(output.file)){
        # remove the RData extension 
        output.var <- gsub(".RData$", "", output.file, ignore.case=T, perl=T)
        output.var <- gsub(".RDat$", "", output.var, ignore.case=T, perl=T)
        output.var <- gsub(".RDa$", "", output.var, ignore.case=T, perl=T)
    
        do.call(assign, list(output.var, x))
        save(list=output.var, file=output.file)
    
        if(file.exists(output.file)){
            message(sprintf("An object of S3 class 'HIS' has been built and saved into '%s'.", file.path(getwd(),output.file)), appendLF=T)
        }
        
    }
    ####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=T)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)


    invisible(x)
}