#' Function to reconstruct ancestral discrete states using maximum parsimony algorithm
#'
#' \code{dcAncestralMP} is supposed to reconstruct ancestral discrete states using a maximum parsimony-modified Fitch algorithm. In a from-tip-to-root manner, ancestral state for an internal node is determined if a state is shared in a majority by all its children. If two or more states in a majority are equally shared, this internal node is temporarily marked as an unknown tie, which is further resolved in a from-root-to-tip manner: always being the same state as its direct parent holds. If the ties also occur at the root, the state at the root is set to the last state in ties (for example, usually being 'present' for 'present'-'absent' two states).
#'
#' @param data an input data matrix/frame storing discrete states for tips (in rows) X characters (in columns). The rows in the matrix are for tips. If the row names do not exist, then addumedly they have the same order as in the tree tips. More wisely, users provide row names which can be matched to the tip labels of the tree. The row names can be more than found in the tree labels, and they should contain all those in the tree labels
#' @param phy an object of class 'phylo'
#' @param output.detail logical to indicate whether the output is returned as a detailed list. If TRUE, a nested list is returned: a list of characters (corresponding to columns of input data matrix), in which each element is a list consisting of three components ("states", "transition" and "relative"). If FALSE, a matrix is returned: the columns respond to the input data columns, and rows responding to all node index in the phylo-formatted tree
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("foreach","doMC"))}. If not yet installed, this option will be disabled
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' It depends on the 'output.detail'.
#' If FALSE (by default), a matrix is returned, with the columns responding to the input data columns, and rows responding to node index in the phylo-formatted tree.
#' If TRUE, a nested list is returned. Outer-most list is for characters (corresponding to columns of input data matrix), in which each elemenl is a list (inner-most) consisting of three components ("states", "transition" and "relative"):
#' \itemize{
#'  \item{\code{states}: a named vector storing states (extant and ancestral states)}
#'  \item{\code{transition}: a posterior transition matrix between states}
#'  \item{\code{relative}: a matrix of nodes X states, storing relative probability}
#' }
#' @note
#' This maximum parsimony algorithm for ancestral discrete state reconstruction is attributable to the basic idea as described in \url{http://sysbio.oxfordjournals.org/content/20/4/406.short}
#' @export
#' @seealso \code{\link{dcAncestralML}}, \code{\link{dcTreeConnectivity}}, \code{\link{dcDuplicated}}
#' @include dcAncestralMP.r
#' @examples
#' # 1) a newick tree that is imported as a phylo-formatted tree
#' tree <- "(((t1:5,t2:5):2,(t3:4,t4:4):3):2,(t5:4,t6:4):6);"
#' phy <- ape::read.tree(text=tree)
#'
#' # 2) an input data matrix storing discrete states for tips (in rows) X four characters (in columns)
#' data1 <- matrix(c(0,rep(1,3),rep(0,2)), ncol=1)
#' data2 <- matrix(c(rep(0,4),rep(1,2)), ncol=1)
#' data <- cbind(data1, data1, data1, data2)
#' colnames(data) <- c("C1", "C2", "C3", "C4")
#' ## reconstruct ancestral states, without detailed output
#' res <- dcAncestralMP(data, phy, parallel=FALSE)
#' res
#'
#' # 3) an input data matrix storing discrete states for tips (in rows) X only one character
#' data <- matrix(c(0,rep(1,3),rep(0,2)), ncol=1)
#' ## reconstruct ancestral states, with detailed output
#' res <- dcAncestralMP(data, phy, parallel=FALSE, output.detail=TRUE)
#' res
#' ## get the inner-most list
#' res <- res[[1]]
#' ## visualise the tree with ancestral states and their conditional probability
#' Ntip <- ape::Ntip(phy)
#' Nnode <- ape::Nnode(phy)
#' color <- c("white","gray")
#' ## visualise main tree
#' ape::plot.phylo(phy, type="p", use.edge.length=TRUE, label.offset=1, show.tip.label=TRUE, show.node.label=FALSE)
#' ## visualise tips (state 1 in gray, state 0 in white)
#' x <- data[,1]
#' ape::tiplabels(pch=22, bg=color[as.numeric(x)+1], cex=2, adj=1)
#' ## visualise internal nodes
#' ### thermo bar to illustrate relative probability (state 1 in gray, state 0 in white)
#' ape::nodelabels(thermo=res$relative[Ntip+1:Nnode,2:1], piecol=color[2:1], cex=0.75)
#' ### labeling reconstructed ancestral states
#' ape::nodelabels(text=res$states[Ntip+1:Nnode], node=Ntip+1:Nnode, frame="none", col="red", bg="transparent", cex=0.75)

dcAncestralMP <- function(data, phy, output.detail=F, parallel=T, multicores=NULL, verbose=T)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################
    
    if (class(phy) != "phylo"){
        stop("The input 'phy' must belong to the class 'phylo'!")
    }
    
    Ntip <- ape::Ntip(phy)
    Nnode <- ape::Nnode(phy)
    Ntot <- Ntip+Nnode

    # In the "postorder" order, the rows are arranged so that postorder tree traversal can be done by descending along these rows
    phy <- ape::reorder.phylo(phy, "postorder")
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    
    ## calculate the sparse connectivity matrix between parents and children
    connectivity <- suppressMessages(dcTreeConnectivity(phy, verbose=verbose))
    
    if (Nnode != Ntip-1){
        stop("The input 'phy' is not binary and rooted!")
    }
    
    if(is.vector(data)){
        tmp_data <- matrix(data, ncol=1)
        if(!is.null(names(data))){
            rownames(tmp_data) <- names(data)
        }
        data <- tmp_data
    }
    if(is.data.frame(data)){
        data <- as.matrix(data)
    }
    
    if (!is.null(rownames(data))) {
        
        ind <- match(rownames(data), phy$tip.label)
        data <- data[!is.na(ind),]
        
        if(nrow(data) != Ntip){
            stop(message(sprintf("The row names of input 'data' do not contain all of the tip labels of the input 'phy': %d NOT FOUND!", Ntip-nrow(data)), appendLF=T))
        }
        
        ind <- match(phy$tip.label, rownames(data))
        data <- data[ind,]
    }else{
        if(nrow(data) != Ntip){
            stop(message(sprintf("The row number of input 'data' do not equal the tip number of the input 'phy'!"), appendLF=T))
        }
    }

    if(verbose){
        message(sprintf("The input tree has '%d' tips.", Ntip), appendLF=T)
    }

    ######################################################################################
    ## A function to do prediction
    doReconstruct <- function(x, Ntot, Ntip, Nnode, connectivity, e1, e2, output.detail, verbose){
    
        if (!is.factor(x)){
            x_tmp <- base::factor(x)
        }else{
            x_tmp <- x
        }
        nl <- base::nlevels(x_tmp)
        lvls <- base::levels(x_tmp)
        x_tmp <- as.integer(x_tmp)
    
        ##########################
        if(nl==1){

            if(verbose){
                message(sprintf("\tNote, there is only one state '%s' in tips", lvls), appendLF=T)
            }
            
            if(output.detail){
                res <- list()
                res$states <- rep(lvls, Ntot)
            }else{
                res <- rep(lvls, Ntot)
            }
            
            return(invisible(res))
        }
        ########################## 
    
        ################################################################################################
        if(verbose){
            message(sprintf("\tFirst, do maximum parsimony-modified Fitch algorithm in a bottom-up manner (%s) ...", as.character(Sys.time())), appendLF=T)
        }

        ## dimension: #nodes X #states (current node)
        Cx <- matrix(NA, nrow=Ntot, ncol=nl, dimnames=list(1:Ntot, lvls))
    
        ## for tips
        Cx[cbind(1:Ntip, x_tmp)] <- 1
        
        Cx_final <- Cx
        ## for internal nodes (in a postordered manner)
        for (i in seq(from=1, by=2, length.out=Nnode)) {
        
            j <- i + 1L
            cur <- e1[i]
        
            all_children <- which(connectivity[cur-Ntip,]==1)
        
            ## only tips
            if(0){
                ind <- match(1:Ntip, all_children)
                all_children <- all_children[ind[!is.na(ind)]]
            }
        
            ### do calculation
            tmp <- Cx[all_children, ]
            ttmp <- apply(tmp, 2, function(x) sum(x,na.rm=T))
            ind <- which(ttmp==max(ttmp))
            if(length(ind)==1){
                Cx[cur, ind] <- 1
                Cx_final[cur, ind] <- 1
            }else{
                Cx[cur, ind] <- NA
                Cx_final[cur, ind] <- Inf
            }
        
        }
    
        anc <- apply(Cx, 1, function(x){
            tmp <- lvls[which(x==1)]
            if(length(tmp)==0){
                return("tie")
            }else{
                return(tmp)
            }
        })

        ################################################################################################
        if(verbose){
            message(sprintf("\tSecond, resolve unknown states being the same as its direct parent in a top-down manner (%s) ...", as.character(Sys.time())), appendLF=T)
        }
    
        # break ties: the tie always follows the direct parent state (in a preorder)
        anc_final <- anc
        ties <- which(anc_final=='tie')
        if(length(ties) > 0){
            if(verbose){
                message(sprintf("\t\tbreak %d tie(s)", length(ties)), appendLF=T)
            }
            for(i in 1:length(ties)){
                child_ind <- ties[i]
                if(child_ind==Ntip+1){
                    # break the tie at the root
                    ind <- which(is.infinite(Cx_final[child_ind,]))
                    anc_final[child_ind] <- lvls[ind[length(ind)]]
                    Cx_final[child_ind, ind[length(ind)]] <- 1
                }else{
                    child <- names(child_ind)
                    parent <- e1[match(child, e2)]
                    parent_ind <- match(parent, names(anc))
                    anc_final[child_ind] <- anc_final[parent_ind]
                    Cx_final[child_ind, match(anc_final[parent_ind],lvls)] <- 1
                }
            }
        }
        Cx_final[is.infinite(Cx_final)] <- NA
    
        ####################################################################################
        
        ## A summary of changes between states
        p2c <- cbind(anc_final[e1], anc_final[e2])
        p2c_final <- p2c
        ### for the root: being present
        if(nl==2){
            if(anc_final[e1[2*Nnode-1]] == lvls[nl]){
                p2c_final <- rbind(p2c_final, lvls)
            }
        }
        
        if(verbose){
            all <- paste("\t", p2c_final[,1], "->", p2c_final[,2], sep='')
            changes <- sapply(unique(all), function(x){
                sum(x==all)
            })
            changes <- sort(changes)
            msg <- paste(names(changes),changes, sep=": ", collapse="\n")
    
            message(sprintf("\tIn summary, the number of between-state changes:\n%s\n", msg), appendLF=T)
        }
        
        if(output.detail){
            ### Calculate relative probability
            ## dimension: #nodes X #probability (current node)
            Lx <- matrix(0, nrow=Ntot, ncol=nl, dimnames=list(1:Ntot, lvls))
            ## for tips
            Lx[cbind(1:Ntip, x_tmp)] <- 1
            ## for internal nodes (in a postordered manner)
            for (i in seq(from=1, by=2, length.out=Nnode)) {
                j <- i + 1L
                cur <- e1[i]
                all_children <- which(connectivity[cur-Ntip,]==1)
                ## only tips
                if(0){
                    ind <- match(1:Ntip, all_children)
                    all_children <- all_children[ind[!is.na(ind)]]
                }
                ### do calculation
                tmp <- Cx_final[all_children, ]
                ttmp <- apply(tmp, 2, function(x) sum(x,na.rm=T))
                Lx[cur,] <- ttmp
            }
            ## relative probability
            rp <- t(apply(Lx, 1, function(x) x/sum(x)))
    
            ### Calculate posterior transition matrix
            rate <- matrix(0, nl, nl)
            ind <- matrix(match(p2c_final, lvls),ncol=2)
            for(i in 1:nrow(ind)){
                rate[ind[i,1], ind[i,2]] <- rate[ind[i,1], ind[i,2]] + 1
            }
            colnames(rate) <- rownames(rate) <- lvls
        }

        if(output.detail){
            res <- list()
            res$states <- anc_final
            res$transition <- rate
            res$relative <- signif(rp, digits=4)
        }else{
            res <- anc_final
        }
        
        return(invisible(res))
    }
    
    ## A function to indicate the running progress
    progress_indicate <- function(i, B, step, flag=F){
        if(i %% ceiling(B/step) == 0 | i==B | i==1){
            if(flag & verbose){
                message(sprintf("\t%d out of %d (%s)", i, B, as.character(Sys.time())), appendLF=T)
            }
        }
    }
    
    ######################################################################################
    
    integer_vec <- suppressMessages(dcDuplicated(data, pattern.wise="column", verbose=verbose))
    ind_unique <- sort(unique(integer_vec))
    data_unique <- as.matrix(data[, ind_unique], ncol=length(ind_unique))
    
    if(verbose){
        message(sprintf("The input data has %d characters/columns (with %d distinct patterns).", ncol(data), ncol(data_unique)), appendLF=T)
    }
    
    ###### parallel computing
    flag_parallel <- F
    if(parallel){
        flag_parallel <- dnet::dCheckParallel(multicores=multicores, verbose=verbose)
        if(flag_parallel){
            j <- 1
            res_list <- foreach::`%dopar%` (foreach::foreach(j=1:ncol(data_unique), .inorder=T), {
                progress_indicate(i=j, B=ncol(data_unique), 10, flag=T)
                suppressMessages(doReconstruct(x=data_unique[,j], Ntot=Ntot, Ntip=Ntip, Nnode=Nnode, connectivity=connectivity, e1=e1, e2=e2, output.detail=output.detail, verbose=verbose))
            })
            
            ###################################
            ## return back to input data matrix
            ind <- match(integer_vec, ind_unique)
            res_list <- res_list[ind]
            ###################################            
            
            if(!is.null(colnames(data))){
                names(res_list) <- colnames(data)
            }else{
                names(res_list) <- 1:ncol(data)
            }
            
            if(!output.detail){
                res <- do.call(base::cbind, res_list)
                if(is.numeric(data)){
                    res_tmp <- matrix(as.numeric(res), ncol=ncol(res), nrow=nrow(res))
                    if(!is.null(rownames(res))){
                        rownames(res_tmp) <- rownames(res)
                    }
                    if(!is.null(colnames(res))){
                        colnames(res_tmp) <- colnames(res)
                    }
                    res <- res_tmp
                }
            }else{
                res <- res_list
            }
            
        }
    }
    
    ###### non-parallel computing
    if(flag_parallel==F){
        res_list <- lapply(1:ncol(data_unique),function(j) {
            progress_indicate(i=j, B=ncol(data_unique), 10, flag=T)
            suppressMessages(doReconstruct(x=data_unique[,j], Ntot=Ntot, Ntip=Ntip, Nnode=Nnode, connectivity=connectivity, e1=e1, e2=e2, output.detail=output.detail, verbose=verbose))
        })
        
        ###################################
        ## return back to input data matrix
        ind <- match(integer_vec, ind_unique)
        res_list <- res_list[ind]
        ###################################  
        
        if(!is.null(colnames(data))){
            names(res_list) <- colnames(data)
        }else{
            names(res_list) <- 1:ncol(data)
        }
        
        if(!output.detail){
            res <- do.call(base::cbind, res_list)
            if(is.numeric(data)){
                res_tmp <- matrix(as.numeric(res), ncol=ncol(res), nrow=nrow(res))
                if(!is.null(rownames(res))){
                    rownames(res_tmp) <- rownames(res)
                }
                if(!is.null(colnames(res))){
                    colnames(res_tmp) <- colnames(res)
                }
                res <- res_tmp
            }
        }else{
            res <- res_list
        }
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message("", appendLF=T)
        message(paste(c("Finish at ",as.character(endT)), collapse=""), appendLF=T)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)
    
    invisible(res)
}
