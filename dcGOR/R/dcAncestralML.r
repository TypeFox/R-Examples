#' Function to reconstruct ancestral discrete states using fast maximum likelihood algorithm
#'
#' \code{dcAncestralML} is supposed to reconstruct ancestral discrete states using fast maximum likelihood algorithm. It takes inputs both the phylo-formatted tree and discrete states in the tips. The algorithm assumes that state changes can be described by a probablistic reversible model. It first determines transition matrix between states (also considering branch lengths), then uses dynamic programming (from tips to the root) to estimate conditional maximum likelihood, and finally reconstructs the ancestral states (from the root to tips). If the ties occur at the root, the state at the root is set to the last state in ties (for example, usually being 'present' for 'present'-'absent' two states).
#'
#' @param data an input data matrix storing discrete states for tips (in rows) X characters (in columns). The rows in the matrix are for tips. If the row names do not exist, then addumedly they have the same order as in the tree tips. More wisely, users provide row names which can be matched to the tip labels of the tree. The row names can be more than found in the tree labels, and they should contain all those in the tree labels
#' @param phy an object of class 'phylo'
#' @param transition.model a character specifying the transition model. It can be: "different" for all-transition-different model (such as \eqn{matrix(c(0,1,2,0),2)}), "symmetric" for the symmetric model (such as \eqn{matrix(c(0,1,1,0),2)} or \eqn{matrix(c(0,1,2,1,0,3,2,3,0),3)}), "same" for all-transition-same model (such as \eqn{matrix(c(0,1,1,0),2)}), "customised" for the user-customised model (see the next parameter)
#' @param customised.model a matrix customised for the transition model. It can be: \eqn{matrix(c(0,1,1,0),2)}, \eqn{matrix(c(0,1,2,0),2)}, or \eqn{matrix(c(0,1,2,1,0,3,2,3,0),3)}
#' @param edge.length.power a non-negative value giving the exponent transformation of the branch lengths. It is useful when determining transition matrix between states
#' @param initial.estimate the initial value used for the maximum likelihood estimation
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
#'  \item{\code{transition}: an estimated transition matrix between states}
#'  \item{\code{relative}: a matrix of nodes X states, storing conditional maximum likelihood being relative to each state}
#' }
#' @note
#' This fast dynamic programming for ancestral discrete state reconstruction is partially inspired by a joint estimation procedure as described in \url{http://mbe.oxfordjournals.org/content/17/6/890.full}
#' @export
#' @seealso \code{\link{dcAncestralMP}}, \code{\link{dcDuplicated}}
#' @include dcAncestralML.r
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
#' res <- dcAncestralML(data, phy, parallel=FALSE)
#' res
#'
#' # 3) an input data matrix storing discrete states for tips (in rows) X only one character
#' data <- matrix(c(0,rep(1,3),rep(0,2)), ncol=1)
#' ## reconstruct ancestral states, with detailed output
#' res <- dcAncestralML(data, phy, parallel=FALSE, output.detail=TRUE)
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

dcAncestralML <- function(data, phy, transition.model=c("different","symmetric","same","customised"), customised.model=NULL, edge.length.power=1, initial.estimate=0.1, output.detail=F, parallel=T, multicores=NULL, verbose=T)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=T)
        message("", appendLF=T)
    }
    ####################################################################################

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    transition.model <- match.arg(transition.model)
    
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

    if (is.null(phy$edge.length)){
        stop("tree has no branch lengths")
    }
    EL <- phy$edge.length
    
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
    doReconstruct <- function(x, Ntot, Ntip, Nnode, E, e1, e2, output.detail, verbose){

        if (!is.factor(x)){
            x_tmp <- factor(x)
        }else{
            x_tmp <- x
        }
        nl <- nlevels(x_tmp)
        lvls <- levels(x_tmp)
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
    
    
        if(transition.model != 'customised'){
            rate <- matrix(NA, nl, nl)
            switch(transition.model, 
                same = np <- rate[] <- 1, 
                different = {
                    np <- nl * (nl - 1)
                    rate[col(rate) != row(rate)] <- 1:np
                },
                symmetric = {
                    np <- nl * (nl - 1)/2
                    sel <- col(rate) < row(rate)
                    rate[sel] <- 1:np
                    rate <- t(rate)
                    rate[sel] <- 1:np
                }
            )
        }else{
            if(is.null(customised.model)){
                stop("Customised transition model/matrix is not given!")
            }else{
                if (ncol(customised.model) != nrow(customised.model)) stop("Customised transition model/matrix is not square!")
                if (ncol(customised.model) != nl) stop("Customised transition model/matrix must have as many rows as the number of categories in 'x'")
                rate <- customised.model
                np <- max(rate)
            }
        }
    
        ################################################################################################
        if(verbose){
            message(sprintf("First, estimate the transition matrix based on the '%s' model and initial estimate '%.2f' (%s) ...", transition.model, initial.estimate, as.character(Sys.time())), appendLF=T)
        }
    
        ## determine transition matrix
        Q <- matrix(0, nl, nl)
    
        rate[cbind(1:nl, 1:nl)] <- 0
        rate[rate == 0] <- np + 1
        liks <- matrix(0, Ntip + Nnode, nl)
        TIPS <- 1:Ntip
        liks[cbind(TIPS, x_tmp)] <- 1
    
        ## the cost function (maximum likelihood of the whole tree) to optimise
        if(1){
            dev <- function(p) {
                if (any(is.nan(p)) || any(is.infinite(p))) return(1e+50)
                comp <- numeric(Ntip + Nnode)
                Q[] <- c(p, 0)[rate]
                diag(Q) <- -rowSums(Q)
                decompo <- eigen(Q)
                lambda <- decompo$values
                GAMMA <- decompo$vectors
                invGAMMA <- solve(GAMMA)
                for (i in seq(from=1, by=2, length.out=Nnode)) {
                    j <- i + 1L
                    anc <- e1[i]
                    des1 <- e2[i]
                    des2 <- e2[j]
                    v.l <- GAMMA %*% diag(exp(lambda * EL[i])) %*% invGAMMA %*% liks[des1, ]
                    v.r <- GAMMA %*% diag(exp(lambda * EL[j])) %*% invGAMMA %*% liks[des2, ]
                    v <- v.l * v.r
                    comp[anc] <- sum(v)
                    liks[anc, ] <- v/comp[anc]
                }

                dev <- -2 * sum(log(comp[-TIPS]))
                if (is.na(dev)){
                    Inf
                }else{
                    dev
                }
            }
        }else{
            dev <- function(p) {
                if (any(is.nan(p)) || any(is.infinite(p))) return(1e+50)
                comp <- numeric(Ntip + Nnode)
                Q[] <- c(p, 0)[rate]
                diag(Q) <- -rowSums(Q)
                for (i in seq(from=1, by=2, length.out=Nnode)) {
                    j <- i + 1L
                    anc <- e1[i]
                    des1 <- e2[i]
                    des2 <- e2[j]
                    v.l <- Matrix::expm(Q * EL[i]) %*% liks[des1, ]
                    v.r <- Matrix::expm(Q * EL[j]) %*% liks[des2, ]
                    v <- as.matrix(v.l * v.r)
                    comp[anc] <- sum(v)
                    liks[anc, ] <- v/comp[anc]
                }
                dev <- -2 * sum(log(comp[-TIPS]))
                if (is.na(dev)){
                    Inf
                }else{
                    dev
                }
            }
        }
    
        out <- stats::nlminb(rep(initial.estimate,length.out=np), function(p) dev(p), lower=rep(0, np), upper=rep(1/nl, np))
    
        ############################################################
        ## make sure the estimated transition rate no more than 1/nl
        ind <- which(out$par>=1/nl)
        out$par[ind] <- min(out$par)
        ############################################################
    
        ## transition matrix: from (in rows) -> to (in columns)
        t_matrix <- matrix(out$par[rate], nl, nl)
        diag(t_matrix) <- apply(t_matrix, 1, function(x) 1-sum(x,na.rm=T))
        colnames(t_matrix) <- rownames(t_matrix) <- lvls
    
        ################################################################################################
        if(verbose){
            message(sprintf("Second, estimate conditional maximum likelihood in a bottom-up manner (%s) ...", as.character(Sys.time())), appendLF=T)
        }
    
        ## Lx(s): the likelihood of the best reconstruction of the subtree (under the node x) on the condition that the direct parent of the node x is assigned character state s
        ## Cx(s): the character state assigned to the node x in this optimal condition reconstruction
    
        ## dimension: #nodes X #states (direct parent)
        Cxs <- Lxs <- matrix(NA, nrow=Ntot, ncol=nl, dimnames=list(1:Ntot, lvls))
    
        ## for tips
        ind <- match(x, colnames(t_matrix)) 
        Lxs[1:Ntip, ] <- t(t_matrix[, ind])
        Cxs[1:Ntip, ] <- matrix(rep(x, nl), ncol=nl)

        ## for internal nodes (in a postordered manner)
        for (i in seq(from=1, by=2, length.out=Nnode)) {
        
            j <- i + 1L
            cur <- e1[i]
            des1 <- e2[i]
            des2 <- e2[j]
        
            ### do calculation
            for(k in 1:nl){
                tmp <- t_matrix[k,] * (Lxs[des1,] * Lxs[des2,])
                Lxs[cur, k] <- max(tmp)
                ## for ties, always use the last state in ties
                ind <- which(tmp==max(tmp))
                Cxs[cur, k] <- names(ind[length(ind)])
            }
        }

        ## (conditional) maximum likelihood
        cml <- t(apply(Lxs, 1, function(x) x/sum(x)))

        ################################################################################################
        if(verbose){
            message(sprintf("Finally, reconstruct ancestral states in a top-down manner (%s) ...", as.character(Sys.time())), appendLF=T)
        }
    
        ## reconstruct the ancestral states in a top-down manner
        anc <- rep(NA, Ntot)
        names(anc) <- 1:Ntot
    
        ## in a preordered manner
        for (i in seq(to=1, by=-2, length.out=Nnode)) {
        
            ### for the root
            if(i == 2*Nnode-1){
                cur <- e1[i]
                ind <- which(Lxs[cur,]==max(Lxs[cur,]))
            
                if(length(ind)==1){
                    anc[cur] <- colnames(Lxs)[ind]
                }else{
                    ### in ties, always the last state in ties
                    anc[cur] <- colnames(Lxs)[ind[length(ind)]]
                }
            }
        
            j <- i + 1L
            cur <- e1[i]
            des1 <- e2[i]
            des2 <- e2[j]
        
            ind <- match(anc[cur], colnames(Cxs))
            anc[des1] <- Cxs[des1,ind]
            anc[des2] <- Cxs[des2,ind]
        }
    
        ####################################################################################
        if(verbose){
            ## A summary of changes between states
            p2c <- cbind(anc[phy$edge[,1]], anc[phy$edge[,2]])
            ### for the root: being present
            if(nl==2){
                if(anc[e1[2*Nnode-1]] == lvls[2]){
                    p2c <- rbind(p2c, lvls)
                }
            }
            all <- paste(p2c[,1], "->", p2c[,2], sep='')
            changes <- sapply(unique(all), function(x){
                sum(x==all)
            })
            changes <- sort(changes)
            msg <- paste(names(changes),changes, sep=": ", collapse="\n")
    
            message(sprintf("In summary, the number of between-state changes:\n%s", msg), appendLF=T)
        }
    
        if(output.detail){
            res <- list()
            res$states <- anc
            res$transition <- signif(t_matrix, digits=4)
            res$relative <- signif(cml, digits=4)
        }else{
            res <- anc
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
                suppressMessages(doReconstruct(x=data_unique[,j], Ntot=Ntot, Ntip=Ntip, Nnode=Nnode, E=E, e1=e1, e2=e2, output.detail=output.detail, verbose=verbose))
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
            suppressMessages(doReconstruct(x=data_unique[,j], Ntot=Ntot, Ntip=Ntip, Nnode=Nnode, E=E, e1=e1, e2=e2, output.detail=output.detail, verbose=verbose))
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
