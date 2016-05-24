###############################################################
# Chong Wu
# email: wuxx0845@umn.edu
###############################################################

# Find the corresponding leaves
correspLeaves <- function(tree) {
    tip.label <- tree$tip.label
    ntip <- length(tip.label)
    nbr <- nrow(tree$edge)
    edge = tree$edge
    edge2 = edge[, 2]
    
    res = vector("list", length = dim(edge)[1])
    for(i in 1:ntip) {
        node.loc <- which(edge2 == i)

        while (length(node.loc)) {
            res[[node.loc]] = c(res[[node.loc]] ,i)
            node <- edge[node.loc, 1]
            node.loc <- which(edge2 == node)
        }
    }
    res
}

#g.taxon.index = 1: weighted generalized taxon proportion
#g.taxon.index = 2: unweighted one
ranking <- function(y, X,tree,cov = NULL,gamma, g.taxon.index,model = "binomial") {
    
    ## calcualte cum and br.len (branch length)
    GuniF.cum = GUniFrac_cum(X,tree)
    cum = GuniF.cum$cum
    br.len = GuniF.cum$br.len
    br.len = as.matrix(br.len)
    
    # create weighted generalized taxa proportions
    if(g.taxon.index == 1) {
        cum2 = cum *matrix(rep(br.len,each = dim(cum)[2]),nrow = dim(cum)[1],ncol = dim(cum)[2],byrow = TRUE)
    } else {
        tmp.cum = cum
        tmp.cum[tmp.cum != 0] = 1
        # create unweighted generalized taxa proportions
        cum2 = tmp.cum *matrix(rep(br.len,each = dim(cum)[2]),nrow = dim(cum)[1],ncol = dim(cum)[2],byrow = TRUE)
    }
    
    Y = y
    X = t(cum2)

    n <- length(Y)
    k <- ncol(X)

    #### Score vector:
    if (is.null(cov)){
        ## NO nuisance parameters:
        r<-Y-mean(Y)
        U<-as.vector(t(X) %*% r)
    } else {
        tdat1 <-data.frame(trait=Y, cov)
        fit1 <-glm(trait~.,family=model,data=tdat1)
        pis <-fitted.values(fit1)
        r<-Y - pis
        
        U<-t(X) %*% r
    }
    U.score = abs(U)^gamma
    sum.U = sum(U.score)
    U.score2 = U.score/sum.U
}







