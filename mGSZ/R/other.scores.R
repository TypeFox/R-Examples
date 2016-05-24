KS.score <-
function(expr.data, gene.sets){
    ord <- order(expr.data,decreasing=FALSE)
    expr.data <- expr.data[ord]
    gene.sets <- gene.sets[ord,]
    set.dim <- dim(gene.sets)
    result1 <- rep(0,set.dim[2])
    result2 <- result1
    profile2 <- matrix(0,set.dim[2],set.dim[1])
    expr.data.len <- length(expr.data)
    for(i in 1:set.dim[2]){
        po1 <- which(gene.sets[,i]==1)
        po0 <- which(gene.sets[,i]==0)
        tmp1 <- expr.data
        tmp1[po0] <- 0
        tmp0 <- expr.data
        tmp0[po1] <- 0
        tmp2 <- cumsum(abs(tmp1))/sum(abs(tmp1))
        tmp2b <- tmp1
        tmp2b[po1] <- 1
        tmp2c <- cumsum(tmp2b)/sum(tmp2b)
        tmp0b <- tmp0
        tmp0b[po0] <- 1
        tmp3  <- cumsum(abs(tmp0))/sum(abs(tmp0))
        tmp3b <- cumsum(abs(tmp0b))/sum(abs(tmp0b))
        result1[i]  = max(abs(tmp2c - tmp3b))
        result2[i] = max(abs(tmp2 - tmp3b))
        #profile2.abs <- abs(tmp2 - tmp3b)
        #profile1.abs <- abs(tmp2c - tmp3b)
        #profile2[i,] <- tmp2 - tmp3b
        # profile1 <- tmp2c - tmp3b
    }
    
    out <- list(KS.org = result1,KS.mod = result2)
    
    
}

#############

WRS.test.score <-
function(expr.data, gene.sets){
    ord <- order(expr.data)
    expr.data <- expr.data[ord]
    rnk <- rank(expr.data)
    gene.sets <- gene.sets[ord,]
    num.classes <- dim(gene.sets)[2]
    WRS.score <- rep(0,num.classes)
    for(i in 1:num.classes){
        pos1 <- which(gene.sets[,i]==1)
        num.cl.members <- length(pos1)
        pos0 <- which(gene.sets[,i]==0)
        num.non.cl.members <- length(pos0)
        R.cl.members <- sum(rnk[pos1])
        
        ##sum of ranks of non-class members
        R.non.cl.members <- sum(rnk[pos0])
        WRS.score[i] <- R.cl.members
       
        
    }
    
    return(WRS.score)
    
}

##############

SS.test.score <-
function(expr.data, gene.sets){
    num.classes <- dim(gene.sets)[2]
    Mu <- mean(expr.data)
    SS.score = rep(0,num.classes)
    for(i in 1:num.classes){
        pos1 <- which(gene.sets[,i]==1)
        ##Calculation of sum of squares
        SS.score[i] <- sum((expr.data[pos1]-Mu)^2)
        
    }
    return(SS.score)
}

################

sumTestscore <-
function(expr.data, gene.sets){
    num.classes <- dim(gene.sets)[2]
    Mu <- mean(expr.data)
    sum.score = rep(0,num.classes)
    for(i in 1:num.classes){
        pos1 <- which(gene.sets[,i]==1)
        ##Calculation of sum of squares
        sum.score[i] <- sum(expr.data[pos1]-Mu)
        
    }
    return(sum.score)
}

##################