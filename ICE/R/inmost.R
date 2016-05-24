"inmost" <-
function(data, eps=0.0001)
{ 
    lr <- data.frame(data)
    names(lr) <- c("l", "r")
    lr$l <- lr$l*(1 + eps)
    lr.stk <- stack(lr, select=c("l","r"))
    lr.order <- lr.stk[order(lr.stk[,1]),]
    n <- length(lr.order[,1])
    lr1 <- data.frame(lr.order[-n,1], lr.order[-1,1], paste(lr.order[-n,2],
    lr.order[-1,2]))
    names(lr1) <- c("q", "p", "label")
    innermost <- lr1[lr1$label=="l r",]
    innermost[,1] <- innermost[,1]/(1+eps)
    innermost
}

