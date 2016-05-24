epi.nomogram <- function(se, sp, lr, pre.pos, verbose = FALSE){
    # If likelihood ratios are known:
    if(is.na(se) & is.na(sp) & !is.na(lr[1])& !is.na(lr[2])){
       lr.pos <- lr[1]
       lr.neg <- lr[2]
       }

    # If likelihood ratios are not known:
    if(!is.na(se) & !is.na(sp) & is.na(lr[1]) & is.na(lr[2])){
       # se <- ifelse(se == 1.0, 1 - 1E-04, se)
       # sp <- ifelse(sp == 1.0, 1 - 1E-04, sp)
       lr.pos <- se / (1 - sp)
       lr.neg <- (1 - se) / sp
       }
       
    pre.odds <- pre.pos / (1 - pre.pos)
    post.odds.pos <- pre.odds * lr.pos
    post.odds.neg <- pre.odds * lr.neg
   
    post.prob.pos <- post.odds.pos / (1 + post.odds.pos)
    post.prob.neg <- post.odds.neg / (1 + post.odds.neg)

    lr <- as.data.frame(cbind(pos = lr.pos, neg = lr.neg))
    prob <- as.data.frame(cbind(pre.pos = pre.pos, post.pos = post.prob.pos, post.neg = post.prob.neg))
    rval <- list(lr = lr, prob = prob)
         
   if(verbose == TRUE){
     return(rval)
     }
   
   if(verbose == FALSE){
     post.prob.pos <- ifelse(post.prob.pos < 0.01, round(post.prob.pos, digits = 4), round(post.prob.pos, digits = 2))
     post.prob.neg <- ifelse(post.prob.neg < 0.01, round(post.prob.neg, digits = 4), round(post.prob.neg, digits = 2))
     cat("The post-test probability of being disease positive is", post.prob.pos, "\n")
     cat("The post-test probability of being disease negative is", post.prob.neg, "\n") 
     }  
}