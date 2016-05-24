resample<-function(A, n, fvar=1.5, ...)
{
   # if(nrow(A) != length(n)){ rep(n, length=nrow(A)) } 
   # stages <- rownames(A)
   a <- splitA(A, ...)
   # add dead fates to transitions so columns sum to 1
   a1 <- rbind(a$T, dead=1-colSums(a$T))
   #  add stage vector to transtions (2 extra rows )
   ## this makes it  easier to apply the number in stage 1  to the
   ## transitions in column 1 and so on
   a1 <- rbind(a1, n)
   r  <- nrow(A) + 2
   # sample transitions using sample sizes in stage vector (last row)
   # from  a multinomial distribtion  
   t1 <- apply(a1, 2, function(x) rmultinom(1, x[r], prob=x[-r]) /x[r])
   # remove dead fates in last row
   t1 <- t1[-nrow(t1),]
   ## get fertilities (remove zeros) and then
   f1 <- a$F[a$F > 0]
   ## fvar can be a vector, if so, then it should match the number of fertilities
   # sample fertilities from a log normal distribution
   # OCT 2011 - fixed bug noted by Tali Vardi - some elements may have both T and F values and were not added
   #t1[which(a$F > 0)] <- lnorms(1, f1, fvar)
    t1[which(a$F > 0)] <- t1[which(a$F > 0)] + lnorms(1, f1, fvar)
   # dimnames(t1) <- list(stages, stages)
   t1
}
