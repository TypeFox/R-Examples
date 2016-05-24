shrinkBeta <- function(Y, A, JJs, q){

mini <- min(unlist(JJs))
kj <- unlist(lapply(JJs, length)) + 1

## Generate set of all possible inequalities
## (".col" stands for "collapsed", 
##  ".exp" for "explicit")
JJs.exp <- matrix(ncol = 3)
for (u in 1:length(JJs)){
     JJ <- JJs[[u]]
     JJs.exp <- rbind(JJs.exp, c(u, NA, JJ[1]))
     
     if (length(JJ) > 1){for (v in 1:(length(JJ) - 1)){JJs.exp <- rbind(JJs.exp, c(u, JJ[v:(v + 1)]))}}
}

JJs.exp <- matrix(JJs.exp[-1, ], ncol = 3)
JJs.exp <- cbind(1:length(JJs.exp[, 1]), JJs.exp)

act <- rep(0, q)
act[A] <- 1
JJs.A <- matrix(cbind(JJs.exp, act), ncol = 5)

JJs.exp <- matrix(JJs.A[JJs.A[, 5] == 1, -5], ncol = 4)
dimnames(JJs.exp) <- list(NULL, c("A", "fac", "", ""))
dimnames(JJs.A) <- list(NULL, c("A", "fac", "", "", "active"))

facs <- sort(unique(JJs.exp[, 2]))
facs <- 1:length(kj)
rems <- NULL  # collect columns to be removed
sums <- NULL  # collect columns to be merged

## go through all possible cases for each factor
## record actions
## compute which columns of Y must be merged

for (u in facs){

    JJs.A.u <- matrix(JJs.A[JJs.A[, 2] == u, ], ncol = 5)
    tmp <- matrix(JJs.exp[JJs.exp[, 2] == u, -1], ncol = 3)

    ndummy <- kj[u] - 1
    nactive <- sum(JJs.A.u[, 5])

    if (nactive >= 1){

    ## -------------------------------------------
    ## the current factor has only two levels 
    ## (totally, parametrized by only one dummy variable)
    ## and this single constraint is active
    if (ndummy == 1){rems <- rbind(rems, tmp)} 
    ## -------------------------------------------

    if (ndummy > 1){
    ## -------------------------------------------
    ## beta_{j, 2} >= 0
    if (is.na(tmp[1, 2]) == FALSE){
        sum1 <- matrix(tmp[1, ], ncol = 3)        
                
        if (nactive >=2){
           for (i in 2:length(tmp[, 1])){
               if (tmp[i, 2] == sum1[length(sum1[, 1]), 3]){
                   sum1[length(sum1[, 1]), 3] <- tmp[i, 3]} else {
                   sum1 <- rbind(sum1, tmp[i, ])} # end if
           } # end for i  
        } # end nactive >= 2   
          
        sums <- rbind(sums, sum1)

    } # end if    
    ## -------------------------------------------
        
    ## -------------------------------------------
    ## beta_{j, 2} < 0, beta_{j, 3} > beta_{j, 2}     
    if ((JJs.A.u[1, 5] == 1) && (JJs.A.u[2, 5] == 0)){
        rem2 <- matrix(tmp[1, ], ncol = 3)  
        rems <- rbind(rems, rem2)
                
        if ((ndummy >= 3) && (nactive >= 2)){  
            sum2 <- matrix(tmp[2, ], ncol = 3)    

            if (nactive >= 3){        
                for (i in 3:length(tmp[, 1])){
                    if (tmp[i, 2] == sum2[length(sum2[, 1]), 3]){
                        sum2[length(sum2[, 1]), 3] <- tmp[i, 3]} else {
                        sum2 <- rbind(sum2, tmp[i, ])} # end if
                } # end for i    
            }
            
            sums <- rbind(sums, sum2)  
         }
            
    } # end if           
    ## -------------------------------------------
    
    ## -------------------------------------------        
    ## beta_{j, 2} < 0, beta_{j, 2} > beta_{j, 3}
    if ((JJs.A.u[1, 5] == 1) && (JJs.A.u[2, 5] == 1)){
        rem3 <- NULL 
        sum3 <- NULL
           
        ## extract all beta_{j, l}'s that have to be set to 0
        tmp3 <- cbind(tmp, tmp[1, 3] + (1:length(tmp[, 1])) - 1)
        rem3 <- cbind(rem3, tmp3[tmp3[, 3] == tmp3[, 4], -4])       
        
        ## extract beta_{j, l}'s that have to be summed
        max.rem <- max(rem3[, 3]) 
        max.q <- max(JJs.A[, 4])
        tmp3.sum <- matrix(tmp3[tmp3[, 3] != tmp3[, 4], -4], ncol = 3)
        
        if ((nrow(tmp3.sum) > 0) && (max.q > max.rem) && (sum(JJs.A[JJs.A[, 3] > max.rem, 5], na.rm = TRUE) > 0)){
            sum3 <- matrix(tmp3.sum[1, ], ncol = 3) 
            
            if (nrow(tmp3.sum) > 1){
                for (i in 2:length(tmp3.sum[, 1])){
                    if (tmp3.sum[i, 2] == sum3[length(sum3[, 1]), 3]){
                        sum3[length(sum3[, 1]), 3] <- tmp3.sum[i, 3]} else {
                        sum3 <- rbind(sum3, tmp3.sum[i, ])} # end if
                } # end for i                    
            }  # end if length(tmp3.sum)

            sums <- rbind(sums, sum3)            
        } # end 
        
        rems <- rbind(rems, rem3)
    } # end if                   
    ## -------------------------------------------        

    } # end nactive >= 1

    } # end ndummy > 1
    
} # end for u


## complement sums
for (u in facs){
JJu <- JJs[[u]]
if (length(sums[sums[, 1] == u, ]) == 0){
    sums <- rbind(sums, cbind(u, JJu, JJu))} else {   
    tmp <- sort(as.vector(sums[sums[, 1] == u, -1]))
    
    miss <- (1 - JJu %in% tmp) * JJu
    miss <- miss[miss > 0]   
    
    if (length(miss) > 0){
        sums.u <- matrix(sums[sums[, 1] == u, ], ncol = 3)
        miss.tmp <- NULL
        for (h in 1:length(miss)){
            if (max((miss[h] >= sums.u[, 2]) * (miss[h] <= sums.u[, 3])) == 0){miss.tmp <- c(miss.tmp, miss[h])}
            } # end for h

    if (is.null(miss.tmp) == 0){
        sums <- rbind(sums, cbind(u, miss.tmp, miss.tmp))
        } # end if
    } # end if length(miss) > 0
    
} # end if else
} #end for u

## sort sums according to first two columns
sums <- matrix(sums[order(sums[, 1], sums[, 2]), ], ncol = 3)

for (u in facs){
## remove columns to be removed
if (is.null(rems) == FALSE){
    rem.col <- unique(as.vector(rems[rems[, 1] == u, 2:3]))
    rem.col <- sort(rem.col[is.na(rem.col) == FALSE])
    sums2 <- NULL     
    for (o in 1:length(sums[, 1])){
        if ((sums[o, 2] %in% rem.col) == FALSE){sums2 <- rbind(sums2, sums[o, ])}
    }
    sums <- sums2
    
    } # end if        
} # end for u

## now generate new design matrix
if (is.null(rems) == FALSE){
    rems.col <- unique(as.vector(rems[, 2:3]))
    rems.col <- rems.col[is.na(rems.col) == FALSE]
    }

if (mini == 1){Y.col <- Y[, 1] * NA} else {Y.col <- Y[, 1:(mini - 1)]}
for (s in 1:length(JJs)){
     if (s %in% sums[, 1]){
         sums.s <- matrix(sums[sums[, 1] == s, ], ncol = 3)
         for (d in 1:length(sums.s[, 1])){
              Y.col <- cbind(Y.col, apply(matrix(Y[, sums.s[d, 2]:sums.s[d, 3]], nrow = dim(Y)[1]), 1, sum))} # end for d
         } else { # end if
         
         add.col <- setminus(JJs[[s]], rems.col)
         if (length(add.col) > 0){
            Y.col <- cbind(Y.col, Y[, JJs[[s]]])}
         }
} # end for s
if(mini == 1){Y.col <- Y.col[, -1]} 
Y.col <- as.matrix(Y.col)

return(list("Y.col" = Y.col, "sums" = sums, "rems" = rems, "JJs.A" = JJs.A))
}













#
