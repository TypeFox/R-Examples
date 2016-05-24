
########################################################################
##get random s2n  matrix and observed s2n matrix
S2N <- function(A, class.labels, miR.labels, nperm) { 
     A <- A + 0.00000001

     N <- length(A[,1])
     Ns <- length(A[1,])

     subset.mask <- matrix(0, nrow=Ns, ncol=nperm)
     reshuffled.class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
     reshuffled.class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
     class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
     class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)

     order.matrix <- matrix(0, nrow = N, ncol = nperm)
     obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
     s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
     obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)

	 M1 <- matrix(0, nrow = N, ncol = nperm)
     M2 <- matrix(0, nrow = N, ncol = nperm)
     S1 <- matrix(0, nrow = N, ncol = nperm)
     S2 <- matrix(0, nrow = N, ncol = nperm)

     gc()

     C <- split(class.labels, class.labels)
     class1.size <- length(C[[1]])
     class2.size <- length(C[[2]])
     class1.index <- seq(1, class1.size, 1)
     class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)

     for (r in 1:nperm) {
        class1.subset <- sample(class1.index, size = ceiling(class1.size))
        class2.subset <- sample(class2.index, size = ceiling(class2.size))
        class1.subset.size <- length(class1.subset)
        class2.subset.size <- length(class2.subset)
        subset.class1 <- rep(0, class1.size)
        for (i in 1:class1.size) {
            if (is.element(class1.index[i], class1.subset)) {
                subset.class1[i] <- 1
            }
        }
        subset.class2 <- rep(0, class2.size)
        for (i in 1:class2.size) {
            if (is.element(class2.index[i], class2.subset)) {
                subset.class2[i] <- 1
            }
        }
        subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
        class1 <- class1.size/Ns
        class2 <- class2.size/Ns

        
        full.subset <- c(class1.subset, class2.subset)
        label1.subset <- sample(full.subset, size = Ns * class1)
        reshuffled.class.labels1[, r] <- rep(0, Ns)
        reshuffled.class.labels2[, r] <- rep(0, Ns)
        class.labels1[, r] <- rep(0, Ns)
        class.labels2[, r] <- rep(0, Ns)
        for (i in 1:Ns) {
            m1 <- sum(!is.na(match(label1.subset, i)))
            m2 <- sum(!is.na(match(full.subset, i)))
            reshuffled.class.labels1[i, r] <- m1
            reshuffled.class.labels2[i, r] <- m2 - m1
            if (i <= class1.size) {
                class.labels1[i, r] <- m2
                class.labels2[i, r] <- 0
           } else {
                class.labels1[i, r] <- 0
                class.labels2[i, r] <- m2
            }
        }
    }

# compute S2N for the random permutation matrix
     
     P <- reshuffled.class.labels1 * subset.mask
     n1 <- sum(P[,1])         
     M1 <- A %*% P
     M1 <- M1/n1      
     gc()
     A2 <- A*A        
     S1 <- A2 %*% P   
     S1 <- S1/n1 - M1*M1    
     S1 <- sqrt(abs((n1/(n1-1)) * S1))   
     gc()
     P <- reshuffled.class.labels2 * subset.mask
     n2 <- sum(P[,1])           
     M2 <- A %*% P           
     M2 <- M2/n2          
     gc()
     A2 <- A*A           
     S2 <- A2 %*% P      
     S2 <- S2/n2 - M2*M2 
     S2 <- sqrt(abs((n2/(n2-1)) * S2))
     rm(P)
     rm(A2)
     gc()

       
     S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
     S2 <- ifelse(S2 == 0, 0.2, S2)
     S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
     S1 <- ifelse(S1 == 0, 0.2, S1)
     gc()
     

     M1 <- M1 - M2
     rm(M2)
     gc()
     S1 <- S1 + S2
     rm(S2)
     gc()

     s2n.matrix <- M1/S1

   

    
# compute S2N for the "observed" permutation matrix

     P <- class.labels1 * subset.mask
     n1 <- sum(P[,1])         
     M1 <- A %*% P
     M1 <- M1/n1      
     gc()
     A2 <- A*A        
     S1 <- A2 %*% P   
     S1 <- S1/n1 - M1*M1    
     S1 <- sqrt(abs((n1/(n1-1)) * S1))   
     gc()
     P <- class.labels2 * subset.mask
     n2 <- sum(P[,1])           
     M2 <- A %*% P           
     M2 <- M2/n2          
     gc()
     A2 <- A*A           
     S2 <- A2 %*% P      
     S2 <- S2/n2 - M2*M2 
     S2 <- sqrt(abs((n2/(n2-1)) * S2))
     rm(P)
     rm(A2)
     gc()

      
     S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
     S2 <- ifelse(S2 == 0, 0.2, S2)
     S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
     S1 <- ifelse(S1 == 0, 0.2, S1)
     gc()
     

     M1 <- M1 - M2
     rm(M2)
     gc()
     S1 <- S1 + S2
     rm(S2)
     gc()

     obs.s2n.matrix <- M1/S1
     gc()

   
    return(list(s2n.matrix = s2n.matrix, 
                 obs.s2n.matrix = obs.s2n.matrix))
}