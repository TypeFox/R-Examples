Altscreen<-function(nfac, randomize=FALSE) {
 if (nfac<6) {stop("Alternate screening designs exist only for 6, 7 or 8 factors") }
 if (nfac>8) {stop("Alternate screening designs exist only for 6, 7 or 8 factors") }
 if (nfac==6) {v <- c(4,  6, 11, 13, 19, 24, 25, 30, 34, 37, 44, 47, 49, 55, 58, 64)
               Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1))
               AS <- Full[v, ]
               v2 <- c(16,  1,  4, 13,  6, 11, 10,  7,  8,  9, 12,  5,  2, 15, 14,  3)
               AS <- AS[v2, ]
               rownames(AS)<- c(1:16)
               if (randomize==TRUE) {AS <- AS[sample(1:16), ]}
               return(AS)        
             }
 if (nfac==7) {v <- c(8,  13,  19,  28,  33,  46,  50,  63,  71,  74,  86,  89, 100, 107, 117, 128)
               Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1),G=c(-1,1))
               AS <- Full[v, ]
               v2 <- c(16,  1,  4, 13,  6, 11, 10,  7,  8,  9, 14,  3,  2, 15, 12,  5)
               AS <- AS[v2, ]
               rownames(AS)<- c(1:16)
               if (randomize==TRUE) {AS <- AS[sample(1:16), ]}
               return(AS)    
             }
 if (nfac==8) {v <- c(16,  23,  52,  57,  75,  86, 101, 106, 141, 154, 163, 166, 196, 209, 255, 256)
               Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1),G=c(-1,1),H=c(-1,1))
               AS <- Full[v, ]
               v2 <- c(16,  1,  3, 13,  6, 12, 10,  8, 15,  2,  5, 11,  9,  7,  4, 14)
               AS <- AS[v2, ]
               rownames(AS)<- c(1:16)
               if (randomize==TRUE) {AS <- AS[sample(1:16), ]}
               return(AS)        
             }
}



