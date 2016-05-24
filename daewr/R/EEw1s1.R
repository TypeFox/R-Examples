EEw1s1 <- function(des='',randomize=FALSE) {
  if (des=='')  {
cat(" ", "\n")
cat("Catalog of D-efficient Estimation Equivalent RS","\n")
cat("  Designs for (1 wp factor and  1 sp factor)  ","\n")
cat(" ", "\n")
cat("   Jones and Goos, JQT(2012) pp. 363-374","\n")
cat(" ", "\n")
cat(format("Design Name",width=11),format("whole plots",width=11),format("sub-plots/whole plot",width=21),"\n")
cat("----------------------------------------","\n")
cat(format("EE8R4WP", width=11),format("   4",width=11),format("           2",width=21),"\n")
cat(format("EE10R5WP", width=11),format("   5",width=11),format("           2",width=21),"\n")
cat(format("EE12R4WP", width=11),format("   4",width=11),format("           3",width=21),"\n")
cat(format("EE12R6WP", width=11),format("   6",width=11),format("           2",width=21),"\n")
cat(format("EE14R7WP", width=11),format("   7",width=11),format("           2",width=21),"\n")
cat(format("EE15R5WP", width=11),format("   5",width=11),format("           3",width=21),"\n")
cat(format("EE16R4WP", width=11),format("   4",width=11),format("           4",width=21),"\n")
cat(format("EE18R6WP", width=11),format("   6",width=11),format("           3",width=21),"\n")
cat(format("EE20R4WP", width=11),format("   4",width=11),format("           5",width=21),"\n")
cat(format("EE20R5WP", width=11),format("   5",width=11),format("           4",width=21),"\n")
cat(format("EE21R7WP", width=11),format("   7",width=11),format("           3",width=21),"\n")
cat(format("EE24R4WP", width=11),format("   4",width=11),format("           6",width=21),"\n")
cat(format("EE24R6WP", width=11),format("   6",width=11),format("           4",width=21),"\n")
cat(format("EE25R5WP", width=11),format("   5",width=11),format("           5",width=21),"\n")
cat(format("EE28R7WP", width=11),format("   7",width=11),format("           4",width=21),"\n")
cat(format("EE30R5WP", width=11),format("   5",width=11),format("           6",width=21),"\n")
cat(format("EE30R6WP", width=11),format("   6",width=11),format("           5",width=21),"\n")
cat(format("EE35R7WP", width=11),format("   7",width=11),format("           5",width=21),"\n")
cat(format("EE36R6WP", width=11),format("   6",width=11),format("           6",width=21),"\n")
cat(format("EE42R7WP", width=11),format("   7",width=11),format("           6",width=21),"\n")
cat(" ","\n")
cat("==> to retrieve a design type EEw1s1('EE10R5WP') etc.","\n")

         } else if(des=='EE8R4WP') {v <- c( 1, 25, 34, 10, 31, 19, 28,  4)
                                     Full <-expand.grid(WP=c(1:4),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:8)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:4),each=2))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:8)
                                                           #this randomizes subplots
                                                           s<-2
                                                           w<-4
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                      return(EE)
         } else if(des=='EE10R5WP') {v <- c(41, 11, 42, 12, 33,  3,  9, 24, 35,  5)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:10)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:5),each=2))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:10)
                                                           #this randomizes subplots
                                                           s<-2
                                                           w<-5
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                      return(EE)

         } else if(des=='EE12R4WP') {v <- c(17,  5, 29, 26, 14,  2, 11, 35, 23, 12, 36, 24)
                                     Full <-expand.grid(WP=c(1:4),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:12)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:4),each=3))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:12)
                                                           #this randomizes subplots
                                                           s<-3
                                                           w<-4
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                      return(EE)
         } else if (des=='EE12R6WP') {v <- c(13, 49, 38,  2,  3, 39, 10, 28, 17, 53, 30, 12)
                                     Full <-expand.grid(WP=c(1:6),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:12)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:6),each=2))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:12)
                                                           #this randomizes subplots
                                                           s<-2
                                                           w<-6
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)
         } else if(des=='EE14R7WP') { v <- c(43,  1, 16, 58, 52, 31, 32, 53,  5, 47, 55, 34, 21, 63)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:14)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:7),each=2))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:14)
                                                           #this randomizes subplots
                                                           s<-2
                                                           w<-7
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)
         } else if(des=='EE15R5WP') { v <- c(26, 11, 41, 17,  2, 32, 28, 43, 13, 24,  9, 39, 35, 20,  5)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:15)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:5),each=3))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:15)
                                                           #this randomizes subplots
                                                           s<-3
                                                           w<-5
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)
         } else if(des=='EE16R4WP') { v <- c(25, 13,  1,  1, 18, 30,  6, 18, 11, 35, 11, 35,  4,  4, 28, 16)
                                     Full <-expand.grid(WP=c(1:4),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:16)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:4),each=4))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:16)
                                                           #this randomizes subplots
                                                           s<-4
                                                           w<-4
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)
         } else if(des=='EE18R6WP') { v <- c(49, 31, 13,  8, 26, 44, 39, 21,  3, 34, 16, 52, 41, 23,  5, 12, 30, 48)
                                     Full <-expand.grid(WP=c(1:6),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:18)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:6),each=3))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:18)
                                                           #this randomizes subplots
                                                           s<-3
                                                           w<-6
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)
                     
         } else if(des=='EE20R4WP') {v <- c(1, 25, 13,  1, 25, 10, 34, 10, 22, 34, 31, 31, 19,  7, 19,  4,  4, 28, 28, 16)
                                     Full <-expand.grid(WP=c(1:4),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:20)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:4),each=5))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:20)
                                                           #this randomizes subplots
                                                           s<-5
                                                           w<-4
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }                                    
                                     return(EE) 
         } else if(des=='EE20R5WP') {v <- c(11, 41, 11, 26, 42, 12, 27, 12,  3, 33, 18,  3, 34,  4,  4, 19, 10, 25, 25, 40)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:20)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:5),each=4))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:20)
                                                           #this randomizes subplots
                                                           s<-4
                                                           w<-5
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }                                    
                                     return(EE)   
         } else if(des=='EE21R7WP') {v <- c( 1, 43, 22, 16, 37, 58, 38, 59, 17, 46,  4, 25, 12, 33, 54, 34, 55, 13, 21, 42, 63)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:21)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:7),each=3))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:21)
                                                           #this randomizes subplots
                                                           s<-3
                                                           w<-7
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }                                    
                                     return(EE)                                            
         } else if(des=='EE24R4WP') {v <- c(29, 17,  5, 29, 17,  5, 26,  2, 14,  2, 26,  2, 23, 35, 11, 23, 35, 11, 24, 12, 12, 24, 36, 36)
                                     Full <-expand.grid(WP=c(1:4),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:24)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:4),each=6))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:24)
                                                           #this randomizes subplots
                                                           s<-6
                                                           w<-4
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE24R6WP') {WP<-c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6)
                                     w1<-c(1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1,  0,  0,  0,  0, -1, -1, -1, -1,  0,  0, 0,  0)
                                     s1<-c(-1.00, -1.00,  0.17, 1.00, -1.00,  1.00, -1.00,  0.14, -1.00,  1.00, -1.00,  0.17, -1.00, -0.05, -0.04,
                                            1.00,  1.00, -1.00, -1.00,  0.14, -0.05, -0.04, -1.00,  1.00)
                                     EE<-data.frame(WP,w1,s1)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:6),each=4))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:24)
                                                           #this randomizes subplots
                                                           s<-4
                                                           w<-6
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)            
         } else if(des=='EE25R5WP') {v <- c(31, 16,  1, 31,  1, 27, 42, 12, 12, 42, 43, 43, 13, 13, 28, 34,  4, 19,  4, 34, 40, 25, 10, 25, 40)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:25)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:5),each=5))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:25)
                                                           #this randomizes subplots
                                                           s<-5
                                                           w<-5
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE28R7WP') {v <- c( 1, 43,  1, 43, 16, 16, 58, 37, 17, 59, 17, 38,  4, 46, 46,  4, 61, 19, 40, 19, 13, 34, 55, 34, 56, 14, 35, 35)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:28)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:7),each=4))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:28)
                                                           #this randomizes subplots
                                                           s<-4
                                                           w<-7
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)     
         } else if(des=='EE30R5WP') {v <- c(11, 41, 41, 11, 26, 26, 32,  2, 32,  2, 32, 17,  8, 38,  8, 38, 23, 23, 44, 14, 14, 44, 29, 29,  5, 35,  5, 35, 20, 35)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:30)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:5),each=6))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:30)
                                                           #this randomizes subplots
                                                           s<-6
                                                           w<-5
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)        
         } else if(des=='EE30R6WP') {v <- c(37, 19,  1, 37,  1, 14, 32, 14, 50, 50, 39, 39,  3, 21,  3, 46, 46, 10, 28, 28, 53, 17, 35, 17, 53, 48, 12, 30, 30, 48)
                                     Full <-expand.grid(WP=c(1:6),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:30)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:6),each=5))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:30)
                                                           #this randomizes subplots
                                                           s<-5
                                                           w<-6
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE35R7WP') {v <- c(1, 43,  1, 43, 22, 16, 37, 58, 16, 58, 10, 52, 52, 31, 31, 46,  4,  4, 25, 46, 19, 61, 61, 19, 40, 55, 13, 55, 34, 34,
                                            7, 49, 49,  7, 28)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:35)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:7),each=5))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:35)
                                                           #this randomizes subplots
                                                           s<-5
                                                           w<-7
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)  
         } else if(des=='EE36R6WP') {v <- c( 1, 37, 37,  1, 37, 19, 50, 50, 14, 32, 50, 14, 39, 39,  3, 21, 39,  3, 16, 52, 16, 52, 52, 34, 11, 47, 11, 47, 29, 29,
                                             12, 30, 12, 48, 48, 30)
                                     Full <-expand.grid(WP=c(1:6),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:36)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:6),each=6))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:36)
                                                           #this randomizes subplots
                                                           s<-6
                                                           w<-6
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE42R7WP') {v <- c(57, 15, 57, 15, 36, 36, 23, 44,  2,  2, 44,  2, 17, 17, 59, 38, 38, 59,  4,  4, 46, 46,  4, 25, 19, 61, 19, 40, 61, 40,
                                            13, 55, 13, 34, 55, 34, 56, 14, 35, 35, 56, 14)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:42)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:7),each=6))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:42)
                                                           #this randomizes subplots
                                                           s<-6
                                                           w<-7
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 

          } else cat(" Design name misspelled-Enter EEw1s1( ) to display list of names","\n")
                 }

