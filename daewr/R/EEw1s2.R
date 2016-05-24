EEw1s2 <- function(des='',randomize=FALSE) {
  if (des=='')  {
cat(" ", "\n")
cat("Catalog of D-efficient Estimation Equivalent RS","\n")
cat("  Designs for (1 wp factor and  2 sp factors)  ","\n")
cat(" ", "\n")
cat("   Jones and Goos, JQT(2012) pp. 363-374","\n")
cat(" ", "\n")
cat(format("Design Name",width=11),format("whole plots",width=11),format("sub-plots/whole plot",width=21),"\n")
cat("----------------------------------------","\n")
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
cat("==> to retrieve a design type EEw1s2('EE12R6WP') etc.","\n")

         } else if(des=='EE12R4WP') {WP<-c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
                                     w1<-c(-1, -1, -1,  1,  1,  1, -1, -1, -1,  0,  0,  0)
                                     s1<-c(0.0000,  1.0000, -1.0000,  1.0000,  1.0000, -1.0000, -1.0000,  0.0000,  1.0000, -1.0000, -1.0000,
                                           0.1524)
                                     s2<-c(-1,  1,  0,  1, -1,  0, -1,  1,  0, -1,  1,  0)
                                     EE<-data.frame(WP,w1,s1,s2)
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
         } else if (des=='EE12R6WP') {v <- c(145,   1,  38, 110,  45,  81,  52, 124,  17, 161,  30, 102)
                                     Full <-expand.grid(WP=c(1:6),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE14R7WP') { v <- c(57, 141, 184,  16, 157, 115, 116, 158,   5, 173,  48, 132,  91,  91)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE15R5WP') { v <- c(1,  76, 106,  37,  67,  97, 123,  18,  48,  44, 104, 134, 110,  80,   5)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE16R4WP') { v <- c(45,  33,  93,   9, 102,  78,  54,  30,  99,  27,   3,  75,  84,  24,  72,  12)
                                     Full <-expand.grid(WP=c(1:4),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE18R6WP') { v <- c(157,  67,  49, 134,  98,  26,  39, 147, 147, 100,  28, 136,   5, 113,  77, 138,  30, 102)
                                     Full <-expand.grid(WP=c(1:6),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
                     
         } else if(des=='EE20R4WP') {v <- c(33, 81, 69, 9, 93, 82, 106, 34, 46, 22, 91, 103, 43, 7, 31, 4,  52, 28, 100, 76)
                                     Full <-expand.grid(WP=c(1:4),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE20R5WP') {v <- c(86, 131, 101,  26,  92,  32, 122,   2,   8,  68,  98,  38, 134,  59,  44, 119,   5,  95, 125,  35)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE21R7WP') {v <- c(78, 162,  57,  58,  79, 163,  45, 171,  66, 172,  46,  67, 124,  40, 145,  13,  97, 139,  14, 140,  98)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE24R4WP') {v <- c(97, 1, 25, 25, 73, 49, 106, 34, 82, 10, 70, 22, 35, 11, 107, 35, 95, 47, 80, 20, 8, 104, 68, 56)
                                     Full <-expand.grid(WP=c(1:4),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE24R6WP') {v <- c(1, 37, 109, 145, 50, 14, 140, 104, 51, 159, 33, 69, 148, 40, 112, 4,  53, 161, 71, 35, 12, 102, 84, 120)
                                     Full <-expand.grid(WP=c(1:6),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:24)
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
         } else if(des=='EE25R5WP') {v <- c(131, 116, 11, 56, 41, 22, 82, 97, 112, 52, 33, 3, 63, 93, 123, 4, 124, 94, 34, 64, 90, 135, 15,
                                            105,  30)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE28R7WP') {v <- c(78, 162, 57, 162, 128, 170, 107, 23, 38, 122, 143, 164, 137, 158, 11, 95, 152, 68, 47, 173, 13, 97, 139,
                                            160,  98,  14, 161, 140)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE30R5WP') {v <- c(36, 66, 126, 36, 66, 6, 122, 2, 47, 92, 107, 32, 118, 43, 13, 58, 133, 103, 124, 4, 19, 94, 94,
                                            79, 135, 105,  15,  90,  30, 105)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE30R6WP') {WP<-c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6)
                                     w1<-c(1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000,  1.000, -1.000, -1.000, -1.000,
                                          -1.000, -1.000,  1.000,  1.000,  1.000,  1.000,  1.000, -0.064, -0.064, -0.064, -0.064, -0.064, -1.000,
                                          -1.000, -1.000, -1.000, -1.000)
                                     s1<-c(-1, -1,  0,  1,  1,  1, -1,  1, -1,  0, -1, -1,  0,  1,  1,  1,  1, -1, -1,  0,  0,  1,  1,  0, -1, -1,  1,  1,  0, -1)
                                     s2<-c(1, -1,  0, -1,  1, -1,  1,  1, -1,  0, -1,  1, -1,  0,  1, -1,  1, -1,  1,  0, -1, -1,  0,  1,  0, -1, -1,  1,  1,  0)
                                     EE<-data.frame(WP,w1,s1,s2)
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
         } else if(des=='EE35R7WP') {v <- c(22, 106, 1, 127, 169, 184, 16, 58, 142, 100, 3, 129, 24, 171, 108, 18, 60, 186, 144, 102, 75, 33, 54,
                                            75, 159, 146,  20, 188,  62, 104,  49, 154, 175,  70, 7)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE36R6WP') {v <- c(31, 157, 49, 121, 67, 157, 38, 20, 146, 110, 146, 56, 39, 111, 3, 147, 93, 129, 160, 16, 52, 124, 142,
                                            106,  47, 119,  11, 101,  11,  83,  18, 162,  54, 126, 144, 108)
                                     Full <-expand.grid(WP=c(1:6),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE42R7WP') {v <- c(15,  57, 183, 141,  36,  78,  58, 142, 121,  16, 163,  16,   3, 129, 108, 171,  24, 129,  53,  11, 179, 137,  95,
                                             95,  47, 173,   5, 152,  68, 131, 153,   6, 132,  69, 174,  48,  56,  14, 140, 182,  98,  98)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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

          } else cat(" Design name misspelled-Enter EEw1s2( ) to display list of names","\n")
                 }

EEw1s2( )
EE <- EEw1s2('EE12R4WP')
EE <- EEw1s2('EE12R4WP', randomize=TRUE)
EE <- EEw1s2('EE12R6WP')
EE <- EEw1s2('EE12R6WP', randomize=TRUE)
EE <- EEw1s2('EE14R7WP')
EE <- EEw1s2('EE14R7WP', randomize=TRUE)
EE <- EEw1s2('EE15R5WP')
EE <- EEw1s2('EE15R5WP', randomize=TRUE)
EE <- EEw1s2('EE16R4WP')
EE <- EEw1s2('EE16R4WP', randomize=TRUE)
EE <- EEw1s2('EE18R6WP')
EE <- EEw1s2('EE18R6WP', randomize=TRUE)
EE <- EEw1s2('EE20R4WP')
EE <- EEw1s2('EE20R4WP', randomize=TRUE)
EE <- EEw1s2('EE20R5WP')
EE <- EEw1s2('EE20R5WP', randomize=TRUE)
EE <- EEw1s2('EE21R7WP')
EE <- EEw1s2('EE21R7WP', randomize=TRUE)
EE <- EEw1s2('EE24R4WP')
EE <- EEw1s2('EE24R4WP', randomize=TRUE)
EE <- EEw1s2('EE24R6WP')
EE <- EEw1s2('EE24R6WP', randomize=TRUE)
EE <- EEw1s2('EE25R5WP')
EE <- EEw1s2('EE25R5WP', randomize=TRUE)
EE <- EEw1s2('EE28R7WP')
EE <- EEw1s2('EE28R7WP', randomize=TRUE)
EE <- EEw1s2('EE30R6WP')
EE <- EEw1s2('EE30R6WP', randomize=TRUE)
EE <- EEw1s2('EE30R5WP')
EE <- EEw1s2('EE30R5WP', randomize=TRUE)
EE <- EEw1s2('EE35R7WP')
EE <- EEw1s2('EE35R7WP', randomize=TRUE)
EE <- EEw1s2('EE36R6WP')
EE <- EEw1s2('EE36R6WP', randomize=TRUE)
EE <- EEw1s2('EE42R7WP')
EE <- EEw1s2('EE42R7WP', randomize=TRUE)





