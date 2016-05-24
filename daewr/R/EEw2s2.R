EEw2s2 <- function(des='',randomize=FALSE) {
  if (des=='')  {
cat(" ", "\n")
cat("Catalog of D-efficient Estimation Equivalent RS","\n")
cat("  Designs for (2 wp factors and  2 sp factors)  ","\n")
cat(" ", "\n")
cat("   Jones and Goos, JQT(2012) pp. 363-374","\n")
cat(" ", "\n")
cat(format("Design Name",width=11),format("whole plots",width=11),format("sub-plots/whole plot",width=21),"\n")
cat("----------------------------------------","\n")
cat(format("EE16R8WP", width=11),format("   8",width=11),format("           2",width=21),"\n")
cat(format("EE18R9WP", width=11),format("   9",width=11),format("           2",width=21),"\n")
cat(format("EE20R10WP", width=11),format("  10",width=11),format("           2",width=21),"\n")
cat(format("EE21R7WP", width=11),format("   7",width=11),format("           3",width=21),"\n")
cat(format("EE22R11WP", width=11),format("  11",width=11),format("           2",width=21),"\n")
cat(format("EE24R8WP", width=11),format("   8",width=11),format("           3",width=21),"\n")
cat(format("EE27R9WP", width=11),format("   9",width=11),format("           3",width=21),"\n")
cat(format("EE28R7WP", width=11),format("   7",width=11),format("           4",width=21),"\n")
cat(format("EE30R10WP", width=11),format("  10",width=11),format("           3",width=21),"\n")
cat(format("EE32R8WP", width=11),format("   8",width=11),format("           4",width=21),"\n")
cat(format("EE33R11WP", width=11),format("  11",width=11),format("           3",width=21),"\n")
cat(format("EE35R7WP", width=11),format("   7",width=11),format("           5",width=21),"\n")
cat(format("EE36R9WP", width=11),format("   9",width=11),format("           4",width=21),"\n")
cat(format("EE40R8WP", width=11),format("   8",width=11),format("           5",width=21),"\n")
cat(format("EE40R10WP", width=11),format("  10",width=11),format("           4",width=21),"\n")
cat(format("EE42R7WP", width=11),format("   7",width=11),format("           6",width=21),"\n")
cat(format("EE44R11WP", width=11),format("  11",width=11),format("           4",width=21),"\n")
cat(format("EE45R9WP", width=11),format("   9",width=11),format("           5",width=21),"\n")
cat(format("EE48R8WP", width=11),format("   8",width=11),format("           6",width=21),"\n")
cat(format("EE50R10WP", width=11),format("  10",width=11),format("           5",width=21),"\n")
cat(format("EE54R9WP", width=11),format("   9",width=11),format("           6",width=21),"\n")
cat(format("EE55R11WP", width=11),format("  11",width=11),format("           5",width=21),"\n")
cat(format("EE60R10WP", width=11),format("  10",width=11),format("           6",width=21),"\n")
cat(format("EE66R11WP", width=11),format("  11",width=11),format("           6",width=21),"\n")
cat(" ","\n")
cat("==> to retrieve a design type EEw2s2('EE21R7WP') etc.","\n")

         } else if (des=='EE16R8WP') {v <- c(145, 433, 202, 490, 643,  67,  52, 628, 589,  13, 454, 166, 103, 391, 320, 176)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:16)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:8),each=2))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:16)
                                                           #this randomizes subplots
                                                           s<-2
                                                           w<-8
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)
         } else if(des=='EE18R9WP') { v <- c(181, 505, 164, 488,   3, 651, 544, 220, 608, 284,  24, 672,  79, 727,  62, 710, 504, 342)
                                     Full <-expand.grid(WP=c(1:9),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:18)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:9),each=2))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:18)
                                                           #this randomizes subplots
                                                           s<-2
                                                           w<-9
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)
         } else if(des=='EE20R10WP') { v <- c(1, 721, 622, 262,  23, 743,  64, 784, 185, 545, 566, 206, 287, 647,  88, 808, 609, 249, 590, 410)
                                     Full <-expand.grid(WP=c(1:10),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:20)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:10),each=2))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:20)
                                                           #this randomizes subplots
                                                           s<-2
                                                           w<-10
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)                                   
         } else if(des=='EE21R7WP') {v <- c(456,  15, 141,  65, 506, 380, 143, 458,  17, 123, 438, 564, 404,  26, 341, 237, 489, 174, 308, 560,  56)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE22R11WP') {v <- c(595, 199,   2, 794, 267, 663, 818,  26,  71, 863, 886,  94, 436, 634, 690, 294, 625, 229, 219, 417, 341, 737)
                                     Full <-expand.grid(WP=c(1:11),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:22)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:11),each=2))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:22)
                                                           #this randomizes subplots
                                                           s<-2
                                                           w<-11
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)           
         } else if(des=='EE24R8WP') {v <- c(49, 193, 553, 578, 434,  74, 643, 139, 499, 436, 580,  76, 525, 237, 165, 198,  54, 558,  15, 519, 375, 480, 408,  48)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:24)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:8),each=3))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:24)
                                                           #this randomizes subplots
                                                           s<-3
                                                           w<-8
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE27R9WP') {WP<-c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9)
                                     w1<-c(1.0000,  1.0000,  1.0000,  0.0000,  0.0000,  0.0000,  0.7805,  0.7805,  0.7805,  1.0000,  1.0000,  1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000,
                                           -1.0000, -1.0000, -1.0000,  0.0000,  0.0000,  0.0000,  1.0000,  1.0000,  1.0000)
                                     w2<-c(-1.0000, -1.0000, -1.0000,  1.0000,  1.0000,  1.0000,  0.0000,  0.0000,  0.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000,
                                            1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  0.4752,  0.4752,  0.4752)
                                     s1<-c(-1,  1,  0,  1,  0, -1, -1,  0,  1,  0,  1, -1,  0, -1,  1, -1,  0,  1, -1,  1,  0,  0, -1,  1, -1,  0,  1)
                                     s2<-c(-1.000,  0.000,  1.000,  1.000, -1.000,  1.000, -1.000,  0.000,  0.000,  1.000,  0.000, -1.000,  1.000,  0.000, -1.000,  1.000, -1.000,  0.000, -1.000, -1.000,
                                            1.000, -1.000,  1.000,  1.000, -1.000,  0.827,  0.000)
                                     EE<-data.frame(WP,w1,w2,s1,s2)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:9),each=3))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:27)
                                                           #this randomizes subplots
                                                           s<-3
                                                           w<-9
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE28R7WP') {WP<-c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7)
                                     w1<-c(1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0)
                                     w2<-c(1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0, -1, -1, -1, -1,  1,  1,  1,  1)
                                     s1<-c(1.0000, -1.0000,  0.0000,  1.0000,  1.0000, -1.0000,  1.0000, -1.0000,  1.0000, -1.0000, -1.0000,  1.0000, -1.0000, -1.0000,  1.0000,  0.0000,  1.0000, -1.0000,
                                           1.0000,  0.0000, -1.0000, -1.0000,  1.0000,  0.0000, -1.0000, -1.0000,  0.0555,  1.0000)
                                     s2<-c(-1, -1,  1,  0, -1, -1,  1,  1,  1,  1, -1, -1,  1, -1,  0, -1, -1,  1,  1,  0, -1,  0, -1,  1,  0,  1, -1,  1)
                                     EE<-data.frame(WP,w1,w2,s1,s2)
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
         } else if(des=='EE30R10WP') {WP<-c(1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,  5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9,  9, 10, 10, 10)
                                     w1<-c(0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  0, -1, -1, -1,  1,  1,  1, -1, -1, -1,  0,  0,  0, -1, -1, -1,  1,  1,  1)
                                     w2<-c(0.0000,  0.0000,  0.0000, -1.0000, -1.0000, -1.0000,  1.0000,  1.0000,  1.0000,  0.0000,  0.0000,  0.0000,  1.0000,
                                           1.0000,  1.0000,  1.0000,  1.0000,  1.0000, -1.0000, -1.0000, -1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,
                                           1.0000, -0.2242, -0.2242, -0.2242)
                                     s1<-c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  1,  1,  1, -1, -1, -1, -1, -1, -1,  1)
                                     s2<-c(0,  0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1, -1,  0, -1,  1, -1, -1,  1, -1)
                                     EE<-data.frame(WP,w1,w2,s1,s2)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:10),each=3))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:30)
                                                           #this randomizes subplots
                                                           s<-3
                                                           w<-10
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
    
    
         } else if(des=='EE32R8WP') {v <- c(433, 145,   1, 577, 194, 626, 482,  50, 499, 211, 643,  67, 596, 164, 452,  20, 101, 605, 317, 245, 374, 518, 302,  14, 335, 119,
                                           623, 263,  64, 568, 352, 424)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:32)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:8),each=4))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:32)
                                                           #this randomizes subplots
                                                           s<-4
                                                           w<-8
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE30R10WP') {WP<-c(1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,  5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9,  9, 10, 10, 10)
                                     w1<-c(0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  0, -1, -1, -1,  1,  1,  1, -1, -1, -1,  0,  0,  0, -1, -1, -1,  1,  1,  1)
                                     w2<-c(0.0000,  0.0000,  0.0000, -1.0000, -1.0000, -1.0000,  1.0000,  1.0000,  1.0000,  0.0000,  0.0000,  0.0000,  1.0000,
                                           1.0000,  1.0000,  1.0000,  1.0000,  1.0000, -1.0000, -1.0000, -1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,
                                           1.0000, -0.2242, -0.2242, -0.2242)
                                     s1<-c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  1,  1,  1, -1, -1, -1, -1, -1, -1,  1)
                                     s2<-c(0,  0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1, -1,  0, -1,  1, -1, -1,  1, -1)
                                     EE<-data.frame(WP,w1,w2,s1,s2)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:10),each=3))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:30)
                                                           #this randomizes subplots
                                                           s<-3
                                                           w<-10
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
    

         } else if(des=='EE33R11WP') {WP<-c(1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,  5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9,  9, 10, 10, 10, 11, 11, 11)
                                     w1<-c(0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1, -1, -1, -1,  0,  0,  0,  1,  1,  1,  0,  0,  0,  1,  1,  1,  0,  0,  0,  1,  1, 1)
                                     w2<-c(-1, -1, -1,  0,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0, -1, -1, -1,  0,  0,  0,  1,  1, 1)
                                     s1<-c(-0.24132, -0.74571,  0.85059,  0.77328,  0.15347, -0.94949, -0.97068,  0.26282,  0.57142,  0.66691,
                                            0.92762, -0.99334, -0.02421, -0.50774, -0.51018, -0.93534,  0.13127,  0.66763, -1.00000,  0.99705,
                                            0.99621,  0.63952,  0.33773, -1.00000,  0.68169, -0.99524,  0.91475, -0.51733, -0.50540,  1.00000,
                                           -0.86947, -0.86400,  0.49390)
                                     s2<-c(0.91061, -0.87649, -0.15561, -0.20217, -0.71439,  0.98894, -0.43842,  0.98889, -0.67196, -0.80857,
                                           0.99032,  0.45553, -1.00000,  1.00000,  0.99999,  0.09592, -1.00000,  0.78259, -0.98569, -0.95685,
                                          -0.96119, -0.01310, -0.83059,  0.91607, -0.80435,  0.44163,  1.00000,  0.96132, -0.12044, -0.76850,
                                           1.00000,  1.00000, -1.00000)
                                     EE<-data.frame(WP,w1,w2,s1,s2)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:11),each=3))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:33)
                                                           #this randomizes subplots
                                                           s<-3
                                                           w<-11
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 

         } else if(des=='EE35R7WP') {v <- c(379, 505, 127,   1, 253, 520, 142, 457, 205,  16,  45, 423, 297, 549, 171, 557, 494, 242, 179,  53, 516,  12,
                                            390, 327,  75, 440, 566,  62, 125, 377,  42, 168, 357, 420, 483)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE36R9WP') {v <- c(1, 487, 649, 163, 164, 650, 488,   2, 219, 543, 705,  57, 274, 679, 355, 112, 239, 725,  77, 563, 420,  15,
                                             339, 582,  70, 475, 637, 394, 674,  26, 512, 188, 135, 297, 702, 378)
                                     Full <-expand.grid(WP=c(1:9),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:36)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:9),each=4))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:36)
                                                           #this randomizes subplots
                                                           s<-4
                                                           w<-9
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE40R8WP') {v <- c(593, 161, 449, 305,  17,  50, 626, 482, 194, 338, 579, 147, 219, 435,  75, 620, 404,  44, 116, 476, 637, 277,
                                            61, 565, 205, 646, 142, 214, 286, 502, 231, 519, 159, 591,  15, 608, 464, 104, 392,  32)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:40)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:8),each=5))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:40)
                                                           #this randomizes subplots
                                                           s<-5
                                                           w<-8
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE40R10WP') {WP<-c(1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,
                                            8,  8,  8,  9,  9,  9,  9, 10, 10, 10, 10)
                                     w1<-c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                           -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1)
                                     w2<-c(1,  1,  1,  1,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0, -1, -1, -1, -1,  1,  1,  1,  1, -1,
                                           -1, -1, -1,  0,  0,  0,  0, -1, -1, -1, -1)
                                     s1<-c(1.00000, -0.77966, -0.76995,  1.00000, -0.99675, -0.99944, -0.30072,  0.97375, -0.31686,
                                          -0.65113, -0.12419,  1.00000, -0.31122, -0.06396,  0.97780, -0.69480, -0.99438,  0.97565,
                                          -0.99830, -0.30613,  1.00000, -0.99993, -0.99898,  0.99929,  0.98534,  0.97935, -0.10962,
                                           0.99028, -0.99895,  0.99941,  0.99988, -0.99996,  0.60424,  0.57595,  0.55075,  0.57399,
                                          -0.99975,  1.00000,  0.99929, -0.99916)
                                     s2<-c(-0.93578,  1.00000,  1.00000, -0.94037, -0.98676,  0.41478,  0.85200, -0.30739, -0.94247,
                                           -0.93255,  0.65211,  0.44627, -0.73337, -0.95381,  0.96526, -0.05473, -0.99411,  0.00400,
                                            0.98853, -0.02579, -0.99975,  1.00000, -0.99918,  0.99931, -1.00000, -1.00000,  1.00000,
                                           -0.99382, -0.99897,  0.99943, -0.99996,  0.99988,  0.09254,  0.12663,  0.99971,  1.00000,
                                            1.00000, -0.99992,  0.99931, -0.99901)
                                     EE<-data.frame(WP,w1,w2,s1,s2)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:10),each=4))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:40)
                                                           #this randomizes subplots
                                                           s<-4
                                                           w<-10
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE42R7WP') {v <- c(43, 421, 547, 232, 169, 106,   2, 380,   2, 443, 317, 128, 290,  38, 290, 416, 164, 542, 375, 438,  60, 501,
                                            186,  60, 432, 180, 558, 117,  54, 243, 391,  13, 202, 517, 139,  76, 525, 147, 210, 399,  21,  84)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
         } else if(des=='EE44R11WP') {WP<-c(1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8,  9,  9,  9,
                                            9, 10, 10, 10, 10, 11, 11, 11, 11)
                                     w1<-c(0,  0,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1,  1,  1,  1,  1,  0,  0,  0,
                                            0, -1, -1, -1, -1,  1,  1,  1,  1)
                                     w2<-c(1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,
                                           0,  1,  1,  1,  1, -1, -1, -1, -1)
                                     s1<-c(0.18532, -0.84429,  0.96867, -0.94704,  0.99218,  0.98611, -0.99838, -0.98064,  0.98252,  0.99579, -0.99042,
                                          -0.98862,  0.97823, -0.98030,  1.00000, -0.99866, -0.76197,  1.00000,  0.11848, -0.99385,  0.13043, -0.76756,
                                          -0.99414,  0.99393,  0.98086,  0.99739, -0.97939, -0.99959,  0.49499,  0.70847, -0.94145,  0.45942,  0.29916,
                                          -1.00000,  0.18823,  0.18850, -0.51531,  0.57717,  0.24671,  0.24820, -0.98461,  0.24222,  1.00000, -0.99714)
                                     s2<-c(0.77708,  0.90810,  0.01244, -0.99856,  0.90220, -0.98479, -0.99650,  0.99940,  0.91658, -0.98187, -1.00000,
                                           0.98560,  0.91392,  0.98818, -0.98345, -0.99834,  0.76129, -0.12240,  0.98649, -0.92632, -0.96962,  0.71739,
                                          -0.03351,  0.98480,  0.90914, -0.98216,  0.99276, -0.99943, -0.97279,  0.80294, -0.72137, -0.98347, -0.96515,
                                          -0.51817, -1.00000, -0.99406, -0.99344, -0.83605, -0.99950, -0.99656, -0.30913,  0.90542, -0.99123, -0.25718)
                                     EE<-data.frame(WP,w1,w2,s1,s2)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:11),each=4))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:44)
                                                           #this randomizes subplots
                                                           s<-4
                                                           w<-11
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE45R9WP') {v <- c(271, 595,  28, 676, 190, 218,  56, 542, 704, 380, 183, 669, 345, 507,  21, 616, 211, 697,  49, 292, 482, 563,
                                            644,  77, 239, 150,  69, 717, 555, 474,  97, 664,  16, 421, 502, 224, 710,  62, 386, 548,   9, 495, 414, 576,
                                            171)
                                     Full <-expand.grid(WP=c(1:9),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:45)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:9),each=5))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:45)
                                                           #this randomizes subplots
                                                           s<-5
                                                           w<-9
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE48R8WP') {v <- c(625,  49, 337, 481, 193, 193, 450,  18, 450, 162, 594, 306, 147,   3, 579, 147, 219, 507, 476,  44, 620, 188,
                                            260, 548, 445, 373, 517,  85, 373,  13,  30, 606, 462, 102, 174, 390, 207, 495, 351,  63, 639, 351, 648, 288,
                                            504,  72, 216, 576)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:48)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:8),each=6))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:48)
                                                           #this randomizes subplots
                                                           s<-6
                                                           w<-8
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE50R10WP') {v <- c(61, 421, 241, 601, 781, 482, 662, 212,  32, 572, 643, 283, 733, 193,  13, 134, 584, 764,  44, 494, 565, 475,
                                             205,  25, 655, 726, 186, 546, 276,  96, 787, 247,  67, 607, 427, 268, 628, 448,  88, 808, 119, 569, 749, 299,
                                             209,  90, 630, 270, 810, 450)
                                     Full <-expand.grid(WP=c(1:10),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:50)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:10),each=5))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:50)
                                                           #this randomizes subplots
                                                           s<-5
                                                           w<-10
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 

         } else if(des=='EE54R9WP') {v <- c(505, 667, 181,  19, 505, 343, 704, 137,  56, 299, 542, 218,   3, 570, 246, 165, 651, 165,  13, 661,  94, 499,
                                            418, 499, 545, 626,  59, 221,  59, 464, 213, 618,  51, 375, 699, 294, 565, 241,  79, 727, 565, 241,   8, 413,
                                            494, 170,  89, 656, 540, 702, 135,  54, 459, 378)
                                     Full <-expand.grid(WP=c(1:9),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:54)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:9),each=6))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:54)
                                                           #this randomizes subplots
                                                           s<-6
                                                           w<-9
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE55R11WP') {WP<-c(1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  6,  6,  6,  6,
                                            6,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11)
                                     w1<-c(1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1,
                                          -1,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1)
                                     w2<-c(-1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,
                                            1,  0, 0,  0,  0,  0,  1,  1,  1,  1,  1,  0,  0,  0, 0,  0, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0)
                                     s1<-c(0.98429, -0.82507, -0.98451,  0.99750,  0.04071, -0.87034, -0.99475, -0.01204,  0.98776,
                                           0.99778,  0.99995, -0.91951,  0.09887, -1.00000,  0.92910,  0.43305, -0.94788,  0.72799,
                                          -1.00000,  0.99976, -0.87261,  0.09876,  1.00000,  0.95365, -0.96688, -0.87350, -0.02119,
                                          -0.98713,  0.99174,  0.99849, -0.59677,  0.81325,  0.77006, -1.00000, -1.00000,  0.66729,
                                           1.00000,  1.00000,  0.93112,  0.14840,  0.00000, -1.00000,  1.00000,  1.00000, -1.00000,
                                           0.26815, -0.91270, -0.96677,  0.21395, -0.97235,  0.00000,  1.00000,  1.00000, -1.00000,
                                          -1.00000)
                                     s2<-c(0.91021,  0.97557, -0.98402, -0.99577, -0.02061,  0.88528, -0.94831, -0.02324, -0.99669,
                                           0.99857, -0.33466, -0.74973, -1.00000,  1.00000,  1.00000, -0.98809, -0.95946,  0.91483,
                                           1.00000, -0.08190, -1.00000, -0.12157,  0.99462, -0.92560,  0.93793,  0.87980,  0.02329,
                                          -0.97363, -0.99471,  0.98086, -0.91112, -0.94322, -0.95956,  0.88470,  0.87411,  0.74084,
                                          -0.95482, -0.95997, -1.00000,  0.78796, -1.00000,  0.00000,  1.00000, -1.00000,  1.00000,
                                          -0.89773,  0.88177,  0.89081, -0.88393,  0.88409, -1.00000, -1.00000,  1.00000,  0.00000,
                                           1.00000)
                                     EE<-data.frame(WP,w1,w2,s1,s2)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:11),each=5))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:55)
                                                           #this randomizes subplots
                                                           s<-5
                                                           w<-11
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE60R10WP') {WP<-c(1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  7, 7,  7,  7,  7,  7,  8,  8,  8,  8,  8,
                                             8,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10)
                                     w1<-c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                                            1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0)
                                     w2<-c(0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,
                                           1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0)
                                     s1<-c(-0.25288, -0.72584,  0.83486,  0.36873, -0.98195,  0.90698, -0.31684, -1.00000,  1.00000, -1.00000,  1.00000, -0.04849, -0.80187,  0.93467, -0.93572,
                                            0.99053,  0.15517, -0.99421, -0.92930,  1.00000,  0.99966, -0.98391, -0.99357,  1.00000, -0.73483,  0.10336,  0.99975, -0.97198,  0.94820, -0.99593,
                                            0.97861, -0.11000, -0.97978, -0.01820,  0.99953, -1.00000,  0.99962,  1.00000, -0.92972,  1.00000, -0.97702, -1.00000, -1.00000, -0.54625,  0.08232,
                                            0.48161,  0.98423, -0.95217, -1.00000, -0.20007,  0.85249,  1.00000, -1.00000,  1.00000,  0.96057, -1.00000,  0.99870, -1.00000, -0.09485,  0.00575)
                                     s2<-c(-0.71040, -0.79763,  0.54092,  0.86236,  0.90678, -0.99326, -0.38804, -0.79542,  0.29917,  0.95800, -0.98796,  1.00000,  0.29070,  0.75794, -0.99059,
                                           -0.20564, -1.00000,  0.99095,  0.99724, -0.99715,  0.97250,  0.07376, -0.99933, -1.00000,  0.90590, -0.09903, -0.99837,  0.12662,  0.90729, -0.99905,
                                            0.98608, -0.15295,  0.98263,  0.02622, -0.99723, -1.00000, -0.96872, -1.00000, -1.00000,  1.00000,  1.00000,  0.01574, -0.71426,  0.90832,  0.99987,
                                           -0.90413, -0.94444, -0.30111, -0.71224,  0.94995, -0.95141, -0.42612, -0.88513,  1.00000, -0.11086,  0.97592,  0.99732, -0.01926, -0.99837, -1.00000)
                                     EE<-data.frame(WP,w1,w2,s1,s2)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:10),each=6))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:60)
                                                           #this randomizes subplots
                                                           s<-6
                                                           w<-10
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE66R11WP') {WP<-c(1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,
                                             8,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11)
                                     w1<-c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0, 0,  0,  0,  0,
                                           0, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
                                     w2<-c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,
                                            1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1)
                                     s1<-c(-0.08536,  0.04859, -0.97256,  0.96469, -0.99999,  1.00000, -0.97123, -0.06605,  0.98627, -0.98369, -0.00895,  0.99901, -0.03921,  0.97856, -0.99437,
                                           -0.97307,  0.99560, -0.01215, -0.03334,  0.96042, -0.99434,  0.99813, -0.98866,  0.01316,  0.45629, -0.48804,  0.84232, -0.85490, -1.00000,  0.99969,
                                            0.04630,  0.07748, -0.31119, -0.19387,  0.46002, -0.31104, -1.00000,  0.70299,  1.00000,  0.38519,  1.00000, -0.96757, -0.79257, -0.76546, -0.88832,
                                            0.52015, -0.70245, -0.80991,  1.00000, -0.84982,  0.97308, -0.95245, -0.05020,  1.00000,  0.89930,  0.88143, -0.95362, -1.00000,  0.86108,  0.84376,
                                           -0.82466, -1.00000, -0.82337, -0.98441,  0.67493, -0.98599)
                                     s2<-c(-0.41563,  0.26767, -0.91084, -0.95726,  0.97809,  1.00000,  0.98496, -0.01963, -0.99834, -0.99104, -0.01267,  0.99875, -0.03929, 1.00000, -1.00000,
                                            0.98629, -0.98627,  0.00130, -0.99050, -0.03640, -0.99788,  0.00250,  0.98432,  0.99998, -0.22692,  0.99952,  0.99991, -0.99631, 0.12799, -0.94216,
                                            0.77437,  0.79950, -0.95416, -0.96031, -0.13838, -0.99194, -0.99358, -0.94996, -0.99574,  0.01866,  0.99006,  0.81201, -0.33748, 0.23941, -0.45839,
                                            0.06763, -0.02703, -0.29384, -1.00000,  0.08272,  0.76364, -0.96491,  1.00000, -1.00000,  0.81487,  0.79824, -0.76645, -0.57502,  0.76975, 0.79419,
                                            0.82274, -0.94014,  0.86408, -0.80296,  0.50191, -0.83181)
                                     EE<-data.frame(WP,w1,w2,s1,s2)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:11),each=6))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:66)
                                                           #this randomizes subplots
                                                           s<-6
                                                           w<-11
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 

          } else cat(" Design name misspelled-Enter EEw2s2( ) to display list of names","\n")
                 }


EEw2s2( )
EE<-EEw2s2('EE33R11WP')
EE<-EEw2s2('EE33R11WP', randomize = TRUE)