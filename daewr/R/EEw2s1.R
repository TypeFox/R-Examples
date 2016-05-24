EEw2s1 <- function(des='',randomize=FALSE) {
  if (des=='')  {
cat(" ", "\n")
cat("Catalog of D-efficient Estimation Equivalent RS","\n")
cat("  Designs for (2 wp factors and  1 sp factors)  ","\n")
cat(" ", "\n")
cat("   Jones and Goos, JQT(2012) pp. 363-374","\n")
cat(" ", "\n")
cat(format("Design Name",width=11),format("whole plots",width=11),format("sub-plots/whole plot",width=21),"\n")
cat("----------------------------------------","\n")
cat(format("EE14R7WP", width=11),format("   7",width=11),format("           2",width=21),"\n")
cat(format("EE16R8WP", width=11),format("   8",width=11),format("           2",width=21),"\n")
cat(format("EE18R9WP", width=11),format("   9",width=11),format("           2",width=21),"\n")
cat(format("EE20R10WP", width=11),format("  10",width=11),format("           2",width=21),"\n")
cat(format("EE21R7WP", width=11),format("   7",width=11),format("           3",width=21),"\n")
cat(format("EE22R11WP", width=11),format("  11",width=11),format("           2",width=21),"\n")
cat(format("EE24R8WP", width=11),format("   8",width=11),format("           3",width=21),"\n")
cat(format("EE24R12WP", width=11),format("  12",width=11),format("           2",width=21),"\n")
cat(format("EE27R9WP", width=11),format("   9",width=11),format("           3",width=21),"\n")
cat(format("EE28R7WP", width=11),format("   7",width=11),format("           4",width=21),"\n")
cat(format("EE30R10WP", width=11),format("  10",width=11),format("           3",width=21),"\n")
cat(format("EE32R8WP", width=11),format("   8",width=11),format("           4",width=21),"\n")
cat(format("EE33R11WP", width=11),format("  11",width=11),format("           3",width=21),"\n")
cat(format("EE35R7WP", width=11),format("   7",width=11),format("           5",width=21),"\n")
cat(format("EE36R9WP", width=11),format("   9",width=11),format("           4",width=21),"\n")
cat(format("EE36R12WP", width=11),format("  12",width=11),format("           3",width=21),"\n")
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
cat(format("EE60R12WP", width=11),format("  12",width=11),format("           5",width=21),"\n")
cat(format("EE66R11WP", width=11),format("  11",width=11),format("           6",width=21),"\n")
cat(format("EE72R12WP", width=11),format("  12",width=11),format("           6",width=21),"\n")
cat(" ","\n")
cat("==> to retrieve a design type EEw2s1('EE21R7WP') etc.","\n")
         } else if(des=='EE14R7WP') {WP<-c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
                                     w1<-c(-1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1,  0,  0)
                                     w2<-c(-1, -1,  1,  1, -1, -1,  1,  1,  0,  0,  1,  1, -1, -1)
                                     s1<-c(1.0000, -1.0000, -1.0000,  1.0000, -1.0000,  1.0000,  1.0000, -1.0000, -1.0000,  0.0000, -1.0000,  1.0000, -0.0736,
                                           1.0000)
                                     EE<-data.frame(WP,w1,w2,s1)
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

         } else if (des=='EE16R8WP') {v <- c(145,   1, 186,  42,  67, 211, 164,  20, 173,  29, 198,  54, 183, 111, 184, 112)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
         } else if(des=='EE18R9WP') {WP<-c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9)
                                     w1<-c(1,  1, -1, -1, -1, -1,  0,  0,  1,  1, -1, -1,  1,  1,  0,  0, -1, -1)
                                     w2<-c(1,  1,  1,  1, -1, -1,  1,  1, -1, -1, -1, -1,  0,  0,  1,  1,  1,  1)
                                     s1<-c(-1.0000,  1.0000, 1.0000, -1.0000,  1.0000, -1.0000,  1.0000,  0.0000, -1.0000,  1.0000, -1.0000,  1.0000, -1.0000,
                                            0.2096,  1.0000,  0.0000,  1.0000, -1.0000)
                                     EE<-data.frame(WP,w1,w2,s1)
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
         } else if(des=='EE20R10WP') { v <- c(1, 181, 182,   2,  83, 263,  64, 244,  85, 265,  36, 126,  27, 207, 258, 168, 209,  29, 260, 170)
                                     Full <-expand.grid(WP=c(1:10),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
         } else if(des=='EE21R7WP') {WP<-c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7)
                                     w1<-c(1,  1,  1, -1, -1, -1,  1,  1,  1,  0,  0,  0,  1,  1,  1,  0,  0,  0, -1, -1, -1)
                                     w2<-c(1.0000,  1.0000,  1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000,  1.0000,  1.0000,  1.0000,  0.0000,
                                           0.0000,  0.0000, -0.1586, -0.1586, -0.1586,  1.0000,  1.0000,  1.0000)
                                     s1<-c(1, -1,  0, -1,  0,  1,  1, -1,  0,  1, -1,  0,  1, -1,  0, -1,  1,  0,  1, -1,  0)
                                     EE<-data.frame(WP,w1,w2,s1)
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
         } else if(des=='EE22R11WP') {v <- c(1, 199, 123, 123,  69, 267,  26, 224,  93, 291,   6, 204, 117, 117, 195, 195, 141, 141, 263,  65, 286,  88)
                                     Full <-expand.grid(WP=c(1:11),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
         } else if(des=='EE24R8WP') {v <- c(161,  17,  89, 130,  58, 202, 147,   3,  75, 212,  68, 140, 197,  53, 125,  46, 190, 118,  15, 159,  87, 104,  32, 176)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
         } else if(des=='EE24R12WP') {v <- c(313,  97, 134, 134, 195, 195,  16, 232, 113, 113, 258,  42, 115, 115, 176, 176, 165, 165, 298,  82, 191, 191,  36, 252)
                                     Full <-expand.grid(WP=c(1:12),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:24)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:12),each=2))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:24)
                                                           #this randomizes subplots
                                                           s<-2
                                                           w<-12
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
 
         } else if(des=='EE27R9WP') {v <- c(235,  73, 154, 137,  56, 218, 165,   3,  84, 211,  49, 130, 185,  23, 104,   6, 168,  87,  34, 196, 115,  71, 233, 152,  18, 180, 99)
                                     Full <-expand.grid(WP=c(1:9),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:27)
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
         } else if(des=='EE28R7WP') {v <- c(57, 183, 183,  57, 142,  16, 142,  16,   3, 129,   3, 129, 172,  46, 172,  46,  12,  75,  75, 138, 167, 104, 104,  41,  28, 154,
                                            91,  91)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
         } else if(des=='EE30R10WP') {v <- c(91, 181,   1, 242,  62, 152,  83, 263, 173, 204,  24, 114,  25, 205, 115, 246,  66, 156,  17, 197, 107,  38, 218, 128, 269,  89,
                                             179,  50, 140, 230)
                                     Full <-expand.grid(WP=c(1:10),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:30)
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
    
    
         } else if(des=='EE32R8WP') {v <- c(145,   1, 145,   1, 194,  50,  50, 194,  19, 163, 163,  19, 100,  28, 100, 172, 157,  85,  85,  13, 214,  70, 214,  70, 119,  47, 191, 119, 136, 208, 136,  64)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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

         } else if(des=='EE33R11WP') {v <- c(23, 221, 122,  90, 288, 189,  36, 234, 135, 257,  59, 158,   5, 203, 104, 171, 270,  72, 227,  29, 128, 283,  85, 184, 295,  97,
                                             196,  10, 208, 109,  22, 121, 220)
                                     Full <-expand.grid(WP=c(1:11),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:33)
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

         } else if(des=='EE35R7WP') {WP<-c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7)
                                     w1<-c(1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
                                     w2<-c(1.0000,  1.0000,  1.0000,  1.0000,  1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000,  0.0000,  0.0000,  0.0000,
                                           0.0000,  0.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000,  1.0000,
                                           1.0000,  1.0000,  1.0000,  1.0000, -0.1213, -0.1213, -0.1213, -0.1213, -0.1213)
                                     s1<-c(-1, -1,  1,  0,  1, -1,  1, -1,  0,  1, -1,  1,  1, -1,  0, -1,  1, -1,  0,  1, -1,  1, -1,  0,  1,  1, -1,  1, -1,  0,  1,  1, -1,  0, -1)
                                     EE<-data.frame(WP,w1,w2,s1)
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

         } else if(des=='EE36R9WP') {v <- c(235,  73,  73, 235,  20, 182,  20, 182, 219,  57,  57, 219,   4,   4, 166, 166,  68, 230, 149, 149, 213, 132, 132,  51, 169,   7,
                         7, 169, 179,  98,  17,  98,  36, 117, 198, 117)
                                     Full <-expand.grid(WP=c(1:9),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
         } else if(des=='EE36R12WP') {v <- c(289, 181,  73, 242,  26, 134,  27, 243, 135, 100, 208, 316,  65, 173, 281, 294,  78, 186,   7, 223, 115, 104, 320, 212,   9, 225,
                                             117, 262, 154,  46,  95, 311, 203, 240,  24, 132)
                                     Full <-expand.grid(WP=c(1:12),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:36)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:12),each=3))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:36)
                                                           #this randomizes subplots
                                                           s<-3
                                                           w<-12
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)
         } else if(des=='EE40R8WP') {v <- c(161,  17, 161,  17,  89, 194,  50,  50, 194, 122,  43, 115,  43, 187, 187, 204, 204,  60,  60, 132, 213,  69,  69, 141, 213,  78,
                                            6,   6, 150, 150, 175,  31,  31, 103, 175,  16,  16,  88, 160, 160)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
         } else if(des=='EE40R10WP') {WP<-c(1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8,  9,  9,  9,
                                            9, 10, 10, 10, 10)
                                     w1<-c(-1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1,  0,  0,  0,  0, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1,
                                           -1, -1, -1, -1, -1)
                                     w2<-c(1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1,
                                           -1,  0,  0,  0,  0)
                                     s1<-c(1.0000, -1.0000,  0.0000,  0.0000,  1.0000,  1.0000, -1.0000,  0.0000, -1.0000,  0.0000,  1.0000,  0.0000,  1.0000,
                                           1.0000, -1.0000, -1.0000,  1.0000, -1.0000,  0.0000,  0.3427,  0.0000, -1.0000, 0.0000,  1.0000, 0.0000, -1.0000,
                                           1.0000,  1.0000,  0.0000,  0.0000, -1.0000,  1.0000, -1.0000,  1.0000,  0.0000,  0.0000,  1.0000, -1.0000,  0.0000,
                                           0.0000)
                                     EE<-data.frame(WP,w1,w2,s1)   
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
         } else if(des=='EE42R7WP') {WP<-c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7)
                                     w1<-c(1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,
                                           1,  0,  0,  0,  0,  0,  0)
                                     w2<-c(-1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000, -1.0000,
                                           -1.0000, -1.0000, -1.0000, -1.0000, -1.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -1.0000, -1.0000,
                                           -1.0000, -1.0000, -1.0000, -1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  0.0977,  0.0977,  0.0977,
                                            0.0977,  0.0977,  0.0977)
                                     s1<-c(-1,  1, -1, -1,  1,  1,  0, -1, -1,  1,  1,  0,  0, -1,  1, -1,  1,  0, -1,  1,  1, -1,  0,  0,  1, -1,  0, -1,  0,  1,  1,  1, -1, -1, -1,
                                            1,  1,  0,  1, -1,  0, -1)
                                     EE<-data.frame(WP,w1,w2,s1) 
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
                                     w1<-c(1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1,  0,  0,  0,
                                           0,  1,  1,  1,  1,  0,  0,  0,  0)
                                     w2<-c(1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1, -1,
                                           -1,  1,  1,  1,  1,  1,  1,  1,  1)
                                     s1<-c(1.0000,  1.0000,  0.0000, -1.0000,  0.0000, -1.0000,  1.0000, -1.0000,  0.0000, -1.0000, -1.0000,  1.0000, -1.0000,
                                          -1.0000,  1.0000,  0.0000,  1.0000,  1.0000, -1.0000,  0.0000, -1.0000,  1.0000,  0.0000,  0.0000,  0.0000, -0.6924,
                                          -0.7143,  1.0000, -1.0000,  0.0000, -1.0000,  1.0000, -1.0000,  1.0000,  0.0000,  0.0000,  1.0000,  1.0000,  0.0000,
                                          -1.0000, -1.0000,  0.0000,  0.0000,  1.0000)
                                     EE<-data.frame(WP,w1,w2,s1)
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
         } else if(des=='EE45R9WP') {v <- c(235, 154,  73, 235,  73, 182,  20, 182,  20, 101, 210, 210,  48,  48, 129,   4, 166,  85, 166,   4, 221,  59, 140, 221,  59, 195,
                                            114,  33,  33, 195,  70, 232,  70, 151, 232, 107,  26,  26, 188, 188,  18,  99, 180,  18, 180)
                                     Full <-expand.grid(WP=c(1:9),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
        } else if(des=='EE48R8WP') {v <- c(145,   1,   1, 145,  73,   1,  50, 122, 194,  50, 194,  50, 171,  27, 171,  27,  27,  99,  20, 164,  20, 164,  92, 164, 213,  69,
                                           69, 213, 213, 141, 158,  14,  14, 158,  86,  86,  63, 207, 207,  63, 135, 135,  48, 192, 120, 192, 192,  48)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
        } else if(des=='EE48R12WP') {v <- c(241, 133, 133,  25,  98, 314, 206, 206, 315,  99, 207, 207,   4, 220,   4, 220,  65, 173, 281, 173, 294,  78,  78, 294, 223,   7,
                                            223,   7, 128,  20, 128, 236, 177,  69, 285, 177, 250, 142,  34, 142, 275,  59, 167, 167, 204, 312, 204,  96)
                                     Full <-expand.grid(WP=c(1:12),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:48)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:12),each=4))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:48)
                                                           #this randomizes subplots
                                                           s<-4
                                                           w<-12
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE)
        } else if(des=='EE50R10WP') {v <- c(61, 241, 151, 241,  61, 172, 262,  82,  82, 262,  33,  33, 213, 123, 213,  94,   4,   4, 184, 184,  55, 235, 235, 145,  55,  66,
                                            66, 246, 156, 246,  27,  27, 117, 207, 207, 168,  78, 258,  78, 258, 189, 189,  99,   9,   9,  20,  20, 200, 110, 200)
                                     Full <-expand.grid(WP=c(1:10),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
         } else if(des=='EE54R9WP') {v <- c(1, 163,   1,   1, 163, 163, 236, 236, 236,  74,  74,  74,  48, 210,  48, 210, 129, 129,  22, 184,  22, 184, 103, 103, 221, 221,
                                            59, 140,  59, 140, 195, 195, 114, 114,  33,  33,  25, 106, 106,  25, 187, 187, 179,  17,  98,  17, 179,  98, 234,  72, 234,  72,
                                            153, 153)
                                     Full <-expand.grid(WP=c(1:9),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
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
         } else if(des=='EE55R11WP') {v <- c(254, 254,  56,  56, 155, 266,  68,  68, 266, 167,  14, 113, 212,  14, 212,  92,  92, 290, 290, 191,  27, 126, 225,  27, 225, 226,
                                             226,  28,  28, 127, 282, 282, 183,  84,  84,  96,  96, 195, 294, 294, 108,   9, 207, 207,   9,  76, 274, 175, 274,  76, 242, 242,
                                             143,  44,  44)
                                     Full <-expand.grid(WP=c(1:11),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:55)
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
         } else if(des=='EE60R10WP') {WP<-c(1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  6,
                                            6,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10)
                                     w1<-c(0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1)
                                     w2<-c(1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,
                                           0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0)
                                     s1<-c(-0.98490,  0.92013, -0.99021, -0.99009,  0.92155, -0.98891,  0.91538,  0.99317, -0.81307,  1.00000, -1.00000,
                                           -1.00000,  0.89666,  0.90155,  0.89731, -0.97395, -0.97473,  0.89803,  0.95949, -0.98310,  1.00000, -1.00000,
                                           -0.22019,  0.07504, -0.82216, -0.99146,  0.90910,  1.00000,  1.00000, -1.00000, -0.06885, -0.08836,  0.97621,
                                           -0.98777, -1.00000,  1.00000,  1.00000,  0.00000, -1.00000, -1.00000,  1.00000, -1.00000, -1.00000,  0.00000,
                                            1.00000,  1.00000, -1.00000, -1.00000, -0.87078, -0.88205, -0.87636,  0.87867,  0.87426, -0.87620, -0.99937,
                                            0.95959, -1.00000, -0.14622,  0.01723,  1.00000)
                                     EE<-data.frame(WP,w1,w2,s1)
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
         } else if(des=='EE60R12WP') {v <- c(61, 277, 169,  61, 277,  74,  74, 290, 182, 290, 135,  27, 243, 243,  27,   4, 220, 220,   4, 112, 137, 245, 245,  29,  29, 294,
                                             78, 294,  78, 186, 103, 319, 319, 103, 211, 260,  44,  44, 152, 260,  21,  21, 237, 237, 129, 106, 214, 106, 322, 322,  95,  95,
                                             311, 311, 203,  12,  12, 228, 120, 228)
                                     Full <-expand.grid(WP=c(1:12),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1))
                                     EE <- Full[v, ]
                                     rownames(EE)<- c(1:60)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:12),each=5))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:60)
                                                           #this randomizes subplots
                                                           s<-5
                                                           w<-12
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 
         } else if(des=='EE66R11WP') {WP<-c(1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  6,
                                            6,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11)
                                     w1<-c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,
                                           1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
                                     w2<-c(0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,
                                            1,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1)
                                     s1<-c(-0.99966, -0.00031,  0.00000,  0.99999,  1.00000, -1.00000,  0.00131,  1.00000,  1.00000, -0.00053, -1.00000,
                                           -1.00000,  0.00079,  1.00000,  0.00000, -1.00000,  1.00000, -1.00000, -1.00000,  0.00092,  1.00000,  1.00000,
                                           -0.00014, -1.00000, -0.00001, -0.99970, -0.00030, -0.99996,  0.99999,  1.00000, -0.20031, -0.05055, -0.23939,
                                           -0.48349, -0.38893, -0.27241, -0.99965, -0.00033, -1.00000,  1.00000,  1.00000,  0.00000, -0.31444, -0.89282,
                                           -0.97489,  0.20502,  1.00000,  1.00000,  0.00002,  0.99972, -0.99993, -1.00000,  0.00021,  1.00000, -1.00000,
                                            0.99658, -1.00000, -1.00000, -1.00000,  0.99416,  0.79655,  0.71364, -0.90724, -0.97798, -0.96451, -0.95494)
                                     EE<-data.frame(WP,w1,w2,s1)
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
         } else if(des=='EE72R12WP') {WP<-c(1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  5,
                                            6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10,
                                            11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12)
                                     w1<-c(0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1,
                                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,
                                            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
                                     w2<-c(1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1,
                                           1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,
                                           -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
                                     s1<-c(0.13221,  0.98465, -0.98591, -1.00000, -0.10153,  1.00000, -1.00000, -0.99777, -0.99716, -1.00000,
                                           1.00000, -1.00000, -0.96235,  0.14533, -0.99940,  0.99980,  1.00000, -0.15396,  0.98485,  0.03364,
                                           0.99947, -0.99961, -1.00000,  0.01107,  0.97580, -0.96581,  0.99985,  0.00498, -1.00000, -1.00000,
                                           0.37902, -0.23213, -0.22691,  0.36009, -0.23354, -0.22438, -0.20814, -0.23277,  0.33819,  0.33275,
                                          -0.21581, -0.19911,  1.00000,  0.99858,  0.99813,  1.00000,  0.99933,  1.00000, -0.93947, -1.00000,
                                          -1.00000, -0.04571,  1.00000,  1.00000,  0.02544,  0.00086, -0.99058,  0.99370,  1.00000, -1.00000,
                                          -0.96182,  0.97958,  1.00000, -1.00000, -0.00294, -1.00000, -0.96054, -0.99975, -0.02311,  0.98957,
                                           1.00000, -0.99135)
                                     EE<-data.frame(WP,w1,w2,s1)
                                     if (randomize==TRUE) {#This orders the whole plots randomly
                                                           o<-c(rep(sample(1:12),each=6))
                                                           EE<-EE[order(order(o)), ]
                                                           rownames(EE)<-c(1:72)
                                                           #this randomizes subplots
                                                           s<-6
                                                           w<-12
                                                           rowp<-c(sample(1:s))
                                                           for (i in 1:(w-1)) {rowp<-c(rowp,i*s+sample(1:s))}
                                                           EE<-EE[order(rowp), ]
                                                          }
                                     return(EE) 

          } else cat(" Design name misspelled-Enter EEw2s1( ) to display list of names","\n")
                 }


