EEw1s3 <- function(des='',randomize=FALSE) {
  if (des=='')  {
cat(" ", "\n")
cat("Catalog of D-efficient Estimation Equivalent RS","\n")
cat("  Designs for (1 wp factor and  3 sp factors)  ","\n")
cat(" ", "\n")
cat("   Jones and Goos, JQT(2012) pp. 363-374","\n")
cat(" ", "\n")
cat(format("Design Name",width=11),format("whole plots",width=11),format("sub-plots/whole plot",width=21),"\n")
cat("----------------------------------------","\n")
cat(format("EE15R5WP", width=11),format("   5",width=11),format("           3",width=21),"\n")
cat(format("EE16R4WP", width=11),format("   4",width=11),format("           4",width=21),"\n")
cat(format("EE18R6WP", width=11),format("   6",width=11),format("           3",width=21),"\n")
cat(format("EE20R4WP", width=11),format("   4",width=11),format("           5",width=21),"\n")
cat(format("EE20R5WP", width=11),format("   5",width=11),format("           4",width=21),"\n")
cat(format("EE24R4WP", width=11),format("   4",width=11),format("           6",width=21),"\n")
cat(format("EE24R6WP", width=11),format("   6",width=11),format("           4",width=21),"\n")
cat(format("EE30R6WP", width=11),format("   6",width=11),format("           5",width=21),"\n")
cat(format("EE36R6WP", width=11),format("   6",width=11),format("           6",width=21),"\n")
cat(" ","\n")
cat("==> to retrieve a design type EEw1s3('EE15R5WP') etc.","\n")
         } else if (des=='EE15R5WP') {v <- c(241, 391, 121, 132,  42, 147, 218,  98, 293, 374,  74, 314, 320, 305,   5)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE16R4WP') { v <- c(17, 161, 209, 257,  74, 314,  26, 218, 243,  99, 291,   3, 312, 252, 120,  72)
                                     Full <-expand.grid(WP=c(1:4),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE18R6WP') { v <- c(316, 106, 400, 205, 142,  37, 521, 185, 437, 445, 151, 235, 551,   5, 551,  55, 391, 286)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE20R4WP') {WP<-c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4)
                                     w1<-c(-1.0000, -1.0000, -1.0000, -1.0000, -1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,
                                            1.0000,  1.0000,  1.0000,  1.0000, -0.2478, -0.2478, -0.2478, -0.2478, -0.2478)
                                     s1<-c(-1, -1,  1,  1, -1,  1, -1,  0,  1, -1, -1,  0, -1,  1,  1, -1,  1, -1,  1,  0)
                                     s2<-c(0,  1,  1, -1, -1,  0,  1, -1,  1, -1,  0,  1, -1,  1, -1, -1,  1,  1, -1,  0)
                                     s3<-c(-1,  1, -1,  1,  0, -1,  0, -1,  1,  1,  1,  1, -1, -1,  0,  1,  1, -1, -1,  0)
                                     EE<-data.frame(WP,w1,s1,s2,s3)
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
         } else if(des=='EE20R5WP') {v <- c(391,  91,  31, 271,   2, 362, 122, 302,  88, 238, 298,  13, 399, 204,  39, 279,  15, 330, 120, 180)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE24R4WP') {v <- c(221,  17,  77, 161, 317,  29, 290,   2, 254, 242, 206,  86, 315,  63, 231, 291,  75, 111, 228, 108, 288, 144, 300,  12)
                                     Full <-expand.grid(WP=c(1:4),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE24R6WP') {v <- c(1, 145, 361, 433,  80, 404, 242, 242, 471,  39, 327, 111, 412, 412, 250,  88, 485,  53, 179, 125, 432, 396, 72, 270)
                                     Full <-expand.grid(WP=c(1:6),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE25R5WP') {v <- c(136,  91, 376,  31, 346, 402, 282, 207,  42, 102, 233, 398, 278,  83, 113,  49, 364, 259,  19, 304, 210, 135,
                                            375, 315,  15)
                                     Full <-expand.grid(WP=c(1:5),w1=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE30R6WP') {WP<-c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6 )
                                     w1<-c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0, -1, -1, -1,
                                           -1, -1,  1,  1,  1,  1, 1)
                                     s1<-c(0.1467,  0.4372, -0.9533,  0.9523, -0.8610,  0.7324,  0.6926, -0.8613,  0.1524, -0.9942, -1.0000,
                                          -0.1228,  0.5316, -0.6869,  1.0000,  0.8588, -0.9366,  1.0000,  0.8337,  0.9420, -0.8481,  0.6643,
                                          -0.8878,  0.8998, -0.1064,  1.0000,  0.1705, -0.7848,  0.9997, -0.8476)
                                     s2<-c(0.9970,  0.7162,  0.0875, -0.8638, -0.9352,  0.8870, -0.7992,  0.8500,  0.0558, -0.9919,  0.9949,
                                          -1.0000,  0.5110, -0.8674,  0.3632, -0.6348,  1.0000, -0.8105,  0.9934, -0.8724,  0.8085, -0.8940,
                                          -0.8715,  0.9606, -0.0019,  0.8828,  1.0000,  0.9896, -0.9783,  0.8635)
                                     s3<-c(-0.9659,  0.7165,  1.0000,  0.2234, -0.8596, -0.8207,  0.9175,  0.9220,  0.0277, -0.9321, -0.0126,
                                           -1.0000, -0.7197,  0.8467,  1.0000, -1.0000, -1.0000,  0.9976,  0.5478, -0.9948, -0.8287, -0.8866,
                                            0.9148,  0.9609, -0.0459, -0.9117,  0.9884, -0.8769,  0.9997, -0.9999)
                                     EE<-data.frame(WP,w1,s1,s2,s3)
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
         } else if(des=='EE36R6WP') {WP<-c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6 )
                                     w1<-c(0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1,  0, 0,  0,  0,  0,
                                           0, -1, -1, -1, -1, -1, -1)
                                     s1<-c(0.53908,  0.76183,  0.85120, -0.65988, -0.94260,  0.96642,  0.47852,  0.53728,  0.32607,
                                          -0.81904, -1.00000, -0.92735,  0.50225, -0.19625, -0.60576,  1.00000,  1.00000,  0.08923,
                                          -0.58067, -0.16993, -0.94720, -0.96033,  0.92974,  0.32387,  0.92279, -1.00000,  0.87572,
                                          -0.45357,  0.17111,  1.00000, -0.40728, -0.82531, -1.00000, -0.57218,  0.70431,  0.69594)
                                     s2<-c(-0.79936,  0.89476,  0.42489,  0.63926,  0.88700, -0.37711,  0.88231, -0.76266,  0.16716,
                                           -0.74521, -0.50252,  0.90081, -0.11153,  0.46263,  0.85752, -0.02876,  0.00745,  0.35810,
                                           -0.93614,  1.00000, -0.10915, -0.20044, -0.63579,  0.82141,  0.39575,  0.65611, -1.00000,
                                            1.00000,  0.00885,  0.60873,  0.90481,  0.28202, -0.54116, -0.71167, -0.80632,  0.81221)
                                     s3<-c(-0.92414,  0.88755, -0.57336,  0.66455,  0.25221, -0.03008, -0.24680,  0.10316,  0.80028,
                                           -0.72825,  1.00000, -0.66508, -0.75278,  0.70306, -0.60648,  0.89481,  0.90464,  1.00000,
                                           -0.39641,  0.12134, -0.61963,  0.98790,  0.83953, -0.66942,  0.20694,  0.35728, -0.50765,
                                            1.00000, -1.00000,  0.22016, -0.02109,  0.76701, -0.79858, -0.00717,  1.00000, -0.67686)
                                     EE<-data.frame(WP,w1,s1,s2,s3)
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
          } else cat(" Design name misspelled-Enter EEw1s3( ) to display list of names","\n")
                 }






