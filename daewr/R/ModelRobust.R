ModelRobust <- function(des='',randomize=FALSE) {
  if (des=='')  {
cat(" ", "\n")
cat("       Model Robust Factorial Designs","\n")
cat(" ", "\n")
cat(" Li and Nachtsheim, JQT(2000) p. 345-352","\n")
cat(" ", "\n")
cat(format(" ",width=11),format(" ",width=7),format(" ",width=8),format("    g",width=13),"\n")
cat(format(" ",width=11),format(" ",width=7),format(" ",width=8),format("  Maximum",width=13),"\n")
cat(format(" ",width=11),format(" ",width=7),format("  m",width=8),format("  # of",width=13),"\n")
cat(format("Design Name",width=11),format("  runs",width=7),format("factors",width=8),format("interactions",width=13),"\n")
cat("----------------------------------------","\n")
cat(format("MR8m4g3", width=11),format("   8",width=7),format(" 4",width=8),format("    3",width=13),"\n")
cat(format("MR8m5g2", width=11),format("   8",width=7),format(" 5",width=8),format("    2",width=13),"\n")
cat(format("MR8m6g1", width=11),format("   8",width=7),format(" 6",width=8),format("    1",width=13),"\n")
cat(format("MR12m5g5", width=11),format("   12",width=7),format(" 5",width=8),format("    5",width=13),"\n")
cat(format("MR12m6g5", width=11),format("   12",width=7),format(" 6",width=8),format("    5",width=13),"\n")
cat(format("MR12m7g4", width=11),format("   12",width=7),format(" 7",width=8),format("    4",width=13),"\n")
cat(format("MR12m8g3", width=11),format("   12",width=7),format(" 8",width=8),format("    3",width=13),"\n")
cat(format("MR12m9g2", width=11),format("   12",width=7),format(" 9",width=8),format("    2",width=13),"\n")
cat(format("MR16m7g5", width=11),format("   16",width=7),format(" 7",width=8),format("    5",width=13),"\n")
cat(format("MR16m8g5", width=11),format("   16",width=7),format(" 8",width=8),format("    5",width=13),"\n")
cat(format("MR16m9g5", width=11),format("   16",width=7),format(" 9",width=8),format("    5",width=13),"\n")
cat(format("MR16m10g3", width=11),format("   16",width=7),format("10",width=8),format("    3",width=13),"\n")
cat(" ","\n")
cat("==> to retrieve a design type ModelRobust('MR8m4g3') etc.","\n")
         } else if (des=='MR8m4g3') {v <- c(1,  3,  6,  7, 10, 11, 14, 16)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(4, 5, 8, 1, 3, 2, 6, 7)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:8)
                                     if (randomize==TRUE) {MR <- MR[sample(1:8), ]}
                                     return(MR)
         } else if(des=='MR8m5g2') { v <- c(1,  6, 12, 15, 19, 22, 25, 32)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(4, 1, 5, 8, 3, 7, 2, 6)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:8)
                                     if (randomize==TRUE) {MR <- MR[sample(1:8), ]}
                                     return(MR)
         } else if(des=='MR8m6g1') { v <- c(8, 11, 21, 30, 34, 47, 52, 57)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(5, 4, 7, 3, 1, 2, 8, 6)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:8)
                                     if (randomize==TRUE) {MR <- MR[sample(1:8), ]}
                                     return(MR)                                   
         } else if(des=='MR12m5g5') {v <- c(3,  6, 10, 11, 13, 15, 18, 19, 21, 24, 28, 30)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(10, 6, 4, 9, 1, 5, 7, 8, 2, 3, 12, 11)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:12)
                                     if (randomize==TRUE) {MR <- MR[sample(1:12), ]}
                                     return(MR)                                   
         } else if(des=='MR12m6g5') {v <- c(3, 10, 13, 16, 18, 31, 34, 39, 52, 53, 59, 62)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(1, 12, 3, 4, 11, 6, 9, 7, 2, 10, 8, 5)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:12)
                                     if (randomize==TRUE) {MR <- MR[sample(1:12), ]}
                                     return(MR)           
         } else if(des=='MR12m7g4') {v <- c(3,  14,  18,  38,  44,  63,  74,  84,  85, 111, 119, 121)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1),G=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(4, 10, 8, 7, 11, 3, 1, 2, 9, 5, 12, 6)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:12)
                                     if (randomize==TRUE) {MR <- MR[sample(1:12), ]}
                                     return(MR)     
         } else if(des=='MR12m8g3') {v <- c(32,  55,  77,  81, 106, 127, 129, 140, 150, 180, 232, 233)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1),G=c(-1,1),H=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(10, 5, 1, 2, 3, 9, 6, 12, 7, 8, 4, 11)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:12)
                                     if (randomize==TRUE) {MR <- MR[sample(1:12), ]}
                                     return(MR)  
         } else if(des=='MR12m9g2') {v <- c(31,  65, 104, 173, 204, 242, 262, 308, 345, 425, 448, 471)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1),G=c(-1,1),H=c(-1,1), J=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(11, 3, 8, 4, 1, 6, 2, 10, 9, 5, 7, 12)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:12)
                                     if (randomize==TRUE) {MR <- MR[sample(1:12), ]}
                                     return(MR)  

         } else if(des=='MR16m7g5') {v <- c(4,   8,  15,  25,  43,  45,  54,  55,  65,  75,  85,  90, 101, 116, 121, 128)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1),G=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(8, 5, 4, 16, 11, 12, 1, 3, 14, 9, 7, 10, 13, 6, 2, 15)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:16)
                                     if (randomize==TRUE) {MR <- MR[sample(1:16), ]}
                                     return(MR)  

         } else if(des=='MR16m8g5') {v <- c(6,  17,  32,  46,  77,  88, 116, 123, 141, 162, 171, 183, 195, 218, 240, 241)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1),G=c(-1,1),H=c(-1,1),J=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(8, 3, 15, 1, 16, 14, 11, 2, 9, 4, 13, 7, 6, 5, 10, 12)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:16)
                                     if (randomize==TRUE) {MR <- MR[sample(1:16), ]}
                                     return(MR)  

         } else if(des=='MR16m9g5') {v <- c(34,  47,  62,  76,  93, 116, 151, 230, 269, 323, 394, 425, 436, 465, 472, 511)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1),G=c(-1,1),H=c(-1,1),J=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(14, 13,  1,  6,  5,  9, 11,  3,  2, 12, 10,  7, 16,  4, 15,  8)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:16)
                                     if (randomize==TRUE) {MR <- MR[sample(1:16), ]}
                                     return(MR)  

         } else if(des=='MR16m10g3') {v <- c(5,   28,  223,  249,  291,  334,  344,  491,  539,  613,  656,  756,  822,  898,  942, 1009)
                                     Full <-expand.grid(A=c(-1,1),B=c(-1,1),C=c(-1,1),D=c(-1,1),E=c(-1,1),F=c(-1,1),G=c(-1,1),H=c(-1,1),J=c(-1,1),K=c(-1,1))
                                     MR <- Full[v, ]
                                     v2 <- c(7,  8, 11, 16, 10, 12, 13,  3,  5, 14,  1,  4, 15,  2,  6,  9)
                                     MR <- MR[v2, ]
                                     rownames(MR)<- c(1:16)
                                     if (randomize==TRUE) {MR <- MR[sample(1:16), ]}
                                     return(MR)  

         } else cat(" Design name misspelled-Enter ModelRobust() to display list of names","\n")
                 }


