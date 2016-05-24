EEw2s3 <- function(des='',randomize=FALSE) {
  if (des=='')  {
cat(" ", "\n")
cat("Catalog of D-efficient Estimation Equivalent RS","\n")
cat("  Designs for (2 wp factors and  3 sp factors)  ","\n")
cat(" ", "\n")
cat("   Jones and Goos, JQT(2012) pp. 363-374","\n")
cat(" ", "\n")
cat(format("Design Name",width=11),format("whole plots",width=11),format("sub-plots/whole plot",width=21),"\n")
cat("----------------------------------------","\n")
cat(format("EE21R7WP", width=11),format("   7",width=11),format("           3",width=21),"\n")
cat(format("EE24R8WP", width=11),format("   8",width=11),format("           3",width=21),"\n")
cat(format("EE28R7WP", width=11),format("   7",width=11),format("           4",width=21),"\n")
cat(format("EE32R8WP", width=11),format("   8",width=11),format("           4",width=21),"\n")
cat(format("EE35R7WP", width=11),format("   7",width=11),format("           5",width=21),"\n")
cat(format("EE40R8WP", width=11),format("   8",width=11),format("           5",width=21),"\n")
cat(format("EE42R7WP", width=11),format("   7",width=11),format("           6",width=21),"\n")
cat(format("EE48R8WP", width=11),format("   8",width=11),format("           6",width=21),"\n")
cat(" ","\n")
cat("==> to retrieve a design type EEw2s3('EE21R7WP') etc.","\n")
         } else if (des=='EE21R7WP') {v <- c(1191,  183,  435,  492, 1311,  807,  969,  150, 1158, 1278, 1530,  522,   47, 1118, 1559, 1427, 1679, 41,  637,  196, 1645)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE24R8WP') { v <- c(161, 1889,  449, 1298,  362, 1226,   11, 1451, 1739,  284,  860, 1868,  981, 1341,  621, 1902, 462,  174, 1351, 1351,  631,  144, 1728, 1152)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE28R7WP') { v <- c(1569, 1317,  561,   57,  170, 1682,  422, 1178,   66, 1452,  507,  948, 1579,  193,  697,  508, 1153, 1657,  145,  397, 1287, 1539,  846,   27, 1694,  497,  623,  371)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE32R8WP') {v <- c(1129, 1417,1921, 409, 1442, 2,578, 1730, 1315, 163, 1891, 451,  844, 1564, 1924,  556, 1701, 1269,  117, 1773,  438, 1878,  150, 1302, 1511,  647, 1871,   71, 1360,  208, 1000,  496)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE35R7WP') {v <- c(1639, 1135,  127,  379,  820, 1129, 1381,  184, 1255,  436,  990, 1557,   45,  549, 1305,  508, 1516,
                                            1264,    4,  823, 1489, 1615,  355,  607,  418, 1532, 1280,   83,  209,  524,  245, 1694,  749, 1190,
                                            497 )
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE40R8WP') {v <- c(1433, 1721, 1793,   65,  641, 1490,  122, 1562,  626, 1130, 1875,  147,  651, 1587,  435, 1748, 1460,
                                            884,  596,   20, 1357, 1213, 1933,  277,  205, 1894,  166, 1318,  454,  886,  295, 1231, 1447,    7,
                                           1735, 1832, 1328,  392,  464,  824)
                                     Full <-expand.grid(WP=c(1:8),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE42R7WP') {v <- c(1310, 1562,   50,  932, 1373,  491, 1682,  611, 1241,  548,  422,  170, 1641, 1011,  192,  129, 1137,
                                            1263, 1152, 1656,  396,  333,  711,   81, 1552,   40, 1300, 1426,  796,  544, 1518,  510,    6, 1266,
                                            636, 1455, 1197, 1323, 1701, 1008,  315,  189)
                                     Full <-expand.grid(WP=c(1:7),w1=c(-1,0,1),w2=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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
         } else if(des=='EE48R8WP') {WP<-c(1,1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8,8, 8, 8, 8, 8)
                                     w1<-c(-1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1)
                                     w2<-c(-1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1)
                                     s1<-c(0.56084, -0.06357,  1.00000,  0.51562, -0.06737,  0.70429, -0.11015, -0.95967, -0.20421,
                                          -0.35938,  0.99579, -0.99950,  1.00000, -0.86437, -0.85703, -0.09420, -0.92222,  0.88506,
                                           0.56986, -0.21668, -0.07369,  0.91509,  0.91019,  0.86681,  1.00000,  1.00000,  1.00000,
                                           0.09854,  0.15793,  0.17613, -0.12690, -0.90367, -0.66383,  0.10124,  0.92932, -0.97328,
                                          -0.78345, -1.00000,  0.15928,  0.75899, -0.97078,  0.98320, -0.48317,  0.68254, -0.29762,
                                          -0.43140,  1.00000,  0.57906)
                                     s2<-c(-0.92854,  0.10102, -1.00000, -1.00000,  0.18066, -0.90389, -0.92912, -0.15766, -0.95595,
                                            0.90964,  0.63698, -0.23680,  0.91050, -0.59467, -0.97021, -0.55536,  0.74053, -0.89746,
                                           -1.00000, -1.00000, -0.99894,  0.82650,  0.86108,  1.00000, -1.00000, -1.00000, -1.00000,
                                            0.99876, -1.00000, -1.00000,  0.96027, -0.84218,  0.62294, -0.45932, -0.09132, -0.92330,
                                           -0.78227,  0.84761, -0.99991,  0.73327, -0.92855, -0.23682, -0.42720,  0.98957,  0.54908,
                                           -0.42988,  0.87858,  0.99915)
                                     s3<-c(-0.92486,  0.97572,  0.07996, -0.56142,  0.97940, -0.93836,  0.96591, -0.24995, -0.91156,
                                           -0.89274,  0.99626,  0.99349, -0.99128, -0.72718, -0.97457, -0.32284,  0.90005,  0.75505,
                                           -0.92204,  0.58001,  0.52779, -0.88350, -0.84440,  0.95054, -0.00586,  0.06026,  0.02739,
                                            1.00000,  0.32016,  0.10662, -1.00000,  0.90604,  0.87733,  0.11895,  0.99909, -1.00000,
                                            0.51880, -1.00000, -0.89626,  1.00000, -0.05022, -0.93309,  1.00000, -0.75242,  0.97294,
                                            0.99664, -1.00000, -0.75579)
                                     EE<-data.frame(WP,w1,w2,s1,s2,s3)
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

          } else cat(" Design name misspelled-Enter EEw2s3( ) to display list of names","\n")
                 }






