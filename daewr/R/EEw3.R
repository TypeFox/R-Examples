EEw3 <- function(des='',randomize=FALSE) {
  if (des=='')  {
cat(" ", "\n")
cat("Catalog of D-efficient Estimation Equivalent RS","\n")
cat("      Designs for ( 3 wp factors )  ","\n")
cat(" ", "\n")
cat("   Jones and Goos, JQT(2012) pp. 363-374","\n")
cat(" ", "\n")
cat(format("Design Name",width=11),format("whole plots",width=11),format("sub-plots/whole plot",width=21),format("sub-plot factors",width=15),"\n")
cat("--------------------------------------------------------------","\n")
cat(format("EE22R11WP", width=11),format("   11",width=11),format("          2",width=21),format("           2",width=15),"\n")
cat(format("EE48R12WP", width=11),format("   12",width=11),format("          4",width=21),format("           3",width=15),"\n")
cat(" ","\n")
cat("==> to retrieve a design type EEw3('EE22R11WP') etc.","\n")
         } else if (des=='EE22R11WP') {v <- c(2663,  287, 2444,   68, 2005,  817, 1005, 2193, 1787,  599, 1876,  688,   29, 2405, 2054,  866,  207, 2583, 1539,  351,
                                              2541, 1353)
                                     Full <-expand.grid(WP=c(1:11),w1=c(-1,0,1),w2=c(-1,0,1),w3=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1))
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
     
         } else if(des=='EE48R12WP') {v <- c(241, 7045, 2833, 3805, 2594, 7778, 6482,    2, 5859,  675, 1971, 8451,  724, 8500, 5908, 2020, 6557, 2669,   77, 7853,
                                             8718,  942, 6126, 2238,  319, 8095, 2911, 6799, 8000, 7676, 2492, 3140, 8697,  921, 4485, 2217, 6382, 5734, 8002, 1198,
                                             1475, 8603, 6011, 5039, 7008, 4740,  204, 8304)
                                     Full <-expand.grid(WP=c(1:12),w1=c(-1,0,1),w2=c(-1,0,1),w3=c(-1,0,1),s1=c(-1,0,1),s2=c(-1,0,1),s3=c(-1,0,1))
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

          } else cat(" Design name misspelled-Enter EEw3( ) to display list of names","\n")
                 }





