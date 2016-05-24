# nfca main function
nfca <- function(data, type=0, method='hist', choice=1, n=30, alpha=0.05) {

        # Assign values to variables
        data.type <- type
        method.type <- method
        hist.manual <- choice
        n.CI <- n
        alpha.CI <- alpha
        
        # nFCAmain function: Generating hierarchical clustering using nFCA algorithm
        nFCAmain <- function(data, manual=1, method='hist', cutoff1=3, cutoff2=10, n=30, alpha=0.05) {
                    
                    # datafile <- data
                    cin <- n
                    exec <- 1
                    
                    # Reading the numerical matrix
                    xx <- data
                    names <- colnames( xx )
                    n <- dim(xx)[1]
                    x <- matrix(0, ncol=n, nrow=n)
                    for ( i in 1:n ) {
                        for ( j in 1:n ) {
                            if ( isTRUE( all.equal(i,j) ) ) {
                               x[i,j] <- 1
                            } else {
                               x[i,j] <- xx[i,j][[1]]
                            }
                        }
                    }
                    thresholdlist <- rep( 0, length(names) )
                    ithreshold <- 1
                    
                    # Form the data for Histogram
                    histdata <- rep( 0, (n^2-n)/2 )
                    k <- 1
                    for ( i in 2:n ) {
                        for ( j in 1: (i-1) ) {
                            histdata[ k ] <- x[j,i]
                            k <- k + 1
                        }
                    }
                    
                    # Display Histogram and determining the threshold
                    if ( isTRUE( all.equal(1,manual) ) ) {
                       hist( histdata )
                       ires <- 1
                       while ( !isTRUE( all.equal(ires,0) ) ) {
                             threshold <- readline( "Please input the threshold: ")
                             ires <- thresholdTEST( threshold, names, x )
                             if ( isTRUE( all.equal(1,ires) ) ) {
                                print( "This threshold is not large enough. Please choose another larger threshold!" )
                                flush.console()
                             }
                             if ( isTRUE( all.equal(2,ires) ) ) {
                                print( "The threshold is too large. Please choose a smaller one! " )
                                flush.console()
                             }
                       }
                       thresholdlist[ ithreshold ] <- threshold
                       ithreshold <- ithreshold + 1
                    } else {
                       if ( method=="hist" ) {
                          threshold <- hist( histdata,plot=F )$breaks[ (length(hist(histdata,plot=F)$breaks)-1) ]
                          ires <- thresholdTEST( threshold, names, x )
                          while ( !isTRUE( all.equal(ires,0) ) ) {
                                hgi <- length( histdata ) 
                                threshold <- sort(histdata)[ hgi ]
                                ires <- thresholdTEST( threshold, names, x )
                                while ( !isTRUE( all.equal(ires ,0) ) ) {
                                      if ( isTRUE( all.equal(1,hgi) ) ) {
                                         ires <- 0
                                         print( "Automatically choosing threshold failed!" )
                                         flush.console()
                                         exec <- 0
                                      } else {
                                         hgi <- hgi-1
                                         threshold <- sort(histdata)[ hgi ]
                                         ires <- thresholdTEST( threshold, names, x )
                                      }
                                }
                          }
                          thresholdlist[ ithreshold ] <- threshold
                          ithreshold <- ithreshold + 1                          
                       } else {
                          cihistdata <- sort(histdata)
                          cires <- 111
                          cii <- length( cihistdata )
                          while ( !isTRUE( all.equal(cires , 0) ) ) {
                                if ( isTRUE( all.equal(1,cires) ) ) {
                                   hgi <- length( histdata ) 
                                   threshold <- sort(histdata)[ hgi ]
                                   ires <- thresholdTEST( threshold, names, x )
                                   while ( !isTRUE( all.equal(ires ,0) ) ) {
                                         if ( isTRUE( all.equal(1,hgi) ) ) {
                                            ires <- 0
                                            print( "Automatically choosing threshold failed!" )
                                            flush.console()
                                            exec <- 0
                                         } else {
                                            hgi <- hgi-1
                                            threshold <- sort(histdata)[ hgi ]
                                            ires <- thresholdTEST( threshold, names, x )
                                         }
                                   }
                                   cires <- ires
                                } else {
                                   ciindex <- cii
                                   cir <- cihistdata[ciindex]
                                   ciz <- 0.5*log( (1+cir)/(1-cir) )
                                   cilz <- ciz-qnorm( 1-alpha/2 )*sqrt(1/(cin-3))
                                   cithreshold <- ( exp(2*cilz)-1 )/( exp(2*cilz)+1 )
                                   cires <- thresholdTEST( cithreshold, names, x )
                                   if ( isTRUE(all.equal(0,cires)) ) {
                                      threshold <- cithreshold
                                   } else {
                                      cii <- cii-1
                                   }
                                }
                          }
                          thresholdlist[ ithreshold ] <- threshold
                          ithreshold <- ithreshold + 1        
                       }
                    }
                    if ( isTRUE( all.equal(1,exec) ) ) {
                       
                       # Threshold the data matrix
                       xnames <- names
                       xbin <- x
                       xbin[ which( xbin<threshold ) ] <- 0
                       xbin[ which( xbin>=threshold ) ] <- 1
                       
                       # Apply FCA to construct the concept lattice
                       if ( file.exists( paste( "graph",".txt",sep="") )==TRUE ) {
                          file.remove( paste( "graph",".txt",sep="") )
                       }
                       file.copy( "graph0.txt", "graph.txt" )
                       if ( file.exists( paste( "graph001",".txt",sep="") )==TRUE ) {
                          file.remove( paste( "graph001",".txt",sep="") )
                       }
                       file.copy( "graph00.txt", "graph001.txt" )
                       if ( file.exists( "example.oal" )==TRUE ) {
                          file.remove( "example.oal" )
                       }
                       for ( i in 1:(dim(xbin)[1]) ) {
                           rowname <- paste( names[i], ":", sep="" )
    	                     colname <- ""
	                         for ( j in 1:(dim(xbin)[2]) ) {
		                           if ( isTRUE(all.equal(1,xbin[i,j])) ) {
			                            colname <- paste( colname, names[j], ";", sep="" )
                               }
                           }
	                         write( paste( rowname, colname, sep="" ), file=paste("example",".oal",sep=""), append=T )
                       }
                       system( "ruby multistagefca.rb" )
                       file.rename( "example.dot", "example.txt" )
                       
                       # Define and generating the cluster matrix
                       xcluster <- matrix( 0, ncol=(n+2), nrow=50 )
                       tlevel <- 1
                       cindex <- 1
                       concepts <- as.matrix( read.csv( "example.txt", header=F ) )
                       index1 <- grep( "--", concepts )
                       iorder <- rep( 0, length(index1) )
                       for ( i in 1:length(index1) ) {
                           xbefore <- constr( concepts[ index1[i] ] )
                           iloc <- which( xbefore=="-" )
                           if ( isTRUE(all.equal(1,iloc)) ) {
                              iorder[ i ] <- 0
                           } else {
                              if ( ( isTRUE(all.equal(2,iloc)) ) & (  isTRUE(all.equal(1,( length(xbefore) - iloc ))) ) ) {
                                 iorder[ i ] <- 0
                              } else {
                                 if ( ( iloc-1 ) > ( length( xbefore ) - iloc ) ) {
                                    iorder [i] <- 0
                                 } else {
                                    iorder[i] <- iloc-1
                                 }
                              }
                           }
                       }
                       
                       index2 <- which( iorder != 0 )
                       index1 <- index1[ index2 ]
                       iorder <- iorder[ index2 ]
                       index0 <- index1[ sort( iorder, index.return=T )$ix ]
                       conindex <- rep( 0, length(concepts) )
                       
                       for ( i in index0 ) {
                           if ( isTRUE( all.equal(0,conindex[i]) ) ) {
                              conindex[i] <- 1
                              xbefore <- constr( concepts[i] )
                              iloc <- which( xbefore=="-" )
                              if ( isTRUE( all.equal(2,iloc) ) ) {
                                 ifirst <- xbefore[iloc-1]
                                 xcluster[ cindex, which( xnames==ifirst ) ] <- 1
                                 xcluster[ cindex, (length(xnames)+1) ] <- tlevel
                                 if ( file.exists( paste( "graph",cindex,".txt",sep="") )==TRUE ) {
                                    file.remove( paste( "graph",cindex,".txt",sep="") )
                                 }
                                 write( paste("subgraph", paste( "cluster",cindex, sep=""), "{", sep=" "), length( paste("subgraph", paste("cluster",cindex,sep=""), "{", sep=" ") ), file=paste( "graph",cindex,".txt",sep=""), append=T )
                                 for ( j in index0 ) {
                                     if ( isTRUE( all.equal(0,conindex[j]) ) ) {
                                        newxbefore <- constr( concepts[j] )
                                        newloc <- which( newxbefore=="-" )
                                        if ( (length( which( newxbefore[1:(newloc-1)]==ifirst ) )>0 ) & ( isTRUE(all.equal(newloc,iloc)) ) ) {
                                           conindex[j] <- 1
                                        }
                                        if ( (length( which( newxbefore[1:(newloc-1)]==ifirst ) )>0 ) & ( newloc>iloc ) ) {
                                           conindex[ j ] <- 1
                                           for ( ii in 1: (newloc-1) ) {
                                               if ( !isTRUE( all.equal(newxbefore[ii],ifirst) ) ) {
                                                  xcluster[ cindex, which( xnames==newxbefore[ii] ) ] <- 1
                                                  write( paste( ifirst, "--", newxbefore[ii], paste("[label=", format( x[which(names==ifirst),which(names==newxbefore[ii])],digits=2 ),",","dir=forward];",sep=""), sep=" " ), length( paste( ifirst, "--", newxbefore[ii], paste("[label=", format( x[which(names==ifirst),which(names==newxbefore[ii])],digits=2),",","dir=forward];",sep=""), sep=" " ) ), file=paste( "graph",cindex,".txt",sep=""), append=T )
                                               }
                                           }
                                        }
                                     }
                                 }
                                 write( paste( "}", sep=" " ), length( paste( "}", sep=" " ) ), file=paste( "graph",cindex,".txt",sep="" ), append=T )
                                 xcluster[ cindex, (length(names)+2) ] <- threshold
                                 cindex <- cindex + 1
                              } else {
                                 if ( isTRUE( all.equal((iloc-1),(length(xbefore)-iloc)) ) ) {
                                    ifirst <- xbefore[1:(iloc-1)]
                                    for ( j in 1:(iloc-1) ) {
                                        xcluster[ cindex, which( xnames==ifirst[j] ) ] <- 1
                                    }
                                    xcluster[ cindex, (length(xnames)+1) ] <- tlevel
                                    if ( file.exists( paste( "graph",cindex,".txt",sep="") )==TRUE ) {
                                       file.remove( paste( "graph",cindex,".txt",sep="") )
                                    }
                                    write( paste("subgraph", paste( "cluster",cindex, sep=""), "{", sep=" "), length( paste("subgraph", paste( "cluster",cindex, sep=""), "{", sep=" ") ), file=paste( "graph",cindex,".txt",sep=""), append=T )
                                    for ( j in 1:(iloc-2) ) {
                                        write( paste( ifirst[j], "--", ifirst[j+1], paste("[label=", format( x[which(names==ifirst[j]),which(names==ifirst[j+1])],digits=2 ),"];",sep=""),sep=" " ), length(paste( ifirst[j], "--", ifirst[j+1],paste("[label=", format( x[which(names==ifirst[j]),which(names==ifirst[j+1])],digits=2 ),",","dir=forward];",sep=""), sep=" " )),file=paste( "graph",cindex,".txt",sep=""), append=T )
                                    }
                                    if ( iloc > 3 ) {
                                       write( paste( ifirst[iloc-1], "--", ifirst[1],paste("[label=", format( x[which(names==ifirst[iloc-1]),which(names==ifirst[1])],digits=2 ),"];",sep=""), sep=" " ), length( paste( ifirst[iloc-1], "--", ifirst[1],paste("[label=", format( x[which(names==ifirst[iloc-1]),which(names==ifirst[1])],digits=2 ),",","dir=forward];",sep=""), sep=" " )),file=paste( "graph",cindex,".txt",sep=""), append=T )
                                    }
                                    write( paste( "}", sep=" " ), length( paste( "}", sep=" " ) ), file=paste( "graph",cindex,".txt",sep="" ), append=T )
                                    xcluster[ cindex, (length(names)+2) ] <- threshold
                                    cindex <- cindex + 1
                                 }
                              }
                           }
                       }
                       
                       # Determine the column number of new matrix
                       istatus <- rep( 0, length(names) )
                       for ( i in 1:length(names) ) {
                           if ( !isTRUE( all.equal(length( which( xcluster[,i]==1 ) ),0) ) ) {
                              istatus[i] <- max( which( xcluster[,i]==1 )  )
                           } else {
                              istatus[i] <- 0
                           }
                       }
                       xnames <- ""
                       inames <- 1
                       for ( i in 1:length(names) ) {
                           if ( !isTRUE( all.equal(istatus[i],0) ) ) {
                              if ( isTRUE( all.equal( length( which(xnames==paste( istatus[i],sep="") ) ),0) ) ) {
                                 xnames[inames] <- paste( istatus[i],sep="" )
                                 inames <- inames+1
                              }
                           } else {
                              xnames[inames] <- names[ i ]
                              inames <- inames+1
                           }
                       }
                       icol <- length( xnames )           
                       xx <- matrix( 0, ncol=icol, nrow=icol )
                       xxavg <- matrix( 0, ncol=icol, nrow=icol )
                       for ( i in 1:(icol-1) ) {
                           for ( j in (i+1):icol ) {
                               if ( ( isTRUE( all.equal(length( which(names==xnames[i]) ),0) ) ) & ( isTRUE( all.equal(length( which(names==xnames[j]) ),0) ) ) ) {
                                  iindex <- which( xcluster[ as.numeric( xnames[i] ),1:n]==1 )
                                  jindex <- which( xcluster[ as.numeric( xnames[j] ),1:n]==1 )
                               }
                               if ( ( isTRUE( all.equal(length( which(names==xnames[i]) ),0) ) ) & ( !isTRUE( all.equal( length( which(names==xnames[j]) ),0) ) ) ) {
                                  iindex <- which( xcluster[ as.numeric( xnames[i] ),1:n]==1 )
                                  jindex <- which( names==xnames[j] )
                               }
                               if ( (!isTRUE( all.equal( length( which(names==xnames[i]) ),0) ) ) & (isTRUE( all.equal(length( which(names==xnames[j]) ),0) ) ) ) {
                                  iindex <- which( names==xnames[i] )
                                  jindex <- which( xcluster[ as.numeric( xnames[j] ),1:n]==1 )
                               }
                               if ( (!isTRUE( all.equal( length( which(names==xnames[i]) ),0) ) ) & ( !isTRUE( all.equal(length( which(names==xnames[j]) ),0) ) ) ) {
                                  iindex <- which( names==xnames[i] )
                                  jindex <- which( names==xnames[j] )
                               }
                               temp <- 0
                               tempavg <- 0
                               for ( ii in iindex ) {
                                   for ( jj in jindex ) {
                                       temp <- max( temp, x[ii,jj] )
                                       tempavg <- tempavg + x[ii,jj]
                                   }
                               }
                               xx[i,j] <- temp
                               xxavg[i,j] <- tempavg / ( length(iindex)*length(jindex) )
                           }
                       }
                       for ( i in 1: icol ) {
                           xx[i,i] <- 1
                           xxavg[i,i] <- 1
                       }
                       for ( i in 2: icol ) {
                           for ( j in 1:(i-1) ) {
                               xx[i,j] <- xx[ j,i]
                               xxavg[i,j] <- xxavg[j,i]
                           }
                       }
                       while ( icol>2 ) {
                             
                             # MAIN LOOP STARTS FROM HERE
                             histdata <- rep( 0, (icol^2-icol)/2 )
                             k <- 1
                             for ( i in 2:icol ) {
                                 for ( j in 1: (i-1) ) {
                                     histdata[ k ] <- xx[j,i]
                                     k <- k + 1
                                 }
                             }
                             
                             ###################################################
                             if ( isTRUE( all.equal(manual, 1) ) ) {
                                hist( histdata )
                                ires <- 1
                                while ( !isTRUE( all.equal(ires,0) ) ) {
                                      threshold <- readline( "Please input the threshold: ")
                                      ires <- thresholdTEST( threshold, xnames, xx )
                                      if ( isTRUE( all.equal(ires,1) ) ) {
                                         print( "This threshold is not large enough. Please choose another larger threshold!" )
                                         flush.console()
                                      }
                                      if ( isTRUE( all.equal(ires,2) ) ) {
                                         print( "The threshold is too large. Please choose a smaller one! " )
                                         flush.console()
                                      }
                                }
                                thresholdlist[ ithreshold ] <- threshold
                                ithreshold <- ithreshold + 1
                             } else {
                                if ( method=="hist" ) {
                                   threshold <- hist( histdata,plot=F )$breaks[ (length(hist(histdata,plot=F)$breaks)-1) ]
                                   ires <- thresholdTEST( threshold, xnames, xx )
                                   while ( !isTRUE( all.equal(ires,0) ) ) {
                                         hgi <- length( histdata ) 
                                         threshold <- sort(histdata)[ hgi ]
                                         ires <- thresholdTEST( threshold, xnames, xx )
                                         while ( !isTRUE( all.equal(ires ,0) ) ) {
                                               if ( isTRUE( all.equal(hgi,1) ) ) {
                                                  ires <- 0
                                                  print( "Automatically choosing threshold failed!" )
                                                  flush.console()
                                                  exec <- 0
                                               } else {
                                                  hgi <- hgi-1
                                                  threshold <- sort(histdata)[ hgi ]
                                                  ires <- thresholdTEST( threshold, xnames, xx )
                                               }
                                         }
                                   }
                                   thresholdlist[ ithreshold ] <- threshold
                                   ithreshold <- ithreshold + 1                          
                                } else {
                                   cihistdata <- sort(histdata)
                                   cires <- 111
                                   cii <- length( cihistdata )
                                   while ( !isTRUE( all.equal(cires , 0) ) ) {
                                         if ( isTRUE( all.equal(cires, 1) ) ) {
                                            hgi <- length( histdata ) 
                                            threshold <- sort(histdata)[ hgi ]
                                            ires <- thresholdTEST( threshold, xnames, xx )
                                            while ( !isTRUE( all.equal(ires ,0) ) ) {
                                                  if ( isTRUE( all.equal(hgi,1) ) ) {
                                                     ires <- 0
                                                     print( "Automatically choosing threshold failed!" )
                                                     flush.console()
                                                     exec <- 0
                                                  } else {
                                                     hgi <- hgi-1
                                                     threshold <- sort(histdata)[ hgi ]
                                                     ires <- thresholdTEST( threshold, xnames, xx )
                                                  }
                                            }
                                            cires <- ires
                                         } else {
                                            ciindex <- cii
                                            cir <- cihistdata[ciindex]
                                            ciz <- 0.5*log( (1+cir)/(1-cir) )
                                            cilz <- ciz-qnorm( 1-alpha/2 )*sqrt(1/(cin-3))
                                            cithreshold <- ( exp(2*cilz)-1 )/( exp(2*cilz)+1 )
                                            cires <- thresholdTEST( cithreshold, xnames, xx )
                                            if ( isTRUE( all.equal(cires, 0) ) ) {
                                               threshold <- cithreshold
                                            } else {
                                               cii <- cii-1
                                            }
                                         }
                                   }
                                   thresholdlist[ ithreshold ] <- threshold
                                   ithreshold <- ithreshold + 1        
                                }
                             }
                             
                             ###################################################
                             
                             # Threshold the data matrix
                             xbin <- xx
                             xbin[ which( xbin<threshold ) ] <- 0
                             xbin[ which( xbin>=threshold ) ] <- 1
                             tlevel <- tlevel + 1
                             
                             # Apply FCA to construct the concept lattice
                             if ( file.exists( "example.oal" )==TRUE ) {
                                file.remove( "example.oal" )
                             }
                             for ( i in 1:(dim(xbin)[1]) ) {
                                 rowname <- paste( xnames[i], ":", sep="" )
    	                           colname <- ""
	                               for ( j in 1:(dim(xbin)[2]) ) {
		                                 if ( isTRUE( all.equal(xbin[i,j], 1) ) ) {
			                                  colname <- paste( colname, xnames[j], ";", sep="" )
                                     }
                                 }
	                               write( paste( rowname, colname, sep="" ), file=paste("example",".oal",sep=""), append=T )
                             }
                             system( "ruby multistagefca.rb" )
                             file.rename( "example.dot", "example.txt" )
                             
                             # Take all the concepts and put them in increasing order
                             concepts <- as.matrix( read.csv( "example.txt", header=F ) )
                             index1 <- grep( "--", concepts )
                             iorder <- rep( 0, length(index1) )
                             for ( i in 1:length(index1) ) {
                                 xbefore <- constr( concepts[ index1[i] ] )
                                 iloc <- which( xbefore=="-" )
                                 if ( isTRUE( all.equal(iloc,1) ) ) {
                                    iorder[ i ] <- 0
                                 } else {
                                    if ( (isTRUE( all.equal(iloc, 2) )) & ( isTRUE( all.equal(( length( xbefore ) - iloc ),1) ) ) ) {
                                       iorder[ i ] <- 0
                                    } else {
                                       if ( ( iloc-1 ) > ( length( xbefore ) - iloc ) ) {
                                          iorder [i] <- 0
                                       } else {
                                          iorder[i] <- iloc-1
                                       }
                                    }
                                 }
                             }
                             
                             index2 <- which( iorder != 0 )
                             index1 <- index1[ index2 ]
                             iorder <- iorder[ index2 ]
                             index0 <- index1[ sort( iorder, index.return=T )$ix ]
                             conindex <- rep( 0, length(concepts) )
                             for ( i in index0 ) {
                                 if ( isTRUE( all.equal(conindex[i],0) ) ) {
                                    conindex[i] <- 1
                                    xbefore <- constr( concepts[i] )
                                    iloc <- which( xbefore=="-" )
                                    if ( isTRUE( all.equal(iloc,2) ) ) {
                                       ifirst <- xbefore[iloc-1]
                                       if ( isTRUE( all.equal(length( which( names==ifirst ) ) , 0) ) ) {
                                          xcluster[ cindex, which( xcluster[as.numeric(ifirst),]==1 ) ] <- 1
                                          xcluster[ cindex, (length(names)+1) ] <- tlevel
                                       } else {
                                          xcluster[ cindex, which( names==ifirst ) ] <- 1
                                          xcluster[ cindex, (length(names)+1) ] <- tlevel
                                       }
                                       if ( file.exists( paste( "graph",cindex,".txt",sep="") )==TRUE ) {
                                          file.remove( paste( "graph",cindex,".txt",sep="") )
                                       }
                                       write( paste("subgraph", paste("cluster",cindex,sep=""), "{", sep=" "), length( paste("subgraph", paste("cluster",cindex,sep=""), "{", sep=" ") ), file=paste( "graph",cindex,".txt",sep=""), append=T )
                                       for ( j in index0 ) {
                                           if ( isTRUE( all.equal(conindex[j],0) ) ) {
                                              newxbefore <- constr( concepts[j] )
                                              newloc <- which( newxbefore=="-" )
                                              if ( (length( which( newxbefore[1:(newloc-1)]==ifirst ) )>0 ) & ( isTRUE(all.equal(newloc,iloc)) ) ) {
                                                 conindex[j] <- 1
                                              }
                                              if ( (length( which( newxbefore[1:(newloc-1)]==ifirst ) )>0 ) & ( newloc>iloc ) ) {
                                                 conindex[ j ] <- 1
                                                 for ( ii in 1: (newloc-1) ) {
                                                     if ( !isTRUE(all.equal(newxbefore[ii] , ifirst)) ) {
                                                        if ( isTRUE(all.equal(length( which( names==newxbefore[ii] ) ) , 0)) ) {
                                                           xcluster[ cindex, which( xcluster[as.numeric(newxbefore[ii]),]==1 ) ] <- 1
                                                           xcluster[ cindex, (length(names)+1) ] <- tlevel
                                                        } else {
                                                           xcluster[ cindex, which( names==newxbefore[ii] ) ] <- 1
                                                           xcluster[ cindex, (length(names)+1) ] <- tlevel
                                                        }
                                                        write( paste( ifirst, "--", newxbefore[ii],paste("[label=", format( xxavg[which(xnames==ifirst),which(xnames==newxbefore[ii])],digits=2 ),",","dir=forward];",sep=""), sep=" " ), length( paste( ifirst, "--", newxbefore[ii],paste("[label=", format( xxavg[which(xnames==ifirst),which(xnames==newxbefore[ii])],digits=2 ),",","dir=forward];",sep=""), sep=" " ) ), file=paste( "graph",cindex,".txt",sep=""), append=T )
                                                     }
                                                 }
                                              }
                                           }
                                       }
                                       write( paste( "}", sep=" " ), length( paste( "}", sep=" " ) ), file=paste( "graph",cindex,".txt",sep="" ), append=T )
                                       xcluster[ cindex, (length(names)+2) ] <- threshold
                                       cindex <- cindex + 1
                                    }
                                    if ( isTRUE( all.equal((iloc-1),(length(xbefore)-iloc)) ) ) {
                                       ifirst <- xbefore[1:(iloc-1)]
                                       for ( j in 1:(iloc-1) ) {
                                           if ( isTRUE( all.equal(length( which( names==ifirst[j] ) ) , 0) ) ) {
                                              xcluster[ cindex, which( xcluster[as.numeric(ifirst[j]),1:n]==1 ) ] <- 1
                                              xcluster[ cindex, (length(names)+1) ] <- tlevel
                                           } else {
                                              xcluster[ cindex, which( names==ifirst[j] ) ] <- 1
                                              xcluster[ cindex, (length(names)+1) ] <- tlevel
                                           }
                                       }
                                       xcluster[ cindex, (length(names)+1) ] <- tlevel
                                       if ( file.exists( paste( "graph",cindex,".txt",sep="") )==TRUE ) {
                                          file.remove( paste( "graph",cindex,".txt",sep="") )
                                       }
                                       write( paste("subgraph", paste("cluster",cindex,sep=""), "{", sep=" "), length( paste("subgraph", paste("cluster",cindex,sep=""), "{", sep=" ") ), file=paste( "graph",cindex,".txt",sep=""), append=T )
                                       for ( j in 1:(iloc-2) ) {
                                           write( paste( ifirst[j], "--", ifirst[j+1], paste("[label=", format( xxavg[which(xnames==ifirst[j]),which(xnames==ifirst[j+1])],digits=2 ),"];",sep=""),sep=" " ), length(paste( ifirst[j], "--", ifirst[j+1], paste("[label=", format( xxavg[which(xnames==ifirst[j]),which(xnames==ifirst[j+1])],digits=2 ),",","dir=forward];",sep=""),sep=" " )),file=paste( "graph",cindex,".txt",sep=""), append=T )
                                       }
                                       if ( iloc > 3 ) {
                                          write( paste( ifirst[iloc-1], "--", ifirst[1], paste("[label=", format( xxavg[which(xnames==ifirst[iloc-1]),which(xnames==ifirst[1])],digits=2 ),"];",sep="" ), sep=" "), length(paste( ifirst[iloc-1], "--", ifirst[1], paste("[label=", format( xxavg[which(xnames==ifirst[iloc-1]),which(xnames==ifirst[1])],digits=2 ),",","dir=forward];",sep=""),sep=" " )),file=paste( "graph",cindex,".txt",sep=""), append=T )
                                       }
                                       write( paste( "}", sep=" " ), length( paste( "}", sep=" " ) ), file=paste( "graph",cindex,".txt",sep="" ), append=T )
                                       xcluster[ cindex, (length(names)+2) ] <- threshold
                                       cindex <- cindex + 1
                                    }
                                 }
                             }
                             istatus <- rep( 0, length(names) )
                             for ( i in 1:length(names) ) {
                                 if ( !isTRUE(all.equal(length( which( xcluster[,i]==1 ) ) ,0)) ) {
                                    istatus[i] <- max( which( xcluster[,i]==1 )  )
                                 } else {
                                    istatus[i] <- 0
                                 }
                             }
                             xnames <- ""
                             inames <- 1
                             for ( i in 1:length(names) ) {
                                 if ( !isTRUE(all.equal(istatus[i],0)) ) {
                                    if ( isTRUE(all.equal(length( which(xnames==paste( istatus[i],sep="") ) ),0)) ) {
                                       xnames[inames] <- paste( istatus[i],sep="" )
                                       inames <- inames+1
                                    }
                                 } else {
                                    xnames[inames] <- names[ i ]
                                    inames <- inames+1
                                 }
                             }
                             icol <- length( xnames )           
                             xx <- matrix( 0, ncol=icol, nrow=icol )
                             xxavg <- matrix( 0, ncol=icol, nrow=icol )
                             
                             ####################################################
                             if ( icol >1 ) {
                                for ( i in 1:(icol-1) ) {
                                    for ( j in (i+1):icol ) {
                                        if ( (isTRUE(all.equal(length( which(names==xnames[i]) ),0)) ) & (isTRUE(all.equal(length( which(names==xnames[j]) ),0)) ) ) {
                                           iindex <- which( xcluster[ as.numeric( xnames[i] ),1:n]==1 )
                                           jindex <- which( xcluster[ as.numeric( xnames[j] ),1:n]==1 )
                                        }
                                        if ( (isTRUE(all.equal(length( which(names==xnames[i]) ),0)) ) & (!isTRUE(all.equal( length( which(names==xnames[j]) ),0)) ) ) {
                                           iindex <- which( xcluster[ as.numeric( xnames[i] ),1:n]==1 )
                                           jindex <- which( names==xnames[j] )
                                        }
                                        if ( (!isTRUE(all.equal(length( which(names==xnames[i]) ),0)) ) & (isTRUE(all.equal(length( which(names==xnames[j]) ),0)) ) ) {
                                           iindex <- which( names==xnames[i] )
                                           jindex <- which( xcluster[ as.numeric( xnames[j] ),1:n]==1 )
                                        }
                                        if ( (!isTRUE(all.equal(length( which(names==xnames[i]) ),0)) ) & (!isTRUE(all.equal(length( which(names==xnames[j]) ),0)) ) ) {
                                           iindex <- which( names==xnames[i] )
                                           jindex <- which( names==xnames[j] )
                                        }
                                        temp <- 0
                                        tempavg <- 0
                                        for ( ii in iindex ) {
                                            for ( jj in jindex ) {
                                                temp <- max( temp, x[ii,jj] )
                                                tempavg <- tempavg + x[ii,jj]
                                            }
                                        }
                                        xx[i,j] <- temp
                                        xxavg[i,j] <- tempavg / ( length(iindex)*length(jindex) )
                                    }
                                }   
                             }
                             for ( i in 1: icol ) {
                                 xx[i,i] <- 1
                                 xxavg[i,i] <- 1
                             }
                             if ( icol >1 ) {
                                for ( i in 2: icol ) {
                                    for ( j in 1:(i-1) ) {
                                        xx[i,j] <- xx[ j,i]
                                        xxavg[i,j] <- xxavg[j,i]
                                    }
                                }
                             }
                             if ( isTRUE(all.equal(icol,2)) ) {
                                if ( length(which(xcluster[(cindex-1),]==0))>1 ) {
                                   lastpair <- paste( paste("cluster",cindex-1,sep=""), "--", paste( "cluster",which( xcluster[,which(xcluster[(cindex-1),]==0)[1]]==1 ),sep=""), paste("[label=", format( xxavg[1,2],digits=2 ),",","dir=forward];",sep=""),sep=" " )
                                } else {     
                                   lastpair <- paste( paste("cluster",cindex-1,sep=""), "--", names[which(xcluster[(cindex-1),]==0)], paste("[label=", format( xxavg[1,2],digits=2 ),",","dir=forward];",sep=""),sep=" " )
                                }
                             }
                       }
                       
                       # Generating the scripts for the resulting graph
                       if ( length(which(xcluster[(cindex-1),]==0))>1 ) {
                          gScript( which( xcluster[,which(xcluster[(cindex-1),]==0)[1]]==1 ), names )
                       } else {
                          if ( isTRUE(all.equal(length(which(xcluster[(cindex-1),]==0)),1)) ) {
                             write( paste( names[ which(xcluster[(cindex-1),]==0) ], sep=""), length(  paste( names[ which(xcluster[(cindex-1),]==0) ], sep="") ),file= paste( "graph",".txt",sep=""), append=T )
                          }
                       }     
                       gScript( cindex-1, names )
                       gScript0( cindex-1, names, xcluster, x )
                       if ( isTRUE(all.equal(icol,2)) ) {
                          write( lastpair, length( lastpair ), file=paste( "graph",".txt",sep=""), append=T )
                       }
                       cMatrix <- gMatrix( cindex-1, names, xcluster, x )
                       write( paste( "}",sep="" ), length( paste( "}",sep="" ) ), file=paste( "graph",".txt",sep=""), append=T )
                       if ( isTRUE(all.equal(icol,2)) ) {
                          write( paste( names[ which( x[ ,which(xcluster[(cindex-1),1:length(names)]==0) ]==max(x[ -which(xcluster[(cindex-1),]==0),which(xcluster[(cindex-1),]==0) ]) )[1] ], "->", names[ which(xcluster[(cindex-1),1:length(names)]==0) ], "[len=", paste( 3*(1.5-max(x[ -which(xcluster[(cindex-1),]==0),which(xcluster[(cindex-1),]==0) ])),sep="" ),",label=", paste(max(x[ -which(xcluster[(cindex-1),]==0),which(xcluster[(cindex-1),]==0) ]),sep=""),"];",sep="" ), length(paste( names[ which( x[ ,which(xcluster[(cindex-1),]==0) ]==max(x[ -which(xcluster[(cindex-1),]==0),which(xcluster[(cindex-1),]==0) ]) )[1] ], "->", names[ which(xcluster[(cindex-1),]==0) ], "[len=", paste( 3*(1.5-max(x[ -which(xcluster[(cindex-1),]==0),which(xcluster[(cindex-1),]==0) ])),sep="" ),",label=", paste(max(x[ -which(xcluster[(cindex-1),]==0),which(xcluster[(cindex-1),]==0) ]),sep=""),"];",sep="" )),file=paste( "graph001",".txt",sep=""), append=T ) 
                       }
                       write( paste( "}",sep="" ), length( paste( "}",sep="" ) ), file=paste( "graph001",".txt",sep=""), append=T ) 
                       thresholdlist <- thresholdlist[ 1: length( which( thresholdlist != 0 )  ) ]
                       hc <- rep( 0,3)
                       hc[1] <- cor( as.dist(1-x), cophenetic( hclust( as.dist(1-x), method="average" ) ) )
                       hc[2] <- cor( as.dist(1-x), cophenetic( hclust( as.dist(1-x), method="single" ) ) )
                       hc[3] <- cor( as.dist(1-x), cophenetic( hclust( as.dist(1-x), method="complete" ) ) )
                       rest <- list( 0, cor( as.dist(x), as.dist(cMatrix) ), hc, cMatrix )
                       file.rename( "graph.txt", "Hgraph.txt" )
                       file.rename( "graph001.txt", "Igraph.txt" )
                       return( rest )
                    } else {
                       return( 1 )
                    }
        }
        
        # constr function: convert a concept to a string
        constr <- function ( concept ) {
                  k <- 1
                  xbefore <- ""
                  xtemp <- substr( concept, k, k )
                  if ( xtemp == ":" ) {
                     xbefore <- "-"
                  } else {
                     while ( xtemp != ":" ) {
                           xtemp1 <- ""
                           while ( xtemp != " " ) {
                                 k <- k+1
                                 xtemp1 <- paste( xtemp1, xtemp, sep="" )
                                 xtemp <- substr( concept, k, k )
                           }
                           if ( xbefore[1] == "" ) {
                              xbefore <- xtemp1
                           } else {
                              if ( xtemp1 != "" ) {
                                 xbefore <- c( xbefore, xtemp1 )
                              }
                           }
                           k <- k+1
                           xtemp <- substr( concept, k, k )
                     }
                  }
                  if ( xbefore[1] != "-" ) {
                     xbefore <- c( xbefore, "-" )
                  }
                  k <- k+1
                  xtemp <- substr( concept, k, k )
                  while ( xtemp != "-" ) {
                        xtemp1 <- ""
                        while ( xtemp != " " ) {
                              k <- k+1
                              xtemp1 <- paste( xtemp1, xtemp, sep="" )
                              xtemp <- substr( concept, k, k )
                        }
                        if ( xbefore[1] == "" ) {
                           xbefore <- xtemp1
                        } else {
                           if ( xtemp1 != "" ) {
                              xbefore <- c( xbefore, xtemp1 )
                           }
                        }
                        k <- k+1
                        xtemp <- substr( concept, k, k )
                  }
                  return( xbefore )
        }
        
        # constr0 function: convert first part of a concept to a string
        constr0 <- function ( concept ) {
                   k <- 1
                   xbefore <- ""
                   xxtemp <- ""
                   xtemp <- substr( concept, k, k )
                   while ( xtemp != "[" ) {
                         xxtemp <- ""
                         while ( xtemp != " " ) {
                               k <- k+1
                               xxtemp <- paste( xxtemp,xtemp,sep="" )
                               xtemp <- substr( concept, k, k )
                         }
                         xbefore <- c( xbefore, xxtemp )
                         k <- k+1
                         xtemp <- substr( concept, k, k )
                   }
                   return( xbefore[2:length(xbefore)] )
        }
        
        # thresholdTEST function: test whether a given threshold is appropriate
        thresholdTEST <- function ( threshold, names, x ) {
                         res <- 0
                         # Threshold the data matrix
                         xbin <- x
                         xbin[ which( xbin<threshold ) ] <- 0
                         xbin[ which( xbin>=threshold ) ] <- 1
                         # Apply FCA to construct the concept lattice
                         if ( file.exists( "example.oal" )==TRUE ) {
                            file.rename( "example.oal", "exampleoriginal.oal" )
                         }
                         if ( file.exists( "example.txt" )==TRUE ) {
                            file.rename( "example.txt", "exampleoriginal.txt" )
                         }
                         for ( i in 1:(dim(xbin)[1]) ) {
                             rowname <- paste( names[i], ":", sep="" )
    	                       colname <- ""
	                           for ( j in 1:(dim(xbin)[2]) ) {
		                             if ( isTRUE(all.equal(xbin[i,j] , 1)) ) {
                                    colname <- paste( colname, names[j], ";", sep="" )
                                 }
                             }
	                           write( paste( rowname, colname, sep="" ), file=paste("example",".oal",sep=""), append=T )
                         }
                         system( "ruby multistagefca.rb" )
                         file.rename( "example.dot", "example.txt" )
                         # Define the cluster informative matrix
                         xcluster <- matrix( 0, ncol=length(names)+1, nrow=50 )
                         tlevel <- 1
                         cindex <- 1
                         # Take all the concepts and put them in increasing order
                         concepts <- as.matrix( read.csv( "example.txt", header=F ) )
                         index1 <- grep( "--", concepts )
                         iorder <- rep( 0, length(index1) )
                         for ( i in 1:length(index1) ) {
                             xbefore <- constr( concepts[ index1[i] ] )
                             iloc <- which( xbefore=="-" )
                             if ( isTRUE(all.equal(iloc,1)) ) {
                                iorder[ i ] <- 0
                             } else {
                                if ( (isTRUE(all.equal(iloc , 2))) & ( isTRUE(all.equal(( length( xbefore ) - iloc ),1)) ) ) {
                                   iorder[ i ] <- 0
                                } else {
                                   if ( ( iloc-1 ) > ( length( xbefore ) - iloc ) ) {
                                      iorder [i] <- 0
                                   } else {
                                      iorder[i] <- iloc-1
                                   }
                                }
                             }
                         }
                         index2 <- which( iorder != 0 )
                         index1 <- index1[ index2 ]
                         iorder <- iorder[ index2 ]
                         index0 <- index1[ sort( iorder, index.return=T )$ix ]
                         conindex <- rep( 0, length(concepts) )
                         # Given threshold is too large     
                         if ( isTRUE(all.equal(length(index0),0)) ) {
                            res <- 2
                            return( res )
                         }
                         # Test if it is too small
                         for ( i in index0 ) {
                             if ( isTRUE(all.equal(conindex[i] ,0)) ) {
                                conindex[i] <- 1
                                xbefore <- constr( concepts[i] )
                                iloc <- which( xbefore=="-" )
                                if ( isTRUE(all.equal(iloc,2)) ) {
                                   ifirst <- xbefore[iloc-1]
                                   for ( j in index0 ) {
                                       if ( isTRUE(all.equal(conindex[j],0)) ) {
                                          newxbefore <- constr( concepts[j] )
                                          newloc <- which( newxbefore=="-" )
                                          if ( (length( which( newxbefore[1:(newloc-1)]==ifirst ) )>0 ) & ( isTRUE(all.equal(newloc,iloc)) ) ) {
                                             conindex[j] <- 1
                                          }
                                          if ( (length( which( newxbefore[1:(newloc-1)]==ifirst ) )>0 ) & ( newloc>iloc ) ) {
                                             conindex[ j ] <- 1
                                          }
                                       }
                                   }
                                } else {
                                   if ( isTRUE(all.equal((iloc-1),(length(xbefore)-iloc))) ) {
                                      iloc <- iloc
                                   } else {
                                      res <- 1
                                      return( res )
                                   }
                                }
                             }
                         }
                         for ( i in index0 ) {
                             for ( j in index0 ) {
                                 if ( !isTRUE(all.equal(i ,j)) ) {
                                    ibefore <- constr( concepts[i] )
                                    jbefore <- constr( concepts[j] )
                                    iloc <- which( ibefore=="-" )
                                    jloc <- which( jbefore=="-" )
                                    isame <- 0
                                    for ( ii in 1: (iloc-1) ) {
                                        if ( isTRUE(all.equal(length( which( jbefore[1:(jloc-1)]==ibefore[ii] ) ), 0)) ) {
                                           isame <- 1
                                        }
                                    }
                                    jsame <- 0
                                    for ( jj in 1: (jloc-1) ) {
                                        if ( isTRUE(all.equal(length( which( ibefore[1:(iloc-1)]==jbefore[jj] ) ), 0)) ) {
                                           jsame <- 1
                                        }
                                    }
                                    if ( (isTRUE(all.equal(isame ,1))) & (isTRUE(all.equal( jsame,1))) ) {
                                       for ( kk in max( which( index0==i ), which( index0==j ) ): length(index0) ) {
                                           k <- index0[kk]
                                           kbefore <- constr( concepts[k] )
                                           kloc <- which( kbefore=="-" )
                                           if ( kloc> max(iloc,jloc-1) ) {
                                              iiisame <- 0
                                              for ( iii in 1:(iloc-1) ) {
                                                  if ( !isTRUE(all.equal(length( which( kbefore[1:(kloc-1)]==ibefore[iii] ) ) , 0)) ) {
                                                     iiisame <- iiisame + 1
                                                  }
                                              }
                                              jjjsame <- 0
                                              for ( jjj in 1:(jloc-1) ) {                                            
                                                  if ( !isTRUE(all.equal(length( which( kbefore[1:(kloc-1)]==jbefore[jjj] ) ) , 0)) ) {
                                                     jjjsame <- jjjsame + 1
                                                  }
                                              }
                                              if ( (isTRUE(all.equal(iiisame,(iloc-1)))) & (isTRUE(all.equal( jjjsame,(jloc-1))) ) ) {
                                                 res <- 1
                                                 return( res )
                                              }
                                           }
                                       }
                                    }
                                 }
                             }
                         }
                         for ( i in index0 ) {
                             ibefore <- constr( concepts[i] )
                             iloc <- which( ibefore=="-" )
                             if ( isTRUE(all.equal(iloc,2)) ) {
                                for ( j in index0 ) {
                                    jbefore <- constr( concepts[j] )
                                    jloc <- which( jbefore=="-" )
                                    if ( jloc > iloc ) {
                                       if ( !isTRUE(all.equal(length( which( jbefore[1:(jloc-1)]==ibefore[1] ) ), 0)) ) {
                                          for ( ii in index0 ) {
                                              if (  (!isTRUE(all.equal(ii ,i))) & (!isTRUE(all.equal(ii , j))) ) {
                                                 iibefore <- constr( concepts[ii] )
                                                 iiloc <- which( iibefore == "-" )
                                                 if ( (iiloc>iloc) & ( iiloc>jloc ) & (!isTRUE(all.equal(length(which( iibefore[1:(iiloc-1)]==ibefore[1] ) ),0)) ) ) {
                                                    iiicount <- 0
                                                    for ( jjj in 1:(jloc-1) ) {
                                                        if ( !isTRUE(all.equal(jjj , which( jbefore[1:(jloc-1)]==ibefore[1] )))  ) {
                                                           if ( isTRUE(all.equal(length( which( iibefore[1:(iiloc-1)]==jbefore[jjj] ) ),0)) ) {
                                                              iiicount <- 1
                                                           }
                                                        }
                                                    }
                                                    if ( isTRUE(all.equal(iiicount , 0)) ) {
                                                       res <- 1
                                                       return( res ) 
                                                    } 
                                                 }
                                              }
                                          }
                                       }
                                    }
                                }
                             }
                         }
                         if ( file.exists( "example.oal" )==TRUE ) {
                            file.remove( "example.oal" )
                         }
                         if ( file.exists( "example.txt" )==TRUE ) {
                            file.remove( "example.txt" )
                         }
                         if ( file.exists( "exampleoriginal.oal" )==TRUE ) {
                            file.rename( "exampleoriginal.oal", "example.oal" )
                         }
                         if ( file.exists( "exampleoriginal.txt" )==TRUE ) {
                            file.rename( "exampleoriginal.txt", "example.txt" )
                         }
                         return( res )
        }
        
        # gScript function: generating the script based on the thresholding results
        gScript <- function( ccindex, names ) {
                   ggindex <- rep( 0, ccindex )
                   filename <- paste( "graph",ccindex,".txt",sep="" )
                   xconcepts <- as.matrix( read.csv( paste(filename,sep=""), header=F ) )
                   concepts <- rep( 0, dim(xconcepts)[1] )
                   clusterinfo <- ""
                   if ( dim(xconcepts)[[2]]>1 ) {
                      for ( ii in 2: (length( concepts )-1) ) {
                          concepts[ii] <- paste( xconcepts[ii,1], ",", xconcepts[ii,2], sep="" )
                      }
                      concepts[1] <- paste( xconcepts[1,1], xconcepts[1,2], sep="" )
                      concepts[length( concepts )] <- paste( xconcepts[length( concepts ),1], xconcepts[length( concepts ),2], sep="" )
                   } else {
                      concepts <- xconcepts
                   }
                   index0 <- grep( "--", concepts )
                   write( concepts[1], length( concepts[1]), file=paste( "graph",".txt",sep=""), append=T )
                   for ( ii in index0 ) {
                       xbefore <- constr0( concepts[ii] )
                       if ( !isTRUE(all.equal(length( which( names==xbefore[1] ) ),0)) ) {
                          if ( !(!isTRUE(all.equal(length( which( names==xbefore[3] ) ),0))) ) {
                             write( paste( xbefore[1], sep="" ), length( paste( xbefore[1], sep="" ) ), file=paste( "graph",".txt",sep=""), append=T )
                          } else {
                             write( paste( concepts[ii], sep="" ), length( paste( concepts[ii], sep="" ) ), file=paste( "graph",".txt",sep=""), append=T ) 
                          }
                       } else {
                          if ( (!isTRUE(all.equal(length( which( names==xbefore[3] ) ),0))) ) {
                             write( paste( xbefore[3], sep="" ), length( paste( xbefore[3], sep="" ) ), file=paste( "graph",".txt",sep=""), append=T ) 
                          }    
                       }
                   }
                   for ( ii in index0 ) {
                       xbefore <- constr0( concepts[ii] )
                       if ( (isTRUE(all.equal(length( which( names==xbefore[1] ) ),0))) ) {
                          if ( isTRUE(all.equal(ggindex[ as.numeric( xbefore[1] ) ],0)) ) {
                             gScript( as.numeric( xbefore[1] ), names )
                             ggindex[ as.numeric( xbefore[1] ) ] <- 1       
                          }
                       }
                       if ( (isTRUE(all.equal(length( which( names==xbefore[3] ) ),0))) ) {
                          if ( isTRUE(all.equal(ggindex[ as.numeric( xbefore[3] ) ],0)) ) {                 
                             gScript( as.numeric( xbefore[3] ), names )
                             ggindex[ as.numeric( xbefore[3] ) ] <- 1   
                          }
                       }
                   }
                   write( paste( "label=Cluster", ccindex, sep="" ), length(paste( "label=Cluster", ccindex, sep="" ) ) , file=paste( "graph",".txt",sep=""), append=T )            
                   write( paste( "}",sep="" ), length( paste( "}",sep="" ) ), file=paste( "graph",".txt",sep=""), append=T ) 
                   return( 0 )          
        }
        
        # gScript0 function: generating addition information based on the thresholding results
        gScript0 <- function( ccindex, names, xcluster, x ) {
                    ggindex <- rep( 0, ccindex )
                    filename <- paste( "graph",ccindex,".txt",sep="" )
                    xconcepts <- as.matrix( read.csv( paste(filename,sep=""), header=F ) )
                    concepts <- rep( 0, dim(xconcepts)[1] )
                    clusterinfo <- ""
                    if ( dim(xconcepts)[[2]]>1 ) {
                       for ( ii in 2: (length( concepts )-1) ) {
                           concepts[ii] <- paste( xconcepts[ii,1], ",", xconcepts[ii,2], sep="" )
                       }
                       concepts[1] <- paste( xconcepts[1,1], xconcepts[1,2], sep="" )
                       concepts[length( concepts )] <- paste( xconcepts[length( concepts ),1], xconcepts[length( concepts ),2], sep="" )
                    } else {
                       concepts <- xconcepts
                    }
                    index0 <- grep( "--", concepts )
                    for ( ii in index0 ) {
                        i <- 1
                        while ( substr( concepts[ii],i,i)!="[" ) {
                              i <- i + 1
                        }
                        j <- 1
                        while ( substr( concepts[ii],j,j)!=";" ) {
                              j <- j + 1
                        }
                        xbefore <- constr0( concepts[ii] )
                        if ( !isTRUE(all.equal(length( which( names==xbefore[1] ) ),0)) ) {
                           if ( !(!isTRUE(all.equal(length( which( names==xbefore[3] ) ),0))) ) {
                              iii <- paste( xbefore[1], xbefore[2], paste( "cluster",xbefore[3],sep=""), substr( concepts[ii] , i,j ), sep=" " )
                              write( iii, length( iii ), file=paste( "graph",".txt",sep=""), append=T )
                              iii00 <- which( names==xbefore[1] )
                              jjj00 <- which( xcluster[ as.numeric( xbefore[3] ), 1:length(names) ]!=0 )
                              write( paste( names[iii00], "->", names[ which( x[iii00,]==max(x[iii00,jjj00]) )[1] ],"[len=", paste(3*(1.5-max(x[iii00,jjj00])),sep=""), ",", "label=", paste(max(x[iii00,jjj00]),sep=""),",dir=back];",sep="" ), length(paste( names[iii00], "->", names[ which( x[iii00,]==max(x[iii00,jjj00]) )[1] ],"[len=", paste(3*(1.5-max(x[iii00,jjj00])),sep=""), ",", "label=", paste(max(x[iii00,jjj00]),sep=""),",dir=back];",sep="" )),file=paste( "graph001",".txt",sep=""), append=T )                               
                           } else {
                              iii00 <- which( names==xbefore[1] )
                              jjj00 <- which( names==xbefore[3] )
                              if ( isTRUE(all.equal(length(index0),1)) ) {
                                 write( paste( names[iii00], "->", names[jjj00],"[len=", paste(3*(1.5-x[iii00,jjj00]),sep=""), ",", "label=", paste(x[iii00,jjj00],sep=""),",dir=both];",sep="" ), length( paste( names[iii00], "->", names[jjj00],"[len=", paste(3*(1.5-x[iii00,jjj00]),sep=""), ",", "label=", paste(x[iii00,jjj00],sep=""),",dir=both];",sep="" ) ),file=paste( "graph001",".txt",sep=""), append=T )  
                              } else {
                                 write( paste( names[iii00], "->", names[jjj00],"[len=", paste(3*(1.5-x[iii00,jjj00]),sep=""), ",", "label=", paste(x[iii00,jjj00],sep=""),"];",sep="" ), length( paste( names[iii00], "->", names[jjj00],"[len=", paste(3*(1.5-x[iii00,jjj00]),sep=""), ",", "label=", paste(x[iii00,jjj00],sep=""),"];",sep="" ) ),file=paste( "graph001",".txt",sep=""), append=T )  
                              }
                           }
                        } else {
                           if ( !(!isTRUE(all.equal(length( which( names==xbefore[3] ) ),0))) ) {                 
                              iii <- paste( paste("cluster",xbefore[1],sep=""), xbefore[2], paste( "cluster",xbefore[3],sep=""), substr( concepts[ii] , i,j ), sep=" " )
                              write( iii, length( iii ), file=paste( "graph",".txt",sep=""), append=T )
                              iii00 <- which( xcluster[ as.numeric( xbefore[1] ), 1:length(names) ]!=0 )
                              jjj00 <- which( xcluster[ as.numeric( xbefore[3] ), 1:length(names) ]!=0 )
                              for ( ij00 in jjj00 ) {
                                  write( paste( names[ which( x[,ij00]==max(x[iii00,ij00]) )[1] ], "->", names[ij00],"[len=", paste(3*(1.5-x[which( x[,ij00]==max(x[iii00,ij00]) )[1],ij00]),sep=""), ",", "label=", paste(x[which( x[,ij00]==max(x[iii00,ij00]) )[1],ij00],sep=""),"];",sep="" ), length( paste( names[which( x[,ij00]==max(x[iii00,ij00]) )[1]], "->", names[ij00],"[len=", paste(3*(1.5-x[which( x[,ij00]==max(x[iii00,ij00]) )[1],ij00]),sep=""), ",", "label=", paste(x[which( x[,ij00]==max(x[iii00,ij00]) )[1],ij00],sep=""),"];",sep="" ) ),file=paste( "graph001",".txt",sep=""), append=T )  
                              }
                           } else {
                              iii <- paste( paste("cluster",xbefore[1],sep=""), xbefore[2], xbefore[3], substr( concepts[ii] , i,j ), sep=" " )
                              write( iii, length( iii ), file=paste( "graph",".txt",sep=""), append=T )
                              jjj00 <- which( names==xbefore[3] )
                              iii00 <- which( xcluster[ as.numeric( xbefore[1] ), 1:length(names) ]!=0 )
                              write( paste( names[which( x[,jjj00]==max(x[iii00,jjj00]) )[1]], "->", names[jjj00],"[len=", paste(3*(1.5-x[which( x[,jjj00]==max(x[iii00,jjj00]) )[1],jjj00]),sep=""), ",", "label=", paste(x[which( x[,jjj00]==max(x[iii00,jjj00]) )[1],jjj00],sep=""),"];",sep="" ), length( paste( names[which( x[,jjj00]==max(x[iii00,jjj00]) )[1]], "->", names[jjj00],"[len=", paste(3*(1.5-x[which( x[,jjj00]==max(x[iii00,jjj00]) )[1],jjj00]),sep=""), ",", "label=", paste(x[which( x[,jjj00]==max(x[iii00,jjj00]) )[1],jjj00],sep=""),"];",sep="" ) ),file=paste( "graph001",".txt",sep=""), append=T )  
                           }
                        }
                    }
                    for ( ii in index0 ) {
                        xbefore <- constr0( concepts[ii] )
                        if ( (isTRUE(all.equal(length( which( names==xbefore[1] ) ),0))) ) {
                           if ( isTRUE(all.equal(ggindex[ as.numeric( xbefore[1] ) ],0)) ) {
                              gScript0( as.numeric( xbefore[1] ), names, xcluster, x )
                              ggindex[ as.numeric( xbefore[1] ) ] <- 1       
                           }
                        }
                        if ( (isTRUE(all.equal(length( which( names==xbefore[3] ) ),0))) ) {
                           if ( isTRUE(all.equal(ggindex[ as.numeric( xbefore[3] ) ],0)) ) {                 
                              gScript0( as.numeric( xbefore[3] ), names, xcluster, x )
                              ggindex[ as.numeric( xbefore[3] ) ] <- 1   
                           }
                        }
                    }
                    return( 0 )          
        }
        
        # gMatrix function
        gMatrix <- function(ccindex, names, xcluster, x) {
                   cMatrix <- matrix(0, ncol=dim(x)[1], nrow=dim(x)[1])
                   for ( iindex in 1:ccindex ) {
                       filename <- paste( "graph",iindex,".txt",sep="" )
                       xconcepts <- as.matrix( read.csv( paste(filename,sep=""), header=F ) )
                       concepts <- rep( 0, dim(xconcepts)[1] )
                       if ( dim(xconcepts)[[2]]>1 ) {
                          for ( ii in 2: (length( concepts )-1) ) {
                              concepts[ii] <- paste( xconcepts[ii,1], ",", xconcepts[ii,2], sep="" )
                          }
                          concepts[1] <- paste( xconcepts[1,1], xconcepts[1,2], sep="" )
                          concepts[length( concepts )] <- paste( xconcepts[length( concepts ),1], xconcepts[length( concepts ),2], sep="" )
                       } else {
                          concepts <- xconcepts
                       }
                       index0 <- grep( "--", concepts )
                       for ( ii in index0 ) {
                           xbefore <- constr0( concepts[ii] )
                           if ( !isTRUE( all.equal(length( which( names==xbefore[1] ) ),0) ) ) {
                              if ( !(!isTRUE( all.equal(length( which( names==xbefore[3] ) ),0) )) ) {
                                 rowindex <- which( names == xbefore[1] )
                                 colindex <- which( xcluster[ as.numeric( xbefore[3] ), 1:(length(names)) ]==1 )
                                 iiitemp <- 0
                                 for ( iii in colindex ) {
                                     iiitemp <- iiitemp + x[ rowindex, iii ]
                                 }
                                 cMatrix[ rowindex, colindex ] <- iiitemp/length(colindex)
                                 for ( jjj in colindex ) {
                                     cMatrix[ rowindex,jjj ] <- x[ rowindex,jjj ]
                                 } 
                              } else {
                                 rowindex <- which( names == xbefore[1] )
                                 colindex <- which( names == xbefore[3] )
                                 cMatrix[ rowindex, colindex ] <- x[ rowindex, colindex ]
                              }
                           } else {
                              if ( !(!isTRUE( all.equal(length( which( names==xbefore[3] ) ),0) )) ) {
                                 rowindex <- which( xcluster[ as.numeric( xbefore[1] ), 1:(length(names)) ]==1 )
                                 colindex <- which( xcluster[ as.numeric( xbefore[3] ), 1:(length(names)) ]==1 )
                                 iiitemp <- 0
                                 for ( iii in rowindex ) {
                                     for ( jjj in colindex ) {
                                         iiitemp <- iiitemp + x[ iii , jjj ]
                                     }
                                 }
                                 cMatrix[ rowindex, colindex ] <- iiitemp/ ( length(rowindex)*length(colindex) )                                 
                                 for ( jjj in colindex ) {
                                     cMatrix[ which( x[,jjj]==max(x[rowindex,jjj]) ), jjj ] <- max( x[rowindex,jjj] )
                                 }
                              } else {
                                 rowindex <- which( xcluster[ as.numeric( xbefore[1] ), 1:(length(names)) ]==1 )
                                 colindex <- which( names == xbefore[3] )
                                 iiitemp <- 0
                                 for ( iii in rowindex ) {
                                     iiitemp <- iiitemp + x[ iii , colindex ]
                                 }
                                 cMatrix[ rowindex, colindex ] <- iiitemp/length(rowindex)
                                 cMatrix[ which( x[,colindex]==max(x[rowindex,colindex]) ), colindex ] <- max( x[rowindex,colindex] )
                              }
                           }
                       }
                       fixvec <- rep( 0, 1 )
                       for ( fixi in index0 ) {
                           fixibefore <- constr0( concepts[fixi] )
                           if ( isTRUE( all.equal(length( which( fixvec==fixibefore[1] ) ) , 0) ) ) {
                              fixvec <- c( fixvec, fixibefore[1] )
                           }
                           if ( isTRUE( all.equal(length( which( fixvec==fixibefore[3] ) ) , 0) ) ) {
                              fixvec <- c( fixvec, fixibefore[3] )
                           } 
                       }
                       if ( !isTRUE( all.equal((length( fixvec)-1) , length( index0 )) ) ) {
                          if ( length(index0) >1 ) {
                             for ( ii in index0 ) {
                                 for ( jj in index0 ) {
                                     if ( !isTRUE( all.equal(ii , jj) ) ) {
                                        iibefore <- constr0( concepts[ii] )
                                        jjbefore <- constr0( concepts[jj] )
                                        ii3 <- iibefore[3]
                                        jj3 <- jjbefore[3]
                                        if ( !isTRUE( all.equal(length( which( names==ii3 ) ),0) ) ) {
                                           if ( !(!isTRUE( all.equal(length( which( names==jj3 ) ),0) )) ) {
                                              rowindex <- which( names == ii3 )
                                              colindex <- which( xcluster[ as.numeric( jj3 ), 1:(length(names)) ]==1 )
                                              iiitemp <- 0
                                              for ( iii in colindex ) {
                                                  iiitemp <- iiitemp + x[ rowindex, iii ]
                                              }
                                              if ( ii<jj ) {
                                                 write( paste( ii3, iibefore[2], paste( "cluster",jj3,sep=""),"[label=",paste(format(iiitemp/length(colindex),digits=2),",",sep=""),"color=red,fontcolor=red];",sep="" ), length( paste( ii3, iibefore[2], paste( "cluster",jj3,sep=""),"[label=",paste(format(iiitemp/length(colindex),digits=2),",",sep=""),"color=red,fontcolor=red];",sep="" ) ), file=paste( "graph",".txt",sep=""), append=T ) 
                                              }
                                              cMatrix[ rowindex, colindex ] <- iiitemp/length(colindex)
                                           } else {
                                              rowindex <- which( names == ii3 )
                                              colindex <- which( names == jj3 )
                                              if ( ii<jj ) {
                                                 write( paste( ii3, iibefore[2], jj3,"[label=",paste(format(x[ rowindex, colindex ],digits=2),",",sep=""),"color=red,fontcolor=red];",sep="" ), length( paste( ii3, iibefore[2], jj3,"[label=",paste(format(x[ rowindex, colindex ],digits=2),",",sep=""),"color=red,fontcolor=red];",sep="" ) ), file=paste( "graph",".txt",sep=""), append=T ) 
                                              }                                                 
                                              cMatrix[ rowindex, colindex ] <- x[ rowindex, colindex ]
                                           }
                                        } else {
                                           if ( !(!isTRUE( all.equal(length( which( names==jj3 ) ),0) )) ) {
                                              rowindex <- which( xcluster[ as.numeric( ii3 ), 1:(length(names)) ]==1 )
                                              colindex <- which( xcluster[ as.numeric( jj3 ), 1:(length(names)) ]==1 )
                                              iiitemp <- 0
                                              for ( iii in rowindex ) {
                                                  for ( jjj in colindex ) {
                                                      iiitemp <- iiitemp + x[ iii , jjj ]
                                                  }
                                              }
                                              if ( ii<jj ) {
                                                 write( paste( paste( "cluster",ii3,sep=""), iibefore[2], paste( "cluster",jj3,sep=""),"[label=",paste(format(iiitemp/ ( length(rowindex)*length(colindex) ),digits=2),",",sep=""),"color=red,fontcolor=red];",sep="" ), length( paste( paste( "cluster",ii3,sep=""), iibefore[2], paste( "cluster",jj3,sep=""),"[label=",paste(format(iiitemp/ ( length(rowindex)*length(colindex) ),digits=2),",",sep=""),"color=red,fontcolor=red];",sep="" ) ), file=paste( "graph",".txt",sep=""), append=T ) 
                                              }                                                 
                                              cMatrix[ rowindex, colindex ] <- iiitemp/ ( length(rowindex)*length(colindex) )
                                           } else {
                                              rowindex <- which( xcluster[ as.numeric( ii3 ), 1:(length(names)) ]==1 )
                                              colindex <- which( names == jj3 )
                                              iiitemp <- 0
                                              for ( iii in rowindex ) {
                                                  iiitemp <- iiitemp + x[ iii , colindex ]
                                              }
                                              if ( ii<jj ) {
                                                 write( paste( paste( "cluster",ii3,sep=""), iibefore[2], jj3,"[label=",paste(format(iiitemp/length(rowindex),digits=2),",",sep=""),"color=red,fontcolor=red];",sep="" ), length( paste( paste( "cluster",ii3,sep=""), iibefore[2], jj3,"[label=",paste(format(iiitemp/length(rowindex),digits=2),",",sep=""),"color=red,fontcolor=red];",sep="" ) ), file=paste( "graph",".txt",sep=""), append=T ) 
                                              }                                                 
                                              cMatrix[ rowindex, colindex ] <- iiitemp/length(rowindex)
                                           }
                                        }
                                     }
                                 }
                             }
                          }  
                       }
                   }
                   lastone <- which( xcluster[ccindex,1:length(names)]==0 )
                   if ( !isTRUE( all.equal(length(lastone),0) ) ) {
                      if ( length( lastone )>1 ) {
                         cMatrix[ lastone[1], lastone[2] ] <- x[ lastone[1], lastone[2] ]
                         iiitemp <- 0
                         for ( iii in lastone ) {
                             for ( jjj in which( xcluster[ccindex,1:length(names)]!=0 ) ) {
                                 iiitemp <- iiitemp + x[iii,jjj]
                             }
                         }
                         cMatrix[ lastone, which( xcluster[ccindex,1:length(names)]!=0 ) ] <- iiitemp / ( length(lastone)*length(which( xcluster[ccindex,1:length(names)]!=0 )) )
                      } else {
                         iiitemp <- 0
                         for ( i in 1:length(names) ) {
                             if ( !isTRUE( all.equal(i , lastone) ) ) {
                                iiitemp <- iiitemp + x[ lastone, i ] 
                             }
                         }
                         cMatrix[ lastone, which( names != names[lastone] ) ] <- iiitemp / length( which( names != names[lastone] ) )
                         cMatrix[ lastone, which( x[,lastone]==max( x[-lastone,lastone] ) ) ] <-  max( x[-lastone,lastone] )
                      }
                   }
                   for ( i in 1:(dim(x)[1]) ) {
                       for ( j in 1:(dim(x)[1]) ) {
                           if ( isTRUE( all.equal(i,j) ) ) {
                              cMatrix[ i,j ] <- 1
                           }
                           if ( isTRUE( all.equal(cMatrix[ i,j ] , 0) ) ) {
                              cMatrix[ i,j ] <- cMatrix[ j,i ]
                           }
                       }
                   }
                   return( cMatrix )
        }

        ################
        ################
        cluster.cat <- function(string.graph,space.index,iline.index,file.graph){
                       return.index <- rep(0,3)
                       subgraph.locate <- regexpr("subgraph",string.graph,fixed=T)[1]
                       label.locate <- regexpr(" [",string.graph,fixed=T)[1]
                       start.label.locate <- regexpr("label=Cluster",string.graph,fixed=T)[1]
                       if (label.locate>1) {
                          dir.index <- regexpr("dir=",string.graph,fixed=T)
                          temp.locate <- regexpr(" -- ",string.graph,fixed=T)
                          if (dir.index>1) {
                             cat(paste(replicate(space.index,"   ")));cat(paste(substring(string.graph,1,temp.locate-1),"-->",substring(string.graph,temp.locate+4,label.locate-1),sep=""))
                          } else {
                             cat(paste(replicate(space.index,"   ")));cat(paste(substring(string.graph,1,temp.locate-1),"<->",substring(string.graph,temp.locate+4,label.locate-1),sep=""))
                          }
                          temp1.index <- regexpr(",",string.graph,fixed=T)
                          temp2.index <- regexpr("];",string.graph,fixed=T)
                          if (temp1.index>1) {
                             cat(paste("(",substring(string.graph,label.locate+8,temp1.index-1),")",sep=""))
                          } else {
                             cat(paste("(",substring(string.graph,label.locate+8,temp2.index-1),")",sep=""))
                          }
                          cat("\n")
                       } else {
                          if ( isTRUE(all.equal(subgraph.locate,1)) ) {
                             return.index[1] <- 1
                             return.index[3] <- 1
                             temp.index <- regexpr("cluster",string.graph,fixed=T)
                             temp1.index <- regexpr(" {",string.graph,fixed=T)
                             cat(paste(replicate(space.index+1,"---")));cat(paste(substring(string.graph,temp.index,temp1.index-1),sep=""))
                             cat("\n")
                          } else {
                             if ( string.graph=='}' ) {
                                return.index[2] <- 1
                                return.index[3] <- -1
                             } else {
                                if ( start.label.locate!=1 ) {
                                   cat(paste(replicate(space.index,"   "),sep=""))
                                   cat(paste(string.graph,sep=""))
                                   for ( ii in iline.index:length(file.graph) ) {
                                       icluster.temp <- regexpr("cluster",file.graph[ii],fixed=T)
                                       if ( isTRUE( all.equal(icluster.temp,1) ) ) {
                                          ileft.temp <- regexpr(" [",file.graph[ii],fixed=T)
                                          ieq.temp <- regexpr(" -- ",file.graph[ii],fixed=T)
                                          if ( substring(file.graph[ii],ieq.temp+4,ileft.temp-1) == string.graph ) {
                                             ilabel.temp <- regexpr(",",file.graph[ii],fixed=T)
#                                            cat(paste("(",substring(file.graph[ii],ileft.temp+8,ilabel.temp-1),")",sep=""))
                                          }
                                       }
                                   }
                                   cat("\n")
                                }
                             }
                          }
                       }
                       return(return.index)
        }
        
        Igraph.list <- function(string.graph) {
                       left.brace.locate <- regexpr("[",string.graph,fixed=T)[1]
                       arrow.locate <- regexpr("->",string.graph,fixed=T)[1]
                       label.locate <- regexpr("label",string.graph,fixed=T)[1]
                       dir.locate <- regexpr("dir",string.graph,fixed=T)[1]
                       if ( isTRUE( all.equal(dir.locate,-1) ) ) {
                          right.brace.locate <- regexpr("]",string.graph,fixed=T)[1]
                          cat(paste(substring(string.graph,1,arrow.locate-1)," ----> ", substring(string.graph,arrow.locate+2,left.brace.locate-1), " (", substring(string.graph,label.locate+6,right.brace.locate-1), ")" ,sep=""),file="Iscript.txt",append=T)
                       } else {
                          cat(paste(substring(string.graph,1,arrow.locate-1)," <----> ", substring(string.graph,arrow.locate+2,left.brace.locate-1), " (", substring(string.graph,label.locate+6,dir.locate-2), ")" ,sep=""),file="Iscript.txt",append=T)
                       }
                       cat("\n",file="Iscript.txt",append=T)
                       return()
        }   
                          
        ##########################################################################################

        # copy graph0.txt, graph00.txt and multistagefca.rb to current directory
        file.copy( system.file("extdata/graph0.txt",package="nFCA"), "graph0.txt", overwrite=T )
        file.copy( system.file("extdata/graph00.txt",package="nFCA"), "graph00.txt", overwrite=T )
        file.copy( system.file("exec/multistagefca.rb",package="nFCA"), "multistagefca.rb", overwrite=T )
        
        # parameters checking
        if ( !(data.type %in% c(0,1,2)) ) {
           stop( "the specified 'type' is not valid" )
        }
        
        if ( !(method.type %in% c("hist","CI")) ) {
           stop( "the specified 'method' is not valid" )
        }
        
        if ( !(hist.manual %in% c(0,1)) ) {
           stop( "the specified 'choice' is not valid" )
        }
        
        if ( !missing(n.CI) ) {
           if ( n.CI<0 ) {
              stop( "the specified 'n' is not valid" )
           } else {
              n.CI <- ceiling(n.CI)
           }
        }
        
        if ( !missing(alpha.CI) ) {
           if ( (alpha.CI<=0) | (alpha.CI>=1) ) {
              stop( "the specified 'alpha' is not valid" )
           }
        }
        
        if ( !missing(data) ) {
           if ( is.data.frame(data) ) {
               data.ri <- data
           } else {
               stop("input data does not exist or is not a data frame")
           }
        }
        
        if ( isTRUE( all.equal(hist.manual,0) ) ) {
           hist.manual <- 1
        } else {
           hist.manual <- 0
        }
        
        if ( method.type!="hist" ) {
            hist.manual <- 0
        }

        # nFCA is now working
        if ( isTRUE( all.equal(data.type,0) ) ) {
           nFCAmain(data=data.ri, manual=hist.manual, method=method.type, n=n.CI, alpha=alpha.CI)
        }
        
        if ( isTRUE( all.equal(data.type,1) ) ) {
           nFCAmain(data=data.ri, manual=hist.manual)
        }
        
        if ( isTRUE( all.equal(data.type,2) ) ) {
           data.p <- as.data.frame(cor(data.ri))
           nFCAmain(data=data.p, manual=hist.manual, method=method.type, n=nrow(data.ri), alpha=alpha.CI)
        }
        
        # remove graph0.txt, graph00.txt, multistagefca.rb and other interstate files from current directory
        unlink("graph*.txt")
        unlink("example.*")
        unlink("multistagefca.rb")
        
        file.rename("Hgraph.txt", "Hgraph.dot")
        file.rename("Igraph.txt", "Igraph.dot")
        
        ##################################################################################################
        cat("\n")
        cat("Hgraph.dot and Igraph.dot are generated in the current R working directory.\n")
        cat("You can go outside of R and use Graphviz to visualize high quality H- and I-graphs.\n")
        cat("\n")
        cat("If you can not visualize the high quality H and I graphs using Graphviz,\n")
        cat("here is the digitalized presentation of these graphs.\n") 
        cat("\n")
        cat("Explanation of the digitized H-graph:\n")
        cat("\n")
        cat("Each '---' is one depth further into the hierarchical center from the\n")
        cat("outer most boundary, e.g.\n")
        cat("--- 1-depth into the center from the outside,\n")
        cat("--- --- 2-depth into the center,\n")
        cat("if not starting with '---', it is in the most outside layer.\n")
        cat("\n")
        cat("Style at the same depth:\n")
        cat("\n")
        cat("--- Cluster#\n")
        cat("    A member, its relationship to other members (relation strength value)\n")
        cat("    A member of the same cluster\n")
        cat("\n")
        cat("Actual Presentation:\n")
        cat("\n")
        cat("Digital presentation of clustering results from H-graph:\n")
        cat("\n")
        file.copy("Hgraph.dot", "Hprint.dot", overwrite=T)
        file.graph <- readLines("Hprint.dot")
        file.graph <- file.graph[3:(length(file.graph)-1)]
        iright.count <- 0
        for ( i in 1:length(file.graph) ) {
            if (file.graph[i]=='}') {
               iright.count <- iright.count + 1
            }
        }
        ileft.index <- 0
        iright.index <- 0
        iline.index <- 1
        space.index <- 0
        while ( iright.index < iright.count ) {
              i.index <- cluster.cat(file.graph[iline.index],space.index,iline.index,file.graph)
              ileft.index <- ileft.index + i.index[1]
              iright.index <- iright.index + i.index[2]
              space.index <- space.index + i.index[3]
              iline.index <- iline.index + 1
        }
        unlink("Hprint.dot")
        
        cat("\n")
        cat("Digital inherent structure results from I-graph:\n")
        cat("\n")
        file.copy("Igraph.dot", "Iprint.dot", overwrite=T)
        file.graph <- readLines("Iprint.dot")
        file.graph <- file.graph[3:(length(file.graph)-1)]
        for ( i in 1:length(file.graph) ) {
             Igraph.list(file.graph[i])
        }
        Iscript.graph <- readLines("Iscript.txt")
        inumber <- rep(0,length(Iscript.graph))
        for ( i in 1:length(Iscript.graph) ) {
            ileft_locate <- regexpr("(",Iscript.graph[i],fixed=T)
            iright_locate <- regexpr(")",Iscript.graph[i],fixed=T)
            inumber[i] <- substring(Iscript.graph[i],ileft_locate+1,iright_locate-1)
        }
        inumber <- as.numeric(inumber)
        isort <- (sort(inumber,decreasing=T,index.return=T))$ix
        for ( i in 1:length(Iscript.graph) ) {
            cat(Iscript.graph[isort[i]])
            cat("\n")
        }
        cat("\n")
        unlink("Iprint.dot")
        unlink("Iscript.txt")
        ##################################################################################################
}