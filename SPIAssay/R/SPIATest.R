############################################################################
#
# Compute SPIA test
#  Input: 1. x: a matrix with a column for each cell line and a row for each SNP
#            in the SPIA format (use toSPIAData before)
#         2. row.names: specify if the fisrt column contains the name of the SNPs
#         3. test.prob: specify if probabilistic test has to be performed
#         4. test.param: specify the parameters of the test
#            - test.param$Pmm: SNP probability of mismatch in a matching population (dafault 0.1)
#            - test.param$nsigma: area limit for Pmm
#            - test.param$Pmm_nonM: SNP probability of mismach in a non matching population (dafault 0.6)
#            - test.param$nsigma_nonM: area limit for Pmm_nonM
#            - test.param$PercValidCall: percentage of valid call to consider the test valid
#         5. verbose: print verbose information
#  Output: a matric with a line for each cell line and with columns with the
#          informationss about distances 
#
#############################################################################

SPIATest<-function (x, row.names = TRUE, test.prob = TRUE,  
                    test.param = list(Pmm = 0.1, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, PercValidCall=0.9), verbose = FALSE) 
{  
  
  #Print welcome message  
  message("SPIA: start analysis:")

  
  #If row.names is set, the first column is removed from x
  if (row.names){
    rowNames <- x[,1]
    x <- x[,c(2:dim(x)[2])]
  }
  
  ##Nr is the number of SNP
  Nr <- nrow(x);  
  ##Nc is the number of cell lines
  Nc <- ncol(x);  
  ##cuples gives the number of pairs: if there is 4 cell lines  
  ##       a total of (4 * 3 )/2 = 6 combinations are possible
  couples<-(Nc * (Nc - 1)/2)
  
  
  if (verbose){
    #Print the statistics about the data
    message(paste("Number of samples:",Nc));
    message(paste("Number of pairs:",couples));
    message(paste("Number of SNPs:",Nr));
  }
  
  if (test.prob){
    #Create the that store the results of the distance computing
    result<-array(,c(couples,13));  
    colnames(result)<-c("Cell_1","Cell_2","Distance","SPIA_Score","SNP_available","Total_SNP","One_SNP_NA","Bot_SNP_NA","Diff_AvsB_or_BvsA","Diff_AorBvsAB_or_vic","DiffABvsAorB","counterBothHomoz","counterBothHeter");
    }  else  {
    #Create the that store the results of the distance computing with the result of the statistical test
    result<-array(,c(couples,12));  
    colnames(result)<-c("Cell_1","Cell_2","Distance","SNP_available","Total_SNP","One_SNP_NA","Bot_SNP_NA","Diff_AvsB_or_BvsA","Diff_AorBvsAB_or_vic","DiffABvsAorB","counterBothHomoz","counterBothHeter");
  }
  
  #Verify if cell line names are available
  if (!is.null(colnames(x))){
    CLnames <- colnames(x)
  } else if (! is.null(names(x)) ) {
        CLnames <- names(x)
    } else {
        CLnames <- c(1:dim(x)[2])
      }      
  #Print dots when computing (each 20 pairs)
  if (couples < 20){
    steps <- 1
  } else {
    steps <- as.integer(couples/20) + 1
  }
  #The counter is used to trace the current pair of cell line processed by the for cycle
  counter<-0;
  for(i in (1:(Nc-1))){  
    for(j in ((i+1):Nc)){                 
        #Compute distance and other information for the cell i vs cell j
        dist<-ComputingSPIADistance(x[,i],x[,j]);
        
        #Save the result in result and increment counter        
        result[(counter+1),1] <- CLnames[i]
        result[(counter+1),2] <- CLnames[j]
        
        #If the probabilistic test is not eables
        if (!test.prob){
          #save the current value of the distance
          result[(counter+1),c(3:12)]<-dist  
        }  else  {
          #save the current value of the distance and verify the probabilist test
          result[(counter+1),3]<-dist[1]
          #call the probabilistic test
          limits <- getSPIALimits(dist[2],test.param$Pmm,test.param$nsigma,test.param$Pmm_nonM,test.param$nsigma_nonM)
          #verify if there are problems with the limits defined by the binomial distributions
          if (limits$err)
          {
            message("SPIA error: the parameters of the probabilistic probabilistic test define ")
            message("  a negative 'similar' region.")
            message(paste("  The limit for the mismatch in a matching population is ",limits$liminf,sep=""))
            message(paste("  The limit for the mismatch in a mismatching population is ",limits$limsup,sep=""))
            return(-1)
          }
          
          #save the result of the probabilistic test
          if ( dist[2] / dist[3] < test.param$PercValidCall){
            result[(counter+1),4] <-"Not Valid"
          } else 
          if ( dist[1] < limits$liminf ){
            result[(counter+1),4] <-"Similar"
          } else 
          if ( dist[1] > limits$limsup ){
            result[(counter+1),4] <-"Different"
          } else {
            result[(counter+1),4] <-"Uncertain"
          }                    
          result[(counter+1),c(5:13)]<-dist[c(2:10)]
        }
        
        #print one dot each n steps
        if ((counter %% steps) == 0) { cat(".") }
        counter<-counter+1
      }
    }
  message("SPIA: analysis completed.")  
  return(list(SPIAresult = result, parameters = test.param, input.param = list(N_samples = Nc, N_SNPs = Nr, testDone = test.prob)))
}
