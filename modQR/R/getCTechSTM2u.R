getCTechSTM2u <- function(){
#getCTechSTM2u <- function(), output: CTechST
#this function sets the options for computing all the directional (regression) quantiles by means of the main program (= compContourM2u)
CTechST <- list()
#universal activity indicators
CTechST$ReportI       <- 0           #1/0 ... if some information (such as the progress of computation) is displayed on the screen (1) or not (0)
                                     #- (1) may be useful for large problems with M > 2
CTechST$OutSaveI      <- 0           #1/0 ... if the detailed output is stored in file(s) on the disk (1) or not (0)
                                     #- (0) might be useful for novice users, for processing large/high-dimensional data sets or for simulation studies
CTechST$D2SpecI       <- 1           #1/0 ... if, for M = 2, the cones are visited counter-clockwise (1) or as in any other dimension (0)
                                     #- (1) leads to a more precise computation than (0)
CTechST$BriefOutputI  <- 1           #1/0 ... if the brief (1) or verbose (0) version of the output is provided by the program
                                     #- (1) leads to the (file) output rows of length (1+1+M+P*M+M):
                                     #        c(ConeID, Nu, UVec, vec(ACOMat), MuBRow)
                                     #- (0) leads to the (file) output rows of length (1+1+M+P*M+M+P+P):
                                     #        c(ConeID, Nu, UVec, vec(ACOMat), MuBRow, MuR0Row, IZ)
                                     #  where
                                     #  ConeID  ..  the number/order of the cone investigated
                                     #  Nu      ... the number of negative residuals
                                     #  UVec    ... an L2 or max normalized vector of the cone investigated (depending on CTechST$D2SpecI)
                                     #              (= BVec' after L2 normalization)
                                     #  ACOMat  ... the matrix for computing AVec, AVec = ACOMat%*%UVec (if ||UVec|| = 1)
                                     #  MuBRow  ... the row multiplier vector associated with the constraint BVec = UVec
                                     #  MuR0Row ... the row multiplier vector associated with zero residuals
                                     #  IZ      ... the indices of observations with the zero residuals (i.e. on the tau- and UVec-quantile hyperplane)

#relevant only if (M > 2) or !CTechST$D2SpecI
CTechST$CubRegWiseI   <- 1            #1/0 ... if the directional space is divided into orthants investigated separately (1) or not (0)
                                      #- (1) splits the problem to 2^M smaller ones but generates some artificial cones with a facet in the orthant borders
CTechST$ArchAllFI     <- 1            #1/0 ... if all the past cone facet identifiers (1) or only those from the last few layers (0) are stored
                                      #- (1) makes the computation more likely to terminate successfully than (0)
                                      #- (0) is faster and less memory demanding than (1)
                                      #- if M > 3 then CTechST$ArchAllFI == 1 is set by the program
CTechST$SkipRedI      <- 0            #1/0 ... if the information should be skipped (1) or stored (0) also from the cones with all non-artificial facets
                                      # already known (such cones are redundant/irrelevant for the quantile contour computation by this method)
                                      #- (1) makes the output smaller but maybe also slightly less reliable

#user-defined strings (not checked for errors by the main program compContourM2u!)
CTechST$OutFilePrefS <- 'DQOutputM2_' # ... the prefix of possible output files
CTechST$getCharST    <- getCharSTM2u  # ... the user-defined output function
return(CTechST)
}