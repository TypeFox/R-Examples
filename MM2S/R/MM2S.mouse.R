# Deena M.A. Gendoo
# October 29, 2014
# Code to conduct MM2S predictions on mouse samples

# DISCLAIMER:
# MM2S package (and its code components) is provided "AS-IS" and without any warranty of any kind. 
# In no event shall the University Health Network (UHN) or the authors be liable for any consequential damage of any kind, 
# or any damages resulting from the use of MM2S.

#################################################################################
#################################################################################

if(getRversion() >= "2.15.1")  utils::globalVariables(c("MouseGMT", "genesetHuman","en_ES_Rank_Matrix","MB_SampleInfo"))

MM2S.mouse<-function(InputMatrix,xls_output,parallelize)
{
  set.seed(12345)
  options(warn=-1)
  ###################################
  ## Parameter Checks
  ###################################
  ## Check & (re)format Input Matrix as needed
  #mouseData <- read.csv(file=FILES$jobs_data$tmp_name, sep=SEP, stringsAsFactors=FALSE, header=FALSE)
  mouseData<-InputMatrix
  mdm <- NULL
  if (is.na(as.numeric(mouseData[1,1])))
  {
    mdm <- mouseData[-1,-1, drop=FALSE]
    colnames(mdm) <- mouseData[1,][-1]
    rownames(mdm) <- mouseData[,1][-1]
    mouseData <- mdm
  }
  if (is.na(as.numeric(mouseData[1,ncol(mouseData)])))
  {
    mdm <- mouseData[-1,-1, drop=FALSE]
    colnames(mdm) <- mouseData[1,][-ncol(mouseData)]
    rownames(mdm) <- mouseData[,1][-1]
    mouseData <- mdm
  } 
  # Check as.numeric
  ExpressionMatrixMouse <- as.matrix(mouseData)
  ExpressionMatrixMouse <- apply(ExpressionMatrixMouse, c(1,2), as.numeric)
  rownames(ExpressionMatrixMouse) <- rownames(mouseData)
  
  ## Check boolean CSV
  if(!is.logical(xls_output))
  {
    message("TRUE or FALSE needed for XLS output")
    stop()
  }
  
  ###################################
  ## Perform ssGSEA & get Rank Matrix
  ###################################
  availcore=1
  if(!is.numeric(parallelize))
  {
    message("Number of Cores needed")
    stop()
  } else{availcore=parallelize}
  
  ## Call the data
  MouseData<-ExpressionMatrixMouse
  #Estimate ssGSEA scores using GSVA function
  # GSVA -> S4 method for signature 'matrix,list,character'
  set.seed(12345)
  MouseGSVA<-gsva(MouseData, MouseGMT$genesets,method="ssgsea", ssgsea.norm=FALSE, min.sz=20,max.sz=100, parallel.sz=availcore,verbose=FALSE)
  genesetMouse<-rownames(MouseGSVA)
  
  #Find common geneset between the Mouse and Human, 
  # especially because GSVA filtering reduced the sets
  # Now the Rank matrix will be reduced based on the commonSet
  commonSet<-intersect(genesetHuman,genesetMouse)
  message("There are ",length(commonSet)," common genesets between Human MB and the Test Data.")
  
  MouseGSVA<-MouseGSVA[commonSet,, drop=FALSE]
  MouseGSVA<-t(MouseGSVA)
  
  ## FEATURE SELECTION FOR GENESETS TO USE IN THE RANKING
  #First get the subset of the data pertaining to just one group
  GenesetStatNormal<-GenesetStatNormal[commonSet]
  GenesetStatGroup3<-GenesetStatGroup3[commonSet]
  GenesetStatGroup4<-GenesetStatGroup4[commonSet]
  GenesetStatWNT<-GenesetStatWNT[commonSet]
  GenesetStatSHH<-GenesetStatSHH[commonSet] 
  
  #Pick the features/genesets that best differentiate each subtype, get their names
  #Pick from the FROZEN Matrix of Wilcoxon P-values for all genesets, across subtypes, for Human MB
  geneset<-commonSet
  FeatureSelection<-c(names(sort(GenesetStatSHH,decreasing=FALSE))[1:24],
                      names(sort(GenesetStatNormal,decreasing=FALSE))[1:24],
                      names(sort(GenesetStatGroup4,decreasing=FALSE))[1:24],
                      names(sort(GenesetStatGroup3,decreasing=FALSE))[1:24],
                      names(sort(GenesetStatWNT,decreasing=FALSE))[1:24])   
  #Remove the redundant genesets 
  NorthcottFeatures<-unique(FeatureSelection)
  message("Of these, ", length(NorthcottFeatures)," feature-selected genesets are being used for classification")
 
  MouseGSVA<-MouseGSVA[,NorthcottFeatures, drop=FALSE]
  MouseGSVA<-t(MouseGSVA)
  genesetMouse<-rownames(MouseGSVA)
  
  #Get Ranking Matrix for Human Data
  Matrix_RANK_Human<-data.frame(NorthcottFeatures)
  for(sample in 1:ncol(Frozen_ES_Rank_Matrix))
  {
    TempRankHuman<-Frozen_ES_Rank_Matrix[,sample]
    TempRankHuman<-TempRankHuman[which(TempRankHuman %in% NorthcottFeatures)]
    Matrix_RANK_Human[,(colnames(Frozen_ES_Rank_Matrix)[sample])]<-match(factor(NorthcottFeatures),factor(TempRankHuman))
  }
  rownames(Matrix_RANK_Human)<-NorthcottFeatures
  Human_GSVA_Matrix<-Matrix_RANK_Human[,-1]
  
  #Get Ranking Matrix for Mouse Data
  Matrix_RANK_Mouse<-data.frame(genesetMouse)
  for(sample in 1:ncol(MouseGSVA))
  {
    TempRankMouse<-sort(MouseGSVA[,sample],decreasing=TRUE)
    Matrix_RANK_Mouse[,(colnames(MouseGSVA)[sample])]<-match(Matrix_RANK_Mouse$geneset,names(TempRankMouse))
  }
  rownames(Matrix_RANK_Mouse)<-Matrix_RANK_Mouse$geneset #Rename the matrix by sample names
  Mouse_GSVA_Matrix<-Matrix_RANK_Mouse[,-1, drop=FALSE] #Remove the geneset column
  
  Human_GSVA_Matrix80<-Human_GSVA_Matrix
  Mouse_GSVA_Matrix80<-Mouse_GSVA_Matrix
  
  ###################################  
  ## PART B - Generate Subtype Predictions
  ###################################
  
  geneset<-rownames(Human_GSVA_Matrix80) #List of all the genesets that were used in ssGSEA
  
  # Transpose the matrix, so that Samples are in rows and genesets in columns
  Northcott<-t(Human_GSVA_Matrix80) ##CHANGE DEPENDING ON WHETHER YOU WANT THE FULL MATRIX OR THE SUBSET!!
  Mouse<-t(Mouse_GSVA_Matrix80)
  
  # Add the group prediction to the Northcott data (WNT, SHH, Group3, Group4)
  # Groups: WNT, SHH, Group3, Group4, NA
  Northcott<-as.data.frame(Northcott)
  Northcott[,"Group"]<- MB_SampleInfo$subtype[match((rownames(Northcott)),MB_SampleInfo$Sample_ID)]
  
  # Add a column to the dataset, containing the Group labels per sample (Only Normal vs Tumour in this case)
  Mouse<-as.data.frame(Mouse)
  Mouse[,"Group"]<- "MouseSamples"
  Mouse$Group<-as.factor(Mouse$Group)
  
  TrainSet<-Northcott
  TestSet<-Mouse
  
  TrainSet<-TrainSet[,-ncol(TrainSet)]
  ## Generate the predictions
  set.seed(12345)
  
  TestKKNN<-kknn(formula = Northcott$Group ~ ., TrainSet, TestSet, na.action = na.omit(),k = 5, distance = 1, kernel = "rectangular", scale=TRUE)
  
  # ConfusionMatrixTest<-table(Mouse$Group,TestKKNN$fitted.values)
  MM2S_Prediction<-as.character(TestKKNN$fitted.values)
  
  RESULTS<-(cbind(rownames(Mouse),MM2S_Prediction,TestKKNN$prob*100,TestKKNN$CL))
  listOfCols<-c("SampleName","MM2S_Prediction","Gr3_Confidence","Gr4_Confidence","Normal_Confidence","SHH_Confidence","WNT_Confidence","Neighbor1","Neighbor2","Neighbor3","Neighbor4","Neighbor5")
  
  ## PRODUCE OUTPUT WITH THE MM2S PREDICTIONS PER SAMPLE
  colnames(RESULTS) <- listOfCols
  
  message("\n")
  message("OUTPUT OF MM2S:","\n")
  
  print.table(RESULTS) 
  
  if(xls_output==TRUE)
  {
    write.table(RESULTS,file="MM2S_Predictions.xls",sep = "\t",col.names=listOfCols,row.names=FALSE)
  }
  
  FINAL<-TestKKNN$prob*100
  colnames(FINAL)<-c("Group3","Group4","Normal","SHH","WNT")
  rownames(FINAL)<-rownames(Mouse)
  
  return(list(RankMatrixTesting=t(Mouse_GSVA_Matrix80), RankMatrixTraining=t(Human_GSVA_Matrix80), Predictions=FINAL,MM2S_Subtype=RESULTS[,1:2]))
  
}   # End of function


