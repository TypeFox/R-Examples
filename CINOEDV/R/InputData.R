InputData <-
function(FileName){
  # Input SNP data from a file with (.mat) format
  #
  # input
  #   file name with (.mat) format that saves SNP data
  #   the file has two variables, i.e., pts and class.
  #   pts: Row -> Sample, Column -> SNP
  #      1 -> AA
  #      2 -> Aa
  #      3 -> aa
  #   class: Row -> 1, Column -> class label
  #      1 -> case
  #      2 -> control
  #
  # output
  #   pts
  #   class
  #
  # for example
  #   Data <- InputData("test.mat")
  #
  # Junliang Shang
  # 3.25/2014
  
  #library(R.matlab)
  
  cat(" The input file is :", FileName,"\n")
  num=nchar(FileName)
  
  # Check file format (.mat).
  if(length(grep(substr(FileName,num-3,num),".mat"))==1){
    # Check whether the file is exist.
    if (file.exists(FileName)){
      data <- readMat.default(FileName)
      # Check the variables in the file
      if (("pts" %in% names(data)) && ("class" %in% names(data))){
        pts <- data$pts
        class <- data$class
        RowColNum <- dim(pts)
        # Check values in pts and class
        if ((max(pts)==3)&&(min(pts)==1)&&(max(class)==2)&&(min(class)==1)&&
              (nrow(pts)==ncol(class))&&(nrow(class)==1)){
          colnames(class) <- paste("Sample",1:RowColNum[1],sep="")
          rownames(pts) <- paste("Sample",1:RowColNum[1],sep="")
          colnames(pts) <- paste("SNP",1:RowColNum[2],sep="")
          cat("    Samples:",nrow(pts),"\n")
          cat("      Cases:",length(class[class==1]),"\n")
          cat("      Controls:",length(class[class==2]),"\n")
          cat("    SNPs:",ncol(pts),"\n\n")
          list(pts=pts,class=class)
        }
        else
        {
          stop(" Values in pts or class is error !\n")
        } 
      }else
      {
        stop(" No 'pts' or 'class' variables in the file !\n")
      }
    }else
    {
      stop("  ",FileName," is not exist !\n")
    } 
  }else
  {
    stop("  Error! File format should be (.mat).\n")
  } 
}
