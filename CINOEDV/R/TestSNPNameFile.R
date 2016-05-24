TestSNPNameFile <-
function(RowNum,SNPNameFileName=NA){
  # Check the parameter 'SNPNameFileName'.
  #
  # If there are real SNP names which will be used for constructing graphs and further
  # analysis, the name of file that saves real SNP names should be provided.
  #
  # input
  #     RowNum: the number of SNPs considered.
  #     SNPNameFileName: the name of file that saves real SNP names.
  #         The format of file is (.mat)
  #         It has only one variable, i.e., Name.
  #         Name: Row -> 1, Column -> SNP Name, and the length of Name is RowNum.
  #
  #     If not exist such file, SNPNameFileName -> NA
  #
  # output
  #     SNPNames
  #
  # For example, test_Name.mat
  #
  # Junliang Shang
  # 3.29/2014
  
  #library(R.matlab)
  
  if (is.na(SNPNameFileName)){
    SNPNames <- NA
    cat("    The SNP Name File :", SNPNameFileName,"\n\n")
  }else
  {
    num=nchar(SNPNameFileName)
    if(length(grep(substr(SNPNameFileName,num-3,num),".mat"))==1){
      if (file.exists(SNPNameFileName)){
        data <- readMat.default(SNPNameFileName)
        if ("Name" %in% names(data)){
          SNPNames1 <- data$Name
          if (length(SNPNames1)==RowNum){
            SNPNames <- sapply(SNPNames1,function(x) return(x[[1]]))
            # names(SNPNames) <- seq(1,RowNum)
            cat("    The SNP Name File :", SNPNameFileName,"\n\n")
          }else
          {
            stop("    The Number of SNP names is not equal to the number of SNPs !\n")
          }          
        }else
        {
          stop("    No 'Name' variable in the file !\n")
        }
      }else
      {
        stop("   ",SNPNameFileName," is not exist !\n")
      }
    }else
    {
      stop("   Error! File format should be (.mat).\n")
    }
  }
  
  # return
  list(SNPNames=SNPNames)
}
