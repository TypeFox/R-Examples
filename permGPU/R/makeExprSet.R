makeExprSet <-
function(exprdat,phenodat,anno="custom")
  {
    if(ncol(exprdat)==nrow(phenodat))
      {
        colnames(exprdat)=rownames(phenodat)
        datobj=new("ExpressionSet",exprs=as.matrix(exprdat))
        pData(datobj)=phenodat
        annotation(datobj)=anno
      }
    else
      {
        print("Number of samples of the expression and phenotype data sets do not match")
        datobj=NULL
      }
    return(datobj)
  }

