"df2kasc" <- function(df, index, x)
  {
      ## Verifications
      if (!inherits(df,"data.frame")) stop("non convenient data type")
      if ((!inherits(x,"kasc"))&(!inherits(x,"mapattr")))
          stop("non convenient data type")

      ## prepare the data
      class(x)<-"data.frame"
      N<-attr(x, "nrow")*attr(x, "ncol")
      indw<-c(1:N) ## a vector of the indices of the rows of the output kasc
      n1<-nrow(df)

      ## The missing values to be inserted ("compl"ementary)
      compl<-as.data.frame(matrix(NA, nrow=N-n1, ncol=ncol(df)))
      names(compl) <- names(df)

      ## we bind these NA to the DF
      output<-rbind.data.frame(df, compl)

      ## The indices of the pixels with NA in the output kasc
      indcompl<-indw[is.na(match(indw, index))]

      ## concatenate the vectors of indices of the pixels with values (index)
      ## and the vectors of indices of the pixels with NAs (indcompl)
      indtot<-c(index, indcompl)

      ## sort the table according to this resulting vector
      output<-output[sort(indtot, index.return=TRUE)$ix,]

      ## Output
      class(output)<-c("kasc","data.frame")
      attr(output, "nrow")<-attr(x, "nrow")
      attr(output, "ncol")<-attr(x, "ncol")
      attr(output, "xll")<-attr(x, "xll")
      attr(output, "yll")<-attr(x, "yll")
      attr(output, "cellsize")<-attr(x, "cellsize")

      return(output)
  }

