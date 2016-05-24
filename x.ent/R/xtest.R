xtest <- function(v1,v2)
{
  if((length(v1) < 1) || (length(v2) < 2))
  {
    print("This vector must be greater than 2 elements")
  }
  else
  {
    #get the result as a data frame
    dt <- xdata()
    rows <- nrow(dt)
    #build a data frame for containing the result
    df <- NULL
    for(i in 1:length(v1))
    {
      for(j in 1:length(v2))
      {
        reg <- paste(v1[i],":",v2[j],sep="")
        df1 <- data.frame(entity = reg,stringsAsFactors=FALSE)
        reg <- paste(reg,":",sep="")
        for(r in 1:rows)
        {
          name_col <- paste("d",r,sep="")
          df1[,name_col] <- 0
          #search each row of data frame
          cols <- grep(reg,dt[r,])
          if(length(cols) > 0)
          {
            df1[1,name_col] <- 1
          }
        }
        if(is.null(df))
        {
          df= df1
        }
        else
        {
          df <- rbind(df,df1) 
        }
      }
    }
    #tranpose the  data frame
    names <- df$entity
    df <- as.data.frame(t(df[,-1]))
    colnames(df) <- names
    n <- ncol(df)
    col_names <- names(df)
    #create a data frame for store every test of relation
    result <- NULL
    if(n >= 2)
    {
     #pair-wised comparaisons
      for(i in 1:(n-1))
      {
        for(j in (i+1):n)
        {
          mat <- NULL
          v1 = df[,i]
          v2 = df[,j]
          mat <- data.matrix(v1)
          mat <- cbind(mat,v2)
          #print labels
          label = paste(col_names[i],"/",col_names[j],sep= "")
          count1 <- 0
          count2 <- 0
          v <- c()
          if(sum(mat[,1]) > 0 && sum(mat[,2]) > 0)
          {
            rs1 <- data.frame(relation = label,stringsAsFactors=FALSE)
            #TEST du KOLMOGOROV
            test <- ks.test(mat[,1], mat[,2])
            rs1[,"KOLMOGOROV"] <- test$p.value  
            #v <- c(v,test$p.value)
            #print("TEST du KOLMOGOROV:",test$p.value,sep=""))
            #TEST de WILCOXON
            test = wilcox.test(mat[,1], mat[,2])
            rs1[,"WILCOXON"] <- test$p.value 
            #v <- c(v,test$p.value)
            #print(paste("TEST de WILCOXON:",test$p.value,sep=""))
            #TEST de STUDENT
            test = t.test(mat[,1], mat[,2])
            rs1[,"STUDENT"] <- test$p.value 
            #v <- c(v,test$p.value)
            #print(paste("TEST de STUDENT:",test$p.value,sep=""))
            rs1[,"GrowthCurves"] <- NA
            tryCatch({
              #TEST de comparaison de courbe de croissance
              test = compareGrowthCurves(mat[,1], as.matrix(mat[,2]))
              rs1[1,"GrowthCurves"] <- test$P.Value
              #v <- c(v,test$P.Value)
              #print(paste("TEST de comparaison de courbe de croissance:",test$P.Value,sep=""))
              #count1 = count1 + 1
            },
            error=function(cond) {
              message("Less than 2 groups to compare!")
            },
            warning=function(cond) {
              message("Less than 2 groups to compare!")
            },
            finally={
              #rm(list=ls())
            })
            ##############
            if(is.null(result))
            {
              result= rs1
            }
            else
            {
              result <- rbind(result,rs1) 
            }
          }
        }
      }
    }
    path = paste(.libPaths()[1],"x.ent/www/statistics.html",sep="/")
    html= print(xtable(result),"html",file=path)
    html = paste("<meta http-equiv=Content-Type content=text/html; charset=utf-8>",html,sep="")
    write(html, file=path)
    browseURL(path)
    return(df)
  }
}

