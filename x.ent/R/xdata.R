xdata <- function(e=NULL)
{
  tryCatch(
  {
    #create a data frame
    conf = fromJSON(paste(.libPaths()[1], "x.ent/www/config/ini.json", sep='/'))
    path = conf$result$file;
    lst_f <- xfile(sep=":")
    #get all tags in the file config
    lst_tag <- xentity();
    lst_tag <-add_unique(lst_tag,xrelation())
    dta <- data.frame()
    reg_ent = ":\\$:"
    reg_rel = ":\\$\\$:"
    dta <- data.frame(file=lst_f)
    #read file output.txt
    text <- readLines(path)
    if(is.null(e))
    {
      for(i in 1:length(lst_tag))
      {
        dta[,lst_tag[i]] <- "N/A" 
      }
      for(i in 1:length(text))
      {
        #get name of file
        f <- unlist(strsplit(text[i],":"))[1]
        #find entity
        if(grepl(pattern=reg_ent,x=text[i]))#entity
        {
          data_ele <- ""
          eles <- unlist(strsplit(text[i],":"))
          if(length(eles) == 4)
          {
            data_ele <- eles[4]
          }
          else if(length(eles) >= 5)
          {
            data_ele <- eles[4]
            for(j in 5:length(eles))
            {
              data_ele <- paste(data_ele,eles[j],sep="; ")
            } 
          }
          dta[dta$file == f,eles[2]] <- data_ele
        }
        #relation
        if(grepl(pattern=reg_rel,x=text[i]))
        {
          eles <- unlist(strsplit(text[i],reg_rel))
          col = sub(pattern = paste(f,":",sep=""), replacement = "",x = eles[1])
          col = gsub(":$", "",col, perl=TRUE)#delete ":" at the end of sentence
          if(col %in% names(dta))
          {
            if(dta[dta$file == f,col] == "N/A")
            {
              dta[dta$file == f,col] <- eles[2] 
            }
            else
            {
              dta[dta$file == f,col] <- paste(dta[dta$file == f,col],eles[2],sep =";")  
            }
          }
        }
      }
    }
    else
    {
      for(i in 1:length(e))
      {
        dta[,e[i]] <- "N/A" 
      }
      for(i in 1:length(text))
      {
        #get name of file
        f <- unlist(strsplit(text[i],":"))[1]
        for(j in 1:length(e))
        {
          #entity
          reg = paste(f,":",e[j],reg_ent,sep="")
          #find all
          if(grepl(pattern=reg,x=text[i]))#entity
          {
            data_ele <- ""
            eles <- unlist(strsplit(text[i],":"))
            if(length(eles) == 4)
            {
              data_ele <- eles[4]
            }
            else if(length(eles) >= 5)
            {
              data_ele <- eles[4]
              for(j in 5:length(eles))
              {
                data_ele <- paste(data_ele,eles[j],sep="; ")
              } 
            }
            dta[dta$file == f,eles[2]] <- data_ele
          }
          #relation
          reg = paste(f,":",e[j],reg_rel,sep="")
          if(grepl(pattern=reg,x=text[i]))#relation
          {
            result = gsub(reg, "",text[i], perl=TRUE)
            if(e[j] %in% names(dta))
            {
              if(dta[dta$file == f,e[j]] == "N/A")
              {
                dta[dta$file == f,e[j]] <- result 
              }
              else
              {
                dta[dta$file == f,e[j]] <- paste(dta[dta$file == f,e[j]],result,sep =";")  
              }
            }
          }
        }
      }
    }
    return(dta)
  },
  error=function(cond) {
    message("Parameters are incorrect or there are problems in paths, please check your parameters!")
  },
  warning=function(cond) {
    message("Parameters are incorrect or there are problems in paths, please check your parameters!")
  },
  finally={
    rm(list=ls())
  })
}