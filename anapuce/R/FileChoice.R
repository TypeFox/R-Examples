FileChoice <-
function(pattern=NULL){
    fic1     <- file.choose()
    datapath <- dirname(fic1)
    setwd(datapath)
    files <- list.files(pattern=pattern)

 cat("\n")
        files.ana <- NULL
        for( i in 1:length(files)) {
        cat("\t[",i,"]\t",files[i],"\n", sep="")
        repeat{
        answer <- readline("\tValidate (y/n)? ")
        if(answer=="n") break
        if(answer=="y")
        {
        files.ana <- c(files.ana,files[i])
        break           
        }
              }
     } 
     files <- files.ana
     if(length(files)>1){
        cat("\n\n*** Change slides order \n\n")


                # Display files names in order
  cat("\n")
  for( i in 1:length(files)){
    cat("\t[",i,"]\t",files[i],"\n", sep="")
  }


    repeat{
      ord <- readline("\n\tChange order [ ex : 1,3,2,..] or type [enter] if no change : ")
      


      if(ord=="")
        ord <- "-1"
         else{
        ord <- as.numeric( unlist(strsplit(ord,",")) )
        ord <- ord[!is.na(ord)]
      }
      if( ord[1] == -1)
        break


      if(length(ord) != length(files) )
        cat("Error! a file has been omitted.")
      else{
        if(length(ord)!=length(unique(ord)) ){
          cat("Error! a file has been set more than once.")
        } else{
          files <- files[ord]
          cat("\n")
          for( i in 1:length(files)){
            cat("\t[",i,"] ",files[i],"\n", sep="")
          }
          repeat{
            answer <- readline("\tValidate (y/n)? ")
            if(answer=="n" | answer=="y")
              break
          }
          if(answer=="y")
            break
        }
      }
    }
    
  }
  cat("\n")
  invisible(files)
  # (c) 2007 Institut National de la Recherche Agronomique

  }

