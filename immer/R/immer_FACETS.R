
#####################################################
# Wrapper to FACDOS (Linacre, 1999)
immer_FACETS <- function(
  title=NULL,
  convergence=NULL,
  totalscore=NULL,
  facets=NULL,
  noncenter=NULL,
  arrange=NULL,
  entered_in_data=NULL,
  models=NULL,
  inter_rater=NULL,
  pt_biserial=NULL,
  faire_score=NULL,
  unexpected=NULL,
  usort=NULL,
  positive=NULL,
  labels=NULL,
  fileinput=NULL,
  data=NULL,
  path.dosbox = NULL,
  path.facets = "",
  model.name=NULL,
  facetsEXE=NULL
){
  
  
  # ------------------------------------
  os_system <- Sys.info()["sysname"]
  user <- Sys.info()['user']
  options(digits.secs=6)									
  time <- base::format(Sys.time(),"_%y%m%d_%H%M%S%OS")
  # ------------------------------------
 if(!is.null(fileinput)){
   # Check if the specification is either NULL or a character vector
  checkInput <- function(x) {is.numeric(x) || is.character(x) || is.null(x)}
  stopifnot(checkInput(title))
  stopifnot(checkInput(totalscore))
  stopifnot(checkInput(facets))
  stopifnot(checkInput(noncenter))
  stopifnot(checkInput(arrange))
  stopifnot(checkInput(entered_in_data))
  stopifnot(checkInput(models))
  stopifnot(checkInput(inter_rater))
  stopifnot(checkInput(pt_biserial))
  stopifnot(checkInput(faire_score))
  stopifnot(checkInput(unexpected))
  stopifnot(checkInput(usort))
  stopifnot(checkInput(positive))
  stopifnot(checkInput(labels))
 }
  # Check if the DosBox is on its place
  if(!is.null(path.dosbox)){
    if(!any(list.files(path.dosbox)=="DOSBox.exe")){
      stop("The Path of DosBox.exe is not correct")
    }  
  }
  if(is.null(facetsEXE)){
    if( length(grep("Facets.exe",path.facets,ignore.case=TRUE))>0){
      stop("The Path of Facets is not correct")
    }
  }
  # ------------------------------------
  # Save the Input-File in the Facets-Folder
  if(is.data.frame(data)){
    data <- rbind("data=",matrix(apply(data,1,paste,collapse=","),ncol=1))
  }
  if(!is.matrix(data) & !is.null(data)){
    data <- try(readLines(data))
  }
  # try to read the specified Inputfile
  if(is.null(fileinput)){
    
  
  models <- rbind("models=",matrix(models,ncol=1),"*")
  labels <- tapply(labels,rep(1:(length(labels)/2),each=2),function(x) c(x,"*"))
  labels <- unlist(labels)
  labels <- rbind("Labels=",matrix(labels,ncol=1))
  # "ISFILE" = paste0("Category",time,".txt"),
  input <- list(
    "title" = title,
    "Output" = paste0("Output",time,".txt"),
    "Scorefile" = paste0("Score",time,".txt"),
    "Residualfile" = paste0("Residual",time,".txt"),
    "Headin lines" = "YES",
    "QM" = "Double",
    "CSV" = "Commas",
    "convergence" = convergence,
    "totalscore" = totalscore,
    "facets" = facets,
    "noncenter" = noncenter,
    "arrange" = arrange,
    "entered in data" = entered_in_data,
    "inter rater" = inter_rater,
    "pt biserial" = pt_biserial,
    "Faire Score" = faire_score,
    "unexpected" = unexpected,
    "usort" = usort,
    "positive" = positive
    )

  
  input.data <- do.call(rbind,input)
  inputfile.1 <- paste0(row.names(input.data)," = ",input.data," ;")
  inputfile <- rbind(matrix(inputfile.1,ncol=1),models,labels,data)
  }else{
    inputfile <-fileinput
  }
  # ------------------------------------------------------------------------------------------------
  # Checks START
  # Check length of specified modelname
  if(!is.null(path.dosbox)){ model.name <- bit8(model.name)}
  filenames <- c()
  woOut <- grep("Output",inputfile,ignore.case=TRUE)
    Out <- inputfile[woOut] ; Out <- grepInput(Out) 
      if(!is.null(path.dosbox)){ 
        inputfile[woOut] <- paste0("Output = ",bit8(Out)[1],"; Name of Output file")
        filenames <- c(filenames,"Outputfile"=bit8(Out)[1]) 
      }
   filenames <- c(filenames,"Outputfile"=Out)
  
  woScor <- grep("Scorefile",inputfile,ignore.case=TRUE)
    Scor <- inputfile[woScor] ; Scor <- grepInput(Scor) 
      if(!is.null(path.dosbox)){ 
        inputfile[woScor] <- paste0("Scorefile = ",bit8(Scor)[1],"; Name of Scorefile file")
        filenames <- c(filenames,"Scorefile"=bit8(Scor)[1]) 
      }
    filenames <- c(filenames,"Scorefile"=Scor)
  
  woResid <- grep("Residualfile",inputfile,ignore.case=TRUE)
    Resid <- inputfile[woResid] ; Resid <- grepInput(Resid) 
      if(!is.null(path.dosbox)){ 
        inputfile[woResid] <- paste0("Residualfile = ",bit8(Resid)[1],"; Name of Residualfile file")
        filenames <- c(filenames,"Residualfile"=bit8(Resid)[1]) 
      }
    filenames <- c(filenames,"Residualfile"=Resid)
    
 # woCate <- grep("ISFILE",inputfile,ignore.case=TRUE)
#    Cate <- inputfile[woCate] ; Resid <- grepInput(Cate) 
#    if(!is.null(path.dosbox)){ 
#      inputfile[woCate] <- paste0("IFILE = ",bit8(Cate)[1],"; Name of Residualfile file")
#      filenames <- c(filenames,"IFILE"=bit8(Cate)[1]) 
#    }
#    filenames <- c(filenames,"Categoryfile"=Cate)  
  
  #woXfil <- grep("Xfile",inputfile,ignore.case=TRUE)
  #if(length(woXfil)>0){
  #  Xfil <- inputfile[woXfil] ; Xfil <- grepInput(Xfil) 
  #  inputfile[woXfil] <- paste0("Xfile = ",bit8(Xfil)[1],"; Name of Xfile file")
  #  # filenames <- c(filenames,"Xfile"=bit8(Xfil)[1])
  #  filenames <- c(filenames,Out)
  #}else{
  #  c(
  #    inputfile[c(1:5),1],
  #    c("Xfile = XFILE.txt; Name of Xfile file"),
  #    inputfile[c(5:nrow(inputfile)),1],
  #    filenames <- c(filenames,"Xfile"="XFILE.txt")
  #  )
  #}
  
  # Checks END
  # ------------------------------------------------------------------------------------------------
  # Check if CSV is specified
  csv <- grep("CSV",inputfile,ignore.case=TRUE)
  if(length(csv)==0){
    c(
      inputfile[c(1:6),1],
      c("CSV=Commas ;"),
      inputfile[c(6:nrow(inputfile)),1]
    )
  }
  # END Check if CSV is specified
  
  # ---------------------------------------------------------------------------------------------------
  
  cat(paste0("writing the inputfile to \n --> ",path.facets," <-- \n"))
  writeLines(inputfile,paste0(path.facets,"\\",model.name))
  # ------------------------------------
  
  # Ausgabe fuer den User
  switch(os_system,
         Windows = 
         {
           cat("I'm a Windows PC. writing batch-job \n")
           # it is very important to quote the first path to the DOSbox (" \" ")
           if(!is.null(path.dosbox)){
             start.dosbox <- c( paste0( "echo start the DOSbox: " , Sys.time() ), 
                                paste( paste0("\"" ,path.dosbox,"\\DOSbox.exe"),"\" -c mount ", path.facets ,"\\mymodel_",time,".bat",sep="")
             )
             
             writeLines( start.dosbox , paste0("C:\\Users\\",user,"\\AppData\\Local\\Temp\\startDOSbox.bat") )
             
             # KEINE BLANKS IN DEN INPUT_FILE (Dateinamen)beim FACDOS!!
             
             
             start.facets <- c( paste( "echo change directory in DOSbox" , Sys.time(),sep=" " ), 
                                # paste( "cd", paste(path.facets,"/",sep="",collapse="") , sep=" ") ,
                                paste( "FACETS BATCH=YES ",model.name,  sep="") ,
                                paste( "echo shutdown DOSbox",  sep="") ,
                                paste( "exit",  sep="")
             )
             writeLines( start.facets , paste0(path.facets,"\\mymodel_",time,".bat") )
             
             shell(paste0("C:\\Users\\",user,"\\AppData\\Local\\Temp\\startDOSbox.bat"),wait=TRUE)
             invisible(file.remove(paste0("C:\\Users\\",user,"\\AppData\\Local\\Temp\\startDOSbox.bat")))
           }else{
             if(!is.null(facetsEXE)){
               Exe <- facetsEXE
             }else{
               Exe <- "Facets.exe"
             }
             start.facets <- c( paste( "echo change directory in DOSbox" , Sys.time(),sep=" " ), 
                                paste( "cd", paste(path.facets,"\\",sep="",collapse="") , sep=" ") ,
                                paste( substring(path.facets,1,2),  sep="") ,
                                paste( "START /w ",Exe," BATCH=YES ",model.name,  sep="") ,
                                #paste( "timeout /t 20 >nul 2>nul",  sep="") ,
                                #paste( "taskkill /im Facets.exe >nul 2>nul",  sep="") ,
                                # paste( "echo shutdown DOSbox",  sep="") ,
                                 paste( ":: exit",  sep="")
             )
             writeLines( start.facets , paste0(path.facets,"\\mymodel_",time,".bat") )
             
             shell(paste0(path.facets,"\\mymodel_",time,".bat"),wait=TRUE,translate=TRUE)
             # system2(paste0(path.facets,"\\mymodel_",time,".bat"),wait=TRUE,invisible = TRUE)
           }
           cat("removing temp-files \n")
           invisible(file.remove(paste0(path.facets,"\\mymodel_",time,".bat")))
           
         },
         
         Linux   = 
         {
           cat("I'm a penguin. Install Facets and Dosbox for Linux \n")
           cat("currently your system is not supported \n")
         },
         Darwin  = 
         {
           cat("I'm a Mac. Install Facets and Dosbox for Mac \n")
           cat("currently your system is not supported \n" )
           
         }
  )
  # output on screen: saved files by FACDOS:
 # cat(paste0("\n The estimation in FACETS is finished, the results ['",inputfile[[grep("Output",inputfile)]],"''] are stored here:\n --> ",path.facets," <-- \n"))
#  cat("Outputs are:\n",
#      paste0(
#        "\n", grep("Output",inputfile,ignore.case=TRUE,value=TRUE),
#        "\n", grep("Score",inputfile,ignore.case=TRUE,value=TRUE),
#        "\n", grep("Xfile",inputfile,ignore.case=TRUE,value=TRUE),
#        "\n", grep("Residual",inputfile,ignore.case=TRUE,value=TRUE))
#  )
  # read files created from FACDOS and store them in a list
  files <- file.info(list.files(path.facets))

  filenames <- sapply(filenames,function(x) substr(x,1,nchar(x)-4))
  fileListe <- sapply(filenames,function(x) grep(x,row.names(files),value = TRUE))
  # lapply(fileListe,FUN=function(x) sapply(x,function(x)readLines(x,skipNul = TRUE,)))
  # namScorefile <- c("T.Score","T.Count","Obs.Avge","Fair.Avge","Measure","S.E.","InfitMS","InfitZ","OutfitMS","OutfitZ","PtMea","PtMeExp","Discrim","Displace","Status","Group","Weight","Lable","Teams","F-Number","F-Label" )
  score <- try(
    lapply(fileListe$Scorefile,function(x){
      all_content = readLines(paste0(path.facets,"/",x))
      skip_second = all_content[-1]
      dat <- read.csv(textConnection(skip_second), header = TRUE, stringsAsFactors = FALSE)  
      return(dat)
    })
  )
  names(score) <- fileListe$Scorefile
  residualfile <- try(
    lapply(fileListe$Residualfile,function(x){
      all_content = readLines(paste0(path.facets,"/",x))
      skip_second = all_content
      dat <- read.csv(textConnection(skip_second), header = TRUE, stringsAsFactors = FALSE)  
      return(dat)
    })
  )
  names(residualfile) <- fileListe$residualfile
  
  names(cate) <- fileListe$categoryfile
  categoryfile <- try(
    lapply(fileListe$Categoryfile,function(x){
      all_content = readLines(paste0(path.facets,"/",x))
      skip_second = all_content
      dat <- read.csv(textConnection(skip_second), header = TRUE, stringsAsFactors = FALSE)  
      return(dat)
    })
  )
  names(categoryfile) <- fileListe$categoryfile
  
  output <- list(
    "Inputfile" = try(readLines(paste0(path.facets,"/",model.name),skipNul = TRUE)),
    "Outputfile" = try(readLines(paste0(path.facets,"/",fileListe$Outputfile),skipNul = TRUE)),
    "Scorefile" = score ,
    "Residualfile" = residualfile[[1]]
  )
  # "Categoryfile" = categoryfile
  invisible( try( file.remove( paste0(path.facets,"/",c(unlist(fileListe),model.name))),silent = TRUE ) )
  
  return(output)
}

#############################################################################			