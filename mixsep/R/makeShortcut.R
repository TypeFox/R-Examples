makeShortcut <- function(dbfile="",wang2desktop=FALSE){
  if(grepl("win",tolower(.Platform$OS.type))){ ## Windows
    mixsepdir <- find.package("mixsep")
    ## Create runmixsep.R script
    if(dbfile!=""){ ## If a file of DB connector information is given:
      dblines <- readLines(file(dbfile,'r')) 
      comments <- unlist(lapply(dblines,grepl,pattern="#"))
      ## removes comments from lines/or entire comment lines
      if(any(comments)){
        dblines[comments] <- unlist(lapply(strsplit(dblines[comments],"#"),function(x) x[1]))
      }
      ## removes leading and tailing white space
      dblines <- gsub("^\\s+|\\s+$", "", dblines)
      dblines <- dblines[dblines!=""]
      if(length(dblines)!=4){
        warning("Something went wrong when setting database parameters")
        cat(" library(mixsep)\n mixsep()",file=paste(mixsepdir,"R","runmixsep.R",sep=.Platform$file.sep))
      }
      else{
        ## Writes information to runmixsep.R
        dbcols <- grepl("dbcols",dblines)
        dbcol <- unlist(lapply(strsplit(dblines[dbcols],"="),function(x) x[2]))
        dbcol <- gsub("^\\s+|\\s+$","",unlist(strsplit(gsub("\"", "", dbcol),";")))
        dblines[dbcols] <- paste("dbcols = c(\"",paste(gsub("^\\s+|\\s+$","",unlist(strsplit(gsub("\"", "", dbcol),";"))),collapse="\",\"",sep=""),"\")",sep="",collapse="")
        cat(" library(mixsep)\n options(mixsep=list(",paste(dblines,collapse=","),"))\n mixsep()",file=paste(mixsepdir,"R","runmixsep.R",sep=.Platform$file.sep))
      }
    }
    ## If no DB connector is supplied - This disables [DB] in GUI
    else cat(" library(mixsep)\n mixsep()",file=paste(mixsepdir,"R","runmixsep.R",sep=.Platform$file.sep))
    ## Change from 0.1-1 to 0.1-2
    ## Changes due to difference between R2.11 and R2.12 where Rterm is moved to i386/ or x64/ (architechture dependent).
    ## Hence, rather than using "Rterm.exe -q -f ../runmixsep.R" the call is now "Rscript.exe ../runmixsep.R"
    ## Create .bat-file: 
    cat(paste(' @echo off \n start /MIN ', paste(R.home(),"bin","Rscript.exe",sep=.Platform$file.sep), ' ',
               paste(mixsepdir,"R","runmixsep.R",sep=.Platform$file.sep), '\n'),
        file=paste(mixsepdir, "R", "mixsep.bat",sep=.Platform$file.sep))
    ## Create Visual Basic Script for shortcut creation:
    cat(paste('Set objShell = WScript.CreateObject("WScript.Shell") \n
             strDesktopFld = objShell.SpecialFolders("Desktop") \n
             strMyDocumentFld = objShell.SpecialFolders("MyDocuments") \n
             Set objMixsepShortcut = objShell.CreateShortcut(strDesktopFld & "\\mixsep.lnk") \n
             objMixsepShortcut.TargetPath = "', mixsepdir, '/R/mixsep.bat" \n
             objMixsepShortcut.WorkingDirectory = strMyDocumentFld \n
             objMixsepShortcut.Save \n', sep=""), file=paste(mixsepdir,"R","shortcut.vbs",sep=.Platform$file.sep))
    ## Check if user wants path/to/mixsep/package/data/wang.csv copied to the desktop:
    if(wang2desktop){
      wangfile <- paste(mixsepdir,"data","wang.csv",sep=.Platform$file.sep)
      cat(paste('\n\n dim filesys \n
      set filesys=CreateObject("Scripting.FileSystemObject") \n
      If filesys.FileExists(\"',wangfile,'\") Then \n
      filesys.CopyFile \"',wangfile,'\", strDesktopFld + "',.Platform$file.sep,'", true \n
      End If',sep=""), file=paste(mixsepdir,"R","shortcut.vbs",sep=.Platform$file.sep), append=TRUE)
    }
    ## Create shortcut and/or copy wang.csv to desktop by calling VB script:
    shell(paste("CSCRIPT",paste(mixsepdir,"R","shortcut.vbs",sep=.Platform$file.sep)))
    ## Remove VB script:
    file.remove(paste(mixsepdir,"R","shortcut.vbs",sep=.Platform$file.sep))
  }
  else stop("Creation of shortcut is only possible for Windows operating systems")
}
