# SPECTRINO V1.6.0 - binary package distribution
# variable types
#
# Grp,Spc: string (the name) or integer (the index)
# working with names is saver, with indexes - easier
# if Grp=0 then Grp=(Active Group Index); if Spc=0 then Spc=(Active Spec Index)
# Grp="*" or Spc="*" masks all of them
#
# OnlyChecked, Checked, InclOpt, NewGrp: boolean;
# *Filename, GDirectory: string  "
# .. in the beginning will be replaced by Spectrino's path
# All the filenames or directories must be with forward slashes e.g. D:/Prime/Data/Test.txt

# FOR INTERNAL USE ONLY
.spnDebug <- FALSE
.spnVersion <- "spnVer160"
.spnConnectError <- "Error: The connection to Spectrino app is not active or Spectrino server has been disrupted."

.sprnApp <- list(
  call = function(req) {
    wsUrl = paste(sep='',
                  '"',
                  "ws://",
                  ifelse(is.null(req$HTTP_HOST), req$SERVER_NAME, req$HTTP_HOST),
                  '"')
    
    list(
      status = 200L,
      headers = list(
        'Content-Type' = 'text/html'
      ),
      body = 
        'Hello Spectrino\r\n'
    )
  },
  onWSOpen = function(ws) {
    assign(".wsSprn", ws, envir = .spnEnv) 
    assign(".wsConnected", TRUE, envir = .spnEnv)    
    ws$onMessage(function(binary, message) {  
      assign(".wsSprnFlag", FALSE, envir = .spnEnv)
      assign(".wsSprnText", message, envir = .spnEnv)                                         
    })
    ws$onClose(function() { 
      assign(".wsConnected", FALSE, envir = .spnEnv) 
      print("closing websocket connection to Spectrino app.")
    })
  }
)

.cpos <- function(str,sub,start=1) {
  # find the first position of string sub in string str, starting from position start
  lstr  <- nchar(str)
  lsub1 <- nchar(sub)-1
  if (start+lsub1 > lstr) return(NA)
  else {
    str <- substring(str,start,lstr)
    str <- substring(str, 1:(lstr-lsub1), (1+lsub1):lstr)
    p <- ((start:lstr)[str==sub])[1]
    if (is.na(p>0)) return(NA)
    else return(p)
  }
}

.spnIsExec <- function () {
  spnFile <- paste(path.package("spectrino"),"/exec/spectrino.exe",sep="")
  return(file.exists(spnFile))
}

.spnUnpack <- function (msg) {
  # B - boolean; I - integer; D - double; S - string; E - error  
  # if number after the letter - array of that lenght  
  i = .cpos(msg," ",start = 1)
  st = substr(msg,1,1)
  if(st == "E") return(paste("Error:",substr(msg,3,1000))) 
  sb <- substr(msg,2,i-1) 
  if(st == 'M') {
    j <- .cpos(sb,":",start = 1)
    sn <- as.integer(substr(sb,1,j-1))
    sn2 <- as.integer(substr(sb,j+1,100))
    }
    else sn = as.integer(sb)
  ss = substr(msg,i+1,nchar(msg))
  if(.spnDebug) print(paste("st=",st," sn=",sn," ss=",ss))  
  if(st == 'M') bb <- FALSE
    else bb <- (sn == 1)
  if(bb) {
    switch(st,
           B = as.logical(ss),
           I = as.integer(ss),
           D = as.double(ss),
           S = as.character(ss)
          ) 
  }
  else {
    sa <- fromJSON(ss)
    if(!(st == 'M'))
      if(!(length(sa) == sn)) return('Error: wrong number of elements in array')
    return(sa)
  }
}

.spnCommand <- function(command) { 
  if (!get(".wsConnected", envir = .spnEnv)) return(list(stts = FALSE, rslt = "Error: no connection to spectrino app"))
  if(.spnDebug) print(paste("sent:",command))
  get(".wsSprn", envir = .spnEnv)$send(command)  
  # follow by the command itself - command(prm1, prm2, ...)
  assign(".wsSprnFlag", TRUE, envir = .spnEnv); i = -1
  while (get(".wsSprnFlag", envir = .spnEnv)) {  i <- i + 1   
    service();  Sys.sleep(0.001); 
    if (i > get(".sprnOpt", envir = .spnEnv)$timeOut) break 
  } 
  if (i > get(".sprnOpt", envir = .spnEnv)$timeOut) return(list(stts = FALSE, rslt = "Error: Cannot connect to Spectrino app: TIME OUT"))
  msg <- get(".wsSprnText", envir = .spnEnv)
  if(.spnDebug) print(paste("back:", msg))
  if(!is.character(msg)) return(list(stts = FALSE, rslt = "Error: Wrong result format from Spectrino app"))
  status = !(substr(msg,1,1) == "E")    
  list(stts = status, rslt = .spnUnpack(msg))
}

# Check for spectrino object and (optionally) spectrino application
spnCheck <- function (inclApp = FALSE) {
  bb <- TRUE
  if(inclApp) { cmd = .spnCommand(.spnVersion)
    bb <- (cmd$stts) && (cmd$rslt)
  }  
  for(i in 1:3) { service(); Sys.sleep(0.001); }; Sys.sleep(0.1)
  return(exists(".wsSprn", envir = .spnEnv) &&  get(".wsConnected", envir = .spnEnv) && bb) 
}

# Check if the result of command is Error
spnIsError <- function (rslt) {
  res <- FALSE
  if(is.list(rslt)) res <- (substr(rslt$rslt,1,6) == "Error:")
  if(is.character(rslt)) res <- (substr(rslt,1,6) == "Error:")
  return(res)
}

spnInstallApp <- function (zipFile = "") {
  if(.spnIsExec()) {
    cat("The application Spectrino is already installed. \n")
    q <- readline(prompt="Enter I to reInstall it; OR any other key to skip: ") 
    if(!((q == 'i') || (q == 'I'))) return(TRUE)
  }  
  if(is.null(zipFile)) return("void installation")
  spnDir <- path.package("spectrino")
  if(zipFile == "") {
    temp <- tempfile()
    dw <- try(download.file("http://www.spectrino.com/down/exec.zip",temp), silent = TRUE) 
    if (dw == 0) {
       unzip (temp, exdir = spnDir)
       unlink(temp) 
    } else {
        cat("Cannot download Spectrino application from spectrino.com, check your i-net connection. \n")
        return(FALSE)
      }
  } else {
    if(! file.exists(zipFile)) {
      cat("Cannot find the local zip file, see the help.")
      return(FALSE)
    }
    unzip (zipFile, exdir = spnDir)
  }
  return(.spnIsExec())
}

# Initialization - opens and initializes Spectrino application
spnNew <- function(TimeOut = 2, Host = "127.0.0.1", Port = 9876) { 
  if(! .spnIsExec()) {
    cat("The application Spectrino is missing. For installation use spnInstallApp (see help).")
    q <- readline(prompt="Enter I to run spnInstallApp(); OR any other key to quit: ")
    if((q == 'i') || (q == 'I')) { 
      TimeOut = 5
      assign(".sprnOpt", list(timeOut = 100*TimeOut, host = Host, port = Port), envir = .spnEnv)       
      spnInstallApp()
    }                               
    if(! .spnIsExec()) return(FALSE)
  }
  # clean-up left-overs
  if (get(".wsConnected", envir = .spnEnv)) return("Error: Connection to Spectrino app is still active. Close it with <spnFree(T)>.") 
  # the good stuf 
  spnFile <- paste(path.package("spectrino"),"/exec/spectrino.exe",sep="")
  shell.exec(spnFile)
  assign(".sprnOpt", list(timeOut = 100*TimeOut, host = Host, port = Port), envir = .spnEnv) 
  if(.spnDebug) print(get(".sprnOpt", envir = .spnEnv))
  assign(".wsSprnHandle", startServer(Host, Port, .sprnApp), envir = .spnEnv)
  assign(".wsSprnFlag", TRUE, envir = .spnEnv); i = -1; 
  while (get(".wsSprnFlag", envir = .spnEnv)) {  i <- i + 1 
    service();  Sys.sleep(0.001); 
    if (i > (2*get(".sprnOpt", envir = .spnEnv)$timeOut)) break 
  }
  if (i > (2*get(".sprnOpt", envir = .spnEnv)$timeOut)) return("Error: Cannot start Spectrino app: Time out")
  cmd <- .spnCommand(.spnVersion)
  if(.spnDebug) print(cmd)
  if ((!cmd$stts) || (!cmd$rslt)) return("Error: Failed handshake with Spectrino app") 
    else return( TRUE ) 
}

# Get number of the groups loaded
spnGetGrpCount <- function() {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand("spnGetGrpCount()")$rslt
}

# Get number spec in Grp group
spnGetSpcCount <- function(OnlyChecked = FALSE, Grp = 0) {
  if (!spnCheck()) return(.spnConnectError)
 .spnCommand(paste("spnGetSpcCount(",OnlyChecked,",",Grp,")",sep=""))$rslt 
}

# Get group name ; if GrpIdx="*" gets back a list of all
spnGetGrpName <- function(GrpIdx = "*") {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnGetGrpName(",GrpIdx,")",sep=""))$rslt
}

# Get spec name ; if SpcIdx="*" gets back a list of all
spnGetSpcName <- function(Grp = 0,SpcIdx = "*")  {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnGetSpcName(",Grp,",",SpcIdx,")",sep=""))$rslt
}

# Get reference X set of values (vector)
spnGetRefer <- function() {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand("spnGetRefer")$rslt
}

# Get one spectrum (vector)
spnGetSpc <- function(Grp = 0, Spc = 0) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnGetSpc(",Grp,",",Spc,")",sep=""))$rslt
}

# Get spectra from one group (matrix)
spnGetGrp <- function(OnlyChecked = FALSE, Grp = 0) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnGetGrp(",OnlyChecked,",",Grp,")",sep=""))$rslt
}

# Get spectra from all the groups (matrix); WITHOUT "Unknows" group
spnGetTree <- function(OnlyChecked = FALSE) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnGetTree(",OnlyChecked,")",sep=""))$rslt
}

# Get the boolean vector of the state of checking boxes of Grp group
spnGetSpcChecked <- function(Grp = 0) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnGetSpcChecked(",Grp,")",sep=""))$rslt
}

# Set Spc spectrum of Grp group checkbox to checked/unchecked state
spnSetSpcChecked <- function(Grp = 0, Spc = 0, Checked = FALSE) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnSetSpcChecked(",Grp,",",Spc,",",Checked,")",sep=""))$rslt
}

# Open SFilename in Grp group. result = spnGetSpcCount(false,Grp)
spnOpenSpc <- function(Grp = 0, SFilename) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnOpenSpc(",Grp,",",SFilename,")",sep=""))$rslt
}

# Open a group from GFilename file. result = spnGetGrpCount
spnOpenGrp <- function(GFilename,NewGrp = FALSE) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnOpenGrp(",GFilename,",",NewGrp,")",sep=""))$rslt
}

# Open spec-tree from TFilename. result = spnGetGrpCount
# InclOpt = 0 -> factory setting (no preproc.); 1 -> last used; 2 -> from LFilename
spnOpenTree <- function(TFilename, InclOpt = 0) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnOpenTree(",TFilename,",",InclOpt,")",sep=""))$rslt
}

# Save Grp group as GFilename file (the path is ignored)
# if GFilename="" then (use its proper name) else (rename and save)
# if Grp=* then save all groups (GFilename is ignored, but must be pressent ;)
spnSaveGrp <- function(Grp = 0, GFilename = "") {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnSaveGrp(",Grp,",",GFilename,")",sep=""))$rslt
}

# Save spec-tree along with the preprocessing options
# if LFilename='' then use proper name
spnSaveTree <- function(TFilename = "") {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnSaveTree(",TFilename,")",sep=""))$rslt
}

# Delete Spc spectrum from Grp group. result = spnGetSpcCount(false,Grp)
# if Spc="*" then all of them
spnDelSpc <- function(Grp = 0, Spc = 0) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnDelSpc(",Grp,",",Spc,")",sep=""))$rslt
}

# Delete Grp group. result = spnGetGrpCount
# if Grp="*" then all of them
spnDelGrp <- function(Grp = 0) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnDelGrp(",Grp,")",sep=""))$rslt
}

# Show/Hide Spectrino
spnSetVis <- function(Visible = TRUE) {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnSetVis(",Visible,")",sep=""))$rslt
}

# Get/Set active group; if Grp=0 only get Active group; else set one
spnActGrp <- function(Grp = 0){
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste("spnActGrp(",Grp,")",sep=""))$rslt
}

# Set Spectrino pre-processing options
spnSetPPOpt <- function(OptionList = "") {
  if (!spnCheck()) return(.spnConnectError)
  .spnCommand(paste('spnSetPPOpt("',OptionList,'")',sep=""))$rslt
}

# Validation (not conclusive, only the most common functions and modes)
spnValidation <- function() {
  if (!spnCheck()) return(.spnConnectError)
  if (spnDelGrp("*")!=0) return("Error: #1")
  if (spnOpenTree("<test>",0)!=3) return("Error #2")
  if (sum(as.numeric(spnGetRefer()))!=20100) return("Error: #3")
  if (spnGetSpcCount(FALSE,"Test1")!=3) return("Error #4")
  if (sum(as.numeric(spnGetSpc("Test1",2)))!=1274) return("Error: #5")
  if (!(spnSetSpcChecked("Test2","*",TRUE))) return("Error: #6")
  if (!(spnSetSpcChecked("Test2",3,FALSE))) return("Error: #7")
  if (sum(as.numeric(spnGetSpcChecked("Test2")))!=2) return("Error: #8")
  if (sum(as.numeric(spnGetGrp(FALSE,"Test2")))!=6346) return("Error: #9")
  if (sum(as.numeric(spnGetGrp(TRUE,"Test2")))!=4442) return("Error: #10")
  if (sum(as.numeric(spnGetTree(FALSE)))!=18002) return("Error: #11")
  if (sum(as.numeric(spnGetTree(TRUE)))!=8250) return("Error: #12")
  if (spnDelSpc("Test2",2)!=2) return("Error: #13")
  if (sum(as.numeric(spnGetGrp(TRUE,"Test2")))!=1904) return("Error: #14")
  if (spnDelGrp("Test2")!=2) return("Error: #15")
  if (sum(as.numeric(spnGetTree(FALSE)))!=11656) return("Error: #16")
  if (spnOpenGrp("Test2")!=3) return("Error: #17")
  if (sum(as.numeric(spnGetTree(FALSE)))!=18002) return("Error: #18")

  cat(" Active group: ",spnGetGrpName(0),"\t"," All groups: ",spnGetGrpName("*"),"\n")
  cat(" Specs in the active group: \n",spnGetSpcName(0,"*"),"\n\n" )
  cat("Validation confirmed. \n") 
}

#Closes Spectrino object and application
spnFree <- function(inclApp = FALSE) {
  if (!spnCheck()) return(.spnConnectError)
  if (inclApp) get(".wsSprn", envir = .spnEnv)$send("spnQuit()")
  stopServer(get(".wsSprnHandle", envir = .spnEnv))
  Sys.sleep(1)
}
