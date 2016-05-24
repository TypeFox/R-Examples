"getWinEdt" <-
function(){
    WinEdtRegistryKey <- file.path("SOFTWARE", "WinEdt 8", fsep="\\")
    WinEdtReg <- try(readRegistry(WinEdtRegistryKey, hive = "HCU",
                                  maxdepth = 1), silent = TRUE)
    WinEdtVersion <- 8
    if(inherits(WinEdtReg, "try-error")){
        WinEdtRegistryKey <- file.path("SOFTWARE", "WinEdt 7", fsep="\\")
        WinEdtReg <- try(readRegistry(WinEdtRegistryKey, hive = "HCU",
                                      maxdepth = 1), silent = TRUE)
        WinEdtVersion <- 7
    }
    if(inherits(WinEdtReg, "try-error")){    
        WinEdtRegistryKey <- file.path("SOFTWARE", "WinEdt 6", fsep="\\")
            WinEdtReg <- try(readRegistry(WinEdtRegistryKey, hive = "HCU",
                                          maxdepth = 1), silent = TRUE)
        WinEdtVersion <- 6                                      
    }                                  
    if(inherits(WinEdtReg, "try-error")){
      ## WinEdt 5.x
        WinEdtVersion <- 5
        WinEdtRegistryKey <- file.path("SOFTWARE", "WinEdt", fsep="\\")
        WinEdtReg <- try(readRegistry(WinEdtRegistryKey, hive = "HCU",
                                      maxdepth = 1), silent = TRUE)
        if(!inherits(WinEdtReg, "try-error")){
            InstallRoot <- WinEdtReg[["Install Root"]]
            ApplData <- WinEdtReg[["ApplData"]]
        } else{
            temp <- try(readRegistry(WinEdtRegistryKey, hive = "HLM",
                                     maxdepth = 1, view = "32-bit"), silent = TRUE)
            temp2 <- try(readRegistry(file.path("SOFTWARE", "Team WinEdt", fsep="\\"),
                                      hive = "HLM", maxdepth = 1, view = "32-bit"), silent = TRUE)
            if(!(inherits(temp, "try-error") && inherits(temp2, "try-error"))){
                cat("\nStart WinEdt manually once before first usage.\n")
                stop("Start WinEdt manually once before first usage.")
            } else
                stop("WinEdt is not installed properly.", "\n",
                    "Either reinstall WinEdt or install R-WinEdt manually as described in the ReadMe")
        }
        if(!is.null(ApplData)){
            ApplData <- file.path(ApplData, "WinEdt", fsep="\\")
        } else
            ApplData <- InstallRoot
        RWinEdtInstalled <- RWinEdtVersion <- file.exists(file.path(InstallRoot, "R.ver", fsep = "\\"))
        if(RWinEdtVersion)
            RWinEdtVersion <- scan(file.path(InstallRoot, "R.ver", fsep = "\\"), quiet = TRUE)
    } else {
      ## WinEdt 6-8
        InstallRoot <- WinEdtReg[["Install Root"]]
        ApplData <- WinEdtReg[["AppData"]]
        if(is.null(ApplData)){
            stop("\n", "Before the first usage, you have to start WinEdt manually once or twice and create a personal profile")
        }
        RWinEdtInstalled <- RWinEdtVersion <- file.exists(file.path(ApplData, "R.ver", fsep = "\\"))
        if(RWinEdtVersion)
            RWinEdtVersion <- scan(file.path(ApplData, "R.ver", fsep = "\\"), quiet = TRUE)
    }
    return(list(InstallRoot = InstallRoot, RWinEdtInstalled = RWinEdtInstalled,
                RWinEdtVersion = RWinEdtVersion, ApplData = ApplData,
                WinEdtVersion = WinEdtVersion))
}
