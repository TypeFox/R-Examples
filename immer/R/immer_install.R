
#####################################################
# installation of the free FACETS DOS_Version
# for this the DOSBOX is needed
immer_install <- function(DosBox_path=NULL,Facets_path=NULL){
  os_system <- Sys.info()["sysname"]
  user <- Sys.info()['user']
  # Ausgabe fuer den User
  switch(os_system,
         Windows = 
{cat("I'm a Windows PC. Install Facets and Dosbox for Windows \n")},
Linux   = 
{cat("I'm a penguin. Install Facets and Dosbox for Linux \n")},
Darwin  = 
{cat("I'm a Mac. Install Facets and Dosbox for Mac \n")}
  )

if(os_system=="Windows"){
  # Link fuer Facets
  link_Facets <- "http://www.winsteps.com/a/facdos.zip"
  # wo hinein
  destination_facets <- 
    paste0("C:\\Users\\",user,"\\Downloads\\facdos.zip")
  # win DOSbox
  link_DosBox <- 
    "http://sourceforge.net/projects/dosbox/files/dosbox/0.74/DOSBox0.74-win32-installer.exe/download"
  # wohinein
  destination_dosBox <- 
    paste0("C:\\Users\\",user,"\\Downloads\\DOSBox0.74-win32-installer.exe")
  installpath <- 
    paste0("C:\\Users\\",user,"\\Documents\\facets")
  
  
  # -----------------------------------------
  # Files herunterladen
  error_facets <- tryCatch(
    download.file(
      url=link_Facets,
      destfile=destination_facets,
      method="internal"
    )
  )
  error_DosBox <- tryCatch(
    download.file(
      url=link_DosBox,
      destfile=destination_dosBox,
      method="internal"
    )
  )
  # -----------------------------------------
  
  
  # Den Admin des Computers herausfinden: und Installation von DosBox
  # -----------------------------------------
    if(error_DosBox[[1]]!=0){
      cat("Attention, there was an error while downloading the DosBox, 
          please try again or install die DosBox manually \n")
      cat(paste0("for the manual installation pleas go to: \n",link_DosBox,"\n",
                 "after the download process finished we recomand to install 
                 the DosBox into \n--> \\Users\\yourUser\\Documents <--"))
    }
    if(error_DosBox[[1]]==0){
      cat("starting installation process of DosBox \n")
      admin <- system("net localgroup",intern=TRUE)
      findeAdmin <- paste0("net localgroup ",gsub("\\*","",admin[grep("adm|Adm",admin)]))
      test <- system(findeAdmin,intern=TRUE)
      administrator <- test[grep("--",test)+1]
      system(
        paste0("runas /noprofile /user:",administrator," DOSBox0.74-win32-installer.exe"),
        invisible=FALSE,
        wait = FALSE)
      # -----------------------------------------
      # Edit the configFile to speed up the process [if the installation is successful]
      if(!is.na(match("DOSBox",list.files(paste0("C:\\Users\\",user,"\\AppData\\Local\\"))))){
        config <- readLines(paste0("C:\\Users\\",user,"\\AppData\\Local\\DOSBox\\dosbox-0.74.conf"))
        config[grep("cycles=",config)] <- "cycles=max"
        writeLines(config,paste0("C:\\Users\\",user,"\\AppData\\Local\\DOSBox\\dosbox-0.74.conf"))
        }
    }
    if(error_facets[[1]]!=0){
      cat("Attention, there was an error while downloading facets, 
          please try again or install die DosBox manually \n")
      cat(paste0("for the manual installation pleas go to: \n",link_Facets,"\n",
                 "after the download process finished we recomand to unzip the Folder and 
                 move the content to \n-->",installpath,"<--"))
    }
    if(error_facets[[1]]==0){
      cat("unzip facets\n")  
      utils::unzip (destination_facets, exdir = installpath)
      cat("move facets to ",installpath,"\n") 
      # ----------------------------------------- 
    }
  } # End Windows

  # LINUX
  if(os_system=="Linux"){
    # Link fuer Facets
    link_Facets <- "http://www.winsteps.com/a/facdos.zip"
    # wo hinein
    destination_facets <- paste0("/home/",user,"/Downloads/facdos.zip")
    # win DOSbox
    link_DosBox <- "http://sourceforge.net/projects/dosbox/files/dosbox/0.74/dosbox-0.74.tar.gz/download"
    # wohinein
    destination_dosBox <- paste0("C:\\Users\\",user,"\\Downloads\\dosbox-0.74.tar.gz")  
    
    cat("\n install automake for dosbox \n")
    system("sudo -kS apt-get -y install build-essential autoconf automake",input=readline("Enter your password: "))
  
    cat("\n install dosbox \n")
    system("sudo -kS apt-get -y install build-dep dosbox",input=readline("Enter your password: "))
    cat("\n DONE \n")
    
    cat("\n unzip facdos to ",paste0("/home/",user,"/facdos \n"))
    system(paste0("unzip -o ", destination_facets," -d /home/",user,"/facdos"))
    
    cat("\n now you will find dosbox and facdos in:",paste0("/home/",user,"/facdos"))
  } # End Linux

  # MAC-OS X
  if(os_system=="Darwin"){
    # Link fuer Facets
    link_Facets <- "http://www.winsteps.com/a/facdos.zip"
    # wo hinein
    destination_facets <- paste0("C:\\Users\\",user,"\\Downloads\\facdos.zip")
    # win DOSbox
    link_DosBox <- "http://sourceforge.net/projects/dosbox/files/dosbox/0.74/DOSBox-0.74-1_Universal.dmg/download"
    # wohinein
    destination_dosBox <- paste0("C:\\Users\\",user,"\\Downloads\\DOSBox-0.74-1_Universal.dmg")  
    
    # -----------------------------------------
    # Files herunterladen
    error_facets <- tryCatch(
      download.file(
        url=link_Facets,
        destfile=destination_facets,
        method="internal"
      )
    )
    error_DosBox <- tryCatch(
      download.file(
        url=link_DosBox,
        destfile=destination_dosBox,
        method="internal"
      )
    )
    # -----------------------------------------
  } # End Mac OS X

}
#############################################################################  		