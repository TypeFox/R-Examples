#####EXAMPLES:

##note: to avoid misunderstandings: 'SMHandler' stands for /S/tartup/M/essage/Handler/

#mySMHandler <- function(c) {
#    pkg <- startupPackage(c)
#    npkg <- nchar(pkg)
#    linestarter <- paste(":",pkg,"> ", sep ="")                              
#    linestarterN <- paste("\n",linestarter, sep ="")
#    linestarterE <- paste(linestarterN,"$",sep="")
#    writeLines(paste(linestarter, sub(linestarterE,"\n", 
#                     gsub("\n", linestarterN, 
#                          conditionMessage(c))),sep=""),stderr())                                               
#}


mystartupMessage <- function(..., domain = NULL, pkg = "", type = "version", 
                             SMHandler = mySMHandler, endline = FALSE){
 withRestarts(withCallingHandlers(
                     startupMessage(..., domain = domain, 
                                    pkg = pkg, type=type, endline = endline), 
                     StartupMessage=function(m)
                            {signalCondition(m)
                             invokeRestart("custom",c=m,f=SMHandler)}
                                 ),
                                               #as suggested by Seth Falcon:
              onlytypeMessage = function(c0,atypes)
                              {if(startupType(c0) %in% atypes) 
                                      SMHandler(c=c0)    
                              },                                
              #as suggested by Seth Falcon:
              custom = function(c,f) f(c),
              muffleMessage = function() NULL )
 invisible(NULL) 
}

buildStartupMessage <- function(..., pkg, library=NULL, domain=NULL,
                                packageHelp=FALSE, MANUAL = NULL, 
                                VIGNETTE = NULL,
                                SMHandler=mySMHandler){
#
tit.vers <- readVersionInformation(pkg,library)
if((!getOption("StartupBanner")=="off")||is.null(getOption("StartupBanner"))) 
       mystartupMessage(tit.vers$"title", " (version ", tit.vers$"ver", ")", 
                        domain = domain, pkg = pkg, type="version", 
                        SMHandler = SMHandler)
###
if((getOption("StartupBanner")=="complete")||
    is.null(getOption("StartupBanner"))){ 
     llist <- length(list(...))
     ### checks as to existence of URL- NEWS- and MANUAL-information
     #
     URL <- readURLInformation(pkg,library)
     NEWS <- pointertoNEWS(pkg,library) 
     #
     if ( packageHelp)  packageHelpS <- c("?\"", pkg, "\"")
               else     packageHelpS <- ""
     if (!is.null(NEWS))  NEWSS <- NEWS
                 else     NEWSS <- ""
     if (!is.null(URL))   URLS <- c("\n  ",URL)
                 else     URLS <- ""
     
     ## MANUALL : is there a MANUAL entry?
     MANUALL <- FALSE
     MANUALS <- ""
     if(!is.null(MANUAL))
        {if (all(substr(as.character(MANUAL),1,7)=="http://")) 
               {MANUALL <- TRUE
                MANUALS <- c("\n  ",MANUAL)}
         else  {MANUAL1 <- paste(MANUAL, 
                                 sep = .Platform$file.sep,
                                 collapse = .Platform$file.sep)
                MANUALpath <- file.path(system.file(package = pkg),
                                        MANUAL1, collapse = "")
                if (file.exists(MANUALpath)) 
                              { MANUALL <- TRUE
                                 MANUALS <- c("\n  ",MANUALpath)}  
                }
         }                          
     VIGNETTES = ifelse(!is.null(VIGNETTE), 
                         paste("\n",VIGNETTE, sep = "", collapse = ""), "")

     ## are there any info-lines?
     L <- sum(!is.null(URL), packageHelp , !is.null(NEWS) , MANUALL, 
              !is.null(VIGNETTE))
     
     ##determining the separators:
     seps <- character(3)
     seps[1] <- ifelse(packageHelp&&L>1,", ","")
     seps[2] <- ifelse(!is.null(NEWS)&&
                        sum(!is.null(NEWS) , MANUALL, !is.null(URL))>1,
                        gettext(", as well as", domain = domain),
                        "")
     seps[3] <- ifelse(MANUALL && sum(MANUALL, !is.null(URL))>1,
                       ", ", "")
     if( (MANUALL|| !is.null(URL)) && is.null(NEWS)) 
          seps[1] <- gettext(", as well as", domain = domain) 
     #
     if (L>0){
          if (llist > 0)
          mystartupMessage(..., domain=domain, pkg=pkg, type="notabene", 
                          SMHandler=SMHandler)

          mystartupMessage("For more information see ", 
                         packageHelpS, seps[1], NEWSS, seps[2], URLS, seps[3], 
                         MANUALS, VIGNETTES, "\n",  
                         domain = domain, pkg = pkg, type = "information", 
                         SMHandler = SMHandler, endline = TRUE)
    }
    else{
          if (llist > 0)
          mystartupMessage(..., domain=domain, pkg=pkg, type="notabene", 
                          SMHandler=SMHandler, endline = TRUE)
    } 
   }
}
   
########### end Examples
