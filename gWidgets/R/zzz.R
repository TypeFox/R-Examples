##' function to bypass call to require to keep things quiet
##'
##' @param package name
##' @return result of call to require(pkg)
.bypassRequire <- function(pkg) {
  do.call(sprintf("%s", "require"), list(pkg))
}

## needed to find methods
.onLoad <- function(lib, pkg) {

##   ## Bad practice here (see ?.onAttach), but don't want to require tcltk although likely it isn't too much to ask for
##   doRequire <- function(pkg) do.call("require",list(pkg))
##   popup <- function() {
##     ## popup install message if not presnet
##     ## already checked that no gWidgetsXXX is available
##     all <- utils:::installed.packages()
##     pkgs <- rownames(all)
    
##     title <- "gWidgets needs a toolkit package"
##     msg <- "gWidgets needs a toolkit package installed."

##     if("tcltk" %in% pkgs) {
##       if("RGtk2" %in% pkgs) {
##         msg <- paste(msg, "Try installing gWidgetsRGtk2.", sep="\n")
##       } else {
##         msg <- paste(msg, "Try installing gWidgetstcltk.", sep="\n")
##       }

##       installing_gWidgets_toolkits()
##       doRequire("tcltk")
##       w <- tktoplevel()
##       tkdialog(w, title,msg,"",0,"close")
      
##     } else {
##       msg <- paste(msg,
##                    "You must install either RGtk2 or tcltk,",
##                    "then a toolkit package.",
##                    sep="\n")
##       packageStartupMessage(msg, "\n")

##       installing_gWidgets_toolkits()
##     }

##   }
  
  
  
##   ## check that a toolkit package is loaded 
## #  tmp = installed.packages() 
## #  choices = tmp[grep("^gWidgets.",tmp[,1]),1]
## #  if(length(choices) == 0) {
##     ### popup() ## disabled -- had issues iwth automatic builds

##     msg <- gettext("gWidgets requires a toolkit implementation to be\n installed, for instance gWidgetsRGtk2 or gWidgetstcltk.")
##     warning(paste("\n\n**",msg,"**\n\n", sep=""))
## #  }

  
}

.onAttach <- function(...) {

}
