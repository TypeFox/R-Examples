data.collection.obj <- new("RzDataCollection")
rzSettings <- new("RzSettings")
rzTools    <-new("RzTools")
column.definition <- c(index=0, select=1, vars=2, var.labs=3, msr=4, msr.image=5, val.labs=6, missing=7)

#      .Rz.path <- "/home/masahiro/Documents/R/Rz/Rz/inst"  # for debug
#      .Rz.path <- "/media/sf_Dropbox/Documents/R/Rz/Rz/inst"  # for debug

.onAttach <- function(lib, pkg){
  if(grepl("mingw", R.Version()$os)){
    temp<-try(winMenuAdd("Rz"),silent=TRUE)
    if(class(temp)!="try-error"){
      winMenuAddItem("Rz", gettext("Start"), "Rz()")
    }
  }
  rzSettings$load()
  gtkAccelMapLoad(file.path(rzSettings$getRzPath(), "themes", "Default", "gtk-2.0-key", "key.txt"))
  
  txt1 <- "################################ Rz ################################"
  txt2 <- gettext("Excute Rz() to start,")
  txt3 <- gettext("or you can start from menu bar if you use R on Windows.")
  txt4 <- "####################################################################"
  txt  <- format(c(txt1,txt2,txt3,txt4), justify="centre")
  txt  <- paste(txt, collapse="\n")
  packageStartupMessage(txt)
  if (interactive()) {
    Rz()
    checkConfDir()
  }
}

