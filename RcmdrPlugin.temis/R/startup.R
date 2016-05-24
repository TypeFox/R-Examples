.onAttach <- function(libname, pkgname){
  # Work around a bug with Java on Mac OS blocking Snowball and RWeka
  if(Sys.info()["sysname"] == "Darwin")  {
      Sys.setenv(NOAWT="true")
  }

  # Work around a bug in JGR: restore the correct locale
  # https://www.rforge.net/bugzilla/show_bug.cgi?id=244
  if(Sys.info()["sysname"] == "Darwin" && Sys.getlocale() == "C")  {
      loc <- system("defaults read -g AppleLocale", intern=TRUE)
      Sys.setlocale("LC_ALL", loc)
  }

  if (!interactive()) return()


  Rcmdr <- options()$Rcmdr
  plugins <- Rcmdr$plugins

  if (!pkgname %in% plugins) {
      Rcmdr$plugins <- c(plugins, pkgname)
      Rcmdr$ask.on.exit <- FALSE
      options(Rcmdr=Rcmdr)

      if("package:Rcmdr" %in% search()) {
          if(!getRcmdr("autoRestart")) {
              closeCommander(ask=FALSE, ask.save=TRUE)
              Commander()
          }
      }
      else {
          if(packageVersion("Rcmdr") < "2.0-0")
              stop(.gettext("This package requires Rcmdr version 2.0-0 or higher."))

          Commander()
      }
  }

  if (getRcmdr("use.markdown")) {
      putRcmdr("startNewCommandBlock", FALSE)
      beginRmdBlock()
  }
  if (getRcmdr("use.knitr")) {
      putRcmdr("startNewKnitrCommandBlock", FALSE)
      beginRnwBlock()
  }

  # HTML.matrix() does not allow passing scientific separately,
  # and vocabulary summary tables often end up printed in scientific notation
  doItAndPrint(.gettext("# Prefer fixed to scientific notation"))
  doItAndPrint('options(scipen=5)')

  doItAndPrint("")
  doItAndPrint(.gettext("# Print numbers with two significant digits"))
  doItAndPrint('options(digits=2)')
  doItAndPrint('options(R2HTML.format.digits=2)')

  # Keep in sync with disableBlackAndWhite()
  # We can stop specifying region when latticeExtra uses RColorBrewer:: for its default value:
  # https://r-forge.r-project.org/tracker/index.php?func=detail&aid=4853&group_id=232&atid=942
  doItAndPrint("")
  doItAndPrint(.gettext("# Set a nice color palette for plots"))
  doItAndPrint('lattice.options(default.theme=latticeExtra::custom.theme(symbol=RColorBrewer::brewer.pal(8, "Set1")[c(2:1, 3:5, 7:9)], fill=RColorBrewer::brewer.pal(8, "Set1")[c(2:1, 3:5, 7:9)], region=RColorBrewer::brewer.pal(n=11, name="Spectral")))')

  doItAndPrint("")

  if (getRcmdr("use.markdown")){
      removeNullRmdBlocks()
      putRcmdr("startNewCommandBlock", TRUE)
      if (getRcmdr("rmd.generated")) {
          endRmdBlock()
          putRcmdr("rmd.generated", FALSE)
      }
      removeNullRmdBlocks()
  }
  if (getRcmdr("use.knitr")){
      removeNullRnwBlocks()
      putRcmdr("startNewKnitrCommandBlock", TRUE)
      if (getRcmdr("rnw.generated")) {
         endRnwBlock()
         putRcmdr("rnw.generated", FALSE)
      }
      removeNullRnwBlocks()
  }
}

