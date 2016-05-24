# Rattle: A GUI for Data Mining in R
#
# BASE FUNCTIONS
#
# Time-stamp: <2015-09-21 22:21:48 gjw>
#
# Copyright (c) 2009-2015 Togaware Pty Ltd
#
# This files is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.

# 120704 Avoid "no visible binding for global variable" warnings on a
# check. However, this then requires R >= 2.15.1, so only do this
# conditionally, particularly that a lot of users are not in the
# upgrade habit as yet, and Revolution R is not up to 2.15 yet.

if(getRversion() >= "2.15.1")
  utils::globalVariables(c("rattle.entered.dataset",
                           "gladeXMLNew",
                           "gladeXMLSignalAutoconnect",
                           "biocLite",
                           "Caseload",
                           "Risk",
                           "Precision",
                           "pos",
                           "ticks",
                           "target",
                           "ignore",
                           "digit", "variable",
                           "split.labels",
                           "rbin",                    
                           "pacc",                    
                           "x",
                           "y",
                           "lbl",
                           "hj",
                           "vj",
                           "score",
                           "low",
                           "high"
                           ))

# The function paste0() was introduced in 2.15.0

if (! exists("paste0")) paste0 <- function(...) paste(..., sep="")

Rtxt <- function(...)
{
  # 100130 Currently, on Windows we are waiting for 2.12.17 of  RGtk2 with
  # rgtk2_bindtextdomain().

#  if (.Platform$OS.type == "windows")
#    paste(...)
#  else
    gettext(paste(...), domain="R-rattle")
}

# This is used to avoid the string being identified as a translation, as in
# RtxtNT(paste(vals ...))

RtxtNT <- Rtxt

VERSION <- "4.1.0"
DATE <- "2016-01-26"

# 091223 Rtxt does not work until the rattle GUI has started, perhaps?
COPYRIGHT <- paste(Rtxt("Copyright"), "(C) 2006-2015 Togaware Pty Ltd.")

# Acknowledgements: Frank Lu has provided much feedback and has
# extensively tested early versions of Rattle. Many colleagues at the
# Australian Taxation Office have used Rattle and made many and
# varied suggestions. These include Anthony Nolan, Stuart Hamilton,
# Liyin Zue, Weiqiang Lin, Robert Williams, Shawn Wicks, Ray Lindsay.

# LICENSE
#
# This files is part of Rattle.
#
# Rattle is open source software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.

# STYLE GUIDE
#
#    Use the "_" convention only for Glade variables and functions.
#    Use capitalised verbs for own functions: displayPlotAgain
#    Use dot separated words for variables: list.of.frames, lib.cmd
#    RGtk2 uses the capitalised word convention.
#    Use same names in R code as for the Glade objects.
#    Hide global variables, all capitalised, by including in crv$

# INTERFACE STYLE
#
#    080427 For options like a button to display a model once it has been
#    built or which model builders are available given the nature of
#    the data, we generally toggle the Sensistivity of the widgets
#    appropraitely.
#
#    In general, show all available widgets at any time, but grey out
#    those that can not yet be used because, for example, a model has
#    not yet been built.
#
#    If the functionality is not yet implemented, full stop, then have
#    the interface item(s) not present. This is better than having
#    them greyed out as the expectation is that perhaps there is some
#    way within the interface of getting it not to be greyed out! But
#    displaying future functionality also encourages those with an
#    interest in the greyed out bits to either complain (i.e., I get
#    to know what is wanted) or else help implement them!
#
#    If the functionality is not appropriate in a particular
#    circumstance then don't provide the button. Simply check, in the
#    one place in the code (e.g., when the button is pushed) and pop
#    up an error dialogue.
#
#    This doesn't always work, as in the case of sample where you do
#    want greyed out functionality, but you don't want it to mean not
#    yet implemented.

# BUGS
#
#   Tooltips used to have issues on GNU/Linux. Just fine on
#   MS/Windows.
#
#   The RGtk2 author, Michael Lawrence, notes that most of the GUI
#   functionality in Gnome (i.e., libgnome and libgnomeui) will soon
#   be merged into GTK. At that time, that functionality will be part
#   of RGtk2.

# GLOBALS
#
#   Original design placed state variables into the crs list and
#   global constants into . variables that then moved into the crv
#   list instead, after R CMD check started complaining about possibly
#   unbound variables. The real solution seems to be
#   environments. This was implemented temporarily simply by replacing
#   crv and crs with environments. The list notation then continued to
#   work for them! 090316 Finally removed all <<- assignments into the
#   environments, since, as Chambers (2008) page 124 points out a
#   reference to the environemt ralways refers to the same
#   environment.
#
#   Be aware that the trick of doing
#
# 	crs <- crs
#
#   within functions only works if we <<- assign to crs and don't make
#   use of the value in crs after it might change within the function
#   (or a sub function)! Probably not a good thing to do.

########################################################################
#
# INITIALISATIONS

## overwritePackageFunction <- function(fname, fun, pkg)
## {
##   # 090207 This allows a plugin to easily overwrite any Rattle funtion
##   # with their own functionality. Simply define your own FUN that is
##   # to overwrite the Rattle defined function FNAME. 090517 We do it
##   # this way rather than having to export the function to be
##   # overridden. Note that the override only happens within the
##   # namespace of the package. Thus it does not make sense to use this
##   # overwrite function to overwrite an exported function, since the
##   # overwrite will not be seen externally to the package. 120117
##   # Remove this for now since it could be harmful. Kurt has suggested
##   # only allowing overwriting when 're' is asNamespace('rattle') to
##   # reduce risk of malicious use by other packages.

##   re <- eval(parse(text=sprintf("environment(%s)", pkg)))
##   if (re == asNamespace('rattle')) # NOT RIGHT
##   {
##     unlockBinding(fname, re)
##     assign(fname, fun, re)
##     lockBinding(fname, re)
##   }
## }

toga <- function() browseURL("http://rattle.togaware.com")

########################################################################
# RATTLE Version 2

rattle <- function(csvname=NULL, dataset=NULL, useGtkBuilder=NULL)
{
  # 101113 Add the useGtkBuilder argument so that a user can override
  # the automatic determination of which one to use: libglade versus
  # GtkBuilder. If NULL then automatically determine.
  
  # 090517 Require pmml. Now that there is an indication on the Data
  # tab as to whether the varaiable (i.e., a transformed variable) can
  # be exported to PMML we need pmml to be loaded. Thus pmml is now a
  # "Depends:" in the DESCRIPTION file.

  # If crv$tooltiphack is TRUE then gtkMain is called on focus,
  # blocking the R console, but at least tooltips work. On losing
  # focus gtkMainQuit is called, and thus the console is no longer
  # blocked!  A bit ugly, but seems to work. This was suggested by
  # Felix Andrew, 080705. I notice that to load the supplied audit
  # dataset I need to change focus out of Rattle.

  # 080906 If crv$close="quit" then when the window close is pressed, we
  # also quit R.

  # 080319 Create global crv and crs to avoid many "no visible
  # binding" messages from "R CMD check" by adding all hidden
  # variables to crs and crv. Previously they all began with "." as in
  # crv$ADA used to be .ADA. "R CMD check" complained a lot, once for
  # each of these, so putting them all into crv means only one
  # complaint each time! Then defining crv in .onLoad removes the
  # NOTE altogether.

  # 090303 Make sure crv has been defined. This was necessitated
  # because CHECK does not run .onLoad in checking.

  if (! exists("crv"))
  {
    .onLoad()
    .onAttach()
  }

  # 090309 Reset the environment, crs, which stores the curret Rattle
  # state and used extensively throughout Rattle as a global
  # state. Not ideal for functional programming and only a hopefully
  # small deviation from Chamber's (2008) Prime Directive principle,
  # and similar to the "option" exception to the Prime Directive!

  # crs <<- new.env()
  sapply(ls(crs), function(x) assign(x, NULL, envir=crs))
  
  # crv$tooltiphack <<- tooltiphack # Record the value globally

  # 090525 Move to having the Setting option work on Linux. Thus
  # remove all this tooltip stuff.

  # if (crv$tooltiphack) crv$load.tooltips <- TRUE

  crv$.gtkMain <- FALSE # Initially gtkMain is not running.

  # 150712 No longer required as the package DEPENDS on RGtk2 now
  #if (packageIsAvailable("RGtk2", Rtxt("display the Rattle GUI")))
  #  suppressPackageStartupMessages(library("RGtk2", quietly=TRUE))
  #else
  #  stop(sprintf(Rtxt("The RGtk2 package is not available but is required",
  #                    "for the %s GUI."), crv$appname))

  # 101113 Use GtkBUilder or LibGlade?
  
  # 101009 We need to handle the case of an old install of Gtk (e.g.,
  # 2.12.9 on MS/Windows or GNU/Linux) where GtkBuilder does not
  # recognise the 'requires' element. We construct a string for the
  # xml and try to test this situation, and if the result from
  # gtkBuilderAddFromString has $error$message of "Unhandled tag:
  # 'requires'" then set crv$useGtkBuilder to FALSE.

  if (missing(useGtkBuilder))
  {
    op <- options(warn=-1)
    g <- RGtk2::gtkBuilderNew()
    res <- g$addFromString('<requires/>', 20)
    options(op)

    if (! res$retval && res$error$message[1] == "Unhandled tag: 'requires'")
      crv$useGtkBuilder <- FALSE
    else if (.Platform$OS.type=="windows" && version$major<="2" && version$minor<"12")
      # 101009 Always use glade for old installs of R on MS/Windows
      # rather than trying to figure out when it might work with
      # GtkBuilder.
      crv$useGtkBuilder <- FALSE
    else
      crv$useGtkBuilder <- TRUE
  }
  else
    crv$useGtkBuilder <- useGtkBuilder
  
  # Check to make sure libglade is available.

  if (! crv$useGtkBuilder)
    if (! exists("gladeXMLNew"))
      stop(Rtxt("The RGtk2 package did not find libglade installed.",
                "Please install it."))

  # On the Macintosh (and when using GtkBuilder 100821) we seem to
  # need to initialise all of the types for the GTK widgets. So do
  # that here.

  # 101127 No longer needed if (crv$useGtkBuilder || Sys.info()["sysname"] == "Darwin")
  # 111203 Is this still needed????? Try removing it.
  # 130412 Remove for now????
  if (isMac())
    fixMacAndGtkBuilderTypes()
 
  # Ensure the About dialog will respond to the Quit button.

  #on_aboutdialog_response <<- gtkWidgetDestroy

  # When an error is reported to the R Console, include a time
  # stamp. 120122 Remove the error timestamp for now. The message
  # remains after Rattle and users then attribute errors to Rattle.

#  options(error=function()
#          cat(sprintf(Rtxt("%s timestamp (for the message above):"), crv$appname),
#              sprintf("%s\n%s\n", Sys.time(), paste(rep("^", 72), collapse=""))))

  # Keep the loading of Hmisc quiet.

  options(Hverbose=FALSE)

  # Try firstly to load the glade file from the installed rattle
  # package, if it exists. Otherwise, look locally.

  if (crv$useGtkBuilder)
  {
    crv$rattleGUI <- RGtk2::gtkBuilderNew()
    crv$rattleGUI$setTranslationDomain("R-rattle")
  }
  
  result <- try(etc <- file.path(path.package(package="rattle")[1], "etc"),
                silent=TRUE)
  if (inherits(result, "try-error"))
    if (crv$useGtkBuilder)
      crv$rattleGUI$addFromFile(crv$rattleUI)
    else
      crv$rattleGUI <- gladeXMLNew("rattle.glade",
                                root="rattle_window", domain="R-rattle")
  else
    if (crv$useGtkBuilder)
      crv$rattleGUI$addFromFile(file.path(etc, crv$rattleUI))
    else
      crv$rattleGUI <- gladeXMLNew(file.path(etc,"rattle.glade"),
                                root="rattle_window", domain="R-rattle")

  if (crv$useGtkBuilder)

    # 101009 This sometimes gives an error on older GNU/Linux,
    # complaining that the element "require" is an unhandled tag. I
    # should be able to test this programatically in .onAttach and
    # then set crv$useGtkBuilder to FALSE in that case so we don't get
    # here.
    
    crv$rattleGUI$getObject("rattle_window")$show()
  
  # Really need an second untouched crv$rattleGUI

  #121212 DO WE NEED THIS NOW? Global_.rattleGUI <-crv$rattleGUI

  set.cursor("watch")
  on.exit(set.cursor())

  # 100821 As a temporary fix whilst Michael Lawrence gets theses
  # fixed.

  if (crv$useGtkBuilder) fixGtkBuilderAdjustments()
  
  # 090206 Tune the interface to suit needs, and in particular allow
  # packages to overwrite these functions so that the interface can be
  # tuned to suit plugins.

  setMainTitle()
  configureGUI()
  setDefaultsGUI()
  # 101008 Show toolbar text under the icons, if option is set.
  if (crv$toolbar.text) theWidget("toolbar")$setStyle("GTK_TOOLBAR_BOTH")
  
  # 100120 A temporary fix for MS/Windows where translations of stock
  # items by RGtk2 don't seem to be happening. It works just fine for
  # GNU/Linux. We probably only want to do this if we have a foreign
  # locale.
  
  if (isWindows())
  {
    fixTranslations()
    translateComboBoxes()
    translateMenus()
  }

  if (crv$load.tooltips) loadTooltips()

  # 120121 Don't show timestamps any more.
#  if (not.null(crv$show.timestamp) && crv$show.timestamp)
#    cat(crv$appname, "timestamp:",
#        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

  # 090708 Set the icon for the current window, and then make it the
  # default for all other windows. We do it here rather than earlier
  # in case configureGUI is overriddent to not change the icon.

  theWidget("rattle_window")$setIcon(crv$icon)
  if (! is.null(crv$icon)) RGtk2::gtkWindowSetDefaultIcon(crv$icon)

  # 080511 Record the current options and set the scientific penalty
  # to be 5 so we generally get numerics pinted using fixed rather
  # than exponential notation. We reset all options to what they were
  # at the startup of Rattle on closing Rattle. Not necessarily a good
  # idea since the knowing user may actually also change options
  # whilst Rattle is running.

  crv$options <- options(scipen=5)

  # 080924 Load of a supplied data file occurs here, but may take time
  # and whilst the UI is not fully set up yet, we see the Welcome
  # screen in Rattle displayed in plugins for 30 seconds or so. So
  # perhaps move it to later in the process.

  # Load data from the file identified by the csvname supplied in the
  # call to Rattle, or from the environment variable RATTLE_DATA if
  # defined, or from the variable RATTLE.DATA (as might be defined in
  # a .Rattle file), or else don't load any data by default.

  # First, always execute any .Rattle file in the current working
  # directory.

  # When reading the .Rattle file and identifying a dataset to load,
  # for some reason the stats package will not have been loaded at
  # this stage. The symptom is that median is not defined. So make
  # sure it is always available.

  # suppressPackageStartupMessages(library(stats, quietly=TRUE))

  if (file.exists(".Rattle")) source(".Rattle")

  # 110531 Some error checking first.

  if (!missing(csvname) && !missing(dataset))
  {
    errorDialog(Rtxt("Both a csvname= and a dataset= were specified in",
                     "the call to rattle(). At most, only one is allowed.",
                     "We will continue as if neither were specified."))
    csvname <- NULL
    dataset <- NULL
  }

  if (!missing(dataset) && class(dataset) != "character")
  {
    errorDialog(Rtxt("An actual dataset rather than the name of the",
                     "dataset has been supplied as the argument to",
                     "dataset=. Please supply the dataset name as a",
                     "character string.",
                     "We will continue as if no dataset was specified.",
                     "You can load it through the Data tab."))
    dataset <- NULL
  }
  
  if (is.null(csvname) && is.null(dataset))
  {
    # Use the .Rattle settings first, but these might be overriden if
    # the environment variable is defined.

    if (exists("RATTLE.DATA") && ! is.null(RATTLE.DATA)) csvname <- RATTLE.DATA

    # Obtain the value of the RATTLE_DATA environment variable and if
    # it is defined then use that at the csvname.

    if ((.rattle.data <- Sys.getenv("RATTLE_DATA")) != "")
      csvname <- .rattle.data
  }

  # Tidy up the csvname. TODO Is there an R command to do this, or
  # else put this into a function as I want to do it in a couple of
  # places (like further below in using RATTLE.SCORE.IN).

  if (not.null(csvname) && substr(csvname, 1, 4) == "http")
  {
    errorDialog(sprintf(Rtxt("URLs for the CSV filename are not currently supported.",
                             "\n\nWe found %s.",
                             "\n\nWe will continue but you will need to choose a",
                             "data file to load using the Filename button."),
                        csvname))
    csvname <- NULL
  }

  if (not.null(csvname))
  {
    csvname <- path.expand(csvname)

    # If it does not look like an absolute path then add in the
    # current location to make it absolute.

    if (substr(csvname, 1, 1) %notin% c("\\", "/")
        && substr(csvname, 2, 2) != ":")
      csvname <- file.path(getwd(), csvname)
    if (! file.exists(csvname))
    {
      infoDialog(sprintf(Rtxt("The supplied CSV filename '%s' does not exist."),
                         csvname))
      csvname <- NULL
    }
    else
    {
      # 081020 gjw If the csvname is supplied then prefix it with
      # file:/// to make it conform to the filename obtained from the
      # file chooser button. Without doing this crs$dwd does not
      # include file:/// and when compared in dataNeedsLoading to the
      # filename obtained with getUri they don't match, and hence the
      # data is reloaded! Take care of MS/Windows where the csvname
      # will be prefixed by the drive, so we add three slashes in
      # front.

      if (isWindows())
        csvname <- paste0("file:///", csvname)
      else
        csvname <- paste0("file://", csvname)
    }
  }

  # Keep the loading of Hmisc quiet.

  options(Hverbose=FALSE)

  # Load the Rattle GUI specification. The three commands here
  # represent an attempt to be independent of where R is running and
  # where rattle.R is located by finding out from the system calls the
  # actual call to source rattle.R, and then point to this same
  # location for finding rattle.glade. Assumes the call to source is
  # something like: source("abc/def/rattle.R"). The better alternative
  # might be to tell people to use the chdir=TRUE option in source.

  ##s <- as.character(sys.calls())
  ##n <- grep("source", s)
  ##p <- gsub("\.R..$", ".glade", gsub("source..", "", s[n]))

  # Constants: I would like these available within this package, but
  # not outside? Do I use assign in some way? That is, how to keep
  # these constants within the package only.

  # TODO Put these constants into the top level of this file, defined
  # as NULL. Then keep these double arrow assignments here. I think
  # then that they will stay with the package, but not be in
  # .GlobalEnv because the package scope will be found before the top
  # level.

  ########################################################################
  # PACKAGE GLOBAL CONSTANTS
  #
  # These are double arrow assigned here to place them in
  # .GlobalEnv. I couldn't figure out an easy way to keep them scoped
  # locally. TODO Needs cleaning up.
  #
  # Various Treeview Columns

  crv$COLUMN <- c(number = 0, variable = 1, type = 2, input = 3,
                  target = 4, risk = 5, ident = 6, ignore = 7,
                  weight = 8, comment = 9)
  crv$COLUMNstart <- crv$COLUMN[["input"]]
  crv$COLUMNend <- crv$COLUMN[["weight"]]
  
  crv$IMPUTE <- c(number=0, variable=1, comment=2)

  crv$CATEGORICAL <- c(number = 0, variable = 1, barplot = 2,
                       dotplot = 3, mosplot = 4, paiplot=5, comment = 6)

  crv$CONTINUOUS <-  c(number = 0, variable = 1, boxplot = 2,
                       hisplot = 3, cumplot = 4, benplot = 5, paiplot=6, comment = 7)

  # Create constants naming DESCRIBE (i.e., the descriptive model
  # builders) and PREDICT (i.e., the predictive model builders). Note
  # that these are migrating into the crv variable, but not all are
  # done yet. 081227 Also note that kmeans, hclust and apriori will
  # also be migrating into being treated as first class models.

  crv$KMEANS 	<- "kmeans"
  crv$EWKM 	<- "ewkm"
  crv$CLARA 	<- "clara"
  crv$PAM  	<- "pam"
  crv$HCLUST 	<- "hclust"
  crv$BICLUST 	<- "biclust"
  crv$APRIORI 	<- "apriori"

  # 091218 Not yet - avoid issues with RStat release.
  # crv$DESCRIBE <- c(crv$KMEANS, crv$CLARA, crv$PAM, crv$HCLUST, crv$BICLUST, crv$APRIORI)
  crv$DESCRIBE <- c(crv$KMEANS, crv$HCLUST, crv$APRIORI)

  crv$GLM   	<- "glm"
  crv$RPART 	<- "rpart"
  #GBM <- "gbm"
  crv$ADA   	<- "ada"
  crv$RF    	<- "rf"
  crv$SVM   	<- "svm"
  crv$KSVM  	<- "ksvm"
  crv$NNET  	<- "nnet"
  crv$SURVIVAL <- "survival"

  crv$PREDICT <- c(crv$RPART, crv$ADA, crv$RF, crv$KSVM, crv$GLM,
                     crv$NNET, crv$SURVIVAL)

  # PACKAGE STATE VARIABLE

  # 090309 The following is now taken care of in .onLoad as defined in
  # zzz.R.

  ## if (TRUE)
  ##   crs <<- new.env()
  ## else
  ##   crs <<- list(dataset=NULL,
  ##              dataname=NULL,
  ##              dwd=NULL, 	# Data Working Directory
  ##              mtime=NULL,	# Modification time of file
  ##              pwd=NULL,	# Project Working Directory
  ##              input=NULL,
  ##              target=NULL,
  ##              weights=NULL,
  ##              risk=NULL,
  ##              ident=NULL,
  ##              ignore=NULL,
  ##              nontargets=NULL, # 080426 Started but not yet implemented
  ##              sample=NULL,
  ##              sample.seed=NULL,
  ##              kmeans=NULL,
  ##              kmeans.seed=NULL,
  ##              hclust=NULL,
  ##              page="",
  ##              smodel=NULL, # Record whether the sample has been modelled
  ##              glm=NULL,
  ##              rpart=NULL,
  ##              ada=NULL,
  ##              apriori=NULL,
  ##              rf=NULL,
  ##              svm=NULL,
  ##              ksvm=NULL,
  ##              perf=NULL,
  ##              eval=NULL,
  ##              testset=NULL,
  ##              testname=NULL,
  ##              alog=NULL,	# Record of interaction - useful?
  ##              transforms=NULL  # Record of variable transforms for inclusion in PMML
  ##              )

  # Main notebook related constants and widgets.  Track the widgets
  # that are needed for removing and inserting tabs in the notebook,
  # depending on the selected paradigm. TODO Paradigms have gone as of
  # 080519 so we may not need all this machinery now!

  crv$NOTEBOOK <- theWidget("notebook")

  # 100122 The Rtxt is required for these since Glade will translate
  # these labels. These labels are for tabs that are visible in the
  # GUI.

  crv$NOTEBOOK.DATA.NAME <- Rtxt("Data")

  crv$NOTEBOOK.TEST.NAME <- Rtxt("Test")

  crv$NOTEBOOK.EXPLORE.NAME <- Rtxt("Explore")

  crv$NOTEBOOK.TRANSFORM.NAME <- Rtxt("Transform")

  crv$NOTEBOOK.CLUSTER.NAME    <- Rtxt("Cluster")
  crv$NOTEBOOK.CLUSTER.WIDGET <- theWidget("cluster_tab_widget")
  crv$NOTEBOOK.CLUSTER.LABEL  <- theWidget("cluster_tab_label")

  crv$NOTEBOOK.ASSOCIATE.NAME    <- Rtxt("Associate")
  crv$NOTEBOOK.ASSOCIATE.WIDGET <- theWidget("associate_tab_widget")
  crv$NOTEBOOK.ASSOCIATE.LABEL  <- theWidget("associate_tab_label")

  # 100716 Revert to using Model rather than Predictive.... Model fits
  # the other tabs better.
  crv$NOTEBOOK.MODEL.NAME    <- Rtxt("Model")
  # crv$NOTEBOOK.MODEL.NAME <- theWidget("model_tab_label")$getLabel()
#  if (is.null(crv$NOTEBOOK.MODEL.NAME)) # 100423 Fix for RStat using Model
#    crv$NOTEBOOK.MODEL.NAME <- Rtxt("Predictive")
  crv$NOTEBOOK.MODEL.WIDGET  <- theWidget("model_tab_widget")
  crv$NOTEBOOK.MODEL.LABEL   <- theWidget("model_tab_label")

  crv$NOTEBOOK.EVALUATE.NAME    <- Rtxt("Evaluate")
  crv$NOTEBOOK.EVALUATE.WIDGET <- theWidget("evaluate_tab_widget")
  crv$NOTEBOOK.EVALUATE.LABEL  <- theWidget("evaluate_tab_label")

  crv$NOTEBOOK.LOG.NAME       <- Rtxt("Log")

  # 100122 Every call to getNotebookPage used to need the second
  # argument wrapped with an Rtxt. Glade translates these on
  # loading. 100416 But that was causing issues. Let's instead ensure
  # these tabs, which are never visible, remain in English, not
  # translated, and we use them directly as English throughout the
  # Rattle code.

  # 080921 Define the DATA tab pages

  crv$DATA.NOTEBOOK 	<- theWidget("data_notebook")
  crv$DATA.CORPUS.TAB      <- getNotebookPage(crv$DATA.NOTEBOOK, "corpus")
  crv$DATA.CSV.TAB         <- getNotebookPage(crv$DATA.NOTEBOOK, "csv")

  crv$DATA.DISPLAY.NOTEBOOK     <- theWidget("data_display_notebook")
  crv$DATA.DISPLAY.TREEVIEW.TAB <- getNotebookPage(crv$DATA.DISPLAY.NOTEBOOK,
                                                   "treeview")
  crv$DATA.DISPLAY.WELCOME.TAB  <- getNotebookPage(crv$DATA.DISPLAY.NOTEBOOK,
                                                   "welcome")
  if (isJapanese())
  {
    # 100227 For some reason the following is not working properly:
    #   nb <- rattle:::theWidget("notebook")
    #   nb$getTabLabelText(nb$getNthPage(0))
    # The result should be the same as
    #   rattle:::Rtxt("Data")
    # It appears the UTF is being interpreted as Shift-JIS
    # So hardcode these (perhaps a growing list)
    
    crv$DATA.DISPLAY.TREEVIEW.TAB <- 0
    crv$DATA.DISPLAY.WELCOME.TAB  <- 1
  }
  

  # Define the TRANSFORM tab pages

  crv$TRANSFORM               <- theWidget("transform_notebook")
  # TODO 080423 Change to RESCALE
  crv$TRANSFORM.NORMALISE.TAB <- getNotebookPage(crv$TRANSFORM, "normalise")
  crv$TRANSFORM.IMPUTE.TAB    <- getNotebookPage(crv$TRANSFORM, "impute")
  crv$TRANSFORM.REMAP.TAB     <- getNotebookPage(crv$TRANSFORM, "remap")
  crv$TRANSFORM.OUTLIER.TAB   <- getNotebookPage(crv$TRANSFORM, "outlier")
  crv$TRANSFORM.CLEANUP.TAB   <- getNotebookPage(crv$TRANSFORM, "cleanup")

  crv$EXPLORE                 <- theWidget("explore_notebook")
  crv$EXPLORE.SUMMARY.TAB     <- getNotebookPage(crv$EXPLORE, "summary")
  crv$EXPLORE.PLOT.TAB        <- getNotebookPage(crv$EXPLORE, "explot")
  crv$EXPLORE.CORRELATION.TAB <- getNotebookPage(crv$EXPLORE, "correlation")
  crv$EXPLORE.PRCOMP.TAB      <- getNotebookPage(crv$EXPLORE, "prcomp")
  crv$EXPLORE.INTERACTIVE.TAB <- getNotebookPage(crv$EXPLORE, "interactive")

  crv$CLUSTER             <- theWidget("cluster_notebook")
  crv$CLUSTER.KMEANS.TAB  <- getNotebookPage(crv$CLUSTER, "kmeans")
  crv$CLUSTER.EWKM.TAB    <- getNotebookPage(crv$CLUSTER, "ewkm")
  crv$CLUSTER.CLARA.TAB   <- getNotebookPage(crv$CLUSTER, "clara")
  crv$CLUSTER.PAM.TAB     <- getNotebookPage(crv$CLUSTER, "pam")
  crv$CLUSTER.HCLUST.TAB  <- getNotebookPage(crv$CLUSTER, "hclust")
  crv$CLUSTER.BICLUST.TAB <- getNotebookPage(crv$CLUSTER, "biclust")

  crv$MODEL           <- theWidget("model_notebook")
  crv$MODEL.RPART.TAB <- getNotebookPage(crv$MODEL, "rpart")
  crv$MODEL.GLM.TAB   <- getNotebookPage(crv$MODEL, "glm")
  crv$MODEL.ADA.TAB   <- getNotebookPage(crv$MODEL, "ada")
  ## crv$MODEL.GBM.TAB   <- getNotebookPage(crv$MODEL, .GBM)
  crv$MODEL.RF.TAB    <- getNotebookPage(crv$MODEL, "rf")
  crv$MODEL.SVM.TAB   <- getNotebookPage(crv$MODEL, "svm")
  crv$MODEL.NNET.TAB   <- getNotebookPage(crv$MODEL, "nnet")
  crv$MODEL.SURVIVAL.TAB <- getNotebookPage(crv$MODEL, "survival")

  crv$SVMNB           <- theWidget("svm_notebook")
  crv$SVMNB.ESVM.TAB  <- getNotebookPage(crv$SVMNB, "esvm")
  crv$SVMNB.KSVM.TAB  <- getNotebookPage(crv$SVMNB, "ksvm")

  crv$EVALUATE                 <- theWidget("evaluate_notebook")
  crv$EVALUATE.CONFUSION.TAB   <- getNotebookPage(crv$EVALUATE, "confusion")
  crv$EVALUATE.RISK.TAB        <- getNotebookPage(crv$EVALUATE, "risk")
  crv$EVALUATE.LIFT.TAB        <- getNotebookPage(crv$EVALUATE, "lift")
  crv$EVALUATE.ROC.TAB         <- getNotebookPage(crv$EVALUATE, "roc")
  crv$EVALUATE.PRECISION.TAB   <- getNotebookPage(crv$EVALUATE, "precision")
  crv$EVALUATE.SENSITIVITY.TAB <- getNotebookPage(crv$EVALUATE, "sensitivity")
  crv$EVALUATE.COSTCURVE.TAB   <- getNotebookPage(crv$EVALUATE, "costcurve")
  crv$EVALUATE.PVO.TAB         <- getNotebookPage(crv$EVALUATE, "pvo")
  crv$EVALUATE.SCORE.TAB       <- getNotebookPage(crv$EVALUATE, "score")

  # Turn off the sub-notebook tabs.

  # Sys.sleep(5) 080924 to test delays....

  crv$DATA.NOTEBOOK$setShowTabs(FALSE)
  crv$DATA.DISPLAY.NOTEBOOK$setShowTabs(FALSE)
  crv$EXPLORE$setShowTabs(FALSE)
  crv$TRANSFORM$setShowTabs(FALSE)
  crv$CLUSTER$setShowTabs(FALSE)
  crv$MODEL$setShowTabs(FALSE)
  crv$EVALUATE$setShowTabs(FALSE)

  ########################################################################
  # Connect the callbacks.

  if (crv$useGtkBuilder)
    crv$rattleGUI$connectSignals()
  else
    gladeXMLSignalAutoconnect(crv$rattleGUI)

  # Enable the tooltips Settings option on GNU/Linux. Under MS/Windows
  # tooltips have always worked so this option is not relevant. 110409
  # Tooltips seem to be on by default, even on GNU/Linux now, so I
  # changed the FALSE to TRUE here to reflect that. However, it seems
  # that we can't actually turn tooltips off from the Settings menu.

  if (isLinux() && crv$load.tooltips)
  {
    theWidget("tooltips_menuitem")$show()
    theWidget("tooltips_menuitem")$setActive(TRUE)
  }

  ########################################################################
  # User interface initialisations.

  initialiseVariableViews()

  # Ensure the filechooserbutton by default will filter CSVs.

  updateFilenameFilters("data_filechooserbutton", "CSV")

  # Do not enable ARFF option for versions before 2.5.0 where it was
  # not included in the foreign package.

  if (!exists("getRversion", baseenv()) || getRversion() <= "2.4.0")
    theWidget("arff_radiobutton")$hide()

  theWidget("model_tree_include_missing_checkbutton")$setActive(FALSE)
  #theWidget("glm_family_comboboxentry")$setActive(0)
  theWidget("svm_kernel_combobox")$setActive(0)

  ## Check if some external applications are available and if not
  ## de-sensitise their functionality.

  # How to test if ggobi is actually available?

  # If the cairoDevice package is not available then turn off the
  # option in the settings menu and make it insensitive.

  # 100706 The asCairo is failing:
  # Error in asCairoDevice(da) : Graphics API version mismatch
  # 111111 This usually can be solved with a reinstall of the package:
  # > install.packages("cairoDevice")
  
  if (! packageIsAvailable("cairoDevice", Rtxt("enable the cairo device option")))
  {
    theWidget("use_cairo_graphics_device")$setActive(FALSE)
    theWidget("use_cairo_graphics_device")$hide()
  }

  # 110810 On MS/Windows the CairoDevice seems to drop some graphics
  # elements whe ndrawing multiple plots, so by default, on Windows,
  # turn it off for now. See further comments in newPlot(). The
  # problem is exhibited in Figs 2.8 and 15.5 of the Rattle book.

  if (isWindows())
    theWidget("use_cairo_graphics_device")$setActive(FALSE)

  # Tell MS/Windows to use 2GB (TODO - What's needed under Win64?)
  #
  # Brian D. Ripley 15 Jul 2007 07:57:49 +0100 requested the memory mod
  # be removed:
  #
  # First, because you should not be setting the limit high if the
  # user has only limited memory: the defaults are chosen with a lot
  # of care.  Second, because the default can be as high as 2.5Gb on a
  # machine with 4Gb RAM and the /3GB switch set (the case here).
  #
  # The correct way to refer to things in packages on which you have
  # not declared a dependence is utils::memory.limit.

  # if (isWindows()) utils::memory.limit(2073)

  ## By default the CLUSTER page is not showing.

  ## Don't turn this off until we move away from using the numeric tab
  ## variables above, since a Execute on the Model tab runs the
  ## Cluster tab :-)

##  crv$NOTEBOOK$removePage(getNotebookPage(crv$NOTEBOOK, crv$NOTEBOOK.CLUSTER.NAME))
##  crv$NOTEBOOK$removePage(getNotebookPage(crv$NOTEBOOK, crv$NOTEBOOK.ASSOCIATE.NAME))

##  while (gtkEventsPending()) gtkMainIteration() # Make sure window is displayed

   # Tooltips work when gtkMain is called, but the console is blocked
   # and need gtkMainQuit.

  # if (tooltiphack) gtkMain()

  # TODO Add a console into Rattle to interact with R.

  # 080510 Display a relevant welcome message in the textview.

  displayWelcomeTabMessage()

  initiateLog()

  # Make sure the text is shown on startup.

  while (RGtk2::gtkEventsPending()) RGtk2::gtkMainIterationDo(blocking=FALSE)

  # Now deal with any arguments to rattle.

  if (not.null(dataset))
  {
    theWidget("data_rdataset_radiobutton")$setActive(TRUE)
    
    # 110531 TODO Get list of available data frames from the combobox,
    # choose the right one, and then Execute. How to get the list of
    # current values in the combobox? Instead, for now do the same
    # listing of the data frames, and assume we get the same
    # list. TODO This takes some time and so not really the right
    # thing to do.

    dl <- unlist(sapply(ls(sys.frame(0)),
                        function(x)
                        {
                          cmd <- sprintf(paste("is.data.frame(%s) ||",
                                               'inherits(%s,',
                                               '"sqlite.data.frame")'), x, x)
                          var <- try(ifelse(eval(parse(text=cmd), sys.frame(0)),
                                            x, NULL), silent=TRUE)
                          if (inherits(var, "try-error"))
                            var <- NULL
                          return(var)
                        }))
    theWidget("data_name_combobox")$setActive(which(dataset == dl)[1]-1)
    # Make sure GUI updates
    while (RGtk2::gtkEventsPending()) RGtk2::gtkMainIterationDo(blocking=FALSE)
    executeDataTab()
  }      
  else if (not.null(csvname))
  {
    if (!theWidget("data_filechooserbutton")$setUri(csvname))
      infoDialog(Rtxt("The setting of the filename box failed."), crv$support.msg)
    # Make sure GUI updates
    while (RGtk2::gtkEventsPending()) RGtk2::gtkMainIterationDo(blocking=FALSE)
    executeDataTab(csvname)
  }

  ## theWidget("csv_filechooserbutton")$setFilename("audi.csv")

  # Call resetRattle to ensure all textviews get their default texts

  resetRattle(FALSE)

  invisible()
}

rattleReport <- function()
{
  result <- paste0("Rattle Report: Summary of the Current Model(s)\n\n",
                   "Date\n",
                   "\t", Sys.time(), "\n",
                   "Project Name\n",
                   "\t", crs$dataname, "\n",
                   "Data Miner\n",
                   "\t", Sys.info()["user"], "\n",
                   "\n",
                   "Input Variables\n",
                   paste0("\t", crs$input, collapse="\n"), "\n",
                   "Target Variable\n",
                   "\t", crs$target, "\n",
                   "Risk Variable\n",
                   "\t", crs$risk, "\n",
                   "Identifiers\n",
                   paste0("\t", crs$ident, collapse="\n"), "\n",
                   "Ignored Variables\n",
                   paste0("\t", crs$ignore, collapse="\n"), "\n",
                   "\n",
                   "Models\n",
                   if (not.null(crs$rpart))
                   {
                     "\tTree\n"
                     # Add in the function call
                   })
  return(result)
}

########################################################################
# Configurable functions - these are here because plugins may want to
# overwrite them.

configureGUI <- function()
{
  # Toolbar

  theWidget("report_toolbutton")$show()

  id.string <- paste0('<span foreground="blue">',
                      '<i>', crv$appname, '</i> ',
                      '<i>', Rtxt("Version"), ' ', VERSION, '</i> ',
                      # 100115 It was found that crv$version we not
                      # being updated so use VERSION instead. Not sure
                      # why.
                      '<i><span underline="single">togaware.com</span></i>',
                      '</span>')

  rattle.menu <- theWidget("rattle_menu")
  rattle.menu$SetRightJustified(TRUE)
  rattle.menu$getChild()$setMarkup(id.string)

  # Icon 090705 Set the icon to be the R logo. Save the pixbuf in
  # crv$icon so that plots can also set the icon appropriately. How to
  # get all windows to inherit this icon?

  crv$icon <- system.file("etc/Rlogo.png", package="rattle")
  if (crv$icon == "" && file.exists("./Rlogo.png"))
    crv$icon <- "./Rlogo.png"
  if (crv$icon == "")
    crv$icon <- NULL
  else
    crv$icon <- RGtk2::gdkPixbufNewFromFile(crv$icon)$retval

  # 150921 Change the Connect-R button to be the Connect-R logo.

  connectr.logo <- system.file("etc/ConnectRlogo.png", package="rattle")
  connectr.pixbuf <- RGtk2::gdkPixbufNewFromFile(connectr.logo)$retval
  connectr.icon <- RGtk2::gtkImageNewFromPixbuf(connectr.pixbuf)
  connectr.button <- theWidget("connectr_toolbutton")
  RGtk2::gtkToolButtonSetIconWidget(connectr.button, connectr.icon)
  
  # 101202 Remove the By Group button and instead if a rescale has a
  # categoric selected then do by group. TODO.
  
  # theWidget("normalise_interval_radiobutton")$hide()


  # 110911 Although this function is deprecated, it works to ensure
  # that a Maximize, Un-Maximize returns to the original
  # size. Otherwise it miscalculates that the minimum width is
  # actaully quite wide, and so we end up with a very wide window -
  # Ugly and also difficult to shrink it. We suppress warnings to
  # avoid seeing:
  #
  # Warning message:
  # 'method' is deprecated.
  # Use 'gtkWindowSetResizable' instead.
  # See help("Deprecated") and help("RGtk2-deprecated"). 
  #
  # setResizable(TRUE) is the default but we stillget this problem.
  
   suppressWarnings(crv$rattleGUI$getObject("rattle_window")$setPolicy(TRUE, TRUE, TRUE))

}

setDefaultsGUI <- function()
{
  # 100315 Handle CSV defaults typical in Europe, as suggested by
  # Denis Brion.
  
  decimal <- Sys.localeconv()["decimal_point"]
  if (decimal == ",")
  {
    theWidget("data_separator_entry")$setText(";")
    theWidget("data_decimal_entry")$setText(",")
  }
}  

fixMacAndGtkBuilderTypes <- function()
{
  # 100821 This is required for the max, and also for the GtkBuilder
  # part of Rgtk2 until Michael Lawrence gets a chance to fix it. The
  # GtkBulder stuff added 100821 on the move from libglade2. Note that
  # it may not be needed for Mac (Sys.info()["sysname"] == "Darwin")
  # any more.
  
  # Use the following to extract all widgets from the glade file:
  #
  # $ grep '<widget' rattle.glade | sed 's|^.*widget class="||' |\
  #   sed 's|".*$||' | sort -u | sed 's|^Gtk|gtk|' |\
  #   awk '{printf("%sGetType()\n", $1)}'

  RGtk2::gtkAboutDialogGetType()
  RGtk2::gtkAlignmentGetType()
  RGtk2::gtkButtonGetType()
  RGtk2::gtkCheckButtonGetType()
  RGtk2::gtkCheckMenuItemGetType()
  RGtk2::gtkComboBoxEntryGetType()
  RGtk2::gtkComboBoxGetType()
  RGtk2::gtkDrawingAreaGetType()
  RGtk2::gtkEntryGetType()
  RGtk2::gtkFileChooserButtonGetType()
  RGtk2::gtkFileChooserDialogGetType()
  RGtk2::gtkFrameGetType()
  RGtk2::gtkHBoxGetType()
  RGtk2::gtkHButtonBoxGetType()
  RGtk2::gtkHSeparatorGetType()
  RGtk2::gtkHandleBoxGetType()
  RGtk2::gtkImageGetType()
  RGtk2::gtkImageMenuItemGetType()
  RGtk2::gtkLabelGetType()
  RGtk2::gtkListStoreGetType()
  RGtk2::gtkMenuBarGetType()
  RGtk2::gtkMenuGetType()
  RGtk2::gtkMenuItemGetType()
  RGtk2::gtkMiscGetType()
  RGtk2::gtkNotebookGetType()
  RGtk2::gtkRadioButtonGetType()
  RGtk2::gtkScrolledWindowGetType()
  RGtk2::gtkSeparatorMenuItemGetType()
  RGtk2::gtkSeparatorToolItemGetType()
  RGtk2::gtkSpinButtonGetType()
  RGtk2::gtkStatusbarGetType()
  RGtk2::gtkTableGetType()
  RGtk2::gtkTextViewGetType()
  RGtk2::gtkToolButtonGetType()
  RGtk2::gtkToolItemGetType()
  RGtk2::gtkToolbarGetType()
  RGtk2::gtkTreeViewGetType()
  RGtk2::gtkVBoxGetType()
  RGtk2::gtkVSeparatorGetType()
  RGtk2::gtkWidgetGetType()
  RGtk2::gtkWindowGetType()
}

fixGtkBuilderAdjustments <- function()
{
  # 100821 As a temporary fix whilst Michael Lawrence gets theses
  # fixed.

  wid <- theWidget("sample_seed_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, NULL, 100000000, 1, 100, 0)
  wid$setAdjustment(nad)
  wid$setValue(42)
  
  wid <- theWidget("data_odbc_limit_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 100000000, 1, 100, 0)
  wid$setAdjustment(nad)
  wid$setValue(0)
  
  wid <- theWidget("sample_percentage_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 100, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(70)
  
  wid <- theWidget("sample_count_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 100000000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(0)
  
  wid <- theWidget("plots_per_page_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 1, 9, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(4)
  
  wid <- theWidget("benford_digits_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 1, 9, 1, 2, 0)
  wid$setAdjustment(nad)
  wid$setValue(1)
  
  wid <- theWidget("normalise_interval_numgroups_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 100000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(100)

  wid <- theWidget("remap_bins_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 100, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(4)
  
  wid <- theWidget("kmeans_clusters_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 2, 100000000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(10)
  
  wid <- theWidget("kmeans_seed_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 100000000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(42)
  
  wid <- theWidget("kmeans_runs_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 1, 1000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(1)
  
  wid <- theWidget("hclust_nbproc_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 1, 100, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(1)
  
  wid <- theWidget("hclust_clusters_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 2, 100000000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(10)
  
  wid <- theWidget("associate_support_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 1, 0.01, 0.1, 0)
  wid$setAdjustment(nad)
  wid$setValue(0.1)
  
  wid <- theWidget("associate_confidence_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 1, 0.01, 0.1, 0)
  wid$setAdjustment(nad)
  wid$setValue(0.1)
  
  wid <- theWidget("associate_lift_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 100, 0.1, 0.5, 0)
  wid$setAdjustment(nad)
  wid$setValue(1.5)
  
  wid <- theWidget("rpart_minsplit_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 100000000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(20)
  
  wid <- theWidget("rpart_minbucket_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 1, 1000000000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(7)
  
  wid <- theWidget("rpart_maxdepth_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 1, 30, 1, 5, 0)
  wid$setAdjustment(nad)
  wid$setValue(20)
  
  wid <- theWidget("model_tree_cp_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0.00001, 1, 0.0001, 0.001, 0)
  wid$setAdjustment(nad)
  wid$setValue(0.01)
  
  wid <- theWidget("ada_ntree_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 10000, 10, 50, 0)
  wid$setAdjustment(nad)
  wid$setValue(50)
  
  wid <- theWidget("ada_maxdepth_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 1, 30, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(30)
  
  wid <- theWidget("ada_minsplit_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 10000000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(20)
  
  wid <- theWidget("ada_cp_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, -1, 1, 0.00001, 0.001, 0)
  wid$setAdjustment(nad)
  wid$setValue(0.01)
  
  wid <- theWidget("ada_xval_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 100, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(10)
  
  wid <- theWidget("ada_draw_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 1000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(1)
  
  wid <- theWidget("rf_ntree_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 10000, 10, 50, 0)
  wid$setAdjustment(nad)
  wid$setValue(500)
  
  wid <- theWidget("rf_mtry_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 1, 1000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(10)
  
  wid <- theWidget("rf_print_tree_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 1000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(1)
  
  wid <- theWidget("svm_poly_degree_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 1, 10, 1, 2, 0)
  wid$setAdjustment(nad)
  wid$setValue(1)
  
  wid <- theWidget("nnet_hidden_nodes_spinbutton")
  nad <- RGtk2::gtkAdjustmentNew(NULL, 0, 10000, 1, 10, 0)
  wid$setAdjustment(nad)
  wid$setValue(10)
}  

fixTranslations <- function(w=theWidget("rattle_window"))
{
  # 101111 If the widget does not have a name, ignore it. This was
  # needed for MS/Windows R 2.12.0 for some reason. 101127 But now the
  # children widgets are not getting translated! I guess previously
  # getName() rutnerned an empty string, but is now returning NULL.
  
  ## if (! length(w$getName())) return()
  
  # Ignore these since they are already translated and we end up with
  # a corrupted string passing through to Rtxt again. generally they
  # are Stock Items.

  if (length(w$getName()) &&
      w$getName() %in% c("execute_button", "new_button", "open_button",
                         "save_button", "stop_button", "quit_button",
                         "data_filechooserbutton",
                         "continuous_clear_button", "categorical_clear_button",
                         "execute_menu",
                         "print_textview_menu", "about_menu",
                         "evaluate_filechooserbutton"))
    return()

  # 100410 The following should be translated, unless we are in RStat
  # where they are named Regression rather than Linear, or are not
  # used, or otherwise differently handled.
  
  if (crv$appname == "RStat" && length(w$getName()) && w$getName() %in%
      c("data_sample_checkbutton",
        "data_script_radiobutton",
        "model_linear_radiobutton",
        "glm_linear_radiobutton",
        "evaluate_glm_checkbutton"))
    return()
  
  if ("GtkLabel" %in% class(w))
    w$setLabel(Rtxt(w$getLabel()))
  else if ("GtkNotebook" %in% class(w))
    lapply(RGtk2::gtkChildren(w),
           function(wc)
             w$getTabLabel(wc)$setLabel(Rtxt(w$getTabLabelText(wc))))

  #  if ("GtkLabel" %in% class(w)) w$setLabel("Fred")
  if ("GtkContainer" %in% class(w))
    lapply(RGtk2::gtkChildren(w), fixTranslations)
  
  return()
}

translateMenus <- function()
{
  # 100328 The menus were not getting fixed, since we need to
  # specifically traverse them it seems.
  
  menus <- c("tools_menu", "settings_menu", "help_menu",
             "help_data_menu", "help_explore_menu", "help_test_menu",
             "help_transform_menu", "help_transform_rescale_menu",
             "help_transform_impute_menu", "help_transform_remap_menu",
             "help_transform_cleanup_menu", "help_cluster_menu", "help_model_menu",
             "help_evaluate_menu")
  sapply(sapply(menus, theWidget), fixTranslations)
}


translateComboBoxes <- function()
{
  # 100313 We do this in the code when we are running MS/Windows
  # because the list is not translated using GTK+.

  combos <- c("data_odbc_table_combobox",
              "explore_correlation_method_combobox",
              "svm_kernel_combobox", "hclust_distance_combobox",
              "hclust_link_combobox")
  
  printNode <- function(model, path, iter, data)
    {vals <<- c(vals, model$getValue(iter, 0)$value); integer(1)}
  for (cb in combos)
  {
    # Iterate over the current entries and get label, then set label
    # to Rtxt value.

    # Get the actual object.
    
    cbw <- theWidget(cb)

    # Retrieve the current entries for the combobox.
    
    vals <- NULL
    cbw$getModel()$foreach(printNode)

    # Clear the current entries
    
    cbw$getModel()$clear()

    # Add the translated entries. Note that for entries defined in
    # glade, the actual string that is translated is made up of all of
    # the entries concatenated, with "\n" separating them. So we need
    # to reconstructt this string, translate, then split, then
    # appendText for each one.
    
    sapply(strsplit(RtxtNT(paste(vals, collapse="\n")), "\n")[[1]], cbw$appendText)

    # Reset default choice. Assume to be 0.
    
    cbw$setActive(0)
  }
}

displayWelcomeTabMessage <- function()
{
  crv$DATA.DISPLAY.NOTEBOOK$setCurrentPage(crv$DATA.DISPLAY.WELCOME.TAB)
  resetTextview("rattle_welcome_textview",
                paste0(Rtxt("Welcome to Rattle (rattle.togaware.com)."),
                       "\n\n",
                       Rtxt("Rattle is a free graphical user",
                            "interface for Data Mining, developed using R.",
                            "R is a free software environment",
                            "for statistical computing and graphics.",
                            "Together they provide a sophisticated",
                            "environments for data mining,",
                            "statistical analyses, and data visualisation."),
                       "\n\n",
                       Rtxt("See the Help menu for extensive support in",
                            "using Rattle.",
                            "The book Data Mining with Rattle and R is available from",
                            "Amazon.",
                            "The Togaware Desktop Data Mining Survival Guide",
                            "includes Rattle documentation",
                            "and is available from",
                            "datamining.togaware.com"),
                       "\n\n",
                       Rtxt("Rattle is licensed under the",
                            "GNU General Public License, Version 2.",
                            "Rattle comes with ABSOLUTELY NO WARRANTY.",
                            "See Help -> About for details."),
                       "\n\n",
                       sprintf(Rtxt("Rattle Version %s.",
                                    "Copyright 2006-2015 Togaware Pty Ltd."),
                               VERSION),
#LOG_LICENSE
                       "\n",
                       Rtxt("Rattle is a registered trademark of Togaware Pty Ltd."),
                       "\n",
                       Rtxt("Rattle was created and implemented by Graham Williams.")),
                tvsep=FALSE)
}

writeCSV <- function(x, file="", ...)
{
  write.csv(x, file=file, row.names=FALSE, ...)
}

rattleTodo <- function(...) cat("Rattle TODO:", ..., "\n")

#-----------------------------------------------------------------------
# MAINLOOP ITERATION
#
# Attempt to get tooltips working forGNU/Linux by starting up gtkMain
# on the window getting focus, and stopping it when it loses
# focus. Based on idea from Felix Andrews.

gtkmain_handler <- function(widget, event)
{
  # 090525 Can't get this one working yet - to be able to turn
  # tooltips on and off. playwith does it?

  #if (! theWidget("tooltip_menuitem")$getActive())
  #  return(gtkmainquit_handler(widget, event))

  # Switch to GTK event loop while the window is in focus (for tooltips)

  if (! crv$.gtkMain)
  {
    crv$.gtkMain <- TRUE
    RGtk2::gtkMain()
  }
  return(FALSE)
}

gtkmainquit_handler <- function(widget, event)
{
  if (crv$.gtkMain)
  {
    crv$.gtkMain <- FALSE
    RGtk2::gtkMainQuit()
  }
  return(FALSE)
}

#-----------------------------------------------------------------------
# RESET RATTLE

resetRattle <- function(new.dataset=TRUE)
{
  # Cleanup various bits of Rattle, as when a new dataset is loaded or
  # a project is loaded. Might also be useful for the New button. If
  # new.dataset is FALSE then just reset various textviews and default
  # options.

  if (new.dataset) setMainTitle()

  if (new.dataset)
  {
    # Initialise CRS

    crs$dataset  <- NULL
    crs$dataname <- NULL
    # crs$dwd      <- NULL
    crs$mtime    <- NULL
    crs$input    <- NULL
    crs$target   <- NULL
    crs$weights  <- NULL
    crs$risk     <- NULL
    crs$ident    <- NULL
    crs$ignore   <- NULL
    crs$nontargets <- NULL # 080426 Started but not yet implemented.
    crs$sample   <- NULL
    crs$sample.on <- TRUE
    crs$sample.seed <- NULL
    crs$tain <- NULL # 100110 For now use crs$sample for the sample until migrate rstat
    crs$validate <- NULL
    crs$test <- NULL
    crs$testset  <- NULL
    crs$testname <- NULL
    crs$transforms <- NULL
    crs$projname <- NULL # 101115
    crs$filename <- NULL # 101115
  }

  # Clear out all current models.

  crs$kmeans   <- NULL
  crs$kmeans.seed <- NULL
  crs$clara    <- NULL
  crs$pam      <- NULL
  crs$hclust   <- NULL
  crs$biclust  <- NULL
  crs$apriori  <- NULL
  crs$page     <- ""
  crs$smodel   <- NULL
  crs$glm      <- NULL
  crs$rpart    <- NULL
  crs$ada      <- NULL
  crs$rf       <- NULL
  crs$svm      <- NULL
  crs$ksvm     <- NULL
  crs$nnet     <- NULL
  crs$survival <- NULL
  crs$perf     <- NULL
  crs$eval     <- NULL

  # Clear all now outdated text views

  resetTextviews()

  # Set all sub tabs back to the default tab page and reflect this in
  # the appropriate radio button.

  # TODO 080423 Change name to RESCALE
  crv$TRANSFORM$setCurrentPage(crv$TRANSFORM.NORMALISE.TAB)
  theWidget("normalise_radiobutton")$setActive(TRUE)
  theWidget("normalise_recenter_radiobutton")$setActive(TRUE)
  theWidget("impute_zero_radiobutton")$setActive(TRUE)
  theWidget("impute_constant_entry")$setText("")
  theWidget("remap_quantiles_radiobutton")$setActive(TRUE)
  theWidget("delete_ignored_radiobutton")$setActive(TRUE)

  crv$EXPLORE$setCurrentPage(crv$EXPLORE.SUMMARY.TAB)
  theWidget("summary_radiobutton")$setActive(TRUE)

  crv$CLUSTER$setCurrentPage(crv$CLUSTER.KMEANS.TAB)
  theWidget("kmeans_radiobutton")$setActive(TRUE)

  crv$MODEL$setCurrentPage(crv$MODEL.RPART.TAB)
  theWidget("rpart_radiobutton")$setActive(TRUE)
  #theWidget("all_models_radiobutton")$setActive(TRUE)

  crv$EVALUATE$setCurrentPage(crv$EVALUATE.CONFUSION.TAB)
  theWidget("evaluate_confusion_radiobutton")$setActive(TRUE)
  theWidget("score_class_radiobutton")$setActive(TRUE)
  theWidget("score_class_radiobutton")$setLabel(Rtxt("Class"))
  theWidget("score_probability_radiobutton")$setLabel(Rtxt("Probability"))

  # Reset the DATA tab. But we don't want to do this because
  # resetRattle is called on loading a database table, and this ends
  # up clearing all the widgets!

  if (new.dataset)
  {
    theWidget("sample_count_spinbutton")$setValue(0)
    theWidget("data_sample_checkbutton")$setActive(FALSE)
    theWidget("data_target_auto_radiobutton")$setActive(TRUE)
    theWidget("data_sample_entry")$setText("70/15/15")
  }

  # 080520 Don't turn these off - it makes sense to allow the user to
  # set these options even before the dataset is loaded.

  # theWidget("target_type_radiobutton")$setSensitive(FALSE)
  # theWidget("data_target_categoric_radiobutton")$setSensitive(FALSE)
  # theWidget("data_target_numeric_radiobutton")$setSensitive(FALSE)

##   theWidget("odbc_dsn_entry")$setText("")
##   theWidget("odbc_combobox")$setActive(-1)
##   theWidget("odbc_limit_spinbutton")$setValue(0)
##   theWidget("odbc_believeNRows_checkbutton")$setActive(FALSE)

  if (packageIsAvailable("SnowballC"))
    theWidget("data_corpus_stem_checkbutton")$setActive(TRUE)

  if (new.dataset)
  {
    # Clear the treeviews.

    theWidget("select_treeview")$getModel()$clear()
    theWidget("impute_treeview")$getModel()$clear()
    theWidget("categorical_treeview")$getModel()$clear()
    theWidget("continuous_treeview")$getModel()$clear()

    theWidget("weight_entry")$setText("")
    theWidget("model_tree_rpart_weights_label")$
    setText("")

    # Data -> Corpus

    theWidget("data_corpus_location_filechooserbutton")$setCurrentFolder(getwd())
    
    # Reset Test

    theWidget("test_distr_radiobutton")$setActive(TRUE)
    theWidget("test_vars1_combobox")$getModel()$clear()
    theWidget("test_vars2_combobox")$getModel()$clear()
    #theWidget("test_vars1_combobox")$setActive(-1)
    #theWidget("test_vars2_combobox")$setActive(-1)
    theWidget("test_groupby_checkbutton")$setActive(TRUE)
    theWidget("test_groupby_target_label")$setText(Rtxt("No Target"))
    theWidget("test_groupby_checkbutton")$setSensitive(TRUE)
    theWidget("test_groupby_target_label")$setSensitive(TRUE)

    # Reset Describe -> Cluster -> KMeans

    theWidget("kmeans_clusters_spinbutton")$setValue(10)
    theWidget("kmeans_seed_spinbutton")$setValue(crv$seed)
    theWidget("kmeans_runs_spinbutton")$setValue(1)
    theWidget("kmeans_stats_button")$setSensitive(FALSE)
    theWidget("kmeans_data_plot_button")$setSensitive(FALSE)
    theWidget("kmeans_discriminant_plot_button")$setSensitive(FALSE)

    # Reset Describe -> Cluster -> Clara

    # Reset Describe -> Cluster -> PAM

    # Reset Describe -> Cluster -> HClust

    theWidget("hclust_clusters_spinbutton")$setValue(10)
    theWidget("hclust_nbproc_spinbutton")$setValue(1)
    theWidget("hclust_dendrogram_button")$setSensitive(FALSE)
    theWidget("hclust_stats_button")$setSensitive(FALSE)
    theWidget("hclust_data_plot_button")$setSensitive(FALSE)
    theWidget("hclust_discriminant_plot_button")$setSensitive(FALSE)

    # Reset Describe -> Cluster -> Biclust

    # Reset Predict -> Tree -> RPart

    theWidget("model_tree_priors_entry")$setText("")
    theWidget("model_tree_loss_entry")$setText("")
    theWidget("rpart_minsplit_spinbutton")$setValue(crv$rpart.minsplit.default)
    theWidget("rpart_maxdepth_spinbutton")$setValue(crv$rpart.maxdepth.default)
    theWidget("model_tree_cp_spinbutton")$setValue(crv$rpart.cp.default)
    theWidget("rpart_minbucket_spinbutton")$setValue(crv$rpart.minbucket.default)
    theWidget("model_tree_include_missing_checkbutton")$setActive(FALSE)
    theWidget("model_tree_rpart_radiobutton")$setActive(TRUE)

    # Reset Predict -> ADA

    showModelAdaExists()
    setGuiDefaultsAda()

    # Reset Predict -> RF

    showModelRFExists()

    # Reset Predict -> SVM

    setGuiDefaultsSVM()

    # Reset Predict -> Survival

    setGuiDefaultsSurvival()

    # Update EXPLORE, MODEL and EVALUATE targets

    theWidget("explot_target_label")$setText(Rtxt("No Target"))
    theWidget("explot_annotate_checkbutton")$setActive(FALSE)
    theWidget("summary_find_entry")$setText("")
    theWidget("benford_bars_checkbutton")$setActive(FALSE)
    theWidget("benford_abs_radiobutton")$setActive(TRUE)
    theWidget("benford_digits_spinbutton")$setValue(1)
    theWidget("explore_correlation_method_combobox")$setActive(0)
    theWidget("pairs_color_combobox")$getModel()$clear()

    theWidget("glm_target_label")$setText(Rtxt("No Target"))
    theWidget("rpart_target_label")$setText(Rtxt("No Target"))
    ##theWidget("gbm_target_label")$setText("No Target")
    theWidget("rf_target_label")$setText(Rtxt("No Target"))
    theWidget("svm_target_label")$setText(Rtxt("No Target"))
    theWidget("nnet_target_label")$setText(Rtxt("No Target"))

    theWidget("evaluate_risk_label")$setText(Rtxt("No risk variable selected"))

    theWidget("evaluate_training_radiobutton")$setActive(TRUE)
    theWidget("evaluate_filechooserbutton")$setFilename("")
    theWidget("evaluate_rdataset_combobox")$setActive(-1)

    # If there is a RATTLE.SCORE.IN defined, as might be from a .Rattle
    # file, then use that for the filename of the CSV evaluate option.

    if (exists("RATTLE.SCORE.IN") && ! is.null(RATTLE.SCORE.IN))
    {
      scorename <- RATTLE.SCORE.IN
      if (not.null(scorename))
      {
        scorename <- path.expand(scorename)

        # If it does not look like an absolute path then add in the
        # current location to make it absolute.

        if (substr(scorename, 1, 1) %notin% c("\\", "/")
            && substr(scorename, 2, 2) != ":")
          scorename <- file.path(getwd(), scorename)
        if (! file.exists(scorename))
        {
          errorDialog(sprintf(Rtxt("The specified SCORE file '%s'",
                                   "(sourced from the .Rattle file through the",
                                   "RATTLE.SCORE.IN variable)",
                                   "does not exist. We will continue",
                                   "as if it had not been speficied."),
                              scorename))

          # Remove the variable (from the global environment where the
          # source command will have plade the bindings) so the rest of
          # the code continues to work on the assumption that it has not
          # been supplied.

          RATTLE.SCORE.IN <- NULL
        }
        else
        {
          theWidget("evaluate_filechooserbutton")$setFilename(scorename)
          theWidget("evaluate_csv_radiobutton")$setActive(TRUE)
        }
      }
    }
  }

  # 100224 Things to do irrespective of whether it is a new dataset.

  showModelRPartExists()

  #091112 resetEvaluateTab("all_inactive")
  #091112 resetEvaluateTab("all_insensitive")
  resetEvaluateTab()

  #theWidget("rpart_evaluate_checkbutton")$hide()
  #theWidget("rf_evaluate_checkbutton")$hide()
  #theWidget("ksvm_evaluate_checkbutton")$hide()
  #theWidget("glm_evaluate_checkbutton")$hide()
  #theWidget("ada_evaluate_checkbutton")$hide()

  ## Update CLUSTER tab

  theWidget("kmeans_hclust_centers_checkbutton")$setActive(FALSE)
  theWidget("hclust_distance_combobox")$setActive(FALSE)
  theWidget("hclust_link_combobox")$setActive(1)
  theWidget("hclust_dendrogram_button")$setSensitive(FALSE)
  theWidget("hclust_clusters_label")$setSensitive(FALSE)
  theWidget("hclust_clusters_spinbutton")$setSensitive(FALSE)
  theWidget("hclust_stats_button")$setSensitive(FALSE)
  theWidget("hclust_data_plot_button")$setSensitive(FALSE)
  theWidget("hclust_discriminant_plot_button")$setSensitive(FALSE)
  if (! isMac()) theWidget("associate_sort_comboboxtext")$setActive(0)
  
  setStatusBar(Rtxt("To Begin: Choose the data source,",
                    "specify the details,",
                    "then click the Execute button."))

}

########################################################################
# UTILITIES

"%notin%" <- function(x,y) ! x %in% y

not.null <- function(x) ! is.null(x)

uri2file <- function(u)
{
  sub("^file://", "", u)
}

listVersions <- function(file="", ...)
{
  result <- installed.packages()[,c("Package", "Version")]
  row.names(result) <- NULL
  write.csv(result, file=file, ...)
  invisible(result)
}

########################################################################
## Common Dialogs

debugDialog <- function(...)
{
  dialog <- RGtk2::gtkMessageDialogNew(NULL, "destroy-with-parent", "info", "ok",
                                "Debug Message:", ...)
  RGtk2::connectSignal(dialog, "response", RGtk2::gtkWidgetDestroy)
}

infoDialog <- function(...)
{
  # If the RGtk2 package's functions are not available, then just
  # issue a warning instead of a popup.

  if (exists("gtkMessageDialogNew"))
  {
    dialog <- RGtk2::gtkMessageDialogNew(NULL, "destroy-with-parent", "info", "close",
                                  ...)
    RGtk2::connectSignal(dialog, "response", RGtk2::gtkWidgetDestroy)
  }
  else
    # 080706 This fails the MS/Windows check with "crv" not defined?????
    if (! isWindows()) warning(...)
}

warnDialog <- function(...)
{
  dialog <- RGtk2::gtkMessageDialogNew(NULL, "destroy-with-parent", "warn", "close",
                                ...)
  RGtk2::connectSignal(dialog, "response", RGtk2::gtkWidgetDestroy)
}

errorDialog <- function(...)
{
  # 110320 Note that this is a non-blocking dialog. Thus it could
  # actually remain active. At times this is useful as the error
  # dialogue contains instructions on "fixing" the error and you can
  # keep the dialogue visible whilst fixing the error.
  
  dialog <- RGtk2::gtkMessageDialogNew(NULL, "destroy-with-parent", "error", "close",
                                ...,
                                sprintf("\n\n%s %s",
                                        crv$appname, crv$version))
  RGtk2::connectSignal(dialog, "response", RGtk2::gtkWidgetDestroy)
  return(FALSE)
}

questionDialog <- function(...)
{
  if (package.installed("RGtk2"))
  {
    dialog <- RGtk2::gtkMessageDialogNew(NULL, "destroy-with-parent", "question",
                                  "yes-no",
                                  ...)
    result <- dialog$run()
    dialog$destroy()
    answer <- result == RGtk2::GtkResponseType["yes"]
  }
  else
  {
    cat(paste(strwrap(...), collapse="\n"))
    answer <- tolower(readline(" (yes/NO) ")) %in% c("yes", "y")
  }
  return(answer)
}

notImplemented <- function(action, window)
{
  ## Popup a little information window for non-implemented functions.

  aname <- action$getName()
  result <- try(atype <- action$typeName(), silent=TRUE)
  if (inherits(result, "try-error")) atype <- NULL

  infoDialog(sprintf(Rtxt("The function you activated (via %s)",
                          "%s is not yet implemented."),
                     aname,
                     ifelse(is.null(atype), "", sprintf("of type %s", atype))))
#  infoDialog(sprintf(paste("The function you activated (via %s)",
#                            "of type %s is not yet implemented."),
#                      action$getName(), action$typeName()))
}

noDatasetLoaded <- function()
{
  # Popup an error dialog if no dataset has been loaded, and return
  # TRUE, otherwise return FALSE.

  if (is.null(crs$dataset))
  {
    errorDialog(Rtxt("No dataset has been loaded at this time.",
                     "\n\nAt a minimum, please load a dataset from the Data tab",
                     "before attempting any other operation.",
                     "\n\nBe sure to Execute the Data tab once the",
                     "data source has been specified."))
    return(TRUE)
  }
  else
    return(FALSE)
}

variablesHaveChanged <- function(action)
{
  # PARAMETERS
  #
  # action: a string that is displayed in the error dialogue.

  if (length(crs$ignore) != length(getSelectedVariables("ignore")) ||
      length(crs$ident) != length(getSelectedVariables("ident")) ||
      length(crs$input) != length(getSelectedVariables("input")))
  {
    errorDialog(sprintf(Rtxt("It appears that there have been changes made",
                             "to the variables in the",
                             "Data tab without the tab being Executed.",
                             "\n\nPlease click Execute on the Data tab before",
                             "%s."),
                        action))
    return(TRUE)
  }
  else
    return(FALSE)
}

# 110703: I used to test if the package name was in the result from
# installed.packages(), but as Brian Ripley points out and from the
# man page for the function, installed.packages() is very slow on
# MS/Windows and on networked file systems as it touches a couple of
# files for each package, and with over a thousand packages installed
# that will be a lot of files. So simply check for the package using
# system.file().

package.installed <- function(package) nchar(system.file(package=package)) > 0
  
packageIsAvailable <- function(pkg, msg=NULL)
{
  appname <- ifelse(exists("crv") && ! is.null(crv$appname), crv$appname, "Rattle")
  localmsg <- sprintf(Rtxt("The package '%s' is required to %s.",
                           "It does not appear to be installed.",
                           "This package (and its dependencies) can be installed",
                           "using the following R command:",
                           "\n\ninstall.packages('%s')",
                           "\n\nThis one-time install will allow access to",
                           "the full functionality of %s.",
                           "\n\nWould you like %s to install the package now?"),
                      pkg, msg, pkg, appname, appname)
  if (! package.installed(pkg))
  {
    if (not.null(msg))
      if (questionDialog(localmsg))
      {
        install.packages(pkg)
        return(TRUE)
      }
    return(FALSE)
  }
  else
    return(TRUE)
}

sampleNeedsExecute <- function()
{
  # Popup an error dialog if sampling needs to be executed and return
  # TRUE.

  # If sampling is active, make sure there is a sample.

  if (theWidget("data_sample_checkbutton")$getActive()
      && is.null(crs$sample))
  {
    errorDialog(Rtxt("Sampling is active but has not been Executed.",
                     "Either ensure you Execute the sampling by clicking",
                     "the Execute button on the Transform tab,",
                     "or else de-activate Sampling on the Data tab."))
    return(TRUE)
  }

  # If sampling is inactive, make sure there is no sample. 080601 Why
  # would I need this test?

###   if (! theWidget("data_sample_checkbutton")$getActive()
###       && not.null(crs$sample))
###   {
###     errorDialog("Sampling is inactive but has not been Executed",
###                  "since being made inactive.",
###                  "Please ensure you Execute the Transform tab",
###                  "after de-activating the Sampling on the Transform tab.")
###         return(TRUE)
###   }

  return(FALSE)
}

errorMessageFun <- function(call, result)
{
  # 100109 Generate a message reporting an error in a function call.

  return(sprintf(Rtxt("An error occured in the call to '%s'.",
                      "The error message was:\n\n%s\n\n%s"),
                 call, result, crv$support.msg))
}

errorMessageCmd <- function(call, result)
{
  # 100109 Generate a message reporting an error in a command line.

  return(sprintf(Rtxt("An error occured in the following command:\n\n%s.",
                      "\n\nThe error message was:\n\n%s\n\n%s"),
                 call, result, crv$support.msg))
}

errorReport <- function(cmd, result)
{
  # A standard command error report that is not being captured by
  # Rattle. Eventually, all of these should be identified by Rattle
  # and a sugggestion given as to how to avoid the error.

  errorDialog(errorMessageCmd(cmd, result))
}

########################################################################
#
# Simplify updates to status bar
#

setMainTitle <- function(title=NULL)
{
  standard <- Rtxt("R Data Miner - [Rattle]")
  if (is.null(title))
    theWidget("rattle_window")$setTitle(standard)
  else
  {
    Encoding(title) <- "UTF-8" # 100415 Just in case? Japanese window title not proper in RStat
    theWidget("rattle_window")$setTitle(sub("]",
                                            sprintf(" (%s)]", title),
                                            standard))
  }
}

setStatusBar <- function(..., sep=" ")
{
  msg <- paste(sep=sep, ...)
  if (length(msg) == 0) msg <-""
  theWidget("statusbar")$push(1, msg)
  while (RGtk2::gtkEventsPending()) RGtk2::gtkMainIterationDo(blocking=FALSE) # Refresh status/windows
  invisible(NULL)
}

reportTimeTaken <- function(tv, time.taken, model, msg)
{
  # 091224 This is called after building a model to report on how long
  # the build took in the text view, to append the time taken to the
  # log for information purposes, and to update the status bar. At
  # least one of and only one of model or msg must be supplied.

  if (missing(model) && missing(msg) || (!missing(model) && !missing(msg)))
    stop(Rtxt("rattle: reportTimeTaken:",
              "one and only one of model/msg must be supplied."))

  time.msg <- sprintf(Rtxt("Time taken: %0.2f %s"),
                      time.taken, Rtxt (attr(time.taken, "units")))

  # Rtxt("secs") Rtxt("mins") So that the above units gets
  # translated. Note also the gap after Rtxt above to avoid it being
  # picked up as a string to be translated.
  
  addTextview(tv, "\n", time.msg, textviewSeparator())
  appendLog(time.msg)

  if (missing(msg))
    msg <- sprintf(Rtxt("The %s model has been built."), model)

  setStatusBar(msg, time.msg)
}


collectOutput <- function(command, use.print=FALSE, use.cat=FALSE,
                          width=getOption("width"), envir=parent.frame())
{
  # TODO Should this use cat or print? Cat translates the \n to a
  # newline and doesn't precede the output by [1].  For pretty output
  # with sprintf() you probably want cat(), but if you have a vector
  # of formatted text and you want to look at it (as data), print()
  # would be better.

  owidth <- getOption("width")
  options(width=width)
  if (use.print)
    command <- paste("print(", command, ")", sep="")
  else if (use.cat)
    command <- paste("cat(", command, ")", sep="")

  # 080829 - Let's try out capture.output as a simpler way of doing
  # sink. Seems to work okay!

  if (FALSE)
  {
    zz <- textConnection("commandsink", "w", TRUE)
    sink(zz)
    result <- try(eval(parse(text=command)))#121212, envir=envir))
    sink()
    close(zz)
  }
  else
  {
    result <- try(commandsink <- capture.output(eval(parse(text=command))))#121212, envir=envir)))
  }

  if (inherits(result, "try-error"))
  {
    if (any(grep("cannot allocate vector", result)) ||
        any(grep("vector size specified is too large", result)))
      errorDialog(Rtxt("E141: The dataset is too large for this operation.",
                       "It is terminating now without any output.",
                       "The R Console may contain further information."))
    else
      errorDialog(sprintf(Rtxt("E142: A command has failed\n\n%s\n\n",
                               "The action you requested has not been completed.",
                               "Refer to the R Console for details."),
                          command))
    commandsink <- Rtxt("No output generated.")
  }
  options(width=owidth)
  return(paste(commandsink, collapse="\n"))
}

########################################################################
##
## Miscellaneous Support
##

theWidget <- function(widget)
{
  #crv$rattleGUI <- Global_.rattleGUI # Global - to avoid a "NOTE" from "R CMD check"

  if (crv$useGtkBuilder)
    return(crv$rattleGUI$getObject(widget))
  else
    return(crv$rattleGUI$getWidget(widget))
}

getNotebookPageLabel <- function(nb, page)
{
  # Given a notebook object and a numeric page (from 0 to npages-1),
  # return the label on the tab for that page.
  
  # 100301 Japanese on MS/Windows returns what might be a Shift-JIS
  # string from nb$getTabLabelText(nb$getNthPage(nb$getCurrentPage()))
  # rather than UTF-8, and so the tab name comparisons fail. For now
  # we assume the tab ordering, and so get the tab page number and
  # then map that to the tab label.

  # 100408 Remove the special code for Japanese - instead, we simply
  # need to ensure the encoding of the string returned from GTK is
  # UTF-8 rather than "unknown". That seems to fix the problem.

  # TODO - Remove the commented code.
  
  ## if (! isJapanese()) # Test this first to avoid too much testing otherwise.
    label <- nb$getTabLabelText(nb$getNthPage(page))
    Encoding(label) <- "UTF-8"
  ## else if (nb == crv$NOTEBOOK)
  ##   label <- switch(page+1,
  ##                   Rtxt ("Data"),
  ##                   Rtxt ("Explore"),
  ##                   Rtxt ("Test"),
  ##                   Rtxt ("Transform"),
  ##                   Rtxt ("Cluster"),
  ##                   Rtxt ("Associate"),
  ##                   Rtxt ("Predictive"),
  ##                   Rtxt ("Evaluate"),
  ##                   Rtxt ("Log"))
  ## else if (nb == crv$EXPLORE)
  ##   label <- switch(page+1,
  ##                   Rtxt ("summary"),
  ##                   Rtxt ("explot"),
  ##                   Rtxt ("correlation"),
  ##                   Rtxt ("prcomp"),
  ##                   Rtxt ("interactive"))
  ## else if (nb == crv$TRANSFORM)
  ##   label <- switch(page+1,
  ##                   Rtxt ("normalise"),
  ##                   Rtxt ("impute"),
  ##                   Rtxt ("remap"),
  ##                   Rtxt ("outliers"),
  ##                   Rtxt ("cleanup"))
  ## else if (nb == crv$CLUSTER)
  ##   label <- switch(page+1,
  ##                   Rtxt ("kmeans"),
  ##                   Rtxt ("clara"),
  ##                   Rtxt ("hclust"),
  ##                   Rtxt ("biclust"))
  ## else if (nb == crv$MODEL)
  ##   label <- switch(page+1,
  ##                   Rtxt ("rpart"),
  ##                   Rtxt ("ada"),
  ##                   Rtxt ("rf"),
  ##                   Rtxt ("svm"),
  ##                   Rtxt ("glm"),
  ##                   Rtxt ("nnet"),
  ##                   Rtxt ("gbm"),
  ##                   Rtxt ("survival"))
  ## else if (nb == crv$EVALUATE)
  ##   label <- switch(page+1,
  ##                   Rtxt ("confusion"),
  ##                   Rtxt ("lift"),
  ##                   Rtxt ("roc"),
  ##                   Rtxt ("precision"),
  ##                   Rtxt ("sensitivity"),
  ##                   Rtxt ("risk"),
  ##                   Rtxt ("pvo"),
  ##                   Rtxt ("score"),
  ##                   Rtxt ("costcurve"))
  ## else
  ##   # Fall through to the default.
  ##   label <- nb$getTabLabelText(nb$getNthPage(page))
  
  return(label)
}

getNotebookPage <- function(notebook, label)
{
  # Obtain the notebook page number given its tab's label's text
  # (already translated using Rtxt when it is passed in.  Return NULL
  # if the label is not found.

  for (i in 0:(notebook$getNPages()-1))
   if (getNotebookPageLabel(notebook, i) == label) return(i)
  return(NULL)
}

getCurrentPageLabel <- function(nb)
{
  return(getNotebookPageLabel(nb, nb$getCurrentPage()))
}

isMac <- function()
{
  # 140307 Added to check for GUI tings not migrated back into the Mac
  # GUI XML.
    
  # Perhaps should use .Platform$OS.type as below for isWindows.
  return(Sys.info()["sysname"] == "Darwin")
}


isWindows <- function()
{
  # The use of .Platform$OS.type is as recommended in the R.version
  # manual page.
  return(.Platform$OS.type == "windows")
}

fixWindowsSlash <- function(s)
{
  if (isWindows()) s <- gsub('\\\\', '/', s)
  return(s)
}

isLinux <- function()
{
  return(.Platform$OS.type == "unix")
}

isJapanese <- function()
{
  # 091222 For plots and pdf export under MS/Windows. Tested by
  # acken_sakakibara@ibi.com

  return(isWindows() && Sys.getlocale("LC_CTYPE") == "Japanese_Japan.932")
}

listBuiltModels <- function(exclude=NULL)
{
  # Build a list of models that have been built.
  models <- c()
  for (m in setdiff(c(crv$PREDICT, crv$DESCRIBE), exclude))
    if (not.null(eval(parse(text=sprintf("crs$%s", m)))))
      models <- c(models, m)
  return(models)
}

########################################################################
## PLOTTING
##
## Callbacks

on_plot_save_button_clicked <- function(action)
{
  # To know which window we are called from we extract the plot
  # number from the window title!!!. This then ensures we save the
  # right device.
  #
  # Also, export to pdf (from Cairo) is not too good it seems. Gets a
  # grey rather than white background. PNG and JPEG look just fine.
  # This is being fixed by Michael Lawrence.

  ttl <- action$getParent()$getParent()$getParent()$getParent()$getTitle()
  savePlotToFileGui(dev.num(ttl))
}

on_plot_copy_button_clicked <- function(action)
{
  ttl <- action$getParent()$getParent()$getParent()$getParent()$getTitle()
  dn <- dev.num(ttl)
  startLog(Rtxt("Copy the plot to the clipboard."))
  appendLog(sprintf(Rtxt("Copy the plot on device %d to the clipboard."), dn),
            sprintf('copyPlotToClipboard(%s)', dn))
  copyPlotToClipboard(dn)
  setStatusBar(sprintf(Rtxt("Plot %d has been copied to the clipboard",
                            "using the PNG format."),
                       dn))
}

on_plot_print_button_clicked <- function(action)
{
  ## To know which window we are called from we extract the plot
  ## number from the window title!!!. This then ensures we save the
  ## right device.

  ttl <- action$getParent()$getParent()$getParent()$getParent()$getTitle()
  dn <- dev.num(ttl)
  startLog(Rtxt("Print the plot."))
  appendLog(sprintf(Rtxt("Send the plot on device %d to the printer."), dn),
            sprintf('printPlot(%s)', dn))
  printPlot(dn)
  setStatusBar(sprintf(Rtxt("Plot %d has been sent to the printer",
                             "using the command: %s."),
                       dn, options("printcmd")))
}

on_plot_close_button_clicked <- function(action)
{
  ttl <- action$getParent()$getParent()$getParent()$getParent()$getTitle()
  dn <- dev.num(ttl)
  dev.off(dn)

  pw <- action$getParentWindow()
  
  # 100830 "destroy" causes R to crash. So try hide - but does that
  # not release the object and hence accumulates memory usage. Does
  # withdraw do any better?
  
  # pw$destroy()
  # pw$hide()
  pw$withdraw()
}

dev.num <- function(title)
{
  # 100408 Return the device number for the device with the given
  # title. This was needed because Japanes on MS/Windows was returning
  # a title in some encoding that was not the original, and sub(".* ",
  # "", ttl) was failing.

  Encoding(title) <- "UTF-8"
  return(as.integer(sub(".* ", "", title)))
}
  

########################################################################

newPlot <- function(pcnt=1)
{
  # Create a new device into which the plot is to go.

  # Trial the use of the Cairo device. This was the only place I
  # needed to change to switch over to the Cairo device. As backup,
  # revert to the x11() or windows() device. TODO Under Windows
  # (R2.13.1/Rattle2.6.9/Gtk2.20.17) the plot in Figure 2.8 of the
  # Rattle book does not show the box plot in the top right plot -
  # only the stars. Seems to be an issue with CairoDevice? For
  # Windows, for now, do not use Cairo by default.

  if (!exists("RStudioGD"))
  {
    if (theWidget("use_cairo_graphics_device")$getActive() &&
        packageIsAvailable("cairoDevice", Rtxt("display plots")))
    {
      if (crv$useGtkBuilder)
      {
        plotGUI <- RGtk2::gtkBuilderNew()
        plotGUI$setTranslationDomain("R-rattle")
      }

      result <- try(etc <- file.path(path.package(package="rattle")[1], "etc"),
                    silent=TRUE)
      if (inherits(result, "try-error"))
        if (crv$useGtkBuilder)
          plotGUI$addFromFile(crv$rattleUI)
        else
          plotGUI <- gladeXMLNew("rattle.glade", root="plot_window", domain="R-rattle")
      else
        if (crv$useGtkBuilder)
          plotGUI$addFromFile(file.path(etc, crv$rattleUI))
        else
          plotGUI <- gladeXMLNew(file.path(etc,"rattle.glade"),
                                 root="plot_window", domain="R-rattle")
      if (crv$useGtkBuilder)
      {
        plotGUI$getObject("plot_window")$show()
        plotGUI$connectSignals()
        da <- plotGUI$getObject("drawingarea")
      }
      else
      {
        gladeXMLSignalAutoconnect(plotGUI)
        da <- plotGUI$getWidget("drawingarea")
      }
      
      cairoDevice::asCairoDevice(da)
      if (isJapanese())
      {
        # 091222 Use a font that MS/Windows can display Japanese
        # characters. Would like to use opar to record old value, but
        # not easy to know where the end of this scope is.
        
        fnt.cmd <- 'par(family=windowsFont("MS Gothic"))'
        appendLog(Rtxt("Use a Japanese font for the plots."), fnt.cmd)
        eval(parse(text=fnt.cmd))
      }
      
      if (crv$useGtkBuilder)
        plotGUI$getObject("plot_window")$setTitle(paste(crv$appname, ": ",
                                                        Rtxt("Plot"), " ",
                                                        dev.cur(), sep=""))
      else
        plotGUI$getWidget("plot_window")$setTitle(paste(crv$appname, ": ",
                                                        Rtxt("Plot"), " ",
                                                        dev.cur(), sep=""))
    }
  }

  if (pcnt==1)
    layout(matrix(c(1), 1, 1, byrow=TRUE))
  else if (pcnt==2)
    layout(matrix(c(1,2), 2, 1, byrow=TRUE))
  else if (pcnt==3)
    layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE))
  else if (pcnt==4)
    layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
  else if (pcnt==5)
    layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow=TRUE))
  else if (pcnt==6)
    layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow=TRUE))
  else if (pcnt==7)
    layout(matrix(c(1,1,2,3,3,4,5,6,7), 3, 3, byrow=TRUE))
  else if (pcnt==8)
    layout(matrix(c(1,1,2,3,4,5,6,7,8), 3, 3, byrow=TRUE))
  else if (pcnt==9)
    layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow=TRUE))
}

########################################################################

copyPlotToClipboard <- function(dev.num=dev.cur())
{
  # This is designed to be called from the Gtk window that displays
  # the Cairo device, to copy the plot displayed there into the
  # Clipboard. It has not been tested on non-Cairo devices.
  #
  # We can place a GdkPixbuf image into the CLIPBOARD using
  # GtkClipboardSetImage. I've not figured out yet how to get the
  # image directly from the Cairo device as a GdkPixbuf. So instead I
  # save to PNG file then load that file as a GdkPixmap then copy that
  # to the clipboard.
  #
  # This works for GNU/Linux and more recent MS/Windows (e.g., on my
  # recent Dell laptop but not on ATOnet computers). It has not been
  # tested on Mac/OSX. Perhaps it is a bug and needs to be reported to
  # Michael Lawrence. Michael has also mentioned a new version of
  # cairoDevice supporting cairo backends for PDF, PS, and PNG to
  # output in those formats directly (070406).
  #
  # Note that in oodraw, for example, you can select an object, then
  # grab the selection and have it available in R:
  #
  # im <- gtkClipboardGet("CLIPBOARD")$waitForImage()
  #
  # Of course you can also load the image from file:
  #
  # im <- gdkPixbufNewFromFile("audit_auto_plot3.png")$retval
  #
  # Once we have the image:
  #
  # gtkClipboardGet("CLIPBOARD")$setImage(im)
  #
  # Which can then be pasted into oowriter, for example.

  temp.name <- paste(tempfile(), ".png", sep="")
  savePlotToFile(temp.name, dev.num)
  im <- RGtk2::gdkPixbufNewFromFile(temp.name)$retval
  RGtk2::gtkClipboardGet("CLIPBOARD")$setImage(im)
  file.remove(temp.name)
}

savePlotToFileGui <- function(dev.num=dev.cur(), name="plot")
{

  if (is.null(dev.list()))
  {
    warnDialog(Rtxt("There are currently no active graphics devices.",
                    "So there is nothing to export!",
                    "Please click the Execute button (or press F2)",
                    "to obtain a plot to export."))
    return()
  }

  # Obtain a filename to save to. Ideally, this would also prompt for
  # the device to export, and the fontsize, etc.

  dialog <- RGtk2::gtkFileChooserDialog(paste(Rtxt("Export Graphics"),
                                       " (.pdf, .png, .jpg",
                                       ifelse(isWindows(), ", wmf", ""),
                                       ")", sep=""),
                                 NULL, "save",
                                 "gtk-cancel", RGtk2::GtkResponseType["cancel"],
                                 "gtk-save", RGtk2::GtkResponseType["accept"])
  dialog$setDoOverwriteConfirmation(TRUE)

  if(not.null(crs$dataname))
    dialog$setCurrentName(paste(get.stem(crs$dataname),
                                "_", name, ".pdf", sep=""))

  ff <- RGtk2::gtkFileFilterNew()
  if (isWindows())
    ff$setName(paste(Rtxt("Graphics Files"), "(pdf png jpg wmf)"))
  else
    ff$setName(paste(Rtxt("Graphics Files"), "(pdf png jpg)"))
  ff$addPattern("*.pdf")
  ff$addPattern("*.png")
  ff$addPattern("*.jpg")
  if (isWindows()) ff$addPattern("*.wmf")
  dialog$addFilter(ff)

  ff <- RGtk2::gtkFileFilterNew()
  ff$setName(Rtxt("All Files"))
  ff$addPattern("*")
  dialog$addFilter(ff)

  if (dialog$run() == RGtk2::GtkResponseType["accept"])
  {
    save.name <- dialog$getFilename()
    dialog$destroy()
  }
  else
  {
    dialog$destroy()
    return()
  }

#  if (get.extension(save.name) == "")
#    save.name <- sprintf("%s.pdf", save.name)

  startLog(Rtxt("Save the plot to a file."))
  appendLog(sprintf(Rtxt("Save the plot on device %d to a file."), dev.num),
            ifelse(packageIsAvailable("cairoDevice") &&
                   theWidget("use_cairo_graphics_device")$getActive(),
                   sprintf('library(cairoDevice)\n'), ''),
            sprintf('savePlotToFile("%s", %s)',
                    fixWindowsSlash(save.name), dev.num))

  if (savePlotToFile(save.name, dev.num))
    setStatusBar(sprintf(Rtxt("Plot %d has been exported to the file %s."),
                         dev.num, save.name))
}

savePlotToFile <- function(file.name, dev.num=dev.cur())
{
  cur <- dev.cur()
  dev.set(dev.num)
  ext <- get.extension(file.name)
  if (ext == "pdf")
    # Set version to 1.4 since dev.copy from a Cairo device needs
    # this.  It is done automatically with a warning anyhow, but might
    # as well avoid the warning so as not to worry anyone.  091222 Add
    # the test for Japanese to add the family option so we get
    # Japanese fonts. This also kind of works on GNU/Linux but the
    # viewer compains about missing fonts. Cairo_pdf works just fine
    # on GNU/Linux, and if it works also on MS/Windows in Japanese the
    # we will go with that.
    #if (isJapanese())
    #  dev.copy(pdf, file=file.name, width=10, height=10, version="1.4", family="Japan1")
    #else
    dev.copy(cairoDevice::Cairo_pdf, file=file.name, width=10, height=10)
  else if (ext == "png")
    dev.copy(png, file=file.name, width=1000, height=1000)
  else if (ext == "jpg")
    dev.copy(jpeg, file=file.name, width=1000, height=1000)
  else if (ext == "wmf")
    eval(parse(text=sprintf("dev.copy(win.metafile, file='%s', width=10, height=10)", file.name)))
  else
  {
    infoDialog(sprintf(Rtxt("The supplied filename extension '%s'",
                            "is not supported."), ext))
    return(FALSE)
  }
  dev.off()
  dev.set(cur)
  return(TRUE)
}

printPlot <- function(dev.num=dev.cur())
{
  cur <- dev.cur()
  dev.set(dev.num)
  if (isWindows())
    eval(parse(text="dev.print(win.print)"))
  else
    dev.print()
  dev.set(cur)
}

########################################################################

genPlotTitleCmd <- function(..., vector=FALSE)
{
  # 080817 Use month name rather than number - less ambiguous.

  if (! exists("crv"))
  {
    crv <- list()
    crv$appname <- "Rattle"
    crv$verbose <- TRUE
    crv$show.timestamp <- TRUE
  }

  main = paste(...)
  if(vector)
  {
    if (! crv$verbose)
      sub <- ""
    else if (crv$show.timestamp)
      sub <- sprintf("%s %s %s", crv$appname,
                     format(Sys.time(), "%Y-%b-%d %H:%M:%S"), Sys.info()["user"])
    else
      sub <- sprintf(Rtxt("Generated by %s"), crv$appname)
    return(c(main, sub))
  }
  else
  {
    if (! crv$verbose)
      sub <- ""
    else if (crv$show.timestamp)
      sub <- sprintf(paste('paste("%s", format(Sys.time(),',
                           '"%%Y-%%b-%%d %%H:%%M:%%S"), Sys.info()["user"])'),
                     crv$appname)
    else
      sub <- sprintf('paste("%s")', sprintf(Rtxt("Generated by %s"), crv$appname))
    
    return(sprintf('title(main="%s",\n    sub=%s)', main, sub))
  }
}

set.cursor <- function(cursor="left-ptr", message=NULL)
{
  if (! is.null(message)) setStatusBar(message)
  theWidget("rattle_window")$getWindow()$
  setCursor(RGtk2::gdkCursorNew(cursor))

  # 091106 For now, set cursor specifically on the textview
  # windows. Under Ubuntu it is not needed, but is on Vista. Is this a
  # GTK+ issue? Remove this once MS/Windows no longer has this problem.

  # 091106 The first approach, lapply, did not work! Whlist all the
  # textview widgets do exist, the getWindow returned NULL unless the
  # textview had been visited. So, instead, loop through the
  # textviews.

  # lapply(allTextviews(), function(x) theWidget(x)$
  #            getWindow("GTK_TEXT_WINDOW_TEXT")$
  #            setCursor(gdkCursorNew(cursor)))

  # 111203 On Mac this started causing attmpt to apply non-funciton
  # errors, since the textviews are not yet defined on starting up
  # Rattle.  This started happening with R 2.14.0 on Mac after failing
  # to properly load rattle.ui.  I could get rid of thes for now and
  # test if it works okay on Linux/Windows/Mac for textviews, but on
  # Mac at least, some textviews were not changing cursor. I should test
  # if theWdiget(tv) is NULL then don't proceed.
  
  for (tv in allTextviews())
  {
    win <- theWidget(tv)$getWindow("GTK_TEXT_WINDOW_TEXT")
    if (! is.null(win)) win$setCursor(RGtk2::gdkCursorNew(cursor))
  }
}

simplifyNumberList <- function(nums)
{
  ## Convert 3 4 6 7 8 9 10 12 14 16 17 18 19 21 to
  ## "3:4, 6:10, 12, 14, 16:19, 21"

  if (length(nums) == 1)
    return(sprintf("%s", nums))
  else if (is.null(nums) || length(nums) == 0)
    return(NULL)

  result <- ""
  start <- nums[1]
  len <- 1

  for (i in 2:length(nums))
  {
    if (nums[i] != start + len)
    {
      if (len == 1)
        result <- sprintf("%s, %d", result, start)
      else
        result <- sprintf("%s, %d:%d", result, start, nums[i-1])
      start <- nums[i]
      len <- 1
    }
    else
      len <- len + 1
  }

  if (len == 1)
    result <- sprintf("%s, %d", result, start)
  else
    result <- sprintf("%s, %d:%d", result, start, nums[i])

  result <- sub('c\\(, ', 'c(', sprintf("c(%s)", result))
  return(result)
}

get.extension <- function(path)
{
  ## Extract and return the extension part of a filename

  parts <- strsplit(path, "\\.")[[1]]
  if (length(parts) > 1)
    last <- parts[length(parts)]
  else
    last <- ""
  last
}

get.stem <- function(path)
{
  # Given a filename PATH extract the basename, and from this, the
  # name without an extension.  090718 If the PATH supplied is a
  # string with no extension than just return the PATH.

  parts <- strsplit(basename(path), "\\.")[[1]]
  if (length(parts) > 1)
    last <- paste(parts[seq_len(length(parts)-1)], collapse=".")
  else
    last <- parts
  last
}

########################################################################
#
# Shared callbacks
#

on_rattle_window_delete_event <- function(action, window)
{
  if (crv$close %in% c("quit", "ask"))
  {
    msg <-sprintf(Rtxt("Do you want to terminate %s?"), crv$appname)
    if (!questionDialog(msg))
      return(TRUE)
    else
      if (crv$close == "quit")
        quit(save="no")
      else
        return(FALSE)
  }
}

close_rattle <- function(action, window)
{
  # 090401 This callback seems to be called after the window is
  # destroyed!!!  So the question serves no purpose... Not clear how
  # to fix that.

  closeRattle()
}

quit_rattle <- function(action, window)
{
  # 080815 This function used to return NULL or "yes" and I always
  # tested whether it's results was NULL. But why not return a
  # logical? Start doing that now, by returning TRUE instead of "yes",
  # and look to return FALSE instead of NULL on a negative response to
  # the question.

  closeRattle(TRUE)
}

closeRattle <- function(ask=FALSE)
{
  if (ask || crv$close %in% c("quit", "ask"))
  {
    msg <- sprintf(Rtxt("Do you want to terminate %s?"), crv$appname)
    if (!questionDialog(msg)) return(FALSE)
  }

  # Don't remove the graphics for now. In moving to the Cairo device,
  # this blanks the device, but does not destroy the containing
  # window. I wonder if there is some way to get a list of the plot
  # windows, and destroy each one?

  # graphics.off() # for (i in dev.list()) dev.off(i)

  # 080523 When this is called as a callback from the destroy signal
  # of the GtkObject, the window has already been destroyed, so no
  # need to try again.

  rw <- theWidget("rattle_window")
  if (not.null(rw)) rw$destroy()

  # Communicate to R that Rattle has finished. This is used by the
  # rattle script on GNU/Linux using the littler package which allows
  # one to use R as a scripting language. But rattle dispatches
  # itself from R, and so normally the script immediately
  # terminates. Instead we can have a loop that checks if crv$rattleGUI
  # is NULL, and when it is we finish! Seems to work.

  crv$rattleGUI <- NULL

  # 080511 Restore options to how they were when Rattle was started.

  options(crv$options)

  # if (crv$tooltiphack) gtkMainQuit() # Only needed if gtkMain is run.

  # 080906 Deal with R not finishing up when rattle is called from
  # littler or R CMD BATCH and we close rather than quit.

  if (crv$close == "quit") quit(save="no")

}

interrupt_rattle <- function(action, window)
{
  # The multicore or fork packages may provide some hope under
  # GNU/Linux, but not MS/Wdinwos. Under MS the Esc seems to send a
  # SIGBREAK to the R process. How to do that?

  infoDialog(Rtxt("This operation is not yet functional."))
}

########################################################################

## General Menu Callbacks

on_rattle_menu_activate <- function(action, window)
{
  browseURL("http://rattle.togaware.com")
}

on_delete_menu_activate <- notImplemented

on_connectr_toolbutton_clicked <- function(action, window)
{
  browseURL("http://connect-r.com/posting.php?mode=post&f=2")
}

## Map the unchanged glade defaults

on_cut1_activate <- notImplemented

on_about_menu_activate <-  function(action, window)
{
  result <- try(etc <- file.path(path.package(package="rattle")[1], "etc"),
                silent=TRUE)
  if (crv$useGtkBuilder)
  {
    about <<- RGtk2::gtkBuilderNew()
    about$setTranslationDomain("R-rattle")
  }
  
  if (inherits(result, "try-error"))
    if (crv$useGtkBuilder)
      about$addFromFile(crv$rattleUI)
    else
    about <- gladeXMLNew("rattle.glade", root="aboutdialog", domain="R-rattle")
  else
    if (crv$useGtkBuilder)
      about$addFromFile(file.path(etc, crv$rattleUI))
    else
      about <- gladeXMLNew(file.path(etc, "rattle.glade"),
                           root="aboutdialog", domain="R-rattle")

  if (crv$useGtkBuilder)
  {
    ab <- about$getObject("aboutdialog")
    ab$show()
  }
  else
    ab <- about$getWidget("aboutdialog")

  ab$setVersion(paste0(crv$version,
#LICENSE
                       ""))

  configureAbout(ab)

  if (crv$useGtkBuilder)
    about$connectSignals()
  else
    gladeXMLSignalAutoconnect(about)
}

configureAbout <- function(ab)
{
  ab["program-name"] <- "Rattle"
  ab$setCopyright(paste(DATE, "\n\n", COPYRIGHT, "\n" ,
                        Rtxt("All rights reserved.")))
  
#XX#  ab$setLicense(paste("This program (Rattle) is copyright software, owned by Togaware Pty Ltd.",
#XX#                      "\n\nThis program is licensed and distributed by Togaware Pty Ltd",
#XX#                      "\nto XLICENSEEX for up to XNUSERSX users until XEXPIREX.",
#XX#                      "\nThe license number is XSNX.",
#XX#                      "\n\nThe program is made available under the terms of the GNU",
#XX#                      "\nGeneral Public License as published by the Free",
#XX#                      "\nSoftware Foundation; either version 2 of the License, or (at your",
#XX#                      "\noption) any later version. See the file gpl-license in the",
#XX#                      "\ndistribution and at http://www.gnu.org/copyleft/gpl.html for details.",
#XX#                      "\n\nThis program is distributed without any warranty; without even the",
#XX#                      "\nimplied warranty of merchantability or fitness for a particular purpose.",
#XX#                      "\nPlease see the GNU General Public License for more details.",
#XX#                      "\n\nBy using this program the licensee acknowledges that they have",
#XX#                      "\nevaluated the program and accept the program as is."))

}


on_paste1_activate <- notImplemented
on_copy1_activate <- notImplemented

on_tooltips_activate <- function(action, window)
{

  ## infoDialog("Currently this functionality is not implemented.",
  ##             "It is awaiting some insight into how to get hold of",
  ##             "the glade GtkTooltips group, which can then be",
  ##             "disabled or enabled as requested.")

  if(action$getActive())
  {
    myWin <- theWidget("rattle_window")
    myWin$addEvents(RGtk2::GdkEventMask["focus-change-mask"])
    RGtk2::gSignalConnect(myWin, "focus-in-event", gtkmain_handler)
    RGtk2::gSignalConnect(myWin, "focus-out-event", gtkmainquit_handler)
    RGtk2::gSignalConnect(myWin, "delete-event", gtkmainquit_handler)
  }
  ## else
  ## {
  ##   infoDialog("Currently the functionality to turn tooltips off",
  ##              "is not implemented.")
  ## }
}

on_verbose_menuitem_toggled <- function(action, window)
{
  crv$verbose <- theWidget("verbose_menuitem")$getActive()
}

##----------------------------------------------------------------------

## Miscellaneous callbacks

on_notebook_switch_page <- function(notebook, window, page)
{
  ## notebook is the GtkNotebook object.
  ## window is ??.
  ## page is the index of the page switched to.

  #ct <- current_(page)

  switchToPage(page)
}

on_tools_data_activate <- function(action, window)
{
  crv$NOTEBOOK$setCurrentPage(getNotebookPage(crv$NOTEBOOK,
                                              crv$NOTEBOOK.DATA.NAME))
  switchToPage(crv$NOTEBOOK.DATA.NAME)
}

on_tools_test_activate <- function(action, window)
{
  crv$NOTEBOOK$setCurrentPage(getNotebookPage(crv$NOTEBOOK,
                                              crv$NOTEBOOK.TEST.NAME))
  switchToPage(crv$NOTEBOOK.TEST.NAME)
}

on_tools_transform_activate <- function(action, window)
{
  crv$NOTEBOOK$setCurrentPage(getNotebookPage(crv$NOTEBOOK,
                                              crv$NOTEBOOK.TRANSFORM.NAME))
  switchToPage(crv$NOTEBOOK.TRANSFORM.NAME)
}

on_tools_explore_activate <- function(action, window)
{
  crv$NOTEBOOK$setCurrentPage(getNotebookPage(crv$NOTEBOOK,
                                              crv$NOTEBOOK.EXPLORE.NAME))
  switchToPage(crv$NOTEBOOK.EXPLORE.NAME)
}

on_tools_cluster_activate <- function(action, window)
{
  crv$NOTEBOOK$setCurrentPage(getNotebookPage(crv$NOTEBOOK,
                                              crv$NOTEBOOK.CLUSTER.NAME))
  switchToPage(crv$NOTEBOOK.CLUSTER.NAME)
}

on_tools_model_activate <- function(action, window)
{
  crv$NOTEBOOK$setCurrentPage(getNotebookPage(crv$NOTEBOOK,
                                              crv$NOTEBOOK.MODEL.NAME))
  switchToPage(crv$NOTEBOOK.MODEL.NAME)
}

on_tools_evaluate_activate <- function(action, window)
{
  crv$NOTEBOOK$setCurrentPage(getNotebookPage(crv$NOTEBOOK,
                                              crv$NOTEBOOK.EVALUATE.NAME))
  switchToPage(crv$NOTEBOOK.EVALUATE.NAME)
}

on_tools_log_activate <- function(action, window)
{
  crv$NOTEBOOK$setCurrentPage(getNotebookPage(crv$NOTEBOOK,
                                              crv$NOTEBOOK.LOG.NAME))
  switchToPage(crv$NOTEBOOK.LOG.NAME)
}

switchToPage <- function(page)
{

  # Blank the status bar whenever we change pages

  setStatusBar()

  # This function used to accept numeric pages, so check for that and
  # convert to the page name rather than the now changing page number
  # (page numbers used to be fixed).

  if (is.numeric(page))
    page <- getNotebookPageLabel(crv$NOTEBOOK, page)

  # Record the current page so when we change we know which was last.

  crs$page <- page

}
