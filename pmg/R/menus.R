### Define the main menu here



pmg.menu = list()
pmg.menu$File$"Source file.."$handler =
  function(h,...) gfile("Source file",type="open", action="source")
pmg.menu$File$"Source file.."$icon="file"
pmg.menu$File$"Save Workspace..."$handler =
  function(h,...) gfile("Save workspace",type="save", action="save.image")
pmg.menu$File$"Save Workspace..."$icon = "save"
pmg.menu$File$"Restore Workspace"$handler =
  function(h, ...) gfile("Restore workspace",type="open", action="load")
pmg.menu$File$"Restore Workspace"$icon = "revert-to-saved"
pmg.menu$File$"Load package..."$handler =
  function(h,...) pmg.loadPackages()
pmg.menu$File$"Install CRAN package..."$handler =
  function(h,...) pmg.installCRANPackage()
pmg.menu$File$"Install CRAN package..."$icon = "network"
pmg.menu$File$"Install local package..."$handler =
  function(h,...) {
    old = options("repos")$repos; options("repos"=NULL);
    gfile("Select a package file...","open",action="install.packages",
          filter=   list(
            "tar.gz files"=list(
              patterns=c("*.tgz","*.tar.gz")
              ),
            "zip files"=list(
              patterns=c("*.zip")
              ),
            "All files"=list(
              patterns=c("*")
              )
            )
          )
    options("repos"=old)
  }
pmg.menu$File$"Install local package..."$icon = "file"
pmg.menu$File$"Set working directory..."$handler =
  function(h,...) gfile("Select a directory","selectdir",action="setwd")
pmg.menu$File$"Set working directory..."$icon = "directory"

pmg.menu$File$"pmg options..."$handler =
  function(h,...) pmg.options()
pmg.menu$File$"pmg options..."$icon = "preferences"
pmg.menu$File$"View window list"$handler = function(h,...) pmgWC$show()
pmg.menu$File$"Exit pmg"$handler =
  function(h,...)  {
    dispose(pmg.dialogs.window)
    assignInNamespace("pmg.dialogs.window", NULL,"pmg")
    pmg.closeAll()
  }
pmg.menu$File$"Exit pmg"$icon ="quit"
##
##
pmg.menu$Data$browseEnv$handler =
  function(h,...) browseEnv()
pmg.menu$Data$"Load data set..."$handler =
  function(h,...) pmg.viewDataSets()
pmg.menu$Data$"Import data set..."$handler =
  function(h,...) pmg.specifyFileForImport();
pmg.menu$Data$"Write data as CSV file..."$handler =
  function(h,...) pmg.gw(write.csv.list)
#pmg.menu$Data$"Import data set..."$"import table..."$handler =
#  function(h,...) pmg.gw(read.table.list)
#pmg.menu$Data$"Import data set..."$"import csv file..."$handler =
#  function(h,...) pmg.gw(read.csv.list)
#pmg.menu$Data$"Import data set..."$"import fwf file..."$handler =
#  function(h,...) pmg.gw(read.fwf.list)
## dynamic
pmg.menu$Data$"Dynamic summaries"$handler =
  function(h,...) dSummaryDialog()
pmg.menu$Data$"Dynamic summaries"$icon = "execute"
pmg.menu$Data$"Univariate summaries"$table$handler =
  function(h,...) pmg.gw(table.list)
pmg.menu$Data$"Univariate summaries"$"stem and leaf"$handler =
  function(h,...) pmg.gw(stem.list)
pmg.menu$Data$"Univariate summaries"$"summary"$handler =
  function(h,...) pmg.gw(summary.list)
pmg.menu$Data$"Univariate summaries"$"summary"$icon = "info"
pmg.menu$Data$"Univariate summaries"$mean$handler =
  function(h,...) pmg.gw(mean.list)
pmg.menu$Data$"Univariate summaries"$median$handler =
  function(h,...) pmg.gw(median.list)
pmg.menu$Data$"Univariate summaries"$"standard deviation"$handler =
  function(h,...) pmg.gw(sd.list)
pmg.menu$Data$"Univariate summaries"$IQR$handler =
  function(h,...) pmg.gw(IQR.list)
pmg.menu$Data$"Univariate summaries"$mad$handler =
  function(h,...) pmg.gw(mad.list)
pmg.menu$Data$"Univariate summaries"$quantiles$handler =
  function(h,...) pmg.add(quantileWidget(),label="quantile()") # add here
pmg.menu$Data$"Univariate summaries"$skewness$handler =
  function(h,...) pmg.gw(skewnessList)
pmg.menu$Data$"Univariate summaries"$kurtosis$handler =
  function(h,...) pmg.gw(kurtosisList)
##
pmg.menu$Data$"Bivariate summaries"$correlation$handler =
  function(h,...) pmg.gw(cor.list)
pmg.menu$Data$"Bivariate summaries"$"Cross tabulation"$handler =
  function(h,...) pmg.gw(xtabs.list)
##
pmg.menu$Data$"Random data"$"Cumulative Probabilities"$handler = function(h,...) 
  add(pmg.dialog.notebook,dpqrfuncs(type="p"),label="p funcs")
pmg.menu$Data$"Random data"$"Probabilities"$handler = function(h,...) 
  add(pmg.dialog.notebook,dpqrfuncs(type="d"),label="d funcs")
pmg.menu$Data$"Random data"$Quantiles$handler = function(h,...) 
  add(pmg.dialog.notebook,dpqrfuncs(type="q"),label="q funcs")
pmg.menu$Data$"Random data"$"Random samples"$handler = function(h,...) 
  add(pmg.dialog.notebook,dpqrfuncs(type="r"),label="r funcs")


pmg.menu$Data$"Random data"$"Sample"$handler =
  function(h,...) pmg.gw(sample.list)
##

##########
## Simulation. What else?
pmg.menu$Data$Simulation$"Repeat trials"$handler = function(h,...) {
  add(pmg.dialog.notebook, repeatTrialsGUI(), label = "Repeat trials")
}




## Manipulate
if("reshape" %in% .packages(TRUE)) {
  pmg.menu$Data$"Manipulate"$reshape$melt$handler = function(h,...) pmg.meltGUI()
  pmg.menu$Data$"Manipulate"$reshape$cast$handler = function(h,...) pmg.castGUI()
}
#pmg.menu$Data$"Manipulate"$"subset"$handler =
#  function(h,...) pmg.gw(subset.list)
pmg.menu$Data$"Manipulate"$"subset"$handler =
  function(h,...) add(pmg.dialog.notebook,pmg.subset.dialog(),label="subset()")
pmg.menu$Data$"Manipulate"$"subset"$icon = "subset"
pmg.menu$Data$"Manipulate"$"subset"$handler =
  function(h,...) add(pmg.dialog.notebook,pmg.subset.dialog(),label="subset()")
pmg.menu$Data$"Manipulate"$"stack"$handler =
  function(h,...) pmg.gw(stack.list)
pmg.menu$Data$"Manipulate"$"unstack"$handler =
  function(h,...) pmg.gw(unstack.list)
pmg.menu$Data$"Manipulate"$"Edit data frame properties"$handler =
  function(h,...) add(pmg.dialog.notebook,pmg.edit.dataframe.properties.dialog(),label="edit properties")
pmg.menu$Data$"Manipulate"$"Edit data frame properties"$icon = "properties"
##
pmg.menu$Data$"Coerce"$"as.numeric"$handler =
  function(h,...) pmg.gw(as.numeric.list)
pmg.menu$Data$"Coerce"$"as.character"$handler =
  function(h,...) pmg.gw(as.character.list)
pmg.menu$Data$"Coerce"$"as.data.frame"$handler =
  function(h,...) pmg.gw(as.data.frame.list)
pmg.menu$Data$"Coerce"$"as.matrix"$handler =
  function(h,...) pmg.gw(as.matrix.list)
pmg.menu$Data$"Coerce"$"matrix"$handler =
  function(h,...) pmg.gw(matrix.list)
pmg.menu$Data$"Coerce"$"groupedData"$handler =
  function(h,...) pmg.gw(groupedData.list)
pmg.menu$Data$"Coerce"$"factor"$handler =
  function(h,...) pmg.gw(factor.list)
##
## Plots
## Dynamic widget
pmg.menu$Plots$"Lattice explorer"$handler = function(h,...) {
  dLatticeExplorer(container=pmgWC$new("Lattice explorer", v=T))
}
pmg.menu$Plots$"Lattice explorer"$icon = "execute"
###
pmg.menu$Plots$"Set plot parameters"$Setup$handler =
  function(h,...) pmg.gw(par.setup.list)
pmg.menu$Plots$"Set plot parameters"$Setup$icon = "preferences"
pmg.menu$Plots$"Set plot parameters"$Axes$handler =
  function(h,...) pmg.gw(par.axes.list)
pmg.menu$Plots$"Set plot parameters"$Axes$icon = "preferences"
pmg.menu$Plots$"Set plot parameters"$Colors$handler =
  function(h,...) pmg.gw(par.colors.list)
pmg.menu$Plots$"Set plot parameters"$Colors$icon = "preferences"
pmg.menu$Plots$"Set plot parameters"$Fonts$handler =
  function(h,...) pmg.gw(par.fonts.list)
pmg.menu$Plots$"Set plot parameters"$Fonts$icon = "preferences"
pmg.menu$Plots$"Set plot parameters"$"Number of figures"$handler =
  function(h,...) pmg.gw(par.nofigures.list)
pmg.menu$Plots$"Set plot parameters"$"Number of figures"$icon = "preferences"

##
pmg.menu$Plots$univariate$"barplot"$handler = 
  function(h,...) pmg.gw(barplot.list)
pmg.menu$Plots$univariate$"barplot"$icon="barplot"
pmg.menu$Plots$univariate$"piechart"$handler = 
  function(h,...) pmg.gw(piechart.list)
pmg.menu$Plots$univariate$"boxplot"$handler = 
  function(h,...) pmg.gw(univariate.boxplot.list)
pmg.menu$Plots$univariate$"boxplot"$icon = "boxplot"
pmg.menu$Plots$univariate$"histogram"$handler = 
  function(h,...) pmg.gw(hist.list)
pmg.menu$Plots$univariate$"histogram"$icon ="hist"
pmg.menu$Plots$univariate$"density plot"$handler = 
  function(h,...) pmg.gw(densityplot.list)
pmg.menu$Plots$univariate$"quantile-normal plot"$handler = 
  function(h,...) pmg.gw(qqnorm.list)
pmg.menu$Plots$univariate$"stripchart"$handler = 
  function(h,...) pmg.gw(stripchart.list)
pmg.menu$Plots$univariate$"dotchart"$handler = 
  function(h,...) pmg.gw(dotchart.list)
pmg.menu$Plots$univariate$"ecdf"$handler = 
  function(h,...) pmg.gw(ecdf.list)
##
pmg.menu$Plots$bivariate$"boxplot"$handler = 
  function(h,...) pmg.gw(bivariate.boxplot.list)
pmg.menu$Plots$bivariate$"boxplot"$icon = "boxplot"
pmg.menu$Plots$bivariate$"scatterplot"$handler = 
  function(h,...) pmg.gw(scatterplot.list)
pmg.menu$Plots$bivariate$"scatterplot"$icon = "points"
pmg.menu$Plots$bivariate$"sunflower plot"$handler = 
  function(h,...) pmg.gw(sunflower.list)
pmg.menu$Plots$bivariate$"quantile-quantile plot"$handler = 
  function(h,...) pmg.gw(qqplot.list)
##
pmg.menu$Plots$multivariate$"plot"$handler = 
  function(h,...) pmg.gw(scatterplot.model.list)
pmg.menu$Plots$multivariate$"plot"$icon = "plot"
pmg.menu$Plots$multivariate$"boxplot"$handler = 
  function(h,...) pmg.gw(model.boxplot.list)
pmg.menu$Plots$multivariate$"boxplot"$icon = "boxplot"
pmg.menu$Plots$multivariate$"pairs plot"$handler = 
  function(h,...) pmg.gw(pairs.list)
##
pmg.menu$Plots$"Lattice graphics"$"xyplot"$handler = 
  function(h,...) pmg.gw(xyplot.list)
pmg.menu$Plots$"Lattice graphics"$"dotplot"$handler = 
  function(h,...) pmg.gw(dotplot.list)
pmg.menu$Plots$"Lattice graphics"$"barchart"$handler = 
  function(h,...) pmg.gw(barchart.list)
pmg.menu$Plots$"Lattice graphics"$"stripplot"$handler = 
  function(h,...) pmg.gw(stripplot.list)
pmg.menu$Plots$"Lattice graphics"$"bwplot"$handler = 
  function(h,...) pmg.gw(bwplot.list)
##
pmg.menu$Plots$"Add to graphic"$"points"$handler = 
  function(h,...) pmg.gw(add.points.list)
pmg.menu$Plots$"Add to graphic"$"points"$icon = "points"
pmg.menu$Plots$"Add to graphic"$"lines"$handler = 
  function(h,...) pmg.gw(add.lines.list)
pmg.menu$Plots$"Add to graphic"$"lines"$icon = "lines"
pmg.menu$Plots$"Add to graphic"$"density"$handler = 
  function(h,...) pmg.gw(add.density.list)
pmg.menu$Plots$"Add to graphic"$"curve"$handler = 
  function(h,...) pmg.gw(add.curve.list)
pmg.menu$Plots$"Add to graphic"$"curve"$icon = "curve"
pmg.menu$Plots$"Add to graphic"$"rug"$handler = 
  function(h,...) pmg.gw(rug.list)
pmg.menu$Plots$"Add to graphic"$"title"$handler = 
  function(h,...) pmg.gw(add.title.list)


pmg.menu$Plots$"Teaching demos"$handler =
  function(h,...) pmg.teachingDemos()
###
### tests
pmg.menu$Tests$"Dynamic tests"$handler = function(h,...) {
  dTestsDialog()
}
pmg.menu$Tests$"Dynamic tests"$icon="execute"
##
pmg.menu$Tests$centers$"t.test"$handler =
  function(h,...) pmg.gw(t.test.list)
pmg.menu$Tests$centers$"t.test (summarized data)"$handler =
  function(h,...) pmg.gw(t.test.summaries.list)
pmg.menu$Tests$centers$"wilcox.test"$handler =
  function(h,...) pmg.gw(wilcox.test.list)
pmg.menu$Tests$centers$"oneway.test"$handler =
  function(h,...) pmg.gw(oneway.test.list)
pmg.menu$Tests$centers$"kruskal.test"$handler =
  function(h,...) pmg.gw(kruskal.test.list)
#
pmg.menu$Tests$scales$"var.test"$handler =
  function(h,...) pmg.gw(var.test.list)
pmg.menu$Tests$scales$"ansari.test"$handler =
  function(h,...) pmg.gw(ansari.test.list)
pmg.menu$Tests$scales$"bartlett.test"$handler =
  function(h,...) pmg.gw(bartlett.test.list)
pmg.menu$Tests$scales$"fligner.test"$handler =
  function(h,...) pmg.gw(fligner.test.list)
#
pmg.menu$Tests$shape$"ks.test"$handler =
  function(h,...) pmg.gw(ks.test.list)
pmg.menu$Tests$shape$"shapiro.test"$handler =
  function(h,...) pmg.gw(shapiro.test.list)
#
pmg.menu$Tests$proportion$"prop.test"$handler =
  function(h,...) pmg.gw(prop.test.list)
pmg.menu$Tests$proportion$"binom.test"$handler =
  function(h,...) pmg.gw(binom.test.list)
#
pmg.menu$Tests$counts$"chisq.test"$handler =
  function(h,...) pmg.gw(chisq.test.list)
pmg.menu$Tests$counts$"mantelhaen.test"$handler =
  function(h,...) pmg.gw(mantelhaen.test.list)
pmg.menu$Tests$counts$"mcnemar.test"$handler =
  function(h,...) pmg.gw(mcnemar.test.list)
#
pmg.menu$Tests$correlation$"cor.test"$handler =
  function(h,...) pmg.gw(cor.test.list)
#

###
pmg.menu$Models$"Dynamic models"$handler = function(h,...) {
  dModelsDialog()
}
pmg.menu$Models$"Dynamic models"$icon="execute"
##
pmg.menu$Models$Regression$"lm"$handler =
  function(h,...) pmg.gw(lm.list)
pmg.menu$Models$Regression$"lqs"$handler =
  function(h,...) pmg.gw(lqs.list)
pmg.menu$Models$Regression$"glm"$handler =
  function(h,...) pmg.gw(glm.list)
#
pmg.menu$Models$ANOVA$"aov"$handler =
  function(h,...) pmg.gw(aov.list)
pmg.menu$Models$ANOVA$"anova"$handler =
  function(h,...) pmg.gw(anova.list)
#
pmg.menu$Models$"Mixed effects"$"gls"$handler =
  function(h,...) pmg.gw(gls.list)
pmg.menu$Models$"Mixed effects"$"lmList"$handler =
  function(h,...) pmg.gw(lmList.list)
pmg.menu$Models$"Mixed effects"$"lme"$handler =
  function(h,...) pmg.gw(lme.list)
#
pmg.menu$Models$Diagnostics$"plot.lm"$handler =
  function(h,...) pmg.gw(lm.diagnostics.list)
pmg.menu$Models$Diagnostics$"plot.lme"$handler =
  function(h,...) pmg.gw(plot.lme.diagnostics.list)
pmg.menu$Models$Diagnostics$"qqnorm.lme"$handler =
  function(h,...) pmg.gw(qqnorm.lme.diagnostics.list)
pmg.menu$Models$Diagnostics$"pairs.lme"$handler =
  function(h,...) pmg.gw(pairs.lme.diagnostics.list)
#
pmg.menu$Models$"Model selection"$"anova"$handler =
  function(h,...) pmg.gw(anova.list)
pmg.menu$Models$"Model selection"$"stepAIC"$handler =
  function(h,...) pmg.gw(stepAIC.list)



## help menu
## help menu is separate
## Question: do I want popup dialogs here or integrated in framework?
help.menu=list()
help.menu$Help$"About R"$handler =
  function(h,...) add(pmg.dialog.notebook,pmg.about.R(),label = "About R")
help.menu$Help$"About PMG"$handler =
  function(h,...) add(pmg.dialog.notebook,pmg.about(),label = "About P M G")
##pmg.about(container=pmgWC$new(v=TRUE))
help.menu$Help$"About PMG"$icon="about"
## help.menu$Help$"R helpbrowser"$handler =
##   function(h,...) {
##     if(is.null(pmg.helpBrowser.window) ||
##        is.invalid(pmg.helpBrowser.window)) {
##       assignInNamespace("pmg.helpBrowser.window", ghelpbrowser(),"pmg")
##     } else {
##       focus(pmg.helpBrowser.window) <- TRUE # will this work
##     }
##   }
help.menu$Help$"R Site Search"$handler = function(h,...) RSiteSearch.Dialog()
help.menu$Help$"View vignettes"$handler = function(h,...) viewVignettes.Dialog()
help.menu$Help$"View demos"$handler = function(h,...) viewDemos.Dialog()
help.menu$Help$"View P M G vignette"$handler = function(h,...) print(vignette("pmg",package="pmg"))
## help.menu$Help$"PMG manual"$handler =
##   function(h,...) vignette("pmg")
## help.menu$Help$"Help on topic..."$handler =
##   function(h,...) pmg.helpBrowser()
## help.menu$Help$"Run examples from..."$handler =
##   function(h,...) pmg.examplesBrowser()
## help.menu$Help$"Run demos from..."$handler =
##   function(h,...) pmg.demosBrowser()
## help.menu$Help$"Search web for answers..."$handler =
##   function(h, ...) pmg.RSiteSearch()



#pmg.gw = function(lst) {
#  widget = ggenericwidget(lst, container=NULL, cli=cli)
#  win = gwindow(lst$title, v=T)
#  group = ggroup(container=win)
#  gvarbrowser(container=group)
#  add(group, widget, expand=TRUE)
#}


##################################################
#win = gwindow("P M G", v=T)
#group = ggroup(horizontal=FALSE, container=win)
#m = gmenu(menu, container=group)
#cli = icli(container=group)
