## ---- eval=FALSE---------------------------------------------------------
#  install.packages("MonoPhy")

## ---- eval=FALSE---------------------------------------------------------
#  # 1. Install package 'devtools' (if not already installed):
#  install.packages("devtools")
#  
#  # 2. Load 'devtools' and temporarily install and load 'MonoPhy' in developers mode:
#  library(devtools)
#  dev_mode(on=T)
#  install_github("oschwery/MonoPhy")  # install the package from GitHub
#  library(MonoPhy)  # load the package
#  
#  #3. Leave developers mode after done trying out 'MonoPhy':
#  dev_mode(on=F)
#  # the package will not stay on your system permanently like that, which make sense in case of the
#  # development version, as opposed to the stable version

## ---- results="hide", warning=FALSE, message=FALSE-----------------------
library(MonoPhy)

## ---- results="hide"-----------------------------------------------------
#load data
data(Ericactree)
data(Ericactribes)
data(Ericacsubfams)

## ------------------------------------------------------------------------
#check data
Ericactree
head(Ericactribes)
head(Ericacsubfams)

## ---- eval=FALSE---------------------------------------------------------
#  phy <- read.tree(file="/your_path/your_phylogeny.tre")

## ---- eval=FALSE---------------------------------------------------------
#  your_clades <- read.csv(file="/your_path/your_taxonomy_file.csv", header=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  #check data
#  phy
#  head(your_clades)

## ------------------------------------------------------------------------
solution0 <- AssessMonophyly(Ericactree)

## ------------------------------------------------------------------------
solution1 <- AssessMonophyly(Ericactree, Ericacsubfams)

## ---- eval=FALSE---------------------------------------------------------
#  AssessMonophyly <-
#  function(tree, taxonomy=NULL, verbosity=5, outliercheck=TRUE, outlierlevel=0.5,
#  taxizelevel= NULL, taxizedb='ncbi', taxizepref='ncbi', taxask=FALSE, taxverbose=FALSE)
#  

## ------------------------------------------------------------------------
GetSummaryMonophyly(solution0)  # pull out summary table

## ------------------------------------------------------------------------
GetResultMonophyly(solution1, taxlevels=2)  # pull out summary only for taxlevel 2, subfamilies in this case

## ------------------------------------------------------------------------
GetIntruderTaxa(solution0, taxa="Erica")  # get list of genera that cause monophyly-issues for Erica

## ------------------------------------------------------------------------
GetIntruderTips(solution1, taxa="Ericoideae", taxlevels=2)  # get list of species causing monophyly-issues for the subfamily Ericoideae

## ------------------------------------------------------------------------
GetAncNodes(solution1, taxa = NULL, taxlevels='ALL')

## ------------------------------------------------------------------------
GetOutlierTaxa(solution0, taxlevels='ALL')
GetOutlierTips(solution1, taxa = NULL, taxlevels='ALL')

## ---- fig.width=7, fig.height=14-----------------------------------------
PlotMonophyly(solution0, Ericactree, plot.type='monophyly', ladderize=TRUE, cex=0.5)

## ---- fig.width=7, fig.height=14-----------------------------------------
PlotMonophyly(solution1, Ericactree, taxlevels=2, plot.type='taxonomy', cex=0.5)

## ---- fig.width=7, fig.height=14-----------------------------------------
PlotMonophyly(solution1, Ericactree, taxlevels=2, plot.type='monoVStax', cex=0.4, label.offset=18)

## ---- fig.width=7, fig.height=7------------------------------------------
PlotMonophyly(solution1, Ericactree, taxlevels=1, plot.type='intruders', monocoll=TRUE, cex=0.5, label.offset=5)

## ---- eval=FALSE---------------------------------------------------------
#  PlotMonophyly(solution0, Ericactree, taxlevels=1, plot.type='monophyly', monocoll=TRUE,
#  ladderize=TRUE, PDF=TRUE, PDF_filename='Monophylyplot_Ericac_Example.pdf')

## ---- echo=FALSE, fig.width=7, fig.height=9------------------------------
PlotMonophyly(solution0, Ericactree, taxlevels=1, plot.type='monophyly', monocoll=TRUE, ladderize=TRUE, cex=0.5, label.offset=5)

## ---- eval=FALSE---------------------------------------------------------
#  PlotMonophyly(solution1, Ericactree, taxlevels=2, plot.type='intruders', ladderize=TRUE, PDF=TRUE,
#  PDF_filename='Intruderplot_Ericac_Example.pdf', type="fan")

## ---- echo=FALSE, fig.width=7.5, fig.height=7.5--------------------------
PlotMonophyly(solution1, Ericactree, taxlevels=2, plot.type='intruders', ladderize=TRUE, type="fan", cex=0.4)

