## ----table1, echo = FALSE, results = 'asis'------------------------------
toy<-data.frame(Name=c('EC$eset','EC$genelist','EC$david','EC$genes','EC$process'),Description=c('Data frame of normalized expression values of brain and heart endothelial cells (3 replicates)','Data frame of differentially expressed genes (adjusted p-value < 0.05)','Data frame of results from a functional analysis of the differentially expressed genes performed with DAVID','Data frame of selected genes with logFC','Character vector of selected enriched biological processes'),Dimension=c('20644 x 7','2039 x 7','174 x 5','37 x 2','7'))
knitr::kable(toy, colnames=c('Name','Description','Dimension (row, col)'))

## ----glimpse, warning = FALSE, message = FALSE---------------------------
library(GOplot)
# Load the dataset
data(EC)
# Get a glimpse of the data format of the results of the functional analysis... 
head(EC$david)
# ...and of the data frame of selected genes
head(EC$genelist)

## ----circ_object, warning = FALSE, message = FALSE-----------------------
# Generate the plotting object
circ <- circle_dat(EC$david, EC$genelist)

## ----GOBar, warning = FALSE, message = FALSE, fig.width = 8.3, fig.height = 6----
# Generate a simple barplot
GOBar(subset(circ, category == 'BP'))

## ----GOBar2, eval = FALSE, warning = FALSE, message = FALSE--------------
#  # Facet the barplot according to the categories of the terms
#  GOBar(circ, display = 'multiple')

## ----GOBar3, eval = FALSE, warning = FALSE, message = FALSE--------------
#  # Facet the barplot, add a title and change the colour scale for the z-score
#  GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))

## ----GOBubble1, warning = FALSE, message = FALSE, fig.keep = 'none'------
# Generate the bubble plot with a label threshold of 3
GOBubble(circ, labels = 3)

## ----GOBubble2, warning = FALSE, message = FALSE, fig.keep = 'none'------
# Add a title, change the colour of the circles, facet the plot according to the categories and change the label threshold
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)

## ----GOBubble3, warning = FALSE, message = FALSE, fig.keep = 'none'------
# Colour the background according to the category
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)

## ----GOBubble4, warning = FALSE, message = FALSE, fig.keep = 'none', eval = FALSE----
#  # Reduce redundant terms with a gene overlap >= 0.75...
#  reduced_circ <- reduce_overlap(circ, overlap = 0.75)
#  # ...and plot it
#  GOBubble(reduced_circ, labels = 2.8)

## ----GOCircle1, warning = FALSE, message = FALSE, fig.keep = 'none'------
# Generate a circular visualization of the results of gene- annotation enrichment analysis
GOCircle(circ)

## ----GOCircle2, eval = FALSE---------------------------------------------
#  # Generate a circular visualization of selected terms
#  IDs <- c('GO:0007507', 'GO:0001568', 'GO:0001944', 'GO:0048729', 'GO:0048514', 'GO:0005886', 'GO:0008092', 'GO:0008047')
#  GOCircle(circ, nsub = IDs)

## ----GOCircle3, eval = FALSE---------------------------------------------
#  # Generate a circular visualization for 10 terms
#  GOCircle(circ, nsub = 10)

## ----GOChord1, warning = FALSE, message = FALSE--------------------------
# Define a list of genes which you think are interesting to look at. The item EC$genes of the toy 
# sample contains the data frame of selected genes and their logFC. Have a look...
head(EC$genes)
# Since we have a lot of significantly enriched processes we selected some specific ones (EC$process)
EC$process
# Now it is time to generate the binary matrix
chord <- chord_dat(circ, EC$genes, EC$process)
head(chord)

## ----GOChord2, eval=FALSE, warning = FALSE, message = FALSE--------------
#  # Generate the matrix with a list of selected genes
#  chord <- chord_dat(data = circ, genes = EC$genes)
#  # Generate the matrix with selected processes
#  chord <- chord_dat(data = circ, process = EC$process)

## ----GOChord3, warning = FALSE, message = FALSE, fig.keep = 'none'-------
# Create the plot
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)

## ----GOChord4, warning = FALSE, message = FALSE, fig.keep = 'none'-------
# Display only genes which are assigned to at least three processes
GOChord(chord, limit = c(3, 0), gene.order = 'logFC')

## ----GOHeat1, warning = FALSE, message = FALSE, fig.keep = 'none'--------
# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord[,-8], nlfc = 0)

## ----GOHeat2, warning = FALSE, message = FALSE, fig.keep = 'none'--------
# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))

## ----GOCluster, warning=FALSE, eval=FALSE, message=FALSE, fig.keep='none'----
#  GOCluster(circ, EC$process, clust.by = 'logFC', term.width = 2)

## ----GOCluster2, warning=FALSE, eval=FALSE, message=FALSE, fig.keep='none'----
#  GOCluster(circ, EC$process, clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))

## ----GOVenn, warning=FALSE, message=FALSE, fig.keep='none'---------------
l1 <- subset(circ, term == 'heart development', c(genes,logFC))
l2 <- subset(circ, term == 'plasma membrane', c(genes,logFC))
l3 <- subset(circ, term == 'tissue morphogenesis', c(genes,logFC))
GOVenn(l1,l2,l3, label = c('heart development', 'plasma membrane', 'tissue morphogenesis'))

