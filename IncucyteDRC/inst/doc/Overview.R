## ----message=FALSE, echo=FALSE, warnings=FALSE---------------------------
library(IncucyteDRC)
# note - how to export as pdf rmarkdown::render('vignettes/Overview.Rmd', 'pdf_document', output_dir='~/Desktop')

## ------------------------------------------------------------------------
test_pm <- importPlatemapXML(system.file(file='extdata/example.PlateMap', package='IncucyteDRC'))
head(test_pm)

## ----fig.width=10, fig.height=6, dpi=70----------------------------------
plotPlatemap(test_pm)

## ------------------------------------------------------------------------
test_data <- importIncucyteData(system.file(file='extdata/example_data.txt', package='IncucyteDRC'), metric='pc')
names(test_data)
head(test_data$data)

## ----message=FALSE, warning=FALSE----------------------------------------
test_list <- splitIncucyteDRCPlateData(test_pm, test_data, group_columns='growthcondition')
test_idrcset <- test_list[[2]]
class(test_list)
class(test_idrcset)
names(test_idrcset)
length(test_list)
for (t in test_list) print(t$metadata)

## ----message=FALSE, warning=FALSE----------------------------------------
test_pm <- importPlatemapXML(system.file(file='extdata/example2.PlateMap', package='IncucyteDRC'))
test_pm$celltype <- 'T47D'  #modify contents of platemap
test_data <- importIncucyteData(system.file(file='extdata/example_data2.txt', package='IncucyteDRC'), metric='pc')
test_idrcset <- splitIncucyteDRCPlateData(test_pm, test_data, group_columns='celltype')
class(test_idrcset)

## ----message=FALSE, warning=FALSE----------------------------------------
test_idrcset <- fitGrowthCurvesGrouped(test_idrcset)
test_idrcset <- fitGrowthCurvesIndividual(test_idrcset)
class(test_idrcset)

## ----fig.width=6, fig.height=4-------------------------------------------
plotIncucyteDRCSet(test_idrcset, grouped=FALSE)
plotIncucyteDRCSet(test_idrcset, grouped=TRUE)


## ----fig.width=6, fig.height=4, message=FALSE, warning=FALSE-------------
test_idrcset <- calculateDRCData(test_idrcset, cut_time = 170)
plotIncucyteDRCSet(test_idrcset, grouped=FALSE)
exportDRCDataToDataFrame(test_idrcset)[1:5,]
names(test_idrcset)


## ------------------------------------------------------------------------
exportDRCDataToPRISM(test_idrcset)[,1:5]

## ------------------------------------------------------------------------
exportDRCDataToDataFrame(test_idrcset, include_control = TRUE)[1:5,]
exportDRCDataToDataFrame(test_idrcset, add_metadata = TRUE)[1:5,]

## ----fig.width=4.5, fig.height=3, message=FALSE, warning=FALSE-----------
test_idrcset <- fitDoseResponseCurve(test_idrcset, include_control = TRUE)
names(test_idrcset)
test_idrcset <- calculateEC50(test_idrcset)
names(test_idrcset)
exportEC50Data(test_idrcset)
plotDoseResponseCurve(test_idrcset, 'PDD00017252')


## ----fig.width=6, fig.height=4, message=FALSE, warning=FALSE-------------
test_idrcset <- calculateCutTimeForIDRCSet(test_idrcset, baseline_time = 5, no_doublings = 3.5, max_val = 80)
names(test_idrcset)
test_idrcset$cut_plot
test_idrcset$calculated_cut
test_idrcset <- calculateDRCData(test_idrcset)
plotIncucyteDRCSet(test_idrcset, grouped=TRUE)


## ----message=FALSE, warning=FALSE----------------------------------------
library(dplyr)
test_pm <- importPlatemapXML(system.file(file='extdata/example2.PlateMap', package='IncucyteDRC'))
test_pm$celltype <- 'T47D'  #modify contents of platemap
test_data <- importIncucyteData(system.file(file='extdata/example_data2.txt', package='IncucyteDRC'), metric='pc')
test_idrcset <- splitIncucyteDRCPlateData(test_pm, test_data, group_columns='celltype')
test_idrcset %>% 
    fitGrowthCurvesGrouped() %>% 
    fitGrowthCurvesIndividual() %>%
    calculateCutTimeForIDRCSet(baseline_time=24, no_doublings=4, max_val=80) %>%
    calculateDRCData() %>%
    fitDoseResponseCurve(include_control = TRUE) %>%
    calculateEC50() %>%
    exportEC50Data()

## ----message=FALSE, warning=FALSE----------------------------------------
test_pm <- importPlatemapXML(system.file(file='extdata/example.PlateMap', package='IncucyteDRC'))
test_data <- importIncucyteData(system.file(file='extdata/example_data.txt', package='IncucyteDRC'), metric='pc')
test_list <- splitIncucyteDRCPlateData(test_pm, test_data, group_columns='growthcondition') %>%
    lapply(fitGrowthCurvesGrouped) %>%
    lapply(fitGrowthCurvesIndividual) %>%
    lapply(calculateCutTimeForIDRCSet, baseline_time=24, no_doublings=4, max_val=80) %>%
    lapply(calculateDRCData) %>%
    lapply(fitDoseResponseCurve, include_control=TRUE) %>%
    lapply(calculateEC50)
    
#export EC50 data and combine
lapply(test_list, exportEC50Data, add_metadata=TRUE) %>% dplyr::bind_rows()


## ----eval=FALSE----------------------------------------------------------
#  shinyVisApp()

## ----eval=FALSE----------------------------------------------------------
#  #save code as app.R
#  library(IncucyteDRC)
#  library(shiny)
#  
#  options(warn=-1) #disable warnings
#  
#  #initiate the shiny app
#  shinyApp(
#      ui = IncucyteDRC::shinyVisUI(),
#      server = function(input, output) {
#          IncucyteDRC::shinyVisServer(input, output)
#      }
#  )

