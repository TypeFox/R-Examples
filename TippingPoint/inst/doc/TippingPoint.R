## ----set, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------
# Change the width of html file
options(width = 150)


## ----data-------------------------------------------------------------------------------------------------------------------------------------------

# If you haven't install the package, you can download it from cran

# install.packages("TipingPoint")

library(TippingPoint)

# Load the dataset

data(tippingdata)

# Show the first 6 rows of the data

head(tippingdata)


## ----basic plot, fig.width=8,fig.height=6-----------------------------------------------------------------------------------------------------------

## for binary outcome


# Using `estimate`
TippingPoint(outcome=tippingdata$binary,
             treat= tippingdata$treat,group.infor=TRUE,
             plot.type = "estimate",ind.values = TRUE,
             impValuesT  = NA,  impValuesC = NA,
             summary.type = "density", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(0.38,0.4), HistMeanC =  c(0.2,0.55))

# Using `p.value` with formula class
TippingPoint(binary~treat, data=tippingdata,
             plot.type = "p.value",ind.values = TRUE,
             impValuesT  = NA,  impValuesC = NA,
             summary.type = "density", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(0.38,0.4), HistMeanC =  c(0.2,0.55))

# Using `both` 
TippingPoint(outcome=tippingdata$binary,treat= tippingdata$treat,
             plot.type = "both",ind.values = TRUE,
             impValuesT  = NA,  impValuesC = NA,
             summary.type = "density", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(0.38,0.4), HistMeanC =  c(0.2,0.55))



# for continuous outcome
TippingPoint(continuous~treat, data=tippingdata,
             group.infor=TRUE, plot.type = "estimate",ind.values = TRUE,
             impValuesT  = NA,  impValuesC = NA,
             summary.type = "density", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(120), HistMeanC =  c(131,137))

TippingPoint(outcome=tippingdata$continuous,treat= tippingdata$treat,
             plot.type = "p.value",ind.values = TRUE,
             impValuesT  = NA,  impValuesC = NA,
             summary.type = "density", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(120), HistMeanC =  c(131,137))

TippingPoint(outcome=tippingdata$continuous,treat= tippingdata$treat,
             plot.type = "both",ind.values = TRUE,
             impValuesT  = NA,  impValuesC = NA,
             summary.type = "density", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(120), HistMeanC =  c(131,137))


## ----impute plot, fig.width=8,fig.height=6-----

# Load the imputed dataset

data(imputedata)

# Show the first 6 rows of the data

head(imputedata)


## for binary outcome

TippingPoint(outcome=tippingdata$binary,
             treat= tippingdata$treat, group.infor=TRUE,
             plot.type = "estimate",ind.values = TRUE,
             impValuesT  = imputedata[,c("MAR_T2","MCAR_T2")],  
             impValuesC = imputedata[,c("MAR_C2","MCAR_C2")],
             summary.type = "density", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(0.38,0.4), HistMeanC =  c(0.2,0.55))


# User-defined colors
TippingPoint(outcome=tippingdata$binary,treat= tippingdata$treat,
             plot.type = "p.value",ind.values = TRUE,
             impValuesT  = imputedata[,c("MAR_T2","MCAR_T2")],  
             impValuesC = imputedata[,c("MAR_C2","MCAR_C2")],
             impValuesColor = RColorBrewer::brewer.pal(8,"Accent")[5:6],
             summary.type = "credible.region", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(0.38,0.4), HistMeanC =  c(0.2,0.55))

# Using `point.size` and `point.shape` to control the points.
TippingPoint(outcome=tippingdata$binary,treat= tippingdata$treat,
             plot.type = "both",ind.values = TRUE,
             impValuesT  = imputedata[,c("MAR_T2","MCAR_T2")],  
             impValuesC = imputedata[,c("MAR_C2","MCAR_C2")],
             impValuesColor =c("red","blue"),
             point.size=0.8,point.shape = 15,
             summary.type = "convex.hull", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(0.38,0.4), HistMeanC =  c(0.2,0.55))



## for continuous outcome
TippingPoint(outcome=tippingdata$continuous,
             treat= tippingdata$treat, group.infor=TRUE,
             plot.type = "p.value",ind.values = TRUE,
             impValuesT  = imputedata[,c("MAR_T1","MCAR_T1")],  
             impValuesC = imputedata[,c("MAR_C1","MCAR_C1")],
             summary.type = "density", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(120), HistMeanC =  c(131,137))


# Using `point.size` and `point.shape` to control the points.
TippingPoint(outcome=tippingdata$continuous,treat= tippingdata$treat,
             plot.type = "p.value",ind.values = TRUE,
             impValuesT  = imputedata[,c("MAR_T1","MCAR_T1")],  
             impValuesC = imputedata[,c("MAR_C1","MCAR_C1")],
             point.size = 0.8, point.shape = 15,
             summary.type = "credible.region", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(120), HistMeanC =  c(131,137))

TippingPoint(outcome=tippingdata$continuous,treat= tippingdata$treat,
             plot.type = "both",ind.values = TRUE,
             impValuesT  = imputedata[,c("MAR_T1","MCAR_T1")],  
             impValuesC = imputedata[,c("MAR_C1","MCAR_C1")],
             summary.type = "convex.hull", alpha = 0.95, S=1.5, n.grid = 100,
             HistMeanT = c(120), HistMeanC =  c(131,137))



