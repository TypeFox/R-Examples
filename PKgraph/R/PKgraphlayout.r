#############################################################################################
## Project: PKgraph
## File: PKgraphlayout.R
## Author: Xiaoyong Sun
## Date: 08/19/2009
## Goal: PKgraph
##        - interface
## Notes:
#############################################################################################



################################################################################
## Main layout: menubar, toolbar
################################################################################
mbl = list(

    Project = list(

      
        "Open data" = list(handler= openDataHandler, icon="open"),
        "Save data" = list(handler=saveHandler, icon="save-as"),
        "Save workspace" = list(handler=saveProjectHandler, icon="save"),
        Sep = list(separator=TRUE),
        "Clean data" = list(handler= cleanDataHandler, icon="open"),
        "Restore old workspace" = list(handler=restoreHandler, icon="revert-to-saved"),
        "Exit" =  list(handler=exitHandler, icon="quit")
    ),
    
    Configure = list(
        "Set data type" = list(handler= function(h,...) ggobiPlotType()),
        "Set working directory" = list(handler= configDirHandler),
        "Set saving format" = list(handler= configFormatHandler),
        "Set figure configuration" = list(handler= configFigureHandler)

    ),

    "Data management" = list(

        Subset = list(handler= function(h,...) pk.subset()),
        "Factor" = list(handler= function(h,...) pk.factor())
    ),

    "Exploratory data analysis" = list(
        Univariates = list(handler= function(h,...) pk.uni.menu()),
        Bivariates = list(handler=function(h,...) pk.bi.menu()),
        #Trivariates = list(handler=function(h,...) pk.tri.menu()),
        Sep = list(separator=TRUE),
        #Timeplots = list(handler=function(h,...) pk.time.menu()),
        "Parallel coordinates plot" = list(handler=function(h,...) pk.para.menu()),
        #Heatmap = list(handler=function(h,...) pk.heat.menu()) ,
        "Scatterplot matrix" = list(handler=function(h,...) pk.matrix.menu())
    ),
    "PK Models" = list(
        "Configure model result" = list(handler= configDataHandler),
         Sep = list(separator=TRUE),
        "Individual plots" = list(handler= function(h,...) pk.model.ind()),
        "Basic goodness of fit plots" = list(handler= function(h,...) pk.model.gof()),
        "Parameters" = list(handler= function(h,...) pk.model.para()),
        "Random effects" = list(handler= function(h,...) pk.model.random()),
        "Structural model" = list(handler= function(h,...) pk.model.struct()),
        "Residual error model" = list(handler= function(h,...) pk.model.resid()),
        "Covariate model" = list(handler= function(h,...) pk.model.cov())

        
    ),
    "Model validation" = list(
        "Influence analysis summary (PsN)" = list(handler= function(h,...) pk.outlier.psn()),
        "Visualization for influence analysis" = list(handler= function(h,...) pk.outlier.vis()),
        "Bootstrap summary (PsN)" = list(handler= function(h,...) pk.boot.sum()) ,
        "Visualization for bootstrap" = list(handler= function(h,...) pk.boot.vis()) #,

        #"Simulation" = list(handler= function(h,...) pk.sim()) #,
        #"Numerical predicative check" = list(handler= gofHandler),
        #"Visual predicative check" = list(handler= gofHandler)
    ),
    "Model comparison" = list(
        "Select datasets" = list(handler= function(h,...) pk.com.config()),
        "Configure mapping" = list(handler= function(h,...) pk.com.map()),
         Sep = list(separator=TRUE),
        "Histogram comparison" = list(handler=function(h,...) pk.com.hist()),
        "Scatter plot comparison" = list(handler=function(h,...) pk.com.scatter()),
        "Transform comparison" = list(handler=function(h,...) pk.com.time())
    ),
    "Interactive diagnostics" = list(
        "Select datasets" = list(handler= function(h,...) ggobi.data()),
        "Configure mapping" = list(handler= function(h,...) ggobi.map()),
         Sep = list(separator=TRUE),
        "Diagnostics" = list(handler= function(h,...) ggobi.diagose())
    )

    ##,
    #Report = list(
     #   "Report for exploratory Data Analysis" = list(handler= gofHandler),
     #   "Report for PK models" = list(handler= gofHandler),
     #   "Report for model validation" = list(handler= gofHandler),
      #  "Report for model comparison" = list(handler= gofHandler)
    #)


)

## tool bar
tbl = list(
    open = list(handler=openDataHandler, icon="open"),
    preferences=list(handler=configFormatHandler, icon="preferences"),
    subset=list(handler=function(h,...) pk.subset(), icon="subset"),
    clean = list(handler=cleanDataHandler, icon="clear"),    
     save = list(handler=saveHandler, icon="save"),
     help = list(handler=helpHandler, icon="help"),
    quit = list(handler=exitHandler, icon="quit")

)
