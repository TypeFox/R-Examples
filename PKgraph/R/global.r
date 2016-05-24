#############################################################################################
## Project: PKgraph
## File: global.R
## Author: Xiaoyong Sun
## Date: 08/19/2009
## Goal: PKgraph
##        - interface
## Notes:
#############################################################################################

.pk <- local(
{
    .term <- data.frame()
    .datasets <- list()
    #.dataSpeialPlot <- list()
    .validateData <- data.frame()

    ## for save...splom, abs...
    .dataSpecialPlot <- list()
    ## multiple page saving
    .dataLayoutPlot <- NULL
    
    #.datasetsno <- 0
    .dataType <- list()

    ## model comparison
    .comDataName <- NULL
    .comMap <- data.frame()
    
    ## for interactive
    .itDataName <- NULL
    .itMap <- NULL
    
    # 600
    .subheight <- 600
    .subwidth <- 600*1.6

    .pkcode <- list()
    .pkcodenote <- list()
    
    .pkggobi <- list()
    
    # saving format
    .saveFormat <- list()
    .figureConfig <- list(col="royalblue")

    # ggobi data type
    .ggobiPlotType <- list()

    
    list(
            getTerm = function() return(.term),
            setTerm = function(term.df)
            {
                .term <<- term.df
            },

            getValidateData = function() return(.validateData),
            setValidateData = function(vdata) .validateData <<- vdata,

            ## for abs plot, get special data set instead of default data: getCurrentData()
            getNameDataSpecialPlot = function() return(names(.dataSpecialPlot)),
            getDataSpecialPlot = function(i) return(.dataSpecialPlot[[i]]),
            setDataSpecialPlot = function(tdata, tname) .dataSpecialPlot[[tname]] <<- tdata,
            cleanDataSpecialPlot = function() .dataSpecialPlot <<- list(),

            ## for multiple layout
            # cleanDataLayoutPlot always put below $show... syntax set.
            getDataLayoutPlot = function() return(.dataLayoutPlot),
            setDataLayoutPlot = function(tdata) .dataLayoutPlot <<- c(.dataLayoutPlot, tdata),
            cleanDataLayoutPlot = function() .dataLayoutPlot <<- NULL,

            getDatasets = function() return(.datasets),
            setDatasets = function(dataset, dataname)
            {
                .datasets[[dataname]] <<- dataset
                #names(.datasets)[[thisDatano]] <<- thisDatano
                #.datasetsno <<- .datasetsno + 1
            },

            getCurrentData = function(currentNo)
            {
                if (missing(currentNo))
                {
                    if(length(.datasets)!= 0) 
		{
		    return(.datasets[[length(.datasets)]])
		}
                    else return(NULL)
                }
                else
                {
                    #if(currentNo > 0 && currentNo <= length(.datasets))
                    #{
                       return(.datasets[[currentNo]])
                    #}
                    #else return(NULL)
                }

            },
            getTotalDataLen = function() return(length(.datasets)),
            # setTotalDataLen = function(thisDatano) .datasetsno <<- thisDatano,
            getTotalDataName = function() return(names(.datasets)),

            getCurrentDataType = function(currentNo)
            {
                if (missing(currentNo))
                {
                    if(length(.dataType)!= 0) return(.dataType[[length(.dataType)]])
                    else return(NULL)
                }
                else
                {
                    if(currentNo > 0 && currentNo <= length(.dataType))
                    {
                       return(.dataType[[currentNo]])
                    }
                    else return(NULL)
                }

            },
            setCurrentDataType = function(thisDataType, dataname) .dataType[[dataname]] <<- thisDataType,

            getItDataName = function() return(.itDataName),
            setItDataName = function(itname) .itDataName <<- itname,
            
            getItMap = function() return(.itMap),
            setItMap = function(key) .itMap <<- key,

            getComDataName = function() return(.comDataName),
            setComDataName = function(comname) .comDataName <<- comname,
            
            getComMap = function() return(.comMap),
            setComMap = function(key.df) .comMap <<- key.df,

            getSubHeight = function() return(.subheight),
            getSubWidth = function() return(.subwidth),

            getPKCode = function(i) return(.pkcode[[i]]),
            getPKCodeLen = function(i) return(length(.pkcode)),
            setPKCode = function(newlist)
            {
                 newlen <- length(.pkcode)
                 .pkcode[[newlen+1]] <<- newlist
            },
            cleanPKCode = function() .pkcode <<- list(),
 
            getPKCodeNote = function(i) return(.pkcodenote[[i]]),
            getAllPKCodeNote = function(i) return(.pkcodenote),
            setPKCodeNote = function(newlist)
            {
                 newlen <- length(.pkcodenote)
                 .pkcodenote[[newlen+1]] <<- newlist
            },
            cleanPKCodeNote = function() .pkcodenote <<- list(),
            {

            },
            
            # set default x, y for ggobi. a list of x,y name
            setPKGGobi = function(newxy)
            {
                newlen <- length(.pkggobi)
                .pkggobi[[newlen+1]] <<- newxy
            },
            getPKGGobi = function(i) return(.pkggobi[[i]]),
            cleanPKGGobi = function() .pkggobi <<- list(),
            
            setSaveFormat = function(newformat) .saveFormat <<- newformat,
            getSaveFormat = function() return(.saveFormat),
            
            setFigConfig = function(newconfig) .figureConfig <<- newconfig,
            getFigConfig = function() return(.figureConfig),

            ## ggobi time series plot requirement
            getGGobiPlotType = function(currentNo)   
            {
	      

                if (missing(currentNo))
                {
                    if(length(.ggobiPlotType)!= 0) return(.ggobiPlotType[[length(.ggobiPlotType)]])
                    else return(NULL)
                }
                else
                {
                      if (currentNo > length(.ggobiPlotType)) return(NULL) 
                      else return(.ggobiPlotType[[currentNo]])
                }

            },

            setGGobiPlotType = function(typelist, dataname)  
            {
                .ggobiPlotType[[dataname]] <<- typelist
            },
            
            cleanAll = function()
            {
                .term <<- data.frame()
                .datasets <<- list()
                .dataType <<- list()    
                .ggobiPlotType <<- list()    
                    
                .validateData <<- data.frame()
            
                .dataSpecialPlot <<- list()
                .dataLayoutPlot <<- NULL
            
                ## model comparison
                .comDataName <<- NULL
                .comMap <<- data.frame()
                
                ## for interactive
                .itDataName <<- NULL
                .itMap <<- NULL 
                
                .pkcode <<- list()
                .pkcodenote <<- list()
                
                .pkggobi <<- list()                               
            }            

    )
    

})

## mainGUI
PKW = NULL
pmg.dialog.notebook = NULL
pmg.dialog.notebook2 = NULL
pmg.statusBar=NULL
pk.dirname = NULL
pk.dir = NULL

## current data
global.data <- NULL

## for specific data type
requiredDataType.PKmodel <- "PK data"
modelComType <- "ModelComparison"
