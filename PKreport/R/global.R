#############################################################################################
## File: global.R
## Author: Xiaoyong Sun
## Date: 10/12/2009
## Goal: global variables
## Notes:
##      -
#############################################################################################

.pkplot <- local(
{
    .globalConfig <- list(
                        save.format="png",
                        width = 480, height = 480,
                        package=2,
                        univar.dir="univar",
                        bivar.dir="bivar",
                        ind.dir="ind",
                        gof.dir="gof",
                        struct.dir="struct",
                        resid.dir="resid",
                        para.dir="para",
                        cov.dir="cov",
                        eta.dir="eta"
                    )
                    
    .histgraph <- list(
                      lattice = list(type=c("count"), layout=c(1,1)),
                      ggplot = list(geom=c("histogram"), layout=c(1,1)),
                      others = list(ind.layout=c(5,5))
                  )
    .scattergraph <- list(
                      lattice = list(type=c("p", "smooth"), layout=c(1,1)),
                      ggplot = list(geom=c("point", "smooth"), se=FALSE, layout=c(1,1)),
                      others = list(ind.layout=c(5,5))
                  )

    .term <- list(
                    TIME="Time", CONC="conc", ID="Subject"
                  )
                  
    .pkcode <- list()
    .pkcodenote <- list()
    .pkdata <- NULL
    .figno <- 0
    .dataname <- NULL
    
    list(
            getGlobalConfig = function() return(.globalConfig),
            getConfig = function(argsname) return(.globalConfig[[argsname]]),
            setGlobalConfig = function(argsname, value)
            {
                .globalConfig[[argsname]] <<- value
            },
            
            getHistGraph = function(argsname) return(.histgraph[[argsname]]),
            setHistGraph = function(one.list, argsname)
            {
                .histgraph[[argsname]] <<- one.list
            },
            getScatterGraph = function(argsname) return(.scattergraph[[argsname]]),
            setScatterGraph = function(one.list, argsname)
            {
                .scattergraph[[argsname]] <<- one.list
            },
            getTerm = function() return(.term),
            setTerm = function(term.list)
            {
                .term <<- term.list
            },
            
            getPKCode = function(i) return(.pkcode[[i]]),
            getPKCodeLen = function(i) return(length(.pkcode)),
            setPKCode = function(newlist)
            {
                newlen <- length(.pkcode)
                .pkcode[[newlen+1]] <<- newlist
            },
            cleanPKCode = function()
            {
                .pkcode <<- list()
            },
            
            getPKCodeNote = function(i) return(.pkcodenote[[i]]),
            getAllPKCodeNote = function(i) return(.pkcodenote),
            setPKCodeNote = function(newlist)
            {
                newlen <- length(.pkcodenote)
                .pkcodenote[[newlen+1]] <<- newlist
            },
            cleanPKCodeNote = function()
            {
                .pkcodenote <<- list()
            },
            
            getPKData = function() return(.pkdata),
            setPKData = function(dataset)
            {
                .pkdata <<- dataset
            },
            
            getDataName = function() return(.dataname),
            setDataName = function(dataset)
            {
                .dataname <<- dataset
            },
            
            getFigNo = function() return(.figno),
            setFigNo = function(no)
            {
                .figno <<- no
            }
            
        )

})

