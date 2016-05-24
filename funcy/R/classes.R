#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
#

setClassUnion("numNULL", c("numeric", "NULL"))
setClassUnion("listNULL", c("list", "NULL"))
setOldClass("proc_time", "proc_time") #convert S3 class


setClass("fctName", representation(
    name="character"
)
         )

##fpcCtrl control object for functional principal components analysis
setClass("fpcCtrl",
         representation(
             select="character",
             h1Dim="numeric",
             h2Dim="numeric",
             sm1Dim="character",
             sm2Dim="character",
             coeffsCalc="character",
             nrMaxTime="numeric",
             average="logical"
             ),
         prototype(
             select="automatic",
             sm1Dim="sm.regression",
             sm2Dim="sm.regression",
             coeffsCalc="estimate",
             nrMaxTime=15,
             average=TRUE
             ))

setAs("list", "fpcCtrl",
      function(from, to){
          z <- list2object(from, to)
          z@select <- match.arg(z@select,
                                c("automatic", "manual"))

          z@sm1Dim <- match.arg(z@sm1Dim, c("sm1", "sm.regression"))
          z@sm2Dim <- match.arg(z@sm2Dim, c("sm2", "sm.regression"))
        
          z@coeffsCalc <- match.arg(z@coeffsCalc,
                                c("estimate", "integrate"))
          z
      })


##**********************************************************
setClass("funcyCtrl",
         representation(
             baseType="character",
             dimBase="numeric",
             flexDim="logical",
             init="character",
             nrep="numeric",
             seed="numNULL",
             thd="numeric",
             redDim="numeric"
         ),
         prototype(
             baseType="splines",
             dimBase=4,
             flexDim=FALSE,
             init="kmeans",
             nrep=3,
             seed=2804,
             thd=0.01,
             redDim=2
         )
         )


setAs("list", "funcyCtrl",
      function(from, to){
          z <- list2object(from, to)
          z@baseType <- match.arg(z@baseType,
                                  c("eigenbasis", "splines",
                                    "exponential", "fourier",
                                    "power"))
          z@init <- match.arg(z@init,
                              c("kmeans", "random", "hclust"))
          z
      }
      )

setClass("funcyCtrlMbc",
         contains="funcyCtrl",
         representation(
             eps="numeric",
             maxit="numeric",
             hard="logical"     
         ),
         prototype(
             eps=0.01,
             maxit=10,
             hard=FALSE
         )
         )


setAs("list", "funcyCtrlMbc",
      function(from, to){
          z <- list2object(from, to)
      }
      )



###*********************************************************
setClass("sampleFuncy",
         representation(
             data="matrix",
             clusters="numeric",
             reg="logical"
         )
         )
              

setClass("funcyOut",
         representation(
             control="funcyCtrl",
             methodName="character",
             kOut="numeric",
             dimBaseOut="numeric",
             time="numeric",
             cluster="numeric",
             centers="matrix",
             props="numeric",
             dist2centers="matrix",
             cldist="matrix",
             calcTime="proc_time",
             plotParams="listNULL",
             correctCl="numeric"
         )
         )


setClass("funcyOut-iterSubspace",
         contains="funcyOut",
         representation(
             groupDimBase="numeric",
             prms="list",
             nrIter="numeric"
         )
         )

setClass("funcyOutMbc",
         contains="funcyOut-iterSubspace",
         representation(
             control="funcyCtrlMbc",
             probs="matrix",
             AIC="numeric",
             BIC="numeric",
             logLik="numeric"          
         )
         )


setClass("funcyOutMbc-fitfclust",
         contains="funcyOutMbc",
         representation(
             fit="list"                    
         )
    )

setClass("funcyOutMbc-fscm",
         contains="funcyOutMbc",
         representation(
             trends="list",
             location="data.frame"
         )
    )

##Out list**************************************************
setClass("funcyOutList",
         representation(
             call="call",
             models="list",
             data="matrix",
             timeNr="data.frame",
             reg="logical",
             k="numeric",
             methodName="character",
             allClusters="matrix",
             randIndex="matrix",
             votedCluster="numeric",
             accordance="numeric"
         ))



