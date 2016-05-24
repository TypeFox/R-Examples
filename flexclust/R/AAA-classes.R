#
#  Copyright (C) 2005-2009 Friedrich Leisch
#  $Id: AAA-classes.R 14 2013-07-02 09:56:24Z leisch $
#

setClass("flexclustControl",
         representation(iter.max="numeric",
                        tolerance="numeric",
                        verbose="numeric",
                        classify="character",
                        initcent="character",
                        gamma="numeric",          # for kcca
                        simann="numeric",         # for kcca
                        ntry="numeric",           # for qtclust
                        min.size="numeric",       # for qtclust and ockc
                        subsampling="numeric"     # for ockc
                        ),
         prototype(iter.max=200,
                   tolerance=10e-7,
                   verbose=0,
                   classify="auto",
                   initcent="randomcent",
                   gamma=1,
                   simann=c(0.3, 0.95, 10),
                   ntry=5,
                   min.size=2,
                   subsampling=1))


setAs("list", "flexclustControl",
function(from, to){
    z <- list2object(from, to)
    z@classify <- match.arg(z@classify,
                            c("auto", "weighted", "hard", "simann"))
    z
})

setAs("NULL", "flexclustControl",
function(from, to){
    new(to)
})

###**********************************************************

setClass("cclustControl",
         contains="flexclustControl",
         representation(pol.rate="numeric",
                        exp.rate="numeric",
                        ng.rate="numeric",
                        method="character"),
         prototype(pol.rate=c(1,0),
                   exp.rate=c(0.1, 0.0001),
                   ng.rate=c(0.5, 0.005, 10, 0.01),
                   method="polynomial"))


setAs("list", "cclustControl",
function(from, to){
    z <- list2object(from, to)
    z@method <- match.arg(z@method,
                          c("polynomial", "exponential"))
    z
})

setAs("NULL", "cclustControl",
function(from, to){
    new(to)
})


###**********************************************************
###**********************************************************

setClass("kccaFamily",
         representation(name="character",
                        dist="function",
                        cent="function",
                        allcent="function",
                        wcent="function",
                        weighted="logical",
                        cluster="function",
                        preproc="function",
                        groupFun="function"),
         prototype(weighted=FALSE,
                   preproc=function(x) x))



###**********************************************************

setClass("flexclust",
         representation(k="integer",
                        cluster="integer",
                        iter="integer",
                        converged="logical",
                        clusinfo="data.frame",
                        index="numeric",
                        call="call",
                        control="flexclustControl",
                        data="ModelEnv"))

setClass("kccasimple",
         contains="flexclust",
         representation(centers="ANY",
                        family="kccaFamily",
                        cldist="matrix"))       

setClass("kcca",
         contains="kccasimple",
         representation(second="integer",
                        xrange="ANY",           # range of all data
                        xcent="ANY",            # centroid of all data
                        totaldist="numeric",    # total dist data<->xcent
                        clsim="matrix"))

## Achtung: Bei Aenderung der Klassen auch dorelabel() entsprechend anpassen

                        
###**********************************************************

setClass("stepFlexclust",
         representation(models="list",
                        k="integer",
                        nrep="integer",
                        call="call",
                        xcent="ANY",            # centroid of all data
                        totaldist="numeric"     # total dist data<->xcent
                        ))

setClass("ockc",
         contains = "stepFlexclust",
         representation(order = "integer"))

         
###**********************************************************

setClass("bootFlexclust",
         representation(k="integer",
                        centers1="list",
                        centers2="list",
                        cluster1="array",
                        cluster2="array",
                        index1="matrix",
                        index2="matrix",
                        rand="matrix",
                        call="call"))

setClass("resampleFlexclust",
         representation(k="integer",
                        centers1="list",
                        centers2="list",
                        cluster1="list",
                        cluster2="list",
                        index1="list",
                        index2="list",
                        indextest="list",
                        validation="array",
                        call="call"))

setClass("resampleScheme",
         representation(traintest="function",
                        validate="function",
                        valname="character"))



###**********************************************************

setOldClass("trellis")
setClass("propBarchart",
         representation(chart="trellis",
                        gprop="matrix",
                        tprop="numeric",
                        p.value="numeric",
                        table="matrix",
                        call="call"))
         
###**********************************************************




