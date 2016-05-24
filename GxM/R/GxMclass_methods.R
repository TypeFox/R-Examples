GxMclass = setClass(Class="GxMclass",
                    representation(loglikelihood="numeric", BIC="numeric", par="numeric", 
                                   hess="matrix", gradient="numeric",
                                   modelname="character",
                                   zeroset="character", closedform="logical", 
                                   K="numeric", coreNumber="numeric"));

setMethod("show", "GxMclass", definition=function(object) print.GxMclass(object) );


summaryGxMclass = setClass(Class="summaryGxMclass",
         representation(loglikelihood="numeric", BIC="numeric", par="numeric", stderr="numeric"));

setMethod("show", "summaryGxMclass", definition=function(object) print.summaryGxMclass(object) );
