
GetPathwayData<-function(pathway){

if(!exists("envData")) envData<-initializeMiRSEA()

if(pathway=="kegg"){
keggpathway<-get("kegg",envir=envData)
return(keggpathway)
}

if(pathway=="biocarta"){
biocartapathway<-get("biocarta",envir=envData)
return(biocartapathway)
}

if(pathway=="reactome"){
reactomepathway<-get("reactome",envir=envData)
return(reactomepathway)
}

}