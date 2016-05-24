

GetExampleData<-function(exampleData){

if(!exists("envData")) envData<-initialize()




if (exampleData=="code")
{
code<-get("code",envir=envData)

return(code)
}



if(exampleData=="compound")
{

compound<-get("compound",envir=envData)

return(compound)
}


if(exampleData=="g2")
{

g2<-get("g2",envir=envData)

return(g2)
}

if(exampleData=="gene")
{

gene<-get("gene",envir=envData)

return(gene)
}



if(exampleData=="gene2path")
{

gene2path<-get("gene2path",envir=envData)

return(gene2path)
}

if(exampleData=="gene2symbol")
{

gene2symbol<-get("gene2symbol",envir=envData)

return(gene2symbol)
}

if(exampleData=="GeneExp")
{

GeneExp<-get("GeneExp",envir=envData)

return(GeneExp)
}

if(exampleData=="keggGene2gene")
{

keggGene2gene<-get("keggGene2gene",envir=envData)

return(keggGene2gene)
}


if(exampleData=="lnc2Name")
{

lnc2Name<-get("lnc2Name",envir=envData)

return(lnc2Name)
}

if(exampleData=="lncBackground")
{

lncBackground<-get("lncBackground",envir=envData)

return(lncBackground)
}


if(exampleData=="LncExp")
{

LncExp<-get("LncExp",envir=envData)

return(LncExp)
}


if(exampleData=="LncGenePairs")
{

LncGenePairs<-get("LncGenePairs",envir=envData)

return(LncGenePairs)
}


if(exampleData=="mart")
{

mart<-get("mart",envir=envData)

return(mart)
}


if(exampleData=="nocode")
{

nocode<-get("nocode",envir=envData)

return(nocode)
}


if(exampleData=="pp")
{

pp<-get("pp",envir=envData)

return(pp)
}


if(exampleData=="resultT")
{

resultT<-get("resultT",envir=envData)

return(resultT)
}


if(exampleData=="sub")
{

sub<-get("sub",envir=envData)

return(sub)
}


if(exampleData=="SubcodeLncResult")
{

SubcodeLncResult<-get("SubcodeLncResult",envir=envData)

return(SubcodeLncResult)
}

if(exampleData=="geneLnc")
{

geneLnc<-get("geneLnc",envir=envData)

return(geneLnc)
}


}