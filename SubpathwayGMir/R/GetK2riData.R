GetK2riData <-
function(K2riData){

if(!exists("k2ri")) k2ri<-initializeK2ri()
if (K2riData=="expMir2Tar")
{
expMir2Tar <- get("expMir2Tar",envir=k2ri)
return(expMir2Tar)
}

if (K2riData=="miRNA2Org")
{
miRNA2Org <- get("miRNA2Org",envir=k2ri)
return(miRNA2Org)
}

if (length(grep("MetabolicGEGEEMGraph",K2riData))>0)
{

MetabolicGEGEEMGraph <- get(K2riData,envir=k2ri)
return(MetabolicGEGEEMGraph)
}

if (length(grep("MetabolicGEGEUEMGraph",K2riData))>0)
{
MetabolicGEGEUEMGraph <- get(K2riData,envir=k2ri)
return(MetabolicGEGEUEMGraph)
}

if (K2riData=="BGMiRNA")
{
BGMiRNA <- get("BGMiRNA",envir=k2ri)
return(BGMiRNA)
}

if (K2riData=="BGGene")
{
BGGene <- get("BGGene",envir=k2ri)
return(BGGene)
}

if (K2riData=="gene2symbol")
{
gene2symbol <- get("gene2symbol",envir=k2ri)
return(gene2symbol)
}

if (K2riData=="gene2path")
{
gene2path <- get("gene2path",envir=k2ri)
return(gene2path)
}
}
