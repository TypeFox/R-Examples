setClass("haplotype",representation(snp = "integer",qtl = "list",
                                    hID = "numeric",phID0 = "numeric",phID1 = "numeric"),
         prototype = list(hID = 0,phID0 = NaN,phID1 = NaN))

validhaploListObject <- function(object)
  {
    if(class(object)!="haploList")
      FALSE
    else if(all(sapply(object,function(x)class(x)=="haplotype")))
      TRUE
    else
      FALSE
  }

setClass("haploList",contains = "list",representation(genDist = "numeric",nDec = "integer",nChrom = "integer"),
         prototype(genDist = 1,nDec = as.integer(3),nChrom = as.integer(1)),         
         validity = validhaploListObject)





         
                                    
