MeanBySpot <-
function (fileIN, n = 3, name.A = "A", 
    name.M = "M.norm", by.var = "ID", na.rm=TRUE) 
{  

    Temp <- aggregate(fileIN[, c(grep(name.A, names(fileIN)),grep(name.M, names(fileIN)))], list(by.var = fileIN[, match(by.var, names(fileIN))]), 
        mean,na.rm=na.rm)
    fileIN1=fileIN[-which(duplicated(fileIN[,match(by.var,names(fileIN))])),1:n]
    Resultat1 <- merge(fileIN1, Temp, by.x = by.var, by.y = "by.var",all=TRUE)    
    temp <- match(names(fileIN), names(Resultat1))
    Resultat <- Resultat1[, temp[is.na(temp) == FALSE]]
    Resultat
    # (c) 2010 Institut National de la Recherche Agronomique

}

