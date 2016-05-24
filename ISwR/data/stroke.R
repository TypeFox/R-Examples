stroke <-  read.csv2("stroke.csv", na.strings=".")
names(stroke) <- tolower(names(stroke))
stroke <-  within(stroke,{
    sex <- factor(sex,levels=0:1,labels=c("Female","Male"))
    dgn <- factor(dgn)
    coma <- factor(coma, levels=0:1, labels=c("No","Yes"))
    minf <- factor(minf, levels=0:1, labels=c("No","Yes"))
    diab <- factor(diab, levels=0:1, labels=c("No","Yes"))
    han <- factor(han, levels=0:1, labels=c("No","Yes"))
    died <- as.Date(died, format="%d.%m.%Y")
    end <- pmin(died, as.Date("1996-01-01"), na.rm=TRUE)
    dstr <- as.Date(dstr,format="%d.%m.%Y")
    obsmonths <- as.numeric(end-dstr, "days")/30.6
    obsmonths[obsmonths==0] <- 0.1
    dead <- !is.na(died) & died < as.Date("1996-01-01")
    died[!dead] <- NA
    rm(end)
})
