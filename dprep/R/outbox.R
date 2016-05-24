outbox <-
structure(function (data, nclass) 
{
    ncols = dim(data)[2]
    out1 <- NULL
    datatempo = data[data[, ncols] == nclass, 1:(ncols - 1)]
    blim = boxplot(as.data.frame(datatempo))$stats
    for (i in 1:(ncols - 1)) {
        b1 = as.numeric(rownames(rbind(datatempo[datatempo[, 
            i] < blim[1, i], ], datatempo[datatempo[, i] > blim[5, 
            i], ])))
        out1 = c(out1, b1)
    }
    sort(table(out1))
}, source = c("function(data,nclass)", "{#******************************************************************", 
"#This function detects univariate outliers simultaneously using boxplots", 
"#data: name of the dataset", "#nclass: class number", "#********************************************************************", 
"ncols=dim(data)[2]", "out1<-NULL", "datatempo=data[data[,ncols]==nclass,1:(ncols-1)]", 
"blim=boxplot(as.data.frame(datatempo))$stats", "for(i in 1:(ncols-1)){", 
"b1=as.numeric(rownames(rbind(datatempo[datatempo[,i]<blim[1,i],],datatempo[datatempo[,i]>blim[5,i],])))", 
"out1=c(out1,b1)}", "sort(table(out1))", "}"))
