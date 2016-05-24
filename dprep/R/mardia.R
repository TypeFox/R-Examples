mardia <-
function (data) 
{
    p = dim(data)[2]
    f = p - 1
    data[,p]=as.numeric(factor(data[,p]))
    clases = length(table(data[, p]))
    print(clases)
    for (i in 1:clases) {
        data1 = data[data[, p] == i, -p]
        ndat = dim(data1)[1]
        mo3 = mo3(data1)
        mard1 = ndat * mo3/6
        cat("Mardia's test for class", i, "\n")
        cat("mard1=", mard1, "\n")
        p1 = 1 - pchisq(mard1, df = f * (f + 1) * (f + 2)/6)
        cat("pvalue for m3=", p1, "\n")
        mo4 = mo4(data1)
        mard2 = (mo4 - f * (f + 2))/sqrt(8 * f * (f + 2)/ndat)
        cat("mard2=", mard2, "\n")
        p2 = 2 * (1 - pnorm(abs(mard2)))
        cat("p-value for m4=", p2, "\n")
        if (p1 < 0.05 || p2 < 0.05) 
            cat("There is not statistical evidence for normality in class", 
                i, "\n")
        else cat("There is statistical evidence for normality", 
            "\n")
    }
}
