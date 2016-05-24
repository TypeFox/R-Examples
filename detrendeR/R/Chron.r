#บบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบ
Chron = function (x, prefix = "CRONO", biweight = TRUE, prewhiten = FALSE,
                  stc=c(5,2,1), order.max.prewhiten = 3)
{
    prefix = as.character(prefix)

    if (order.max.prewhiten < 1) prewhiten = FALSE
    #if (nchar(prefix) > 5)  prefix=left(prefix, 5)                              
       
    SAMPLE_DEPTH = sampleDepth(x, stc=stc);                                     
    
    {if (!biweight)
          std = apply(x, 1, mean, na.rm = TRUE)
     else
          std = apply(x, 1, tbrm, C = 9)
    }

     out = data.frame(std, SAMPLE_DEPTH)                                        #FC
     colnames(out) = c(prefix,  "NoTrees", "NoCores")                           #FC

    if (prewhiten) {
       PREFIX = paste(prefix, "prewhiten",sep="-")

       
        x.ar = apply(x, 2, ar.func, order.max=order.max.prewhiten)
        {if (!biweight)
            res = apply(x.ar, 1, mean, na.rm = TRUE)
        else
        res = apply(x.ar, 1, tbrm, C = 9)}
        res[is.nan(res)] = NA
        out = data.frame(std, res,   SAMPLE_DEPTH)
        colnames(out) = c(prefix, PREFIX,  "NoTrees", "NoCores")           #FC
    }
    return(out)
}
#บบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบบ
# Chron(rwl, prewhiten =T, order.max.prewhiten = 3)