sva2 <-
function (dat, mod, mod0 = NULL, n.sv = NULL, method = c("irw",  
    "two-step"), vfilter = NULL, B = 5, numSVmethod = "be")  
{ 
    if (is.null(n.sv)) { 
        n.sv = num.sv2(dat, mod, method = numSVmethod, vfilter = vfilter) 
    } 
    if (!is.null(vfilter)) { 
        if (vfilter < 100 | vfilter > dim(dat)[1]) { 
            stop(paste("The number of genes used in the analysis must be between 100 and",  
                dim(dat)[1], "\n")) 
        } 
        tmpv = rowVars(dat) 
        ind = which(rank(-tmpv) < vfilter) 
        dat = dat[ind, ] 
    } 
    if (n.sv > 0) { 
        cat(paste("Number of significant surrogate variables is: ",  
            n.sv, "\n")) 
        method <- match.arg(method) 
        if (method == "two-step") { 
            return(twostepsva.build(dat = dat, mod = mod, n.sv = n.sv)) 
        } 
        if (method == "irw") { 
            return(irwsva.build2(dat = dat, mod = mod, mod0 = mod0, n.sv = n.sv, B = B)) 
        } 
    } 
    else { 
        print("No significant surrogate variables") 
        return(list(sv = 0, pprob.gam = 0, pprob.b = 0, n.sv = 0)) 
    } 
}
