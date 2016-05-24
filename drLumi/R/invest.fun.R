invest.fun <- function(model,bkg.method, in.fun, yy, parameters, log10.bkgmean){
    if(bkg.method!="constraint"){
        if(in.fun=="SSl4"){
            if(length(parameters)!=4){
                stop("Number of parameters and 'in.fun' model does not match")
            }
            inv <- lapply(yy, 
                    function(x) do.call(inSSl4,c(x,as.list(coef(model)))))
        } 
        if(in.fun=="SSl5"){
            if(length(parameters)!=5){
                stop("Number of parameters and 'in.fun' model does not match")
            }
            inv <- lapply(yy, 
                        function(x) do.call(inSSl5,c(x,as.list(coef(model)))))
        } 
        if(in.fun=="SSexp"){
            if(length(parameters)!=2){
                stop("Number of parameters and 'in.fun' model does not match")
            }
            inv <- lapply(yy, 
                        function(x) do.call(inSSexp,c(x,as.list(coef(model)))))
        }      
        form <- lapply(yy, 
                function(x) as.formula(gsub("fx", x, inv[[1]]$formtext)))
    } else {
        if(in.fun=="SSl4"){
            if(length(parameters)!=3){
                stop("Number of parameters and 'in.fun' model does not match")
            }
            inv <- lapply(yy, function(x) {
            do.call(inSSl4,list(x,coef(model)[1], 
                    log10.bkgmean,coef(model)[2],coef(model)[3]))
                })
        }      
        if(in.fun=="SSl5"){
            if(length(parameters)!=4){
                stop("Number of parameters and 'in.fun' model does not match")
            }
            inv <- lapply(yy, function(x){
                do.call(inSSl5,list(x,coef(model)[1], 
                        log10.bkgmean,coef(model)[2],
                        coef(model)[3],coef(model)[4]))
            })
        } 
        if(in.fun=="SSexp"){
            if(length(parameters)!=1){
                stop("Number of parameters and 'in.fun' model does not match")
            }
            inv <- lapply(yy, function(x){
                    do.call(inSSexp,list(x,coef(model)[1], log10.bkgmean))
                    })        
        }        
        form <- lapply(yy, function(x){
            as.formula(gsub("fx", x, inv[[1]]$formtext.cons))
            })
        form <- gsub("cons", log10.bkgmean, form)
        form <- lapply(form, function(x) as.formula(x))
    }
    ans <- list(form=form, inv=inv)
    return(ans)
}

