alr = function(form, data, group = NULL){

    if(!is.null(group)){
        f = formula(paste(group,"~."),sep="")
        mf = model.frame(f, data)
        g = model.response(mf)
        names(g) = group
    }
    mf = model.frame(form, data)
    denom = model.response(mf)
    variables = mf[,-1]

    if(sum(variables==0) > 0)
        warning("Data has zeros which will make transformed data unusable")

    variables = log(variables/denom)

    if(!is.null(group))
        return(data.frame(g, variables))

    return(variables)
}
