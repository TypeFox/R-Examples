clr = function(data, group = NULL){
    geoMean = function(x){
        p = length(x)

        return(prod(x)^(1/p))
    }

    if(!is.null(group)){
        f = formula(paste(group,"~."),sep="")
        mf = model.frame(f, data)
        g = model.response(mf)
        names(g) = group
        data = mf[,-1]
    }

    gms = apply(data, 1, geoMean)

    if(any(gms==0))
        warning("Data has zeros which means the transformation will be unusable")

    data = log(data/gms)

    if(!is.null(group))
        return(data.frame(group, data))

    return(data)
}
