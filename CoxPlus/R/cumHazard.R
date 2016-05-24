# Doesn't support general case weight, frailty
cumHazard = function(head, formula, par, data=NULL){
    begin = Sys.time()
    if(is.null(par$receiver) || is.null(par$beta)) stop("receiver and beta need to be specified to estimate cumulative hazard!")
    # check on start and stop time
    if(any(head[,'start']>=head[,'stop'])){
        print(head[head[,'start']>=head[,'stop'],])
        stop('Error: start time larger than or equal to stop time!!!')
    }
    # order, should be used on all vectors of length n in par!!!!
    if("strata" %in% colnames(head)) {
        ord = order(head[,'strata'], head[,'event'], -head[,'stop'], -head[,'start'], decreasing=T)
    }else ord = order(head[,'event'], -head[,'stop'], -head[,'start'], decreasing=T)
    head = head[ord,]
    par$receiver=par$receiver[ord]
    
    # offset
    if(!is.null(par$offset)) par$offset = par$offset[ord]
    
    # strata
    if("strata" %in% colnames(head)){
        s = factor(head[,'strata'])
        nStrata = length(unique(s))
        strats = sort(levels(s), decreasing=T) # consistent with ord
        counts = table(s)
        offset = 0;
        strataBounds = matrix(0, nStrata, 2)
        for(i in 1:nStrata){
            strataBounds[i,]=c(offset,offset+counts[strats[i]]-1)
            offset = offset + counts[strats[i]]
        }
        par$strataBounds = strataBounds
    }
    
    
    # model.matrix
    if(is.null(data)) {
        X = model.matrix(formula)
    }else X = model.matrix(formula, data)
    if(nrow(X) != length(ord)) stop("nrow(X)!=length(ord), there might be missing values in the regressors")
    X = X[ord,]
    X = X[, -1, drop=F] # drop intercept
    # deleted not so useful demean because it may affect interaction terms 
    
    # figure out the appropriate names for beta  
    fnames = colnames(X)
    realnames = fnames    
    if(!is.null(par$dropCols)){
        par$dropFlag = as.integer(fnames %in% c(par$dropCols, par$fixedCol))
        realnames = fnames[!fnames %in% c(par$dropCols, par$fixedCol)]
    }
    # time varying cols
    ix = 1:ncol(X)
    names(ix) = fnames    
    if(!is.null(par$fixedCol)) par$fixedCoefIndex = ix[par$fixedCol] - 1
    if(!is.null(par$timeVarCols) && length(intersect(par$timeVarCols,fnames))>0){
        par$timeVarCols = intersect(par$timeVarCols,fnames)
        par$timeVarIndex = ix[par$timeVarCols] -1 
        # print(par$timeVarIndex)
        if("age" %in% par$timeVarCols) par$ageIndex=which(par$timeVarCols=="age")-1
        if("logDiggNum" %in% par$timeVarCols) par$diggNumIndex=which(par$timeVarCols=="logDiggNum")-1
        # construction a map for interaction terms
        subnames = strsplit(fnames, split=":")
        nInteractions = 0
        for(i in 1:ncol(X)){
            snames = unlist(subnames[i])
            if(length(snames)>1 && any(snames %in% par$timeVarCols)){
                nInteractions = nInteractions+1
            }
        }   
        # does not support three way interactions at this moment
        interMap = matrix(rep(-1, 3*nInteractions), ncol=3)
        j = 1
        for(i in 1:ncol(X)){
            snames = unlist(subnames[i])
            if(length(snames)>1 && any(snames %in% par$timeVarCols)){
                interMap[j,] = c(i, ix[snames]) - 1
                j = j+1
            }
        }
        par$interMap = interMap
    }
    
    # verbose
    if(is.null(par$verbose)) par$verbose = 1   
    if(nrow(X)>100000) par$verbose = max(2, par$verbose)   
    # show ties or not distinguishable stop times
    stoptime = head[,'stop'][head[,'event']>0]
    time_diff = abs(diff(stoptime))
    isDistinguishable = time_diff < 1e-9 * pmin(stoptime[1:(length(stoptime)-1)], stoptime[2:length(stoptime)])
    if(!is.null(par$showties) && any(isDistinguishable)) {
        print(paste('Warning: The difference between ', sum(isDistinguishable), ' stop times are not distinguishable!!!'))
        false_ix = which(isDistinguishable)
        row1 = stoptime[false_ix]
        row2 = stoptime[1+false_ix]
        print(cbind(false_ix, row1, row2, abs(row2-row1)/pmin(row1,row2)))
    }
    
    # estimate
    mod=Module("cox_module", dyn.load("E:/Dropbox/FastCox/src/FastCox.dll"))
    X = t(X)   
    head = as.matrix(head[,c("start","stop","event","weight")])
    
    # estimate survival curve
    lvs = levels(factor(par$receiver))
    par$receiver = as.integer(factor(par$receiver)) - 1
    if(!is.null(par$savePath)) save(head,X,par,file=par$savePath)
    inst = new(mod$CoxReg,head,X,par)
    suv = inst$survivalNormal()
    names(suv) = c(strats, "ind","indvar")
    colnames(suv$ind) = strats
    rownames(suv$ind) = lvs
    colnames(suv$indvar) = strats
    rownames(suv$indvar) = lvs
    return (suv)
}