rav <- function(data, subset=NULL, mean=FALSE, lev, s.range=c(NA,NA), w.range=exp(c(-5,5)), I0=FALSE,
    par.fixed=NULL, all=FALSE, IC.diff=c(2,2), Dt=0.1, IC.break=FALSE, t.par=FALSE, verbose=FALSE,
    title=NULL, names=NULL, method="BFGS", start=c(s=NA,w=exp(0)), lower=NULL, upper=NULL, control=list())
{
    # --------------------------------------
    # Inizializzazione variabili e controlli
    # --------------------------------------
    
    # lev deve essere numerico:
    if(!is.numeric(lev)) stop("Mistaken definition of the lev argument: check the values of the vector.\n")
    IC.diff <- abs(IC.diff)
    fact <- length(lev) # numero di fattori
    sumlev <- as.integer(sum(lev)) # numero totale di livelli
    sumlev <- c(sumlev,2+sumlev,2+2*sumlev)
    
    # ----------------------------------
    # Verifiche sul dataset
    # ----------------------------------
    # Se richiesto, si estrae un subset di dati:
    if(!is.null(subset)) {
        set <- FALSE
        for(i in 1:length(subset)) set <- set|(data[,1]==subset[i])
        # Si elimina la colonna che codifica per il subset:
        data <- data[set,-1]
        rm(set)
    }
    dim.data <- dim(data)
    if(is.null(dim.data)) {
        data <- t(as.matrix(data))
        dim.data <- dim(data)
    }
    # Numero corretto di colonne della matrice di dati:
    all.resp <- respdim(lev)
    # Se c'e', si elimina la colonna di codifica, e i codici diventano i nomi di riga:
    if(dim.data[2]!=all.resp) {
        if(dim.data[2] == (all.resp+1)) {
            attr(data,"row.names") <- paste(1:dim.data[1],as.character(data[,1]),sep=".")
            data <- data[,-1]
            dim.data[2] <- dim.data[2]-1
        } else stop("Error occurred in columns check")
    }
    # A questo punto, se data e' un data frame lo si trasforma in matrice:
    data <- as.matrix(data)
    if(ncol(data)==1) {
        data <- t(data)
        dim.data <- dim(data)
    }
    # Se l'utente vuole lavorare sulle medie e non sui dati grezzi:
    if(mean==TRUE) {
        data <- t(as.matrix(colMeans(data)))
        dim.data <- dim(data)
    }
    
    N <- sum(!is.na(data)) # numero di R osservati
    TSS <- sum((data-mean(data,na.rm=TRUE))^2,na.rm=TRUE) # Total Sum of Squares
    if(TSS==0) stop("No variability found in data (TSS = 0)")
    
    # Attribuzione dei valori di default a names:
    if(is.null(names)) names <- paste("Factor",LETTERS[1:fact])
    else names <- as.character(names[1:fact])
    
    # Indici per l'estrazione dei parametri s e w dai vettori:
    pos <- list(
        slast = seq_len(sumlev[2]), # sequenza che va da 1 alla posizione dell'ultimo s
        spos = seq_len(sumlev[1])+2, # posizioni degli s (eccetto s0)    
        spos.fact = list(first=3+cumsum(c(0,lev[-fact])), last=2+cumsum(lev)), # posizioni del primo e dell'ultimo s di ogni fattore
        wfirst = as.integer(cumsum(c(1,lev[-fact]))+sumlev[2]), # posizioni dei primi pesi di ogni fattore
        wpos = (sumlev[2]+1):sumlev[3], # posizioni di tutti i pesi eccetto w0
        w0wpos = c(2,(sumlev[2]+1):sumlev[3]), # posizioni di tutti i pesi incluso w0
        fixed = NULL # posizioni dei valori fissati
    )
    pos$s0spos <- c(1,pos$spos) # posizioni degli s considerando s0
    pos$swfirst <- c(pos$slast,pos$wfirst) # Posizioni di: s0, w0, tutti gli s, primo w di ogni fattore
    
    # Organizzazione del vettore di parametri fissi
    fixed <- rep.int(NaN,sumlev[3]) # Vettore dei parametri fissi
    if(is.logical(I0)) {
        if(!I0) {
            fixed <- c(0,1e-10)
            # Eventuale correzione di fiixed[2] e w.range sulla base di lower e upper:
            if(!is.null(lower)) {
                fixed[2] <- lower[2]
                w.range[1] <- min(lower[pos$wpos])
            }
            if(!is.null(upper))
                w.range[2] <- max(lower[pos$wpos])
        } else {
            # Eventuale correzione di w.range sulla base di lower e upper:
            if(!is.null(lower))
                w.range[1] <- min(lower[pos$w0wpos])
            if(!is.null(upper))
                w.range[2] <- max(lower[pos$w0wpos])
        }
    } else {
        fixed[1:2] <- c(0,I0)
        I0 <- FALSE
    }
    
    # -----------------------------------------------
    # Fissaggio dei valori di scala e dei pesi
    # -----------------------------------------------
    s.fixed <- w.fixed <- FALSE
    if(!is.null(par.fixed)) {
        if(is.list(par.fixed)) {
            val.fixed <- par.fixed # Valori dei parametri da fissare
            par.fixed <- names(par.fixed) # Nomi dei parametri da fissare
            # Gli s sono da fissare?
            if(any(par.fixed=="s"))
                s.fixed <- val.fixed$s
            # Quale tipo di parametro tra w e t e' da fissare?
            type.par <- c(any(par.fixed=="w"),any(par.fixed=="t"))
            if(sum(type.par)==2) # entrambi
                stop("Incorrect par.fixed definition: you must choose between 'w' and 't'\n")
            if(sum(type.par)==1) { # uno e' da fissare
                type.par <- c("w","t")[type.par]
                if(type.par=="w")
                    w.fixed <- val.fixed$w
                else
                    w.fixed <- val.fixed$t
            }
        } else {
            if(is.vector(par.fixed)) {
                # (vettore di caratteri che indicano i nomi dei parametri)
                # Gli s sono da fissare?
                s.fixed <- any(par.fixed=="s")
                if(isTRUE(s.fixed))
                    s.fixed <- fixparam(lev,"s",names) # apre l'interfaccia per fissare gli s
                # Quale tipo di parametro tra w e t e' da fissare?
                type.par <- c(any(par.fixed=="w"),any(par.fixed=="t"))
                if(sum(type.par)==2) # entrambi
                    stop("Incorrect par.fixed definition: you must choose between 'w' and 't'\n")
                if(sum(type.par)==1) { # uno e' da fissare
                    type.par <- c("w","t")[type.par]
                    w.fixed <- fixparam(lev,type.par,names) # apre l'interfaccia per fissare i w
                }
            }
        }
        if(is.list(s.fixed)) {
            fixed[pos$spos] <- unlist(s.fixed)
            s.fixed <- TRUE
        }
        if(is.list(w.fixed)) {
            fixed[pos$wpos] <- unlist(w.fixed)
            if(type.par=="w")
                fixed[pos$wpos][which(fixed[pos$wpos]==0)] <- w.range[1]
            if(type.par=="t")
                fixed[pos$wpos] <- exp(fixed[pos$wpos]) # I t vanno trasformati in w
        }
    }
    pos$fixed <- which(!is.na(fixed))
    pos$spos.fix <- which(!is.na(fixed[pos$spos]))+2
    
    # Specifiche dei pesi fissati:
    nwfix <- list(
        num = NULL, # quanti pesi fissi
        pos = which(!is.na(fixed[pos$wpos])), # posizioni
        nwval = fixed[pos$wpos][which(!is.na(fixed[pos$wpos]))] # valori
    )
    nwfix$num <- length(nwfix$nwval)
    
    # Specifiche degli s fissati:
    nsfix <- list(
        num = NULL,
        pos = which(!is.na(fixed[pos$spos])), # posizioni
        nsval = fixed[pos$spos][which(!is.na(fixed[pos$spos]))] # valori
    )
    nsfix$num <- length(nsfix$nsval)
    
    # Trasformazione dei w in t:
    fixed[pos$w0wpos] <- log(fixed[pos$w0wpos])
    
    # Numero di parametri fissi da sottrarre a n.pars:
    if(nwfix$num==0)
        nwfix$n.pars <- 0
    else
        nwfix$n.pars <- numpar(fixed[pos$wpos][nwfix$pos],nwfix$num)
    
    # Non ci devono essere pesi fissi uguali entro il delta:
    if(nwfix$num>1) {
        w.fixed.diff <- NULL
        for(i in 1:nwfix$num)
            w.fixed.diff <- c(w.fixed.diff,nwfix$nwval[i]-nwfix$nwval[-c(1:i)])
        w.fixed.diff <- abs(w.fixed.diff)
        if(sum(w.fixed.diff > 1e-7 & w.fixed.diff<=Dt))
            stop("Some of the fixed weights are different within the Dt")
    }
    
    par.base <- 2*I0+sumlev[1]-nsfix$num
    paramEmpty <- c(0,-6, rep.int(NA,sumlev[1]), rep.int(0,sumlev[1]))
    
    # ----------------------------------
    # Start, lower, upper
    # ----------------------------------
    if(is.na(start[1]))
        start[1] <- mean(data,na.rm=TRUE)
    # Start -------
    startNull <- start[1] # start s per modello null00
    start <- c(NA,NA,start)
    start[4] <- log(start[4])
    if(I0) {
        start[1] <- start[3]
        start[2] <- start[4]
    }
    start[pos$fixed] <- fixed[pos$fixed]
    
    if(method=="L-BFGS-B") {
        if(is.na(s.range[1])) {
            s.range[1] <- min(data,na.rm=TRUE)
            s.range[1] <- floor(s.range[1])
        }
        if(is.na(s.range[2])) {
            s.range[2] <- max(data,na.rm=TRUE)
            s.range[2] <- ceiling(s.range[2])
        }
        w.range <- log(w.range)
        # Upper -------
        if(is.null(upper)) {
            upper <- list(I0=c(mean(s.range),log(2)),s=NULL,t=NULL)
            for(i in 1:fact) {
                upper$s <- c(upper$s,rep.int(s.range[2],lev[i]))
                upper$t <- c(upper$t,rep.int(w.range[2],lev[i]))
            }
            upper <- unlist(upper)
        }
        upperNull <- list(
            null00 = s.range[2],
            null01 = rep.int(s.range[2],fact),
            null02 = rep.int(s.range[2],sumlev[1])
        )
        # Lower -------
        if(is.null(lower)) {
            lower <- list(I0=fixed[1:2],s=NULL,t=NULL)
            for(i in 1:fact) {
                lower$s <- c(lower$s,rep.int(s.range[1],lev[i]))
                lower$t <- c(lower$t,rep.int(w.range[1],lev[i]))
            }
            if(I0) lower$I0 <- start[1:2]
            lower <- unlist(lower)
        }
        lowerNull <- list(
            null00 = s.range[1],
            null01 = rep.int(s.range[1],fact),
            null02 = rep.int(s.range[1],sumlev[1])
        )
        if(I0) {
            upperNull$null01 <- c(upper[1],upperNull$null01)
            upperNull$null02 <- c(upper[1],upperNull$null02)
            lowerNull$null00 <- c(lower[1],lowerNull$null00)
            lowerNull$null01 <- c(lower[1],lowerNull$null01)
            lowerNull$null02 <- c(lower[1],lowerNull$null02)
        }
        # EAM: di tutti i pesi passati come start bisogna mantenerne solo uno per fattore:
        upper.eq <- upper[pos$swfirst]
        lower.eq <- lower[pos$swfirst]
    } else {
        upper <- (+Inf)
        lower <- (-Inf)
        upperNull <- list(null00=Inf,null01=Inf,null02=Inf)
        lowerNull <- list(null00=(-Inf),null01=(-Inf),null02=(-Inf))
        upper.eq <- (+Inf)
        lower.eq <- (-Inf)
    }
    rm(all.resp,s.range,w.range)
    
    # --------------------------------
    # NULL MODEL
    # --------------------------------
    output <- list(
        par = startNull,
        value = TSS,
        convergence = 0,
        message = "",
        n.pars = 1,
        change = FALSE
    )
    paramNull <- paramEmpty
    paramNull[pos$spos] <- output$par
    output$par <- paramNull
    output$n.pars <- 1
    
    # Istanza di classe:
    case.null <- fit.indices(output=output, data=data, lev=lev, fact=fact, sumlev=sumlev,
        N=N, I0=I0, dim.data=dim.data, TSS=TSS, names=names, model=1, start=startNull,
        upper=upperNull$null00, lower=lowerNull$null00, t.par=t.par)
    
    # --------------------------------
    # EQUAL SCALE VALUES MODEL
    # --------------------------------
    
    startNull <- case.null@param[pos$spos.fact$first]
    if(I0) {
        startNull <- c(case.null@param[1],startNull)
        paramEmpty[2] <- 0
    }
    # Stima dei parametri:
    output <- optim(par=startNull, fn=Residual.null, fixed=fixed, model=2, I0=I0, data=data,
        lev=lev, fact=fact, sumlev=sumlev, sposFirst=pos$spos.fact$first-1, dim.data=dim.data,
        method=method, lower=lowerNull$null01, upper=upperNull$null01, control=control)
    
    if(!is.null(output$message)) {
        if(output$message == "ERROR: NO FEASIBLE SOLUTION") {
            # Se optim fornisce questo messaggio di errore, l'RSS in output non corrisponde col valore
            # della funzione obiettivo (forse si tratta di un bug). Bisogna quindi correggere l'output:
            output$value <- Residual.null(output$par, fixed=fixed[pos$slast], model=2, I0, data, lev,
                fact, sumlev, sposFirst=pos$spos.fact$first-1, dim.data)
        }
    }
    
    paramNull <- paramEmpty
    if(I0) {
        paramNull[1] <- output$par[1]
        output$par <- output$par[-1]
    }
    for(i in seq_len(fact))
        paramNull[pos$spos.fact$first[i]:pos$spos.fact$last[i]] <- output$par[i]
    output$par <- paramNull
    
    # Numero di parametri del modello:
    output$n.pars <- fact+2*I0
    
    # Istanza di classe:
    case.ESM <- fit.indices(output=output, data=data, lev=lev, fact=fact, sumlev=sumlev,
        N=N, I0=I0, dim.data=dim.data, TSS=TSS, names=names, model=2, start=startNull,
        upper=upperNull$null01, lower=lowerNull$null01, t.par=t.par)
    
    # Confronto col modello precedente:
    selectedModel <- 1
    case.baseline <- case.null
    if((case.ESM@BIC+IC.diff[1]) < case.baseline@BIC) {
        case.baseline <- case.ESM
        selectedModel <- 2
    } else {
        if((case.ESM@AIC+IC.diff[2]) < case.baseline@AIC) {
            case.baseline <- case.ESM
            selectedModel <- 2
        }
    }
    
    # --------------------------------
    # SIMPLE AVERAGING MODEL
    # --------------------------------
    
    if(!I0) startNull <- NULL
    else startNull <- case.baseline@param[1]
    for(i in 1:fact)
        startNull <- c(startNull,rep.int(case.baseline@param[pos$spos.fact$first][i],lev[i]))
    if(s.fixed) {
        if(!I0)
            startNull[pos$spos.fix-2] <- fixed[pos$spos.fix]
        else
            startNull[-1][pos$spos.fix-2] <- fixed[pos$spos.fix]
    }
    
    # Stima dei parametri:
    output <- optim(par=startNull, fn=Residual.null, fixed=fixed[pos$slast], model=3,
        I0=I0, data=data, lev=lev, fact=fact, sumlev=sumlev, sposFirst=pos$spos.fact$first-1,
        dim.data=dim.data, method=method, lower=lowerNull$null02, upper=upperNull$null02, control=control)
    
    if(!is.null(output$message)) {
        if(output$message == "ERROR: NO FEASIBLE SOLUTION") {
            # Se optim fornisce questo messaggio di errore, l'RSS in output non corrisponde al valore
            # della funzione obiettivo (forse si tratta di un bug). Bisogna quindi correggere l'output:
            output$value <- Residual.null(output$par, fixed=fixed[pos$slast], model=3, I0, data, lev,
                fact, sumlev, sposFirst=pos$spos.fact$first-1, dim.data)
        }
    }
    
    paramNull <- paramEmpty
    if(I0) {
        paramNull[1] <- output$par[1]
        output$par <- output$par[-1]
    }
    paramNull[pos$spos] <- output$par
    output$par <- paramNull
    rm(paramNull,paramEmpty)
    
    # Numero di parametri del modello:
    output$n.pars <- sumlev[1]-nsfix$num+2*I0
    
    # Istanza di classe:
    case.SAM <- fit.indices(output=output, data=data, lev=lev, fact=fact, sumlev=sumlev,
        N=N, I0=I0, dim.data=dim.data, TSS=TSS, names=names, model=3, start=startNull,
        upper=upperNull$null02, lower=lowerNull$null02, t.par=t.par)
    
    # Confronto col modello precedente:
    if((case.SAM@BIC+IC.diff[1]) < case.baseline@BIC) {
        case.baseline <- case.SAM
        selectedModel <- 3
    } else {
        if((case.SAM@AIC+IC.diff[2]) < case.baseline@AIC) {
            case.baseline <- case.SAM
            selectedModel <- 3
        }
    }
    
    # --------------------------------
    # STARTING VALUES PER EQUAL E DIFF
    # --------------------------------
    # Equal weight case:
    start.eq <- case.baseline@param[pos$swfirst]
    if(s.fixed) start.eq[pos$spos.fix] <- fixed[pos$spos.fix]
    start.eq[2]  <- start[2]
    start.eq[(sumlev[2]+1):(sumlev[2]+fact)] <- start[4]
    # Differential weight case:
    start.diff <- case.baseline@param
    start.diff[2] <- start[2]
    start.diff[pos$wpos] <- start[4]
    start.diff[pos$fixed] <- fixed[pos$fixed]
    
    # --------------------------------
    # EQUAL-WEIGHTS MODEL
    # --------------------------------
    
    # Stima dei parametri:
    output <- optim(par=start.eq, fn=Residual.eq, fixed=fixed, data=data, lev=lev, fact=fact,
        sumlev=sumlev, dim.data=dim.data, method=method, lower=lower.eq, upper=upper.eq, control=control)
    
    # Replica dei t per i livelli dei fattori:
    param.equal <- output$par[pos$slast]
    for(i in 1:fact)
        param.equal <- c(param.equal,rep.int(output$par[pos$wpos[i]],lev[i]))
    output$par <- param.equal
    
    # Si eguagliano i pesi uguali entro delta:
    output$par <- parmeanlast(output$par,fixed,sumlev,Dt,nwfix=list(num=0,nwval=0))
    
    # Si scalano i t:
    if(!I0)
        output$par[pos$wpos] <- output$par[pos$wpos]-mean(output$par[pos$wpos])
    else
        output$par[pos$w0wpos] <- output$par[pos$w0wpos]-mean(output$par[pos$w0wpos])
    
    # Dato che i parametri potrebbero essere stati modificati si ricalcola l'RSS:
    output$value <- Residual(output$par,fixed,data,lev,fact,sumlev,dim.data,Dt,nwfix)
    
    # Calcolo del numero di parametri del modello:
    output$n.pars <- par.base+numpar(output$par[pos$wfirst],fact)-1
    
    # Istanza di classe:
    case.EAM <- fit.indices(output=output, data=data, lev=lev, fact=fact, sumlev=sumlev,
        N=N, I0=I0, dim.data=dim.data, TSS=TSS, names=names, model=4, start=start.eq,
        upper=upper.eq, lower=lower.eq, t.par=t.par)
    
    # Confronto col modello precedente:
    if((case.EAM@BIC+IC.diff[1]) < case.baseline@BIC) {
        case.baseline <- case.EAM
        selectedModel <- 4
    } else {
        if((case.EAM@AIC+IC.diff[2]) < case.baseline@AIC) {
            case.baseline <- case.EAM
            selectedModel <- 4
        }
    }
    
    # --------------------------------
    # DIFFERENTIAL-WEIGHTS MODEL
    # --------------------------------
    
    # Stima dei parametri:
    output <- optim(par=start.diff, fn=Residual, fixed=fixed, data=data, lev=lev, fact=fact,
        sumlev=sumlev, dim.data=dim.data, Dt=Dt, nwfix=nwfix[c(1,3)],
        method=method, lower=lower, upper=upper, control=control)
    
    # Si eguagliano i t uguali entro il delta:
    # output$par <- parmeanlast(output$par,fixed,sumlev,Dt,nwfix)
    
    # Si scalano i t:
    if(nwfix$num==0) {
        if(!I0)
            output$par[pos$wpos] <- output$par[pos$wpos]-mean(output$par[pos$wpos])
        else
            output$par[pos$w0wpos] <- output$par[pos$w0wpos]-mean(output$par[pos$w0wpos])
    }
    
    # Calcolo del numero di parametri del modello:
    if(nwfix$num==0)
        output$n.pars <- par.base+numpar(output$par[pos$wpos],sumlev[1])-1
    else
        output$n.pars <- par.base+numpar(output$par[pos$wpos[-nwfix$pos]],sumlev[1]-nwfix$num)-1
    
    # Dato che i parametri potrebbero essere stati modificati si ricalcola l'RSS:
    output$value <- Residual(output$par,fixed,data,lev,fact,sumlev,dim.data,Dt,nwfix)
    
    # Istanza di classe:
    case.DAM <- fit.indices(output=output, data=data, lev=lev, fact=fact, sumlev=sumlev,
        N=N, I0=I0, dim.data=dim.data, TSS=TSS, names=names, model=5, start=start.diff,
        upper=upper, lower=lower, t.par=t.par)
    
    # --------------------------------
    # INFORMATION CRITERIA
    # --------------------------------
    
    case.start <- case.baseline
    if(s.fixed | nwfix$num) {
        # Si sostituiscono i parametri fissi a quelli del modello di start:
        case.start@param[pos$fixed] <- fixed[pos$fixed]
        # Calcolo del nuovo RSS:
        case.start@RSS <- Residual(case.start@param,fixed,data,lev,fact,sumlev,dim.data,Dt,nwfix)
        # Nuova istanza di classe:
        output$par <- case.start@param
        output$value <- case.start@RSS
        output$convergence <- 0
        output$message <- ""
        output$n.pars <- case.start@n.pars
        case.start <- fit.indices(output=output, data=data, lev=lev, fact=fact, sumlev=sumlev,
            I0=I0, N=N, dim.data=dim.data, TSS=TSS, names=names, model=4, t.par=t.par)
    }
    # Stima dei parametri:
    if(IC.break==FALSE) {
        output <- optimization.IC(data=data, fact=fact, lev=lev, sumlev=sumlev, pos=pos, N=N,
            dim.data=dim.data, model.start=case.start, par.base=par.base, nwfix=nwfix,
            fixed=fixed, I0=I0, Dt=Dt, IC.diff=IC.diff, all=all, verbose=verbose, IC.break=IC.break,
            lower=lower, upper=upper, method=method, control=control, change=FALSE
        )
    } else {
        output <- list(
            par = case.start@param,
            value = case.start@RSS,
            convergence = 2,
            message = "MNA", # Model Not Available
            n.pars = case.start@n.pars,
            change = FALSE
        )
    }
    
    # Istanza di classe:
    case.IC <- fit.indices(output=output, data=data, lev=lev, fact=fact, sumlev=sumlev,
        N=N, I0=I0, dim.data=dim.data, TSS=TSS, names=names, model=6, start=case.start@param,
        upper=upper, lower=lower, t.par=t.par)
    
    # --------------------------------
    # MODEL SELECTION
    # --------------------------------
    selectedModel <- c(selectedModel,5,6)
    
    # Selezione in base al BIC
    Bic <- c(case.baseline@BIC,case.DAM@BIC,case.IC@BIC)  
    if(IC.break | !output$change) {
        Bic <- Bic[-3]
        selectedModel <- selectedModel[-3]
    }
    selection <- (Bic-min(Bic))<IC.diff[1]
    # Se e' stato selezionato piu' di un modello:
    if(sum(selection)>1) {
        # Selezione in base all'AIC
        Aic <- c(case.baseline@AIC,case.DAM@AIC,case.IC@AIC)
        if(IC.break | !output$change)
            Aic <- Aic[-3]
        selection <- (Aic-min(Aic))<IC.diff[2]
    }
    # selectedModel conteneva il modello scelto dalla EAM selection,
    # ora verra' aggiornato con questa ultima selezione:
    selectedModel <- selectedModel[selection]
    
    if(!t.par) {
        case.null@param[pos$w0wpos] <- exp(case.null@param[pos$w0wpos])
        case.ESM@param[pos$w0wpos] <- exp(case.ESM@param[pos$w0wpos])
        case.SAM@param[pos$w0wpos] <- exp(case.SAM@param[pos$w0wpos])
        case.EAM@param[pos$w0wpos] <- exp(case.EAM@param[pos$w0wpos])
        case.DAM@param[pos$w0wpos] <- exp(case.DAM@param[pos$w0wpos])
        case.IC@param[pos$w0wpos] <- exp(case.IC@param[pos$w0wpos])
    }
    if(!I0) {
        case.null@param[1:2] <- NA
        case.ESM@param[1:2] <- NA
        case.SAM@param[1:2] <- NA
        case.EAM@param[1:2] <- NA
        case.DAM@param[1:2] <- NA
        case.IC@param[1:2] <- NA
    }
    
    new("rav", observed=data, factors=fact, levels=lev, title=as.character(title),
        names=names, method=method, start=start, lower=lower, upper=upper, selection=selectedModel,
        null=case.null, ESM=case.ESM, SAM=case.SAM, EAM=case.EAM, DAM=case.DAM, IC=case.IC)
}
