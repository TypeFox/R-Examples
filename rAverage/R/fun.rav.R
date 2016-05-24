# --------------------------------------
# Internal functions
# --------------------------------------

# Calcola il numero di risposte osservate per un dato disegno fattoriale:
respdim <- function(lev)
{
    Rlen <-
        .C("respdim",
            lev = as.integer(lev),
            fact = as.integer(length(lev)),
            dim = as.integer(0),
            PACKAGE = "rAverage"
        )$dim
    return(Rlen)
}

# Ricompone gli R a partire dai parametri in input. E' una
# funzione interna che non puo' essere richiamata dall'utente.
averaging <- function(param,lev,fact,sumlev)
{
    # La funzione riceve i t, che saranno trasformati in w nella averaging
    Rlen <- respdim(lev)
    R <-
        .C("averaging",
            param = as.double(param),
            lev = as.integer(lev),
            fact = as.integer(fact),
            sumlev = as.integer(sumlev),
            R = numeric(Rlen),
            PACKAGE = "rAverage"
        )$R
    return(R)
}

# Crea una matrice contenente tutte le possibili k combinazioni
# di n elementi. Corrisponde alla funzione combn, nativa di R,
# ma e' molto piu' veloce.
combin <- function(n,k)
{
    ncomb <- choose(n,k)
    x <-
    .C("combin",
        n = as.integer(n),
        k = as.integer(k),
        out = as.integer(numeric(k*ncomb)),
        PACKAGE = "rAverage"
    )$out
    return(matrix(x,nrow=ncomb,byrow=TRUE))
}

# Crea tutte le combinazioni dei livelli dei fattori
# specificati, creando la griglia di un disegno fattoriale.
# Fa un lavoro simile alla funzione expand.grid, nativa di
# R, ma e' molto piu' veloce. N.B.: attualmente questa
# funzione e' utilizzata solo da 'datgen'.
design.grid <- function(lev)
{
    fact <- length(lev)
    x <-
    .C("grid",
        lev = as.integer(lev),
        fact = as.integer(fact),
        out = as.integer(numeric(fact*prod(lev))),
        PACKAGE = "rAverage"
    )
    x <- matrix(x$out,ncol=fact)
    return(x)
}

# Crea un vettore di nomi per le colonne di una matrice di dati:
data.names <- function(lev, names=NULL)
{
    fact <- length(lev)
    if(is.null(names))
        names <- LETTERS[1:fact]
    R.names <- NULL
    # Disegni a 1 via
    for(i in 1:fact)
        R.names <- c(R.names,paste(names[i],1:lev[i],sep=""))
    # Disegni a k vie
    for(k in 2:fact) {
        cc <- t(combin(fact,k))
        for(j in 1:ncol(cc)) {
            f <- names[cc[,j]]
            g <- design.grid(lev[cc[,j]])
            for(i in 1:nrow(g)) {
                each <- paste(f,g[i,],sep="")
                collapse <- paste(each[1],each[2],sep="")
                h <- 3
                while(h<=k) {
                    collapse <- paste(collapse,each[h],sep="")
                    h <- h+1
                }
                R.names <- c(R.names,collapse)
            }
        }
    }
    return(R.names)
}

# Calcola il numero di parametri (utilizzata per contare i pesi).
numpar <- function(param,len)
{
    .C("numpar",
        param = as.double(param),
        len = as.integer(len),
        out = as.integer(len),
        PACKAGE = "rAverage"
    )$out
}

# Individua i pesi uguali entro il delta e li eguaglia alla loro media.
parmeanlast <- function(param, fixed, sumlev, Dt, nwfix)
{
    .C("parmeanlast",
        param = as.double(param),
        fixed = as.double(fixed),
        sumlev = as.integer(sumlev),
        deltaweights = as.double(Dt),
        numfix = as.integer(nwfix$num),
        valfix = as.double(nwfix$nwval),
        NAOK = TRUE,
        PACKAGE = "rAverage"
    )$param
}

# --------------------------------------
# External functions
# --------------------------------------

# Costruisce, per un dato disegno sperimentale, una matrice
# vuota (NA) da riempire con i dati da analizzare con rav.
# da analizzare con rav.
rav.grid <- function(lev,trials=1,subset=FALSE,names=NULL)
{
    all.resp <- respdim(lev)
    mat <- data.frame(matrix(NA,nrow=trials,ncol=all.resp+subset))
    col.names <- data.names(lev,names)
    if(subset) col.names <- c("subset",col.names)
    colnames(mat) <- col.names
    return(mat)
}

# Crea una istanza della classe 'indices'.
fit.indices <- function(output,data,lev,fact,sumlev,N,I0,dim.data,TSS,names,model,
    start=c(0,0),upper=c(0,0),lower=c(0,0),t.par)
{
    title <- modelNames[model]
    Bic <- N*log(output$value/N)+output$n.pars*log(N)
    Aic <- N*log(output$value/N)+2*output$n.pars
    if(N/output$n.pars < 40) {
        Bic <- Bic+(output$n.pars*log(N)*(output$n.pars+1))/(N-output$n.pars-1)
        Aic <- Aic+(2*output$n.pars*(output$n.pars+1))/(N-output$n.pars-1)
    }
    fitted <- averaging(output$par,lev,fact,sumlev)
    names(fitted) <- colnames(data)
    estim <- matrix(rep.int(fitted,dim.data[1]),nrow=dim.data[1],byrow=TRUE)
    # Indici di adattamento:
    R2 <- 1-output$value/TSS
    adjR2 <- (1-(1-R2)*((N-1)/(N-output$n.pars-1)))
    if(length(output$message)==0)
        output$message <- ""
    if(output$message!="MNA") {
        msg <- ""
        if(output$convergence!=0) {
            if(output$convergence==1)
                msg <- "Convergence error: the iteration limit maxit had been reached"
            if(output$convergence==2)
                msg <- output$message
            if(output$convergence==10)
                msg <- "Convergence error: degeneracy of the Nelder-Mead simplex"
            if(output$convergence==51)
                msg <- paste("Convergence warning:",output$message)
            if(output$convergence==52)
                msg <- paste("Convergence error:",output$message)
        }
        if(output$value>TSS) {
            if(msg=="") msg <- "FATAL ERROR: the model does not fit the data"
            else msg <- paste("FATAL ERROR: the model does not fit the data\n",msg,sep="")
        }
    } else { msg <- output$message }
    new("indices",
        param=output$par, I0=I0, start=start, upper=upper, lower=lower,
        levels=lev, fitted=fitted, residuals=data-estim, AIC=Aic, BIC=Bic,
        RSS=output$value, TSS=TSS, R2=R2, adjR2=adjR2, n.pars=output$n.pars,
        t.par=t.par, title=c(title,"\n"), names=names, message=msg, model=model)
}

# Esegue lo stesso lavoro di fit.indices, pero' non viene invocata da rav()
# ma da una chiamata diretta dell'utente.
rav.indices <- function(param,lev,data,t.par=FALSE,subset=NULL,n.pars=NULL,names=NULL,title=NULL)
{
    # Verifiche su names e title
    if(is.null(names)) names <- paste("Factor",LETTERS[1:length(lev)])
    if(is.null(title)) title <- as.character()
    fact <- length(lev)
    sumlev <- as.integer(sum(lev))
    sumlev <- c(sumlev,2+sumlev,2+2*sumlev)
    
    # subset e data
    if(!is.null(subset)) {
        set <- FALSE
        for(i in 1:length(subset)) set <- set|(data[,1]==subset[i])
        # Si elimina la colonna che codifica per i soggetti:
        data <- data[set,-1]
        rm(set)
    }
    data <- as.matrix(data)
    dim.data <- dim(data) 
    if(dim.data[2]==1) data <- t(data)
    all.resp <- respdim(lev)
    # Se c'e' la colonna di raggruppamento la si elimina:
    if(dim.data[2]!=all.resp) {
        if(dim.data[2] == (all.resp+1)) {
            data <- data[,-1]
            dim.data[2] <- dim.data[2]-1
        } else stop("Error occurred in columns check\n")
    }
    # Costruzione finto modello di output
    sumlev <- as.integer(sum(lev)) # numero totale di livelli
    sumlev <- c(sumlev,2+sumlev,2+2*sumlev)
    wpos <- (sumlev[2]+1):sumlev[3]
    if(is.na(param[1])) {
        I0 <- FALSE
        param[1] <- 0
        param[2] <- log(1e-10)
    } else
        I0 <- TRUE
    estim <- averaging(param,lev,fact,sumlev)
    estim <- matrix(rep.int(estim,dim.data[1]),nrow=dim.data[1],byrow=TRUE)
    output <- list(par=param,value=NULL,n.pars=NULL,convergence=0,message="")
    output$value <- sum((data-estim)^2,na.rm=TRUE)
    if(is.null(n.pars))
        output$n.pars <- 2*I0+sumlev[1]+numpar(param[wpos],sumlev[1])
    else
        output$n.pars <- n.pars
    N <- sum(!is.na(data)) # numero di R osservati
    TSS <- sum(data^2,na.rm=TRUE)-(sum(data,na.rm=TRUE)^2)/N # Total Sum of Squares
    if(!t.par) {
        if(I0) wpos <- c(2,wpos)
        output$par[wpos] <- log(output$par[wpos])
    }
    model <- fit.indices(output,data,lev,fact,sumlev,N,I0,dim.data,TSS,names,model=5,
        start=c(s=0,w=0),upper=c(s=0,w=0),lower=c(s=0,w=0),t.par=t.par)
    if(!t.par)
        model@param[wpos] <- exp(model@param[wpos])
    model@title <- c(title,"\n")
    return(model)
}

# Genera parametri random per il modello averaging.
pargen <- function(lev,s.range=c(0,20),w.range=exp(c(-5,5)),I0=FALSE,t.par=FALSE,digits=2)
{
    w.range <- log(w.range)
    if(I0) {
        I0[1] <- runif(1,min=s.range[1],max=s.range[2]/2)
        I0[2] <- runif(1,min=w.range[1],max=w.range[2]/4)
    } else
        I0 <- c(NA,NA) # c(0,0.1^digits)
    type.par <- ifelse(t.par,"t","w")
    type.par <- c(paste(type.par,"0",sep=""),type.par)
    s.val <- t.val <- s.labels <- w.labels <- NULL
    for(i in 1:length(lev)) {
        s.val <- c(s.val,runif(lev[i],min=s.range[1],max=s.range[2]))
        t.val <- c(t.val,runif(lev[i],min=w.range[1],max=w.range[2]))
        s.labels <- c(s.labels,paste("s",LETTERS[i],1:lev[i],sep=""))
        w.labels <- c(w.labels,paste(type.par[2],LETTERS[i],1:lev[i],sep=""))
    }
    if(is.na(I0[1]))
        mt <- mean(t.val)
    else {
        mt <- mean(c(I0[2],t.val))
        I0[2] <- I0[2]-mt
    }   
    t.val <- t.val-mt
    if(!t.par) {
        I0[2] <- exp(I0[2])
        t.val <- exp(t.val)
    }
    param <- round(c(I0,s.val,t.val),digits)
    names(param) <- c("s0",type.par[1],s.labels,w.labels)
    return(param)
}

# Combina un set di parametri generando le risposte R. Se sd>0 e
# trials>1 viene costruita una matrice di R con errore random.
datgen <- function(param, lev, t.par=FALSE, trials=1, sd=0, range=NULL)
{
    fact <- length(lev)
    sumlev <- as.integer(sum(lev))
    sumlev <- c(sumlev,2+sumlev,2+2*sumlev)
    wpos <- c(2,(sumlev[2]+1):sumlev[3])
    if(!t.par)
        param[wpos] <- log(param[wpos])
    if(is.na(param[2])) { # I0=FALSE
        param[1] <- 0
        param[2] <- log(1e-10)
    }
    R <- averaging(param,lev,fact,sumlev)
    R <- matrix(rep.int(R,trials), nrow=trials, byrow=TRUE)
    # Aggiunta di errore alle colonne
    if(sd>0 & trials>1) {
        numcol <- ncol(R)
        j <- 1 # j-esima colonna
        repeat {
            R.err <- R[,j] + rnorm(trials,mean=0,sd=sd)
            check <- TRUE
            if(!is.null(range)) {
                if(any(R.err<range[1]) | any(R.err>range[2]))
                    check <- FALSE
            }
            if(check) {
                R[,j] <- R.err
                j <- j+1
            }
            if(j>numcol) break
        }
    }
    # Attribuzione dei nomi alle colonne
    R.names <- data.names(lev)
    colnames(R) <- R.names
    return(R)
}

# Analisi con rav per soggetto singolo.
rav.single <- function(data,...)
{
    # La funzione lavora solo con data frame:
    if(class(data) == "matrix")
        data <- data.frame(data)
    if(class(data) != "data.frame")
        stop("Data must be matrix or data.frame")
    # Numero corretto di colonne della matrice:
    lev <- list(...)$lev
    all.resp <- respdim(lev)
    # Numero reale di colonne della matrice:
    numCol <- ncol(data)
    # Se manca, aggiungo la colonna per i soggetti:
    if(numCol == all.resp)
        data <- data.frame(subj=seq_len(nrow(data)),data)
    # Se la colonna dei soggetti non e' fattore, la trasformo:
    if(!is.factor(data[,1]))
        data[,1] <- factor(as.character(data[,1]),levels=unique(data[,1]))
    # Analisi per soggetto singolo:
    subj <- levels(data[,1])
    N <- length(subj)
    output <- vector("list",N)
    for(i in 1:N) {
        cat("Subject:",i,"/",N,"\n")
        output[[i]] <- rav(data=data, subset=subj[i], ...)
    }
    names(output) <- subj
    return(output)
}

rav.fitted <- function(object, whichModel=NULL) {
    if(is.null(whichModel))
        whichModel <- slotNames(object)[-(1:10)][object@selection][1]
    return(getElement(object,whichModel)@fitted)
}

standardized <- function(resid)
{
    dim.data <- dim(resid)
    if(dim.data[1]==1)
        stop("It's impossible to calculate standardized residuals: row(data) == 1")
    m <- apply(resid,2,mean,na.rm=TRUE)
    s <- apply(resid,2,sd,na.rm=TRUE)
    m <- matrix(m,nrow=dim.data[1],ncol=dim.data[2],byrow=TRUE)
    s <- matrix(s,nrow=dim.data[1],ncol=dim.data[2],byrow=TRUE)
    resid <- (resid-m)/s
    return(resid)
}

rav.resid <- function(object, whichModel=NULL, standard=FALSE)
{
    if(is.null(whichModel))
        whichModel <- c("null","ESM","SAM","EAM","DAM","IC")[object@selection[1]]
    object <- getElement(object,whichModel)
    resid <- object@residuals
    if(standard)
        resid <- standardized(resid)
    resid <- list(
        matrix = resid,
        rows = apply(resid^2,1,mean,na.rm=TRUE),
        columns = apply(resid^2,2,mean,na.rm=TRUE)
    )
    if(is.null(names(resid$rows)))
        names(resid$rows) <- paste("row",1:length(resid$rows),sep="")
    if(is.null(names(resid$columns)))
        names(resid$columns) <- paste("col",1:length(resid$columns),sep="")
    return(resid)
}

rav.param <- function(object, whichModel=NULL) {
    if(is.null(whichModel))
        whichModel <- c("null","ESM","SAM","EAM","DAM","IC")[object@selection[1]]
    model <- getElement(object,whichModel)
    param <- model@param
    s.labels <- w.labels <- NULL
    if(model@t.par)
        type.par <- c("t0","t")
    else
        type.par <- c("w0","w")
    for(k in 1:object@factors) {
        s.labels <- c(s.labels,paste("s",LETTERS[k],1:object@levels[k],sep=""))
        w.labels <- c(w.labels,paste(type.par[2],LETTERS[k],1:object@levels[k],sep=""))
    }
    names(param) <- c("s0",type.par[1],s.labels,w.labels)
    return(param)
}

rav.AIC <- function(object, whichModel=NULL) {
    if(is.null(whichModel))
        whichModel <- slotNames(object)[-(1:10)][object@selection][1]
    model <- slot(object,whichModel)
    return(model@AIC)
}

rav.BIC <- function(object, whichModel=NULL) {
    if(is.null(whichModel))
        whichModel <- slotNames(object)[-(1:10)][object@selection][1]
    model <- slot(object,whichModel)
    return(model@BIC)
}

outlier.replace <- function(object, whichModel=NULL, alpha=0.05, value=NA)
{
    if(is.null(whichModel))
        whichModel <- slotNames(object)[-(1:10)][object@selection][1]
    resid <- getElement(object,whichModel)@residuals
    resid <- standardized(resid)
    # Calcolo z critico
    N <- sum(!is.na(object@observed))
    cutoff <- qnorm(1-alpha/2)
    # Rimozione outliers
    outliers <- abs(resid) > cutoff
    if(!is.function(value))
        object@observed[which(outliers)] <- value
    else {
        value <- apply(object@observed,2,value,na.rm=TRUE)
        numCol <- ncol(object@observed)
        for(i in 1:numCol)
            object@observed[outliers[,i],i] <- value[i]
    }
    return(object@observed)
}

rav2file <- function(object,what,whichModel=NULL,file=file.choose(),sep=",",dec=".")
{
    if(!is.list(object)) {
        # nel caso object non sia risultato di rav.single
        all.resp <- respdim(object@levels)
        S <- "1" # codice soggetto
        object <- list(subj=object)
    } else {
        S <- names(object)
        all.resp <- respdim(object[[1]]@levels)
    }
    N <- length(object)
    subject <- NULL
    if(is.null(what))
        stop("you must specify the what argument")
    if(!is.null(whichModel))
        model <- whichModel
    if(what=="resid") {
        whichSlot <- "residuals"
        lab <- colnames(object[[1]]@observed)
        output <- matrix(nrow=0,ncol=all.resp)
    }
    if(what=="param") {
        whichSlot <- what
        # Nomi dei parametri
        s.labels <- w.labels <- NULL
        if(object[[1]]@ESM@t.par)
            type.par <- c("t0","t")
        else
            type.par <- c("w0","w")
        for(k in 1:object[[1]]@factors) {
            s.labels <- c(s.labels,paste("s",LETTERS[k],1:object[[1]]@levels[k],sep=""))
            w.labels <- c(w.labels,paste(type.par[2],LETTERS[k],1:object[[1]]@levels[k],sep=""))
        }
        lab <- c("s0",type.par[1],s.labels,w.labels)
        output <- matrix(nrow=0,ncol=2+sum(2*object[[1]]@levels))
    }
    for(i in 1:N) {
        if(is.null(whichModel))
            model <- slotNames(object[[i]])[-(1:10)][object[[i]]@selection][1]
        values <- slot(slot(object[[i]],model),whichSlot)
        if(is.vector(values))
            values <- matrix(values,nrow=1)
        output <- rbind(output,values)
        subject <- c(subject,rep.int(S[i],nrow(values)))
    }
    output <- data.frame(subject=subject,output)
    colnames(output)[-1] <- lab
    write.table(output,file=file,quote=FALSE,sep=sep,na="",dec=dec,row.names=FALSE)
}
