Residual.null <- function(param,fixed=0,model,I0,data,lev,fact,sumlev,sposFirst,dim.data)
{
    .C("residual_null",
        param = as.double(param),
        fixed = as.double(fixed),
        lev = as.integer(lev),
        fact = as.integer(fact),
        sumlev = as.integer(sumlev),
        spos_first = as.integer(sposFirst),
        in_st = as.integer(I0),
        observed = as.double(data),
        dimdata = as.integer(dim.data),
        model = as.integer(model),
        RSS = as.double(0),
        NAOK = TRUE,
        PACKAGE = "rAverage"
    )$RSS
}

Residual.eq <- function(param,fixed,data,lev,fact,sumlev,dim.data)
{
    # Riceve i pesi in t (la averaging li trasformera' in w)
    .C("residual_eq",
        param = as.double(param),
        fixed = as.double(fixed),
        lev = as.integer(lev),
        fact = as.integer(fact),
        sumlev = as.integer(sumlev),
        observed = as.double(data),
        dimdata = as.integer(dim.data),
        RSS = as.double(0),
        NAOK = TRUE,
        PACKAGE = "rAverage"
    )$RSS
}

Residual <- function(param,fixed,data,lev,fact,sumlev,dim.data,Dt,nwfix)
{
    # Riceve i pesi in t (la averaging li trasformera' in w)
    .C("residual",
        param = as.double(param),
        fixed = as.double(fixed),
        lev = as.integer(lev),
        fact = as.integer(fact),
        sumlev = as.integer(sumlev),
        observed = as.double(data),
        dimdata = as.integer(dim.data),
        deltaweights = as.double(Dt),
        numfix = as.integer(nwfix$num),
        valfix = as.double(nwfix$nwval),		
        RSS = as.double(0),
        NAOK = TRUE,
        PACKAGE = "rAverage"
    )$RSS
}

## ###################################################################################

optimization.IC <- function(
    data, fact, lev, sumlev, pos, N, dim.data,  # dati del disegno sperimentale
    model.start, par.base, nwfix, fixed,        # modello di partenza, parametri e fissi
    I0, Dt, IC.diff, all, verbose, IC.break,    # opzioni per la procedura di stima
    lower, upper, method, control, change)      # opzioni per optim (change indica se lo start e' stato modificato)
{
    # Contiene i parametri del modello di partenza per la routine di minimizzazione:
    START <- list(
         param = model.start@param,
           RSS = model.start@RSS,
           IC  = c(model.start@BIC,model.start@AIC),
        n.pars = model.start@n.pars,
          conv = 0,
           msg = ""
    )
    # Conterra' tutti i parametri dell'ultimo modello accettato:
    BEST <- list(
         param = NULL,
           RSS = NULL,
            IC = NULL,
        n.pars = NULL,
          comb = NULL,
          conv = NULL,
           msg = NULL
    )
    # Traccia le combinazioni di pesi selezionate:
    selected <- rep.int(FALSE,sumlev[1])
    # Indica la presenza dei bounds:
    bounds <- method=="L-BFGS-B"
    logN <- log(N) # Serve per il calcolo del BIC

    if(verbose) {
        cat("Model selection by information criterion","\n\n")
        cat("-> start \t",
            "RSS:",sprintf("%.2f",START$RSS),"\t",
            "BIC:",sprintf("%.2f",START$IC[1]),"\t",
            "AIC:",sprintf("%.2f",START$IC[2]),"\n"
        )
    }

    for(i in 1:sumlev[1]) {
        # if(i > 1 & IC.break == TRUE & change == FALSE) break;
        if(verbose) cat("#",i,"free weights\n")
        START$param[pos$fixed] <- fixed[pos$fixed]
        BEST$IC <- START$IC
        accepted <- FALSE # Diventa TRUE se uno dei nuovi modelli sara' accettato
        decision <- " Refused" # Serve per verbose=T
        cc <- combin(sumlev[1],i) # Combinazioni di pesi
        dim.cc <- dim(cc)
        # --------------------------------------
        # Fissaggio dei pesi
        # --------------------------------------
        # Se alcuni pesi sono stati fissati non possono essere modificati dalla
        # routine. Si devono quindi eliminare da cc le righe con i loro indici.
        if(nwfix$num) {
            # Matrice che conterra' le combinazioni di pesi che saranno selezionate
            # (ovvero le combinazioni che non contengono pesi fissi):
            cc.sel <- matrix(FALSE,dim.cc[1],dim.cc[2])
            # Si confronta la matrice cc con il vettore nwfix$pos. cc.sel e' una matrice di
            # valori logici con le stesse dimensioni di cc. I TRUE identificheranno
            # i punti di cc che contengono gli indici dei pesi fissati e poi eliminati:
            for(j in 1:nwfix$num)
                cc.sel <- cc.sel + (cc==nwfix$pos[j])
            cc.sel <- which(rowSums(cc.sel)==0)
            cc <- cc[cc.sel,]
            # Se cc e' diventato un vettore bisogna riconvertirlo in matrice. Attenzione:
            # nel caso i=1, cc deve essere un vettore colonna, in tutti gli altri casi
            # un vettore riga.
            if(is.vector(cc)) {
                cc <- t(cc)
                if(i==1) cc <- t(cc)
            }
            dim.cc[1] <- nrow(cc)
            # Se tutte le combinazioni sono state eliminate, o se le rimanenti hanno 
            # meno parametri di quanto indicato dalla i si salta al successivo ciclo.
            if(length(cc)==0 | dim.cc[2]<i) next
        }
        
        # --------------------------------------
        # Selezione delle combinazioni
        # --------------------------------------
        # Quando all=F la cc contiene solo le combinazioni con i pesi che hanno portato
        # miglioramenti nel fit. nrow(cc)==1 e' il differential weight case e si omette.
        if(all==FALSE & (i>1 & i<sumlev[1]-1) & dim.cc[1]>1) {
            if(i>1 & dim.cc[1]>1) {
                cc.sel <- FALSE
                elem <- length(vet<-(1:sumlev[1])[selected])
                if(elem!=0) {
                    for(k in 1:elem) cc.sel <- cc.sel+(cc[,]==vet[k])
                    cc <- as.matrix(cc[rowSums(as.matrix(cc.sel))==elem,])
                }
            }
            dim.cc[1] <- nrow(cc)
        }
        
        # --------------------------------------
        # Stima per le date combinazioni
        # --------------------------------------
        for(j in 1:dim.cc[1]) {
            # Fissaggio dei pesi:
            eachfixed <- fixed
            eachfixed[pos$wpos[-cc[j,]]] <- START$param[pos$wpos[-cc[j,]]]
            
            # Stima dei parametri:
            output <- optim(par=START$param, fn=Residual, fixed=eachfixed, data=data,
                lev=lev, fact=fact, sumlev=sumlev, dim.data=dim.data, Dt=Dt,
                nwfix=nwfix[c(1,3)], method=method, lower=lower, upper=upper, control=control
            )
            
            # Si scalano i t
            # Attenzione: scalare i t in questo punto potrebbe risultare critico. Se alcuni t vengono
            # compressi contro il limite superiore e altri contro il limite inferiore, potrebbe
            # accadere che la correzione t = t-mean(t) faccia oltrepassare i bound ai t. Se usati
            # come start nel ciclo successivo, questi t andranno a genererare un errore nella optim.
            if(!I0)
                output$par[pos$wpos] <- output$par[pos$wpos]-mean(output$par[pos$wpos])
            else
                output$par[pos$w0wpos] <- output$par[pos$w0wpos]-mean(output$par[pos$w0wpos])
            
            # Sostituzione dei parametri fissi a quelli in output:
            output$par[pos$fixed] <- fixed[pos$fixed]
            
            # I pesi uguali entro un certo delta vengono eguagliati:
            output$par <- parmeanlast(output$par,fixed,sumlev,Dt,nwfix)
            
            # Calcolo del numero di parametri:
            if(nwfix$num==0)
                output$n.pars <- par.base+numpar(output$par[pos$wpos],sumlev[1])-1
            else
                output$n.pars <- par.base+numpar(output$par[pos$wpos[-nwfix$pos]],sumlev[1]-nwfix$num)-1
            
            # Dato che i parametri potrebbero essere stati modificati si ricalcola l'RSS:
            output$value <- Residual(output$par,fixed,data,lev,fact,sumlev,dim.data,Dt,nwfix)
            
            # Calcolo degli indici AIC e BIC
            output$IC <- c(
                N*log(output$value/N)+output$n.pars*logN, # BIC
                N*log(output$value/N)+2*output$n.pars # AIC
            )
            if(N/output$n.pars < 40) {
                output$IC[1] <- output$IC[1]+(output$n.pars*logN*(output$n.pars+1))/(N-output$n.pars-1) # BICc
                output$IC[2] <- output$IC[2]+(2*output$n.pars*(output$n.pars+1))/(N-output$n.pars-1) # AICc
            }
            
            # Valutazione del modello stimato:
            if((output$IC[1]+IC.diff[1]<BEST$IC[1]) |
            ((output$IC[1]<BEST$IC[1]) & (output$IC[2]+IC.diff[2]<BEST$IC[2]))) {
                # Il nuovo modello viene accettato come migliore:
                BEST$param  <- output$par
                BEST$RSS    <- output$value
                BEST$IC     <- output$IC
                BEST$n.pars <- output$n.pars
                BEST$comb   <- cc[j,]
                BEST$conv   <- output$convergence
                BEST$msg    <- output$message
                accepted    <- change <- TRUE
                decision <- "Accepted"
            }
            
            if(verbose) {
                if(output$convergence==0)
                    convergence <- "converged"
                else
                    convergence <- "not converged"
                cat("Free W:",cc[j,],"\t",
                    "RSS:",sprintf("%.2f",output$value),"\t",
                    "BIC:",sprintf("%.2f",output$IC[1]),"\t",
                    "AIC:",sprintf("%.2f",output$IC[2]),"\t",
                    decision,"\t", convergence,"\n"
                )
                decision <- " Refused"
            }
        }
        if(accepted) {
            # Ultimata la combinazione di pesi j, il miglior modello diventa
            # il modello di partenza per la successiva combinazione:
            START <- BEST
            selected[BEST$comb] <- TRUE
        }
    }
    if(verbose)
        cat("-> select","\t",
            "RSS:",round(START$RSS,2),"\t",
            "BIC:",round(START$IC[1],2),"\t",
            "AIC:",round(START$IC[2],2),"\n\n"
        )
    return(
        list(
            par = START$param,
            value = START$RSS,
            convergence = START$conv,
            message = START$msg,
            n.pars = START$n.pars,
            change = change
        )
    )
}
