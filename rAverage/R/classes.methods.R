setClass("indices",
    slots = c(
        param="numeric",
        I0="vector",
        levels="numeric",
        start="numeric",
        upper="numeric",
        lower="numeric",
        fitted="vector",
        residuals="matrix",
        AIC="numeric",
        BIC="numeric",
        RSS="numeric",
        TSS="numeric",
        R2="numeric",
        adjR2="numeric",
        n.pars="numeric",
        t.par="logical",
        title="character",
        names="character",
        message="character",
        model="numeric"
    )
)

setClass("rav",
    slots = c(
        observed="matrix",
        factors="numeric",
        levels="numeric",
        title="character",
        names="character",
        method="character",
        start="numeric",
        lower="numeric",
        upper="numeric",
        selection="numeric",
        null="indices",
        ESM="indices",
        SAM="indices",
        EAM="indices",
        DAM="indices",
        IC="indices"
    )
)

setMethod("show","indices",
    function(object) {
        cat(object@title)
        if(length(object@message>0)) {
            if(object@message=="MNA") {
                # Significa che il modello IC non e' stato stimato
                cat("Model Not Available\n")
                return()
            }
        }
        fact <- length(object@levels)
        if(!object@t.par)
            wlab <- c("  w0:","w:","w"," w")
        else
            wlab <- c("  t0:","t:","t"," t")
        if(object@I0 & object@model!=1) {
            cat("Initial state:")
            cat("  s0:",sprintf("%.2f",object@param[1]))
            cat(wlab[1],sprintf("%.2f",object@param[2]))
            cat("\n")
        }
        if(object@model==1) {
            cat("s:", sprintf("%.2f",object@param[3]),"\t")
            cat(wlab[2],sprintf("%.2f",object@param[3+sum(object@levels)]),"\n")
        } else {
            param <- list(I0=object@param[1:2],s=NULL,w=NULL)
            param$s <- object@param[3:(2+sum(object@levels))]
            param$w <- object@param[(3+sum(object@levels)):length(object@param)]
            if(object@model<=4) {
                if(object@model==2) {
                    s <- param$s[1+cumsum(c(0,object@levels[-fact]))] # Il primo s di ogni fattore
                    w <- param$w[1+cumsum(c(0,object@levels[-fact]))] # Il primo w di ogni fattore
                } else {
                    s <- list(param$s[1:object@levels[1]]) # Tutti gli s del primo fattore
                    w <- list(param$w[1]) # Il primo w del primo fattore
                    for(i in 2:fact) { # Si creano degli slot per i parametri degli altri fattori
                        s[[i]]<-param$s[(sum(object@levels[1:(i-1)])+1):sum(object@levels[1:i])]
                        w[[i]]<-param$w[(sum(object@levels[1:(i-1)])+1)]
                    }
                }
            } else {
                s <- list(param$s[1:object@levels[1]]) # Tutti gli s del primo fattore
                w <- list(param$w[1:object@levels[1]]) # Tutti i w del primo fattore
                for(i in 2:fact) {
                    s[[i]]<-param$s[(sum(object@levels[1:(i-1)])+1):sum(object@levels[1:i])]
                    w[[i]]<-param$w[(sum(object@levels[1:(i-1)])+1):sum(object@levels[1:i])]
                }
            }
            # Le stringhe dei nomi dei fattori vengono rese a lunghezza uguale.
            # Questo viene fatto aggiungendo spazi vuoti alle stringhe piu' corte.
            empty <- paste(rep(" ",50),collapse=" ") # emorme stringa vuota
            # numero di caratteri di ogni stringa:
            charNames <- nchar(object@names)
            # massimo numero di caratteri in una stringa:
            maxChar <- max(charNames)
            # posizioni occupate dalle stringhe piu' lunghe:
            posMaxChar <- which(charNames==maxChar)
            # posizioni occupate dalle stringhe piu' corte:
            posMinChar <- (1:fact)[-posMaxChar]
            # Le stringhe, se gia' non lo sono, vengono rese a lunghezza uguale:
            if(length(posMinChar)!=0) {
                for(i in 1:length(posMinChar)) {
                    # numero di caratteri vuoti da aggiungere:
                    ladd <- maxChar-charNames[posMinChar[i]]
                    # Amplia le stringhe piu' corte:
                    object@names[posMinChar[i]]<-
                        paste(object@names[posMinChar[i]],substr(empty,1,ladd),sep="")
                }
            }
            # Nomi dei parametri:
            maxLev <- max(object@levels)
            labels <- list(s=NULL,w=NULL)
            # I parametri, per essere stampati a video, devono essere trasformati in stringhe
            # a uguale numero di caratteri. Questo è il numero di caratteri che deve possedere
            # ogni stringa:
            charStr <- max(nchar(as.integer(object@param))+3)
            # (+3 perche' ci sono da aggiungere i caratteri del punto e delle due cifre decimali)
            # Costruzione etichette:
            ladd <- substr(empty,1,charStr-2)
            if(object@model==2) {
                labels$s <- paste(ladd," s",sep="")
                labels$w <- paste(ladd,wlab[3],sep="")
            } else {
                labels$s <- c(paste(ladd,"s",1,sep=""))
                if(object@model<=4)
                    labels$w <- c(paste(ladd,wlab[4],sep=""))
                else
                    labels$w <- c(paste(ladd,wlab[3],1,sep=""))
                for(i in 2:maxLev) {
                    labels$s <- c(labels$s,c(paste(ladd,"s",i,sep="")))
                    if(object@model>4)
                        labels$w <- c(labels$w,c(paste(ladd,wlab[3],i,sep="")))
                }
            }
            
            # Se i fattori hanno numeri diversi di livelli, in corrispondenza
            # dei valori mancanti bisogna mettere degli spazi vuoti (NA):
            if(object@model>2) {
                for(i in 1:fact) {
                    if(object@levels[i]<maxLev) {
                        s[[i]]<-c(s[[i]],rep.int(NA,maxLev-object@levels[i]))
                        if(object@model>4)
                            w[[i]]<-c(w[[i]],rep.int(NA,maxLev-object@levels[i]))
                    }
                }
            }
            
            # Il comando sprintf lavorara in modo diverso a seconda che il parametro
            # sia un numero oppure un NA. Bisogna utilizzare due diverse regole di stampa:
            printRuleNA <- paste("%",charStr,"s",sep="")
            printRule <- paste("%",charStr,".2f",sep="")
            
            # Stampa le etichette contenenti i nomi dei parametri:
            spacer <- paste(substr(empty,1,maxChar),"  ",sep="")
            cat(spacer,labels$s,"",labels$w,"\n")
            
            # Questo ciclo stampa a video i parametri, inserendo una stringa vuota
            # in corripondenza dei valori 'NA':
            parlen <- max(object@levels)
            if(object@model==2) {
                for(i in 1:fact) {
                    cat(object@names[i],"  ")
                    cat(sprintf(printRule,c(s[i],w[i])),"\n")
                }
            } else {
                for (i in 1:fact) {
                    par <- unlist(c(s[i],w[i]))
                    cat(object@names[i],"  ")
                    for(j in 1:(2*parlen)) {
                        if(!is.na(par[j]))
                            cat(sprintf(printRule,par[j]),"")
                        else
                            cat(sprintf(printRuleNA,""), "")
                        if(j==parlen) cat(" ")
                        if(j==(2*parlen)) cat("\n")
                    }
                }
            }
        }
        
        # Stampa a vieo gli indici di fit:
        
        Rsq <- c(object@R2,object@adjR2)
        Rsq[which(Rsq<0)] <- NA
        cat("---\n")
        cat("AIC:", sprintf("%.2f",object@AIC),"",
            "BIC:", sprintf("%.2f",object@BIC),"",
            "n.pars:",object@n.pars, "\n")
        cat("R-squared:", sprintf("%.4f",Rsq[1]),"",
            "Adjusted R-squared:", sprintf("%.4f",Rsq[2]), "\n")
        cat("Residual Sum Sq:", sprintf("%.2f",sum(object@RSS)),"",
            "Total Sum Sq:", sprintf("%.2f",sum(object@TSS)), "\n")
        cat("\n")
        if(object@message!="")
            cat(object@message,"\n")
    }
)

setMethod("show","rav", 
    function(object) {
        cat("Parameter Estimation for Averaging Model\n")
        cat("Optimization algorithm:",object@method,"\n")
        design <- NULL
        for(i in 1:object@factors)
            design <- c(design,paste("(",object@levels[i]," lev",")",sep=""))
        design <- paste(object@names,design)
        for(i in 2:object@factors)
            design[1] <- paste(design[1],design[i],sep=" x ")
        cat(design[1],"\n")
        if(length(object@title)!=0) cat(object@title,"\n")
        cat("\nBEST MODEL(S) SELECTED:\n\n")
        for(i in 1:length(object@selection)) {
            if (object@selection[i]==1) print(object@null)
            if (object@selection[i]==2) print(object@ESM)
            if (object@selection[i]==3) print(object@SAM)
            if (object@selection[i]==4) print(object@EAM)
            if (object@selection[i]==5) print(object@DAM)
            if (object@selection[i]==6) print(object@IC)
        }
    }
)

setMethod("summary","rav", 
    function(object) {
        cat("Parameter Estimation for Averaging Model\n")
        cat("Minimization algorithm:",object@method,"\n")
        design <- NULL
        for(i in 1:object@factors)
            design <- c(design,paste("(",object@levels[i]," levels",")",sep=""))
        design <- paste(object@names,design)
        for(i in 2:object@factors)
            design[1] <- paste(design[1],design[i],sep=" x ")
        cat(design[1],"\n")
        if(length(object@title)!=0) cat(object@title,"\n")
        cat("\n")
        cat("FITTED MODELS\n\n")
        print(object@null);cat("\n")
        print(object@ESM); cat("\n")
        print(object@SAM); cat("\n")
        print(object@EAM); cat("\n")
        print(object@DAM); cat("\n")
        print(object@IC);  cat("\n")
        cat("Best Model(s):", modelNames[object@selection],"\n")
    }
)

modelNames <- c("Null Model", "Equal scale values model", "Simple averaging model",
    "Equal-weights averaging model", "Differential-weights averaging model", "Information criteria")
