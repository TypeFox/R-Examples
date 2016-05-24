fixparam <- function(lev,type.par="w",names=NULL)
{
    # Variabile che contiene il nome dell'environment:
    envirFun <- environment()
    # Variabile che conterra' i pesi fissi:
    par.fixed <- FALSE
    # (e' necessario definirla subito perche' alla fine la variabile
    # e' richiesta dal return: se l'utente ha chiuso la finestra senza
    # aver premuto sul pulsante 'Fix' la funzione deve comunque ritornare
    # qualcosa - NULL in questo caso).
    if(is.null(names)) {
        names <- NULL
        for(i in 1:length(lev))
            names <- c(names,paste("Factor",LETTERS[1:lev[i]]))
    }
    # Crea una finestra toplevel:
    top <- tktoplevel()
    tktitle(top) <- paste(type.par,"-fixing",sep="")
    # Frame per le celle:
    frameCell <- tkframe(top)
    tkgrid(frameCell)
    # Frame per il pulsante:
    frameBut <- tkframe(top,pady=10)
    tkgrid(frameBut)
    # Mostra sul toplevel i nomi dei pesi:
    parLabel <- ""
    for(i in 1:(max(lev)+1)) {
        label <- tklabel(frameCell,text=parLabel)
        tkgrid(label,column=i,row=1)
        parLabel <- paste(type.par,i,sep="")
    }
    # Mostra sul toplevel i nomi dei fattori:
    for(i in 2:(length(lev)+1)) {
        label <- tklabel(frameCell,text=names[i-1])
        tkgrid(label,row=i,column=1)
    }
    # Creazione e display delle celle. N.B.: a ogni cella deve corrispondere
    # una variabile tcl. Queste variabili sono generate dinamicamente, fondendo 
    # il nome 'cellVar' con un numero prograssivo. Alla fine del ciclo la lista
    # entryVar conterra'  tutti i nomi delle variabili; i valori di ogni variabile 
    # potranno essere estratti col comando tclvalue (questo verra' fatto dalla
    # funzione SendFixing).
    entryVar <- entryCell <- list()
    fact <- 1
    col <- row <- 2
    for(i in 1:sum(lev)) {
        entryVar[[i]] <- paste("var",i,sep="")
        entryCell[[i]] <- paste("cell",i,sep="")
        assign( entryVar[[i]], tclVar() )
        cellVar <- eval(parse(text=entryVar[[i]]))
        assign( entryCell[[i]],tkentry(frameCell,textvariable=cellVar,width=5) )
        cell <- eval(parse(text=entryCell[[i]]))
        tkgrid(cell,row=row,column=col)
        if((lev[fact]+1) == col) {
            fact <- fact+1
            row <- row+1
            col <- 2
        } else col <- col+1
    }
    # Funzione per l'estrazione dei valori delle variabili:
    sendFixing <- function() {
        par.fixed <- pairlist()
        values <- NULL
        fact <- i <- 1
        repeat {
            values <- c(values,tclvalue(eval(parse(text=entryVar[[i]]))))
            if(i == cumsum(lev)[fact]) {
                values <- suppressWarnings(as.numeric(values))
                values[which(is.na(values))] <- NA
                par.fixed[[fact]] <- values
                values <- NULL
                fact <- fact+1
                }
            if(i == cumsum(lev)[length(lev)]) break
            i <- i+1
        }
        tkdestroy(top)
        assign("par.fixed",par.fixed,pos=envirFun)
    }
    # Crea il pulsante:
    fixBut <- tkbutton(frameBut,text="Fix",command=sendFixing)
    # Mostra a video il pulsante:
    tkgrid(fixBut,sticky="e")
    tkwait.window(top)
    # Se la lista e' tutta di NA/NaN si rimette par.fixed a FALSE
    if(is.list(par.fixed) & sum(!is.na(unlist(par.fixed)))==0)
        par.fixed <- FALSE
    return(par.fixed)
}
