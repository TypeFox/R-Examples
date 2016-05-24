"factorWW" <- function(p, WW, spalten = NULL, andere = NULL, ngroups=NULL) {
#cat("sp",spalten,"ng:","\n")
#for (i in 1:length(ngroups)) cat(ngroups[[i]],"\n")

# Diese Funktion berechnet den Faktor fuer eine Modellpermutation,
# also z.B. wie haeufig die erste Variable an zweiter Stelle mit
# vorangegeangener zweiter Variable auftritt.
# Alle Bedingungen muessen so gestellt werden, das die Hoerwertigen
# zuletzt kommen, also z.B. alle zweifach WW vor den einfach WW.
# WW = Wechselwirkung, HK = Hauptkompunenten

## bei Gruppen muss man die Effektnamen 
## der richtigen Position zuordnen

 if (!is.null(ngroups)){
## Aufgabe von ngroups ist, einen Bezug auf die richtigen Positionen 
##    in den Spalten der Matrizen in indices herzustellen. 
##    Dort sind die Gruppen immer vorn (auch wenn aufgrund von Faktoren und WWen
##    entstanden.
## ngroups enthaelt an erster Position den Effekt, der der ersten Gruppe entspricht
##    usw., und erst nach den mehrspaltigen Effekten folgen die einspaltigen Effekte
##    in ihrer natuerlichen Reihenfolge.
## Dies funktioniert derzeit nur, wenn alle Elemente von ngroups Zahlen 
## (keine Vektoren) sind.
## Falls echte Gruppen auftreten, liegen Vektoren vor; was dann?
##    --> bisher durch Abbruchbedingung abgefangen.

    if (!is.null(spalten)) spalten <- list2vec(ngroups[spalten])-1
    if (!is.null(andere)) andere <- list2vec(ngroups[andere])-1
 }

 # An erster Stelle (spalten) darf niemals eine WW stehen          
 # wird benoetigt, weil davor und danach leer sein koennten
 if (any(spalten[1] == WW[[1]][,1])) {
    erg <- 0
    div.z <- 1
    div.n <- 1
    davor <- NULL
    danach <- NULL
 } else {
    davor <- andere ## sort necessary because of reordering
                    ## so that change of comb.WW works
    if (!is.null(davor)) davor <- sort(davor,method="quick")
    if (is.null(spalten)) {
        danach <- 1:p
    } else {    
        danach <- (1:p)[-spalten]
        if(length(danach)==0) danach <- NULL
    } 

    WW.sa  <- list("davor" = davor, "danach" = danach)
    # Bestimmung der Zusammenhaenge (Farbkombination)
    WW2 <- WW
    new.order <- 1:p
    for (i in 1:nrow(WW2$WW)){
        if (new.order[WW2$WW[i,2]] == WW2$WW[i,2]) {
            new.order[WW2$WW[i,2]] <- WW2$WW[i,1]
        } else {
            new.order[new.order==WW2$WW[i,1]] <- new.order[WW2$WW[i,2]]
            WW2$WW[,1][WW2$WW[,1] == WW2$WW[i,1]] <- new.order[WW2$WW[i,2]]
        }
     }
    new.order.v <- sort(unique(new.order))

    ## vorlaeufige Werte
    erg <- factorial(length(davor)) * factorial(length(danach))
    div.z <- c(1,1)
    div.n <- c(1,1)

    # kombinationen davor und danach bearbeiten, wenn Elemente enthalten
    for(i in 1:2) {

         ## erstes Element von WW.sa: davor, zweites danach
        if (is.null(WW.sa[[i]])) next
        comb <- new.order[WW.sa[[i]]]  ## davor-Elemente bzw. danach-Elemente

        comb.table <- table(comb)
        # vorlaeufiger Teiler (innerhalb der Gruppen nur eine festgelegte Reihenfolge)
        div.z[i] <- prod(factorial(as.numeric(comb.table)))
        
        # wird wieder verkleinert durch Wegnahme interner Reihenfolgen
        # mehrere Elemente einer Gruppe sind vorhanden
        if (any(comb.table > 1)) {
            comb.gto <- as.numeric(names(comb.table)[comb.table > 1])
            ## mehrfach vorhandene Elemente
            for (j in comb.gto) {
                ## in dieser Gruppe vorkommende Effekte
                comb.s <- WW.sa[[i]][comb == j]
                jetzt <- t(permutations(length(comb.s)))
                comb.WW <- WW$WW
                comb.s.WW <- comb.s[which(comb.s %in% comb.WW[,1])]
                comb.s.HW <- comb.s[which(comb.s %in% comb.WW[,2])]
                ## Nullsetzen der fuer diese Gruppe nicht relevanten Effekte
                comb.WW[!(comb.WW %in% comb.s)] <- 0
                for (k in 1:length(comb.s)){
                    ## comb.WW darf nur maximal-Index=Anzahl verschiedener Indizes haben
                    ## WWen brauchen hoehere Indizes als Hauptwirkungen
                    ## Werte in comb.WW durch ihre Position in comb.s ersetzt
                    if (comb.s[k] %in% comb.s.HW) comb.WW[comb.WW==comb.s[k]] <- 
                          which(comb.s.HW==comb.s[k])
                    else comb.WW[comb.WW==comb.s[k]] <-
                          length(comb.s.HW) + which(comb.s.WW==comb.s[k])
                    }
                comb.WW <- matrix(comb.WW[apply(comb.WW != 0,1,all),], ncol=2)
                ## Nullzeilen geloescht
                if (length(comb.WW) > 0)
                    div.n[i] <- div.n[i]*sum(apply(jetzt,2,checkWW,WW = comb.WW,ord=TRUE))
                else div.n[i] <- div.n[i]*factorial(length(comb.s))
                    ## Multiplikator: how many possibilities for orders of equally-ranked effects
                    ## (before = davor and after = danach)
            }
        }
    }
 } 
 return((erg*prod(div.n)/prod(div.z)))
}