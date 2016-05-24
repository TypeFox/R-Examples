intercept<-function (ogg, parm, gap=TRUE, rev.sgn = FALSE, var.diff = FALSE, 
    digits = max(3, getOption("digits") - 3)){
#corregge in caso di no model intercept -- CHE VOLEVO DIRE?? #forse che adesso funziona se nel modello non c'e' l'interc.
#--
        f.U<-function(nomiU, term=NULL){
        #trasforma i nomi dei coeff U (o V) nei nomi delle variabili corrispondenti
        #and if 'term' is provided (i.e. it differs from NULL) the index of nomiU matching term are returned
            k<-length(nomiU)
            nomiUsenzaU<-strsplit(nomiU, "\\.")
            nomiU.ok<-vector(length=k)
            for(i in 1:k){
                nomi.i<-nomiUsenzaU[[i]][-1]
                if(length(nomi.i)>1) nomi.i<-paste(nomi.i,collapse=".")
                nomiU.ok[i]<-nomi.i
                }
          if(!is.null(term)) nomiU.ok<-(1:k)[nomiU.ok%in%term]
          return(nomiU.ok)
        }
#--        
    #if (!"segmented" %in% class(ogg)) stop("A segmented model is needed")
    if (var.diff && length(ogg$nameUV$Z) > 1) {
        var.diff <- FALSE
        warning("var.diff set to FALSE with multiple segmented variables",
            call. = FALSE)
    }
    nomepsi <- rownames(ogg$psi)
    nomeU <- ogg$nameUV[[1]]
    nomeZ <- ogg$nameUV[[3]]
    if (missing(parm)) {
        nomeZ <- ogg$nameUV[[3]]
        if (length(rev.sgn) == 1)
            rev.sgn <- rep(rev.sgn, length(nomeZ))
    }
    else {
        if (!all(parm %in% ogg$nameUV[[3]])) {
            stop("invalid parm")
        }
        else {
            nomeZ <- parm
        }
    }
    if (length(rev.sgn) != length(nomeZ))
        rev.sgn <- rep(rev.sgn, length.out = length(nomeZ))
    nomi <- names(coef(ogg))
    nomi <- nomi[-match(nomepsi, nomi)]
    Allpsi <- index <- vector(mode = "list", length = length(nomeZ))
    gapCoef<-summary.segmented(ogg)$gap
    Ris <- list()
    rev.sgn <- rep(rev.sgn, length.out = length(nomeZ))
    if("(Intercept)"%in%names(coef(ogg))){
        alpha0 <- alpha00 <- coef(ogg)["(Intercept)"]} else {alpha0 <- alpha00 <-0}
#per ogni variabile segmented...
    for (i in 1:length(nomeZ)) {
#        id.cof.U <- grep(paste("\\.", nomeZ[i], "$", sep = ""), nomi, value = FALSE)
#        psii <- ogg$psi[grep(paste("\\.", nomeZ[i], "$", sep = ""), rownames(ogg$psi), value = FALSE), 2]
        id.cof.U <- f.U(ogg$nameUV$U, nomeZ[i]) + (match(ogg$nameUV$U[1], nomi)-1)
        psii<- ogg$psi[f.U(ogg$nameUV$V, nomeZ[i]) , "Est."] 
        Allpsi[[i]] <- sort(psii, decreasing = FALSE) 
        id.cof.U <- id.cof.U[order(psii)]
        index[[i]] <- id.cof.U
        alpha0<-if("(Intercept)"%in%names(coef(ogg))) coef(ogg)["(Intercept)"] else 0
        ind <- as.numeric(na.omit(unlist(index[[i]])))
        cof <- coef(ogg)[ind]
        alpha <- vector(length = length(ind))
        #gapCoef.i<-gapCoef[grep(paste("\\.",nomeZ[i],"$",sep=""), rownames(gapCoef), value = FALSE),"Est."]
        gapCoef.i<-gapCoef[f.U(rownames(gapCoef), nomeZ[i]) ,"Est."]
        for (j in 1:length(cof)) {
            alpha[j] <- alpha0 - Allpsi[[i]][j] * cof[j]
            if(gap) alpha[j] <- alpha[j] - gapCoef.i[j]
            alpha0 <- alpha[j]

        }
        #if(gap) alpha<-alpha -gapCoef[grep(paste("\\.",nomeZ[i],"$",sep=""), rownames(gapCoef), value = FALSE),"Est."]
        cof.out <- c(alpha00, alpha)
        if(rev.sgn[i]) cof.out <- cof.out[length(cof.out):1]
        ris <- matrix(cof.out)
        dimnames(ris) <- list(paste("intercept", 1:nrow(ris), sep = ""), "Est.")
        Ris[[nomeZ[i]]] <- signif(ris, digits)
    }
    Ris
}
