boot.relimp.lm <- function(object, type = "lmg", groups = NULL, groupnames=NULL, always = NULL, ..., b=1000){
    ## error checking through calc.relimp
    hilf <- calc.relimp(object, type=type, groups=groups, always=always, ...)
    rm(hilf)

    lm<-object
    weights <- lm$weights
    terms <- lm$terms
    resp <- attr(terms,"response")
    if ((!is.null(groups)) && !is.list(groups)) 
       groups<-list(groups)
       ogroups <- groups
       ngroups <- groups

        ## interim value of p, since numbers in always and groups would refer to these columns
    ## in case of factors, this is not identical to the later p
    p <- length(labels(terms))
    
    if (is.null(always)) numalw <- integer(0)
        if (is.character(always)){ 
            numalw <- which(labels(terms) %in% always) + 1
            }
        if (is.numeric(always)) {
              numalw<-always
              always <- labels(terms)[always-1]
            }
        ## now always consists of names of model terms or is NULL

        WW <- NULL   
    ## WW defined in case of no WW'en
    # handle linear model with higher order terms
    if (max(attr(terms,"order")) >= 2) {
      if (sum(attr(lm$terms,"order") > 1) > 1) {
        WW.mat <- attr(lm$terms,"factors")[-1,attr(lm$terms,"order") > 1]
      } 
      else {
        WW.mat <- as.matrix(attr(lm$terms,"factors")[-1,attr(lm$terms,"order") > 1])
        colnames(WW.mat) <- colnames(attr(lm$terms,"factors"))[attr(lm$terms,"order") > 1]
      }
      WW.HK <- sum(attr(lm$terms,"order") == 1)
      WW <- cbind("posWW" = col(WW.mat)[WW.mat == 1]+WW.HK,
          "posHK" = row(WW.mat)[WW.mat == 1])
      order <- attr(lm$terms, "order")
      WW <- list("WW" = WW, "Order" = order) 
      
      WW$WW <- WW$WW[,c(1,2,2)]
      DATA = cbind(lm$model[,1], model.matrix(lm)[,-1])  ## is a matrix
      colnames(DATA)[1] <- colnames(lm$model)[1]
      
      # Positionen der Faktoren
      g_facnames <- character(0)
      factor.var <- which(as.vector(attr(lm$terms, "dataClasses") == "factor"))-1
         ## groups is list or NULL, make it numeric list
        if (!is.null(groups)) {
           groups <- lapply(groups, function(obj){
               if (is.character(obj)) which(labels(terms) %in% obj) + 1 else obj
               })
           ngroups <- groups}

      # Wenn Faktor als Hauptwirkung auftritt 
      # oder bei WWen irgendwelche Faktoren mit mehr als 2 Stufen auftreten, 
      # bearbeiten...
      if (any(factor.var %in% WW[[1]][,2]) || length(lm$assign) > length(WW[[2]])+1) {
        if (length(lm$assign) > length(WW[[2]])+1){
          hilf <- matrix(0,0,3)
          for (we in unique(WW$WW[,1])){
                hes <- which(lm$assign %in% WW$WW[WW$WW[,1]==we,2]) - 1
                hess <- lm$assign[hes+1]
                hilf <- rbind(hilf,matrix(c(rep(we,length(hes)),hes,hess),length(hes),3,byrow=FALSE))
             }
             WW$WW <- hilf
          }
        y.modelm <- model.matrix(lm)
        # Gruppenvariablen hinzufugen, wenn noetig
        g_names <- attr(y.modelm, "dimnames")[[2]]
        g_assign <- attr(y.modelm, "assign")
        g_facpos <- as.data.frame(table(g_assign))
        ## make ogroups contain the correct column numbers
        ## allowing also to group factors
        if (!is.null(groups)) {
           ogroups <- lapply(groups, function(obj){
               c(which(g_assign %in% (obj-1)))
               })
           groups <- lapply(ogroups, function(obj){
               g_names[obj]
               })
           }
        for (k in 2:length(g_facpos[,1])) {
          ##k=1 is intercept
          fac.name <- labels(terms)[as.numeric(g_facpos[k,1])-1]
          fac.varn <- g_names[g_assign == g_facpos[k,1]]
          if (g_facpos[k,2]==1) {
              colnames(DATA)[g_assign == g_facpos[k,1]]<-fac.name
              ## make groups findable in case of two levels only
              if (!is.null(groups)) groups <- lapply(groups, function(obj){
                    if (fac.varn %in% obj) obj[obj==fac.varn]<-fac.name
                    obj
                 })
              next
              }
          if (any(always %in% fac.name)) {
            always <- c(always[always != fac.name], fac.varn)
                 ## darf nachher nicht mehr als Faktor auftauchen
            next
          } 
          ## take care of the case with groups including a factor with more than two levels
           if (!is.null(groups)) {
                   schonda <- FALSE
                   for (obj in groups) {
                        if (all(fac.varn %in% obj)) schonda <- TRUE
                        }
                   if (schonda) next
              }
          g_facnames <- c(g_facnames,fac.name)
           ## groups ist NULL oder character-liste
          if (is.null(groups)) {
              groups <- list(fac.varn)
              ngroups <- list(which(labels(terms) %in% fac.name)+1)
          }
              else {groups[[length(groups)+1]] <- fac.varn
              ngroups[[length(ngroups)+1]] <- which(labels(terms) %in% fac.name)+1
              }
            }       ## end of for loop
            
          if (!is.null(groups)) groups <- lapply(groups, function(obj){which(colnames(DATA) %in% obj)})
          # groups could be NULL if only two-level-factors
          # now character variables
            if (is.null(ogroups)) groupnames <- g_facnames
            else {if(is.null(groupnames)) 
                   groupnames <- c(paste("G", 1:length(ogroups), sep=""), g_facnames)
                 else groupnames <- c(groupnames, g_facnames)}
      }
      hilf <- which(colnames(DATA) %in% always)-1 
           ## numeric version of always compatible with second column of WW
      ## remove all lower order effects from WW that are in always

      if (any(hilf %in% WW$WW[,2])) {
         for (we in intersect(hilf,unique(WW$WW[,2])))
            WW$WW <- WW$WW[which(!WW$WW[,2]==we),]
            if (!(is.matrix(WW$WW))) WW$WW <- matrix(WW$WW,1,3)  
                 ## einzeilig ist es keine Matrix,
                 ## nullzeilig schon
         if (nrow(WW[[1]])==0) WW <- NULL
         }
      ## give error if remaining interactions in WW that are in always
         if (!is.null(WW)){
         if (any(hilf %in% which(lm$assign %in% WW$WW[,1])))
             stop("Interaction must not be in always, \n 
                 if not all corresponding lower level effects are also in always.")
         }
      ## Einzeleffekte wieder aus WW$WW entfernen
      WW$WW <- WW$WW[,c(1,3)]
      if (!(is.matrix(WW$WW))) WW$WW <- matrix(WW$WW,1,2)  
      
      ## make repeated rows unique
      for (we in sort(unique(WW$WW[,1]))){
         hilf <- unique(WW$WW[which(WW$WW[,1]==we),2])
         WW$WW <- rbind(WW$WW[which(WW$WW[,1]<we),],
                        cbind(rep(we,length(hilf)),hilf),
                        WW$WW[which(WW$WW[,1]>we),])
      }

    if (!is.null(ngroups)) {
       ngroups <- append(ngroups, as.list(setdiff(2:(p+1),c(list2vec(ngroups),numalw))))
      ## groups is made longer elsewhere
       ngroups <- lapply(ngroups, function(obj){
           if (length(obj)>1) {
             if (any ((obj-1) %in% WW$WW)) 
               stop ("Groups must not contain effects involved in interactions.")
             else obj <- obj[1] ## replace multi-effect group with one of the effects
                                ## for enabling correct count in factorWW
             }
             obj
             })
     }

# remove blank positions between effect numbers in case of always
       if (length(numalw)>0 && !is.null(ngroups))
       ngroups <- lapply(ngroups, function(obj){
           obj - rowSums(matrix(obj,length(obj),length(numalw),byrow=F) >  
                  matrix(numalw,length(obj),length(numalw),byrow=T))})

    # make 0 length ngroups NULL
    if (!is.null(ngroups)) 
    {if (length(ngroups)==0) ngroups <- NULL}


    bt <- do.call("boot.relimp.default.intern", list(DATA, weights=weights, 
             groups = groups, groupnames = groupnames, WW = WW, always = always, ngroups=ngroups, ..., b=b))
    }
    else {
    ## Berechnung ohne WW
    
    ## selection of columns from model needed because of e.g. lm(y~x1+x2+x3-x2)
    ## selection of columns from model based on formula below 
    ##               works even in case of multi-column terms such as poly(x2,3)

    ## as.data.frame sure that things work correctly even for multi-column effects like polynomials
    ## however, these cannot be grouped as far as I see (would have to be done manually
 
    DATA <- as.data.frame(lm$model[,c(resp,which(rowSums(attr(terms,"factors"))>0))])

    ## change UG 1.3: added default to call, in order to make call output work
    ## added weights to call
    bt <- do.call("boot.relimp.default", list(DATA, weights=weights, type=type, groups=groups, 
          groupnames = groupnames, always=always, ..., b=b))
    }
    bt
}
