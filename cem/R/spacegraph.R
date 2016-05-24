`spacegraph` <-
function (treatment = NULL, data = NULL, R = list(cem = 50, psm = 0, 
    mdm = 0, matchit = 0), grouping = NULL, drop = NULL, L1.breaks = NULL, 
    L1.grouping = NULL, fixed = NULL, minimal = 1, maximal = 15, 
    M = 100, raw.profile = NULL, keep.weights = FALSE, progress = TRUE, 
    rgrouping = FALSE, groups = NULL, psmpoly = 1, mdmpoly = 1, 
    other.matches = NULL, heuristic = FALSE, linear.pscore = FALSE) 
{
    if (nrow(na.omit(data)) != nrow(data)) 
        print("Missing values may lead to errors!")
    if (!treatment %in% names(data)) 
        stop("Treatment is not in the specified data.  Check spelling.")
    if (sum(drop %in% names(data)) != length(drop)) 
        stop("Variables to drop do not appear in the data. Check spelling.")
    if (length(names(data)[!names(data) %in% c(drop, treatment)]) <= 
        3) {
        if ("psm" %in% names(R) & (R$psm > 0)) {
            stop("Must match with 4 or more variables for PSM, MDM, or MatchIt.")
        }
        if ("mdm" %in% names(R) & (R$mdm > 0)) {
            stop("Must match with 4 or more variables for PSM, MDM, or MatchIt.")
        }
        if ("matchit" %in% names(R) & (R$matchit > 0)) {
            stop("Must match with 4 or more variables for PSM, MDM, or MatchIt.")
        }
    }
    if (!is.null(raw.profile) & class(raw.profile) != "L1profile") 
        stop("raw.profile must be of class `L1profile'")
    if (!is.list(R)) 
        stop("R must be supplied as a list.  Ex: R = list(cem=100)")
    if (sum(names(R) %in% c("cem", "psm", "mdm", "matchit") == 
        F) > 0) 
        stop("names for R must be cem, mdm, psm, or matchit. Ex: R = list(cem=5, mdm=5, psm=5)")
    balance.metric <- "all"
    gn <- NULL
    if (!is.null(grouping) & !is.null(names(grouping))) {
        gn <- names(grouping)
        n.gn <- length(gn)
        for (g in 1:n.gn) {
            if (!is.null(data)) 
                data[[gn[g]]] <- group.var(data[[gn[g]]], grouping[[g]])
        }
    }
    if (!is.null(other.matches)) {
        if (!is.list(other.matches)) {
            stop("other.matches must be a list of data frames")
        }
        for (i in 1:length(other.matches)) {
            mm <- other.matches[[i]]
            if (!is.data.frame(mm)) {
                stop("other.matches must be a list of data frames")
            }
            if (sum(names(mm) == c("id", "weight", "method")) != 
                3) {
                stop("other.matches must be a list of dataframes, with names: id, weigth, method.  See the help file for an example.")
            }
        }
    }
    if (sum(names(table(data[[treatment]])) == as.character(c(0, 
        1))) != 2) {
        stop("Treatment must be 0, 1 for now: issue is rcaliper draw in psmspace and mdmspace")
    }
    rnames <- c("cem", "psm", "mdm", "matchit")
    for (i in rnames) {
        if (is.null(R[[i]])) {
            R[[i]] <- 0
        }
    }
    original.data <- data
    imb0 <- NULL
    medianL1 <- NULL
    medianCP <- NULL
    drop <- unique(drop)
    dropped <- match(drop, colnames(data))
    dropped <- dropped[!is.na(dropped)]
    if (length(dropped) > 0) 
        data <- data[-dropped]
    vnames <- colnames(data)
    if (!is.null(treatment)) {
        groups <- as.factor(data[[treatment]])
        idx <- match(treatment, colnames(data))
        if (length(idx) > 0) 
            vnames <- vnames[-idx]
    }
    if (is.null(raw.profile)) {
        if (is.null(L1.breaks)) {
            cat("\nCalculating L1 profile for the raw data...\n")
            imb0 <- L1.profile(groups, data, drop = treatment, 
                M = M, plot = FALSE, progress = progress)
            medianL1 <- median(imb0$L1)
            medianCP <- imb0$CP[[which(imb0$L1 >= medianL1)[1]]]
            medianGR <- imb0$GR[[which(imb0$L1 >= medianL1)[1]]]
        }
        else {
            medianL1 <- L1.meas(groups, data, drop = treatment, 
                breaks = L1.breaks, grouping = L1.grouping)$L1
            medianCP <- L1.breaks
            medianGR <- L1.grouping
        }
    }
    else {
        imb0 <- raw.profile
        medianL1 <- raw.profile$medianL1
        medianCP <- raw.profile$medianCP
        medianGR <- raw.profile$medianGR
    }
    alldist <- NULL
    if (balance.metric == "mdisc" | balance.metric == "all" & R$mdm>0) {
        ff <- as.formula(paste(treatment, "~", paste(vnames, 
            collapse = "+")))
        try(alldist <- mdistfun(ff, data), silent = T)
        if (exists("alldist") == F) {
            alldist <- mdistfun2(ff, data)
        }
    }
    mnames <- vnames
    if (!is.null(gn)) {
        idx <- match(gn, vnames)
        if (length(idx) > 0) 
            mnames <- mnames[-idx]
    }
    if (!is.null(fixed)) {
        idx <- match(fixed, vnames)
        if (length(idx) > 0) 
            mnames <- mnames[-idx]
    }
    nv <- length(mnames)
    v.num <- 1:nv
    b.seq <- vector(nv, mode = "list")
    names(b.seq) <- mnames
    tmp.min <- 2
    tmp.max <- 7
    if (!is.list(minimal)) {
        tmp.min <- minimal + 1
        minimal <- vector(nv, mode = "list")
        for (i in v.num) minimal[[i]] <- tmp.min
    }
    if (!is.list(maximal)) {
        tmp.max <- maximal + 1
        maximal <- vector(nv, mode = "list")
        for (i in v.num) maximal[[i]] <- tmp.max
    }
    for (i in v.num) {
        vna <- mnames[i]
        min.br <- tmp.min
        max.br <- tmp.max
        nuval <- length(unique(data[[vna]]))
        if ((nuval == 2) | !(is.numeric(data[[vna]]) | is.integer(data[[vna]]) | 
            is.logical(data[[vna]]))) 
            max.br <- nuval + 1
        if (!is.null(minimal[[vna]])) 
            min.br <- minimal[[vna]] + 1
        min.br <- max(tmp.min, min.br)
        if (!is.null(maximal[[vna]])) 
            max.br <- maximal[[vna]] + 1
        max.br <- min(tmp.max, max.br)
        b.seq[[i]] <- min.br:max.br
    }
    relax <- NULL
    g.names <- levels(groups)
    n.groups <- length(g.names)
    tab <- as.data.frame(matrix(NA, R$cem + 1, 2 * n.groups + 
        6))
    colnames(tab) <- c(paste("G", g.names, sep = ""), paste("PercG", 
        g.names, sep = ""), "ML1", "Mdiff", "Mdisc", "Relaxed", 
        "Method", "Call")
    n.coars <- dim(tab)[1]
    coars <- vector(n.coars, mode = "list")
    weights <- NULL
    if (keep.weights) 
        weights <- vector(n.coars, mode = "list")
    tab[1, 1:n.groups] <- as.numeric(table(groups))
    tab[1, (n.groups + 1):(2 * n.groups)] <- 100
    tab$Relaxed[1] <- "<raw>"
    coars[[1]] <- NULL
    tab[1, "ML1"] <- medianL1
    tab[1, "Mdiff"] <- mdiff(caliperdat = data, alldat = data, 
        mvars = names(data)[-which(names(data) %in% treatment)], 
        tvar = treatment, wt = NULL)
    if(!is.null(alldist))
    tab[1, "Mdisc"] <- mdisc(caliperdat = data, adist = alldist, 
        wt = NULL)
    else 
     tab[1, "Mdisc"] <- 0    
    tab[1, "Method"] <- "raw"
    tab[1, "Call"] <- "none"
    newcut <- vector(nv, mode = "list")
    names(newcut) <- mnames
    if (R$cem > 0) {
        cat(sprintf("Executing %d different random CEM solutions\n", 
            R$cem))
        if (progress == T & R$cem > 1) {
            pb <- txtProgressBar(min = 1, max = R$cem, initial = 1, 
                style = 3)
        }
        for (r in 1:R$cem) {
            if (progress == T & R$cem > 1) {
                setTxtProgressBar(pb, r)
            }
            for (i in 1:nv) newcut[[i]] <- sample(b.seq[[i]], 
                1)
            if (!rgrouping) 
                obj <- cem(treatment, data, cutpoints = newcut, eval.imbalance = FALSE)
            if (rgrouping) {
                grps <- do.random.groups2(data, groups, tvar = treatment)$grouping
                obj <- cem(treatment, data, cutpoints = newcut, eval.imbalance = FALSE,
                  grouping = grps)
            }
            cpholder <- rep(NA, length(obj$breaks))
            for (ii in 1:length(obj$breaks)) {
                cpholder[ii] <- paste(names(obj$breaks)[ii], 
                  "=c(", paste(round(obj$breaks[[ii]], 2), collapse = ", "), 
                  ")", sep = "")
            }
            breakslist <- paste("list(", paste(cpholder, collapse = ", "), 
                ")")
            if (!rgrouping) {
                grouplist <- ""
            }
            if (rgrouping) {
                gpholder <- rep(NA, length(grps))
                for (ii in 1:length(grps)) {
                  ttt <- paste(names(grps)[ii], "=list(\"", paste(grps[[ii]], 
                    collapse = "\", \""), "\")", sep = "")
                  ttt <- gsub("\"list(", "c(", ttt, fixed = T)
                  ttt <- gsub("\"c(", "c(", ttt, fixed = T)
                  ttt <- gsub("\")\",", "\"),", ttt, fixed = T)
                  ttt <- gsub("\")\")", "\"))", ttt, fixed = T)
                  gpholder[ii] <- ttt
                }
                grouplist <- paste(", grouping = list(", paste(gpholder, 
                  collapse = ", "), ")")
            }
            tab[r + 1, "Call"] <- sprintf("cem(obj$match$treatment, data=obj$data[,c(obj$match$vars,obj$match$treatment) ], cutpoints=%s, eval.imbalance = FALSE, L1.breaks=obj$raw.profile$medianCP %s)",
                breakslist, grouplist)
            coars[[r + 1]] <- obj$breaks
            if (keep.weights) 
                weights[[r + 1]] <- obj$w
            tab[r + 1, 1:(2 * n.groups)] <- as.numeric(c(obj$tab[2, 
                ], obj$tab[2, ]/obj$tab[1, ] * 100))
            tab$Relaxed[r + 1] <- "random"
            if (balance.metric == "L1" | balance.metric == "all") {
                if (max(obj$tab[2, ] == 0) != 1) {
                  tab[r + 1, "ML1"] <- L1.meas(groups, data, 
                    drop = treatment, breaks = medianCP, weights = obj$w, 
                    grouping = medianGR)$L1
                }
            }
            if (balance.metric == "mdiff" | balance.metric == 
                "all") {
                if (max(obj$tab[2, ] == 0) != 1) {
                  tab[r + 1, "Mdiff"] <- mdiff(caliperdat = data, 
                    alldat = data, mvars = obj$vars, tvar = treatment, 
                    wt = obj$w)
                }
            }
            if (balance.metric == "mdisc" | balance.metric == 
                "all") {
                if (max(obj$tab[2, ] == 0) != 1) {
                  tab[r + 1, "Mdisc"] <- mdisc(caliperdat = data, 
                    adist = alldist, wt = obj$w)
                }
            }
            tab[r + 1, "Method"] <- "cem"
        }
        if (progress == T & R$cem > 1) {
            close(pb)
        }
    }
    if (R$cem == 0) {
        obj <- cem(treatment, data, cutpoints = newcut, eval.imbalance = FALSE)
    }
    if (R$psm > 0) {
        if (length(g.names) > 2) {
            stop("PSM requires a dichotomous treatment")
        }
        cat(sprintf("Executing %d different random PSM solutions\n", 
            R$psm))
        myps <- psmspace(treatment = treatment, data = data, 
            R = R$psm, poly = psmpoly, randomgroups = rgrouping, 
            groups = groups, rawCP = medianCP, progr = progress, 
            heur = heuristic, balance.metric = balance.metric, 
            linear.pscore = linear.pscore, alldist = alldist)
        tab <- rbind(tab, myps$tab)
    }
    if (R$mdm > 0) {
        if (length(g.names) > 2) {
            stop("MDM requires a dichotomous treatment")
        }
        cat(sprintf("Executing %d different random MDM solutions\n", 
            R$mdm))
        mymdm <- mdmspace(treatment = treatment, data = data, 
            R = R$mdm, poly = mdmpoly, randomgroups = rgrouping, 
            groups = groups, rawCP = medianCP, progr = progress, 
            heur = heuristic, balance.metric = balance.metric, 
            alldist = alldist)
        tab <- rbind(tab, mymdm$tab)
    }
    if (R$matchit > 0) {
        cat(sprintf("Executing %d different random MatchIt solutions\n", 
            R$matchit))
        mymatchit <- matchitspace(treatment = treatment, data = data, 
            R = R$matchit, poly = psmpoly, randomgroups = rgrouping, 
            groups = groups, rawCP = medianCP, progr = progress, 
            heur = heuristic)
        tab <- rbind(tab, mymatchit$tab)
    }
    if (!is.null(other.matches)) {
        cat(sprintf("Calculating L1 for %d solutions provided by the user\n", 
            length(other.matches)))
        rtab <- as.data.frame(matrix(NA, length(other.matches), 
            2 * n.groups + 6))
        colnames(rtab) <- c(paste("G", g.names, sep = ""), paste("PercG", 
            g.names, sep = ""), "ML1", "Mdiff", "Mdisc", "Relaxed", 
            "Method", "Call")
        if (progress == T) {
            pb <- txtProgressBar(min = 0, max = length(other.matches), 
                initial = 0, style = 3)
        }
        for (i in 1:length(other.matches)) {
            if (progress == T) {
                setTxtProgressBar(pb, i)
            }
            mm <- other.matches[[i]]
            rdata <- data[as.character(mm$id), ]
            if (nrow(na.omit(rdata)) != nrow(rdata)) {
                stop("other.matches specifies observations that aren't in the dataset.")
            }
            if (balance.metric == "L1" | balance.metric == "all") {
                rtab[i, "ML1"] <- (my.imbalance(rdata[[treatment]], 
                  rdata, drop = c(treatment), weights = mm$weight, 
                  breaks = medianCP))$L1$L1
            }
            if (balance.metric == "mdiff" | balance.metric == 
                "all") {
                rtab[i, "Mdiff"] <- mdiff(caliperdat = rdata, 
                  alldat = data, mvars = names(rdata)[-which(names(rdata) %in% 
                    treatment)], tvar = treatment, wt = mm$weight)
            }
            if (balance.metric == "mdisc" | balance.metric == 
                "all") {
                rtab[i, "Mdisc"] <- mdisc(caliperdat = rdata, 
                  adist = alldist, wt = mm$weight)
            }
            for (j in 1:length(g.names)) {
                rtab[i, paste("G", g.names[j], sep = "")] <- sum(rdata[[treatment]] == 
                  g.names[j] & mm$weight > 0)
                rtab[i, paste("PercG", g.names[j], sep = "")] <- sum(rdata[[treatment]] == 
                  g.names[j] & mm$weight > 0)/sum(rdata[[treatment]] == 
                  g.names[j])
                if (length(unique(as.character(mm$method))) > 
                  1) {
                  warning("More than one method specified.  Using the first.")
                }
            }
            rtab[i, "Relaxed"] <- "user"
            rtab[i, "Method"] <- as.character(mm$method)[1]
            rtab[i, "Call"] <- as.character(mm$method)[1]
        }
        if (progress == T) {
            close(pb)
        }
        tab <- rbind(tab, rtab)
    }
    if (balance.metric == "L1" | balance.metric == "all") {
        idx <- order(tab[, "ML1"])
        tab <- tab[idx, ]
    }
    if (balance.metric == "mdiff") {
        idx <- order(tab[, "Mdiff"])
        tab <- tab[idx, ]
    }
    if (balance.metric == "mdisc") {
        idx <- order(tab[, "Mdisc"])
        tab <- tab[idx, ]
    }
    rownames(tab) <- 1:(dim(tab)[1])
    out <- list(space = tab)
    out$L1breaks <- medianCP
    out$raw.profile <- imb0
    out$tab <- obj$tab
    out$medianCP <- medianCP
    out$medianL1 <- medianL1
    out$coars <- coars[idx]
    if (keep.weights) 
        out$weights <- weights[idx]
    out$n.coars <- n.coars
    tmp.list <- c()
    for (i in 1:length(vnames)) {
        tmp.list[[i]] <- c(0, 0)
    }
    names(tmp.list) <- vnames
    for (i in 1:length(out$coars)) {
        if (is.null(out$coars[[i]])) {
            out$coars[[i]] <- tmp.list
        }
    }
    out$match <- obj
    out$data <- data
    out$alldist <- alldist
    class(out) <- "spacegraph"
    if (exists("myps")) {
        if (nrow(na.omit(myps$tab)) < nrow(myps$tab)) {
            print(paste("WARNING: ", nrow(myps$tab) - nrow(na.omit(myps$tab)), 
                " of ", nrow(myps$tab), " propensity score models did not converge", 
                sep = ""))
        }
    }
    return(invisible(out))
}



`spacegraphold` <-
function (treatment=NULL, data = NULL, R=list("cem"=50,"psm"=0,"mdm"=0,"matchit"=0),
    grouping = NULL, drop=NULL,
    L1.breaks = NULL, L1.grouping=NULL, fixed = NULL, 
    minimal = 1, maximal = 15, M=100, 
    raw.profile=NULL, keep.weights=FALSE, progress=TRUE,
    rgrouping=FALSE, groups=NULL, psmpoly=1, mdmpoly=1,
    other.matches=NULL, heuristic=FALSE, linear.pscore=FALSE) 
{

    
    if (nrow(na.omit(data)) != nrow(data))
       print("Missing values may lead to errors!")

    #if (nrow(na.omit(data)) != nrow(data) & R$matchit > 0) ## chokes if no R$matchit is specified
    #   stop("Cannot use MatchIt with missing data.")

    if (!treatment %in% names(data))
       stop("Treatment is not in the specified data.  Check spelling.")

    if (sum(drop %in% names(data)) != length(drop))
       stop("Variables to drop do not appear in the data. Check spelling.")

    if (length(names(data)[!names(data) %in% c(drop,treatment)]) <=3){
      if("psm" %in% names(R) & (R$psm>0)){
       stop("Must match with 4 or more variables for PSM, MDM, or MatchIt.")
      }
      if("mdm" %in% names(R) & (R$mdm>0)){
       stop("Must match with 4 or more variables for PSM, MDM, or MatchIt.")
      }
      if("matchit" %in% names(R) & (R$matchit>0)){
       stop("Must match with 4 or more variables for PSM, MDM, or MatchIt.")
      }
    }
       
    if (!is.null(raw.profile) & class(raw.profile) != "L1profile") 
	 stop("raw.profile must be of class `L1profile'")

    if (!is.list(R))
         stop("R must be supplied as a list.  Ex: R = list(cem=100)")

    if (sum(names(R) %in% c("cem","psm","mdm","matchit") == F) > 0)
         stop("names for R must be cem, mdm, psm, or matchit. Ex: R = list(cem=5, mdm=5, psm=5)")

    ## As of 13 July 2011, I hardcode the balance metric so that it
    ## always does both.  This could be undone and it could be listed
    ## as an argument taking "L1", "mdiff", "mdisc", or "all", but the problem
    ## was that I then couldn't do na.omit on the table to do the
    ## plotting.  I have to do some na.omitting because some of the
    ## cems have no obs, and thus no L1.
    balance.metric <- "all"

    gn <- NULL
    if (!is.null(grouping) & !is.null(names(grouping))) {
        gn <- names(grouping)
        n.gn <- length(gn)
        for (g in 1:n.gn) {
            if (!is.null(data)) 
			data[[gn[g]]] <- group.var(data[[gn[g]]], grouping[[g]])
        }
    }

    if(!is.null(other.matches)){
      if(!is.list(other.matches)){stop("other.matches must be a list of data frames")}
      for(i in 1:length(other.matches)){
        mm <- other.matches[[i]]
        if(!is.data.frame(mm)){stop("other.matches must be a list of data frames")}
        if(sum(names(mm) ==
         c("id","weight","method")) != 3){stop("other.matches must be a list of dataframes, with names: id, weigth, method.  See the help file for an example.")}
      }
    }

    if(sum(names(table(data[[treatment]])) == as.character(c(0,1))) != 2){
      stop("Treatment must be 0, 1 for now: issue is rcaliper draw in psmspace and mdmspace")
    }


        ## fill in the missing elements of R if not specified
        rnames <- c("cem","psm","mdm","matchit")
        for(i in rnames){
          if(is.null(R[[i]])){R[[i]] <- 0}
        }       

        original.data <- data
	imb0 <- NULL
	medianL1 <- NULL
	medianCP <- NULL	

	drop <- unique(drop)
	dropped <- match(drop, colnames(data))
	dropped <- dropped[!is.na(dropped)]
	
	if(length(dropped)>0) 
	 data <- data[-dropped]
	vnames <- colnames(data)

	if(!is.null(treatment)){
		groups <- as.factor(data[[treatment]])
        idx <- match(treatment, colnames(data))
		if(length(idx)>0)
		 vnames <- vnames[-idx]
	}
	
	if(is.null(raw.profile)){
		if(is.null(L1.breaks)){
		 cat("\nCalculating L1 profile for the raw data...\n")
		 imb0 <- L1.profile(groups, data, drop=treatment, M=M,
                                    plot=FALSE, progress=progress)
		 medianL1 <- median(imb0$L1)
		 medianCP <- imb0$CP[[ which(imb0$L1 >= medianL1)[1] ]]
		 medianGR <- imb0$GR[[ which(imb0$L1 >= medianL1)[1] ]]
		} else {
		 medianL1 <- L1.meas(groups, data, drop=treatment, breaks=L1.breaks, grouping=L1.grouping)$L1
		 medianCP <- L1.breaks
		 medianGR <- L1.grouping	
		}
	} else {
		imb0 <- raw.profile
		medianL1 <- raw.profile$medianL1
		medianCP <- raw.profile$medianCP
		medianGR <- raw.profile$medianGR
	}


      #######################################################
      ## Inserting Mahalanobis matching discrepancy here
      if(balance.metric == "mdisc" | balance.metric=="all"){
        ff <- as.formula(paste(treatment,"~",paste(vnames,collapse="+")))
        try(alldist <- mdistfun(ff,data), silent=T)
        if(exists("alldist")==F){
          alldist <- mdistfun2(ff,data)
        }  
      }
    
      #######################################################

        
	mnames <- vnames
	if(!is.null(gn)){
	 idx <- match(gn, vnames)
	 if(length(idx)>0)
	  mnames <- mnames[-idx] 
	}

	if (!is.null(fixed)) {
        idx <- match(fixed, vnames)
        if(length(idx) > 0) 
		 mnames <- mnames[-idx]
    }
	
    nv <- length(mnames)
    v.num <- 1:nv


    b.seq <- vector(nv, mode = "list")
    names(b.seq) <- mnames
	tmp.min <- 2
	tmp.max <- 7
	if(!is.list(minimal)){
	 tmp.min <- minimal+1
	 minimal <- vector(nv, mode="list")
	 for (i in v.num) 
	  minimal[[i]] <- tmp.min
	}
	if(!is.list(maximal)){
		tmp.max <- maximal+1
		maximal <- vector(nv, mode="list")
		for (i in v.num) 
		maximal[[i]] <- tmp.max
	}
    for (i in v.num) {
		vna <- mnames[i]
		min.br <- tmp.min
		max.br <- tmp.max
		
		nuval <- length(unique(data[[vna]]))

        if( (nuval==2) | !(is.numeric(data[[vna]]) | is.integer(data[[vna]]) | is.logical(data[[vna]])) )
			max.br <- nuval+1
		
		if (!is.null(minimal[[vna]]))  
			min.br <- minimal[[vna]] + 1
		min.br <- max(tmp.min, min.br)
		if (!is.null(maximal[[vna]])) 
			max.br <- maximal[[vna]] + 1
		max.br <- min(tmp.max, max.br)
		b.seq[[i]] <- min.br:max.br
    }
    relax <- NULL

	
	g.names <- levels(groups)
	n.groups <- length(g.names)


    tab <- as.data.frame(matrix(NA, R$cem + 1, 2 * n.groups + 6))

    colnames(tab) <- c(paste("G", g.names, sep = ""), 
    					   paste("PercG", g.names, sep
                                                 = ""), "ML1","Mdiff",
                                           "Mdisc",
                                           "Relaxed", "Method","Call")

	n.coars <- dim(tab)[1]
	coars <- vector(n.coars, mode="list")
	weights <- NULL
	if(keep.weights)
	 weights <- vector(n.coars, mode="list")

	tab[1, 1:n.groups] <- as.numeric(table(groups))
    tab[1, (n.groups + 1):(2 * n.groups)] <- 100
    tab$Relaxed[1] <- "<raw>"
	coars[[1]] <- NULL
	tab[1, "ML1"] <- medianL1
      tab[1, "Mdiff"] <- mdiff(caliperdat=data,alldat=data,mvars=names(data)[-which(names(data) %in% treatment)],tvar=treatment, wt=NULL)
      tab[1, "Mdisc"] <- mdisc(caliperdat=data,adist=alldist, wt=NULL)
        tab[1, "Method"] <- "raw"
        tab[1, "Call"] <- "none"

    
	newcut <- vector(nv, mode="list")
    names(newcut) <- mnames

    if(R$cem>0){
      cat(sprintf("Executing %d different random CEM solutions\n", R$cem))

      if(progress==T & R$cem>1){pb <- txtProgressBar(min = 1, max = R$cem, initial = 1, style = 3)}

        for (r in 1:R$cem) {
		if(progress==T & R$cem>1){setTxtProgressBar(pb, r)}
		for(i in 1:nv)
			newcut[[i]] <- sample(b.seq[[i]], 1) 

            if(!rgrouping)
		    obj <- cem(treatment, data, cutpoints=newcut, eval.imbalance=FALSE)

            if(rgrouping){
                grps <- do.random.groups2(data, groups, tvar=treatment)$grouping
                obj <- cem(treatment, data, 
                           cutpoints=newcut, eval.imbalance=FALSE, grouping = grps)
            } 

                ## reconstruct the call
                cpholder <- rep(NA,length(obj$breaks))
                for(ii in 1:length(obj$breaks)){
                  cpholder[ii] <-paste(names(obj$breaks)[ii],"=c(",paste(round(obj$breaks[[ii]],2),collapse=", "),")",sep="")
                }
                breakslist <-  paste("list(",paste(cpholder,collapse=", "),")")
                ## get the groups if there are any
                if(!rgrouping){grouplist <- ""}
                if(rgrouping){
                  gpholder <- rep(NA,length(grps))
                  for(ii in 1:length(grps)){

                    ttt <- paste(names(grps)[ii],"=list(\"",paste(grps[[ii]],collapse="\", \""),"\")",sep="")
                    ttt <- gsub("\"list(","c(",ttt, fixed=T)
                    ttt <- gsub("\"c(","c(",ttt, fixed=T)
                    ttt <- gsub("\")\",", "\"),",ttt, fixed=T)
                    ttt <- gsub("\")\")", "\"))",ttt, fixed=T)
                    gpholder[ii] <- ttt   
                  }
                  grouplist <- paste(", grouping = list(",paste(gpholder,collapse=", "),")")  
                }
                ## note that this "obj" is the "obj" that "spacegraph.plot"
                ## takes as an input, NOT the "obj" above
                tab[r+1,"Call"] <- sprintf("cem(obj$match$treatment, data=obj$data[,c(obj$match$vars,obj$match$treatment) ], cutpoints=%s, eval.imbalance = FALSE, L1.breaks=obj$raw.profile$medianCP %s)", breakslist, grouplist)
                
                
		
		coars[[r+1]] <- obj$breaks
		if(keep.weights)
		 weights[[r+1]] <- obj$w
		
		tab[r+1, 1:(2 * n.groups)] <- as.numeric(c(obj$tab[2,],  obj$tab[2, ]/obj$tab[1, ] * 100))
		tab$Relaxed[r+1] <- "random"
            ## The L1 calculation seems to choke if there are no matches
            if(balance.metric=="L1" | balance.metric=="all"){
              if(max(obj$tab[2,]==0)!=1){
		    tab[r+1, "ML1"] <- L1.meas(groups, data, drop=treatment, breaks=medianCP, 
								    weights=obj$w, grouping=medianGR)$L1
              }
            }
            if(balance.metric=="mdiff" | balance.metric=="all"){
              if(max(obj$tab[2,]==0)!=1){
		    tab[r+1, "Mdiff"] <- mdiff(caliperdat=data,alldat=data,mvars=obj$vars,tvar=treatment, wt=obj$w)

              }
            }
            if(balance.metric=="mdisc" | balance.metric=="all"){
              if(max(obj$tab[2,]==0)!=1){
		    tab[r+1, "Mdisc"] <- mdisc(caliperdat=data,adist=alldist, wt=obj$w)

              }
            }
              
            tab[r+1, "Method"] <- "cem"
  
        }
      if(progress==T & R$cem>1){close(pb)}
    }  # end if R$cem > 0
    ## Do one cem even if R$cem==0 because you need it to make the
    ## output
    if(R$cem==0){
      obj <- cem(treatment, data, cutpoints=newcut, eval.imbalance=FALSE)
    }


    ## Do PSM here
    if(R$psm > 0){
      if(length(g.names) > 2){stop("PSM requires a dichotomous treatment")}
      cat(sprintf("Executing %d different random PSM solutions\n", R$psm))

      myps <- psmspace(treatment=treatment,data=data,R=R$psm,
                     poly=psmpoly,randomgroups=rgrouping,
                     groups=groups,rawCP=medianCP,progr=progress,
                     heur=heuristic, balance.metric=balance.metric,
                     linear.pscore=linear.pscore,alldist=alldist)
      tab <- rbind(tab,myps$tab)
    }


    ## Do MDM here
    if(R$mdm > 0){
      if(length(g.names) > 2){stop("MDM requires a dichotomous treatment")}
      cat(sprintf("Executing %d different random MDM solutions\n", R$mdm))

      mymdm <- mdmspace(treatment=treatment,data=data,R=R$mdm,
                     poly=mdmpoly,randomgroups=rgrouping,
                     groups=groups,rawCP=medianCP,progr=progress,
                     heur=heuristic, balance.metric=balance.metric,
                     alldist=alldist)
      tab <- rbind(tab,mymdm$tab)
    }



  #I PUT THIS MATCHIT THING IN HERE, BUT IT'S NOT PLOTTING FOR SOME REASON.
   
  ## NOTE, "matchitspace" has not been updated to allow mdiff balance.
    ## Run MatchIt here
    if(R$matchit >0){
      cat(sprintf("Executing %d different random MatchIt solutions\n", R$matchit))
      mymatchit <- matchitspace(treatment=treatment,data=data,R=R$matchit,
                     poly=psmpoly,randomgroups=rgrouping,
                     groups=groups,rawCP=medianCP,progr=progress,
                     heur=heuristic)
      tab <- rbind(tab,mymatchit$tab)
    }


    ## Matched Samples provided by User
    ## I think this will allow multinomial treatments.
    if(!is.null(other.matches)){
      cat(sprintf("Calculating L1 for %d solutions provided by the user\n", length(other.matches)))
      rtab <- as.data.frame(matrix(NA, length(other.matches), 2 *
                                   n.groups + 6))
      #if(sum(g.names == c("0","1")) != 2){stop("other.matches requires a dichotomous treatment")}
      colnames(rtab) <- c(paste("G", g.names, sep = ""), 
    					   paste("PercG", g.names, sep
                                                 = ""), "ML1","Mdiff",
                                           "Mdisc",
                                           "Relaxed", "Method","Call")
      if(progress==T){pb <- txtProgressBar(min = 0, max = length(other.matches), initial = 0, style = 3)}
      for(i in 1:length(other.matches)){
        if(progress==T){setTxtProgressBar(pb, i)}
        mm <- other.matches[[i]]
        ## order the dataset according to the weights...
        rdata <- data[as.character(mm$id),]
        if(nrow(na.omit(rdata)) != nrow(rdata)){stop("other.matches specifies observations that aren't in the dataset.")}
        if(balance.metric=="L1" | balance.metric=="all"){
          rtab[i,"ML1"] <- (my.imbalance(rdata[[treatment]],rdata,drop=c(treatment),
                          weights=mm$weight,breaks=medianCP))$L1$L1
        }
        if(balance.metric=="mdiff" | balance.metric=="all"){
          rtab[i,"Mdiff"] <- mdiff(caliperdat=rdata,alldat=data,mvars=names(rdata)[-which(names(rdata) %in% treatment)] ,tvar=treatment, wt=mm$weight)
        }
        if(balance.metric=="mdisc" | balance.metric=="all"){
          rtab[i,"Mdisc"] <- mdisc(caliperdat=rdata,adist=alldist, wt=mm$weight)
        }
        ## loop through the levels of treatment
        for(j in 1:length(g.names)){
          rtab[i,paste("G",g.names[j],sep="")] <- sum(rdata[[treatment]]==g.names[j] & mm$weight > 0)
          rtab[i,paste("PercG",g.names[j],sep="")] <- sum(rdata[[treatment]]==g.names[j] & mm$weight > 0)/sum(rdata[[treatment]]==g.names[j])
          if(length(unique(as.character(mm$method))) > 1){warning("More than one method specified.  Using the first.")}
        }
        rtab[i,"Relaxed"] <- "user"
        rtab[i,"Method"] <- as.character(mm$method)[1]
        rtab[i,"Call"] <- as.character(mm$method)[1]
      }
      if(progress==T){close(pb)}
      tab <- rbind(tab, rtab)

    }



    
    ## Finalize output
    if(balance.metric=="L1" | balance.metric=="all"){
      idx <- order(tab[, "ML1"])
      tab <- tab[idx, ]
    }
    if(balance.metric=="mdiff"){
      idx <- order(tab[, "Mdiff"])
      tab <- tab[idx, ]
    }
    if(balance.metric=="mdisc"){
      idx <- order(tab[, "Mdisc"])
      tab <- tab[idx, ]
    }
    rownames(tab) <- 1:(dim(tab)[1])
	out <- list(space = tab)
	out$L1breaks <- medianCP
	out$raw.profile <- imb0
	out$tab <- obj$tab  ## this is only the table for the last cem anyway.  And it makes it crash if there are no cems.
    out$medianCP <- medianCP
	out$medianL1 <- medianL1
	out$coars <- coars[idx]
	if(keep.weights)
	 out$weights <- weights[idx]
	out$n.coars <- n.coars
        ## make a filler for now
        tmp.list <- c()
        for(i in 1:length(vnames)){
          tmp.list[[i]] <- c(0,0)
        }
        names(tmp.list) <- vnames
        for(i in 1:length(out$coars)){
          if(is.null(out$coars[[i]])){
            out$coars[[i]] <- tmp.list
          }
        }
	out$match <- obj  ## again, this would just be the latest cem
                        ## but I use it for running CEM from the GUI.
        #out$data <- original.data
        out$data <- data
      ## add the mahalanobis discrepancy so we can add new matches to the plot
      out$alldist <- alldist
    class(out) <- "spacegraph"
	
    ## give a warning about nonconverged pscores
    if(exists("myps")){
      if(nrow(na.omit(myps$tab)) < nrow(myps$tab)){
        print(paste("WARNING: ",nrow(myps$tab)-nrow(na.omit(myps$tab))," of ",nrow(myps$tab),
            " propensity score models did not converge",sep=""))
    }
    }
    #if (plot) 
    #  #plot(out,data=data,explore=interactive())
    #      spacegraph.plot(out,explore=interactive())
	
    return(invisible(out))
}


###################################################################
## End spacegraph function
###################################################################


  ## DEFINE FUNCTIONS

####################################################################
## Data prep function
##
## This expects the following arguments: "treatment", "data", "drop"

dataPrep <- function(treatment, data, drop=NULL, dropNA = T){
  if(is.null(drop)==F){
    dropme <- match(drop,colnames(data))
    dat <- data[,-dropme]
  }
  if(dropNA){
    dat <- na.omit(dat)
  }
  return(dat)
}
####################################################################


####################################################################
## psmspace
##
## This function expects the arguments: "treatment", "data"



psmspace <- function(treatment, data, R = 100, poly=2, randomgroups=FALSE,
                     groups=NULL,rawCP, progr, 
                     heur, balance.metric="L1",linear.pscore=FALSE,
                     alldist){
 
  start.time <- Sys.time()
  ## make some holders
  outN <- matrix(NA,R,2)
  outL1 <- rep(NA,R)
  outMdiff <- rep(NA,R)
  outMdisc <- rep(NA,R)
  outobs <- c()
  mycalls <- rep(NA,R)
  ## start the loop
  if(progr & R>1){pb <- txtProgressBar(min = 1, max = R, initial = 1, style = 3)}
  for(qq in 1:R){
    mvars <- names(data)[names(data)!=treatment]
    ## generate random groupings for the factors
    ## PROBLEM -- If you do this, it gives you false L1s because
    ## it changes the data!!!
    #if(randomgroups & is.null(groups)){
    #  newdat <- do.random.groups(data,groups)
    #  data <- newdat$dat
    #  grouping <- newdat$grouping
    #}
    ## generate all of the interactions
    allIntObj <- enumerateInteractions(mvars, poly)
    ## generate a formula
    if(poly==1){
      if(!allIntObj$isinf$total1 & !heur){
        myform <- genForm(allIntObj,treatment=treatment, data=data, poly=poly)
      } else {
        myform <- heuristicForm(allIntObj,treatment=treatment, data=data, poly=poly)
      }
    }
    if(poly==2){
      if(!allIntObj$isinf$total2 & !allIntObj$isinf$total1 & !heur){
        myform <- genForm(allIntObj,treatment=treatment, data=data, poly=poly)
      } else {
        myform <- heuristicForm(allIntObj,treatment=treatment, data=data, poly=poly)
      }
    }   
    ## draw the caliper
    rcaliper <- sample(seq(0,(sum(data[[treatment]])-2),1), size=1)
    ## do the pscore matching, now throwing out non-converged ps models
    if(exists("match.obj")){rm(match.obj)}
    try(
    match.obj <- psm(formula = myform, data = data, ratio=1,
                     caliper=rcaliper, linear.pscore=linear.pscore,rawCP=rawCP)
    ,silent = T)
    if(exists("match.obj")==F){
      #print("PS model didn't converge")
      if(progr & R>1){setTxtProgressBar(pb, qq)}
      next
    }
    mycalls[qq] <- paste("psm(formula = ",myform," , data = obj$data, caliper = ",rcaliper," , linear.pscore = ",linear.pscore,",rawCP=rawCP)",sep="")
    ## calculate L1
    if(balance.metric=="L1" | balance.metric=="all"){
      myL1 <- (my.imbalance(match.obj$match.dat[[treatment]],match.obj$match.dat,drop=c(treatment,"distance","weights"),breaks=rawCP))$L1$L1
    } else {
      myL1 <- NA
    }
    if(balance.metric=="mdiff" | balance.metric=="all"){
      myMdiff <- mdiff(caliperdat=match.obj$match.dat,alldat=data,
                       mvars=names(match.obj$match.dat)[-which(names(match.obj$match.dat) %in% c(treatment,"distance","weights"))],
                       tvar=treatment, wt=match.obj$match.dat$weights)
    } else {
      myMdiff <- NA
    }
    if(balance.metric=="mdisc" | balance.metric=="all"){
      myMdisc <- mdisc(caliperdat=match.obj$match.dat,adist=alldist,
                       wt=match.obj$match.dat$weights)
    } else {
      myMdisc <- NA
    }
    ## put things in their holders
    outN[qq,] <- match.obj$N
    outL1[qq] <- myL1
    outMdiff[qq] <- myMdiff
    outMdisc[qq] <- myMdisc
    outobs[[qq]] <- rownames(match.obj$match.dat)
    if(progr & R>1){setTxtProgressBar(pb, qq)}
    #print(paste("PSM",qq,"of",R,": time =",round(Sys.time()-start.time,2)))
  }
  if(progr & R>1){close(pb)}
  ## return everything
  tab <- as.data.frame(matrix(NA, R, 2 * 2 + 6))
  g.names <- c("0","1")
  colnames(tab) <- c(paste("G", g.names, sep = ""), 
    					   paste("PercG", g.names, sep
                                                 = ""), "ML1","Mdiff",
                                           "Mdisc",
                                           "Relaxed", "Method","Call")
  tab[,1:2] <- outN
  tab[,3] <- tab[,1]/table(data[[treatment]])[1]
  tab[,4] <- tab[,2]/table(data[[treatment]])[2]
  tab[,5] <- outL1
  tab[,6] <- outMdiff
  tab[,7] <- outMdisc
  tab[,8] <- rep("random",R)
  tab[,9] <- rep("psm",R)
  tab[,10] <- mycalls

  out <- list(tab=tab,matched=outobs)
  return(out)
}

####################################################################



####################################################################
## mdmspace
##
## This function expects the arguments: "treatment", "data"



mdmspace <- function(treatment, data, R = 100, poly=2, randomgroups=F,
                     groups=NULL,rawCP, progr,
                     heur, balance.metric="L1",
                     alldist){
  start.time <- Sys.time()
  ## make some holders
  outN <- matrix(NA,R,2)
  outL1 <- rep(NA,R)
  outMdiff <- rep(NA,R)
  outMdisc <- rep(NA,R)
  outobs <- c()
  mycalls <- rep(NA,R)
  ## start the loop
  if(progr & R>1){pb <- txtProgressBar(min = 1, max = R, initial = 1, style = 3)}
  for(qq in 1:R){
    mvars <- names(data)[names(data)!=treatment]
    ## generate random groupings for the factors
    ## PROBLEM -- If you do this, it gives you false L1s because
    ## it changes the data!!!
    #if(randomgroups){
    #  newdat <- do.random.groups(data,groups)
    #  data <- newdat$dat
    #  grouping <- newdat$grouping
    #}
    ## generate all of the interactions
    allIntObj <- enumerateInteractions(mvars, poly)
    ## generate a formula
    if(poly==1){
      if(!allIntObj$isinf$total1 & !heur){
        myform <- genForm(allIntObj,treatment=treatment, data=data, poly=poly)
      } else {
        myform <- heuristicForm(allIntObj,treatment=treatment, data=data, poly=poly)
      }
    }
    if(poly==2){
      if(!allIntObj$isinf$total2 & !allIntObj$isinf$total1 & !heur){
        myform <- genForm(allIntObj,treatment=treatment, data=data, poly=poly)
      } else {
        myform <- heuristicForm(allIntObj,treatment=treatment, data=data, poly=poly)
      }
    } 
    myform.tmp <- myform
    myform <- as.formula(myform)
    ## draw the caliper
    rcaliper <- sample(seq(0,(sum(data[[treatment]])-2),1), size=1)
    ## do the mdm matching
    if(exists("match.obj")){
      try(rm(match.obj), silent=T)
    }
    try(match.obj <- mdm(formula = myform, data = data, ratio=1, caliper=rcaliper,rawCP=rawCP), silent=T)
    if(exists("match.obj")==F){
      print(paste("MDM failed on run",qq,": probably too many colinear interactions"))
      next
    }
    mycalls[qq] <- paste("mdm(formula = ",myform.tmp,", data = obj$data, caliper=",rcaliper,",rawCP=rawCP)",sep="")
    ## calculate L1
    if(balance.metric=="L1" | balance.metric=="all"){
      myL1 <- (my.imbalance(match.obj$match.dat[[treatment]],match.obj$match.dat,drop=c(treatment,"distance","weights"),breaks=rawCP))$L1$L1
    } else {
      myL1 <- NA
    }
    if(balance.metric=="mdiff" | balance.metric=="all"){
      myMdiff <- mdiff(caliperdat=match.obj$match.dat,alldat=data,
                       mvars=names(match.obj$match.dat)[-which(names(match.obj$match.dat) %in% c(treatment,"distance","weights"))],
                       tvar=treatment, wt=match.obj$match.dat$weights)
    } else {
      myMdiff <- NA
    }
    if(balance.metric=="mdisc" | balance.metric=="all"){
      myMdisc <- mdisc(caliperdat=match.obj$match.dat, adist=alldist,
                       wt=match.obj$match.dat$weights)
    } else {
      myMdisc <- NA
    }
    ## put things in their holders
    outN[qq,] <- match.obj$N
    outL1[qq] <- myL1
    outMdiff[qq] <- myMdiff
    outMdisc[qq] <- myMdisc
    outobs[[qq]] <- rownames(match.obj$match.dat)
    if(progr & R>1){setTxtProgressBar(pb, qq)}
    #print(paste("MDM",qq,"of",R,": time =",round(Sys.time()-start.time,2)))
  }
  if(progr & R>1){close(pb)}
  ## return everything
    ## return everything
  tab <- as.data.frame(matrix(NA, R, 2 * 2 + 6))
  g.names <- c("0","1")
  colnames(tab) <- c(paste("G", g.names, sep = ""), 
    					   paste("PercG", g.names, sep
                                                 = ""), "ML1","Mdiff",
                                           "Mdisc",
                                           "Relaxed", "Method","Call")
  tab[,1:2] <- outN
  tab[,3] <- tab[,1]/table(data[[treatment]])[1]
  tab[,4] <- tab[,2]/table(data[[treatment]])[2]
  tab[,5] <- outL1
  tab[,6] <- outMdiff
  tab[,7] <- outMdisc
  tab[,8] <- rep("random",R)
  tab[,9] <- rep("mdm",R)
  tab[,10] <- mycalls

  out <- list(tab=tab,matched=outobs)
  return(out)
}

####################################################################


####################################################################
## matchitspace
##
## This function expects the arguments: "treatment", "data"


#matchitspace <- function(treatment, data, R = 100, poly=2, randomgroups=F,
#                     groups=NULL,rawCP=rawCP, progr=progress,
#                     heur=heuristic){
matchitspace <- function(treatment, data, R = 100, poly=2, randomgroups=FALSE,
                     groups=NULL,rawCP, progr,
                     heur){

  start.time <- Sys.time()
  ## make some holders
  outN <- matrix(NA,R,2)
  outL1 <- rep(NA,R)
  outobs <- c()
  mycalls <- rep(NA,R)
  ## start the loop
  if(progr & R>1){pb <- txtProgressBar(min = 1, max = R, initial = 1, style = 3)}
  for(qq in 1:R){
    mvars <- names(data)[names(data)!=treatment]
    ## generate random groupings for the factors
    ## PROBLEM -- If you do this, it gives you false L1s because
    ## it changes the data!!!
    #if(randomgroups & is.null(groups)){
    #  newdat <- do.random.groups(data,groups)
    #  data <- newdat$dat
    #  grouping <- newdat$grouping
    #}
    ## generate all of the interactions
    allIntObj <- enumerateInteractions(mvars, poly)
    ## generate a formula
    if(poly==1){
      if(!allIntObj$isinf$total1 & !heur){
        myform <- genForm(allIntObj,treatment=treatment, data=data, poly=poly)
      } else {
        myform <- heuristicForm(allIntObj,treatment=treatment, data=data, poly=poly)
      }
    }
    if(poly==2){
      if(!allIntObj$isinf$total2 & !allIntObj$isinf$total1 & !heur){
        myform <- genForm(allIntObj,treatment=treatment, data=data, poly=poly)
      } else {
        myform <- heuristicForm(allIntObj,treatment=treatment, data=data, poly=poly)
      }
    }
    ## draw the caliper
    rcaliper <- round(runif(0,.5,n=1),5)
    ## do the pscore matching
    dist <- "logit"
    meth <- "nearest"
    mcall <- sprintf("matchit(formula = %s, data=data, caliper = %s, method = \"%s\", distance = \"%s\")",myform,rcaliper,meth,dist)
    m.out <- eval(parse(text = mcall))
    mcall2 <- sprintf("matchit(formula = %s, data=obj$data, caliper = %s, method = \"%s\", distance = \"%s\")",myform,rcaliper,meth,dist)
    ## collect the call
    mycalls[qq] <- mcall2
    ## calculate L1
    myL1 <- (my.imbalance(data[[treatment]],data,drop=c(treatment),breaks=rawCP, weights=m.out$weights))$L1$L1
    ## put things in their holders
    outN[qq,] <- m.out$n["Matched", "Treated"]
    outL1[qq] <- myL1
    outobs[[qq]] <- rownames(na.omit(m.out$match.matrix))
    if(progr & R>1){setTxtProgressBar(pb, qq)}
    #print(paste("PSM",qq,"of",R,": time =",round(Sys.time()-start.time,2)))
  }
  if(progr & R>1){close(pb)}
  ## return everything
  tab <- as.data.frame(matrix(NA, R, 2 * 2 + 4))
  g.names <- c("0","1")
  colnames(tab) <- c(paste("G", g.names, sep = ""), 
    					   paste("PercG", g.names, sep
                                                 = ""), "ML1",
                                           "Relaxed", "Method","Call")
  tab[,1:2] <- outN
  tab[,3] <- tab[,1]/table(data[[treatment]])[1]
  tab[,4] <- tab[,2]/table(data[[treatment]])[2]
  tab[,5] <- outL1
  tab[,6] <- rep("random",R)
  tab[,7] <- rep("matchit",R)
  tab[,8] <- mycalls

  out <- list(tab=tab,matched=outobs)
  return(out)
}

####################################################################






####################################################################
  ## enumerate the possible interactions
## This function takes a list of variable names and enumerates all
## the interactions.
## Note that it only works if there are 3 or more vars.

enumerateInteractions <- function(mvars, poly=1){
    #library(combinat,quietly=TRUE)
  ## interactions
  int <- combn(length(mvars),2)
  ## squared terms
  sq <- t(matrix(1:length(mvars),length(mvars),2))

  if(poly >=2){
    ## enumerate interactions
    interact.holder <- rep(NA,ncol(int))
    for(i in 1:ncol(int)){
      newvar <- paste(mvars[int[1,i]],":",mvars[int[2,i]], sep="")
      interact.holder[i] <- newvar
    }
    ## enumerate squared terms
    sq.holder <- rep(NA,ncol(sq))
    for(i in 1:ncol(sq)){
      newvar <- paste("I(",mvars[sq[1,i]],"*",mvars[sq[2,i]],")", sep="")
      sq.holder[i] <- newvar
    }
  }

  if(poly >=3){
    ## for three way interactions and cubed terms
    int3 <- combn(length(mvars),3)
    ## cubed terms
    cube <- t(matrix(1:length(mvars),length(mvars),3))

    ## Three way interactions
    int3.holder <- rep(NA,ncol(int3))
    for(i in 1:ncol(int3)){
      newvar <- paste(mvars[int3[1,i]],":",mvars[int3[2,i]],":",mvars[int3[3,i]], sep="")
      int3.holder[i] <- newvar
    }
    ## cubed terms
    cube.holder <- rep(NA, ncol(cube))
    for(i in 1:ncol(cube)){
      newvar <- paste("I(",mvars[cube[1,i]],"*",mvars[cube[2,i]],"*",mvars[cube[3,i]],")", sep="")
      cube.holder[i] <- newvar
    }
  }

  ## the list of variables to feed into the sampler
  ## first order interactions
  if(poly >=2){
    allvars <- c(mvars, interact.holder, sq.holder)
    ## just the interactions
    interax <- c(interact.holder, sq.holder)
  } else {
    allvars <- NA
    interax <- NA
  }
  ## second order interactions
  if(poly >=3){
    allvarscube <- c(mvars,interact.holder,int3.holder,sq.holder,cube.holder)
  } else {
    allvarscube <- NA
  }

  ## counting the number of main effects regressions to run
  total1 <- rep(0,length(mvars))
  for (i in 1:length(mvars)) {
    total1[i] <- nCm(length(mvars), i)
  }

  if(poly >=2){
    ## counting the number of first order int regressions to run
    total2 <- rep(0,length(allvars))
    for (i in 1:length(allvars)) {
      total2[i] <- nCm(length(allvars), i)
    }
  } else {
    total2 <- Inf
  }

  if(poly >=3){
    ## counting the number of first order int regressions to run
    total3 <- rep(0,length(allvarscube))
    for (i in 1:length(allvarscube)) {
      total3[i] <- nCm(length(allvarscube), i)
    }
  } else {
    total3 <- Inf
  }

  ## Note if things include Inf
  isinf <- c()
  isinf$total1 <- sum(total1==Inf)>0
  isinf$total2 <- sum(total2==Inf)>0
  isinf$total3 <- sum(total3==Inf)>0

  ## RETURN
  return(list(total1=total1, total2=total2, total3=total3,
              mvars=mvars,allvars=allvars,allvarscube=allvarscube,
              isinf=isinf, interax=interax))
}                               
####################################################################

####################################################################
  ## A function to exclude interactions of indicators from a formula
excludeIndicatorInteractions <- function(mvrs,tofactr,form){
  interacts.to.exclude <- paste("I(",mvrs[tofactr==1],"*",mvrs[tofactr==1],")",sep="")
  mfx <- interacts.to.exclude  ## store this for later
  for(iii in 1:length(mvrs[tofactr==1])){
    interacts.to.exclude <- c(interacts.to.exclude,
                            paste(mvrs[tofactr==1][iii],":",mvrs[tofactr==1],sep=""))
  }
    ## add three way interactions
  interacts.to.exclude <- c(interacts.to.exclude,paste("I(",mvrs[tofactr==1],"*",mvrs[tofactr==1],"*",mvrs[tofactr==1],")",sep=""))
  for(iii in 1:length(mfx)){
    tempstring <- paste(mvrs[tofactr==1][iii],":",mvrs[tofactr==1],sep="")
    for(jjj in 1:length(tempstring)){
      interacts.to.exclude <- c(interacts.to.exclude,paste(tempstring[jjj],":",mvrs[tofactr==1],sep=""))
    }
  }
      ## gsub out the bad interactions
  for(iii in 1:length(interacts.to.exclude)){
    form <- gsub(paste(" + ",interacts.to.exclude[iii],sep=""),"",form, fixed = T)
  }
  return(form)
}
####################################################################

####################################################################
## This takes the list of all possible interactions and creates a
## randomly selected formla for the matching

## the argument "allIntObj" is an object created by
## "enumerateInteractions".  The argument "poly" can be 1, 2, or 3
## and says whether to use main effects only (1), squared terms (2),
## or cubed terms (and lower order terms) (3).  THE OPTION POLY=3
## DOES NOT YET EXIST!!!
## The option poly = 4
## gives you a small number of first order interactions.
## THE OPTION POLY = 4 DOES NOT YET EXIST!!!

genForm <- function(allIntObj, poly=2, treatment,data=data){
  if(poly!=1 & poly!=2){stop
                        print("poly must be 1 or 2")}
  ## start the
  ftmp <- paste(treatment, "~ 1")
  if(poly==1){
    nvars <- which(rmultinom(1,1,(allIntObj$total1/(sum(allIntObj$total1)))) > 0)
    bk <- sample(allIntObj$mvars,nvars,replace=F)
    form <- paste(ftmp, "+",paste(bk,collapse=" + "))
  }
  if(poly==2){
    nvars <- which(rmultinom(1,1,(allIntObj$total2/(sum(allIntObj$total2)))) > 0)
    bk <- sample(allIntObj$allvars,nvars,replace=F)
    form <- paste(ftmp, "+",paste(bk,collapse=" + "))

    tofactor <- rep(0,length(allIntObj$mvars))
    for(i in 1:length(allIntObj$mvars)){
      if(nlevels(as.factor(data[[allIntObj$mvars[i]]]))<3){
        tofactor[i] <- 1
      }
    }
    form <- excludeIndicatorInteractions(form=form, tofactr=tofactor,
                                         mvrs=allIntObj$mvars)
  }
  return(form)
}



####################################################################
 
####################################################################
## This takes the list of all possible interactions and creates a 
## formula for matching heuristically

## the argument "allIntObj" is an object created by

heuristicForm <- function(allIntObj, poly=2, treatment, data=data){
  if(poly!=1 & poly!=2){stop
                        print("poly must be 1 or 2")}

  ## start the
  ftmp <- paste(treatment, "~ 1")

  if(poly==1){
    weights.mfx <- (dpois(1:length(allIntObj$mvars), length(allIntObj$mvars)))/sum(dpois(1:length(allIntObj$mvars),length(allIntObj$mvars)))
    nvars.mfx <- which(rmultinom(1,1,weights.mfx) > 0)
    bk.mfx <- sample(allIntObj$mvars,nvars.mfx,replace=F) 
    form <- paste(ftmp, "+",paste(bk.mfx,collapse=" + "))
  }
  if(poly==2){
    weights.mfx <- (dpois(1:length(allIntObj$mvars), length(allIntObj$mvars)))/sum(dpois(1:length(allIntObj$mvars),length(allIntObj$mvars)))
    nvars.mfx <- which(rmultinom(1,1,weights.mfx) > 0)
    bk.mfx <- sample(allIntObj$mvars,nvars.mfx,replace=F)

    #m <- length(allIntObj$interax)
    #weights.interx <- (dpois(1:(m+1), m/15))/sum(dpois(1:(m+1), m/15))
    weights.interx <- (dpois(1:length(allIntObj$mvars), length(allIntObj$mvars)/8))/sum(dpois(1:length(allIntObj$mvars),length(allIntObj$mvars)/8))
      ## Changing the digit in the line above (twice) will change the number of interactions
      ## that are typically drawn.  Smaller numbers mean more interactions.
    nvars.interx <- which(rmultinom(1,1,weights.interx) > 0)
    bk.interx <- sample(allIntObj$interax,nvars.interx,replace=F)

    form <- paste(ftmp, "+",paste(bk.mfx,collapse=" + "), "+",paste(bk.interx,collapse=" + "))
    
    tofactor <- rep(0,length(allIntObj$mvars))
    for(i in 1:length(allIntObj$mvars)){
      if(nlevels(as.factor(data[[allIntObj$mvars[i]]]))<3){
        tofactor[i] <- 1
      }
    }
    form <- excludeIndicatorInteractions(form=form, tofactr=tofactor,
                                         mvrs=allIntObj$mvars)
  }
  return(form)
}



####################################################################
 

####################################################################
psm <- function(formula, data, ratio=1, caliper=0, linear.pscore=FALSE,rawCP){
    ## estimate the propensity score
  if(nrow(na.omit(data))!=nrow(data)){
    print("Missing values in propensity score model.  Observations with missing data ignored.")
    data <- na.omit(data)
  }
  mymod <- (glm(formula, data, family=(binomial(link="logit"))))
  
  if(mymod$converged==F){break}
  ps <- mymod$fitted.values
  if(linear.pscore==T){
    ps <- log(ps/(1-ps))
  }
    ## set the number of matches, 1-to-k
  k <- ratio

  tvar <-  as.character(as.formula(formula))[2]
  Tnms <- row.names(data)[data[[tvar]]==1]
  Cnms <- row.names(data)[data[[tvar]]==0]

    ## make some holders
  match.matrix <- matrix(NA, length(Tnms), k)
  rownames(match.matrix) <- Tnms
  dist.matrix <- match.matrix

  t0 <- Sys.time()
  ixnay <- c()
    ## loop over the match ratio
  for(j in 1:k){
      ## loop over the treated units
    for(i in 1:length(Tnms)){
      #if(i%%100==0){
      #  print(paste("Match",j,"of",k,", T",i,"of",length(Tnms),":",round(Sys.time()-t0,2)))
      #}
        ## define a vector of the Cnms that are still available
      elig <- Cnms[(Cnms%in%ixnay)==F]
      #if(i%%500==0){print(paste("treated unit",i,"of",length(Tnms)))}
        ## make the distance vector for treated t and all the eligable
        ## controls
      x <-abs(ps[Tnms[i]] - ps[elig])
        ## capture the colnames of the ones that have
        ## minimum distances and take the first
      min.x <- min(x)
      m <- (elig[(x==min.x)==T])[1]
        ## add the things to the output matrices
      match.matrix[i,j] <- m
      dist.matrix[i,j] <- min.x
        ## add the new match to ixnay
      ixnay <- c(ixnay, m)
        ## stop the iterations if there are no more control units
      if(length(ixnay)==length(Cnms)){break}
    }
    ## stop the iterations if there are no more control units
    if(length(ixnay)==length(Cnms)){break}
  }


  ## make the matched data
    ## first, select all the treated obs
  mh.data <- data[rownames(match.matrix)[is.na(match.matrix[,1])==F],]
    ## then add the control obs
  for(i in 1:ncol(match.matrix)){
    mh.data <- rbind(mh.data,data[na.omit(match.matrix[,i]),])
  }


  ## add distances -- mean of all distances for a matched group
  dd <- apply(dist.matrix, MARGIN=1, mean, na.rm=T)
  distance <- matrix(NA,nrow(mh.data),1)
  rownames(distance) <- rownames(mh.data)
  distance[names(na.omit(dist.matrix[,1])),] <- dd[names(na.omit(dist.matrix[,1]))]
  for(i in 1:ncol(dist.matrix)){
    names(dd) <- match.matrix[,i]
    distance[na.omit(match.matrix[,i]),] <- dd[na.omit(match.matrix[,i])]
  }
  mh.data$distance <- distance

    ## add weights to the data
  weights <- matrix(NA,nrow(mh.data),1)
  rownames(weights) <- rownames(mh.data)
  weights[names(na.omit(match.matrix[,1])),]<-1
    ## calculate the weights
  countf <- function(x){sum(is.na(x))}
  wts <- 1/(k-apply(match.matrix, MARGIN=1, countf))
    ## fill in the weights
    ## this has me all mixed up, but I think it works
  for(i in 1:ncol(dist.matrix)){
    names(wts) <- match.matrix[,i]
    weights[na.omit(match.matrix[,i]),] <-wts[na.omit(match.matrix[,i])]
  }
  mh.data$weights <- weights

  ## Begin calipering
  match.obj <- list("match.matrix"=match.matrix, "dist.matrix"=dist.matrix,
              "match.dat"=mh.data)
  rc <- randomCaliper(match.obj, tvar, rawCP=rawCP, caliper=caliper)

  ## Need to make calipered match matrix, data, etc.
  keepobs <- rownames(rc$caliperdat)
  match.matrix <- as.matrix(match.matrix[which(rownames(match.matrix) %in% keepobs),])
  dist.matrix <- as.matrix(dist.matrix[which(rownames(dist.matrix) %in% keepobs),])
  
 
    ## RETURN
  out <- list("match.matrix"=match.matrix, "dist.matrix"=dist.matrix,
              "match.dat"=rc$caliperdat,call=match.call(),
              N=rc$N)
  class(out) <- "psm.match"
  return(out)
}

####################################################################




####################################################################
## This function takes a "match.obj" for either myMH or myPS and
## applies calipers.

  ## 1-to-1 matching ONLY at the moment.
randomCaliper <- function(match.obj, treatment,rawCP, caliper=0){
  mdat <- match.obj$match.dat
     ## pull out the solution
  solution <- match.obj
     ## make the dataset
  t.units <- rownames(solution$match.matrix)
  c.units1 <- solution$match.matrix[,1]
    ## this makes sure the names in the match.matrix get matched to
    ## the data rownames correctly
  matchfind <- function(x,vec){
    if(is.na(x)){return(NA)}
    if(!is.na(x)){return(which(vec==x))}
  }
  t.rows <- apply(as.matrix(t.units),MARGIN=1,FUN=matchfind,vec=rownames(mdat))
  if(is.list(t.rows)){
    t.rows <- unlist(t.rows)
  }
  c.rows <- apply(as.matrix(c.units1),MARGIN=1,FUN=matchfind,vec=rownames(mdat))
    ## Then, make sure to only include treated units that have at
    ## least one match
  matched <- which(is.na(solution$match.matrix[,1])==F)
    ## then combine them to make the treated data matrix
  tdat <- mdat[t.rows[matched],]
    ## the cdat1 will be the same, because it has to be 1-to-1
  cdat1 <- mdat[c.rows[matched],]
  psdiff <- tdat$distance  ## diff from pscore

    ## order the datasets so that they are in order of best
    ## to worst matches
  tdat <- tdat[order(psdiff),]
  cdat1 <- cdat1[order(psdiff),]

    ## calculate the ATT, L1s, etc

  #if(exists(psrunsR[i])){
    if(caliper > (nrow(tdat) - 1)){
      caliper <- (nrow(tdat) - 1)
      print("Caliper exceeds number of matches.  Automatically reset to the maximum possible")
    }
    caliperdat <-  na.omit(rbind(tdat[1:(nrow(tdat)-caliper),],
                     cdat1[1:(nrow(cdat1)-caliper),]))
      ## calculate the ATT, L1, etc
      ## this now calls "my.imbalance" which takes out the chi-sq
      ## calculation that was leading to bugs.
    #L1 <- (imbalance(caliperdat[[treatment]],caliperdat,drop=c(treatment,"distance","weights"),breaks=rawCP))$L1$L1
    #L1 <- (my.imbalance(caliperdat[[treatment]],caliperdat,drop=c(treatment,"distance","weights"),breaks=rawCP))$L1$L1
    N <- table(caliperdat[[treatment]])
    matched.rows <- rownames(caliperdat)

  return(list(N=N,matched.rows=matched.rows, caliperdat=caliperdat))

}



####################################################################


###################################################################

## These next functions do the random grouping

############################################################3
## From "cem.R" in cegg, this function groups variables
group.var <- function(x, groups){
    n <- length(x)
#    print(str(x))
    tmp <- numeric(n)
    ngr <- length(groups)
    all.idx <- NULL
#    cat(sprintf("n=%d, ngr=%d\n", n, ngr))
     for(i in 1:ngr){
      idx <- which(x %in% groups[[i]])
       if(length(idx)>0){
        tmp[idx] <- sprintf("g%.2d",i)
        all.idx <- c(all.idx, idx)
       }
      }
      if(length(all.idx)>0 && length(all.idx)<n)
       tmp[-all.idx] <- x[-all.idx]
    tmp <- factor(tmp)
        return(tmp)
}


## This is a modified chunk from the begining of the cem.R code
group.data <- function(data,grouping){
        gn <- names(grouping)
        n.gn <- length(gn)
        nd <- 0
        for (g in 1:n.gn) {
                data[[gn[g]]] <- group.var(data[[gn[g]]], grouping[[g]])

        }
        return(data)
    }


## Random grouping code
## This generates random groupings by splitting the groups
## at a uniformly generated split point, then repeating for
## the remaining until none are left.
## This mixes up the order, so it assumes they are unordered.
## This function forces there to be at least 2 groups -- the
## random formula generation process will make sure we just don't
## condition on some of the covariates
random.group <- function(x){
  lev <- unique(as.character(x))
  grp <- list()
  counter <- 0
  while(length(lev)>0){
    counter <- counter + 1
    n <- length(lev)
    lev.mixed <- sample(lev,length(lev),replace=F)
    if(counter==1){
      s <- sample(1:(n-1),1)
    } else {
      s <- sample(1:n,1)
    }
    samp <- lev.mixed[1:s]
    grp[[counter]] <- samp
    levnew <- lev[-which(lev%in%samp)]
    lev <- levnew
  }
  return(grp)
}


##############################################
## a function that does random groups with specified
## groups.

random.group2 <- function(x, xgroups=NULL){
  lev <- c(unique(as.character(x)),xgroups)
  grp <- list()
  counter <- 0
  while(length(lev)>0){
    counter <- counter + 1
    ## sample one of the groups
    samp <- unlist(sample(lev,1))
    grp[[counter]] <- samp
    ## then eliminate any other group that has the same categorie(s)
    overlapping <- rep(F,length(lev))
    for(ll in 1:length(lev)){
      tmp <- rep(F,length(samp))
      for(ss in 1:length(samp)){
        tmp[ss] <- samp[ss]%in%unlist(lev[ll])
      }
      overlapping[ll] <- max(tmp)
    } 
    levnew <- lev[-which(overlapping==1)]
    lev <- levnew
  }
  return(grp)
}

########################################

## A function to identify all the factor variables
isfactorvar <- function(dat){
  isfactor <- rep(NA,ncol(dat))
  for(i in 1:ncol(dat)){
    isfactor[i] <- is.factor(dat[,i])
  }
  return(isfactor)
}
 
#############################################

## do random factoring for the factor variables
## and put it in a grouping

## THIS FUNCTION IS NOW REPLACED BY "do.random.groups2()"
do.random.groups <- function(dat,...){
  dotargs=as.list(sys.call())
  mvars <- dotargs$mvars
  grouping <- list()
  counter <- 0
  namevec <- c()
  for(ii in 1:ncol(dat)){
    if(isfactorvar(dat)[ii] & (names(dat)%in%mvars)[ii]){
      counter <- counter+1
      rg <- random.group(dat[[ii]])
      grouping[[counter]] <- rg
      namevec <- c(namevec,names(dat)[ii])
    }
    names(grouping) <- namevec
  }

  ## Then apply the random grouping to the variables
  newdat <- group.data(dat,grouping)
  dat <- newdat
  return(list(dat=dat,grouping=grouping))
}



##################################################
## do random factoring where some or all have only
## certain groups that are ok

do.random.groups2 <- function(dat,groups=NULL,tvar = treatment){
  dotargs=as.list(sys.call())
  mvars <- dotargs$mvars
  treatment <- dotargs$treatment
  grouping <- list()
  counter <- 0
  namevec <- c()
  idx <- match(tvar,colnames(dat))
  mvars <- names(dat)[-idx]

  specialvars <- names(groups)
  for(ii in 1:ncol(dat)){
    if(isfactorvar(dat)[ii] & (names(dat)%in%mvars)[ii] & (names(dat)[ii]%in%names(groups))==F){
      counter <- counter+1
      rg <- random.group(dat[[ii]])
      grouping[[counter]] <- rg
      namevec <- c(namevec,names(dat)[ii])
    }
    if(isfactorvar(dat)[ii] & (names(dat)%in%mvars)[ii] & (names(dat)[ii]%in%names(groups))){
      counter <- counter+1
      rg <- random.group2(dat[[ii]],groups[[names(dat)[ii]]])
      grouping[[counter]] <- rg
      namevec <- c(namevec,names(dat)[ii])
    }
    names(grouping) <- namevec
  }

  ## Then apply the random grouping to the variables
  newdat <- group.data(dat,grouping)
  dat <- newdat
  return(list(dat=dat,grouping=grouping))
}



##########################################################################


##################################################################
  ## a mahalanobis function
  ## myMH function from http://www.stat.lsa.umich.edu/~bbh/optmatch/doc/mahalanobisMatching.pdf

  ## the second stopifnot condition was in ben hansen's code
  ## but the code seems to work fine for even one variable
  ## and I have samples with just one covariate.
myMH <- function(Tnms, Cnms, inv.cov, data) {
 stopifnot(!is.null(dimnames(inv.cov)[[1]]),# dim(inv.cov)[1] > 1,
 all.equal(dimnames(inv.cov)[[1]], dimnames(inv.cov)[[2]]),
 all(dimnames(inv.cov)[[1]] %in% names(data)))
 covars <- dimnames(inv.cov)[[1]]
 xdiffs <- as.matrix(data[Tnms, covars])
 xdiffs <- xdiffs - as.matrix(data[Cnms, covars])
 rowSums((xdiffs %*% inv.cov) * xdiffs)
}



  ## This is the FAST m-dist function for large data
myMH4 <- function(Tnms, Cnms, inv.cov, data,k){
    ## the basic idea now is that we do the distance calculation
    ## and matching all in one step, so that we don't have to
    ## make a huge distance matrix and then cycle through it.
    ## This also means that as we get fewer and fewer available
    ## control units, the distance calculation should get faster.
  if(nrow(na.omit(data))!=nrow(data)){
    print("Missing values ignored by Mahalanobis matching.")
    data <- na.omit(data)
  }

  icv <- inv.cov
  mdist <- matrix(NA,length(Tnms),length(Cnms))
  ixnay <- c()
  t0 <- Sys.time()
    ## make some holders
  match.matrix <- matrix(NA, length(Tnms), k)
  rownames(match.matrix) <- Tnms
  dist.matrix <- match.matrix
    ## loop over the match ratio
  for(j in 1:k){
      ## loop over the treated units
    for(i in 1:length(Tnms)){
      if(i%%100==0){
       # print(paste("Match",j,"of",k,", T",i,"of",length(Tnms),":",round(Sys.time()-t0,2)))
      }
        ## define a vector of the Cnms that are still available
      elig <- Cnms[(Cnms%in%ixnay)==F]
      #if(i%%500==0){print(paste("treated unit",i,"of",length(Tnms)))}
        ## make the distance vector for treated t and all the eligable
        ## controls
      x <- outer(Tnms[i], elig, FUN = myMH, inv.cov = icv, data = data)
        ## capture the colnames of the ones that have
        ## minimum distances and take the first
      min.x <- min(x)
      m <- (elig[(x==min.x)==T])[1]
        ## add the things to the output matrices
      match.matrix[i,j] <- m
      dist.matrix[i,j] <- min.x
        ## add the new match to ixnay
      ixnay <- c(ixnay, m)
        ## stop the iterations if there are no more control units
      if(length(ixnay)==length(Cnms)){break}
    }
    ## stop the iterations if there are no more control units
    if(length(ixnay)==length(Cnms)){break}
  }
  return(list("match.matrix"=match.matrix,
              "dist.matrix"=dist.matrix))
}




mdm <- function(formula, data, ratio=1, caliper=0, rawCP){

  ## I think this function has a dependence on "tvar"

    ## pull out the covariates
  X1 <- data
  X2 <- data.frame(model.frame(formula, data))
  tvar <- as.character(formula)[2]
    ## just the covariates
  X3 <- subset(data.frame(model.matrix(formula,data)),select=-1)
  mvrs <- colnames(X3)
    ## covariates and treatment variable
  X <- cbind(X2[,1],X3)
  names(X) <- c(tvar,mvrs)

    ## set the number of matches, 1-to-k
  k <- ratio

  ## from http://tolstoy.newcastle.edu.au/R/help/05/05/3996.html
  matrix.rank <- function(A, eps=sqrt(.Machine$double.eps)){
    sv. <- abs(svd(A)$d)
    sum((sv./max(sv.))>eps)
  }
  if(dim(X)[2]!=matrix.rank(X)){
    #print("MDM failed: The matrix is not of full rank.")
  }
  stopifnot(dim(X)[2]==matrix.rank(X), silent=T)

    ## if it can't invert, this is where it will fail
  try(icv <- solve(cov(subset(X,select=mvrs))))
  #if(exists("icv")==F){
  #  stop
  #  print("VCOV matrix can't invert with these covariates")
  #}
  stopifnot(exists("icv"))
  trtnms <- row.names(X)[X[[tvar]]==1]
  ctlnms <- row.names(X)[X[[tvar]]==0]
    ## use the data X3 that has just the covariates

    ## run the main internal function
  try(mydist <- myMH4(trtnms, ctlnms, inv.cov=icv, data=X3,k=k))
  stopifnot(exists("mydist"))
  match.matrix <- mydist$match.matrix
  dist.matrix <- mydist$dist.matrix

  ## make the matched data
    ## first, select all the treated obs
  mh.data <- data[rownames(match.matrix)[is.na(match.matrix[,1])==F],]
    ## then add the control obs
  for(i in 1:ncol(match.matrix)){
    mh.data <- rbind(mh.data,data[na.omit(match.matrix[,i]),])
  }

  ## PROBLEM HERE WITH THE DISTANCES -- NEED TO GET THIS RIGHT

  ## add distances -- mean of all distances for a matched group
  dd <- apply(dist.matrix, MARGIN=1, mean, na.rm=T)
  distance <- matrix(NA,nrow(mh.data),1)
  rownames(distance) <- rownames(mh.data)
  distance[names(na.omit(dist.matrix[,1])),] <- dd[names(na.omit(dist.matrix[,1]))]
  for(i in 1:ncol(dist.matrix)){
    names(dd) <- match.matrix[,i]
    distance[na.omit(match.matrix[,i]),] <- dd[na.omit(match.matrix[,i])]
  }
  mh.data$distance <- distance

    ## add weights to the data
  weights <- matrix(NA,nrow(mh.data),1)
  rownames(weights) <- rownames(mh.data)
  weights[names(na.omit(match.matrix[,1])),]<-1
    ## calculate the weights
  countf <- function(x){sum(is.na(x))}
  wts <- 1/(k-apply(match.matrix, MARGIN=1, countf))
    ## fill in the weights
    ## this has me all mixed up, but I think it works
  for(i in 1:ncol(dist.matrix)){
    names(wts) <- match.matrix[,i]
    weights[na.omit(match.matrix[,i]),] <-wts[na.omit(match.matrix[,i])]
  }
  mh.data$weights <- weights

  ## Begin calipering
  match.obj <- list("match.matrix"=match.matrix, "dist.matrix"=dist.matrix,
              "match.dat"=mh.data)
  rc <- randomCaliper(match.obj, tvar, rawCP=rawCP, caliper=caliper)

  ## Need to make calipered match matrix, data, etc.
  keepobs <- rownames(rc$caliperdat)
  match.matrix <- as.matrix(match.matrix[which(rownames(match.matrix) %in% keepobs),])
  dist.matrix <- as.matrix(dist.matrix[which(rownames(dist.matrix) %in% keepobs),])
  
 
    ## RETURN
  out <- list("match.matrix"=match.matrix, "dist.matrix"=dist.matrix,
              "match.dat"=rc$caliperdat,call=match.call(),
              N=rc$N)
  class(out) <- "mdm.match"
  return(out)
}


## End of mahalanobis function
###################################################################

###################################################################
## Define my own imbalance function
## This is a modification of the imbalance function in cem
## It cuts out the chi-square table that both produces errors and
## causes a bug every so often

my.imbalance <- function(group, data, drop=NULL, breaks=NULL, weights, grouping = NULL){
 if (!is.data.frame(data))
        stop("Data must be a dataframe", call. = FALSE)

  if(!is.null(drop)){
   drop <- unique(drop)
   dropped <- match(drop, colnames(data))
   dropped <- dropped[!is.na(dropped)]

   if(length(dropped)>0)
        data <- data[-dropped]
  }
 if(is.null(breaks)){
  vars <- colnames(data)
  nv <- length(vars)
  breaks <- vector(nv, mode="list")
  for(i in 1:nv){
   if(is.numeric(data[[i]]) | is.integer (data[[i]]))
    breaks[[i]] <- pretty(range(data[[i]],na.rm=TRUE), n=nclass.scott(na.omit(data[[i]])), 1)
   names(breaks) <- vars
  }
 } 

 if(!is.null(grouping) & !is.null(names(grouping))){
    gn <- names(grouping)
    n.gn <- length(gn)
    for(g in 1:n.gn){
        if(!is.null(data))
            data[[gn[g]]] <- group.var(data[[gn[g]]], grouping[[g]])
        if(!is.null(breaks[gn[g]]))
            breaks[gn[g]] <- NULL                 
    }
 }
   
 n <- dim(data)[1]
 if(missing(weights)){
  weights <- rep(1, n)
 }
 rem <- which(weights<=0)
 if(length(rem)>0){
  data <- data[-rem,]
  weights <- weights[-rem]
  group <- group[-rem]
 }


 tmp <- reduce.data(data, breaks=breaks)$data
 lv <- unique(na.omit(group))
   
 globalL1 <- L1.meas(group=group, data=data, breaks=breaks, weights=weights)

 out <- list(L1=globalL1)

 ## I cut out a bunch of stuff here, becuase I don't want to make the table

    
 class(out) <- "imbalance"
 return( out )
 
}

## end of my imbalance function
###################################################################


###################################################################
## Reduce data function
## I think this can be deleted once it's in the cem package because
## then the function is available to me.

`reduce.var` <-
function(x, breaks){
	if(is.numeric(x) | is.integer(x)){
	 if(is.null(breaks)){
	  breaks <- "sturges"
	  }
	 if(is.character(breaks)){
       breaks <- match.arg(tolower(breaks), c("sturges", 
                "fd", "scott", "ss"))
            breaks <- switch(breaks, sturges = nclass.Sturges(x), 
                 fd = nclass.FD(x), 
				 scott = nclass.scott(x), 
				 ss = nclass.ss(x),
                stop("unknown 'breaks' algorithm"))
        }
	 if(length(breaks) > 0){
		if(length(breaks)==1){
			rg <- range(x, na.rm=TRUE)
			breaks <- seq(rg[1],rg[2], length = breaks)
		}
		breaks <- unique(breaks)
		if(length(breaks)>1)
	     x <- cut(x, breaks=breaks, include.lowest = TRUE, labels = FALSE)
		else 	
		 x <- as.numeric(x) 
	 }
	} else {
	  x <- as.numeric(x) 
	}
	return(list(x=x, breaks=breaks)) 
}

reduce.data <- function(data, breaks=NULL, collapse=FALSE){
  if (!is.data.frame(data))
        stop("Data must be a dataframe", call. = FALSE)
 vnames <- colnames(data)
 nv <- length(vnames)
 new.breaks <- vector(dim(data)[2], mode="list")
 names(new.breaks) <- vnames
 for (i in 1:nv){
   tmp <- reduce.var(data[[i]], breaks[[vnames[i]]] )
   new.breaks[[vnames[i]]] <- tmp$breaks
   data[[i]] <- tmp$x
  }
 if(collapse)
  return(list(data=collapse.data(data), breaks=new.breaks))
 
 return(list(data=data, breaks=new.breaks))
}

collapse.data <- function(data){
  apply(data,1, function(x) paste(x, collapse="\r"))	
}

##################################################################

##################################################################
## A function that calculates Mahalanobis distances without matching

mdistfun <- function(formula, data){
  ## I think this function has a dependence on "tvar"

  ## a mahalanobis function
  ## myMH function from http://www.stat.lsa.umich.edu/~bbh/optmatch/doc/mahalanobisMatching.pdf

  ## the second stopifnot condition was in ben hansen's code
  ## but the code seems to work fine for even one variable
  ## and I have samples with just one covariate.
  myMH <- function(Tnms, Cnms, inv.cov, data) {
   stopifnot(!is.null(dimnames(inv.cov)[[1]]),# dim(inv.cov)[1] > 1,
   all.equal(dimnames(inv.cov)[[1]], dimnames(inv.cov)[[2]]),
   all(dimnames(inv.cov)[[1]] %in% names(data)))
   covars <- dimnames(inv.cov)[[1]]
   xdiffs <- as.matrix(data[Tnms, covars])
   xdiffs <- xdiffs - as.matrix(data[Cnms, covars])
   rowSums((xdiffs %*% inv.cov) * xdiffs)
  }
    ## pull out the covariates
  X1 <- data
  X2 <- data.frame(model.frame(formula, data))
  tvar <- as.character(formula)[2]
    ## just the covariates
  X3 <- subset(data.frame(model.matrix(formula,data)),select=-1)
  mvrs <- colnames(X3)
    ## covariates and treatment variable
  X <- cbind(X2[,1],X3)
  names(X) <- c(tvar,mvrs)

  ## from http://tolstoy.newcastle.edu.au/R/help/05/05/3996.html
  matrix.rank <- function(A, eps=sqrt(.Machine$double.eps)){
    sv. <- abs(svd(A)$d)
    sum((sv./max(sv.))>eps)
  }
  #if(dim(X)[2]!=matrix.rank(X)){
  #  #print("MDM failed: The matrix is not of full rank.")
  #}
  #stopifnot(dim(X)[2]==matrix.rank(X), silent=T)

    ## if it can't invert, this is where it will fail
  try(icv <- solve(cov(subset(X,select=mvrs))))
  stopifnot(exists("icv"))
  trtnms <- row.names(X)[X[[tvar]]==1]
  ctlnms <- row.names(X)[X[[tvar]]==0]
    ## use the data X3 that has just the covariates

    ## run the main internal function
  mydist <- outer(trtnms, ctlnms, FUN = myMH, inv.cov = icv, data = X3)
  stopifnot(exists("mydist"))
  dimnames(mydist) <- list(trtnms, ctlnms)
  return(mydist)
}

## a version that works for bigger datasets



mdistfun2 <- function(formula, data){
  ## I think this function has a dependence on "tvar"

  ## a mahalanobis function
    myMH <- function(Tnms, Cnms, inv.cov, data) {
   stopifnot(!is.null(dimnames(inv.cov)[[1]]),# dim(inv.cov)[1] > 1,
   all.equal(dimnames(inv.cov)[[1]], dimnames(inv.cov)[[2]]),
   all(dimnames(inv.cov)[[1]] %in% names(data)))
   covars <- dimnames(inv.cov)[[1]]
   xdiffs <- as.matrix(data[Tnms, covars])
   xdiffs <- xdiffs - as.matrix(data[Cnms, covars])
   rowSums((xdiffs %*% inv.cov) * xdiffs)
  }
    ## This is the FAST m-dist function for large data
myMHbig <- function(Tnms, Cnms, inv.cov, data){
  ## This is the FAST m-dist function for large data
  icv <- inv.cov
  mdist <- matrix(NA,length(Tnms),length(Cnms))
  ## loop over the treated units
  for(i in 1:length(Tnms)){
    ## define a vector of the Cnms that are still available
    elig <- Cnms
    ## make the distance vector for treated t and all the eligable
    ## controls
    x <- outer(Tnms[i], elig, FUN = myMH, inv.cov = icv, data = data)
    ## capture the colnames of the ones that have
    ## minimum distances and take the first
    mdist[i,] <- x
  }
  return(mdist)
}
#####################################

    ## pull out the covariates
  X1 <- data
  X2 <- data.frame(model.frame(formula, data))
  tvar <- as.character(formula)[2]
    ## just the covariates
  X3 <- subset(data.frame(model.matrix(formula,data)),select=-1)
  mvrs <- colnames(X3)
    ## covariates and treatment variable
  X <- cbind(X2[,1],X3)
  names(X) <- c(tvar,mvrs)

  ## from http://tolstoy.newcastle.edu.au/R/help/05/05/3996.html
  matrix.rank <- function(A, eps=sqrt(.Machine$double.eps)){
    sv. <- abs(svd(A)$d)
    sum((sv./max(sv.))>eps)
  }
  #if(dim(X)[2]!=matrix.rank(X)){
  #  #print("MDM failed: The matrix is not of full rank.")
  #}
  #stopifnot(dim(X)[2]==matrix.rank(X), silent=T)

    ## if it can't invert, this is where it will fail
  try(icv <- solve(cov(subset(X,select=mvrs))))
  stopifnot(exists("icv"))
  trtnms <- row.names(X)[X[[tvar]]==1]
  ctlnms <- row.names(X)[X[[tvar]]==0]
    ## use the data X3 that has just the covariates

    ## run the main internal function
  #mydist <- outer(trtnms, ctlnms, FUN = myMH, inv.cov = icv, data =
  #X3)
  mydist <- myMHbig(trtnms,ctlnms,icv,X3)
  
  stopifnot(exists("mydist"))
  dimnames(mydist) <- list(trtnms, ctlnms)
  return(mydist)
}
##################################################################

##################################################################
## Plotting function

#plot.spacegraph <- function(...) spacegraph.plot(...)

plot.spacegraph <- function(...){
   ## This is a bit more complicated but it allows me to pass
   ## the name of the object.  Otherwise, the object ends up being
   ## named "..1"
   obj.name <- as.character(match.call())[2]
   ## Grab the other things
   plotcall <- match.call()
   spacegraph.plot(objname=obj.name, plot.call=plotcall,...)
}

spacegraph.plot2 <- function(obj, objname=NULL, group="1", 
                            method.col=list(cem="#0000ff75",mdm="#ff000075",
                                            psm="#00ff0075",matchit="#FF880075"),
                            ...){
	if(class(obj) != "spacegraph")
	  stop("obj must be of class `spacegraph'")

        if(is.null(objname)){
          obj.name <- match.call()["obj"]
          obj.name <- gsub("()","",obj.name, fixed=T)
        } else {
          obj.name <- objname
        }
        
        data <- obj$data
        obj$space <- na.omit(obj$space)  ## this solves a problem
                                         ## when some solutions have NA in the table
	g <- sprintf("G%s",group)
	n <-  obj$space[[g]]
	ML1 <- obj$space$ML1
	Relaxed <- obj$space$Relaxed
        Method <- obj$space$Method
        Call <- obj$space$Call
        ## rename the stuff in the call
        for(i in 1:length(Call)){
          Call[i] <- gsub("obj$",paste(obj.name,"$",sep=""),Call[i],fixed=T)
        }
	class <- rep("relax", length(n))
	id.raw <- which(Relaxed=="<raw>")
	id.start <- which(Relaxed=="<start>") 
	class[ id.raw ] <- "raw"
	if(length(id.start)>0)
		class[ id.start ] <- "start"
	ids <- (1:length(n))[-c(id.raw, id.start)]
	name.vars <- names(obj$coars[[1]])
	n.vars <- length(name.vars)
	tab <- data.frame(n=n, ML1=ML1, class=class)
	main.txt <- sprintf("Total units=%d, ML1=%.3f", n[id.raw], ML1[id.raw])
	
	if(length(id.start)>0)
		main.txt <- sprintf("Initial matched units=%d, ML1=%.3f", n[id.start], ML1[id.start])

        ## make different colors of dots for the different methods
        dotcol <- rep("cyan",length(ids))
        ## get colors for other methods...
        approved.methods <- c("raw","cem","psm","mdm","matchit")
        needcol <- unique(Method)[unique(Method) %in% approved.methods==F]


        ## I decided that it's easier to specify my own colors
        colz <- c("#CC00FF75","darkolivegreen","gray50",
                  "lavenderblush","peru")
        rcol <- colz[1:length(needcol)]
        
        ll <- length(approved.methods)-1
        if(length(needcol)>0){
          for(i in (1:length(needcol) + ll)){
            pnames <- names(method.col)
            method.col[[i]] <- rcol[i-ll]
            names(method.col) <- c(pnames,needcol[i-ll])
          }
        }

        for(i in 1:length(method.col)){
          dotcol[Method==names(method.col)[i]] <-  method.col[[names(method.col)[i]]]
        }
       
	#dotcol <- rgb(0.2,0.2,0.8)
	#selcol <- rgb(0.8,0.8,0.2)
        selcol <- "black"
      ## set the plotting args
      myxlab <- "number matched, scaled as 1/sqrt(matched)"
      myylab <- "median of L1 profile"
      myylim <- c(0,range(ML1)[2])
      myxlim <- range(1/sqrt(n))
      ## get user-specified xlim and ylim values
      #if(match.call()["xlim"] != "NULL()"){
      #  xli <- gsub("()","",match.call()["xlim"],fixed=T)
      #  xli2 <- strsplit(xli,",")[[1]]
      #  lim1 <- as.numeric(gsub(" *","",gsub("c(","",xli2[1], fixed=T)))
      #  lim2 <- as.numeric(gsub(" *","",gsub(")","",xli2[2], fixed=T)))
      #  myxlim <- c(1/sqrt(lim1),1/sqrt(lim2))
      #}
      #if(match.call()["ylim"] != "NULL()"){
      #  xli <- gsub("()","",match.call()["ylim"],fixed=T)
      #  xli2 <- strsplit(xli,",")[[1]]
      #  lim1 <- as.numeric(gsub(" *","",gsub("c(","",xli2[1], fixed=T)))
      #  lim2 <- as.numeric(gsub(" *","",gsub(")","",xli2[2], fixed=T)))
      #  myylim <- c(lim1,lim2)
      #}

      ## collect any user-specified plotting args that conflict:
      #if(match.call()["main"] != "NULL()"){main.txt <- gsub("()","",match.call()["main"],fixed=T)}
      #if(match.call()["xlab"] != "NULL()"){myxlab <- gsub("()","",match.call()["xlab"],fixed=T)}
      #if(match.call()["ylab"] != "NULL()"){myylab <- gsub("()","",match.call()["ylab"],fixed=T)}
 
	plot(1/sqrt(n[ids]), ML1[ids],
         xlab=myxlab, 
		 ylab=myylab, 
         pch=20, col=dotcol[ids], ylim=myylim, xlim=myxlim,
         main=main.txt, axes=FALSE)        
	axis(2)
	#x1 <- pretty( 1/sqrt(n), 5)
      x1 <- pretty( myxlim, 5)
	axis(1, x1, round(1/x1^2))
	box()
	points(1/sqrt(n[id.raw]), ML1[id.raw],   col="black", pch=21,bg="white")
	text(1/sqrt(n[id.raw]), ML1[id.raw], "raw",  col="black",adj=-0.5, cex=0.7)
	if(length(id.start)>0){
		points(1/sqrt(n[id.start]), ML1[id.start],   col="green", pch=20)
		text(1/sqrt(n[id.start]), ML1[id.start], "start",  col="green",adj=-0.5, cex=0.7)
	}


}



print.selected.cem <- function(x, ...){
	print(x$breaks)
}






spacegraph.plot <- function(obj, objname=NULL, group="1", explore=TRUE,
                            method.col=list(cem="#0000ff75",mdm="#ff000075",
                                            psm="#00ff0075",matchit="#FF880075"),
                            plot.call=NULL,balance.metric="L1",
                            N="treated",scale.var=TRUE,...){
	if(!interactive() | !explore){
            if(match.call()["xlim"] != "NULL()" | match.call()["ylim"] != "NULL()")
              cat("cannot specify xlim or ylim with explore=F\n")
            if(match.call()["main"] != "NULL()" | match.call()["xlab"] != "NULL()" | match.call()["ylab"] != "NULL()")
              cat("cannot specify main, xlab, or ylab with explore=F\n")
		spacegraph.plot2(obj, group, method.col=method.col)
		return(NULL)	
	}
	
	
	if(class(obj) != "spacegraph")
	stop("obj must be of class `spacegraph'")

	
	haveTCL <- interactive()
	if(!capabilities("tcltk")){
		haveTCL <- FALSE	
		cat("\ntcltk support is absent")
	}
	
	if(haveTCL){
		have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
		if(have_ttk) {
			tkbutton <- ttkbutton
			tkframe <- ttkframe
			tklabel <- ttklabel
			tkradiobutton <- ttkradiobutton
		}
	}

        if(is.null(objname)){
          obj.name <- match.call()["obj"]
          obj.name <- gsub("()","",obj.name, fixed=T)
        } else {
          obj.name <- objname
        }
        
        data <- obj$data
        obj$space <- na.omit(obj$space)
       # keepme <- rownames(na.omit(obj$space[-which(names(obj$space) %in% c("ML1","Mdiff"))]))
       # obj$space <- obj$space[keepme,]  ## this solves a problem
                                         ## when some solutions have NA in the table
                                         ## but it allows ML1 and Mdiff to have
                                         ##  NAs in the table.


	g <- sprintf("G%s",group)
	n <-  obj$space[[g]]
        ## Allow the x-axis to be different
        ## Note, this is over-riding "group"
        if(N=="treated"){
          group <- "1"
	  g <- sprintf("G%s",group)
	  n <-  obj$space[[g]]
        }
        if(N=="control"){
          group <- "0"
	  g <- sprintf("G%s",group)
	  n <-  obj$space[[g]]
        }
        if(N=="all"){
          group1 <- "1"
	  g1 <- sprintf("G%s",group1)
          group0 <- "0"
	  g0 <- sprintf("G%s",group0)
	  n <-  obj$space[[g1]] + obj$space[[g0]]
        }
	ML1 <- obj$space$ML1
      Mdiff <- obj$space$Mdiff
      Mdisc <- obj$space$Mdisc
      alldist <- obj$alldist
      ## some warnings if you try to specify the wrong balance metric
      if(balance.metric=="L1" & length(na.omit(ML1))!=length(ML1)){
        stop("Missing values for L1")
      }
      if(balance.metric=="mdiff" & length(na.omit(Mdiff))!=length(Mdiff)){
        stop("Missing values for Mdiff")
      }
      if(balance.metric=="mdisc" & length(na.omit(Mdisc))!=length(Mdisc)){
        stop("Missing values for Mdiff")
      }
      ## NOTE: If Mdiff is specified, we name it "ML1" here so that the 
      ##       function will still run through.
      if(balance.metric=="mdiff"){
        ML1 <- obj$space$Mdiff
      }
      if(balance.metric=="mdisc"){
        ML1 <- obj$space$Mdisc
      }
	Relaxed <- obj$space$Relaxed
        Method <- obj$space$Method
        Call <- obj$space$Call
        ## rename the stuff in the call
        for(i in 1:length(Call)){
          Call[i] <- gsub("obj$",paste(obj.name,"$",sep=""),Call[i],fixed=T)
        }
	class <- rep("relax", length(n))
	id.raw <- which(Relaxed=="<raw>")
	id.start <- which(Relaxed=="<start>") 
	class[ id.raw ] <- "raw"
	if(length(id.start)>0)
		class[ id.start ] <- "start"
	ids <- (1:length(n))[-c(id.raw, id.start)]
	name.vars <- names(obj$coars[[1]])
	n.vars <- length(name.vars)
	tab <- data.frame(n=n, ML1=ML1, class=class)
	#main.txt <- sprintf("Total units=%d, ML1=%.3f", n[id.raw], ML1[id.raw])
      main.txt <- sprintf("Total units=%d", n[id.raw])

	
	if(length(id.start)>0){
		main.txt <- sprintf("Initial matched units=%d, ML1=%.3f", n[id.start], ML1[id.start])
         }

        ## make different colors of dots for the different methods
        dotcol <- rep("cyan",length(ids))
        ## get colors for other methods...
        approved.methods <- c("raw","cem","psm","mdm","matchit")
        needcol <- unique(Method)[unique(Method) %in% approved.methods==F]


        ## I decided that it's easier to specify my own colors
        ## UPDATE: 28 june 2011, decided that it's better to just
        ##   make them all the same color.  If the user wants, they
        ##   can make their own plot from the space that distinguishes.
        colz <- c("#CC00FF75","darkolivegreen","gray50",
                  "lavenderblush","peru")
        rcol <- colz[1:length(needcol)]
        
        ll <- length(approved.methods)-1
        if(length(needcol)>0){
          for(i in (1:length(needcol) + ll)){
            pnames <- names(method.col)
            #method.col[[i]] <- rcol[i-ll]
            method.col[[i]] <- "#CC00FF75"
            names(method.col) <- c(pnames,needcol[i-ll])
          }
        }

        for(i in 1:length(method.col)){
          dotcol[Method==names(method.col)[i]] <-  method.col[[names(method.col)[i]]]
        }
       
	#dotcol <- rgb(0.2,0.2,0.8)
	#selcol <- rgb(0.8,0.8,0.2)
        selcol <- "black"
      ## set the plotting args
      if(scale.var==F)
        n <- 1/n^2
      myxlab <- "number matched, scaled as 1/sqrt(matched)"
      myxlim <- range(1/sqrt(n))
      if(scale.var==F)
        myxlim <- rev(myxlim)
      if(balance.metric=="L1"){
        myylab <- "L1"
        myylim <- c(0,range(ML1)[2])
      }
      if(balance.metric=="mdiff"){
        myylab <- "Average Difference in Means"
        myylim <- c(0,range(ML1)[2])
      }
      if(balance.metric=="mdisc"){
        myylab <- "Average Mahalanobis Discrepancy"
        myylim <- c(0,range(ML1)[2])
      }
      

      
      ## get user-specified xlim and ylim values
      if(is.null(plot.call)){
        plot.call <- match.call()
      }
      if(plot.call["xlim"] != "NULL()"){
        xli <- gsub("()","",plot.call["xlim"],fixed=T)
        xli2 <- strsplit(xli,",")[[1]]
        lim1 <- as.numeric(gsub(" *","",gsub("c(","",xli2[1], fixed=T)))
        lim2 <- as.numeric(gsub(" *","",gsub(")","",xli2[2], fixed=T)))
        myxlim <- c(1/sqrt(lim1),1/sqrt(lim2))
      }
      if(plot.call["ylim"] != "NULL()"){
        xli <- gsub("()","",plot.call["ylim"],fixed=T)
        xli2 <- strsplit(xli,",")[[1]]
        lim1 <- as.numeric(gsub(" *","",gsub("c(","",xli2[1], fixed=T)))
        lim2 <- as.numeric(gsub(" *","",gsub(")","",xli2[2], fixed=T)))
        myylim <- c(lim1,lim2)
      }

      ## collect any user-specified plotting args that conflict:
      if(plot.call["main"] != "NULL()"){main.txt <- gsub("()","",plot.call["main"],fixed=T)}
      if(plot.call["xlab"] != "NULL()"){myxlab <- gsub("()","",plot.call["xlab"],fixed=T)}
      if(plot.call["ylab"] != "NULL()"){myylab <- gsub("()","",plot.call["ylab"],fixed=T)}
      
        
	plot(1/sqrt(n[ids]), ML1[ids],
         xlab=myxlab, 
		 ylab=myylab, 
         pch=20, col=dotcol[ids], ylim=myylim, xlim=myxlim,
         main=main.txt, axes=FALSE)        
	axis(2)
	#x1 <- pretty( 1/sqrt(n), 5)
      x1 <- pretty( myxlim, 5)
        if(scale.var==T)
	  axis(1, x1, round(1/x1^2))
        if(scale.var==F)
          axis(1, x1, x1)
	box()
	points(1/sqrt(n[id.raw]), ML1[id.raw],   col="black", pch=21,bg="white")
	text(1/sqrt(n[id.raw]), ML1[id.raw], "raw",  col="black",adj=-0.5, cex=0.7)
	if(length(id.start)>0){
		points(1/sqrt(n[id.start]), ML1[id.start],   col="green", pch=20)
		text(1/sqrt(n[id.start]), ML1[id.start], "start",  col="green",adj=-0.5, cex=0.7)
	}
	idx.sav <- id.raw
	if(length(id.start)>0)
		idx.sav <- id.start

	old.idx <- NULL
	xy <- xy.coords(1/sqrt(n), ML1)
	
 	tmp.br <- obj$coars[[2]]
	if(length(id.start)>0)
		tmp.br <- obj$coars[[id.start]]
	new.br <- tmp.br

	
	goOn <- TRUE
	
	if(haveTCL){	
		tclServiceMode(FALSE)
		tt <- tktoplevel()
		
		tkwm.title(tt,"Modify CEM solution")
		entries <- list()
		tcvars <- list()

               
            if(scale.var==T)   
		  infoText <- tclVar( sprintf("Raw Data: Total units=%d, ML1=%.3f", n[id.raw], ML1[id.raw]) )
            if(scale.var==F)   
		  infoText <- tclVar( sprintf("Raw Data: Total units=%d, ML1=%.3f", round(1/sqrt(n[id.raw]),0), ML1[id.raw]) )
                callInfo <- tclVar(Call[id.raw])
                
                if(length(id.start)>0)
                        if(scale.var==T)
                          infoText <- tclVar( sprintf("Method=%s: Matched units=%d, ML1=%.3f", Method[id.start],n[id.start], ML1[id.start]) )
                        if(scale.var==F)
                          infoText <- tclVar( sprintf("Method=%s: Matched units=%d, ML1=%.3f", Method[id.start],round(1/sqrt(n[id.start]),0), ML1[id.start]) )
                        callInfo <- tclVar(Call[id.start])
                        callText <- Call[id.start]
		tkpack(tklabel(tt, textvariable=infoText))


                ## print the call

            entry.call <- tkentry(tt, width="100", textvariable=callInfo)	 
	      tkpack( tklabel(tt, text=sprintf("Call: ") ))
		tkpack(  entry.call )

                ## Make a button so that you can edit and rerun the
                ## call
                OnCall <- function(){
			other.args <- NULL
			cat("\n... running the call ...\n")
			tmpc <- tclvalue( callInfo )
			valid <- try(eval(parse(text=tmpc)), silent = TRUE)
			if( class(valid) == "try-error"){
				cat(sprintf("\nError in call specification. Exiting.\n", tmpc) )
				break
			}
                  tmp.mat <- eval(parse(text= tclvalue(callInfo)))
                 
			tclServiceMode(TRUE)
                  if(class(tmp.mat)=="cem.match"){
                      if(balance.metric=="L1"){
			      tmp.ML1 <- L1.meas(obj$match$groups, data=data[,obj$match$vars], breaks=obj$medianCP, weights=tmp.mat$w)$L1
			    }
                      if(balance.metric=="mdiff"){
                        tmp.ML1 <- mdiff(caliperdat=data[,c(obj$match$treatment, obj$match$vars)],alldat=data[,c(obj$match$treatment, obj$match$vars)],mvars=tmp.mat$vars,tvar=obj$match$treatment, wt=tmp.mat$w)
                      }
                      if(balance.metric=="mdisc"){
                        tmp.ML1 <- mdisc2(caliperdat=data[,c(obj$match$treatment, obj$match$vars)],alldat=data[,c(obj$match$treatment, obj$match$vars)],mvars=tmp.mat$vars,tvar=obj$match$treatment, wt=tmp.mat$w)
                      }
                      tmp.n <- tmp.mat$tab["Matched", g] 
			    len.c <- length(n)
			    n[len.c+1] <<- tmp.n 	 
			    ML1[len.c+1] <<- tmp.ML1 	 
			    Relaxed[len.c+1] <<- "<new>"
			    Method[len.c+1] <<- "cem"
              	    Call[len.c+1] <<- tclvalue(callInfo)
			    obj$coars[[len.c+1]] <<- new.br	
			    xy <<- xy.coords(1/sqrt(n), ML1)
			    update.imbplot(len.c+1)	
                  }
                  if(class(tmp.mat)=="psm.match"){
                      tvar <- as.character(as.formula(gsub("()","",tmp.mat$call["formula"],fixed=T)))[2]
                      if(balance.metric=="L1"){
                        tmp.ML1 <- L1.meas(tmp.mat$match.dat[[tvar]], data=tmp.mat$match.dat[,obj$match$vars], breaks=obj$medianCP, weights=tmp.mat$match.dat[,"weights"])$L1
                      }
                      if(balance.metric=="mdiff"){
                        tmp.ML1 <- mdiff(caliperdat=tmp.mat$match.dat[,c(tvar,obj$match$vars)],alldat=obj$data[,c(tvar,obj$match$vars)],mvars=obj$match$vars,tvar=tvar, wt=tmp.mat$match.dat[,"weights"])
                      }
                      if(balance.metric=="mdisc"){
                        tmp.ML1 <- mdisc2(caliperdat=tmp.mat$match.dat[,c(tvar,obj$match$vars)],alldat=obj$data[,c(tvar,obj$match$vars)],mvars=obj$match$vars,tvar=tvar, wt=tmp.mat$match.dat[,"weights"])
                      }
                      tmp.n <- tmp.mat$N[group]
                      len.c <- length(n)
			    n[len.c+1] <<- tmp.n 	 
			    ML1[len.c+1] <<- tmp.ML1 	 
			    Relaxed[len.c+1] <<- "<new>"
			    Method[len.c+1] <<- "psm"
              	    Call[len.c+1] <<- tclvalue(callInfo)
                       ## make a filler for now
                       new.br <- c()
                       for(i in 1:length(name.vars)){
                         new.br[[i]] <- c(0,0)
                       }
                       names(new.br) <- name.vars	 
			    obj$coars[[len.c+1]] <<- new.br	
			    xy <<- xy.coords(1/sqrt(n), ML1)
			    update.imbplot(len.c+1)
                  }
                  if(class(tmp.mat)=="mdm.match"){
                      tvar <- as.character(as.formula(gsub("()","",tmp.mat$call["formula"],fixed=T)))[2]
                      if(balance.metric=="L1"){
                        tmp.ML1 <- L1.meas(tmp.mat$match.dat[[tvar]], data=tmp.mat$match.dat[,obj$match$vars], 
                                         breaks=obj$medianCP, weights=tmp.mat$match.dat[,"weights"])$L1
                      }
                      if(balance.metric=="mdiff"){
                        tmp.ML1 <- mdiff(caliperdat=tmp.mat$match.dat[,c(tvar,obj$match$vars)],alldat=obj$data[,c(tvar,obj$match$vars)],mvars=obj$match$vars,tvar=tvar, wt=tmp.mat$match.dat[,"weights"])
                      }
                      if(balance.metric=="mdisc"){
                        tmp.ML1 <- mdisc2(caliperdat=tmp.mat$match.dat[,c(tvar,obj$match$vars)],alldat=obj$data[,c(tvar,obj$match$vars)],mvars=obj$match$vars,tvar=tvar, wt=tmp.mat$match.dat[,"weights"])
                      }
                      tmp.n <- tmp.mat$N[group]
                      len.c <- length(n)
			    n[len.c+1] <<- tmp.n 	 
			    ML1[len.c+1] <<- tmp.ML1 	 
			    Relaxed[len.c+1] <<- "<new>"
			    Method[len.c+1] <<- "mdm"
              	    Call[len.c+1] <<- tclvalue(callInfo)
                       ## make a filler for now
                       new.br <- c()
                       for(i in 1:length(name.vars)){
                         new.br[[i]] <- c(0,0)
                       }
                       names(new.br) <- name.vars	 
			    obj$coars[[len.c+1]] <<- new.br	
			    xy <<- xy.coords(1/sqrt(n), ML1)
			    update.imbplot(len.c+1)             
                  }
                  if(class(tmp.mat)=="matchit"){
                      ## Note, this is not changed to allow mean differences
			    tmp.ML1 <- L1.meas(obj$match$groups, data=data[,obj$match$vars], 
                                         breaks=obj$medianCP, weights=tmp.mat$weights)$L1
			    tmp.n <- tmp.mat$nn["Matched", "Treated"]
                      len.c <- length(n)
			    n[len.c+1] <<- tmp.n 	 
			    ML1[len.c+1] <<- tmp.ML1 	 
			    Relaxed[len.c+1] <<- "<new>"
			    Method[len.c+1] <<- "matchit"
              	    Call[len.c+1] <<- tclvalue(callInfo)
                       ## make a filler for now
                       new.br <- c()
                       for(i in 1:length(name.vars)){
                         new.br[[i]] <- c(0,0)
                       }
                       names(new.br) <- name.vars	 
			    obj$coars[[len.c+1]] <<- new.br	
			    xy <<- xy.coords(1/sqrt(n), ML1)
			    update.imbplot(len.c+1)             
                  }
                }

            Call.but <-tkbutton(tt,text="   Run this Call   ",command=OnCall)
		tkpack(Call.but)

                

                ## I was attempting to make a text box so that you
                ## could see the whole call at once, but I couldn't
                ## get it to work.
#                txt <- tktext(tt, width="100")
#                #entry.text <- tkinsert(entry.call, callInfo)
#                tkinsert(txt,"end", callText)
#	        tkpack( tklabel(tt, text=sprintf("Call: ") ))
#		tkpack(  txt )

                #tkpack(tktext(tt, bg="white",font="courier",width="100"))

		n.tmp.br <- length(tmp.br)
		for(i in 1:n.tmp.br){
			tcvars[[i]] <- tclVar( deparse( round(tmp.br[[i]], 2), width.cutoff=500) )  
			entries[[i]] <- tkentry(tt, width="100", textvariable=tcvars[[i]])	    
			
			tkpack( tklabel(tt, text=sprintf("Variable: %s", names(tmp.br)[i]) ))
			tkpack(  entries[[i]] )
		}
	 
               

		tcvars[[n.tmp.br+1]] <- tclVar(  )  
		entries[[n.tmp.br+1]] <- tkentry(tt, width="100", textvariable=tcvars[[n.tmp.br+1]])	    
		
		tkpack( tklabel(tt, text="Additional CEM args:"))
		tkpack(  entries[[n.tmp.br+1]] )
                #tkpack( entries[[n.tmp.br]])
		
		
		OnOK <- function(){
			other.args <- NULL
			cat("\n... running new cem...\n")
			n.tmp.br <- length(tmp.br)
			for(i in 1:n.tmp.br){
				vv <- names(tmp.br)[i]
				tmpc <- tclvalue( tcvars[[i]] )
				new.br[[i]] <<- try(eval(parse(text=tmpc)), silent = TRUE)
				if( class(new.br[[i]]) == "try-error"){
					cat(sprintf("\nError in settings cutpoints of variable << %s >>:\n\n >> %s <<\n\n Using original ones.\n", vv, tmpc) )
					new.br[[i]] <<- tmp.br[[i]] 
				}
			}
			tmpc <- tclvalue( tcvars[[n.tmp.br+1]] )
			other.args <- try(eval(parse(text=tmpc)), silent = TRUE)
			if( class(other.args) == "try-error"){
				cat(sprintf("\nError in additional CEM arguments specification. Ignoring them.\n", tmpc) )
				other.args <- NULL 
			} else 
			 other.args <- tmpc
			tclServiceMode(FALSE)	
                  ## reconstruct the call
                  cpholder <- rep(NA,length(new.br))
                  for(ii in 1:length(new.br)){
                    cpholder[ii] <-paste(names(new.br)[ii],"=c(",paste(new.br[[ii]],collapse=", "),")",sep="")
                  }
                  breakslist <-  paste("list(",paste(cpholder,collapse=", "),")")

                   cem.call <-  sprintf("cem(obj$match$treatment, data=data[,c(obj$match$vars,obj$match$treatment) ], cutpoints=%s, %s)", breakslist, other.args)
			#if(!is.null(other.args))
			tmp.mat <- eval(parse(text=cem.call))
                  #      else
			#   tmp.mat <- cem(obj$match$treatment, data=data[,c(obj$match$vars,obj$match$treatment) ], cutpoints=new.br, eval.imbalance=FALSE)

                        
			tclServiceMode(TRUE)
			#tmp.ML1 <- L1.meas(obj$match$groups, data=data[,obj$match$vars], breaks=obj$medianCP, weights=tmp.mat$w)$L1
                  if(balance.metric=="L1"){
			  tmp.ML1 <- L1.meas(obj$match$groups, data=data[,obj$match$vars], breaks=obj$medianCP, weights=tmp.mat$w)$L1
			}
                  if(balance.metric=="mdiff"){
                    tmp.ML1 <- mdiff(caliperdat=data[,c(obj$match$treatment, obj$match$vars)],alldat=data[,c(obj$match$treatment, obj$match$vars)],mvars=tmp.mat$vars,tvar=obj$match$treatment, wt=tmp.mat$w)
                  }
                  if(balance.metric=="mdisc"){
                    tmp.ML1 <- mdisc2(caliperdat=data[,c(obj$match$treatment, obj$match$vars)],alldat=data[,c(obj$match$treatment, obj$match$vars)],mvars=tmp.mat$vars,tvar=obj$match$treatment, wt=tmp.mat$w)
                  }
                  tmp.n <- tmp.mat$tab["Matched", g] 
			len.c <- length(n)
			n[len.c+1] <<- tmp.n 	 
			ML1[len.c+1] <<- tmp.ML1 	 
			Relaxed[len.c+1] <<- "<new>"	
			Method[len.c+1] <<- "cem"
              	Call[len.c+1] <<- cem.call 
			obj$coars[[len.c+1]] <<- new.br	
			xy <<- xy.coords(1/sqrt(n), ML1)
			update.imbplot(len.c+1)	 
		}
		
		
		OK.but <-tkbutton(tt,text="   Run CEM with these coarsening   ",command=OnOK)
		
		
		tkpack(OK.but) 
		tclServiceMode(TRUE)

##            ## Make a button to update the spacegraph obj.
## I COULDN'T GET THIS TO WORK EASILY
##		UPDATE <- function(){
##			other.args <- NULL
##			cat("\n... updating spacegraph object ...\n")
##                
##                  update <- paste(obj.name,"$space <- cbind(n,ML1,Relaxed,Method,Call)",sep="")
##                  print(update)
##                  eval(parse(text=update))
##                  cat("\n... update complete ...\n")
##		}
##		
##		
##		UPDATE.but <-tkbutton(tt,text="   Save modified spacegraph info   ",command=UPDATE)
##		
##		tkpack(UPDATE.but) 
##		tclServiceMode(TRUE)
        
	}
	
	firstime <- TRUE
	update.imbplot <- function(mT){
		if(length(mT)==0)
		return(FALSE)
		
		if(scale.var==T)
		  idx <- which(n>=n[mT])
            if(scale.var==F)
              idx <- which(n<=n[mT])
		idx2 <- which(ML1[idx]<= ML1[mT])
		idx <- idx[idx2]

                ## my addition -- will it screw things up?
                idxx <- idx[-which(idx %in% mT)]

             
		
		id.new <- which(Relaxed =="<new>")


                ## white-out
                text(1/sqrt(n[id.raw]), ML1[id.raw], "\u2585", col="white", cex=5, adj=0)
                points(1/sqrt(n[ids]), ML1[ids],   col="white", pch=19,cex=2)
                if(length(id.new))
                  points(1/sqrt(n[id.new]), ML1[id.new],   col="white", pch=19,cex=2)
                box()
                points(1/sqrt(n[ids]), ML1[ids],   col=dotcol[ids], pch=20)
                #points(1/sqrt(n[idx.sav]), ML1[idx.sav],   col="white", pch=19,cex=2.2)
		#points(1/sqrt(n[idx.sav]), ML1[idx.sav],   col=dotcol[ids], pch=20)
		points(1/sqrt(n[idxx]), ML1[idxx],   col=selcol, pch=1,cex=1)
               	#points(1/sqrt(n[mT]), ML1[mT],   col="white", pch=3)
                points(1/sqrt(n[mT]), ML1[mT],   col=dotcol[mT], pch=20)
		points(1/sqrt(n[mT]), ML1[mT],   col="black", pch=1,
                       cex=1.5)
		points(1/sqrt(n[id.raw]), ML1[id.raw],   col="black",
                       pch=21, bg="white")
                 
                text(1/sqrt(n[id.raw]), ML1[id.raw], "raw",
                     col="black",adj=-0.5, cex=0.7)
                
		points(1/sqrt(n[mT]), ML1[mT],   col=dotcol[mT], pch=20)
		if(length(id.start)>0){
			points(1/sqrt(n[id.start]), ML1[id.start],   col="green", pch=20)
			text(1/sqrt(n[id.start]), ML1[id.start], "start",  col="green",adj=-0.5, cex=0.7)
		}
		if(length(id.new)){
               newdotcol <- rep("#00000050",length(id.new))
               newdotcol <- unlist(method.col[Method[id.new]])
		   points(1/sqrt(n[id.new]), ML1[id.new],   col=newdotcol, pch=20)
            }
             points(1/sqrt(n[mT]), ML1[mT],   col="black", pch=1, cex=1.5)
		id.bad <- which(idx == id.raw)
		if(length(id.bad>0))
		 idx <- idx[-id.bad]
		if(length(idx)>0){
 		  y <- lapply(obj$coars[idx], function(a) unlist(lapply(a, length)))
		  x <- matrix(unlist(y), length(y), length(y[[1]]), byrow=TRUE) 
		
		  colnames(x) <- names(y[[1]])
		  tmp <- as.data.frame(x)
				
		  tmp2 <- data.frame(tmp, ML1=ML1[idx])	
		  rownames(tmp2) <- idx
		
		 if(haveTCL & !is.null(obj$coars[[mT]])){
			ttt <- obj$coars[[mT]]

 			for(i in 1:length(tmp.br)){
				tclvalue( tcvars[[i]] )  <-  deparse( round(ttt[[i]], 2), width.cutoff=500)
  	 		}
                        
                        mycall <- Call[mT]


                        tclvalue(callInfo) <- mycall
                        
			       if(scale.var==T)
                           tclvalue(infoText) <- sprintf("Method=%s, Matched units=%d, ML1=%.3f",toupper(Method[mT]), n[mT], ML1[mT]) 
                         if(scale.var==F)
                           tclvalue(infoText) <- sprintf("Method=%s, Matched units=%d, ML1=%.3f",toupper(Method[mT]), round(1/sqrt(n[mT]),0), ML1[mT])
                        tkfocus(tt)	
		 }
		}

		
		old.idx <<- list(breaks = obj$coars[[mT]], n=n[mT], ML1=ML1[mT], medianCP=obj$medianCP)
	
		
		idx.sav <<- idx
		
		return(TRUE)
		
	}
	
	old.idx <- NULL
	update.imbplot( 2 )
	
	while(goOn){
		xy <<- xy.coords(1/sqrt(n), ML1)
		goOn <- update.imbplot(identify(xy, n=1,plot=FALSE)) 
	}
	
    if(haveTCL)
	 tkdestroy(tt)
	
	
	class(old.idx) <- "selected.cem"
  	old.idx
	
}
############################################################

############################################################
## random.matchit

random.matchit <- function( treatment=NULL, 
                           data = NULL, drop=NULL, calipermax=.25, 
                           R=10, method = "unknown", progress=T,...){
  holder <- c()
  callholder <- c()
  if(progress==T & R>1){pb <- txtProgressBar(min = 1, max = R, initial = 1, style = 3)}
  for(r in 1:R){
    if(progress==T & R>1){setTxtProgressBar(pb, r)}
    mvars <- names(data)[names(data)!=treatment]
    ## generate all of the interactions
    allIntObj <- enumerateInteractions(mvars)
    ## generate a formula
    myform <- genForm(allIntObj,treatment=treatment,data=data, poly=1)
    ## draw the caliper
    rcaliper <- round(runif(0,calipermax,n=1),5)
    ## do the pscore matching
    ## check if distance and method were specified in the call
    dist <- NULL
    meth <- NULL
    if(match.call()["distance"] != "NULL()"){dist <- gsub("()","",match.call()["distance"],fixed=T)}
    if(match.call()["method"] != "NULL()"){meth <- gsub("()","",match.call()["method"],fixed=T)}

    if(!is.null(dist)){
      dist.text <- paste(", distance = \"",dist,"\"", sep="")
    } else { 
      dist.text <- ""
    }
    if(!is.null(meth)){
      meth.text <- paste(", method = \"",meth,"\"", sep="")
    } else { 
      meth.text <- ""
    }
    dist.meth <- paste(dist.text, meth.text,sep="")
    mcall <- sprintf("matchit(formula = %s, data=data, caliper = %s %s)",myform,rcaliper,dist.meth)
    m.out <- eval(parse(text = mcall))

    out <- c()
    out$id <- names(m.out$w)
    out$weight <- m.out$w
    out$method <- rep(method,length(m.out$w))
    out <- as.data.frame(out, stringsAsFactors=F)
    holder[[r]] <- out
    callholder[[r]] <- mcall
  }
  if(progress==T & R>1){close(pb)}
  return(list(match.list=holder,calls=callholder))
}
############################################################

############################################################
## A difference in means function

mdiff <- function(caliperdat,alldat,mvars,tvar, wt=NULL){
  if(is.null(wt)){
    wt <- rep(1,nrow(caliperdat))
  }
  covs <- subset(caliperdat, select=mvars)
  covs <- apply(covs,MARGIN=2,as.numeric)
  covs2 <- subset(alldat, select=mvars)
  covs2 <- apply(covs2,MARGIN=2,as.numeric)
  ## I ran into an issue with datasets that only have 2 observations
  ## Can there be a weighted avg with only one obs? I code it as NO.
  if(sum(caliperdat[[tvar]]==1) == 1){
    t.mean <- covs[caliperdat[[tvar]]==1,]
  } else {    
    t.mean <- apply(covs[caliperdat[[tvar]]==1,], MARGIN=2, weighted.mean,w=wt[caliperdat [[tvar]]==1], na.rm=T)
  }
  if(sum(caliperdat[[tvar]]==0) == 1){
    c.mean <- covs[caliperdat[[tvar]]==0,]
  } else {
    c.mean <- apply(covs[caliperdat[[tvar]]==0,], MARGIN=2, weighted.mean,w=wt[caliperdat [[tvar]]==0], na.rm=T)
  }
  t.sd <- apply(covs2[alldat[[tvar]]==1,], MARGIN=2, sd, na.rm=T)
  return(mean(abs((t.mean-c.mean)/t.sd), na.rm=T))
}

#############################################################

#############################################################
## A mahalanobis discrepancy function
## adist is a matrix of mah distances produced by mdistfun()

#mdisc <- function(caliperdat, adist=alldist, wt = NULL){
mdisc <- function(caliperdat, adist, wt = NULL){
  if(is.null(wt)){wt <- rep(1,nrow(caliperdat))}
  matched.names <- rownames(caliperdat)[wt>0]
  mdis <- adist[rownames(adist)[rownames(adist) %in% matched.names],
          colnames(adist)[colnames(adist) %in% matched.names]]
  if(is.matrix(mdis)){
    try(minT <- apply(mdis,MARGIN=1,min), silent=T)
    ## If the vector is too big for apply
    if(exists("minT")==F){
      minT <- rep(NA,nrow(mdis))
      for(jj in 1:nrow(mdis)){minT[jj] <- min(mdis[jj,])}
    }
    ## 
    try(minC <- apply(mdis,MARGIN=2,min), silent=T)
    if(exists("minC")==F){
      minC <- rep(NA,ncol(mdis))
      for(jj in 1:ncol(mdis)){minC[jj] <- min(mdis[,jj])}
    }
  } else { 
    minT <- min(mdis)
    minC <- min(mdis)
  }    
  return(mean(c(minT,minC)))
}


## a second version that doesn't need the distances calculated outside
mdisc2 <- function(caliperdat, alldat,mvars,tvar, wt = NULL){
  ff <- as.formula(paste(tvar,"~",paste(mvars,collapse="+")))
  try(adist <- mdistfun(ff,alldat), silent=T)
  if(exists("adist")==F){
    adist <- mdistfun2(ff,alldat)
  }
  if(is.null(wt)){wt <- rep(1,nrow(caliperdat))}
  matched.names <- rownames(caliperdat)[wt>0]
  mdis <- adist[rownames(adist)[rownames(adist) %in% matched.names],
          colnames(adist)[colnames(adist) %in% matched.names]]
  if(is.matrix(mdis)){
    minT <- apply(mdis,MARGIN=1,min)
    minC <- apply(mdis,MARGIN=2,min)
  } else { 
    minT <- min(mdis)
    minC <- min(mdis)
  }     
  return(mean(c(minT,minC)))
}

#############################################################

#############################################################
## A function for combining spacegraph objects for plotting

combine.spacegraphs <- function(x,y){
      if (class(x) != "spacegraph") 
        stop("x must be of class `spacegraph'")
      if (class(y) != "spacegraph") 
        stop("y must be of class `spacegraph'")
      out <- x
      new <- y$space[-which(y$space[,"Relaxed"] == "<raw>"),]
      out$space <- rbind(x$space,new)
      out$coars <- c(x$coars,y$coars)
      return(out)
}

#############################################################
