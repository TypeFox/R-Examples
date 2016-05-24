predictmeans <- function (model, modelterm, pairwise=FALSE, atvar=NULL, adj="none", Df=NULL, 
                          level=0.05, covariate=NULL, trans = NULL, responsen=NULL, count=FALSE, 
						  plotord=NULL, plottitle=NULL, mplot=TRUE, barplot=FALSE, pplot=TRUE, 
						  bkplot=TRUE, plot=TRUE, jitterv=0, basesz=12, prtnum=TRUE, newwd=TRUE, 
						  permlist=NULL)
{  
  options(scipen=6)
  slevel <- level
  if (class(model)[1]=="aovlist") stop("Plese use model 'lme' instead of 'aov'!")
  if (class(model)[1]=="glm") {
    trans <- model$family$linkinv
    if (model$family$family %in% c("poisson", "quasipoisson")) count=TRUE
  }
  if (class(model)[1] == "glmerMod") {
    trans <- slot(model, "resp")$family$linkinv
#    if (slot(model, "resp")$family$family %in% c("poisson", "quasipoisson", "binomial", "quasibinomial")) count=TRUE
  }
  vars <- unlist(strsplit(modelterm, "\\:"))
  mdf <- model.frame(model)
  if(any(!is.element(vars, names(mdf)[sapply(mdf,is.factor)])))
    stop(sQuote(vars), "all must be factor(s)!")
  
  if (length(vars)==1) atvar <- NULL
  if (!is.null(permlist) && !permlist%in%c("NULL", "")) {pairwise <- TRUE; if (adj=="tukey") stop("The p-value can't be adjusted by Tukey methd!")}                        # option checking
  if (!is.null(atvar) && !atvar%in%c("NULL", "")) pairwise <- TRUE
  if (adj != "none") pairwise <- TRUE
  if (!is.null(plotord) && !plotord%in%c("NULL", "")) plot <- mplot <- TRUE

  ctr.matrix <- Kmatrix(model, modelterm, covariate, prtnum)
  KK <- ctr.matrix$K
  label <- ctr.matrix$fctnames
  rownames(label) <- rownames(KK)
  response <- ctr.matrix$response
  mp <- mymodelparm(model) 
  n.table <- table(mdf[, vars, drop = FALSE])
  ndf <- data.frame(n.table)       ## To obtain info from model

  K <- KK[, mp$estimable, drop = FALSE]          # To match coef names
  if (any(ndf$Freq==0)) {
    rnTrt <- do.call("paste", c(ndf[, vars, drop=FALSE], sep=":"))
    rnTrt <- rnTrt[ndf$Freq!=0]      # To delete any missing level in factor
    K <- K[rnTrt,]
    label <- label[rnTrt,]
  }
  pm <- K %*% mp$coef
  vcovm <- mp$vcov
  # ses <- sqrt(diag(K %*% tcrossprod(vcovm, K)))
  ses <- as.numeric(apply(K, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(mp$vcov, y))}))
  mt <- data.frame(pm, ses, label)
  names(mt)[-1:-2] <- vars
  bkmt <- mt  # for back transformed
  mean.table <- round(xtabs(pm ~ ., mt[, c("pm", vars)], drop.unused.levels = TRUE), 4)
  se.table <- round(xtabs(ses ~ ., mt[, c("ses", vars)], drop.unused.levels = TRUE), 5)
  mean.table[!(n.table)] <- NA
  se.table[!(n.table)] <- NA
  if (length(vars) > 1) {
    varsnlevel <- numeric(0)
    for (i in vars) varsnlevel[i] <- nlevels(mdf[, i])
    tbvars <- names(sort(varsnlevel, decreasing = TRUE))
    mean.table <- ftable(mean.table, row.vars =tbvars[1], col.var=tbvars[-1])
    se.table <- ftable(se.table, row.vars =tbvars[1], col.var=tbvars[-1])
  }
  if (length(na.omit(unique(se.table)))==1) {
    se.table <- min(se.table, na.rm=TRUE)
    names(se.table) <- "All means have the same Stder"
  }
  
  nK <- nrow(K)          # To setup various matrix and row, col names
  rnK <- rownames(K)
  varn1 <- varn2 <- rep(0, nK * (nK - 1)/2)
  
  if (nK == 1) {
    SED.out <- NA
    LSD <- NA
  }else {
    kindx <- 1:nK
    CM <-  matrix(0, nrow=nK * (nK - 1)/2, ncol=nK)
    t <- 1
    for (i in 2:nK) {                      # To construct pairwise comparison K matrix by col order
      for (j in 1:(i-1)) {
        CM[t, ] <- (kindx == j) - (kindx == i)
        varn1[t] <- rnK[i]
        varn2[t] <- rnK[j]
        t <- t+1
      }
    }
    rK <- CM%*%K                    # calculate stats
    cm <- rK%*%mp$coef
    #vcov.contr <- rK %*% tcrossprod(vcovm, rK)
    #dses <- sqrt(diag(vcov.contr))
    dses <- as.numeric(apply(rK, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(vcovm, y))}))
	  if (adj == "bonferroni") level <- level/length(dses)
	  
    SED.out <- c(Max.SED = max(dses), Min.SED = min(dses), Aveg.SED = mean(dses))
    dses.df <- data.frame(matrix(unlist(strsplit(varn1, "\\:")), byrow=T, nrow=length(varn1)),
                            matrix(unlist(strsplit(varn2, "\\:")), byrow=T, nrow=length(varn2)), dses)
							
    if (all(length(vars) > 1, SED.out[1]!=SED.out[2])) {
      dses.m <- matrix(0, nrow=3, ncol=length(vars))
      colnames(dses.m) <- vars
      rownames(dses.m) <- c("Aveg.SED", "Min.SED", "Max.SED")
      for (i in 1:length(vars)) {
        varsindx <- as.character(dses.df[,i])==as.character(dses.df[, length(vars)+i])
        if(any(varsindx))  # in case of A/B
        dses.m[,i] <- summary.default(dses.df[varsindx,"dses"])[c(4,1,6)]  else {dses.m[,i] <- NA; atvar <- NULL}
  	  }
      attr(SED.out, "For the Same Level of Factor") <- dses.m
    }

  if (is.null(permlist) || permlist%in%c("NULL", "")) {
	if (length(Df) == 0) {
      if (class(model)[1] == "lme") {
        Df <- terms(model$fixDF)[modelterm]
      }else if (class(model)[1] == "lmerMod") {
		termlabel <- attr(terms(model),"term.labels")
		for (i in vars) termlabel <- termlabel[grep(i, termlabel)]
		termlabel <- paste(termlabel, collapse="-")
		model.b <- update( model, as.formula(paste(".~. -", termlabel)))
		Df <- getKR(KRmodcomp(model, model.b), "ddf")
      }else Df <- mp$df

	  if (Df==0) stop("You need provide Df for this model!")
    }
      LSD <- round(qt(1 - level/2, df = Df) * SED.out, 5)
      names(LSD) <- c("Max.LSD", "Min.LSD", "Aveg.LSD")
      attr(LSD, "For the Same Level of Factor") <- NULL
  #    if (length(vars) > 1) {
  #      rownames(dses.m) <- c("Aveg.LSD", "Min.LSD", "Max.LSD")
  #      attr(LSD, "For the Same Level of Factor") <- qt(1 - level/2, df = Df) * dses.m
  #    } # end of if LSD
      attr(LSD, "Significant level") <- slevel
      attr(LSD, "Degree of freedom") <- round(Df, 2)
    }else{
      LSD <- round(2 * SED.out[1:3], 5)
      names(LSD) <- c("Max.LSD", "Min.LSD", "Aveg.LSD")
      attr(LSD, "Note") <- "This is a approximate LSD which is 2*SED."          
    }

    if (pairwise) {
      tvm <- t.p.valuem <- LSDm <- Diffm <- matrix(0, ncol = nK, nrow = nK)
      rownames(tvm) <- colnames(tvm) <- rownames(t.p.valuem) <- colnames(t.p.valuem) <- rownames(LSDm) <- colnames(LSDm) <- rnK
      t.v <- cm/dses
	  if (all(is.null(permlist) || permlist%in%c("NULL", ""), adj=="tukey")) p.tukey <- ptukey(sqrt(2)*abs(t.v), nK, Df, lower.tail=FALSE)
      tvm[upper.tri(tvm)] <- t.v
      if (is.null(permlist) || permlist%in%c("NULL", "")) {
        t.p.values <- 2 * pt(-abs(t.v), Df)
      }else{      
      	nsim <- length(permlist[[1]])
        tValue <- function(x, rK){
          cm <- rK%*%x$coef
          vcovm <- x$vcov
          #vcov.contr <- rK %*% tcrossprod(vcovm, rK)
          #ses <- sqrt(diag(vcov.contr))
          ses <- as.numeric(apply(rK, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(vcovm, y))}))
          t.v <- cm/ses
          return(t.v)	
        }
        if (nK==2)	t.p.values <-  (sum(sapply(permlist[[1]], function(x) round(abs(tValue(x, rK)),6) >= round(abs(t.v), 6)))+1)/(nsim+1) else	
          t.p.values <-  (rowSums(sapply(permlist[[1]], function(x) round(abs(tValue(x, rK)),6) >= round(abs(t.v), 6)))+1)/(nsim+1)
      } # end of if (is.null(permlist))
            
      if (is.null(atvar) || atvar%in%c("NULL", "")) {
		    if (adj=="tukey") t.p.valuem[upper.tri(t.p.valuem)] <- p.tukey else
          t.p.valuem[upper.tri(t.p.valuem)] <- p.adjust(t.p.values, adj)
        t.p.valuep <- t.p.valuem    # for plot
        t.p.valuem <- t(t.p.valuem) + tvm
        names(t.p.valuem) <- NULL
        if (is.null(permlist) || permlist%in%c("NULL", "")) {
          attr(t.p.valuem, "Degree of freedom") <- Df
          attr(t.p.valuem, "Note") <- paste("The matrix has t-value above the diagonal, p-value (adjusted by",
                                            sQuote(adj), "method) below the diagonal")
		      LSDm.up <- qt(1-level/2, df = Df)*dses
				  LSDm[upper.tri(LSDm)] <- LSDm.up
				  Diffm[upper.tri(Diffm)] <- cm
          LSDm <- t(LSDm)+Diffm
          names(LSDm) <- NULL
          attr(LSDm,"Significant level") <- slevel
          attr(LSDm,"Degree of freedom") <- Df
          attr(LSDm,"Note") <- paste("LSDs matrix has mean differences (row-col) above the diagonal, LSDs (adjusted by",
                                            sQuote(adj), "method) below the diagonal")
        }else{
          attr(t.p.valuem, "Note") <- paste("The matrix has t-value above the diagonal, and", sQuote(nsim), "times permutation p-value (adjusted by",
                                      sQuote(adj), "method) below the diagonal")       
        } # end of if (is.null(permlist)) 
		multGrp <- data.frame(multcompLetters(t.p.valuem, Letters=LETTERS, threshold=level))
		names(multGrp) <- "Group" 
		attr(t.p.valuem, paste("Letter-based representation of pairwise comparisons at significant level", sQuote(level))) <- multGrp
        if (all(nrow(t.p.valuep) > 3, pplot, plot)) {
          mtitle <- plottitle
          if (is.null(plottitle) || plottitle%in%c("NULL", "")) mtitle <- paste("Level Plot of p-value (adjusted by", sQuote(adj), "method)\n for Pairwise Comparison")		
          PMplot(t(t.p.valuep), level=slevel, legendx=0.69, mtitle=mtitle, newwd=newwd)
		} # end of if (all(nrow(t.p.valuep) > 3, pplot, plot)) 
      }else{
        # atvar <- vars[which(vars %in% atvar)] # ensure a proper order for atvar
        dses.df$tvalue <- t.v
        dses.df$pvalue <- t.p.values 
        atvar.df <- dses.df  
        for (i in which(vars%in%atvar)) { # To find rows relating atvar
          atvar.df <- atvar.df[as.character(atvar.df[, i]) == as.character(atvar.df[, length(vars) + i]), ]
        }
		    if (adj=="tukey") atvar.df$adj.pvalue <- ptukey(sqrt(2)*abs(atvar.df$tvalue), nK, Df, lower.tail=FALSE) else
          atvar.df$adj.pvalue <- p.adjust(atvar.df$pvalue, adj)
        names(atvar.df)[1:length(vars)] <- vars
        for (i in vars) {                     # To ensure factor vars  have the same level as original
          atvar.df[,i] <- factor(atvar.df[,i], levels=levels(mdf[, i]))
        }
        atvar.df <- atvar.df[do.call(order, atvar.df[, atvar, drop=FALSE]),]  # sort data by atvar order

        rnK.df <- as.data.frame(matrix(unlist(strsplit(rnK, "\\:")), byrow=T, nrow=nK)) # To find the suitable names for pm
        colnames(rnK.df) <- vars
        for (i in vars) {       # To ensure factor vars  have the same level as original
          rnK.df[,i] <- factor(rnK.df[,i], levels=levels(mdf[, i]))
        }
        rnK.df  <- rnK.df [do.call(order, rnK.df[, atvar, drop=FALSE]),]
        resvar <- vars[!vars%in%atvar]   # The rest of vars rather than at var
        rnK.df <- rnK.df[, c(atvar, resvar)]   # To ensure the right matrix name later
        atvar.levels <- unique(do.call("paste", c(rnK.df[, atvar, drop=FALSE], sep=" : ")))
        resvar.levels <- unique(do.call("paste", c(rnK.df[, resvar, drop=FALSE], sep=" : "))) 
        rcnplotm <- do.call("paste", c(rnK.df[, , drop=FALSE], sep=" : ")) # row col names of image plot 
                  
        listlength <- length(atvar.levels)
        pmlist <- vector("list", listlength)
        indexlength <- nrow(atvar.df)/listlength
        nrow.pm <- length(resvar.levels)

        for (i in 1:listlength) {           # extract pvalues for each factor level
          atvar.pm <- matrix(0, nrow=nrow.pm, ncol=nrow.pm)
          atvar.pm[upper.tri(atvar.pm)] <- atvar.df$adj.pvalue[(indexlength*i-(indexlength-1)):(indexlength*i)]
          rownames(atvar.pm) <- colnames(atvar.pm) <- resvar.levels
          pmlist[[i]] <- t(atvar.pm)
        } 
        
        if (all(nrow.pm > 2, pplot, plot)) { 
    		  mtitle <- plottitle
          if (is.null(plottitle) || plottitle%in%c("NULL", ""))  mtitle <- paste("Adjusted p-value (by", sQuote(adj), 
            "method)\n for Pairwise Comparison at Each Level of",paste(sQuote(atvar), collapse =" and "))
          PMplot(pmlist, level=slevel, xylabel=rcnplotm, legendx=0.69, mtitle=mtitle, newwd=newwd)
        }    
      } # if (is.null(atvar))
    }# end of if(pairwise)

    if (plot) {
      if (length(vars) > 3)
        cat("\n", "There is no plot for more than three-way interaction! \n\n")
	    plotmt <- na.omit(mt)	  
      yMin <- min(plotmt[, "pm"])
      yMax <- max(plotmt[, "pm"])
      offSet <- 0.25 * (yMax - yMin)
      up <- yMin + LSD[3]
      lsdBar <- cbind(plotmt[, vars, drop=FALSE], up=up, yMin=yMin)
      limits <- aes(ymax = (pm + ses)*(pm > 0) + pmin(pm + ses, 0)*(pm <= 0), ymin=(pm - ses)*(pm < 0) + pmax(pm - ses, 0)*(pm >= 0))
	  
      if (length(vars) == 1) {
		  if (mplot) {
  			if (newwd) dev.new()
  			mtitle <- plottitle
  			if (is.null(plottitle) || plottitle%in%c("NULL", "")) mtitle <- paste("Predicted means for \"", vars, "\" with Aveg.LSD (", slevel * 100, "%) Bar", sep="")
    		p <- qplot(eval(parse(text = vars)), pm, data=plotmt, xlab = paste("\n", vars, sep=""),
    		  ylab = paste(response, "\n", sep=""), main = paste(mtitle, "\n", sep=""), 
          ylim = c(yMin - offSet, max(yMax + offSet, yMin + LSD[3] + offSet)),
          xlim = c("Aveg.LSD", levels(plotmt[, vars])))
        p <- p +geom_point(colour="red")+
          geom_errorbar(data=lsdBar, aes(ymax=up, ymin=yMin, x="Aveg.LSD"), width=0.15, size=0.8, colour="blue") + 
          theme_bw(basesz)
        print(p)		  
	    } # end of if mplot	   
      if (barplot) {
        if (newwd) dev.new()
	      mtitle <- plottitle
        if (is.null(plottitle) || plottitle%in%c("NULL", "")) mtitle <- paste("Predicted means for \"", modelterm, "\" with Stder Bars", sep = "")
        p <- qplot(eval(parse(text = vars)), pm, data=plotmt, ylim = c(min(0, (min(pm) - max(ses)) * 1.2),
                   max(0, (max(pm) + max(ses)) * 1.2)), xlab=paste("\n", vars, sep=""), ylab=paste(response, "\n", sep=""), 
                   main=paste(mtitle, "\n", sep=""))
        dodge <- position_dodge(width=0.9)
        p <- p + geom_bar(position=dodge, stat="identity", fill="papayawhip", colour="darkgreen") +
          geom_errorbar(limits, position=dodge, width=0.25, colour="blue")+theme_bw(basesz)+
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())
		print(p)
      }
      }
      if (length(vars) == 2) {
        if (is.null(plotord) || plotord%in%c("NULL", "")) plotord <- 1:2
          fact1 <- (vars[plotord])[1]
          fact2 <- (vars[plotord])[2]	  
  	    if (mplot) {     
    			if (newwd) dev.new()
    	#		if (is.null(plotord)) plotord <- 1:2
        #  fact1 <- (vars[plotord])[1]
        #  fact2 <- (vars[plotord])[2]
        #  nlvel1 <- nlevels(plotmt[, fact1])
        #  nlvel2 <- nlevels(plotmt[, fact2])
     			mtitle <- plottitle
    			if (is.null(plottitle) || plottitle%in%c("NULL", ""))  mtitle <- paste("Predicted means for \"", fact1, "\" by \"", fact2, "\" with Aveg.LSD (", 
            slevel * 100, "%) Bar", sep = "")
          plotmt[, fact1] <- factor(plotmt[, fact1], levels = c("Aveg.LSD", levels(plotmt[, fact1])))
          p <- qplot(eval(parse(text = fact1)), pm, data=plotmt, group=eval(parse(text = fact2)), 
            col=eval(parse(text = fact2)), main=paste(mtitle, "\n", sep=""), xlab = paste("\n", fact1, sep=""),
      		  ylab = paste(response, "\n", sep=""), xlim = levels(plotmt[, fact1]),
            ylim = c(yMin - offSet, max(yMax + offSet, yMin + LSD[3] + offSet)))
          p <- p+ geom_line(aes(linetype=eval(parse(text = fact2)), col=eval(parse(text = fact2))), size=0.96)+
                   geom_errorbar(data=lsdBar, aes(ymax=up, ymin=yMin, x="Aveg.LSD"), width=0.15, size=0.8, colour="blue")+
                  guides(linetype = guide_legend(title = fact2))+
                  guides(col = guide_legend(title = fact2))+
                  theme_bw(basesz)  
		 print(p)
#                  theme(legend.position = c(0.12-0.01*(max(nlvel1-5, 0)), 0.88-0.02*(max(nlvel2-3,0))), 
#                    legend.background = element_rect(colour = "black")) )
    		} # end if mplot
    		if (barplot) {
          if (newwd) dev.new()
          mtitle <- plottitle
    			if (is.null(plottitle) || plottitle%in%c("NULL", ""))  mtitle <- paste("Predicted means for \"", modelterm, "\" with Stder Bars", sep = "")
		      dodge <- position_dodge(width=0.9)
          p <- qplot(eval(parse(text = fact1)), pm, data=plotmt, geom="bar", stat="identity", 
            position=dodge, fill=eval(parse(text = fact2)), ylim = c(min(0, (min(pm) - max(ses)) * 1.1),
            max(0, (max(pm) + max(ses)) * 1.1)), xlab = paste("\n", fact1, sep=""),
      		  ylab = paste(response, "\n", sep=""), main = mtitle)
          p <- p + geom_errorbar(limits, position=dodge, width=0.25, colour="blue")+theme_bw(basesz)+
            scale_fill_brewer()+
            theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+            
            theme(legend.position = "top")+guides(fill = guide_legend(title = fact2))
			print(p)
        }
      }
      
      if (all(length(vars) == 3, mplot)) {
        if (newwd) dev.new()
		mtitle <- plottitle
  		  if (is.null(plotord) || plotord%in%c("NULL", "")) plotord <- 1:3
  		  fact1 <- (vars[plotord])[1]
  		  fact2 <- (vars[plotord])[2]
  		  fact3 <- (vars[plotord])[3]
  			plotmt[, fact1] <- factor(plotmt[, fact1], levels = c("Aveg.LSD", levels(plotmt[, fact1])))
  			if (is.null(plottitle) || plottitle%in%c("NULL", "")) mtitle <- paste("Predicted means for '", fact1, "' by '", fact2, "' for each '",
                           fact3, "'\n with Aveg.LSD (", slevel * 100, "%) Bar\n", sep = "")
  			p <- qplot(eval(parse(text = fact1)), pm,  data=plotmt, main=mtitle,
          xlab = paste("\n", fact1, sep=""), ylab = paste(response, "\n", sep=""),
          xlim = levels(plotmt[, fact1]), ylim=c(yMin-0.5*offSet, max(yMax+0.5*offSet, yMin+LSD[3]+0.5*offSet)),
          group=factor(eval(parse(text = fact2))), col=factor(eval(parse(text = fact2)))) +
          geom_errorbar(data=lsdBar, aes(ymax=up, ymin=yMin, x="Aveg.LSD"), width=0.15, size=0.8, colour="blue")+
  		    geom_line(aes(linetype=eval(parse(text = fact2)), col=eval(parse(text = fact2))), size=0.8)+
  			  facet_grid(eval(parse(text = paste("~",fact3, sep=""))))+
  			  guides(group = guide_legend(fact2))+
  			  guides(linetype = guide_legend(fact2))+
  			  guides(col = guide_legend(fact2))+
  			  theme_bw(basesz)
  			print(p)
      }

    }
  }
  
  if (!is.null(trans)) {
    Mean <- Trt <- NULL
    bkmt$Mean <- trans(bkmt$pm)
    if (is.null(permlist) || permlist%in%c("NULL", "")) {    
      bkmt$LL <- trans(bkmt$pm - qt(1 - slevel/2, df = Df) * bkmt$ses)
      bkmt$UL <- trans(bkmt$pm + qt(1 - slevel/2, df = Df) * bkmt$ses)
      }else{
      bkmt$LL <- trans(bkmt$pm - 2 * bkmt$ses)
      bkmt$UL <- trans(bkmt$pm + 2 * bkmt$ses)
    }
    bkmt$pm <- bkmt$ses <- NULL
    nc <- ncol(bkmt)
    names(bkmt)[c(nc - 1, nc)] <- c(paste("LL of ", (1 - slevel) * 100, "% CI", sep = ""),
      paste("UL of ", (1 - slevel) * 100, "% CI", sep = ""))
    bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], 4)
    if (count) bkmt[, (nc - 2):nc] <- round(bkmt[, (nc - 2):nc], 0)

    if (plot && bkplot) {
      if (response %in% names(mdf)) {    ## Transformed y before modelling
        if (class(mdf[, response])=="factor"){
          bky <- as.numeric(mdf[, response])-1
        }else{
          if (class(model)[1]=="glm" | class(model)[1]=="glmerMod") {
            bky <- mdf[, response]
            if (!is.null(dim(mdf[, response]))) bky <- mdf[, response][,1]/rowSums(mdf[, response])
          }else{
            bky <- trans(mdf[, response])
          }# end of if glm or glmer
        }# end of if factor
      }else{       ## Transformed y within modelling
        nresponse <- regmatches(response, regexec("\\(([^<]+)\\)", response))[[1]][2]
        if (!nresponse %in% names(mdf)) {
          # ques <- paste("R thinks the back transformed response is", sQuote(nresponse), "but the right one is: ")
          # nresponse <- readline(ques)
		  if (is.null(responsen) || responsen%in%c("NULL", ""))  stop("Please provide suitable name for response variable using option 'responsen'!")
          nresponse <- responsen
        }
        bky <- mdf[, nresponse]
      }

      Trtn <- do.call("paste", c(mdf[, vars, drop=FALSE], sep=":"))
  	  newdata2 <- data.frame(bky=bky, Trtn=Trtn)
  	  bkmt$Trt <- do.call("paste", c(bkmt[, vars, drop=FALSE], sep=":"))
      xMax <- max(max(bkmt[, nc], na.rm=TRUE), bky, na.rm=TRUE)
      xMin <- min(min(bkmt[, nc - 1], na.rm=TRUE), bky, na.rm=TRUE)
      xoffSet <- 0.15 * (xMax - xMin)
      mtitle <- plottitle
		  if (is.null(plottitle) || plottitle%in%c("NULL", "")) mtitle <- paste("Back Transformed Means with ", (1 - slevel) * 100, "% CIs\n for '",
        modelterm, "'", "\n", sep = "")
  	  if (newwd) dev.new()
  	  limits <- aes(xmax = bkmt$"UL of 95% CI", xmin=bkmt$"LL of 95% CI")
  	  xlimv <- c(xMin - xoffSet, xMax + xoffSet)
        p <- qplot(Mean, Trt, main = mtitle, ylab = "", xlab = "", xlim = xlimv, data=bkmt)+
          geom_point(colour="red") + geom_errorbarh(limits, height=0.2, size=0.8, colour="red") +
  			  scale_y_discrete(limits = rev(unique(bkmt$Trt)))+
          geom_point(aes(x=bky, y=Trtn), shape=1, position = position_jitter(width = jitterv, height = jitterv), colour="blue", alpha=0.6, data=newdata2)+
  			  theme_bw(basesz)
        print(p)
    }
    rownames(bkmt) <- NULL
    bkmt$Trt <- NULL
    if (pairwise) {
      if (is.null(atvar) || atvar%in%c("NULL", "")) {
        if (all(is.null(permlist) || permlist%in%c("NULL", ""), adj %in% c("none", "bonferroni"))) {
          return(list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                    "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise LSDs"=round(LSDm,5),
                    "Pairwise p-values" = round(t.p.valuem,4), "Back Transformed Means" = bkmt))
        }else{
          return(list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                    "Standard Error of Differences" = SED.out, LSD = LSD,
                    "Pairwise p-values" = round(t.p.valuem,4), "Back Transformed Means" = bkmt))
        }
      }else{
          outputlist <- vector("list", listlength+5)
          outputlist[[1]] <- mean.table
          outputlist[[2]] <- se.table
          outputlist[[3]] <- SED.out
          outputlist[[4]] <- LSD
          
          for (i in 5: (listlength+4)) {
      		  outtab <- round(as.table(pmlist[[i-4]]), 4)
      		  outtab[outtab < 0.0001] <- "0.0001"
      		  outtab[col(outtab)==row(outtab)] <- 1.0000
      		  outtab[upper.tri(outtab)] <-""
			  outtab <- as.table(cbind(outtab, Group=multcompLetters(pmlist[[i-4]], Letters=LETTERS, threshold=level)))		      
      		  outputlist[[i]] <- outtab
    		  }
			  outputlist[[listlength+5]] <- bkmt
    		  if (is.null(permlist) || permlist%in%c("NULL", "")) {
            names(outputlist) <- c(c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
               "LSD"), paste("Pairwise comparison p-value (adjusted by", sQuote(adj), "method)","\n", "for variable", 
                paste(sQuote(resvar),collapse =" and "), "at level <", atvar.levels, "> of", paste(sQuote(atvar), collapse =" and ")), "Back Transformed Means")                   
          }else{
            names(outputlist) <- c(c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
               "Approximated LSD"), paste("Pairwise", sQuote(nsim), "times permuted p-value (adjusted by", sQuote(adj), "method)","\n", "for variable", 
                paste(sQuote(resvar),collapse =" and "), "at level <", atvar.levels, "> of", paste(sQuote(atvar), collapse =" and ")), "Back Transformed Means with an Approximated 95% CIs (Mean +/- 2*SE)")                            
          }          
          return(outputlist)        
        }
    }
    else {
      return(list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                  "Standard Error of Differences" = SED.out, LSD = LSD,
                  "Back Transformed Means" = bkmt))
    }
  }else {
    if (pairwise) {
      if (is.null(atvar) || atvar%in%c("NULL", "")) {
        if (all(is.null(permlist) || permlist%in%c("NULL", ""), adj %in% c("none", "bonferroni"))) {
          return(list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                      "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise LSDs"=round(LSDm,5), "Pairwise p-value" = round(t.p.valuem, 4)))
        }else{
          return(list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                      "Standard Error of Differences" = SED.out, LSD = LSD, "Pairwise p-value" = round(t.p.valuem, 4)))
        }        
       }else{
        outputlist <- vector("list", listlength+4)
        outputlist[[1]] <- mean.table
        outputlist[[2]] <- se.table
        outputlist[[3]] <- SED.out
        outputlist[[4]] <- LSD
    		for (i in 5: (listlength+4)) {
    		  outtab <- round(as.table(pmlist[[i-4]]), 4)
    		  outtab[outtab < 0.0001] <- "0.0001"
    		  outtab[col(outtab)==row(outtab)] <- 1.0000
    		  outtab[upper.tri(outtab)] <-""
			  outtab <- as.table(cbind(outtab, Group=multcompLetters(pmlist[[i-4]], Letters=LETTERS, threshold=level)))
			  outputlist[[i]] <- outtab
    		}
    		if (is.null(permlist) || permlist%in%c("NULL", "")) {
          names(outputlist) <- c(c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
             "LSD"), paste("Pairwise comparison p-value (adjusted by", sQuote(adj), "method)","\n", "for variable", paste(sQuote(resvar), 
             collapse =" and "), "at level <", atvar.levels, "> of", paste(sQuote(atvar), collapse =" and ")))                   
        }else{
          names(outputlist) <- c(c("Predicted Means", "Standard Error of Means", "Standard Error of Differences",
             "Approximated LSD"), paste("Pairwise", sQuote(nsim), "times permuted p-value (adjusted by", sQuote(adj), "method)","\n", "for variable", paste(sQuote(resvar), 
             collapse =" and "), "at level <", atvar.levels, "> of", paste(sQuote(atvar), collapse =" and ")))                   
        }
        return(outputlist)        
        }
    }
    else {
      return(list("Predicted Means" = mean.table, "Standard Error of Means" = se.table,
                  "Standard Error of Differences" = SED.out, LSD = LSD))
    }
  }
}
