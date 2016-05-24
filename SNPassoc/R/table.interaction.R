`table.interaction` <-
function(var, dep, adj = NULL, int, num.status, level)
{       
	# taula.int(Datos$XRCC1.81, Datos$grupo, Datos[,c("sexo", "rcal.dia")], Datos$n.edad)
	# taula.int(Datos$XRCC1.81, Datos$grupo, NULL, Datos$n.edad)

	if (num.status==0) #Categorical response variable
	{
	    var <- as.factor(var)
	    dep <- as.factor(dep)

	    if (is.null(adj))
	    {
	       
	        m.t  <- glm(dep~ as.numeric(var) + int, family = binomial)
	
	        subset <- 1:length(var)%in%as.numeric(rownames(m.t$model));
	
	        m.b   <- glm(dep~ var + int,         subset = subset, family = binomial)
	        m.int <- glm(dep~ var/int,         subset = subset, family = binomial)
	        m.t.int <- glm(dep~ as.numeric(var) * int, subset = subset, family = binomial)
	
	    }
	    else
	    {
	        m.t  <- glm(dep~. + as.numeric(var) + int, family = binomial, data=adj)
	
	        subset <- 1:length(var)%in%as.numeric(rownames(m.t$model));
	
	        m.b   <- glm(dep~. + var + int,         subset = subset, family = binomial, data=adj)
	        m.int <- glm(dep~. + var/int,         subset = subset, family = binomial, data=adj)
	        m.t.int <- glm(dep~. + as.numeric(var) * int, subset = subset, family = binomial, data=adj)
	        
	    }
	       
		var.int <- factor(paste(levels(var)[var], levels(int)[int]), levels = outer(levels(var), levels(int), paste),
						  exclude = c(paste(levels(var), ""), paste("", levels(int)), paste(" ")))
	
		ta <- table(var.int[subset], dep[subset])

		# Matriu de coeficients i cov
	
		mat.coef <- merge(m.int$coef, summary(m.int)$coef, by=0, all.x=TRUE, sort=FALSE)
		nom.pos <- data.frame(names(m.int$coef), ordre=1:length(m.int$coef))
		mat.ordre <- merge(nom.pos, mat.coef, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
		mat.ordre <- mat.ordre[order(mat.ordre$ordre),]
	
		a <- as.matrix(mat.ordre[,c("Estimate")])
		se <- as.matrix(mat.ordre[,c("Std. Error")])
		mat <- cbind(a, se)
		selec <- dim(mat)[1.] - (length(levels(int)) - 1.) * length(levels(var))
		o <- (selec + 1.):dim(mat)[1.]
		k <- matrix(nrow = length(levels(var)), ncol = 3.)
		k[, 1.] <- 1.
		taula <- cbind(exp(mat[o, 1.]), exp(mat[o, 1.] - 1.96 * mat[o, 2.]), exp(mat[o, 1.] + 1.96 * mat[o, 2.]))
		taula[taula > 999.] <- NA
		ktaula <- rbind(k, round(taula, 2.))

		ktaula <- cbind(ta, ktaula)
	
		i <- 1;
		j <- 1;
		step <- length(levels(var));
		taula.int <- NULL;
		while (i <= nrow(ktaula))
	    {
			aux <- ktaula[i:(i+step-1),];
			colnames(aux)[3] <- levels(int)[j];
			taula.int <- cbind(taula.int, aux);
			i <- i + step;
			j <- j + 1;
		}

		#Check if interaction pvalues are NA
		t1 <- anova(m.b, m.int, test = "Chi")
                pval <- t1[2, grep("^P.*Chi",names(t1))]
		if (is.na(pval))
		{
			pval <- "NA"
		}
		else
		{
			pval <- format.pval(pval)
		}
		
		t2 <- anova(m.t, m.t.int, test = "Chi")
                pval.trend <- t2[2, grep("^P.*Chi",names(t2))]
		if (is.na(pval.trend))
		{
			pval.trend <- "NA"
		}
		else
		{
			pval.trend <- format.pval(pval.trend)
		}
		rownames(taula.int) <- levels(var);
		list(table=taula.int,pval=pval,trend=pval.trend)
	}
	else #Continuous response variable
	{
	    var <- as.factor(var)

	    if (is.null(adj))
	    {
	        m.t  <- glm(dep~as.numeric(var) + int, family = gaussian)

	        subset <- 1:length(var)%in%as.numeric(rownames(m.t$model));

	        m.b   <-   glm(dep~ var + int,             subset = subset, family = gaussian)
	        m.int <-   glm(dep~ var/int,               subset = subset, family = gaussian)
	        m.t.int <- glm(dep~ as.numeric(var) * int, subset = subset, family = gaussian)
	    }
	    else
	    {
	        m.t  <- glm(dep~. + as.numeric(var) + int, family = gaussian, data=adj)

	        subset <- 1:length(var)%in%as.numeric(rownames(m.t$model));

	        m.b <-     glm(dep~. + var + int,             subset = subset, family = gaussian, data=adj)
	        m.int <-   glm(dep~. + var/int,               subset = subset, family = gaussian, data=adj)
	        m.t.int <- glm(dep~. + as.numeric(var) * int, subset = subset, family = gaussian, data=adj)
	    }
		var.int <- factor(paste(levels(var)[var], levels(int)[int]), levels = outer(levels(var), levels(int), paste),
						  exclude = c(paste(levels(var), ""), paste("", levels(int)), paste(" ")))

		# Matriu de coeficients i cov

		mat.coef <- merge(m.int$coef, summary(m.int)$coef, by=0, all.x=TRUE, sort=FALSE)
		nom.pos <- data.frame(names(m.int$coef), ordre=1:length(m.int$coef))
		mat.ordre <- merge(nom.pos, mat.coef, by.x=1, by.y=1, all.x=TRUE, sort=FALSE)
		mat.ordre <- mat.ordre[order(mat.ordre$ordre),]
	
		a <- as.matrix(mat.ordre[,c("Estimate")])
		se <- as.matrix(mat.ordre[,c("Std. Error")])
		mat <- cbind(dif=a, lo=a-(1.96*se), up=a+(1.96*se))
		selec <- dim(mat)[1] - (length(levels(int)) - 1.) * length(levels(var))
		o <- (selec + 1):dim(mat)[1]
		mat <- mat[o,];
		
		i <- 1;
		while (i <= length(levels(var)))
		{
			mat <- rbind(c(0,NA,NA),mat);
			i <- i + 1;
		}
        
	    res <- cbind(Table.mean.se(var.int, dep, subset)$tp, mat);

	    i <- 1;
	    j <- 1;
	    step <- length(levels(var));
	    taula.int <- NULL;
	    while (i <= nrow(res))
	    {
	        aux <- res[i:(i+step-1),];
	        colnames(aux)[3] <- levels(int)[j];
	        taula.int <- cbind(taula.int, aux);
	        i <- i + step;
	        j <- j + 1;
	    }

		pval <- anova(m.b, m.int, test = "F")$"Pr(>F)"[2];
		if (is.na(pval))
		{
			pval <- "NA";
		}
		else
		{
			pval <- format.pval(pval);
		}
		
		pval.trend <- anova(m.t, m.t.int, test = "F")$"Pr(>F)"[2];
		if (is.na(pval.trend))
		{
			pval.trend <- "NA";
		}
		else
		{
			pval.trend <- format.pval(pval.trend);
		}
	    rownames(taula.int) <- levels(var);
		list(table=taula.int,pval=pval,trend=pval.trend);
	}
}

