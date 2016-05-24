snk.test <-
function (object, term, among = NULL, within = NULL) 
{
	tm <- estimates(object)$tm
	n.mean <- tm[1:nrow(tm)-1, , drop = FALSE]
	for (j in 1:ncol(tm)) {
		for (i in 1:nrow(tm)) {
			if (tm[i,j] == 0) n.mean[i,j] = 1
			}
	}
	n.mean <- rowProds(n.mean)
	f.versus <- estimates(object)$f.versus
	anova.res <- gad(object)
	se.num <- anova.res$Mean[which(rownames(anova.res) == f.versus[match(term, rownames(f.versus))])]
	se.den <- n.mean[which(term == rownames(f.versus))]
	se <- sqrt(se.num/se.den)
	if (length(se) == 0) stop ("F-ratio versus = No test \nStandard error for this set of means cannot be calculated!")
	Df <- anova.res$Df[which(rownames(anova.res) == f.versus[match(term, rownames(f.versus))])]
	if (length(among) == 0) {
		set.means <- tapply(object$model[,1], subset(object$model, select = term), mean)
		a <- length(set.means)
		rank.means <- sort(set.means)
		m.means <- matrix(NA, ncol = a, nrow = a-1)
		m.final <- matrix("", ncol = a, nrow = a-1)
		i <- 0; g <- a:2
		for (i in 2:a-1) {
			for (j in 1:i) {
				m.means[i,j] <- rank.means[a-i+j]-rank.means[j]
			}
		}
		Q <- m.means/se
		pvalue <- Q
		for (i in 1:nrow(m.means)) {
			for (j in 1:ncol(m.means)) {
				pvalue[i,j] <- ptukey(Q[i,j], nmeans = g[i], df = Df, lower.tail = FALSE)
			}
		}
		for (i in 2:a-1) {
			for (j in 1:i) {
				if (pvalue[i,j] <= 0.001) m.final[i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "***")
				else if (pvalue[i,j] <= 0.01) m.final[i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "**")
				else if (pvalue[i,j] < 0.05) m.final[i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "*")
				else if (pvalue[i,j] > 0.05) m.final[i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "ns")
			}
		}
		result <- data.frame(m.final, stringsAsFactors = FALSE)
		result <- rbind(as.character(""), result)
	    result <- rbind(round(rank.means, 4), result)
		result <- rbind(as.character(1:ncol(m.final)), result)
		colnames(result) <- names(rank.means)
		rownames(result) <- c("Rank order:", "Ranked means:", "Comparisons:", 1:nrow(m.means))	
	} else {
		within.split <- strsplit(within, ":")
		within.int <- interaction(object$model[, within.split[[1]]], lex.order = TRUE)
		new.model <- cbind(object$model, within.int)
		set <- split(new.model, f = factor(within.int, levels = unique(within.int)))
		set.means <- set
		for (i in 1:length(set)) {
			set.means[[i]] <- tapply(set[[i]][,1], subset(set[[i]], select = among), mean)
		}
		a <- length(set.means[[1]])
		rank.means <- set.means
		for (i in 1:length(set.means)) {
			rank.means[[i]] <- sort(set.means[[i]])
		}
		m.means <- matrix(NA, ncol = a, nrow = a-1)
		m.final <- matrix("", ncol = a, nrow = a-1)
		i <- 0; g <- a:2
		l.means <- rank.means
		for (i in 1:length(l.means)) {
			l.means[[i]] <- m.means
		}
		for (i in 2:a-1) {
			for (j in 1:i) {
				for (k in 1:length(l.means)) {
					l.means[[k]][i,j] <- as.vector(rank.means[[k]])[a-i+j]-as.vector(rank.means[[k]])[j]
				}
			}	
		}
		Q <- l.means
		for (i in 1:nrow(m.means)) {
			for (j in 1:ncol(m.means)) {
				for (k in 1:length(Q)) {
					Q[[k]][i,j] <- l.means[[k]][i,j]/se
				}
			}	
		}
		pvalue <- Q
		for (i in 1:nrow(m.means)) {
			for (j in 1:ncol(m.means)) {
				for (k in 1:length(Q)) {
					pvalue[[k]][i,j] <- ptukey(Q[[k]][i,j], nmeans = g[i], df = Df, lower.tail = FALSE)
				}
			}	
		}
		l.final <- l.means
		for (i in 1:length(l.final)) {
			l.final[[i]] <- m.final
		}
		for (i in 2:a-1) {
			for (j in 1:i) {
				for (k in 1:length(l.final)) {
					if (pvalue[[k]][i,j] <= 0.001) l.final[[k]][i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "***")
					else if (pvalue[[k]][i,j] <= 0.01) l.final[[k]][i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "**")
					else if (pvalue[[k]][i,j] < 0.05) l.final[[k]][i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "*")
					else if (pvalue[[k]][i,j] > 0.05) l.final[[k]][i,j] <- paste(paste(c(a-i+j, j), collapse = "-"), "ns")
				}
			}	
		}
		result <- l.final
		for (i in 1:length(l.final)) {
			result[[i]] <- data.frame(l.final[[i]], stringsAsFactors = FALSE)
			result[[i]] <- rbind(as.character(""), result[[i]])
			result[[i]] <- rbind(round(rank.means[[i]], 4), result[[i]])
			result[[i]] <- rbind(1:ncol(m.means), result[[i]])
			colnames(result[[i]]) <- names(rank.means[[i]])
			rownames(result[[i]]) <- c("Rank order:", "Ranked means:", "Comparisons:", 1:nrow(m.means))
		}
	}
	cat("\nStudent-Newman-Keuls test for:", term, "\n")
	cat("\nStandard error =", round(se, 4))
	cat("\nDf =", Df, "\n")
	if (length(among) == 0) {
		print(result, right = FALSE)
	} else {
		cat("\nPairwise comparisons among levels of:", among, "\nwithin each level of:", within, "\n")
		for (i in 1:length(result)) {
			cat("\nLevel:", names(result)[[i]], "\n")
			print(result[[i]], right = FALSE)
		}
	}
	cat("---")
	cat("\nSignif. codes: <0.001 '***' <0.01 '**' <0.05 '*' >0.05 'ns'\n")
}

