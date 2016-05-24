summary.optimal.cutpoints <-
function(object, ...) {
	opt.criterion.methods <- c("MCT","CB", "MaxSpSe", "MaxProdSpSe", "ROC01", "SpEqualSe", "Youden", "MaxEfficiency", "Minimax", "MaxDOR", "MaxKappa", "PROC01", "NPVEqualPPV", "MaxNPVPPV", "MaxSumNPVPPV", "MaxProdNPVPPV", "MinPvalue", "PrevalenceMatching")
	methods <- object[object$methods]
	ci.fit <- ifelse(is.null(object$call$ci.fit), FALSE, object$call$ci.fit)
	levels.cat <- if(is.null(object$levels.cat)) {"Global"} else {object$levels.cat}
	conf.level <- ifelse(is.null(object$call$conf.level), 0.95, object$call$conf.level)
	p.results <- names(methods[[1]][[1]][["optimal.cutoff"]])
	ci.legend <- paste(paste(conf.level*100, "% CI", sep = ""), c("lower limit", "upper limit"))
	res <- vector("list", length(levels.cat))
	for (i in 1:length(levels.cat)) {
		for(j in 1:length(methods)) {
			aux.criterion <- names(methods)[j] %in% opt.criterion.methods 
			row.names <- if(aux.criterion) {c(p.results,"Optimal criterion")} else {c(p.results)}
			col.names <- if(ci.fit) {
				c("Estimate", ci.legend)
			} else {
				"Estimate"
			}
			res[[i]][[j]] <- vector("list", length(methods[[j]][[i]][["optimal.cutoff"]][[1]]))
			if(length(methods[[j]][[i]][["optimal.cutoff"]][[1]]) != 0) {
				for(k in 1:length(methods[[j]][[i]][["optimal.cutoff"]][[1]])){
					m <- matrix(ncol = ifelse(ci.fit, 3, 1), nrow = length(row.names), dimnames = list(row.names,col.names ))
					m[1,1] <- methods[[j]][[i]][["optimal.cutoff"]][[1]][[k]]		
					for (l in 2:length(p.results)) {
						m[l,] <- methods[[j]][[i]][["optimal.cutoff"]][[l]][k,]
					}
					# Auxiliar criterion
					if(aux.criterion)
						m[length(p.results) + 1,1] <- methods[[j]][[i]][["optimal.criterion"]]				
					res[[i]][[j]][[k]] <- m
				}
			}				
		}
		res[[i]][[length(methods)+1]] <- paste(round(methods[[1]][[i]][["measures.acc"]][["AUC"]][1], 3), " (", round(methods[[1]][[i]][["measures.acc"]][["AUC"]][2], 3),"",", ", round(methods[[1]][[i]][["measures.acc"]][["AUC"]][3], 3),")", sep = "")
		names(res[[i]]) <- c(names(methods), "AUC_CI")
	}
	names(res) <- levels.cat
	object$p.table <- res
	class(object) <- "summary.optimal.cutpoints"
	object
}
