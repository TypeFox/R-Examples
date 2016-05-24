summary.gsym.point <-
function(object, ...) {
	methods <- object[object$methods]
	levels.cat <- if(is.null(object$levels.cat)) {"Global"} else {object$levels.cat}
	conf.level <- ifelse(is.null(object$call$confidence.level), 0.95, object$call$confidence.level)
	p.results <- names(methods[[1]][[1]][["optimal.result"]])
  	ci.legend <- paste(paste(conf.level*100, "% CI", sep = ""), c("lower limit", "upper limit"))
	res <- vector("list", length(levels.cat))
	for (i in 1:length(levels.cat)) {
		for(j in 1:length(methods)) {
			row.names <- c(p.results)
			col.names <- c("Estimate", ci.legend)
			res[[i]][[j]] <- vector("list", length(methods[[j]][[i]][["optimal.result"]][["cutoff"]][[1]]))
			if (length(methods[[j]][[i]][["optimal.result"]][["cutoff"]][[1]]) != 0) {
				for(k in 1:length(methods[[j]][[i]][["optimal.result"]][["cutoff"]][[1]])){
					m <- matrix(ncol = 3, nrow = length(row.names), dimnames = list(row.names,col.names ))
					m[1,] <- methods[[j]][[i]][["optimal.result"]][[1]][k,]
					for (l in 2:length(p.results)) {
						m[l,] <- methods[[j]][[i]][["optimal.result"]][[l]][k,]
					}

					res[[i]][[j]][[k]] <- m
				}
			}
		}
  	res[[i]][[length(methods)+1]] <- round(methods[[1]][[i]][["AUC"]][[1]], 3)
      names(res[[i]]) <- c(names(methods), "AUC")
	}

	names(res) <- levels.cat
	object$p.table <- res
	class(object) <- "summary.gsym.point"
	object
}
