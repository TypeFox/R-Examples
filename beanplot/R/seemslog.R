`seemslog` <-
function(groups) {
    if (min(unlist(groups)) <= 0) {
        logevidence = 1
		nologevidence = 1
    }
    else {
        nologevidence = sum(sapply(groups, function(x) {
            ifelse(length(x) > 2, ifelse(length(x) > 5000, 
            log(max(0, shapiro.test(x[seq(1,5000,floor((length(x)-1)/5000+1))])$p.value)),
            log(max(0, shapiro.test(x)$p.value))), 
                0)
        }))	
        logevidence = sum(sapply(groups, function(x) {
            ifelse(length(x) > 2, ifelse(length(x) > 5000, 
			log(max(0, shapiro.test(log(x[seq(1,5000,floor((length(x)-1)/5000+1))]))$p.value)),
            log(max(0, shapiro.test(log(x))$p.value))), 
                0)
        }))
    }
    (logevidence > nologevidence + 1)
}

