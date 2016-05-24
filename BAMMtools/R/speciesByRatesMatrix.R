speciesByRatesMatrix = function(ephy, nslices, index = NULL, spex = "s") {
	if (!spex %in% c('s', 'e', 'netdiv')) {
		stop("arg spex must be 's', 'e' or 'netdiv'.")
	}
	seq.nod <- .Call("seq_root2tip", ephy$edge, length(ephy$tip.label), ephy$Nnode, PACKAGE = "BAMMtools");
	if (nslices <= 100) {
		tvec <- (seq(0, 1, 0.01)+0.005) * max(ephy$end);
		tvec <- tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		ephy <- dtRates(ephy, 0.01, index, tmat = TRUE);
	}
	else if (nslices > 100 && nslices <= 500) {
		tvec <- (seq(0, 1, 0.002)+0.001) * max(ephy$end);
		tvec <- tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		ephy <- dtRates(ephy, 0.002, index, tmat = TRUE);
	}
	else if (nslices > 500 && nslices <= 1000) {
		tvec <- (seq(0, 1, 0.001)+0.0005) * max(ephy$end);
		tvec <- tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		ephy <- dtRates(ephy, 0.001, index, tmat = TRUE);
	}
	else {
		stop("Max slices (1000) exceeded.  Choose a smaller number of slices");
	}
	ret <- lapply(seq.nod, function(x) {
		path = which(ephy$dtrates$tmat[,1] %in% x);
		ids = sapply(tvec[-length(tvec)], function(y) which(ephy$dtrates$tmat[path,2] <= y & ephy$dtrates$tmat[path,3] > y));
		if (is.list(ids))
			ids = unlist(ids);
		if (ephy$type == "trait") {
			rts = ephy$dtrates$rates[path][ids];
		}
		else {
			if (tolower(spex) == "s") {
				rts = ephy$dtrates$rates[[1]][path][ids];
			}
			else if (tolower(spex) == "e") {
				rts = ephy$dtrates$rates[[2]][path][ids];
			}
			else if (tolower(spex) == "netdiv") {
				rts = ephy$dtrates$rates[[1]][path][ids] - ephy$dtrates$rates[[2]][path][ids];
			}
		}
		if (length(rts) < (length(tvec)-1))
			rts = c(rts, rep(NA, length(tvec)-1-length(rts)));
		rts;
	});
	ret <- do.call(rbind, ret);
	rownames(ret) <- ephy$tip.label;
	return(list(times = tvec[-length(tvec)],rates = ret));	
}
