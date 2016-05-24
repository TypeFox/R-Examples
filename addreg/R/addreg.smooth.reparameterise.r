addreg.smooth.reparameterise <- function(coefficients, interpret, allref, knots, design.knots, design.param,
                                            subset, na.action) {
	coefs.new <- coefficients
	coefs.new.int <- coefs.new[1]

	smthterms <- sapply(interpret$smooth.spec,"[[","term")
    smthnames <- NULL
	for(smth in smthterms) {
		smthlabel <- interpret$smooth.spec[[smth]]$termlabel
		which.smth <- which(substr(names(coefficients),1,nchar(smthlabel))==smthlabel)
		coefs.smth <- coefficients[which.smth]
        smthnames <- c(smthnames, names(coefs.smth))
		smthtype <- class(interpret$smooth.spec[[smth]])
		if(smthtype == "Iso.smooth") {
			coefs.smth.new <- coefs.smth
			names.smth.new <- names(coefs.smth)
		} else if(smthtype == "B.smooth") {
			ref <- allref$allref[[smth]][[as.numeric(design.param[smth])]]
			num.knots <- length(knots[[smth]])
			if(length(ref) == 1) {
				coefs.smth.temp <- append(coefs.smth, 0, after = ref - 1)
				coefs.new.int <- coefs.new.int + coefs.smth.temp[1]
				coefs.smth.new <- (coefs.smth.temp - coefs.smth.temp[1])[-1]
			} else {
				coefs.smth.temp <- c(0, cumsum(coefs.smth))
				coefs.smth.ord <- coefs.smth.temp[order(ref)]
				coefs.new.int <- coefs.new.int + coefs.smth.ord[1]
				coefs.smth.new <- (coefs.smth.ord - coefs.smth.ord[1])[-1]
			}
			names.smth.new <- paste(smthlabel, 2:(num.knots-3), sep = "")
		}
		coefs.new[which.smth] <- coefs.smth.new
		names(coefs.new)[which.smth] <- names.smth.new
	}
	
	coefs.new[1] <- coefs.new.int
    
    dummy.design <- rep(1,length(design.param))
    names(dummy.design) <- names(design.param)
    dummy.allref <- allref
    for(smth in names(allref$allref))
        dummy.allref$allref[[smth]][[1]] <- 1
    modelspec <- addreg.smooth.design(interpret, dummy.allref, design.knots, dummy.design)
    data.new <- modelspec$data
    dummy.frame.call <- call("model.frame", formula = eval(modelspec$formula), data = as.name("data.new"))
    dummy.frame.call$drop.unused.levels <- TRUE
    if (!missing(na.action)) dummy.frame.call$na.action <- na.action
    if (!missing(subset)) dummy.frame.call$subset <- subset
    dummy.frame <- eval(dummy.frame.call)
    dummy.terms <- attr(dummy.frame, "terms")
    
    design <- model.matrix(dummy.terms, dummy.frame)
	
	list(coefs = coefs.new, mf = dummy.frame, design = design, mt = dummy.terms,
            smoothnames = smthnames)
}