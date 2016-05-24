print_info<-FALSE
odds <- function(object, ...)
{
	if(is.numeric(object) && all(object < 1.001) && all(object > -0.001))
		object / (1 - object)
	else
	{ message("unsupported object"); stop("unsupported object")}
}

divergence.Bernoulli <- function(P1, P2, base = exp(1), nsample = numeric(0), independent = !is_nothing(nsample), na.rm = logical(0), verbose = FALSE) # see expRelativeEntropy of odd.s, divergence of MDL.s, and relative.entropy of verisimilitude.s
{
	stopifnot(length(base) == 1)
	stopifnot(length(P1) == length(P2))
	assert.are(list(P1, P2), "numeric")
	stopifnot(is_prob(P1) && is_prob(P2))
	lg <- function(vec, ...)
	{
		vec <- as(vec, "numeric")
		logb(vec, ..., base = base)
	}
	if(!independent && is_nothing(nsample) && is_nothing(na.rm))
	{
		if(verbose)
		{
			message("no independence assumption")
			# browser()
		}
		term <- function(p1, p2)
		{
			p1 <- as(p1, "numeric")
			ifelse(p1 == 0, 0, p1 * (lg(p1) - lg(p2)))
		}
		di <- try(term(P1, P2) + term(1 - P1, 1 - P2))
		if(is_err(di))
		{ message("bad divergence!")}#; browser()
		if(any(di < 0))
		{ message("negative divergence!")}#; browser()
		if(verbose) message("returning divergence")
		di
	}
	else if(independent)
	{
		if(verbose) message("independent case for Bernoulli")
		if(is_nothing(na.rm))
			na.rm <- FALSE
		summary.stat <- function(vec, FUN)
		{
			sca <- FUN(vec, na.rm = na.rm)
			if(sca < 0)
			{
				mess <- paste("divergence of ", sca, " truncated at 0", sep = "")
				if(verbose) message(mess)
				warning(mess)
				sca <- 0
			}
			nScalar(sca)
		}
		if(is_nothing(nsample))
		{
			if(verbose) message("making recursive call")
			divergences <- divergence.Bernoulli(P1 = P1, P2 = P2, base = base, nsample = nsample, independent = FALSE, na.rm = logical(0), verbose = verbose)
			if(any(divergences < 0))
			{ message("negative divergence?")}#; browser()
			if(verbose) message("adding up")
			summary.stat(divergences, FUN = sum) # chain rule of relative entropy
		}
		else if(!is_nothing(nsample) && nsample >= 1 && base == exp(1)) # begin Monte Carlo
		{
			warning("Monte Carlo is probably not needed; try independent = TRUE instead of specifying nsample")
			get.lik <- function(p, success)
			{
				stopifnot(is.logical(success) && length(success) == length(p))
				lik <- sum(lg(ifelse(success, p, 1 - p)))
	#			if(!is.finite(lik))
	#			{ message("bad lik:"); print(stats(lik))}
				lik
			}
			get.success <- function(p)
			{
				ran <- runif(n = length(p))
				stopifnot(length(ran) == length(p) && is_prob(ran) && is_prob(p))
				ran <= p
			}
			lik.ratio <- sapply(1:nsample, function(i)
			{
				success <- get.success(P1)
				get.lik(p = P1, success = success) - get.lik(p = P2, success = success)
			})
			if(any(is.na(lik.ratio)))
			{ message("bad lik.ratio"); print(stats(lik.ratio))}#; browser()
			summary.stat(lik.ratio, FUN = mean)
		} # end Monte Carlo
		else
			stop("bad nsample, independent case")
	} # end if(independent)
	else
		stop("bad nsample")
}
#classes----------

setClass("probs.to.combine", representation("list", method.name = "character"))
setValidity("probs.to.combine", function(object)
{
	len.ok <- length(object@method.name) == length(object)
	prob.ok <- all(sapply(object, is_prob))
	class.ok <- all(sapply(object, is, class2 = "Numeric"))
	nam.ok <- all(sapply(object, function(elem){sameNames(elem, object)}))
	oks <- c(len = len.ok, prob = prob.ok, class = class.ok, nam = nam.ok)
	ok <- all(oks)
	if(!ok)
	{ printInvalid(object);  stop("bad probs.to.combine")}
	ok
})
setClass("plot_probs.to.combine", representation(prob="probs.to.combine", comb.prob = "Numeric"))
setValidity("plot_probs.to.combine", function(object)
{
	prob.ok <- is_nothing(object@comb.prob) || is_prob(object@comb.prob)
	nam.ok <- is_nothing(object@comb.prob) || identical(names(object@prob),names(object@ comb.prob))#sameNames(object, object@comb.prob)
	oks <- c(nam.ok = nam.ok, prob.ok = prob.ok)
	ok <- all(oks)
	if(!ok)
	{ printInvalid(object); stop("bad plot_probs.to.combine")} 
	ok
})
setClass("combination.prob", representation("Numeric", weight = "Numeric",names="character"))#added names by marta April 2013
setValidity("combination.prob", function(object)
{
	prob.ok <- is_prob(object) && is_prob(object@weight)
	names.ok<-length(object@.Data)==length(object@names)
	oks <- c(prob.ok = prob.ok,names.ok=names.ok)#added names by marta April 2013
	ok <- all(oks)
	if(!ok)
	{ printInvalid(object); stop("bad combination.prob")}
	ok
})





#
setMethod("[", signature(x = "combination.prob", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	x0 <- x
	x@.Data <- as(x, "Numeric")[i]
	nam0 <- names(x)
	if(is.character(nam0))
	{
		names(nam0) <- nam0
		nam <- nam0[i]
		stopifnot(length(nam) == length(names(x)))
		names(x@.Data) <- nam
		names(x) <- nam
		nam.ok <- is.character(names(x)) && all(names(x) == if(is.null(i)) nam else i)
		prob.ok <- is_prob(x)
		ok <- nam.ok && prob.ok
		if(!ok){ stop("bad [")}#marta July2014 message
	}
	x
})



setMethod("[", signature(x = "plot_probs.to.combine", i = "ANY", j = "missing"), function(x, i, j, exact) # not working?
{
#	method.name <- x@method.name[i]
	ptc <- as(x, "probs.to.combine")[i]
	stopifnot(validObject(ptc))
	x@.Data <- ptc
	x@comb.prob <- x@comb.prob[i]
	stopifnot(validObject(x))
	x
})
setMethod("[", signature(x = "plot_probs.to.combine", i = "ANY", j = "NULL"), function(x, i, j, exact) # not working?
{
	ptc <- as(x, "probs.to.combine")[i, j]
	stopifnot(validObject(ptc))
	x@.Data <- ptc
	x@method.name <- ptc@method.name # method.name
	ok <- try(validObject(x))
	if(is_err(ok) || !ok){ stop("bad x")}#marta July2014 message
	x
})
setMethod("[", signature(x = "probs.to.combine", i = "ANY", j = "missing"), function(x, i, j, drop)
{
	for(k in 1:length(x))
		x[[k]] <- x[[k]][i]
	assert.is(x, "probs.to.combine")
	stopifnot(validObject(x))
	x
})
setMethod("[", signature(x = "probs.to.combine", i = "ANY", j = "NULL"), function(x, i, j, exact) # not working
{
	nam <- names(x)
	x@.Data <- as(x, "list")[i]
	x@method.name <- x@method.name[i]
	stopifnot(validObject(x))
	if(!sameNames(nam, names(x)))
	{ stop("names bad"); browser()}#marta July2014 message
	x
})
#
setMethod("names", signature(x = "probs.to.combine"), function(x)
{
	names(x@.Data[[1]])
})
setMethod("names", signature(x = "combination.prob"), function(x){x@names})#names(x@.Data)
setMethod("names", signature(x = "plot_probs.to.combine"), function(x)
{
	names(x@ prob)
})

setReplaceMethod("names", signature(x="combination.prob",value='character'),function(x,value){
	if(length(value)==length(x@.Data)){x@names<-value}
	x})#added by marta April 2013



#conversions----------------
setAs(from = "combination.prob", to = "Numeric", function(from){
	nn<-as(from, "numeric")
	nNumeric(nn)})
setAs(from = "combination.prob", to = "numeric", function(from){#added by marta April 2013
	z<-as(from@.Data, "numeric");names(z)<-names(from);
	stopifnot(is_prob(z[is.finite(z)]));z
})



#setAs(from = "combination.prob", to = "numeric", function(from)
#{
#	stopifnot(is_prob(from))
#	vec <- as.numeric(from)
#	names(vec) <- names(from)
#	stopifnot(is_prob(vec))
#	vec
#})
#---------new------------
k.as.list.dcombi<-function(x){
	assert.is(x,c("combination.prob","list"))
	if(is(x,c("combination.prob"))){
		cx<-as(x,'numeric')
		wx<-as(x@ weight,'numeric')
		nms<-x@ names
	}
	else{			
		cx<-MakeNames(x$combined)
		wx<-x$weight
		nms<-names(cx)
	}

	if(length(names(cx))==0){names(cx)<-nms}
	if(length(wx)==length(cx)){
		names(wx)<-nms
		stopifnot(identical(names(cx),names(wx)))}
    
	z1<-list(combined=cx,weight=wx)#,names=nms
	class(z1)<-c("list","combination.prob.list")
	z1
}
setAs(from = "probs.to.combine", to = "list", function(from){#added by marta April 2013
    z1<-lapply(from@.Data,FUN=function(x){x<-x@.Data;names(x)<-names(from);x})
	names(z1)<-from@ method.name;
	stopifnot(is_prob(unlist(z1)))
	#z1<-c(z1,list(method.name=from@ method.name))
	class(z1)<-c("list","probs.to.combine.list");z1
})
setAs(from = "combination.prob",to = "list", function(from){k.as.list.dcombi(from)})
setAs(from = "list", to = "probs.to.combine", function(from)
{
	nx<-k.cleanlist.prob2combine(x=from)	
	new.x<-nx$nx
	method.name<-names(new.x)[1:length(new.x)]
		
	nprobs.to.combine(object=new.x, method.name = method.name)

})
setAs(from = "list", to = "combination.prob", function(from)
{
	stopifnot(all(c("combined","weight")%in%names(from)))
	assert.are(from[c("combined","weight")],"numeric")
	from$combined<-MakeNames(from$combined)
	from$combined<-from$combined[is.finite(from$combined)]
	from$weight<-from$weight[is.finite(from$weight)]
	if(length(from$weight)==length(from$combined)){from$weight<-MakeNames(from$weight)}
	else{stopifnot(length(from$weight)==1)}
	new_combination.prob (object=from$combined, weight=from$weight)
})

setAs(from = "list", to = "plot_probs.to.combine", function(from){
	z1<-as(from[setdiff(names(from),"comb.prob")],"probs.to.combine")
	z2<-from$comb.prob
	new_plot_probs.to.combine (object=z1, comb.prob=z2)})
setAs(from = "plot_probs.to.combine", to = "list", function(from){#added by marta April 2013
	z1<-as(as(from,"probs.to.combine"),"list")
	zaux<-as(from@comb.prob,"numeric");names(zaux)<-names(from@comb.prob)
	z1<-c(z1,list(comb.prob=zaux))
	class(z1)<-c("list","plot_probs.to.combine.list");z1
})
setAs(from = "plot_probs.to.combine", to = "probs.to.combine", function(from){z<-from@prob;validObject(z);z})
#------------------------------
new_plot_probs.to.combine <- function(object, comb.prob){
	
	assert.is(object,"probs.to.combine")
	assert.is(comb.prob,c("combination.prob","numeric"))
	cps<-cpsn<-comb.prob
	if(length(cps)>0){
		cpsn<-rep(NA,length(object[[1]]))
		names(cpsn)<-names(object)
		fin<-intersect(names(object),names(cps))
		cpsn[fin]<-cps[fin]
		stopifnot(identical(names(cpsn),names(object)))}
    new("plot_probs.to.combine", prob=object, comb.prob = nNumeric(cpsn))}


new_probs.to.combine <- function(object, method.name= character(0))
{
	stopifnot(length(object) == length(method.name))
	new("probs.to.combine", object, method.name = method.name)
}
new_combination.prob <- function(object, weight)
{
	vec <- try(nNumeric(object))
	if(is_err(vec)){ stop("bad object")}#marta July2014 { message("bad vec"); browser()}
	nam <- try(names(object))
	if(is_error(nam)){object<-MakeNames(object);nam<-names(object) }#message("bad name assignment"); browser()
	names(vec) <- nam
	cprob <- try(new("combination.prob", vec, weight = nNumeric(weight),names=nam))#added names by marta April 2013
	if(is_err(cprob)){ stop("bad combination.prob")}#marta July2014 { message("bad cprob"); browser()}
	names(cprob) <- names(object)
	stopifnot(sameNames(cprob, object))
	cprob
}


#

nprobs.to.combine <- function(object, method.name = character(0)){
	assert.is(object, c("list","probs.to.combine"))
	if(is(object, "probs.to.combine") && is_nothing(method.name))
		object
	else if(!is(object, "probs.to.combine"))
	{
		if(all(sapply(object, function(elem){is.null(names(elem))})))
			object <- lapply(object, function(elem)
			{
				stopifnot(is_prob(elem))
				elem<-MakeNames(elem,nmvar='Y')
				#names(elem) <- make.names(1:length(elem))#XXX|:atention!!
				elem
			})
		if(is_nothing(method.name))
		{
			method.name <- if(is.character(names(object)))
				names(object)
			else
				make.names(1:length(object))
		}
		get_name <- function(x)
		{
			na <- names(x)
			if(!is(na, "character"))
			{ message("bad na: "); print(na); browser()}
			na
		}
		nam <- get_name(object[[1]])
		if(length(object) >= 2)
		{
			for(i in 2:length(object))
			{
				nam <- try(intersect(nam, get_name(object[[i]])))
				if(is_err(nam))
				{ message("bad nam, very bad"); browser()}
			}
		}
		else
			warning(paste("length(object):", length(object)))
		ord <- order(object[[1]][nam])
		get_p <- function(x)
		{
			#if(is(x, "empiricalNull")){x <- lfdr(x)}#XXX|:atention!!
			proba <- x[nam][ord]
			if(length(proba) == 0 || !is(x, "numeric") || !is_prob(proba))
			{ stop("something wrong")}#marta July2014 { message("something wrong"); browser()}
			stopifnot(length(proba) == length(nam))
			nNumeric(proba)
		}
		object <- lapply(object, get_p)
		if(!are(object, "numeric")){ stop("bad list")}#marta July2014 { message("bad list"); browser()}
		nam <- names(object[[1]])
		stopifnot(all(sapply(object, function(elem){sameNames(nam, names(elem))})))
		if(!all(sort(object[[1]]) == object[[1]]))
		{ message("bad sort"); plot(object[[1]], sort(object[[1]])); stop()}#marta July2014 browser()
		assert.are(object, "numeric")
		new_probs.to.combine(object, method.name = method.name)
	}
	else
		stop("bad args: probs.to.combine")
}




ncombination.prob <- function(object, lower.plausible = 0, upper.plausible = 1, method.name = character(0), nsample = numeric(0), independent = !is_nothing(nsample), tolerance = 0.005, verbose = FALSE, arithmetic = TRUE, call.browser = FALSE,
			     loss.type = "information", ...){
	if(is(object, "numeric") && is_prob(object))
		object <- as.list(object)
	assert.is(object, "list")
	if(!is(object, "probs.to.combine"))
		object <- nprobs.to.combine(object, method.name = method.name)
	stopifnot(is(object, "probs.to.combine"))
	get_plausible <- function(bound)
	{
		stopifnot(is_prob(bound))
		if(length(bound) == 1)
			bound <- rep(bound, length(names(object)))
		if(is.null(names(bound)) && length(names(object)) == length(bound))
			names(bound) <- names(object)
		else if(is.character(names(bound)) && is.character(names(object))) # length(names(bound)) != length(names(object)))
		{
			nam <- intersect(names(bound), names(object))
			bound <- bound[nam] # bound[names(object)]
			object <- object[nam]
#			names(bound) <- names(object)
		}
		else
		{	stop("names(bound) is too weird")}##marta July2014 { message("names(bound) is too weird"); browser()}
		if(!sameNames(object, bound))
		{ stop("!sameNames(object, bound)")}#marta July2014 message ; browser()
		plausible <- try(nNumeric(bound))
		if(is_err(plausible))
		{ stop("bad plausible")}#marta July2014 message ; browser()
		plausible
	}
	
	lower.plausible <- get_plausible(lower.plausible)
	upper.plausible <- get_plausible(upper.plausible)
	stopifnot(all(lower.plausible <= upper.plausible))
	combine <- function(x, y, w, arith)
	{
		stopifnot(length(x) == length(y) && length(w) == 1)
		combined <- function(x, y, w, arithmetic)
		{
			if(arithmetic)
				w * x + (1 - w) * y
			else
			{
				if(!is_prob(c(x, y))){ stop("bad c(x, y) 1")}#marta July2014 message ; browser()
				get_unnorm <- function(x, y, w)
				{
					if(!is_prob(c(x, y))){ stop("bad c(x, y) 2")}#marta July2014 message ; browser()
					exp(combined(x = log(x), y = log(y), w = w, arithmetic = TRUE))
				}
				unnorm1 <- get_unnorm(x = x, y = y, w = w)
				unnorm2 <- get_unnorm(x = 1 - x, y = 1 - y, w = w)
				unnorm1 / (unnorm1 + unnorm2)
			}
		}
		combined(x = x, y = y, w = w, arithmetic = if(missing(arith)) arithmetic else arith)
	}
	vec <- if(length(object) == 2 && independent && arithmetic && loss.type == "information")
	{
		if(verbose) message("independent case")
		plausible <- function(p){all(p >= lower.plausible & p <= upper.plausible)}
		p1 <- object[[1]]
		p2 <- object[[2]]
		if(!plausible(p1) || !plausible(p2)) { stop( "something is not plausible")}##marta July2014 message
		get_mixture.p <- function(w){combine(x = p1, y = p2, w = w)}
		mutual <- function(w)
		{
			mixture.p <- get_mixture.p(w = w)
			div <- function(P1)
			{
				if(verbose) message("computing divergence")
				di <- try(divergence.Bernoulli(P1 = P1, P2 = mixture.p, nsample = nsample, independent = independent, verbose = verbose))
				if(is_err(di)){ stop("bad divergence.Bernoulli")}##marta July2014 message
				stopifnot(length(di) == 1)
				if(verbose) message("getting Scalar")
				sca <- try(nScalar(di))
				if(is_err(sca)){ stop("bad scalar")}#marta July2014 message
				sca
			}
			comb <- try(combine(x = div(p1), y = div(p2), w = w))
			if(is_err(comb))
			{ stop("bad combine")}#marta July2014 message
			comb
		}
		if(verbose) message("optimizing")
		opt <- optimize(f = mutual, lower = 0, upper = 1, maximum = TRUE, ...)
		weight <- opt$maximum
		get_mixture.p(w = weight)
	}
	else if(!independent && is_nothing(nsample))
	{
		comb.p <- weight <- rep(NA, length((names(object))))
		for(i in 1:length(weight))
		{
			p <- sapply(object, function(elem){stopifnot(length(elem) == length(names(object))); elem[i]})
			stopifnot(is_prob(p) && length(p) == length(object))
			min_plausible <- lower.plausible[i]
			max_plausible <- upper.plausible[i]
			plausible.p <- p[p >= min_plausible & p <= max_plausible]
			if(length(plausible.p) == 0)
			{
				if(all(p <= min_plausible))
					min_plausible
				else if(all(p >= max_plausible))
					max_plausible
				else
				{ stop("unusual p")}#marta July2014 message
			}
			else if(length(plausible.p) == 1)
				plausible.p
			else
			{
				min_p <- min(plausible.p)
				max_p <- max(plausible.p)
				get_mixture.p <- function(w){combine(x = min_p, y = max_p, w = w)}
				maximized <- function(w)
				{
					loss.w <- function(P1, w) # for debugging with call.browser = TRUE
					{
						mixture.p <- get_mixture.p(w = w)
						quadratic.loss <- function(P1, P2){(P1 - P2) ^ 2}
						geometric.loss <- function(P1, P2){quadratic.loss(P1 = log(P1), P2 = log(P2))}
						if(loss.type == "information")
							divergence.Bernoulli(P1 = P1, P2 = mixture.p)
						else if(loss.type == "quadratic")
							quadratic.loss(P1 = P1, P2 = mixture.p) # (P1 - mixture.p) ^ 2
						else if(loss.type == "harmonic")
							quadratic.loss(P1 = 1 / P1, P2 = 1 / mixture.p)
						else if(loss.type == "geometric")
							geometric.loss(P1 = P1, P2 = mixture.p)
						else if(loss.type == "logarithmic")
							geometric.loss(P1 = P1, P2 = mixture.p) + geometric.loss(P1 = 1 - P1, P2 = 1 - mixture.p) # resembles information divergence
						else if(loss.type == "evidential")
							geometric.loss(P1 = odds(P1), P2 = odds(mixture.p))
						else
							stop("bad loss.type")
					}
					loss <- function(P1){loss.w(P1 = P1, w = w)}
					if(call.browser)
					{ message("call.browser == TRUE; see loss.w"); browser()}
					if(arithmetic && loss.type == "information")
						combine(x = loss(min_p), y = loss(max_p), w = w, arith = TRUE) # mutual information
					else
					{
						to.minimize <- max(loss(min_p), loss(max_p)) # minimax
						-to.minimize # minimax equals maximin of additive inverse
					}
				}
				opt <- optimize(f = maximized, lower = 0, upper = 1, maximum = TRUE, ...)
				w.plus <- weight[i] <- opt$maximum
				ok <- try(is(w.plus, "numeric"))
				if(is_err(ok) || !ok){ stop("bad w.plus")}#marta July2014 message
				stopifnot(abs(opt$objective - maximized(w = w.plus)) < tolerance && is_prob(w.plus))
				comb.p[i] <- get_mixture.p(w = w.plus)
			} # end inner else
		} # end for (was sapply)
		if(any(is.na(weight))){ stop("bad weight")}#marta July2014 message
		if(any(is.na(comb.p))){ stop("bad comb.p")}#marta July2014 message
		comb.p
	} # end outer else (!independent && is_nothing(nsample))
	else
		stop("not yet implemented for those arguments")
	names.compatible <- length(names(object)) %in% c(0, length(vec))
	if(!names.compatible)
	{ stop("incompatible names")}#marta July2014 message
	names(vec) <- names(object)
	stopifnot(is_prob(vec))
	if(verbose) message("creating new object")
	cprob <- new_combination.prob(vec, weight = weight)
	if(!sameNames(cprob, object))
	{ stop("!sameNames(cprob, object)")}#marta July2014 message
	cprob
} # end combination.prob
combine.prob <- function(x, y, tolerance = 0.005, independent = TRUE, ...)
{
	lis <- list(x, y)
	assert.are(lis, "numeric")
	stopifnot(are_prob(lis, tolerance = tolerance))
	stopifnot(sameNames(x, y))
	ncombination.prob(object = lis, tolerance = tolerance, independent = independent, ...)
}



nplot_probs.to.combine <- function(object, method.name = character(0), lty = 2:length(object), legend.x = "bottomright", type = "l", col = 2:length(object), lwd = if(type == "l") 4 else 1, comb.prob = numeric(0), comb.lty = 1, comb.col = 1, verbose = FALSE, ...)
{
	assert.is(object, c("list"))
	stopifnot(length(object) >= 2)
	if(length(col) == 1)
		col <- rep(col, length(object) - 1)
	stopifnot(length(object) == length(col) + 1)
	stopifnot(length(object) == length(lty) + 1)
	if(!is(object, "probs.to.combine"))
	{
		message("calling probs.to.combine on ", date())
		if(length(method.name)==0){method.name<-names(object)}
		object <- new_probs.to.combine(object, method.name = method.name)
	}
	method.name <- object@method.name
	stopifnot(length(object) == length(method.name))
	if(type == "l")
		object <- object[order(object[[1]])]
	message("calling plot on ", date())
	plot(object[[1]], object[[2]], xlab = method.name[1], ylab = "probability of null hypothesis truth", lty = lty[1], pch = lty[1], col = col[1], type = type, lwd = lwd, ...)
	graph <- if(type == "l")
		function(..., lty)
		{
			lines(..., lty = lty, lwd = lwd)
		}
	else
		function(..., lty)
		{
			points(..., pch = lty)
		}
	if(length(object) >= 3)
	{
		for(i in 3:length(object))
		{
			graph(object[[1]], object[[i]], lty = lty[i - 1], col = col[i - 1], ...)
		}
	}
	if(is_nothing(comb.prob))
		message("no combination")
	else
	{
		stopifnot(is_prob(comb.prob))
		nam.comb.prob <- names(comb.prob)
		comb.prob.0 <- comb.prob
		if(is.null(names(comb.prob)))
			warning("is.null(names(comb.prob))")
		else
			comb.prob <- comb.prob[names(object)]
		comb.prob.ok <- is_prob(comb.prob) && length(comb.prob) == length(names(object))
		if(!comb.prob.ok)
		{
			stop("!comb.prob.ok; c(is_prob(comb.prob), class(nam.comb.prob), class(names(comb.prob)), length(nam.comb.prob), length(names(comb.prob))): ", c(is_prob(comb.prob), class(nam.comb.prob), class(names(comb.prob)), length(nam.comb.prob), length(names(comb.prob))))
			
		}#marta July2014 message
		if(!sameNames(comb.prob, object)){
			warning("!sameNames(comb.prob, object)")
			comb.prob<-comb.prob[names(object)]}#marta July2014 message
		message(length(comb.prob), " combined probabilities:")
		print(stats(comb.prob))
		graph(x = object[[1]], y = comb.prob, lty = comb.lty, col = comb.col, ...)
	}
	leg <- if(type == "l")
		function(..., lty)
		{
			legend(..., lty = lty, lwd = lwd,cex=0.7)
		}
	else
		function(..., lty)
		{
			legend(..., pch = lty,cex=0.7)
		}
	abline(a = 0, b = 1, col = "gray")
	get_vec <- function(vec, add.if.comb)
	{
		if(is_nothing(comb.prob))
			vec
		else
		{
			if(verbose)
			{
				message("adding ", add.if.comb, " to legend")
				print(stats(comb.prob))
			}
			c(vec, add.if.comb)
		}
	}
	legend.legend <- get_vec(method.name[-1], "COMBINATION")
	if(verbose) print(list(legend.legend = legend.legend))
	leg(x = legend.x, lty = get_vec(lty, comb.lty), legend = legend.legend, col = get_vec(col, comb.lty))
	#invisible(new_plot_probs.to.combine(object = object, comb.prob = comb.prob))
} # end plot_probs.to.combine


#-----
setClass("est.combi", representation(combined ="numeric",weight = "numeric",info="list"))
setValidity("est.combi", function(object) {
  dx<-object@combined
  wx<-object@weight
  ix<-object@info
  nvar<-length(dx)
stopifnot(length(names(dx))>0)
stopifnot(length(wx)%in%c(1,nvar))
if(length(wx)==nvar){stopifnot(identical(names(wx),names(dx)))}
if(!are_prob(list(dx[is.finite(dx)],wx[is.finite(wx)]))){
	stop('error: combined or weight out of [0,1]')}

})
#
#----------
k.just.finite<-function(x){
	new.x<-lapply(x,FUN=function(z){z[is.finite(z)]})
	nx<-lapply(new.x,length)
	ind<-which(nx>0)
	new.x[ind]
}
k.indep.max2probs<-function(x,na.rm=F){
	m.new.x<-list2matrix(x)	
	m.new.vx<-apply(m.new.x,2,FUN=function(z){c(min(z,na.rm=na.rm),max(z,na.rm=na.rm))})
	new.x.asif<-matrix2list(m.new.vx);names(new.x.asif)<-c("mins","maxs")
	k.just.finite(new.x.asif)
}
k.makelist.prob2combine<-function(x){#allow NAs
	assert.is(x,c('list','numeric'))
	if(is(x,"numeric")){x<-as.list(x)}
	assert.are(x, "numeric")

	x<-lapply(x,FUN=function(z){MakeNames(z,nmvar='Y')})
	MakeNames(x,nmvar='X')
	}
k.cleanlist.prob2combine<-function(x){#remove NAs
	nx<-k.makelist.prob2combine(x=x)	
	new.x<-k.just.finite(x=nx)
	list(nx=new.x,ox=nx)
}
	
k.checklist.prob2combine<-function(x){#allow NAs
	assert.is(x,c("list"))
	orig.cp<-x
	lok<-unlist(lapply(orig.cp,length))
	nmok<-unique(unlist(lapply(orig.cp,names)))
	nvars<-length(nmok)
	names.probs.to.combine<-names(orig.cp)

	if(nvars==0||max(lok)>nvars){stop("some elements in probs do not have names")}
	if(length(names.probs.to.combine)==0||any(nchar(names.probs.to.combine)==0)){
		stop("some probs do not have names")}
	list(nmok=nmok,nvars=nvars,names.probs.to.combine=names.probs.to.combine)
}

k.in.probs<-function(x,y=NULL,method.name=character(0),independent=F,...){max.nprobs<-2
	assert.is(x,c('list',"probs.to.combine","numeric"))
	assert.is(y,c("NULL","numeric"))
 
	if(is(x,'numeric')&&length(y)==0){x<-as.list(x);y<-NULL}
	else if(are(list(x,y),'numeric')){x<-list(x=x,y=y);y<-NULL}

	if(is(x,"list")){
		z<-k.cleanlist.prob2combine(x=x);y<-NULL
		new.x<-z$nx
		old.x<-z$ox
		if(length(method.name)==0){method.name<-names(new.x)}
		
		zz<-nprobs.to.combine(object=new.x, method.name = method.name)
	}
	else if(is(x,"probs.to.combine")){
		method.name<-x@method.name
		new.x<-old.x<-as(x,"list")
		zz<-x
	}
	else{stop("problem in list of probs to combine")}

	#if(!independent && max.nprobs==2 && length(new.x)>2){
	#	new.x<-k.indep.max2probs(x=new.x,na.rm=T)
	#	names(new.x)<-paste(names(new.x),paste0(method.name, collapse="."),sep="_") 
	#}
	
	if(length(method.name)!=length(new.x)){stop("problem with method.name")}
	stopifnot(is(zz,"probs.to.combine"))
	
#-----
	info<-c(list(independent=independent,method=""),list(...))

	if(length(info)>0){
		method1<-paste('combin',length(old.x),'probs:',paste(method.name,collapse="."),sep='')
		method2<-paste('arithmetic',info$arithmetic,'-',info$loss.type,'-indep',info$independent,sep='')
		method<-c(method1,method2)}
	else if(length(info)==0){
		method<-paste('combin',length(old.x),'probs:',paste(names(old.x),collapse="."),sep='')
	}
	info$method<-method
#-----	
	list(probs.to.combine=zz,
	     names.probs.to.combine=method.name,
	     names.elems.probs.to.combine=names(zz),#ok.listprobs=as(zz,"list"),
		orig.listprobs=old.x,info=info)
}

#----------------------
new_est.combi <- function(x=NULL,probs,independent=F,method=character(0),info=list()){
	assert.is(probs,c("list"))
	assert.is(x,c('NULL','combination.prob','try-error'))
	
	ok<-k.checklist.prob2combine(x=probs)
	nvars<-ok$nvars;nmok<-ok$nmok

	combined<-rep(as.numeric(NA),nvars)
	names(combined)<-nmok
	    
	if(independent) {weight<-as.numeric(NA)}
	else{
		weight<-rep(as.numeric(NA),nvars)
		names(weight)<-nmok}
	
	if (is_error(x)){
		message("Combine failed ",Sys.time())
		method<-paste('FAILED-',method,sep='')
		x<-list(combined=combined,weight=weight)}
    	    
	z<-k.as.list.dcombi(x=x)
	nms<-names(z$combined)

	stopifnot(length(intersect(nms,nmok))>0)
	stopifnot(all(nms%in%nmok))
	combined[nms]<-z$combined[nms]
	stopifnot(identical(names(combined),nmok))
	if(!independent){
		weight[nms]<-z$weight[nms]
		stopifnot(identical(names(weight),names(combined)))}
	else{weight<-z$weight}
	
	assert.are(list(combined,weight),c("numeric"))
	    
	if(!all(is.na(combined))){
	    stopifnot(are_prob(list(combined[is.finite(combined)],weight[is.finite(weight)])))}
  
    info<-uniquex(c(list(method=method),info))
    if(!print_info){info<-list(method=method)}
 

new("est.combi",combined =combined,weight = weight,info=info)


}





#plots----------------

setMethod("plot", signature(x = "probs.to.combine", y = "missing"), function(x, y, ...)
{
	nplot_probs.to.combine(object = x, ...)
})
setMethod("plot", signature(x = "probs.to.combine", y = "numeric"), function(x, y, ...)
{
	nplot_probs.to.combine(object = x, comb.prob = nNumeric(y), ...)
})
setMethod("plot", signature(x = "probs.to.combine", y = "combination.prob"), function(x, y, ...)
{
	nplot_probs.to.combine(object = x, comb.prob = as(y,"Numeric"), ...)
})
setMethod("plot", signature(x = "numeric", y = "probs.to.combine"), function(x, y, x.name = "lowest plausible probability", ...)
{
	stopifnot(is_prob(x))
	stopifnot(length(y) == length(y@method.name))
	p <- nNumeric(x)[names(y)]
	stopifnot(sameNames(p, y))
	object <- c(list(p), as(y, "list"))
	stopifnot(sameNames(object[[1]], object[[2]]))
	method.name <- c(x.name, y@method.name)
	stopifnot(length(object) == length(method.name))
	message("calling plot_probs.to.combine on ", date())
	nplot_probs.to.combine(object = object, method.name = method.name, ...)
})
setMethod("plot", signature(x = "plot_probs.to.combine", y = "missing"), function(x, y, ...)
{
	nplot_probs.to.combine(object = as(x,"probs.to.combine"), comb.prob = x@comb.prob, ...)#x
})

ncomb.prob <- function(x,y=NULL,independent = !is_nothing(nsample),lower.plausible = 0, upper.plausible = 1,method.name = character(0), 
     nsample = numeric(0),tolerance = 0.005,arithmetic = TRUE,loss.type = c("information","quadratic","harmonic","geometric","logarithmic","evidential"),
     plots=FALSE,...){loss.type<-loss.type[1]#independent = !is_nothing(nsample)

	z<-k.in.probs(x=x,y=y,method.name = method.name,independent=independent,arithmetic =arithmetic,loss.type = loss.type )

	zo<-try(ncombination.prob (object=z$probs.to.combine, lower.plausible = lower.plausible, upper.plausible =  upper.plausible,
		method.name = method.name, nsample = nsample, independent = independent, tolerance = tolerance , verbose = F,
		arithmetic =arithmetic, call.browser = FALSE, loss.type = loss.type, ...))
    
	infox<-c(list(lower.plausible = lower.plausible, upper.plausible =  upper.plausible,
		method.name = method.name, nsample = nsample, tolerance = tolerance),z$info,
		z[-which(names(z)%in%c("probs.to.combine","orig.listprobs","info"))])
	zz<-new_est.combi(x=zo,probs=z$orig.listprobs,independent=independent,method=z$info$method,info=infox)
	
	if(plots){
		pb<-z$probs.to.combine
		plot(pb,zo)
	}
	est2list(zz)
	

}








