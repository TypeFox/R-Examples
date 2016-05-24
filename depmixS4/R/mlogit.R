# 
# MULTINOMIAL LINK FUNCTION USED IN FAMILY OBJECT
# 

mlogit <-
function(base=1) {
	# 	matrix formulation is possibly very inefficient?!?!?
	# 	moreover it does not admit of bases being different from 1??!?!?
	# 	is it ever used anywhere???????
	linkfun <- function(p,base) {
		lfun <- function(p,base) {
			p <- p/sum(p)
			beta <- numeric(length(p))
			if(any(p==1)) beta[which(p==1)]=Inf
			else beta[-base] <- log(p[-base]/p[base])
			return(beta)
		}
		if(is.matrix(p)) {
			beta <- t(apply(p,1,lfun,base=base))
		} else {
			beta <- lfun(p,base)
		}
		return(beta)
	}
	linkinv <- function(eta,base) {
		linv <- function(eta,base) {
			pp <- numeric(length(eta))
			if(any(is.infinite(eta)) || any(eta > log(.Machine$double.xmax)) || any(eta < log(.Machine$double.xmin))) {
				pp[which(is.infinite(eta))] <- 1
				pp[which(eta > log(.Machine$double.xmax))] <- 1 # change this to something better!
			} else {
				expb <- exp(eta)
				sumb <- sum(expb)
				pp[base] <- 1/sumb
				pp[-base] <- expb[-base]/sumb
			}
			return(pp)
		}
		if(is.matrix(eta)) {
			if(ncol(eta)==1) {
				pp <- as.matrix(apply(eta,1,linv,base=base)) # fixes problem with column matrix eta
			} else pp <- t(apply(eta,1,linv,base=base)) 	
		} else {
			pp <- linv(eta,base)
		}
		return(pp)
	}
	mu.eta <- function(eta) {
		if(length(eta)==1) return(eta-eta^2)
		if(is.vector(eta)) return(diag(eta)-outer(eta,eta))
	}
	valideta <- function(eta) {
		TRUE # fix me
	}
	
	name <- "mlogit"
	structure(list(linkfun=linkfun,
			linkinv=linkinv,
			mu.eta=mu.eta,
			valideta=valideta,
			name=name,
			base=base),
		class="link-glm")
}

