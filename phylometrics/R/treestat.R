#' Significance test of user defined phylogenetic metric
#'
#' This function conducts the significance test of a phylogenetic metric defined by users.
#' @param phy an object of class 'phylo'.
#' @param stlist a vector of tip labels with trait state '1'.
#' @param state a vector of '0' and '1' for trait state of each tip, in the same order as the tip labels.
#' @param func name of the function that calculates the metric. The first two inputs of the function need to be the state vector and the phylo class.
#' @param par values of additional parameters used in the above function.
#' @param traitevol the null model of trait evolution used in the significance test of the metric. traitevol = "TBM" is the threshold brownian motion with the same trait prevalence as observed. traitevol = "random" is the random reshuffle of tip states. If users only want to calculate the value of the metric or if the null model is defined in the 'func' already, traitevol = NULL.
#' @param a number of traits simulated under the null model of trait evolution.
#' @param alternative the alternative hypothesis for the significant test of the metric. alternative = "greater", if users want to test if the observed trait has significantly larger metric value than expected by the null model of trait evolution. alternative = "less", if users want to test if the observed trait has significantly smaller metric value. alternative = "two.sided" (default), if users have no prior knowledge on how the observed trait may differ from the null model.
#' @param simplify if TRUE, the output is simplified. The full output includes the p.value, metric value of the observed trait, metric values of the simulated traits, a matrix, each row of which is the state vector of a simulated trait. The simplified output includes the p.value.
#' @examples
#' phy <- treesim(pars=c(0.1,0.1,0.05,0.05,0.1,0.1),N0=50,N1=50,sampling.f=c(1,1),max.t=Inf)
#' treestat(phy, func=tars, traitevol="random", a=1000, alternative="two.sided", simplify=TRUE)
#' @details
#' In the example, we conduct significane test on a metric defined by 'tars' function for a simulated phylogeny. The null hypothesis is that the trait is randomly distributed across tips. We don't know if the trait will have larger or smaller metric value than expected, so we apply a "two.sided" test. 1000 traits are simulated to generate the null distribution of the metric.
#' The state vector is included in the phylo class as the output of treesim, so we don't need to input state here.
#'

treestat<-function(phy, stlist=NULL, state=NULL, func, par=NULL, traitevol=NULL, a=NULL, alternative="two.sided", simplify=T) {
	if (!is.null(stlist)) {state <- is.element(phy$tip.label,stlist)}
    if (is.null(stlist)&&is.null(state)) {state <- phy$tip.state}
    m<-sum(state)
    nspecies<-length(state)
    if (!is.null(traitevol)) {
    if (traitevol=="TBM") {
		aa<-order(phy$edge[,1],decreasing=T)
		phy$edge<-phy$edge[aa,]
		phy$edge.length<-phy$edge.length[aa]
		anc<-phy$edge[seq(from=1,by=2,length.out=length(phy$edge[,1])/2),1]
		des<-matrix(phy$edge[,2],ncol=2,byrow=T)
		DES<-matrix(0,phy$Nnode,nspecies)
		for (i in 1:phy$Nnode) {
			tmp<-des[i,]
			offs<-which(anc %in% des[i,])
			while (length(offs)>0) {
				tmp<-c(tmp,des[offs,])
				offs<-which(anc %in% des[offs,])
			}
			tmp<-tmp[tmp<=nspecies]
			DES[i,tmp]<-1
		}
		DES<-rbind(diag(nspecies),DES)
		BL<-phy$edge.length[c(order(phy$edge[,2])[1:nspecies],order(phy$edge[,2],decreasing=T)[-c((nspecies-1):(2*nspecies-2))])]
		BL<-c(BL,0)
		V<-crossprod(DES,DES*BL)
		bmstate<-mvtnorm::rmvnorm(n=a,sigma=V)
		bmthred<-apply(bmstate,1,quantile,m/nspecies)
		bmstate<-sweep(bmstate,1,bmthred,'<')
		simstate<-matrix(as.numeric(bmstate),a,nspecies)
	}
	if (traitevol=="random") {
		rdstate<-replicate(a,sample(state,size=nspecies,replace=F))
		simstate<-t(rdstate)
	}
	}
	if (is.null(traitevol)) {
		if (is.null(par)) {
			obsstat <- func(state,phy)
		} else {
			obsstat <- func(state,phy,par)
		}
	} else {
		if (!is.list(func)) {
			if (is.null(par)) {
				obsstat <- func(state,phy)
			} else {
				obsstat <- func(state,phy,par)
			}
			if (!is.null(traitevol)) {
				if (is.null(par)) {
					simstat <- apply(simstate,1,func,phy)
				} else {
					simstat <- apply(simstate,1,func,phy,par)
				}
				if (alternative=="less") {
					p.value <- sum(simstat<=obsstat)/a
				}
				if (alternative=="greater") {
					p.value <- sum(simstat>=obsstat)/a
				}
				if (alternative=="two.sided") {
					p.value <- 2*ifelse(obsstat<=median(simstat),sum(simstat<=obsstat)/a,sum(simstat>=obsstat)/a)
				}
			}
		} else {
			obsstat <- numeric(length(func))
			simstat <- matrix(NA,length(func),a)
			p.value <- numeric(length(func))
			for (i in 1:length(func)) {
				if (is.null(par[[i]])) {
					obsstat[i] <- func[[i]](state,phy)
				} else {
					obsstat[i] <- func[[i]](state,phy,par[[i]])
				}
				if (!is.null(traitevol)) {
					if (is.null(par[[i]])) {
						simstat[i,] <- apply(simstate,1,func[[i]],phy)
					} else {
						simstat[i,] <- apply(simstate,1,func[[i]],phy,par[[i]])
					}
					if (alternative[i]=="less") {
						p.value[i] <- sum(simstat[i,]<=obsstat[i])/a
					}
					if (alternative[i]=="greater") {
						p.value[i] <- sum(simstat[i,]>=obsstat[i])/a
					}
					if (alternative[i]=="two.sided") {
						p.value[i] <- 2*ifelse(obsstat[i]<=median(simstat[i,]),sum(simstat[i,]<=obsstat[i])/a,sum(simstat[i,]>=obsstat[i])/a)
					}
				}
			}
		}
	}
	if (!is.null(traitevol)) {
		if (simplify) {
			out <- p.value
		} else {
			out <- list(p.value=p.value,obsstat=obsstat,simstat=simstat,simstate=simstate)
		}
	} else {
		out <- obsstat
	}
	out
}
