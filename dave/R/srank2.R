srank2<- function(veg,groups,method,y) {
	nspecies<- ncol(veg)
	nreleves<- nrow(veg)
	sp.names<- names(veg)
# Checking proper classification
	if(nreleves != length(groups)) {
		stop("Classification in 'groups' does not accord with data frame used!")
	}
# ff is either indval or F, pp is p, list with both in outputlist
	ff<- rep(0,nspecies)
	pp<- rep(0,nspecies)
	veg<- as.matrix(veg^y)
	if (method == "indval") {
		cat(method,"\n")
		o.ind<- indval(veg,groups)
		ff<- o.ind$indcls
		pp<- o.ind$pval
		ordered.species<- sp.names[order(-ff)]
		ordered.ff<- ff[order(-ff)]
		ordered.pp<- pp[order(-ff)]
#   output for Table 7.1 in Latex
#    for (i in 1:nspecies) {
#    cat(sprintf("%6.0f %1s %6.0f %1s %-30s %1s %10.3f %1s  %10.3g",i,"&",order(-ff)[i],"&",ordered.species[i],"&",ordered.ff[i],"&",ordered.pp[i]))
#    cat("\n")
#    }
#    end Latex
#    outputlist
		Fvalue<- NULL
		indval<- ff[order(-ff)]
	}
	
	if (method == "jancey") {
# aov, ff is F, pp is p, list with both in jancey.tab
		ff<- rep(0,nspecies)
		pp<- rep(0,nspecies)
		for (i in 1:nspecies) {
			onespecies<- as.matrix(veg[,i])
			model1<- aov(onespecies~groups)
			model2<- anova(model1)
			ff[i]<- model2[[4]][[1]]
			pp[i]<- model2[[5]][[1]]
#     jancey.tab[i,]<- c(names(veg)[i],ff[i],pp[i])
		}
# printed output
		ordered.species<- names(veg)[order(-ff)]
		ordered.ff<- ff[order(-ff)]
		ordered.pp<- pp[order(-ff)]
# output for Table 7.1 in Latex
#  for (i in 1:nspecies) {
#  cat(sprintf("%6.0f %1s %6.0f %1s %-30s %1s %10.3f %1s  %10.3g",i,"&",order(-ff)[i],"&",ordered.species[i],"&",ordered.ff[i],"&",ordered.pp[i]))
#  cat("\n")
#  }
# end Latex
#  }
		indval<- NULL
		Fvalue<- ff[order(-ff)]
	}
	o.srank<- list(method=method,rank=seq(1,nspecies,1),species.no=order(-ff),species=sp.names[order(-ff)],Indval=indval,F_value=Fvalue,error.probability=pp[order(-ff)])
}
