w.s<-function(ui, df, ci=rep(1, length(ui)), uc=sqrt(sum((ci*ui)^2)) ) {
	(uc^4)/sum( ((ci*ui)^4) / df) 
}

welch.satterthwaite<-function( ui, df, ci=rep(1, length(ui)), uc=sqrt(sum((ci*ui)^2)) ) {
	w.s(ui,df,ci,uc)
}

