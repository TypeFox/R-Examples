`obj.haplomax.HS.add` <-
function(genea,perf,CD,assoc,res.structure,pi.hap,cor.pere,cor.mere){

nbre.hap	=	length(pi.hap[[1]])
pi.Qh		=	pi.hap[[1]][assoc]
esp.Qh		=	rep(0,nbre.hap)
esp.Qh[assoc]	=	pi.Qh
proba.DL.pere	=	proba.DL(pi.Qh,esp.Qh,res.structure,pi.hap,cor.pere)
proba.DL.mere	=	proba.DL(pi.Qh,esp.Qh,res.structure,pi.hap,cor.mere)
proba.DL	=	proba.DL.pere+proba.DL.mere
fact.pere	=	factor(genea[,2])

#régression linéaire
res	=	aov(perf~0+fact.pere+proba.DL, weight=CD^2)

res

}

