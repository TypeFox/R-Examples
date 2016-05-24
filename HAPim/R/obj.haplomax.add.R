`obj.haplomax.add` <-
function(perf,CD,assoc,res.structure,pi.hap,cor.pere,cor.mere){

nbre.hap	=	length(pi.hap[[1]])
pi.Qh           =	pi.hap[[1]][assoc]
esp.Qh          =	rep(0,nbre.hap)
esp.Qh[assoc]   =	pi.Qh
proba.DL.pere   =       proba.DL(pi.Qh,esp.Qh,res.structure,pi.hap,cor.pere)
proba.DL.mere   =       proba.DL(pi.Qh,esp.Qh,res.structure,pi.hap,cor.mere)
proba.DL        =       proba.DL.pere+proba.DL.mere

#régression linéaire
res             =       aov(perf~proba.DL, weight=CD^2)

res

}

