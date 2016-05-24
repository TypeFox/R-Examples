ResamplelogLik <-
function(xx, D.full, glm.family, nperm.used, covars, use.permMat, permMat, snowfall.args, snowfall.seed) 
{
	# Anzahl von Patienten/Stichproben               	
	N.Subjects<-ncol(xx)
	# Anzahl von Responsevariablen
	N.vars<-nrow(xx)
	N.covars<-length(covars) 
	ord.perm<-NULL
	res.perm.tmp<-matrix(NA, nrow=N.vars, ncol=N.covars)
	res.perm.tmp2<-vector("list", nperm.used)
	terms<-colnames(D.full)
	
	# Spalten der Designmatrix auf die adjustiert werden soll
	D.var.adjust<-cbind(D.full[, !(terms %in% covars)])
	# Spalten der Designmatrix die permutiert werden sollen
	D.var.perm<-cbind(D.full[, (terms %in% covars)])
	colnames(D.var.perm)<-covars
	
	# paralleles Rechnen
	do.call("sfInit", snowfall.args)
	
	sfwrapper.dev1<-function(i) {
		# Permutation erzeugen
		if(i==nperm.used) {
			ord<-1:N.Subjects
		} else {
			ord<-sample(ord.perm, N.Subjects)
		}
		# <<- Objekt wird global (globale Variable), kann dann beschrieben werden, sonst (mit <-) wird das Objekt kopiert
		D.full.perm<-cbind(D.var.adjust, D.var.perm[ord, covars])
		res.perm.tmp<<-varwiselogLik(xx=xx, D=D.full.perm, glm.family=glm.family)
		return(res.perm.tmp)
	}
		
	sfwrapper.dev2<-function(i) {
    	# alle Permutationen
    	ord<-permMat[, i]
    	D.full.perm<-cbind(D.var.adjust, D.var.perm[ord, covars])
    	res.perm.tmp<<-varwiselogLik(xx=xx, D=D.full.perm, glm.family=glm.family)
   		return(res.perm.tmp)
	}
			
	if(snowfall.args$parallel) {
		# Eingabevariablen und Ausgabevariablen
		sfExport("D.var.adjust", "D.var.perm", "xx", "glm.family", "ord.perm", "N.Subjects", "N.covars", "covars", "res.perm.tmp")
		# Zufallszahlengenerator
		sfClusterSetupRNG(seed=snowfall.seed)
	}

	# zufaellig ausgewaehlte Permutationen werden verwendet
	if(use.permMat == 0) { 	
		ord.perm<-1:N.Subjects
    	res.perm.tmp2<-sfLapply(1:nperm.used, sfwrapper.dev1)
  	}
  	
  	# alle moeglichen Permutationen werden verwendet
  	if(use.permMat == 1) { 	
    	res.perm.tmp2<-sfLapply(1:nperm.used, sfwrapper.dev2) 	
	}
	
    sfStop()
	return(res.perm.tmp2)
}
