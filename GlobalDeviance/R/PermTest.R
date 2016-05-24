PermTest <-
function(xx, formula.full, D.full, D.red, model.dat, glm.family, perm, test.vars, cf, adjust, snowfall.args, snowfall.seed, save.stat.data, file.save, use.save.stat.data, file.use)
{
    # Anzahl von Patienten/Stichproben               	
	N.Subjects<-ncol(xx)
	# Anzahl von Responsevariablen
	N.vars<-nrow(xx)
	# Namen der Responsevariablen
	vars<-rownames(xx)
	# Anzahl von Gruppen von Variablen
    N.tests<-length(test.vars)
	# Anzahl von interessierenden Kovariablen
	N.covars<-ncol(D.full) - ncol(D.red)
	# Name der interessierenden Kovariablen
	terms.full<-colnames(D.full)
	terms.red<-colnames(D.red)
	covars<-terms.full[!(terms.full %in% terms.red)]

	# test.vars: Vektor mit Zeilennamen der zu testenden Variablen
	# test.vars0: Liste/Vektor der zu testenden Variablen
	test.vars0<-test.vars
	if(N.tests == 1) {	# eine gruppe von Variablen
		test.vars<-unlist(test.vars)
		test.vars<-which(rownames(xx) %in% test.vars) - 1
		num.test.vars<-length(test.vars)
	} else {
		# Zeilennummer der Variablen
		test.vars<-lapply(test.vars, function(x) which(rownames(xx) %in% x))
		# Laenge der Gruppen von Variablen
		num.test.vars<-sapply(test.vars, length)
		test.vars<-unlist(test.vars) - 1
	}
	
	# Wie viele Permutationen gibt es/sind moeglich
	nperms<-.nPerms(D.full=D.full, model.dat=model.dat, formula.full=formula.full)
	faccounts<-nperms$counts		# Anzahl der jeweils gleichen Zeilen der Designmatrix
	nperm<-nperms$nPerms			# Anzahl der Permutationen

	use.permMat<-0				
	permMat<-1
	nperm.used<-perm			# Anzahl vorgegebener/durchzufuerender Permutationen
	
	# die vorgegebene Anzahl von Permutationen (perm) ist groeßer als die moegliche Anzahl von Permutationen (nperm)
	if(nperm < perm) {
		use.permMat<-permMat
		nperm.used<-nperm		# alle Permutationen werden verwendet
		
		print(paste("Note: Enumerating all", nperm, "permutations."))
		
		if(is.null(faccounts))  {						# design with continuous covariates
			permMat<-.allperms(nums=1:N.Subjects)
		} else {
			permMat<-.allpermsG(counts=faccounts, grouping=faccounts)	# factorial design
		}
		permMat<-apply(permMat, 2, order)  				# get possible orderings
	}

	if(!use.save.stat.data) {
		# log-Likelihood fuer das reduzierte Modell
		T.beob.red<-matrix(NA, ncol=1, nrow=N.vars, dimnames=list(vars, NULL))
		T.beob.red<-varwiselogLik(xx=xx, D=D.red, glm.family=glm.family)

		# log-Likelihood fuer das volle Modell
		T.perm.full<-vector("list", (nperm.used + 1))
		T.perm.full<-ResamplelogLik(xx=xx, D.full=D.full, glm.family=glm.family, nperm.used=(nperm.used+1), 
		covars=covars, use.permMat=use.permMat, permMat=permMat, snowfall.args=snowfall.args, snowfall.seed=snowfall.seed)

		# a) Teststatistik: T_beob R^{q x 1}
		T.beob<-matrix(NA, ncol=1, nrow=N.vars, dimnames=list(vars, NULL))
		T.beob<-(T.beob.red[, 1] - T.perm.full[[nperm.used+1]])
	
		# b) Teststatistik: T_perm R^{q x B}
		T.perm<-matrix(NA, nrow=N.vars, ncol=nperm.used, dimnames=list(vars, 1:nperm.used))
		for(j in 1:nperm.used) {
			T.perm[, j]<-(T.beob.red[, 1] - T.perm.full[[j]])
		}
	} else {
		load(file=file.use)
		
		if(!exists("T.beob")) {
			stop("T.beob is not find.")
		}
		if(!exists("T.perm")) {
			stop("T.perm is not find.")
		}
		
		if(dim(T.beob) != c(N.vars, N.covars)) {
			stop("T.beob has wrong dimensions.")
		}
		if(dim(T.perm) != c(N.vars, N.covars, nperm.used)) {
			stop("T.perm has wrong dimensions.")
		}
	}
	
	if(save.stat.data) {
		save(list=c("T.beob", "T.perm"), file=file.save)
	}
	
	p.T.beob<-matrix(NA, ncol=1, nrow=N.vars, dimnames=list(vars, NULL))
	p.T.perm<-matrix(NA, ncol=nperm.used, nrow=N.vars, dimnames=list(vars, 1:nperm.used))
	S.beob<-p.global<-stat<-matrix(NA, ncol=1, nrow=N.tests, dimnames=list(names(test.vars0), NULL))
	S.perm<-matrix(NA, ncol=nperm.used, nrow=N.tests, dimnames=list(names(test.vars0), 1:nperm.used))

	# c) p-Werte schaetzen
	p.T.beob[, 1]<-(0.5 +  sapply(1:N.vars, function(k) { 
			rowSums(rbind(T.perm[k, ] >= T.beob[k, ]))
	}) )/(nperm.used + 1)
	
	for(k in 1:N.vars) {
		p.T.perm[k, ]<-(0.5 + colSums( sapply(1:nperm.used, function(i) { 
			cbind(T.perm[k, i] >= T.perm[k, ])
		}) ) )/(nperm.used + 1)
	}

	# d) combining function
	S.beob[, 1]<-sapply(1:N.tests, function(k){
		eval(call(cf, p.T.beob[test.vars0[[k]], ] ))
	})
	
	S.perm<-sapply(1:N.tests, function(k){ sapply(1:nperm.used, function(i){
		eval(call(cf, p.T.perm[test.vars0[[k]], i] ))
	}) })

	# e) globaler p-Wert
	p.global[, 1]<-sapply(1:N.tests, function(i) { 
		sum( S.perm[, i] >= S.beob[i, ] )/nperm.used
	})

	if(adjust) {
		p.global.maxt<-NULL
		p.global.by<-NULL
	} else {
		p.global.maxt<-NULL
		p.global.by<-NULL
	}
	
	df.perm<-rep(N.covars, N.tests)*sapply(test.vars0, length)
	stat[, 1]<-S.beob
	
	return(list("statistic"=stat, "df"=df.perm, "p.value"=list(p.global, p.global.maxt, p.global.by), 
	"data.perm"=list("T.perm"=T.perm, "S.perm"=S.perm), "data.original"=list("T.beob"=T.beob, "S.beob"=S.beob), 
	"perm"=nperm.used))
}
