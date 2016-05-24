expr.dev.test <-
function(xx, formula.full, formula.red=NULL, D.red=NULL, model.dat, test.vars, glm.family, perm=100, 
						method=c("chisqstat", "permutation"), cf="fisher", adjust=FALSE, 
						snowfall.args=list(parallel=FALSE), 
						snowfall.seed, save.stat.data=FALSE, file.save="GlobalDeviance_statistic_Tperm_and_Tobs.Rdata",
						use.save.stat.data=FALSE, file.use="GlobalDeviance_statistic_Tperm_and_Tobs.Rdata")
{
### Vorbereitung der Daten und Fehlermeldungen
	# einzelne Variablen des vollen Modells
    design.terms<-as.character(attr(terms(formula.full), "variables"))[-1]		

    # Test, ob alle Variablen des vollen Modells auch in den Daten vorhanden sind
    if(!all(design.terms %in% colnames(model.dat))) {		
    	stop("Error: The terms in the model 'formula.full' do not match variables in data 'model.dat'.")  
    }

    # Test, ob Gendaten und Phenotypdaten zusammenpassen (gleiche Anzahl von Patienten/Stichproben)
    if(ncol(xx) != nrow(model.dat)) {
    	stop("Error: The number of samples in data matrix 'xx' differs from number of samples in data 'model.dat'.")
    }
    
    # Test, ob nur Daten zu einem Y vorhanden sind
    if(is.vector(xx)) {
    	# Umwandlung in einzelnen Zeilenvektor
        xx<-t(as.matrix(xx))
    }
	# Test, ob in den Daten (Y) Zeilennamen vorhanden sind
    if(is.null(rownames(xx))) {
    	# wenn die Gene/Zeilen keine Namen werden die Gene durchnummeriert
        rownames(xx)<-1:nrow(xx)
    }

    # Test, ob alle Variablen getestet werden sollen
    if(missing(test.vars) || is.null(test.vars)) {
    	# Liste aller Gennamen
        test.vars<-list(rownames(xx))
    }
    # Test, ob Gruppen von Variablen (bei Genen: pathways) getestet werden sollen
    # Test, ob die Struktur der zu testenden Variablen stimmt (Vektor oder Liste)
    if(!(missing(test.vars) || is.vector(test.vars))) {
    	stop("Error: 'test.vars' should be a vector or list of groups of variables.")
    }
    # Test, ob eine einzelne Gruppen von Variablen getestet werden soll
    if(!is.list(test.vars)) {
    	# Vektor in eine Liste umwandeln
        test.vars<-list(test.vars)
    }
    # Test, ob test.vars (eine) Zahlen enthaelt
    if(is.numeric(unlist(test.vars))) {
    	# Zahlen listenweise in Character umwandeln
        test.vars<-lapply(test.vars, as.character)
    }
    # Test, ob Variablennamen mit den Daten zusammenpassen
    if(!all(unlist(test.vars) %in% rownames(xx))) {
    	stop("Error: The variable names in 'test.vars' do not correspond to variable names in 'xx'.")
    }
    
### Datenauswahl und Anzahlen
	# Auswahl der benoetigten Gendaten 
    xx2<-xx[unique(unlist(test.vars)), , drop=FALSE]
    # Anzahl von Patienten/Stichproben
    N.Subjects<-ncol(xx2)
    # Anzahl von Gruppen von Variablen
    N.tests<-length(test.vars)
    
### Designmatrizen
    # Designmatrix des vollen Modells
    D.full<-model.matrix(formula.full, data=model.dat)
    
    # Test, ob die Designmatrix des reduzierten Modells gegeben ist
    if(is.null(D.red)) {
    	# Designmatrix des reduzierten Modells
        D.red<-model.matrix(formula.red, data=model.dat)
    }
    
    # Variablen des vollen und des reduzierten Modells
    terms.full<-colnames(D.full)
    terms.red<-colnames(D.red)

    # Test, ob das reduzierte Modell im vollen Modell enthalten ist
    if(!all(terms.red %in% terms.full)) { 
    	stop("Error: The full model and the reduced model are not nested.")
    }

    # Test, auf fehlende Werte (es muessen mehr Stichproben als Zeilen in der Designmatrix des vollen Modells vorhannden sein)
    if(nrow(D.full) < N.Subjects) {
        stop("Error: Missing values in the model variables.")
    }

    method<-match.arg(method)
### chi-Quadrat-Statistik und theoretischer p-Wert    
    if(method == "chisqstat") {
    	# Anzahl interessierender Kovariablen
  	 	N.covars<-ncol(D.full) - ncol(D.red)
    	# Kovariablen
		covars<-terms.full[!(terms.full %in% terms.red)]
		# Anzahl Gene
		N.vars<-nrow(xx2)
		vars<-rownames(xx2)

    	# Devianzen des reduzierten Modells
    	res.red<-varwiselogLik(xx=xx2, D=D.red, glm.family=glm.family)
		
		# Devianzen des vollen Modells
    	res.full<-matrix(ncol=1, nrow=N.vars, dimnames=list(vars, NULL))
    	res.full<-varwiselogLik(xx=xx2, D=cbind(D.red, D.full[, covars]), glm.family=glm.family)
    	
    	# Differenz der max. log-Likelihoods: 2*(l_voll - l_red) in R: fit.red$deviance - fit.full$deviance
    	dev<-(res.red - res.full)

		# Freiheitsgrade
		df<-rep(N.covars, N.vars)

		# p-Wert
		chisq.p.value<-(1 - pchisq(dev, df))
		
    	test.result<-cbind(dev, df, chisq.p.value)
		colnames(test.result)<-c("deviance", "df", "chi.square.p.value")
		rownames(test.result)<-vars
		
		data.perm<-NULL
		data.original<-res.full
		perm<-0
    }
    
### Permutationstest, p-value
    if(method == "permutation") { 
    	p.perm<-NULL
    	    
		# Permutationstest
		p.perm<-PermTest(xx=xx2, formula.full=formula.full, D.full=D.full, D.red=D.red, model.dat=model.dat, 
        	glm.family=glm.family, perm=perm, test.vars=test.vars, cf=cf, adjust=adjust, snowfall.args=snowfall.args, 
        	snowfall.seed=snowfall.seed, save.stat.data=save.stat.data, file.save=file.save, 
        	use.save.stat.data=use.save.stat.data, file.use=file.use)

		data.perm<-p.perm$data.perm
		data.original<-p.perm$data.original
		perm<-p.perm$perm
		
		if(adjust) {
        	test.result<-cbind(p.perm$statistic, p.perm$df, p.perm$p.value[[1]], p.perm$p.value[[2]], p.perm$p.value[[3]])
        	colnames(test.result)<-c("statistic", "df", "p-value", "adj.p-value.maxT", "adj.p-value.BY")
        } else {
        	test.result<-cbind(p.perm$statistic, p.perm$df, p.perm$p.value[[1]])
        	colnames(test.result)<-c("statistic", "df", "p-value")
        }
		rownames(test.result)<-names(test.vars)
    }
    
   main.result<-list("method"=method, "number.of.variables"=N.Subjects, "number.of.permutations"=perm, 
   "formula.full"=formula.full, "formula.red"=formula.red, "test"=test.result, "data.perm"=data.perm, 
   "data.original"=data.original, "test.vars"=test.vars)
   class(main.result)<-"dev.test"

   return(main.result)
}
