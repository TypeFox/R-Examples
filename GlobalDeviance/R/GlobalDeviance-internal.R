.mchoose <-
function(counts) {
	out<-choose(sum(counts), counts[1])
	if(length(counts) > 2) {
		out<-out * .mchoose(counts[-1])
	}
	out
}
.nPermsG <-
function(counts, grouping) {
	total<-.mchoose(counts)
	if(any(!is.na(grouping))) {
		correction<-prod(factorial(sapply(unique(grouping[!is.na(grouping)]), function(cc) sum(grouping == cc, na.rm=TRUE))))
	} else {
		correction<-1
	}
	total/correction
}
.nPerms <-
function(D.full, model.dat, formula.full) {
	# Anzahl der Stichproben/Patienten
	n<-nrow(D.full)

	# alle Kovariablen aus 'model.dat', die im vollen Modell stehen
	design.terms<-as.character(attr(terms(formula.full), "variables"))[-1]
	design.terms<-intersect(design.terms, colnames(model.dat))

	# bei stetigen Variablen im Modell, sind alle n! Permutationen unterschiedlich
	# logical: continuous covariates
	continuous<-sapply(model.dat[, design.terms, drop=FALSE], is.numeric)
	# numeric: two-group variables may be 'numeric' and not 'factor'
	varlength<-sapply(model.dat[, design.terms, drop=FALSE], function(x) length(unique(x))) 
	# logical
	continuous<-ifelse(continuous & (varlength > 2), TRUE, FALSE)
	
	if(any(continuous)) {	# eine stetige Kovariablen
		out<-ifelse(n <= 100, factorial(n), Inf)
		counts<-NULL
	} else {				# andere Kovariablen
		# alle Zeilen der Designmatrix werden unterschiedlich
		unique.rows<-unique(D.full)
		# Nullvektor mit Laenge der Anzahl der gleichen Zeilen der Designmatrix
		Y<-numeric(nrow(unique.rows))	
		for(i in 1:nrow(unique.rows)) {	
			equal.rows<-t(D.full) == unique.rows[i, ]
			equal.rows<-apply(equal.rows, 2, all)
			Y[equal.rows]<-i
		}
		# Vektor mit Anzahl der jeweils gleichen Zeilen der Designmatrix
		counts<-sapply(unique(Y), function(x) sum(Y == x))
		# Anzahl der Permutationen
		out<-.nPermsG(counts, counts)
	}
	return(list(nPerms=out, counts=counts))
}
.allpermsG <-
function(counts, grouping) {
	n<-sum(counts)
	if(n == 1) {
		app<-which.max(counts)
	} else {
		total<-.nPermsG(counts, grouping)
		app<-matrix(, n, total)
		choosable<-(counts > 0) & (is.na(grouping) | (1:length(counts) %in% match(unique(grouping[!is.na(grouping)]), grouping)))
		choosable<-(1:length(counts))[choosable]
		ix<-0
		for(iy in choosable) {
			countstemp<-counts
			countstemp[iy]<-counts[iy] - 1
			groupingtemp<-grouping
			groupingtemp[iy]<-NA
			size<-.nPermsG(countstemp, groupingtemp)
			app[1,(ix+1):(ix+size)]<-iy
			app[2:n, (ix+1):(ix+size)]<-.allpermsG(countstemp, groupingtemp)
			ix<-ix + size
		}
	}
	app
}
.allperms <-
function(nums) {
	# nums: indices to be permuted
	n<-length(nums)
	if(n == 1) {
		app<-nums
	} else {
		app<-matrix(, n, factorial(n))
		for(ix in 1:length(nums)) {
			range<-1:factorial(n-1) + (ix - 1) * factorial(n-1)
			app[1,range]<-nums[ix]
			app[2:n, range]<-.allperms(nums[-ix])
		}
	}
	app
}
