# ctab: oneway, twoway, multiway percentage tables
# first argument must consist of one or more factors
# or a table object (class table, xtabs, or ftable)
# digits: number of digits after the decimal (default 2)
# type: "n" for counts, "row", "column" or "total"
# for percentages (default "n")
# row.vars:
# col.vars: same usage as ftable, ignored for one- and
# two-way tables
# percentages: FALSE==> proportions are presented rather
# than percentages (default TRUE)

# comments to John Hendrickx <John_Hendrickx@yahoo.com>

ctab<-function(...,dec.places=NULL,digits=NULL,type=NULL,style=NULL,row.vars=NULL,col.vars=NULL,percentages=NULL,addmargins=NULL) {
	mk.pcnt.tbl<-function(tbl,type) {
		a<-length(row.vars)
		b<-length(col.vars)
		mrgn<-switch(type,
			column=c(row.vars[-a],col.vars),
			   row=c(row.vars,col.vars[-b]),
			 total=c(row.vars[-a],col.vars[-b]))
		tbl<-prop.table(tbl,mrgn)
		if (percentages) {tbl<-tbl*100}
		tbl
	}

	# options have default NULL so attributes of a ctab object can be used as default
	# defaults for other classes are assigned below
	
	# 01JUN2010
	# Starting R version 2.11.0, it is no longer possible to use attributes(...)
	# if dot-dot-dot contains more than one object, e.g. for ctab(focc,occ)
	# Place "..." in a list and create a table object if the list has a length
	# > 1. If "..." contains a single object, proceed a before
	args<-list(...)
	if (length(args) > 1) {
		if (!all(sapply(args,is.factor))) stop("If more than one argument is passed then all must be factors")
		tbl<-table(...)
	}
	else {
		if (is.factor(...)) {
			tbl<-table(...)
		}
		else if (is.table(...)) {
			tbl<-eval(...)
		}
		else if (class(...)=="ftable") {
			tbl<-eval(...)
			if (is.null(row.vars) && is.null(col.vars)) {
				row.vars<-names(attr(tbl,"row.vars"))
				col.vars<-names(attr(tbl,"col.vars"))
			}
			tbl<-as.table(tbl)
		}
		else if (class(...)=="ctab") {
			tbl<-eval(...)
			if (is.null(row.vars) && is.null(col.vars)) {
				row.vars<-tbl$row.vars
				col.vars<-tbl$col.vars
			}
			for (opt in c("dec.places","type","style","percentages","addmargins")) if (is.null(get(opt))) assign(opt,eval(parse(text=paste("tbl$",opt,sep=""))))
			tbl<-tbl$table
		}
		else {
			stop("first argument must be either factors or a table object")
		}
	}

	# defaults for options and checks for valid values
	if (!is.null(digits)) dec.places<-digits
	if (is.null(dec.places)) dec.places<-2
	stopifnot(as.integer(dec.places)==dec.places,dec.places>0)
	if (is.null(percentages)) percentages<-TRUE
	stopifnot(is.logical(percentages))
	if (is.null(addmargins)) addmargins<-FALSE
	stopifnot(is.logical(addmargins))

	types<-NULL
	choices<-c("n", "row", "column", "total")
	for(tp in type) types<-c(types,match.arg(tp,choices))
	type<-types

	# one dimensional table,restrict choices to "n" and "total"
	if (length(dim(tbl))==1) {
		if (is.null(type)) {
			type<-c("n","total")
			row.vars<-1
			if (is.null(style)) style<-"wide"
		}
		else type<-ifelse(type=="n","n","total")
	}
	else if (is.null(type)) type<-"n"
	style<-match.arg(style,c("long","wide"))
	if (is.null(style)) style<-"long"

	# use row.vars and col.vars to determine the
	# marginals to use when calculating percentages
	# start by translating names to variable positions
	nms<-names(dimnames(tbl))
	z<-length(nms)
	if (!is.null(row.vars) && !is.numeric(row.vars)) {
		row.vars<-order(match(nms,row.vars),na.last=NA)
	}
	if (!is.null(col.vars) && !is.numeric(col.vars)) {
		col.vars<-order(match(nms,col.vars),na.last=NA)
	}
	# calculate the other if only one is given
	if (!is.null(row.vars) && is.null(col.vars)) {
		col.vars<-(1:z)[-row.vars]
	}
	if (!is.null(col.vars) && is.null(row.vars)) {
		row.vars<-(1:z)[-col.vars]
	}
	# evidently, both row.vars and col.vars were NULL
	# assign the last variable to col.vars, the rest to row.vars
	if (is.null(row.vars) && is.null(col.vars)) {
		col.vars<-z
		row.vars<-(1:z)[-col.vars]
	}

	if (type[1] == "n") ctab <- tbl
	else ctab<-mk.pcnt.tbl(tbl,type[1])

	if (length(type) > 1) {
		# create the (percentage) tables, then convert them to data frames
		# stack the data frames, adding a new variable as percentage type
		tbldat<-as.data.frame.table(ctab)
		z<-length(names(tbldat))+1
		tbldat[z]<-1
		pcntlab<-type
		pcntlab[match("n",type)]<-"Count"
		pcntlab[match("row",type)]<-"Row %"
		pcntlab[match("column",type)]<-"Column %"
		pcntlab[match("total",type)]<-"Total %"
		for (i in 2:length(type)) {
			if (type[i] == "n") ctab <- tbl
			else ctab<-mk.pcnt.tbl(tbl,type[i])
			ctab<-as.data.frame.table(ctab)
			ctab[z]<-i
			tbldat<-rbind(tbldat,ctab)
		}
		tbldat[[z]]<-as.factor(tbldat[[z]])
		levels(tbldat[[z]])<-pcntlab
		ctab<-xtabs(Freq ~ .,data=tbldat)
		names(dimnames(ctab))[z-1]<-""
	}

	result<-NULL
	result$row.vars<-row.vars
	result$col.vars<-col.vars
	result$dec.places<-dec.places
	result$type<-type
	result$style<-style
	result$percentages<-percentages
	result$addmargins<-addmargins
	result$ctab<-ctab
	result$table<-tbl
	class(result)<-"ctab"
	result
}

print.ctab<-function(x,dec.places=x$dec.places,addmargins=x$addmargins,...) {
	if (length(dim(x$ctab))==1) {
		tbl<-x$ctab
		if (addmargins) tbl<-addmargins(tbl)
		if (x$style=="long") {
			tbl<-as.matrix(tbl)
			colnames(tbl)<-names(dimnames(x$ctab))
		}
	}
	else {
		row.vars<-x$row.vars
		col.vars<-x$col.vars
		a=length(row.vars)
		if (length(x$type)>1) {
			z<-length(names(dimnames(x$ctab)))
			if (x$style=="long") row.vars<-c(row.vars,z)
			else col.vars<-c(z,col.vars)
		}
		b=length(col.vars)
		tbl<-x$ctab
		mrgn<-c(row.vars[a],col.vars[b])
		# if the table contains counts and percentages of a factor
		if (length(dim(x$table))==1) mrgn<-1
		if (addmargins) tbl<-addmargins(tbl,margin=mrgn)
		tbl<-ftable(tbl,row.vars=row.vars,col.vars=col.vars)
	}

	if (!all(as.integer(tbl)==as.numeric(tbl))) tbl<-round(tbl,dec.places)
	print(tbl,...)
}

summary.ctab<-function(object,...) {
	summary(object$table,...)
}

