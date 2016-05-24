
# -------------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------------- #

#wijfmat <- matrix( c(
#92, 107, 247,
#255, 89, 89,
#92, 203, 92,
#255, 177, 17,
#170, 97, 187,
#255,255,95,
#255,137,235,
#145, 101, 62,
#193,193,193,
#92, 229, 214,
#201, 255, 135,
#255, 224, 204,
#173,45,92,
#227, 196, 239,
#226, 212, 149,
#204, 241, 255,
#87, 142, 82
#), ncol=3, byrow = TRUE)

#wf <- apply(wijfmat,1,function(z){
#	rgb(z[1]/255,z[2]/255,z[3]/255)
#})
wijffelaars <- c("#5C6BF7", "#FF5959", "#5CCB5C", "#FFB111", "#AA61BB", "#FFFF5F", "#FF89EB", "#91653E", "#C1C1C1", "#5CE5D6", "#C9FF87",
 "#FFE0CC", "#AD2D5C", "#E3C4EF", "#E2D495", "#CCF1FF", "#578E52")

draw = function(H,W,X,Y,alpha = 1,border = "black",bg = "white",vp=NULL, lwd = 1){
	grid.rect(x =  unit(X, "npc"), y = unit(Y, "npc"),
          width = unit(W, "npc"), height = unit(H, "npc"),
          just = c("left","bottom"),
          default.units = "npc", name = NULL,
          gp=gpar(col=border,fill = bg, alpha = alpha, lwd = lwd), draw = TRUE, vp = vp
	)
}


rmb = function(x,...){
	UseMethod("rmb")
}
rmb.table = function(x, col.vars = NULL, spine = FALSE, circular = FALSE, eqwidth = FALSE, cat.ord = NULL, 
                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL, residuals = NULL,
 model.opt = list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,label.opt = list(),vp = NULL,...){
	
	#vn <- names(dimnames(x))
	#if(is.null(vn)){
	#	vn <- LETTERS[1:length(dim(x))]
	#}
	#if(is.numeric(col.vars)){
	#		vn <- c(vn[-col.vars],vn[col.vars])
	#		col.vars <- sort(seq_along(	vn ) %in% col.vars)
	#}
	
	#formula <- as.formula(paste("~",paste(vn,collapse="+"),sep=""))
	
	#data <- as.data.frame(ftable(x))
	nd <- length(dim(x))
	data <- subtable(x,1:nd)
	formula <- as.formula(paste("~",paste(names(data)[1:nd],collapse="+"),sep=""))
	rmb.formula(formula, data, col.vars = col.vars, spine = spine, circular = circular,
			 eqwidth = eqwidth, cat.ord = cat.ord, freq.trans = freq.trans, max.scale = max.scale,
			 use.na = use.na, expected = expected, residuals = residuals, model.opt = model.opt, 
			 gap.prop = gap.prop, gap.mult = gap.mult, col = col,col.opt = col.opt, label = label,
			 label.opt = label.opt, vp = vp)
}

rmb.ftable = function(x, col.vars = NULL, spine = FALSE, circular = FALSE, eqwidth = FALSE, cat.ord = NULL, 
                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL, residuals = NULL,
 model.opt = list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,label.opt = list(),vp = NULL,...){
	
	rv <- names(attr(x,"row.vars"))
	cv <- names(attr(x,"col.vars"))
	col.vars <- NULL
	#if(is.null(col.vars)){
		col.vars <- c(rep(F,length(rv)),rep(T,length(cv)))	
	#}
	formula <- as.formula(paste("~",do.call(paste,c(as.list(c(rv,cv)),sep="+")),sep=""))
	data <- as.data.frame(x)
	rmb.formula(formula, data, col.vars = col.vars, spine = spine, circular = circular,
			 eqwidth = eqwidth, cat.ord = cat.ord, freq.trans = freq.trans, max.scale = max.scale,
			 use.na = use.na, expected = expected, residuals = residuals, model.opt = model.opt, 
			 gap.prop = gap.prop, gap.mult = gap.mult, col = col,col.opt = col.opt, label = label,
			 label.opt = label.opt,vp = vp)
}


rmb.glm = function(x, col.vars = NULL, spine = FALSE, circular = FALSE, eqwidth = FALSE, cat.ord = NULL, 
                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL, residuals = NULL,
 model.opt = list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,label.opt = list(),vp = NULL,...){
	
	model.family <- x$family$family
	stopifnot(model.family %in% c("poisson","binomial"))
	
	if(model.family == "poisson"){
	
		fi <- which(names(x$data) == "Freq")
		
		tms <- terms(x$formula)
		vars <- attr(tms,"term.labels")[ which(attr(tms,"order") == 1)]
		formula <- paste( toString(x$formula[2]),"~", paste(vars,collapse="+"),sep="" ) 
		
		if("resid.type" %in% names(model.opt)){
			rtp <- model.opt$resid.type	
		}else{
			rtp <- "pearson"
		}
		data <- subtable(x$data, match(vars,names(x$data)[-fi]), keep.zero=TRUE)
	
		x <- update(x, data = data)
		residuals <- residuals(x,type=rtp)
		expected <- x$fitted

	}
	if(model.family == "binomial"){
		
		tms <- terms(x$formula)
		vars <- attr(tms,"term.labels")[ which(attr(tms,"order") == 1)]
		vars <- c(vars[-1],vars[1])
		formula <- paste( "Freq ~ ",toString(x$formula[2]), paste(vars,collapse="+"),sep="" ) 
		
		nv <- length(vars)+1
		
		unweighted <- all(x$prior.weights==1)
		data <- subtable(x$data, 1:nv, keep.zero=TRUE, allfactor = TRUE)
		
		if(unweighted){
			x <- update(x, data = data, weights = data$Freq)
		}else{
			x <- update(x, data = data)
		}
		if("resid.type" %in% names(model.opt)){
			rtp <- model.opt$resid.type	
		}else{
			rtp <- "pearson"
		}
		
		residuals <- residuals(x,type=rtp)
		expected <- x$fitted
	}
	
	rmb.formula(formula, data, col.vars = col.vars, spine = spine, circular = circular,
			 eqwidth = eqwidth, cat.ord = cat.ord, freq.trans = freq.trans, max.scale = max.scale,
			 use.na = use.na, expected = expected, residuals = residuals, model.opt = model.opt, 
			 gap.prop = gap.prop, gap.mult = gap.mult, col = col,col.opt = col.opt, label = label,
			 label.opt = label.opt, vp = vp)
}



#rmb.polr = function(x, col.vars = NULL, spine = FALSE, circular = FALSE, eqwidth = FALSE, cat.ord = NULL, 
#                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL, residuals = NULL,
# model.opt <- list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,label.opt = list(),...){
#	
#	model.family <- x$family$family
#	stopifnot(model.family %in% c("poisson","binomial"))
#	
#
#		vars <- dimnames(attr(x$terms,"factors"))[[1]]
#		vars <- c(vars[-1],vars[1])
#		nv <- length(vars)
#		formula <- paste( "Freq~", paste(vars,collapse="+"),sep="" ) 
#		
#		data <- x$model
#		
#		unweighted <- ncol(data) == nv
#		
#		if(unweighted){
#			data <- subtable(data, 1:nv, keep.zero=FALSE, allfactor = TRUE)
#			x <- update(x, data = data, weights = data$Freq)
#		}else{
#			names(data)[which(!(names(data) %in% vars))] <- "Freq"	
#		}
#		
#		residuals <- residuals(x,type=resid.type)
#		expected <- x$fitted
#	
#	rmb.formula(formula, data, col.vars = col.vars, spine = spine, circular = circular,
#			 eqwidth = eqwidth, cat.ord = cat.ord, freq.trans = freq.trans, max.scale = max.scale,
#			 use.na = use.na, expected = expected, residuals = residuals, model.opt = model.opt, 
#			 gap.prop = gap.prop, gap.mult = gap.mult, col = col,col.opt = col.opt, label = label,
#			 label.opt = label.opt)
#}




rmb.formula = function(formula, data, col.vars = NULL, spine = FALSE, circular = FALSE, eqwidth = FALSE, cat.ord = NULL, 
             cut = NULL, innerval = 1,   freq.trans = NULL, num.mode = FALSE, max.scale = 1,  use.na = FALSE, expected = NULL, residuals = NULL,
 model.opt = list(), gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), label = TRUE,label.opt = list(), vp = NULL,...){
	
# use.expected.values = FALSE, mod.type = "poisson", resid.type = "pearson",
# resid.display = "both",max.rat = 2.0, resid.max = NULL, cut.rs = 5
# boxes = TRUE, lab.tv = FALSE, varnames = TRUE, abbrev = FALSE, lab.cex = 1.2, yaxis = TRUE

if(!is.null(vp)){
	if(!inherits(vp,"viewport")){
		stopifnot(length(vp)>1)
		
		if(is.null(current.viewport()$layout)){
					print("No layout specified. Try:")
					print("mat.layout <- grid.layout(nrow, ncol); grid.newpage(); vp.mat <- viewport(layout = mat.layout);pushViewport(vp.mat)")
				}
		
		vp=viewport(layout.pos.row = vp[1], layout.pos.col = vp[2])
	}
	
	pushViewport(vp)
}else{
	grid.newpage()
}


min.alpha <- 0.1

if(!is.null(expected)){
	if(is.logical(expected)){
		if(expected == FALSE){
			expected <- NULL	
		}
	}
}

# ----- get model parameters from model.opt or use defaults ------ #
if( "use.expected.values" %in% names(model.opt) ){
	use.expected.values <- model.opt$use.expected.values
}else{
	use.expected.values <- FALSE
}
if( "mod.type" %in% names(model.opt) ){
	mod.type <- model.opt$mod.type
}else{
	mod.type <- "poisson"
}
if( "resid.type" %in% names(model.opt) ){
	resid.type <- model.opt$resid.type
}else{
	resid.type <- "pearson"
}


if( "resid.display" %in% names(model.opt) ){
	resid.display <- model.opt$resid.display
}else{
	resid.display <- "both"
}
if( "max.rat" %in% names(model.opt) ){
	max.rat <- model.opt$max.rat
}else{
	max.rat <- 2.0
}
if( "resid.max" %in% names(model.opt) ){
	resid.max <- model.opt$resid.max
}else{
	resid.max <- NULL
}
if( "cut.rs" %in% names(model.opt) ){
	cut.rs <- model.opt$cut.rs
}else{
	cut.rs <- 5
}


if( "bg" %in% names(col.opt) ){
	bg <- col.opt$bg
}else{
	bg <- NA #hsv(0,0,0, alpha=0.02)#hcl(240,10,85)
}

if( "col2" %in% names(col.opt) ){
col2 <- col.opt$col2
}else{
col2 <- alpha("grey",0.25)
}

if( "ang" %in% names(col.opt) ){
	ang <- col.opt$ang
}else{
	ang <- 360
}



if( "bgs" %in% names(col.opt) ){
	bgs <- col.opt$bgs
}else{
	bgs <- hsv(0,0,0, alpha=0.02)
}
if( "alpha.r" %in% names(col.opt) ){
	stopifnot( col.opt$alpha.r <= 1 & col.opt$alpha.r > 0 )
	alpha.r <- -log(col.opt$alpha.r)
}else{
	alpha.r <- -log(0.2)
}


if( "line.col" %in% names(col.opt) ){
	line.col <- col.opt$line.col
}else{
	if(!is.null(expected)){
		if(spine | circular){
			 line.col <- NA
		}else{
			line.col <- alpha(1,0.9)
		}
	}else{
		line.col <- ifelse( !spine & !circular, "darkgrey", NA)
	}
}

if( "rect" %in% names(col.opt) ){
rect <- col.opt$rect
}else{
rect <- bgs
}

# ----- get label parameters from label.opt or use defaults ------ #
if( "yaxis" %in% names(label.opt) ){
	yaxis <- label.opt$yaxis
}else{
	yaxis <- TRUE & !num.mode
}
if( "boxes" %in% names(label.opt) ){
	boxes <- label.opt$boxes
}else{
	boxes <- TRUE
}
if( "lab.tv" %in% names(label.opt) ){
	lab.tv <- label.opt$lab.tv
}else{
	lab.tv <- FALSE
}
if( "varnames" %in% names(label.opt) ){
	varnames <- label.opt$varnames
}else{
	varnames <- TRUE
}
if( "abbrev" %in% names(label.opt) ){
	abbrev <- label.opt$abbrev
}else{
	abbrev <- NULL
}
if( "lab.cex" %in% names(label.opt) ){
	lab.cex <- label.opt$lab.cex
}else{
	lab.cex <- 1.2
}


if( "s0" %in% names(label.opt) ){
	s0 <- label.opt$s0
	#s1 <- label.opt$s[2]
	#s2 <- label.opt$s[3]
	#s4 <- label.opt$s[4]
}else{
	s0 <- 0.02 # border
}
	s1 <- 0.04 # yaxis
	s2 <- 0.06 # labs
	s3 <- 0.1  # expected scale


tmp <- which(cut == FALSE)
if(any(tmp)) cut[tmp] <- 0

tmp <- which(cut == TRUE)
if(any(tmp)) cut[tmp] <- 12

# ----- parameter check 1 ----- #
		
		dset <- data.frame(data)
		f <- as.formula(formula)
		cut.rv <- TRUE
		if(max.scale > 1 ){
			stop(simpleError("Wrong max.scale specification! (0,1]"))
		}
				
		if(!( (gap.prop < 1) & (gap.prop > 0) ) ){
			stop(simpleError("Wrong gap.prop specification! (0,1)"))
		}
		if( !((min.alpha <= 1) & (min.alpha >= 0)) ){
			stop(simpleError("Wrong alpha specification! [0,1]"))
		}
		
	
		stopifnot(is.logical(label)|is.numeric(label))
		
		#if( !("base.alpha" %in% ls()) ){
		#	base.alpha <- 1	
		#}
	
	
		
# ----- terms extraction ----- #
	fterms <- attr(terms(f),"term.labels")
	nv <- length(fterms)
	
# ----- cutting numeric variables ----- #

int.vars <- sapply(fterms, function(z){ 

v <- data[,match(z,names(data))]
if(is.numeric(v)){
	all( (v%%1)==0 )
}else{
	FALSE	
}
})

if(length(cut) > 1){
	int.vars[cut>0] <- FALSE	
}

# integer variables will not be cut unless the value of cut at this particular position is positive.
# cut variable i iff cut[i] > 0 and Vi is numeric.

num.vars <- sapply(fterms, function(z) is.numeric(data[,match(z,names(data))]))




cuttf <- (cut > 0) & ( !int.vars   & num.vars )

if(is.null(cut)){
	cut <- rep(12,nv)	
}	
if(length(cut) == 1){
	cut <- rep(cut, nv)
}
if(length(cut) < nv){
	cut <- c(cut, rep(0,nv-length(cut)))
}



for(i in fterms[cuttf]){
	j <- which(names(data) == i)
	if(innerval < 1){
		iv <- innerval(dset[,j],p=innerval)
		too.low <- which(dset[,j] < iv[1])
		too.high <- which(dset[,j] > iv[2])
		#data[c(too.high,too.low),j] 	<- NA
		dset[c(too.high,too.low),j] 	<- NA
	}
	dset[,j] <- cut(dset[,j], breaks = cut[match(i, fterms)], right = FALSE)	
	levels(dset[,j]) <- paste(colsplit(levels(dset[,j]),",",1:2)[,1],")", sep="")
}
#if(innerval < 1 & any(num.vars)){
	#kill <- sapply(dset[,match(fterms,names(data))], is.na)
	#kill <- apply(kill,1,any)
	#if(any(kill)){
	#	dset <- dset[-which(kill),]
		# dims of data and dset are now potentially different!
	#}
#}
if(!use.na){
	kill <- sapply(dset[,match(fterms,names(data))], is.na)
	kill <- apply(kill,1,any)
		#dset <- na.omit(dset)
		if(any(kill)){
			dset <- dset[-which(kill),]
		}
	}else{
		for( i in 1:(ncol(dset)-1) ){
			ind <- is.na(dset[,i])
			if(any(ind)){
				levels(dset[,i]) <- c(levels(dset[,i]),"N/A")	
				dset[which(ind),i] <- "N/A"
			}	
		}	
	}
	
	
	
	if(is.numeric(col.vars)){
			col.vars <- seq_along(	fterms ) %in% col.vars
		}
#row.vars <- seq_along(fterms)[-c(col.vars)]

	
	
	if( !is.null(expected) ){
		both <- resid.display %in% sapply(1:4, function(z) abbreviate("both",z))
		res.val.only <- resid.display %in% c(sapply(1:6, function(z) abbreviate("values",z)),sapply(1:6, function(z) abbreviate("value",z)))
		disp.res <- both | res.val.only
		if( disp.res & all( sapply(expected, function(z) all(z == TRUE))) ){
			if( mod.type == "poisson" ){
				expected <- list(1:(nv-1),nv)	
			}else{
				expected <- as.list(1:(nv-1))	
			}
		}
		if( resid.display %in% sapply(1:6, function(z) abbreviate("values",z)) ){
			col <- "rgb"
			col.opt$alpha <- 0.3 	
		}
	}else{
		res.val.only <- FALSE	
	}
	

	
			
	nv0 <- nv
	if(is.null(col.vars)){
		col.vars <- rep(c(T,F),nv)[1:nv]	
	}
		
	col.vars[nv] <- TRUE
	tv <- fterms[nv]
	
	stopifnot(is.logical(col.vars) & length(col.vars) == nv)
	
	ind <-  match(fterms,names(dset))
	#if(num.mode){
	#	dset.num <- as.data.frame(data)[,ind]	
	#}

# 	>>> a dummy-variable is used to handle the case of no row-variables
	if( sum(!col.vars) < 1 ){
		dset <- data.frame(probability = rep("prob",nrow(dset)),dset)
		nv<-nv+1
		col.vars <- c(F,col.vars)
		ind <- c(1,ind+1)
	}
	if( sum(col.vars) < 2 ){
		dset <- data.frame(distribution = rep(tv,nrow(dset)),dset)
		nv<-nv+1
		col.vars <- c(T,col.vars)
		ind <- c(1,ind+1)
		if("probability" %in% names(dset)){ expected <- NULL }
	}

# 	>>> preparing the labels for plotting including the abbreviation option	
	if(is.numeric(label)){
			lab <- TRUE
			tmp <- rep(FALSE,nv)
			tmp[label] <- TRUE
			label <- tmp
	}else{
			lab <- any(label)	
			if(length(label) == 1 & lab){
				label <- rep(TRUE,nv)	
			}
			if(!lab){
				label <- rep(FALSE,nv)	
			}
	}
	if(length(abbrev) < length(ind)){
		abbrev <- rep(abbrev,length(ind))[seq_along(ind)]
	}
		
	
	
	if(is.null(abbrev)){
		rclabs <- lapply(dset[,ind],function(x) levels(as.factor(x)))
	}else{
		#rclabs <- lapply(dset[,ind],function(x) abbreviate(levels(as.factor(x)),abbrev))
		rclabs <- mapply(function(y,z) abbreviate(levels(y),z),y = dset[,ind], z = as.list(abbrev),SIMPLIFY = FALSE)
	}
	lab.tv <- ifelse(spine,FALSE,lab.tv)

	
	label[nv] <- lab.tv
	if( sum(label) == 0 ){
		lab <- FALSE	
	}
	rlabs <- lapply(which( (!col.vars)*label > 0),function(x) rclabs[[x]])
	clabs <- lapply(which(col.vars*label > 0),function(x) rclabs[[x]])
	
#	orig.labs <- lapply(dset[,ind],function(x) levels(as.factor(x)))
	if( length(terms(f)) > 2 ){
		
		fi <- which( suppressWarnings(names(dset) == terms(f)[[2]]))
		if("Freq" %in% names(dset)[-fi]){
			names(dset)[ which(names(dset) == "Freq") ] <- "Freq2"
		}
		names(dset)[fi] <- "Freq"
	}
	dset <- subtable(dset,ind,keep.zero=FALSE,allfactor=TRUE)
	

	
	
	# ----- predefined expected values ----- #

if( is.numeric(expected) ){
	stopifnot( length(expected) == length(dset$Freq) )
	mod.type <- "gen"
	pred <- expected
}else{
	if(suppressWarnings(max(unlist(expected))) > (nv - (mod.type=="logit")) ){
		stop(simpleError("Wrong expected specification!"))
	}	
}

	if(is.null(cat.ord)){
		ntc <- nlevels(dset[,nv])
		cat.ord <- 1:ntc
	}	
	if(length(cat.ord) > 1 & !spine){
		levels(dset[,nv])[which(! 1:nv %in% cat.ord )] <- NA
		dset[,nv] <- factor(dset[,nv],levels=levels(dset[,nv])[rank(cat.ord)])
		if(lab.tv){
			clabs[[length(clabs)]] <- levels(dset[,nv])	
		}
	}
	
	if(!use.na){
		dset <- na.omit(dset)
	}else{
		for( i in 1:(ncol(dset)-1) ){
			ind <- is.na(dset[,i])
			if(any(ind)){
				levels(dset[,i]) <- c(levels(dset[,i]),"N/A")	
				dset[which(ind),i] <- "N/A"
			}	
		}	
	}
	
	
# ----- descriptive parameters ----- #	
	ntc <- nlevels(dset[,nv])
	ntc0 <- ntc
	#nlvl <- sapply(dset[,(sapply(dset,class)=="factor")], nlevels)
	nlvl <- sapply(dset[,-ncol(dset)], nlevels)
	col.nlvl <- nlvl[which(col.vars)]
	row.nlvl <- nlvl[which(!col.vars)]								

	nc <- prod(col.nlvl)
	nr <- prod(row.nlvl)

	ncv <- length(col.nlvl)
	nrv <- length(row.nlvl)




# ----- computation of the underlying (relative) frequencies ----- #		
	tt1 <- ftable(tapply(dset$Freq,as.list(dset[,1:nv]),sum),col.vars=which(col.vars))
	tt2 <- spread(ftable(tapply(dset$Freq,as.list(dset[,1:(nv-1)]),sum),col.vars=which(col.vars[1:(nv-1)])),ncol=ntc)

	tt1[is.na(tt1)] <- 0 
	tt2[is.na(tt2)] <- 1 

	H0 <- tt1/tt2
	H <- H0

# ----- a few more auxiliary variables ----- #	

	rind <- which( (!col.vars)*label > 0)
	cind <- which(col.vars*label > 0)
	nrl <- length(rind)*lab
	ncl <- length(cind)*lab

num.mode.x <- num.mode
num.mode.y <- num.mode
col.gap.prop <- ifelse(nc > ntc, gap.prop * min(1,ncv/nrv), 0) 
row.gap.prop <- ifelse(nr > 1, gap.prop * min(1,nrv/ncv), 0) 	
	
	if( is.factor(data[, ind[cind[length(cind)]]  ])){
		num.mode.x <- FALSE	
	}
	if( is.factor(data[, ind[rind[length(rind)]]  ]) ){
		num.mode.y <- FALSE	
	}
	
	if(num.mode.x & num.mode.y){
		gap.prop <- 0.001
	}
if(num.mode.y){
	row.gap.prop <- 0.001
}
if(num.mode.x){
	col.gap.prop <- 0.001
}

# ----- model computation in expected mode ----- #	
	if(!is.null(expected) ){
		
		int.dset <- lapply(dset,as.integer)
		dset <- dset[do.call("order",int.dset),]
				
		if(!(mod.type %in% c("poisson","polr", "gen"))){
			stop(simpleError("Wrong mod.type specification!"))
		}
		
				
				
		# ------ logit response residuals ----- #
		if(mod.type == "polr"){ 
			if(! (resid.type %in% c("response","prob","prop","rat"))){
				print(simpleWarning("Argument resid.type ignored. Only response residuals are implemented for polr models."))
				resid.type<-"response"	
			}
			nameslist <- names(dset)[1:(nv-1)]		
			nameslist <- nameslist[which(  !(nameslist %in% c("probability","distribution"))  )]
			single.terms <- paste( nameslist ,collapse = "+")
			interaction.terms <- do.call("paste",c(lapply(expected,function(x) do.call("paste",c(as.list(nameslist[x]),sep="*"))),sep="+"))
			
			full.terms <- paste(nameslist,collapse = "*")
			mod.formula <- as.formula(paste(fterms[nv0],"~",single.terms,"+",interaction.terms,sep=""))
			
			modS.formula <- as.formula(paste(fterms[nv0],"~",full.terms,sep=""))
			if (requireNamespace("MASS", quietly = TRUE)) {	
				modS <- MASS::polr( formula = as.formula(modS.formula),data = dset, weights = dset$Freq, method ="logistic")
			mod <- MASS::polr( formula = as.formula(mod.formula),data = dset, weights = dset$Freq, method ="logistic")
			}
			pred <- predict(mod,newdata = subtable(dset,c(1:(nv-1)),keep.zero=F),type="probs")
			fit.values <- as.vector(t(pred))
			tt1r <- ftable(tapply(fit.values,as.list(dset[,1:nv]),sum),col.vars=which(col.vars))
			resid.mat <- (H - tt1r)
			
						
			p.value<-anova(mod,modS)[2,7]
			
			if( use.expected.values ){
				H0 <- tt1r#/tt2
				H00 <- H
				H <- H0
				resid.mat <- -resid.mat
			}else{
				H00 <- tt1r#/tt2
				# H <- H
			}
			
		}
		# ----- poisson model ----- #	
		if(mod.type == "poisson"){
			nameslist <- names(dset)[1:(nv)]		
			nameslist <- nameslist[which(  !(nameslist %in% c("probability","distribution"))  )]
			single.terms <- do.call("paste",c(as.list( nameslist ),sep="+"))
			interaction.terms <- do.call("paste",c(lapply(expected,function(x) do.call("paste",c(as.list(nameslist[x]),sep="*"))),sep="+"))
						
			mod.formula <- paste("Freq~",single.terms,"+",interaction.terms,sep="")
			
			mod <- glm( as.formula(mod.formula),data = dset,  family = "poisson")
			if(resid.type %in% c("prob","rat","prop")){
				mod.residuals <- residuals(mod,type="response")
				resid.mat <- ftable(xtabs(mod.residuals~.,data=dset[,1:nv]),col.vars=which(col.vars))
				resid.mat <- resid.mat / tt2
				resid.type <- "response"
			}else{
				mod.residuals <- residuals(mod,type=resid.type)
				resid.mat <- ftable(xtabs(mod.residuals~.,data=dset[,1:nv]),col.vars=which(col.vars))
			}
		    
			

			#resid.mat <-  ftable(tapply(mod.residuals,as.list(dset[,1:nv]),sum),col.vars=which(col.vars))
			
			Df  <- mod$df.residual
			Dev <- mod$deviance
			p.value <- 1-pchisq(Dev,Df)
			
			pred <- predict(mod,type="response")
			#pred.mat <-  ftable(tapply(pred,as.list(dset[,1:nv]),sum),col.vars=which(col.vars))
			pred.mat <-  ftable(xtabs(pred~.,data=dset[,1:nv]),col.vars=which(col.vars))

			if( use.expected.values ){
				H0 <- pred.mat/tt2
				H00 <- H
				H <- H0
				resid.mat <- -resid.mat
			}else{
				H00 <-  pred.mat/tt2
			}
		}
	
		
		if(mod.type == "gen"){
			
			p.value <- 1
			
			# pred is already defined by expected = c(...)
			pred.mat <- ftable( xtabs(pred~., data = data[,1:nv]), col.vars = which(col.vars))
			
		if(is.null(residuals)){	
				mod.residuals <- dset$Freq - pred
				resid.type <- "response"
			}else{
				stopifnot(length(residuals) == length(pred))
				mod.residuals <- residuals
			}

			resid.mat <- ftable( xtabs(mod.residuals~., data = data[,1:nv]), col.vars = which(col.vars))

			if( use.expected.values ){
				H0 <- pred.mat/tt2
				H00 <- H
				H <- H0
				resid.mat <- -resid.mat
			}else{
				H00 <-  pred.mat/tt2
			}
		}
		
		
	}else{
		disp.res <- FALSE
		H00 <- H0
	}

# ----- modify parameters for rmbplot mode ----- #	
	if( spine | circular ){
		nlvl[nv] <- 1
		nc <- nc/ntc
		ntc <- 1
	}
# preparing the transformation parameters
if( is.list(freq.trans) ){
	stopifnot(freq.trans[[1]] == "sqrt")
	stopifnot(freq.trans[[2]] > 0)
		tfreq <- "sqrt"
		tfn <- freq.trans[[2]]
	
}else{
	tfreq <- freq.trans
	check <- FALSE
	if(is.null(tfreq)){
		tfreq <-"sqrt"
		tfn <- 1
		check <- TRUE
	}else{
		if( tfreq == "sqrt" ){
			tfn <- 2
			check <- TRUE
			}
		if( tfreq == "log" ){
			tfn <- exp(1)
			check <- TRUE
		}
		if( tfreq %in% c("1","c","const","eq","equal") ){
			tfreq <- "const"
			tfn <- 1
			check <- TRUE
		}
	}
	if(!check){
		stop(simpleError("Wrong freq.trans specification!"))	
	}
}
# ----- computation of rectangle widths, heights and x/y coordinates in matrix form ----- #	

	W <- matrix(1,ncol=nc,nrow=nr)
	S <- space(nlvl,col.vars,gap.prop=gap.prop,gap.mult=gap.mult,last.col.zero = TRUE,last.row.zero = FALSE)
# better usage of space here for the num.mode?

	x <- ( seq(0,1-1/nc,1/nc)   )*(1-col.gap.prop)
	
	if(nc*(nlevels(dset[,nv])^(spine|circular)) > nlevels(dset[,nv])){
		x <- x + c(0,cumsum(S[[1]]))*col.gap.prop/sum(S[[1]])
	}
	if( suppressWarnings(any(S[[2]])) ){
		y <- rev((seq(0,1-1/nr,1/nr))*(1-row.gap.prop) + c(0,cumsum(S[[2]]))*row.gap.prop/sum(S[[2]])) #rev test
	}else{
		y <- 0	
	}
		
	X <- spread(t(matrix( x )),nrow=nr)
	Y0 <- spread(matrix( y ),ncol=nc)
	Y <- Y0
	X2 <- X[,seq(1,nc,ntc)]
	Y2 <- Y[,seq(1,nc,ntc)]

	H2 <- matrix(1,ncol=nc/ntc,nrow=nr)
	tt3 <- ftable(tapply(dset$Freq,as.list(dset[,1:(nv-1)]),sum),col.vars=which(col.vars[1:(nv-1)]))
	
	tt3[is.na(tt3)] <- 0
	W2 <- tt3 / max(tt3)

	if(nc > nlevels(dset[,nv])){
		width.cor <- (ntc-1)*S[[1]][1]/sum(S[[1]])*col.gap.prop
	}else{
		width.cor <- 0
	}
	
	if( eqwidth ){
		W2trans <- W2
		if( tfreq == "const" ){
			W2trans <- 1*(W2>0)
		}
		if( tfreq == "log" ){
			W2trans <- log(tt3+1)/max(log(tt3+1))
		}
		if( tfreq == "sqrt" ){
			W2trans <- tt3^(1/tfn) / max(  tt3^(1/tfn) )
		}
		W3 <- W
		X3 <- X
	}else{
		if( tfreq == "const" ){
			W2trans <- 1*(tt3>0)
		}
		if( tfreq == "log" ){
			W2trans <- log(tt3+1)/max(log(tt3+1))
		}
		if( tfreq == "sqrt" ){
			W2trans <- tt3^(1/tfn) / max(  tt3^(1/tfn) )
		}
		W3 <- spread(W2trans,ncol=ntc)
		X3 <- spread(X[,seq(1,nc,ntc),drop=FALSE],ncol=ntc) + t( t(W3) * c(0:(ntc-1)) )/ncol(X)*(1-col.gap.prop) 
	}

	
# ----- ---------------- ----- #	
# ----- plotting section ----- #	
# ----- ---------------- ----- #

# ----- opening the basic viewport ----- #	
	
	#s0 <- 0.02 # border
	#s1 <- 0.04 # yaxis
	#s2 <- 0.06 # labs
	#s3 <- 0.1  # expected scale
	#dev.new()
	
	
	#grid.rect(gp=gpar(fill="white",col="white"))
	
	#if(exists("rmb.id",.GlobalEnv)){
		#popViewport()
		#grid.rect(gp=gpar(fill="white"))
	if(circular) yaxis <- FALSE
	#}else{
	#	.GlobalEnv$rmb.id <- rpois(1,lambda=10000)	
	#}
	vp0 <- viewport(x = 0.5, y = 0.5, width = 1- 2*s0 , height = 1 - 2*s0, just = "centre", name = "vp0")
	pushViewport(vp0)
	vp1 <- viewport(x = 1-yaxis*s1 - s3*(!is.null(expected) ) , y = 0, width = 1- nrl*s2 - yaxis*2*s1 - s3*(!is.null(expected)), height = 1-ncl*s2, just = c("right", "bottom"), name = "vp1")
	pushViewport(vp1)
	grid.rect(gp=gpar(fill=bg,col=NA))
	# alpha values on a log-scale
	#alpha <- spread(W2,ncol=ntc)*base.alpha
	
	#alpha <- spread(exp(2*(W2-1)),ncol=ntc)*base.alpha
	
	if(!("alpha" %in% labels(col.opt))){
		alpha <- exp(alpha.r*(W2-1))
	}else{
		alpha <- matrix(col.opt$alpha, ncol=ncol(W2), nrow=nrow(W2))				
	}
					
	
	#alpha[alpha < min.alpha] <- min.alpha
	if(!is.null(expected)){
		alpha[alpha < 1] <- 1	
	}

# ----- preparing the color matrix ----- #	

			if(is.null(expected) | res.val.only ){
				if(any(c("hsv","hcl","rgb","seq","sequential","sqn","sqt","div","diverging","diverge",
				"d","s","h","t","heat","ht","ter","terrain","Wijffelaars", "w", "q17") %in% col)){
				
					num.col <- ntc0
					colv <- getcolors(num.col,col,col.opt)			
			
				}else{
					if( length(col < ntc0) ){
						col <- rep(col,ntc0)
					}
					col <- sapply(col, function(z){
						if(is.integer(z)){
							return(scales::alpha(z,1))	
						}else{
							return(z)
						}	
					})
					colv <- col[1:ntc0]	
					num.col <- ntc0^(spine | !eqwidth | (eqwidth & circular))
				}					
				col.Mat <- spread( t(rep(colv , nc/ntc)),nrow=nr)

				if(eqwidth & !res.val.only){
					
					dm <- dim(col.Mat)
					alpha <- spread(alpha, ncol = dm[2]/ncol(alpha))#num.col #ntc0
					
					col.Mat <- mapply( function(x,y){
						z <- col2rgb(x)/255
						if(nchar(x) == 9){
							alpha0 <- c(substr(x,8,8),substr(x,9,9))
							alpha0 <- match(alpha0,c(0:9,LETTERS[1:6]))-1
							if(!any(is.na(alpha0))){
								alpha0 <- (alpha0[1]*15 + alpha0[2])/240	
							}else{
								simpleWarning("could not identify the colors alpha value")
								alpha0 <- 1	
							}	
						}else{
							alpha0 <- 1	
						}
						rgb( red = z[[1]], green = z[[2]], blue = z[[3]], alpha = y*alpha0 )
					}, x = col.Mat, y = alpha)
					
					dim(col.Mat) <- dm
				}			
				col.Mat2 <- col.Mat
				
			}else{
				resid.mat[is.nan(resid.mat)] <- 0
				resid.mat[is.na(resid.mat)] <- 0
				if(is.null(resid.max)){
					resid.max <- max(abs(resid.mat),na.rm=T)	
				}
				resid.mat[resid.mat < -resid.max] <- -resid.max
				resid.mat[resid.mat > resid.max]  <-  resid.max
				rmax <- ifelse(abs(resid.max) < 10^{-6},1,resid.max)
				
				#if( !use.expected.values ){
				#	resid.mat <- -resid.mat	
				#}
		
				alf <- 0.05
				
				col.Mat <- apply(resid.mat,2,function(x){
						
						 if(resid.type == "pearson"){
							qx <- abs(x)
							
							qx = 2*(1 - sapply( qx, pnorm ))
							#print(qx)
							qx[qx < alf] <- alf
							
							val <- log( qx )/ log( alf )
							#print(val)
							return( hsv( h = 2/3*(sign(x) > 0), s = val, v = 0.7 ) )
						
						}else{
							qx <- abs(x/rmax)
							#rgb(sign(x) < 0,0,sign(x) > 0,alpha = exp(2*(qx-1))*base.alpha ) #qx # cut.rs/2
							qx <- qx - qx %% 1/(cut.rs-1)
							alf <- (exp(qx - 1)-exp(-1))/(1-exp(-1))
								 return( rgb(sign(x) < 0,0,sign(x) > 0,alpha = alf ) )#qx # cut.rs/2 # base.alp
						}
					
					})
				colv.grey <- hsv(h=1, s = 0, v = seq(0.9,0.6,-0.3/(ntc0-1))) #spine, eqwidth, circular?
				grey.Mat <- spread( t(rep(colv.grey , nc/ntc)),nrow=nr)
				

				if( both ){
					col.Mat2 <- col.Mat
					col.Mat3 <- grey.Mat
					col.Mat4 <- grey.Mat
					col.Mat2[resid.mat >= 0] <- hsv(0,0,0,alpha=0) # red remains
					col.Mat[resid.mat <= 0] <- hsv(0,0,0,alpha=0) # blue remains
					col.Mat4[resid.mat < 0] <- hsv(0,0,0,alpha=0) # grey for nonreds remains
					col.Mat3[resid.mat > 0] <- hsv(0,0,0,alpha=0) # grey for nonblues remains
				}else{
					col.Mat2 <- col.Mat	
					# col.Mat2[H==0] <- NA
				}
				
			}
				
# ----- plotting itself ----- #
# 	>>> divided into eqwidth <- T | eqwidth <- F >>> tfreq <- Id/sqrt/log
# 	>>> in spine mode, the target categories will be plotted one after another (see the for-loop)
# 	>>> the 3 diff. draw-function add the rectangles for the relative frequencies, the weights and the background respectively

# here were the transformations


if( circular ){
		
		P1 = 0*X
		X3 = X2 + 1/ncol(W)/2*(1-col.gap.prop)
		Y3 = Y2 + 1/nrow(W)/2*(1-row.gap.prop)
		
		#W2sqrt = tt3^(1/2) / max(  tt3^(1/2) )
		#W3 = spread(W2sqrt,ncol=ntc)
		
		arat <- ncol(W3)/nrow(W3)
	draw(H2/nrow(H2)*(1-row.gap.prop),H2/ncol(H2)*(1-col.gap.prop) + width.cor,X2,Y2,alpha = 1, bg = bgs, border = rect)
	
	
	if(prod(dim(H2)) > 1000){
		ang <- 	ang/2
	}
	if(prod(dim(H2)) > 10000){
		ang <- 	ang/2
	}
	
	if( disp.res ){ 
		H0 <- apply(H0,1:2, function(z) max(0.001,z))
		ratMat = sqrt(H00/H0)
		
		if( any( ratMat[H0 > 0.001] > max.rat ) ){
			warning(cat("residual with ratio", max(ratMat[H0 > 0.001],na.rm=TRUE), " > ",max.rat," occurred!"))	
		}
		ratMat <- apply(ratMat,1:2, function(z) min(max.rat, z))
		sc <- dim(W3)
		if(eqwidth){
			R1 = matrix(1,ncol = ncol(W3)*ntc0, nrow = nrow(W3)) # 0?, ncol, ??? 
			mv <- 1
		}else{
			R1 = spread(W3/2*(1-gap.prop),ncol=ntc0)# /ncol(W3)
			
			#maximal value/ratio?
			mv <- 1/2*(1-gap.prop)#/max(sc) #ncol
		}
	
			R2 <- R1 * ratMat
			rscf <- max(R2[ratMat <= max.rat], na.rm = TRUE) * 2/(1-gap.prop)# * max(sc) #ncol(W3)
			if(rscf < 1){
				rscf <- 1	
			}
			
			R2[is.na(R2)] <- 1/2*(1-gap.prop)/max(sc)#ncol
			R2 <- R2/rscf
			R1 <- R1/rscf
			
			R2 <- apply(R2, 1:2, function(z) min(z, mv))
			R1 <- apply(R1, 1:2, function(z) min(z, mv))
			
		for( i in 1:length(cat.ord) ){
				
				xind = seq( cat.ord[i], nc*ntc0, ntc0 )
				
				P2 = P1 + H0[,xind]

			#exp
			#gL1 <-
			 mapply( function(p1,p2,rad,bg,x,y){
				cc = seq( p1, p2, (p2-p1)/ceiling( (p2-p1)*ang ))
					grid.polygon(x = x+c(0,rad*cos(cc*2*pi)/sc[2]),y = y+c(0, rad*sin(cc*2*pi)/sc[1]),#ncol, nrow
					gp=gpar(fill=bg, lwd=1.5, col = line.col))
					}, p1 = P1, p2=P2, rad = R1[,xind], bg = col.Mat[,xind], x = X3, y = Y3 )
			# obs
			#gL2<-
			mapply( function(p1,p2,rad,bg,x,y){
					cc = seq( p1, p2, (p2-p1)/ceiling( (p2-p1)*ang ))
					grid.polygon(x = x+c(0,rad*cos(cc*2*pi)/sc[2]),y = y+c(0, rad*sin(cc*2*pi)/sc[1]),
						gp=gpar(fill=bg, col = line.col)
				)}, p1 = P1, p2=P2, rad = R2[,xind], bg = col.Mat2[,xind], x = X3, y = Y3 )
				
				#class(gL1) <- "gList"
				#class(gL2) <- "gList"
				
				#grid.draw(gL1)
				#grid.draw(gL2)
				
			if( both ){
			#exp	
			#gL3<-
				mapply( function(p1,p2,rad,bg,x,y){
					cc = seq( p1, p2, (p2-p1)/ceiling( (p2-p1)*ang ))
					grid.polygon(x = x+c(0,rad*cos(cc*2*pi)/sc[2]),y = y+c(0, rad*sin(cc*2*pi)/sc[1]),
						gp=gpar(fill=bg, lwd=1.5, col = line.col)
					)}, p1 = P1, p2=P2, rad = R1[,xind], bg = col.Mat3[,xind], x = X3, y = Y3 )
				
			# obs
			#gL4<-
				mapply( function(p1,p2,rad,bg,x,y){
					cc = seq( p1, p2, (p2-p1)/ceiling( (p2-p1)*ang ))
					grid.polygon(x = x+c(0,rad*cos(cc*2*pi)/sc[2]),y = y+c(0, rad*sin(cc*2*pi)/sc[1]),
						gp=gpar(fill=bg,col=line.col)
					)}, p1 = P1, p2=P2, rad = R2[,xind], bg = col.Mat4[,xind], x = X3, y = Y3 )
				
				#class(gL3) <- "gList"
				#class(gL4) <- "gList"
				#grid.draw(gL3)
				#grid.draw(gL4)
			}	
				
			
				P1 = P2
			}
		}else{
			#H0 <- apply(H0,1:2, function(z) max(0.0001,z))
			e1<-environment()
			wedge.id <- 0
			wc <- list()
				
			for( i in 1:length(cat.ord) ){
				xind = seq( cat.ord[i], nc*ntc0, ntc0 )
				P2 = P1 + H0[,xind]
				if( eqwidth ){	
					R = H2*(1-gap.prop)/2 #nrow
					sc <- dim(H2)
				}else{
					R = W3/2*(1-gap.prop) #ncol
					sc <- dim(W3)
				}
			
				
				mapply( function(p1,p2,rad,bg,x,y){
					p1 <- min(p1,1)
					p2 <- min(p2,1)
					p1 <- floor(10000*p1)/10000
					p2 <- floor(10000*p2)/10000
					if(p1 == p2){
						invisible(TRUE)
					}else{
						cc = seq( p1, p2, (p2-p1)/ceiling( (p2-p1)*ang ))
						e1$wedge.id <- e1$wedge.id+1
						if( abs(p2-p1) == 1 ){
							#grid.polygon(x = x+c(rad*cos(cc*2*pi)/sc[2]),y = y+c( rad*sin(cc*2*pi)/sc[1]),
							#gp=gpar(fill=bg, lwd=1.2,col=line.col))
							xv <- x+c(rad*cos(cc*2*pi)/sc[2])
							yv <- y+c( rad*sin(cc*2*pi)/sc[1])
							sp <- length(xv)
							idv <- rep(e1$wedge.id,sp)
							bgv <- rep(bg,sp)[1:sp]
							
						}else{
							#grid.polygon(x = x+c(0,rad*cos(cc*2*pi)/sc[2]),y = y+c(0, rad*sin(cc*2*pi)/sc[1]),
							#gp=gpar(fill=bg, lwd=1.2,col=line.col))
							xv <- x+c(0,rad*cos(cc*2*pi)/sc[2])
							yv <- y+c(0, rad*sin(cc*2*pi)/sc[1])
							sp <- length(xv)
							idv <- rep(e1$wedge.id,sp)
							bgv <- rep(bg,sp)[1:sp]
							
						}	
							e1$wc$x <- c(e1$wc$x,xv)
							e1$wc$y <- c(e1$wc$y,yv)
							e1$wc$bg <- c(e1$wc$bg,bgv)
							e1$wc$id <- c(e1$wc$id,idv)
						
						return(invisible(TRUE))
					}
				}, p1 = P1, p2=P2, rad = R, bg = col.Mat[,xind], x = X3, y = Y3 )
				P1 <- P2

			}
			#print(summary(as.data.frame(wc)))
			#print(table(wc$id,wc$bg))
			ubg <- unique(wc$bg)
			for(i in ubg){
				ii <- which(wc$bg == i)
				grid.polygon(wc$x[ii], wc$y[ii],id=wc$id[ii],gp=gpar(fill=wc$bg[ii], lwd=1.2,col=line.col))
			}
			
				
		}
	
}else{ # not circular
	
	if( disp.res ){
		H0 <- apply(H0,1:2, function(z) max(0.001,z))
		ratMat <- sqrt(H00/H0)
		ratMat <- apply(ratMat,1:2, function(z) min(max.rat, z))
	}
	if( disp.res & spine ){	
		if(eqwidth){
			R1 = matrix(1,ncol = ncol(W3)*ntc0/ntc, nrow = nrow(W3))/ncol(W3)/2*(1-gap.prop) # 0?
		}else{
			R1 = spread(W3/ncol(W3)/2*(1-gap.prop),ncol=ntc0/ntc)
		}
			
			R2 <- R1 * ratMat
			rscf <- max(R2[!is.na(R2)]) * (ncol(W3)*2/(1-gap.prop))
			rscf <- max(rscf,1)
	}else{
		rscf <- 1	
	}	
			#R2[is.na(R2)] <- 1/ncol(W3)/2*(1-gap.prop)
			#R2 <- R2/rscf
			#R1 <- R1/rscf

	col2Mat <- matrix(col2,nrow(H2),ncol(H2))
	col2Mat[H2==0] <- NA
	col2Mat[W2trans==0] <- NA
	draw(H2/nrow(H2)*(1-row.gap.prop),( W2trans/ncol(W2trans)*(1-col.gap.prop) + (round(W2trans,digits=7) > 0)*width.cor )/rscf,X2,Y2,alpha = 1, bg = col2Mat, border = ifelse(eqwidth, NA, line.col))
	draw(H2/nrow(H2)*(1-row.gap.prop),H2/ncol(H2)*(1-col.gap.prop) + width.cor,X2,Y2,alpha = 1, bg = bgs, border = rect)

	if( spine ){
		for( i in 1:length(cat.ord) ){
			if(i > 1){ 
				Y <- Y + H/nrow(H)*(1-row.gap.prop)
			}
			xind <- seq( cat.ord[i], nc*ntc0, ntc0 )
			H <- as.matrix(H0[,xind])
			if( nr < 2 ){
				H <- t(H)	
			}
			CM1 <- col.Mat[,xind]
			CM1[H==0] <- NA
			draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop)/rscf,X3,Y,alpha = 1, bg = CM1, border = NA)
			if( disp.res ){
				#ratMat = sqrt(H00/H0)
				max.rat = max(ratMat)
						
				W3res = W3/ncol(W3)*(1-col.gap.prop)*ratMat[,xind]
				# obs
				CM2 <- col.Mat2[,xind]
				CM2[H==0] <- NA
				draw(H/nrow(H)*(1-row.gap.prop),W3res/rscf,X3,Y,alpha = 1, bg = CM2, border = NA)
					
				if( both ){
						#exp	
						draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop)/rscf,X3,Y,alpha = 1, bg = col.Mat3[,xind], border = NA)
						# obs
						draw(H/nrow(H)*(1-row.gap.prop),W3res/rscf,X3,Y,alpha = 1, bg = col.Mat4[,xind], border = NA)
				}	
			}
			draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop)/rscf,X3,Y,alpha = 1, bg = NA, lwd = 1.5, border = line.col)
			
		}
		Y <- Y0
	
	}else{
		
		H <- apply(H,c(1,2),function(x) return(min(x/max.scale,1)))
		
		if( disp.res ){
				draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = col.Mat, border = NA)
				
				ratMat = sqrt(H00/H0)
				max.rat = max(ratMat)		
				
				W3res = W3/ncol(W3)*(1-col.gap.prop)*ratMat#[,xind]
				# obs
				draw(H00/nrow(H00)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = col.Mat2, border = NA)
					
				if( both ){
						#exp	
						draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = col.Mat3, border = NA)
						# obs
						draw(H00/nrow(H00)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = col.Mat4, border = NA)
				}
				draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = NA, lwd =1.5, border = line.col)	
		}else{
			col.Mat[H==0] <- NA
			draw(H/nrow(H)*(1-row.gap.prop),W3/ncol(W3)*(1-col.gap.prop),X3,Y,alpha = 1, bg = col.Mat, border = line.col)	
		}#border = line.col ???
	}
}
		
	if(yaxis & !circular){
		if(spine){
			rscf <- 1	
		}
		sapply(y,function(z){
			grid.yaxis(at =seq(z,z+1/length(y)*(1-row.gap.prop),1/length(y)*(1-row.gap.prop)/5), label = round(seq(0,max.scale,max.scale/5)/rscf,digits = 3), main = TRUE,
			edits = NULL, name = NULL,
			gp = gpar(fontsize = 10*lab.cex), draw = TRUE, vp = NULL)
		})
		sapply(y,function(z){
			grid.yaxis(at =seq(z,z+1/length(y)*(1-row.gap.prop),1/length(y)*(1-row.gap.prop)/5), label = round(seq(0,max.scale,max.scale/5)/rscf,digits = 3), main = FALSE,
			edits = NULL, name = NULL,
			gp = gpar(fontsize=10*lab.cex), draw = TRUE, vp = NULL)
		})
	}
	
	upViewport()
# ----- ---------------- ----- #	
# ----- labeling section ----- #	
# ----- ---------------- ----- #		
	

	
	
	if(lab){
		lab.cex <- rep(lab.cex,nv)[1:nv]
			if(any(col.vars*label)){
	# ----- labeling the x-axis ----- #	

	vpX <- viewport(x = 1-yaxis*s1 - s3*(!is.null(expected)), y = 1-ncl*s2, width = 1-nrl*s2 - yaxis*2*s1 - s3*(!is.null(expected)), height = ncl*s2, just = c("right", "bottom"), name = "vpX")
	pushViewport(vpX)
	cur <- 1
	zc <- cumprod(col.nlvl)[ which(label[which(col.vars)]) ]
    for(i in 1:ncl){
		vpXX <- viewport(x = 1, y = 1 - i*1/ncl, width = 1, height = 1/ncl, just = c("right", "bottom"), name = "vpXX")
		pushViewport(vpXX)
		
		if( (i > ncl-1) & num.mode.x ){
			min.tick <- floor(1000*min(data[,ind[ cind[i]] ],na.rm=T ))/1000
			max.tick <- floor(1000*max(data[,ind[ cind[i]] ],na.rm=T ))/1000
			
			#vpXXb <- viewport(x = 0.5, y = 0, width = 1, height = 0.1, just = "centre", name = "vpXXb")
		#pushViewport(vpXXb)
			my.grid.axis(x0=0,y0=0.1,len=1,ticks= signif(seq(min.tick,max.tick,(max.tick - min.tick)/cut[cind[i]]),5), rot=0, lab.cex = lab.cex[cind[i]])
			
			#grid.xaxis(at = seq(0,1,1/cut[cind[i]]), label = round(seq(min.tick,max.tick,(max.tick - min.tick)/cut[cind[i]]),digits = 3), main = FALSE,
			#edits = NULL, name = NULL,
			#gp = gpar(), draw = TRUE, vp = NULL)
			#popViewport()
			grid.text(  names(dset)[cind[i]],x = 0.5 , y = 5/6, just = "centre",gp=gpar(fontface = "bold",cex=1.1*lab.cex))	
		}else{
			if( varnames ){
				grid.text(  names(dset)[cind[i]],x = 0.5 , y = 5/6, just = "centre",gp=gpar(fontface = "bold",cex=1.1*lab.cex))	
			}
			#z <- (nlvl[cind])[i]*cur
			z <- zc[i]
			grid.text(  rep( clabs[[i]] ,cur), x = X[1,seq(1,nc,nc/z)] + (X[1,seq(nc/z,nc,nc/z)] - X[1,seq(1,nc,nc/z)] +1/nc*(1-col.gap.prop))/2, y = 1/3, just = "centre",gp=gpar(cex=lab.cex[ cind[i] ])) 
			if( boxes ){
				grid.rect(  y = 1/3 , x = X[1,seq(1,nc,nc/z)] , just = c("left","centre"), width = X[1,seq(nc/z,nc,nc/z)] - X[1,seq(1,nc,nc/z)] +1/nc*(1-col.gap.prop), height = 1/2, gp = gpar(fill="transparent"))
			}
			cur <- z
		}
		popViewport()
	}
		
	upViewport()
	}
		if(any( (!col.vars)*label)){
# ----- labeling the y-axis ----- #			
	vpY <- viewport(x = nrl*s2, y = 0, width = nrl*s2, height = 1-ncl*s2, just = c("right", "bottom"), name = "vpY")
	pushViewport(vpY)
	cur <- 1
	zr <- cumprod(row.nlvl)[ which(label[which(!col.vars)]) ]
	
	for(i in 1:nrl){
		
		vpYY <- viewport(x = 0 + (i-1)*1/nrl, y = 0, width = 1/nrl, height = 1, just = c("left", "bottom"), name = "vpYY")
		pushViewport(vpYY)
		if( (i > nrl-1) & num.mode.y ){
			min.tick <- floor(1000*min(data[,ind[ rind[i]] ],na.rm=T ))/1000
			max.tick <- floor(1000*max(data[,ind[ rind[i]] ],na.rm=T ))/1000
			
			#vpYYb <- viewport(x = 1, y = 0.5, width = 0.1, height = 1, just = "centre", name = "vpYYb")
		#pushViewport(vpYYb)
			
			#grid.yaxis(at = seq(0,1,1/cut[rind[i]]), label = round(seq(min.tick,max.tick,(max.tick - min.tick)/cut[rind[i]]),digits = 3), main = TRUE,
			#edits = NULL, name = NULL,
			#gp = gpar(), draw = TRUE, vp = NULL)
			my.grid.axis(x0=0.9,y0=0,len=1,ticks=signif(seq(min.tick,max.tick,(max.tick - min.tick)/cut[rind[i]]),5), rot=90, lab.cex = lab.cex[rind[i]])
			
			#popViewport()
				grid.text(  names(dset)[rind[i]],x = 1/6 , y = 0.5, just = "centre",rot=90,gp=gpar(fontface = "bold",cex=1.1*lab.cex))	
		}else{
		
			if( varnames & (nr > 1) ){
				grid.text(  names(dset)[rind[i]],x = 1/6 , y = 0.5, just = "centre",rot=90,gp=gpar(fontface = "bold",cex=1.1*lab.cex))	
			}			
			#z <- (nlvl[rind])[i]*cur
			z <- zr[i]
			grid.text(  rep( rlabs[[i]] ,cur), y = Y[seq(1,nr,nr/z),1] + (Y[seq(nr/z,nr,nr/z),1] - Y[seq(1,nr,nr/z),1] +1/nr*(1-row.gap.prop))/2, x = 2/3, just = "centre",rot=90,gp=gpar(cex=lab.cex[ rind[i] ])) 
			if( boxes ){
				grid.rect(  x = 2/3 , y = Y[seq(1,nr,nr/z),1] , just = c("centre","bottom"), height = Y[seq(nr/z,nr,nr/z),1] - Y[seq(1,nr,nr/z),1] +1/nr*(1-row.gap.prop), width = 1/2, gp = gpar(fill="transparent"))
			}
					
			cur <- z
		}
		popViewport()
	}
	upViewport()
	}
	}
	# ----- scale for expected mode ----- #
	if(!is.null(expected) & !res.val.only){

		vpS <- viewport(x = 1, y = 0 , width = s3, height = 0.8, just = c("right","bottom"), name = "vpS")
		pushViewport(vpS)
		
		#alpha values for the residuals on a log-scale
		#alfav <- exp(seq(0,-cut.rs/2+0.5,-0.5))
		if( resid.type == "pearson" ){
			sval <- c(log( (1 - seq(1-alf, alf, -(1-2*alf)/(cut.rs-1) ) ) )/ log( alf ))
			sval <- c(sval,0, rev(sval))
			
			alfy <- qnorm(1-alf/2)
			minv <- min(-alfy,min(resid.mat[resid.mat < 0]))
			maxv <- max(alfy,max(resid.mat[resid.mat > 0]))
			ys <- c( minv, seq(-alfy,alfy,2*alfy/(2*cut.rs-1)), maxv)
			colv <- hsv( c(rep(0,cut.rs),0,rep(2/3,cut.rs)),  sval,  0.8)
			
			grid.rect(  x = 0.1 , y = ((ys - minv)/(maxv - minv))[-length(ys)], just = c("left","bottom"), height = diff(ys)/(maxv-minv), width = 0.3 , gp = gpar(fill=colv))
			grid.text( label=round(ys,digits=2), y = (ys - minv)/(maxv - minv), x = 0.45, just = c("left","centre"),rot=0,gp=gpar(cex=lab.cex))
		
			
		}else{
			alfav <- (exp(seq(0,-1,-1/(cut.rs-1)))-exp(-1))/(1-exp(-1))
			colv <- rgb(rep(c(1,0),each=cut.rs),0,rep(c(0,1),each=cut.rs),alpha = c(alfav,rev(alfav)))#*base.alpha
			grid.rect(  x = 0.1 , y = seq(0,1-1/2/cut.rs,1/2/cut.rs) , just = c("left","bottom"), height = 1/2/cut.rs, width = 0.3 , gp = gpar(fill=colv))
			grid.text( label=round(seq(-resid.max,resid.max,resid.max/min(10,cut.rs)),digits=2), y = seq(0,1,1/2/min(10,cut.rs)), x = 0.45, just = c("left","centre"),rot=0,gp=gpar(cex=lab.cex))
		}
		#grid.rect(  x = 0.1 , y = 0, just = c("left","bottom"), height = 1, width = 0.3 , gp = gpar(fill=hsv(0,0,0,alpha=0.1)))
			
		upViewport()	
		vpS2 <- viewport(x = 1, y = 0.8  , width = 0.1, height = 0.2 , just = c("right","bottom"), name = "vpS2")
		pushViewport(vpS2)
		
		grid.text( label=resid.type, y = 0.5, x = 0.2, just = "centre",rot=90,gp=gpar(cex=lab.cex*1.2))
		grid.text( label="residuals", y = 0.5, x = 0.4, just = "centre",rot=90,gp=gpar(cex=lab.cex*1.2))
		grid.text( label=paste("p.value =",format(p.value,digits=2)), y = 0.5, x = 0.7, just = "centre",rot=90,gp=gpar(cex=lab.cex*1.1))
		upViewport()
	}
	
	#op <- options()
	#options(show.error.messages = FALSE)
	#try( upViewport(), silent = TRUE)
	upViewport()
	if(!is.null(vp)){
		try( upViewport(), silent = TRUE )
	}
	#options(op)
	
	return(invisible(TRUE))
}

# ------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------ #


space <- function(nlvl,col.vars,gap.mult = 1.1,gap.prop = 0.1,last.col.zero = F, last.row.zero = F){
	# 	>>> function to compute the gaps widths for the rmb plot
	# 	>>> might be extended in further versions of the package
	
	# ----- parameter check ----- #
		if( !(  (gap.prop > 0)&(gap.prop < 1) ) ){
			stop(simpleError("Wrong gap.prop specification!"))
		}

	col.nlvl <- nlvl[col.vars]
	row.nlvl <- nlvl[!col.vars]								

	nc <- prod(col.nlvl)
	nr <- prod(row.nlvl)

	ncv <- length(col.nlvl)
	nrv <- length(row.nlvl)

	col.spaces <- rep(gap.prop,nc-1)
	row.spaces <- rep(gap.prop,nr-1)

	if( ncv > 1 ){
		cind <- 1:(nc-1)
		cur <- 1
		for( i in (ncv-1):1 ){ 
			cur <- cur*col.nlvl[i+1]
			curind <- which( (cind %% cur) == 0)
			col.spaces[ curind  ] <- col.spaces[ curind  ] * gap.mult
		}
	}
	if( nrv > 1 ){
		rind <- 1:(nr-1)
		cur<-1
		for( i in (nrv-1):1 ){ 
			cur <- cur*row.nlvl[i+1]
			curind <- which( (rind %% cur) == 0)
			row.spaces[ curind  ] <- row.spaces[ curind  ] * gap.mult
		}
	}
	if(last.col.zero){
		col.spaces[which(col.spaces == gap.prop)] <- 0
		col.spaces <- col.spaces/gap.prop
	}
	if(last.row.zero){
		row.spaces[which(row.spaces == gap.prop)] <- 0
		row.spaces <- row.spaces/gap.prop
	}
	
	return( list(col.spaces,row.spaces))
}

# ------------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------------------ #


spread <- function(M,ncol = 1,nrow = 1){
	# 	>>> this function turns each matrix entry into a nrow x ncol - matrix
	M <- as.matrix(M)
	M2 <- apply(M,2,function(x) rep(x,each = nrow)  )
	if( (nrow(M) == 1)&(nrow == 1) ){
		M2 <- matrix(M2,ncol = ncol(M))
	#dim(M2) <- dim(M)
	}
	
	M3 <- t(apply(t(M2),2,function(x) rep(x,each = ncol)  ))
	if( (ncol(M) == 1)&(ncol == 1) ){
		M3 <- matrix(M3,nrow = nrow(M2))
		#dim(M3) <- dim(M2)
	}
	return(M3)
}


prettyscale = function(breaks){
	#a <- min(breaks)
	#b <- max(breaks)
	#n <- length(breaks)
	
	#ep <- floor(log(abs(breaks))/log(10))
	ep <- sapply(breaks, function(z){
		if( abs(z) < 1 ){
			floor( log(abs(1/z))/log(10) )	
		}else{
			floor( log(abs(z))/log(10) )	
		}	
	})
	ep[abs(breaks)<1] <- - ep[abs(breaks)<1]
	esigns = c("+","-")[(sign(ep)<0)+1]
	breaks2 <- breaks
	small <- which(abs(ep) < 3)
	big <- which(abs(ep) > 2)
	a <- format(round(breaks2[small], 2),nsmall=2)
	b <- paste(round(breaks2[big]/(10^ep[big]), digits=1),paste("e",esigns[big],format(paste(0,abs(ep)[big],sep="")),sep=""),sep="")
	
	breaks2[big] <- b
	breaks2[small] <- a
	return(breaks2)	
}

my.grid.axis = function(x0,y0,rot, ticks, len = 1, ltm = 1/30, clockwise = FALSE, lab.cex = 0.8, keep=7, col.axis = 1, trot = 320){
		
	grid.lines(x=c(x0, x0+len*cos(rot/180*pi)), y = c(y0, y0+len*sin(rot/180*pi)),gp=gpar(col=col.axis))
	n <- length(ticks)-1
	
	sgn <- ifelse(clockwise,-1,1)
	
	x.ticks.1 <- x0+seq(0,1,1/n)*len*cos(rot/180*pi)
	y.ticks.1 <- y0+seq(0,1,1/n)*len*sin(rot/180*pi)
	
	x.ticks.2 <- x.ticks.1 + sgn*ltm*len*cos(rot/180*pi+pi/2)
	y.ticks.2 <- y.ticks.1 + sgn*ltm*len*sin(rot/180*pi+pi/2)
	ss <- seq(1,n,n/floor(keep*4 -3))
	mapply( function(x0,y0,x1,y1){
		grid.lines(c(x0,x1),c(y0,y1),gp=gpar(col=col.axis))
		},x0 = x.ticks.1[ss], x1 = x.ticks.2[ss], y0= y.ticks.1[ss], y1 = y.ticks.2[ss])
		if(rot <= 135 & rot > 45) just = ifelse(clockwise,"left","right")
		if(rot <= 225 & rot > 135) just = ifelse(clockwise,"bottom","top")
		if(rot <= 315 & rot > 225) just = ifelse(clockwise,"right","left")
		if(rot <= 45 | rot > 315) just = ifelse(clockwise,"left","right")#top/bottom
	
#	ii <- round( (n-1)/4 )
#	if(ii > 0){
#		ids <- which( floor(1:(n-2)+ii/2) %% ii != 0)
#		ticks[ids+1] = ""	
#	}
	
#	tmp <- n / c(5,6)
#	tick.step <- trunc(tmp)[which.min( tmp - trunc(tmp)  )]
if(n > keep){
	keep <- round(seq(1,n,n/(keep-1)))
	ticks[-c(1,n+1,keep)] = ""
}		
	grid.text(ticks, x.ticks.2, y.ticks.2, just = just, rot = trot, gp=gpar(fontsize=lab.cex*10,col=col.axis))
}

