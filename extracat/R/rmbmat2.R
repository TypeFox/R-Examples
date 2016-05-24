
#############################################
#### ----- << rmb matrix design >> ----- ####
#############################################






rmbmat = function(x, tv, cut = 20, freqvar = NULL, plot.tv = FALSE,
 num.mode = TRUE, mode ="circular", eqwidth = FALSE, freq.trans="sqrt", innerval = 1, allocation = I,
max.scale = 1, use.na = FALSE, expected = FALSE, model.opt = list(), 
    gap.prop = 0.2, gap.mult = 1.5, col = "hcl", col.opt = list(), 
    label = FALSE, label.opt = list(), diag.opt = list(), lower.opt = list(), upper.opt = list(), rc.opt = list(), factor.opt = list(),...){

spine <- FALSE
circular <- FALSE
if(mode %in% c("circular", "circ", "c", "circles", "pie", "p", "piecharts")){
	circular <- TRUE	
}
if(mode %in% c("spine", "s")){
	spine <- TRUE	
}
if(mode %in% c("rect", "rectangles", "nested.rectangles", "r", "nr")){
	rectangular <- TRUE	
}
if(mode %in% c("nested.circles", "nc", "ncirc")){
	nested.circles <- TRUE	
}
if(is.null(allocation)){
	allocation <- function(s) 1	
}

x <- as.data.frame(x)


if( "Freq" %in% names(x) ){
	freqvar <- "Freq"	
}
if(!is.null(freqvar)){
	fi <- which(names(x)==freqvar)
	nv <- ncol(x) - 1 - !plot.tv	
}else{
	nv <- ncol(x) - !plot.tv	
	fi <- nv+3	
}

if("tv2" %in% names(lower.opt)){
	nv <- ifelse(plot.tv,nv,nv-1)
	tv2 <- lower.opt$tv2
	tv2.lower <- TRUE
}else{
	tv2.lower <- FALSE
	if("tv2" %in% names(upper.opt)){
		nv <- ifelse(plot.tv,nv,nv-1)
		tv2 <- upper.opt$tv2
	}else{
		tv2 <- NA	
	}	
}

###########################################
#### ----- setup option defaults ----- ####
###########################################
if(! "bg" %in% names(col.opt) ){
	if(!circular){
		col.opt$bg <- NA	
	}	
}
if(! "line.col" %in% names(col.opt) ){
	col.opt$line.col <- NA	
}
if(! "lab.cex" %in% names(label.opt) ){
	label.opt$lab.cex <- sqrt(2/nv)	
}
if(! "yaxis" %in% names(label.opt) ){
	label.opt$yaxis <- FALSE	
}

###########################################
#### ----- setup option matrices ----- ####
###########################################

spine <- matrix(spine,ncol = nv, nrow=nv)
circular <- matrix(circular , ncol = nv,  nrow=nv)
eqwidth <- matrix(eqwidth , ncol = nv,  nrow=nv)
col <- matrix(col , ncol = nv,  nrow=nv)
label <- matrix(label , ncol = nv,  nrow=nv)

freq.trans <- rep(list(rep(list(freq.trans),nv)),nv)
col.opt <- rep(list(rep(list(col.opt),nv)),nv)
label.opt <- rep(list(rep(list(label.opt),nv)),nv)
model.opt <- rep(list(rep(list(model.opt),nv)),nv)

num.mode <- matrix(num.mode , ncol = nv,  nrow=nv)
max.scale <- matrix(max.scale , ncol = nv,  nrow=nv)
use.na <- matrix(use.na , ncol = nv,  nrow=nv)
if(!is.null(expected)){
	if(!is.logical(expected)){
		expected <- TRUE
	}	
}
expected <- matrix(expected , ncol = nv,  nrow=nv)
gap.prop <- matrix(gap.prop , ncol = nv,  nrow=nv)
gap.mult <- matrix(gap.mult , ncol = nv,  nrow=nv)
# <- matrix( , ncol = nv,  nrow=nv)


# excluded:
# cat.ord, expected = list(), residuals, col.vars



# ----- diag.opt ----- #

#if( "spine" %in% names(diag.opt) ) diag(spine) <- diag.opt$spine
#if( "circular" %in% names(diag.opt) ) diag(circular) <- diag.opt$circular
if( "mode" %in% names(diag.opt) ){
	if(diag.opt$mode %in% c("circular", "circ", "c", "circles", "pie", "p", "piecharts")){
		diag(circular) <- TRUE	
		diag(spine) <- FALSE
	}
	if(diag.opt$mode %in% c("spine", "s")){
		diag(spine) <- TRUE	
		diag(circular) <- FALSE
	}
	if(diag.opt$mode %in% c("bars", "bar", "b")){
		diag(circular) <- FALSE
		diag(spine) <- FALSE		
	}
}
if( "eqwidth" %in% names(diag.opt) ) diag(eqwidth) <- diag.opt$eqwidth
if( "col" %in% names(diag.opt) ) diag(col) <- diag.opt$col
if( "num.mode" %in% names(diag.opt) ) diag(num.mode) <- diag.opt$num.mode
if( "max.scale" %in% names(diag.opt) ) diag(max.scale) <- diag.opt$max.scale
if( "use.na" %in% names(diag.opt) ) diag(use.na) <- diag.opt$use.na

if( "expected" %in% names(diag.opt) ){
	if(!is.logical(diag.opt$expected)){
		diag.opt$expected <- TRUE
	}	
	diag(expected) <- diag.opt$expected
}

if( "gap.prop" %in% names(diag.opt) ) diag(gap.prop) <- diag.opt$gap.prop
if( "gap.mult" %in% names(diag.opt) ) diag(gap.mult) <- diag.opt$gap.mult
if( "label" %in% names(diag.opt) ) diag(label) <- diag.opt$label

if( "freq.trans" %in% names(diag.opt) ){
	for(i in 1:(nv)) freq.trans[[i]][[i]] <- diag.opt$freq.trans	
}
if( "col.opt" %in% names(diag.opt) ){
	for(i in 1:(nv)) col.opt[[i]][[i]] <- diag.opt$col.opt	
}
if( "label.opt" %in% names(diag.opt) ){
	for(i in 1:(nv)) label.opt[[i]][[i]] <- diag.opt$label.opt	
}
if( "model.opt" %in% names(diag.opt) ){
	for(i in 1:(nv)) model.opt[[i]][[i]] <- diag.opt$model.opt	
}


# ----- lower.opt and upper.opt ----- #

ltr <- lower.tri(spine)
#if( "spine" %in% names(lower.opt) ) spine[ltr] <- lower.opt$spine
#if( "circular" %in% names(lower.opt) ) circular[ltr] <- lower.opt$circular
if( "mode" %in% names(lower.opt) ){
	if(lower.opt$mode %in% c("circular", "circ", "c", "circles", "pie", "p", "piecharts")){
		circular[ltr] <- TRUE	
		spine[ltr] <- FALSE
	}
	if(lower.opt$mode %in% c("spine", "s")){
		spine[ltr] <- TRUE	
		circular[ltr] <- FALSE
	}
	if(lower.opt$mode %in% c("bars", "bar", "b")){
		circular[ltr] <- FALSE
		spine[ltr] <- FALSE		
	}
}
if( "eqwidth" %in% names(lower.opt) ) eqwidth[ltr] <- lower.opt$eqwidth
if( "col" %in% names(lower.opt) ) col[ltr] <- lower.opt$col
if( "num.mode" %in% names(lower.opt) ) num.mode[ltr] <- lower.opt$num.mode
if( "max.scale" %in% names(lower.opt) ) max.scale[ltr] <- lower.opt$max.scale
if( "use.na" %in% names(lower.opt) ) use.na[ltr] <- lower.opt$use.na

if( "expected" %in% names(lower.opt) ){
	if(!is.logical(lower.opt$expected)){
		lower.opt$expected <- TRUE
	}	
	expected[ltr] <- lower.opt$expected
}

if( "gap.prop" %in% names(lower.opt) ) gap.prop[ltr] <- lower.opt$gap.prop
if( "gap.mult" %in% names(lower.opt) ) gap.mult[ltr] <- lower.opt$gap.mult
if( "label" %in% names(lower.opt) ) label[ltr] <- lower.opt$label

if( "freq.trans" %in% names(lower.opt) ){
	for(i in 2:(nv)){
		for(j in 1:(i-1)){
			freq.trans[[i]][[j]] <- lower.opt$freq.trans	
		}
	}
}
if( "col.opt" %in% names(lower.opt) ){
	for(i in 2:(nv)){
		for(j in 1:(i-1)){
			col.opt[[i]][[j]] <- lower.opt$col.opt
		}
	}	
}
if( "label.opt" %in% names(lower.opt) ){
	for(i in 2:(nv)){
		for(j in 1:(i-1)){
			label.opt[[i]][[j]] <- lower.opt$label.opt	
		}
	}	
}
if( "model.opt" %in% names(lower.opt) ){
	for(i in 2:(nv)){
		for(j in 1:(i-1)){
			model.opt[[i]][[j]] <- lower.opt$model.opt
		}
	}	
}

utr <- upper.tri(spine)
#if( "spine" %in% names(upper.opt) ) spine[utr] <- upper.opt$spine
#if( "circular" %in% names(upper.opt) ) circular[utr] <- upper.opt$circular
if( "mode" %in% names(upper.opt) ){
	if(upper.opt$mode %in% c("circular", "circ", "c", "circles", "pie", "p", "piecharts")){
		circular[utr] <- TRUE
		spine[utr] <- FALSE	
	}
	if(upper.opt$mode %in% c("spine", "s")){
		spine[utr] <- TRUE	
		circular[utr] <- FALSE
	}
	if(upper.opt$mode %in% c("bars", "bar", "b")){
		circular[utr] <- FALSE
		spine[utr] <- FALSE		
	}
}
if( "eqwidth" %in% names(upper.opt) ) eqwidth[utr] <- upper.opt$eqwidth
if( "col" %in% names(upper.opt) ) col[utr] <- upper.opt$col
if( "num.mode" %in% names(upper.opt) ) num.mode[utr] <- upper.opt$num.mode
if( "max.scale" %in% names(upper.opt) ) max.scale[utr] <- upper.opt$max.scale
if( "use.na" %in% names(upper.opt) ) use.na[utr] <- upper.opt$use.na

if( "expected" %in% names(upper.opt) ){
	if(!is.logical(upper.opt$expected)){
		upper.opt$expected <- TRUE
	}	
	expected[utr] <- upper.opt$expected
}

if( "gap.prop" %in% names(upper.opt) ) gap.prop[utr] <- upper.opt$gap.prop
if( "gap.mult" %in% names(upper.opt) ) gap.mult[utr] <- upper.opt$gap.mult
if( "label" %in% names(upper.opt) ) label[utr] <- upper.opt$label

if( "freq.trans" %in% names(upper.opt) ){
	for(i in 2:(nv)){
		for(j in 1:(i-1)){
			freq.trans[[i]][[j]] <- upper.opt$freq.trans	
		}
	}
}
if( "col.opt" %in% names(upper.opt) ){
	for(i in 1:(nv-1)){
		for(j in (i+1):(nv)){
			col.opt[[i]][[j]] <- upper.opt$col.opt
		}
	}	
}
if( "label.opt" %in% names(upper.opt) ){
	for(i in 1:(nv-1)){
		for(j in (i+1):(nv)){
			label.opt[[i]][[j]] <- upper.opt$label.opt	
		}
	}	
}
if( "model.opt" %in% names(upper.opt) ){
	for(i in 1:(nv-1)){
		for(j in (i+1):(nv)){
			model.opt[[i]][[j]] <- upper.opt$model.opt
		}
	}	
}


# ----- rows and cells ----- #


specced <- names(rc.opt)
id.list <- lapply(specced, function(z){
	rs <- regexpr(pattern="r",text=z)
	cs <- regexpr(pattern="c",text=z)
	zlen <- nchar(z)
	row.id = NA
	col.id = NA
	if(rs > 0){
		if(rs > cs){
			suppressWarnings(row.id <- as.integer(substr(z,rs+1,zlen)))	
			if( cs > 0 ){
				suppressWarnings(col.id <- as.integer(substr(z,cs+1,rs-1)))	
			}
		}else{
			suppressWarnings(col.id <- as.integer(substr(z,cs+1,zlen)))	
			suppressWarnings(row.id <- as.integer(substr(z,rs+1,cs-1)))	
		}	
	}else{
		if( cs > 0 ){
			suppressWarnings(col.id <- as.integer(substr(z,cs+1,zlen)))	
		}	
	}
	return(c(row.id, col.id))
})
#print(id.list)
wrong.spec <- sapply(id.list, function(z) sum(is.na(z))==2)
#print(wrong.spec)
#if(length(wrong.spec) > 0){
if(!all(wrong.spec	== FALSE)){
	cat("Incorrect specification: rc.opt arguments ",specced[which(wrong.spec)]," ignored.")

if(any(wrong.spec)){
	id.list <- id.list[-which(wrong.spec)]
	rc.opt <- rc.opt[-which(wrong.spec)]
}
}
env <- environment()
nix <- mapply(function(x,y){
	applyspecs(y,env,x)
	return(invisible(TRUE))
	}, x = id.list, y = rc.opt)
rm(nix)


# ----------------------------------------------------------------------------------------------- #

data <- idata.frame(x)

# (fi always exists)
if(plot.tv){
	tvfi <- fi
}else{
	tvfi <- c(fi,tv,ifelse(is.na(tv2),nv+2,tv2))
}
fterms <- names(data)[-tvfi]	

if( length(cut) == 1){
	cut <- rep(cut, nv )	
}

cuttable <- !sapply(x[,-tvfi], is.factor)#c(tv,fi)

will.be.cut <- as.logical( cuttable * (cut > 0))

nlvl <- sapply(x[,-tvfi], function(z) nlevels(as.factor(z)) )#c(tv,fi)
nlvl[will.be.cut] <- cut[will.be.cut]

wh <- (wh<-allocation(nlvl))/sum(wh)

mat.layout <- grid.layout(nrow = nv , ncol = nv , widths = wh, heights=wh)
grid.newpage()

# ----------------------------------------------------------------------------------------------- #


# ----- categorical ----- #

factor.ids <- which(sapply(x,is.factor)[-tvfi] | !will.be.cut  )

factor.ids <- expand.grid(factor.ids,factor.ids)
if(length(factor.opt) > 0){
nix <- apply(factor.ids, 1, function(z){
	applyspecs(factor.opt,env,z)
	return(invisible(TRUE))
})
rm(nix)
}
if("mosaic" %in% names(factor.opt)){
    #require(vcd)
	print("mosaics from vcd not yet implemented")
	mosaic <- matrix(FALSE, ncol= nv, nrow = nv)
	nix <- apply(factor.ids, 1, function(z){
		env$mosaic[z[1],z[2]] <- TRUE
	return(invisible(TRUE))
	})
	rm(nix)	
}

# ----------------------------------------------------------------------------------------------- #


ids <- seq_along(fterms)
vp.base <- viewport(x=0.5,y=0.5,width=0.96,height=0.96)
pushViewport(vp.base)
vp.mat <- viewport(layout = mat.layout)
pushViewport(vp.mat)
sapply(ids,function(i){
	sapply(ids,function(j){
		if(i != j){
			use.dset <- FALSE
			if(!is.na(tv2) & ((tv2.lower & i > j) | (!tv2.lower & i < j))){
				if(plot.tv){
				# check whether tv2 equals i or j
					if( tv2 == j ){
						ftj <- paste(fterms[j],".1",sep="")
						use.dset <- TRUE
						if(tv2 == i ){
							fti <- paste(fterms[i],".2",sep="")
						}else{
							fti <- fterms[i]
						}
						dset <- x[, match(c(fterms[i],fterms[j],names(x)[tv2]), names(x))]
						if(fi < nv+3){
							dset$Freq <- x[,fi]	
						}
					}else{
						ftj <- fterms[j]
						if(tv2 == i ){
							fti <- paste(fterms[i],".1",sep="")
							use.dset <- TRUE
							dset <- x[, match(c(fterms[i],fterms[j],names(x)[tv2]), names(x))]
							if(fi < nv+3){
								dset$Freq <- x[,fi]	
							}
						}else{
							fti <- fterms[i]
						}	
					}
				}else{
					ftj <- fterms[j]
					fti <- fterms[i]
				}
					f <- paste("~",ftj,"+",fti,"+",names(x)[tv2],sep="")
			}else{
				if(plot.tv){
				# check whether tv equals i or j
					if( tv == j ){
						ftj <- paste(fterms[j],".1",sep="")
						use.dset <- TRUE
						if(tv == i ){
							fti <- paste(fterms[i],".2",sep="")
						}else{
							fti <- fterms[i]
						}
						dset <- x[, match(c(fterms[i],fterms[j],names(x)[tv]), names(x))]
						if(fi < nv+3){
							dset$Freq <- x[,fi]	
						}
					}else{
						ftj <- fterms[j]
						if(tv == i ){
							fti <- paste(fterms[i],".1",sep="")
							use.dset <- TRUE
							dset <- x[, match(c(fterms[i],fterms[j],names(x)[tv]), names(x))]
							if(fi < nv+3){
								dset$Freq <- x[,fi]	
							}
						}else{
							fti <- fterms[i]
						}	
					}
				}else{
					ftj <- fterms[j]
					fti <- fterms[i]
				}
					f <- paste("~",ftj,"+",fti,"+",names(x)[tv],sep="")
			}
			#print(f)
			#print(c(i,j,tv,tv2))
			#print("-----------")
			vp1 <- viewport(x=0.5,y=0.5,width=0.96,height=0.96)
			vp0 <- viewport(layout.pos.row = i, layout.pos.col = j)
			pushViewport(vp0)
			#grid.rect()
			#pushViewport(vp1)
			if( use.dset ){
				rmb.formula(as.formula(f), data=dset, cut = cut[c(i,j,tv)], col.vars = NULL, spine = spine[i, j], circular = circular[i, j], 
					eqwidth = eqwidth[i, j], cat.ord = NULL, freq.trans = freq.trans[[i]][[j]], 
					num.mode = num.mode[i, j], innerval = innerval, max.scale = max.scale[i, j], use.na = use.na[i, j], expected = expected[i, j], 
					residuals = NULL, model.opt = model.opt[[i]][[j]], gap.prop = gap.prop[i, j], gap.mult = gap.mult[i, j], 
					col = col[i, j], col.opt = col.opt[[i]][[j]], label = label[i, j], label.opt = label.opt[[i]][[j]], vp = vp1)
				
			}else{
#					if(mosaic[i, j]){
#						f <- paste("Freq",f,sep="")
#						
#						mosaic(as.formula(f),data=data,labeling=label[i, j])
#					}else{
					rmb.formula(as.formula(f), data=data, cut = cut[c(i,j,tv)], col.vars = NULL, spine = spine[i, j], circular = circular[i, j], 
					eqwidth = eqwidth[i, j], cat.ord = NULL, freq.trans = freq.trans[[i]][[j]], 
					num.mode = num.mode[i, j], innerval = innerval, max.scale = max.scale[i, j], use.na = use.na[i, j], expected = expected[i, j], 
					residuals = NULL, model.opt = model.opt[[i]][[j]], gap.prop = gap.prop[i, j], gap.mult = gap.mult[i, j], 
					col = col[i, j], col.opt = col.opt[[i]][[j]], label = label[i, j], label.opt = label.opt[[i]][[j]], vp = vp1)
					#}
			}
			popViewport()
		}else{
			pushViewport(vp.mat)
			vp0 <- viewport(layout.pos.row = i,layout.pos.col=j)
			vp1 <- viewport(x=0.5,y=0.5,width=0.96,height=0.96)
			pushViewport(vp0)
			
			grid.rect()
			#pushViewport(vp1)
			if( i != tv | !plot.tv){
				f <- paste("~",paste(fterms[j],".1",sep=""),"+",fterms[i],"+",names(x)[tv],sep="")
				dset <- x[, c(rep(match(fterms[i],names(x)),2),tv)]
				if(fi < nv+3){
					dset$Freq <- x[,fi]	
				}
			}else{
				f <- paste("~",paste(fterms[j],".1",sep=""),"+",paste(fterms[j],".2",sep=""),"+",names(x)[tv],sep="")
				dset <- x[, c(rep(match(fterms[i],names(x)),2),tv)]	
				if(fi < nv+3){
					dset$Freq <- x[,fi]	
				}
			}
			
			rmb.formula(as.formula(f), data=dset, cut = cut[c(i,j,tv)], col.vars = NULL, spine = spine[i, j], circular = circular[i, j], 
    eqwidth = eqwidth[i, j], cat.ord = NULL, freq.trans = freq.trans[[i]][[j]], 
    num.mode = num.mode[i, j], innerval = innerval, max.scale = max.scale[i, j], use.na = use.na[i, j], expected = expected[i, j], 
    residuals = NULL, model.opt = model.opt[[i]][[j]], gap.prop = gap.prop[i, j], gap.mult = gap.mult[i, j], 
    col = col[i, j], col.opt = col.opt[[i]][[j]], label = label[i, j], label.opt = label.opt[[i]][[j]], vp = vp1)
			grid.text(fterms[i],x=0.5,y=0.5,vp=vp0,gp=gpar(fontsize=10, col=hsv(0,0,0,alpha=0.4)))
			
			popViewport()
		}
	})
})
    return(env)
}






## ---------------------------------------------------------------------------------------------------------------- #

applyspecs = function(opt,env,ids){

dims <- dim(env$spine)
if(is.na(ids[1])){
	ids1 <- 1:dims[1]	
}else{
	ids1 <- ids[1]	
}
if(is.na(ids[2])){
	ids2 <- 1:dims[2]	
}else{
	ids2 <- ids[2]	
}
		
#if( "spine" %in% names(opt) ) env$spine[ids[1],ids[2]] <- opt$spine
#if( "circular" %in% names(opt) ) env$circular[ids[1],ids[2]] <- opt$circular
if( "mode" %in% names(opt) ){
	if(opt$mode %in% c("circular", "circ", "c", "circles", "pie", "p", "piecharts")){
		env$circular[ids1,ids2] <- TRUE
		env$spine[ids1,ids2] <- FALSE
	}
	if(opt$mode %in% c("spine", "s")){
		env$circular[ids1,ids2] <- FALSE
		env$spine[ids1,ids2] <- TRUE	
	}
	if(opt$mode %in% c("bars", "bar", "b")){
		env$spine[ids1,ids2] <- FALSE
		env$spine[ids1,ids2] <- FALSE	
	}
}
if( "eqwidth" %in% names(opt) ) env$eqwidth[ids1,ids2] <- opt$eqwidth
if( "col" %in% names(opt) ) env$col[ids1,ids2] <- opt$col
if( "num.mode" %in% names(opt) ) env$num.mod[ids1,ids2] <- opt$num.mode
if( "max.scale" %in% names(opt) ) env$max.scale[ids1,ids2] <- opt$max.scale
if( "use.na" %in% names(opt) ) env$use.na[ids1,ids2] <- opt$use.na

if( "expected" %in% names(opt) ){
	if(!is.logical(opt$expected)){
		opt$expected <- TRUE
	}	
	env$expected[ids1,ids2] <- opt$expected
}

if( "gap.prop" %in% names(opt) ) env$gap.prop[ids1,ids2] <- opt$gap.prop
if( "gap.mult" %in% names(opt) ) env$gap.mult[ids1,ids2] <- opt$gap.mult
if( "label" %in% names(opt) ) env$label[ids1,ids2] <- opt$label

if( "freq.trans" %in% names(opt) ){
	for(i in ids1){
		for(j in ids2){
			env$freq.trans[[i]][[j]] <- opt$label.opt	
		}
	}		
}
if( "col.opt" %in% names(opt) ){
	for(i in ids1){
		for(j in ids2){
			env$col.opt[[i]][[j]] <- opt$col.opt	
		}
	}
}
if( "label.opt" %in% names(opt) ){
	for(i in ids1){
		for(j in ids2){
			env$label.opt[[i]][[j]] <- opt$label.opt	
		}
	}
}
if( "model.opt" %in% names(opt) ){
	for(i in ids1){
		for(j in ids2){
			env$model.opt[[i]][[j]] <- opt$model.opt	
		}
	}
}
#eof
}




## ---------------------------------------------------------------------------------------------------------------- #




innerval <- function(x, p = 0.95, data.points = TRUE){
	if(p == 1){
		return( range(x) )	
	}
	n <- length(x)
	ni <- ceiling(p*n)
	
	x <- na.omit(x)
	#mm <- mean(data, na.rm = T)
	mm <- median(x)
	a <- (1-p)/2#min(p/2, (1-p)/2)
	min.x <- min( abs(quantile(x,c(a,1-a))-mm))
	max.x <- max( abs(quantile(x,c(a,1-a))-mm))
if(min.x == max.x){
	return( mm + c(-min.x,min.x) )
}
	op <- optimize(opt.innerval, data=x, n = ni, center = mm, interval = c(min.x, max.x))#diff(range(data,na.rm=TRUE))
bd <- mm + c(-op$minimum,op$minimum)
if(data.points){
	bd[1] <- min(x[which(x>bd[1])], na.rm = TRUE)
	bd[2] <- max(x[which(x<bd[2])], na.rm = TRUE)	
}
return( bd )
}

## ---------------------------------------------------------------------------------------------------------------- #

opt.innerval <- function(x, data, n, center){
inner <-sum( data <= center+x) - sum( data < center-x)
return( (inner-n)^2)
}




