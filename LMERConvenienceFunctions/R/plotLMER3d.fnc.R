plotLMER3d.fnc<-function(model=NULL,
			pred,
			intr,
			plot.type="contour", # or "persp" or "persp3d" or "image.plot"
			xlim = range(x, na.rm = TRUE),
			ylim = range(y, na.rm = TRUE),
			zlim=range(z, na.rm = TRUE), 
			xlab=NULL, 
			ylab=NULL, 
			zlab=NULL, 
			main=NULL, 
                        shift=0,
                        scale=1,
			cex=1,
			fun=NA, 
			n=30,
			color="topo",
			alpha=1,
			alpha.rs=0.65,
			alpha.u=1,
			lit=TRUE,
			theta=0,
			phi=0,
			contourstepsize=0.2,
            legend.args=NULL,
			play3d=FALSE,   # or TRUE or a list with first element axis, e.g., c(0,0,1), second element rpm, 
					# e.g., 4, and third element duration, e.g., 20.			
			ref.surf=FALSE,
			underneath=FALSE,
			add.raw=FALSE,
			color.raw="grey",
			alpha.raw=0.5,
			rug=FALSE,
			rug.u=FALSE,
			plot.dat="default",
			path="default",
                        ...

){
	if(is.null(model) && plot.dat=="default"){
		stop("either provide a value to argument ''model'' (an object of class ''mer'') \n 
			or provide a file name containing plotting information \n")
	}
####################################################
#               DEFINE SOME FUNCTIONS              #
####################################################

####################################################
# define function getPos.fnc
getPos.fnc<-function (vec, pos) 
{
    if (pos == "end") 
        return(length(vec))
    else {
        if (pos == "beg") 
            return(1)
        else return(floor(length(vec)/2))
    }
}


####################################################
# define function getRange.fnc
getRange.fnc<-function (lst) 
{
    v = vector()
    for (i in 1:length(lst)) {
        if (is.data.frame(lst[[i]])) {
            if ("lower" %in% colnames(lst[[i]])) {
                v = c(v, as.vector(lst[[i]][, c("Y", "lower", 
                  "upper")]))
            }
            else {
                v = c(v, as.vector(lst[[i]][, "Y"]))
            }
        }
        else {
            for (j in 1:length(lst[[i]])) {
                if ("lower" %in% colnames(lst[[i]][[j]])) {
                  v = c(v, as.vector(lst[[i]][[j]][, c("Y", "lower", 
                    "upper")]))
                }
                else {
                  v = c(v, as.vector(lst[[i]][[j]][, "Y"]))
                }
            }
        }
    }
    return(range(v))
}


####################################################
# define function makeDefaultMatrix.fnc
makeDefaultMatrix.fnc<-function (model, n = 100, conditioningPred = "", conditioningValue = NULL, control = NA) {
    coefs = fixef(model)
    ncoefs = length(coefs)
    X = getME(model,"X")
    if (!is.null(names(fixef(model)))) {
        colnames(X) = names(fixef(model))
    }
    nams = strsplit(names(coefs), ":")
    if (is.character(conditioningValue)) {
        condName = paste(conditioningPred, conditioningValue, sep = "")
    }
    else {
        condName = conditioningPred
    }
    m = matrix(0, n, ncoefs)
    rownames(X) = 1:nrow(X)
    for (i in 1:ncoefs) {
        if (names(coefs[i]) == "(Intercept)") {
            m[, i] = rep(1, n)
        }
        else {
            v = names(table(X[, names(coefs[i])]))
            if (length(v) == 2 & v[1] == "0" & v[2] == "1") {
                if (condName == names(coefs)[i]) 
                  m[, i] = rep(1, length = n)
                else m[, i] = rep(0, length = n)
            }
            else {
                if (length(nams[[i]]) == 1) {
                  if (condName == names(coefs)[i]) {
                    m[, i] = rep(conditioningValue, length = n)
                  }
                  else {
                    if (regexpr("^poly\\(", names(coefs[i])) > 
                      0) {
                      if (regexpr("1$", names(coefs[i])) > 0) {
                        maxval = max(X[X[, i] < median(X[, i]), 
                          ][, i])
                        maxbelowmedianpos = which(X[, i] == maxval)[1]
                      }
                      m[, i] = rep(X[maxbelowmedianpos, names(coefs[i])], 
                        length = n)
                    }
                    else {
                      if (regexpr("^rcs\\(", names(coefs[i])) > 
                        0) {
                        if (regexpr("[^']$", names(coefs[i])) > 
                          0) {
                          maxval = max(X[X[, i] < median(X[, 
                            i]), ][, i])
                          maxbelowmedianpos = which(X[, i] == 
                            maxval)[1]
                        }
                        m[, i] = rep(X[maxbelowmedianpos, names(coefs[i])], 
                          length = n)
                      }
                      else {
                        m[, i] = rep(median(X[, names(coefs[i])]), 
                          length = n)
                      }
                    }
                  }
                }
                else {
                  m[, i] = rep(0, length = n)
                }
            }
        }
    }
    colnames(m) = colnames(X)
    if (!is.na(control)[[1]]) {
        controlPredName = as.vector(control[[1]])
        if (!is.element(controlPredName, colnames(m))) {
            stop(paste("the control predictor name", controlPredName, "is not a valid column name\n", sep = " "))
        }
        else {
            m[, controlPredName] = rep(as.vector(control[[2]]), nrow(m))
        }
    }
    return(m)
}

####################################################
# define function degreesOrKnots.fnc
degreesOrKnots.fnc<-function (name) 
{
    s = strsplit(name, " ")[[1]][2]
    s2 = strsplit(s, "[,)]")
    return(as.numeric(s2[[1]][1]))
}


####################################################
# define function getKnots.fnc
getKnots.fnc<-function (colnms, xlb) 
{
    pos = grep(paste("rcs\\(", xlb, ",", sep = ""), colnms)
    tmp = strsplit(colnms[pos[1]], ")")[[1]][1]
    return(as.numeric(strsplit(tmp, " ")[[1]][2]))
}

####################################################
# define function implementInteractions.fnc
implementInteractions.fnc<-function (m) 
{
    nams = strsplit(colnames(m), ":")
    for (i in 1:length(nams)) {
        if (length(nams[[i]]) > 1) {
            m[, i] = m[, nams[[i]][1]]
            for (j in 2:length(nams[[i]])) {
                m[, i] = m[, i] * m[, nams[[i]][j]]
            }
        }
    }
    return(m)
}


####################################################
# define function transforming.fnc
transforming.fnc<-function (y, fun) 
{
    if (is.function(fun)) {
        return(fun(y))
    }
    else return(y)
}

####################################################
# define function parsePredName.fnc
parsePredName.fnc<-function (name) 
{
    s = strsplit(name, "[\\(\\)]")[[1]][2]
    s2 = strsplit(s, ", ")[[1]]
    return(list(baseName = s2[1], knots = as.numeric(s2[2])))
}


####################################################
# define function makeDefaultMatrix.fnc
preparePredictor.fnc<-function (pred, model, m, ylabel, fun, val, xlabel, ranefs, ...) {
    X = getME(model,"X")
    if (!is.null(names(fixef(model)))) {
        colnames(X) = names(fixef(model))
    }
    polynomial = FALSE
    namesplit = strsplit(pred, ", ")[[1]]
    a = regexpr("poly\\(", namesplit[1])
    if ((a == 1) & (attr(a, "match.length") == 5)) {
        polynomial = TRUE
        degree = degreesOrKnots.fnc(pred)
    }
    rcspline = FALSE
    namesplit = strsplit(pred, ", ")[[1]]
    a = regexpr("rcs\\(", namesplit[1])
    if ((a == 1) & (attr(a, "match.length") == 4)) {
        rcspline = TRUE
        knots = degreesOrKnots.fnc(pred)
    }
    if ((!polynomial) & (!rcspline)) {
        pred2 = paste("rcs\\(", pred, sep = "")
        if (length(grep(pred2, colnames(X))) > 0) {
            rcspline = TRUE
        }
        else {
            pred2 = paste("poly\\(", pred, sep = "")
            if (length(grep(pred2, colnames(X))) > 0) {
                polynomial = TRUE
            }
        }
    }
    isfactor = FALSE
    fixefs = fixef(model)
    if (!is.na(ranefs[[1]])) {
        nm = as.vector(ranefs[[4]])
        if (nm %in% names(fixefs)) {
            blup = ranef(model)[[ranefs[[1]]]][ranefs[[2]], ranefs[[3]]]
            fixefs[nm] = fixefs[nm] + blup
            fixefs = as.numeric(fixefs)
        }
        else stop(paste(nm, "is not a valid predictor name, check 'fixef(model)'\n", 
            sep = " "))
    }
    if ((pred %in% colnames(model@frame)) & polynomial == FALSE & 
        rcspline == FALSE) {
        if (is.numeric(model@frame[, pred])) {
            if (pred %in% colnames(X)) {
                m[, pred] = seq(min(X[, pred]), max(X[, pred]), 
                  length = nrow(m))
                m = implementInteractions.fnc(m)
                vals = m %*% fixefs
                vals = transforming.fnc(vals, fun)
                dfr = data.frame(X = m[, pred], Y = vals)
                dfr$Predictor = rep(xlabel, nrow(dfr))
                dfr$Type = rep(isfactor, nrow(dfr))
                if (is.na(val)) {
                  dfr$Interaction = rep(NA, nrow(dfr))
                }
                else {
                  dfr$Interaction = rep(val, nrow(dfr))
                }
            }
            else {
                stop(paste(pred, " is not plotted (not a fixed effect predictor)\n"))
            }
        }
        else {
            if (is.logical(model@frame[, pred])) 
                model@frame[, pred] = factor(model@frame[, pred])
            if (is.factor(model@frame[, pred])) {
                isfactor = TRUE
                factnames = paste(pred, levels(model@frame[, pred])[-1], sep = "")
                m = m[1:(length(factnames) + 1), ]
                for (i in 1:length(factnames)) {
                  m[i + 1, factnames[i]] = 1
                }
                m = implementInteractions.fnc(m)
                vals = m %*% fixefs
                vals = transforming.fnc(vals, fun)
                x = 1:nrow(m)
                dfr = data.frame(X = x, Y = vals)
                dfr$Predictor = rep(xlabel, nrow(dfr))
                dfr$Type = rep(isfactor, nrow(dfr))
                if (is.na(val)) {
                  dfr$Interaction = rep(FALSE, nrow(dfr))
                }
                else {
                  dfr$Interaction = rep(TRUE, nrow(dfr))
                }
                dfr$Levels = levels(model@frame[, pred])
            }
            else {
                cat("warning: I don't know how to handle ", pred, 
                  "\n")
            }
        }
    }
    else {
        if (!(pred %in% colnames(X))) {
            pos = grep(pred, colnames(X), fixed = TRUE)
            degree = 1
            knots = 1
            if (length(pos) > 0) {
                name = colnames(X)[pos][1]
                namesplit = strsplit(name, ", ")[[1]]
                a = regexpr("poly", namesplit[1])
                if ((a == 1) & (attr(a, "match.length") == 4)) {
                  polynomial = TRUE
                  degree = as.numeric(namesplit[2])
                  xlabel = parsePredName.fnc(pred)[[1]]
                  name = pred
                }
                if (!polynomial) {
                  a = regexpr("rcs", namesplit[1])
                  if ((a == 1) & (attr(a, "match.length") == 
                    3)) {
                    rcspline = TRUE
                    aa = parsePredName.fnc(name)
                    knots = aa[[2]]
                    xlabel = aa[[1]]
                  }
                }
            }
        }
        else {
            namesplit = strsplit(pred, ", ")[[1]]
            name = pred
            arg2 = as.numeric(substr(namesplit[2], 1, nchar(namesplit[2]) - 1))
            cat("DIT ZOU DOOD STUK CODE MOETEN ZIJN\n")
        }
        if (polynomial | rcspline) {
            if (is.na(xlabel)) {
                xlabel = pred
            }
            if (polynomial) {
                hasPoly = FALSE
                if (length(grep("^poly\\(", pred)) > 0) {
                  vec = paste(name, "1", sep = "")
                  hasPoly = TRUE
                }
                else {
                  xlabel = pred
                  vec = paste("poly(", name, ", ", degree, ", raw = TRUE)1", sep = "")
                }
                name1 = vec
                m[, vec] = seq(min(X[, vec]), max(X[, vec]), 
                  length = nrow(m))
                for (i in 2:degree) {
                  if (hasPoly) {
                    vec = c(vec, paste(name, as.character(i), 
                      sep = ""))
                  }
                  else {
                    vec = c(vec, paste("poly(", name, ", ", degree, ", raw = TRUE)", as.character(i), sep = ""))
                  }
                  m[, vec[i]] = m[, vec[i - 1]] * m[, vec[1]]
                }
            }
            else {
                if (length(grep("^rcs\\(", pred)) > 0) {
                  nms = unlist(parsePredName.fnc(pred))
                  basename = nms[1]
                  knots = as.numeric(nms[2])
                  name1 = paste("rcs(", basename, ", ", knots, 
                    ")", basename, sep = "")
                  xlabel = basename
                }
                else {
                  knots = getKnots.fnc(colnames(X), pred)
                  name1 = paste("rcs(", pred, ", ", knots, ")", 
                    pred, sep = "")
                }
                vec = rep(name1, knots - 1)
                vec[2] = paste(vec[1], "'", sep = "")
                if (knots > 3) {
                  for (i in 3:(knots - 1)) {
                    vec[i] = paste(vec[i - 1], "'", sep = "")
                  }
                }
                mtmp = unique(X[, vec])
                if (nrow(mtmp) <= nrow(m)) {
                  m = m[1:nrow(mtmp), ]
                  m[, vec] = mtmp
                }
                else {
                  vecIndices = c(1, sort(sample(2:(nrow(mtmp) - 
                    1), nrow(m) - 2)), nrow(mtmp))
                  m[, vec] = mtmp[vecIndices, ]
                }
                m = m[order(m[, vec[1]]), ]
            }
            m = implementInteractions.fnc(m)
            vals = m %*% fixefs
            vals = transforming.fnc(vals, fun)
            dfr = data.frame(X = m[, vec[1]], Y = vals)
            dfr$Predictor = rep(xlabel, nrow(dfr))
            dfr$Type = rep(isfactor, nrow(dfr))
            if (is.na(val)) {
                dfr$Interaction = rep(FALSE, nrow(dfr))
            }
            else {
                dfr$Interaction = rep(val, nrow(dfr))
            }
        }
        else {
            stop(paste("unknown function used in", pred, "\n"))
        }
    }
    return(dfr)
}

#######################################################################	
# define function plotLMERTweaked
plotLMERTweaked<-function (model, xlabel = NA, xlabs = NA, ylabel = NA, ylimit = NA, ilabel = NA, fun = NA, pred = NA, control = NA, ranefs = NA, n = 100, intr = NA, lockYlim = TRUE, addlines = FALSE, withList = FALSE, cexsize = 0.5, linecolor = 1, addToExistingPlot = FALSE, verbose = TRUE, ...) {
    if (!(is(model, "lmerMod") || (is(model, "glmerMod")))) {
        stop("argument should be a (g)lmerMod object")
    }
    if (!is.na(xlabel[1])) {
        if (!is.character(xlabel)) 
            stop("xlabel should be a string\n")
    }
    if (!is.na(ylabel)) {
        if (!is.character(ylabel)) 
            stop("ylabel should be a string\n")
    }
    if (!is.na(ylimit[1])) {
        if ((!is.numeric(ylimit)) | (length(ylimit) != 2)) 
            stop("ylimit should be a two-element numeric vector\n")
    }
    if (!is.na(intr[1])) {
        if (!is.list(intr)) 
            stop("intr should be a list\n")
    }
    if (!is.numeric(n)) {
        stop("n should be an integer\n")
    }
    if (!is.na(pred)) {
        if (!is.character(pred)) 
            stop("pred should be a string\n")
    }
    if (!is.function(fun)) {
        if (!is.na(fun)) {
            stop("fun should be a function (not the name of a function)\n")
        }
    }
    if ((length(grep("^glmer", as.character(model@call))) == 
        1) & (length(grep("binomial", as.character(model@call))) == 
        1)) {
        if (!is.function(fun)) {
            fun = plogis
            if (verbose == TRUE) 
                cat("log odds are back-transformed to probabilities\n")
        }
    }
    if (is.na(pred)) 
        addToExistingPlot = FALSE
    conditioningPred = ""
    conditioningVals = NULL
    conditioningPos = NA
    conditioningColors = 1
    conditioningLines = 1
    if (!is.na(intr[[1]])) {
        conditioningPred = intr[[1]]
        conditioningVals = intr[[2]]
        conditioningPos = intr[[3]]
        if (length(intr) == 4) {
            conditioningColors = intr[[4]][[1]]
            if (length(conditioningColors) != length(conditioningVals)) {
                stop("number of colors and number of conditioning values mismatch")
            }
            conditioningLines = intr[[4]][[2]]
            if (length(conditioningLines) != length(conditioningLines)) {
                stop("number of line types and number of conditioning values mismatch")
            }
        }
    }
    if (length(ylimit) > 1) {
        lockYlim = FALSE
    }
    if (!is.na(control[[1]])) {
        if (!((length(control) == 2) & is.list(control))) {
            stop("control should be a two-element list\n")
        }
    }
    if (is.na(ylabel)) 
        ylabel = as.character(eval(model@call[2]$formula))[2]
    if (is.na(pred)) {
        predictors = colnames(model@frame)
        ranefnames = unique(names(ranef(model)))
        depvar = as.character(eval(model@call[2]$formula))[2]
        predictors = predictors[1:(which(predictors == ranefnames[1]) - 
            1)]
        predictors = predictors[!predictors %in% c(ranefnames, 
            depvar)]
    }
    else {
        predictors = pred
    }
    if (!is.na(xlabs[1])) {
        if (length(xlabs) != length(predictors)) {
            stop("number of labels in xlabs is not the same as the number of predictors\n")
        }
    }
    plots = list()
    for (i in 1:length(predictors)) {
        if (length(predictors) == 1) 
            xlabelShow = xlabel
        else {
            if (!is.na(xlabs[1])) {
                xlabelShow = xlabs
            }
            else {
                xlabelShow = NA
            }
        }
        if (is.na(xlabel[1]) | length(predictors) > 1) {
            xlabel = predictors[i]
        }
        if ((length(predictors) == 1) & (!is.null(conditioningVals))) {
            if (is.null(conditioningColors)) {
                colors = rep(1, length(conditioningVals))
                lineTypes = rep(1, length(conditioningVals))
            }
            else {
                colors = conditioningColors
                lineTypes = conditioningLines
                if (length(colors) < length(conditioningVals)) {
                  nc = (length(conditioningVals)%%length(colors)) + 
                    1
                  colors = rep(colors, nc)
                }
                if (length(lineTypes) < length(conditioningVals)) {
                  nc = (length(conditioningLines)%%length(lineTypes)) + 
                    1
                  lineTypes = rep(lineTypes, nc)
                }
            }
            val = conditioningVals[1]
            m = makeDefaultMatrix.fnc(model, n, conditioningPred, 
                val, control)
            subplots = list()
            dfr = preparePredictor.fnc(predictors[i], model, m, ylabel, fun, val, xlabel = xlabel, ranefs, lty = 1, col = 0, ...)
            subplots[[1]] = dfr
            if (verbose == TRUE) {
                cat("effect sizes (ranges) for the interaction of ", 
                  predictors[i], " and ", conditioningPred, ":\n")
                cat("   ", conditioningPred, " = ", val, ": ", 
                  max(dfr$Y) - min(dfr$Y), "\n")
            }
            for (j in 2:length(conditioningVals)) {
                val = conditioningVals[j]
                m = makeDefaultMatrix.fnc(model, n, conditioningPred, 
                  val, control)
                dfr = preparePredictor.fnc(predictors[i], model, 
                  m, ylabel, fun, val, ranefs, lty = j, xlabel = xlabel, ...)
                subplots[[j]] = dfr
                if (verbose == TRUE) {
                  cat("   ", conditioningPred, " = ", val, ": ", 
                    max(dfr$Y) - min(dfr$Y), "\n")
                }
            }
            plots[[i]] = subplots
        }
        else {
            lineTypes = 1
            m = makeDefaultMatrix.fnc(model, n, "", NULL, control)
            dfr = preparePredictor.fnc(predictors[i], model, 
                m, ylabel, fun, val = NA, xlabel = xlabel, ranefs, ...)
            plots[[i]] = dfr
            if (verbose == TRUE) {
                cat("effect size (range) for ", predictors[i], 
                  "is ", max(dfr$Y) - min(dfr$Y), "\n")
            }
        }
    }
    names(plots) = predictors
    if (!is.na(ilabel)) {
        intrName = ilabel
    }
    else {
        intrName = conditioningPred
    }
    #plotAll.fnc(plots, sameYrange = lockYlim, ylabel, xlabel = xlabelShow, 
        #intrName = intrName, pos = conditioningPos, ylimit = ylimit, 
        #addlines = addlines, cexsize = cexsize, conditioningVals = conditioningVals, 
        #conditioningColors = colors, conditioningLines = lineTypes, 
        #lineColor = linecolor, addToExistingPlot, ...)
    #if (withList) 
        return(invisible(plots))
}


###########################################################################
#                    rest of function continues here                      #
###########################################################################
	# set labels if NULL
	if(is.null(xlab)){
		xlab=pred
	}

	if(is.null(ylab)){
		ylab=intr
	}

     	options(warn=-1)
	if(is.null(zlab)){
		if(try(is.null(model),silent=TRUE)){
			zlab="Response"
		}else{
			zlab<-as.character(model@call)[2]
			zlab<-gsub(" ","",unlist(strsplit(zlab,"~"))[1])
		}
	}
     	options(warn=0)

	if(is.null(main)){
		if(plot.type=="contour"){
			main=zlab
		}else{
			main=""
		}
	}

	# create file name for saving in temp dir
	if(plot.dat[1]!=FALSE){
		if(path=="default"){
			temp.dir<-tempdir()
		}else{
			temp.dir=path
		}
		if(plot.dat[1]=="default"){
			model.name<-as.character(model@call)
			model.name<-gsub(" ","",model.name)
			model.name<-paste(model.name[1],"___",model.name[2],"___",model.name[3],"___",pred,"_",intr,sep="")
			model.name<-gsub("\\+","__",model.name)
			model.name<-gsub("\\:","_",model.name)
			model.name<-gsub("\\*","_",model.name)
			model.name<-gsub("\\^","_",model.name)
			model.name<-gsub("\\|","_",model.name)
			model.name<-gsub("\\~","_",model.name)
			model.name<-gsub("\\(","WWW",model.name)
			model.name<-gsub("\\)","WWW",model.name)
			model.name<-paste(model.name,".rda",sep="")
		}else{
			model.name=paste("lmer___",plot.dat,".rda",sep="")
		}
	}


	# get previously generated plotting info if it exists
	if(plot.dat[1]!=FALSE){
		if(!model.name%in%list.files(path=temp.dir,pattern="lmer___.*\\.rda$")){
			# create LMER plot and store values
			list1<-plotLMERTweaked(model=model,fun=fun,pred=pred,intr=list(intr,
				quantile(model@frame[,intr],seq(0,1,1/n)),"end",list(seq(1,length(seq(0,1,1/n)),1),
				seq(1,length(seq(0,1,1/n)),1))),n=n,verbose=TRUE)
		
			# get info from stored plotting list
			length.intr<-eval(parse(text=paste("length(list1$",pred,")",sep="")))
			x<-eval(parse(text=paste("list1$",pred,"[[1]]$X",sep="")))
			y=quantile(model@frame[,intr],seq(0,1,1/n))
		
			# create plotting matrix
			for(i in 1:length.intr){
				if(i==1){
					z<-eval(parse(text=paste("list1$",pred,"[[",i,"]]$Y",sep="")))
				}else{
					z<-cbind(z,eval(parse(text=paste("list1$",pred,"[[",i,"]]$Y",sep=""))))
				}	
			}
			z<-as.matrix(z)
		
			# add row and column names to matrix
			rownames(z)<-x
			colnames(z)<-y
		
			# remove identical columns
			rem<-vector("numeric")
			for(i in 2:ncol(z)){
				if(unique(z[,i-1]==z[,i])[1]){
					rem<-c(rem,i)
				}
			}
			if(length(rem)!=0){
				z<-z[,-rem]
			}
			
			save(z,file=file.path(temp.dir,model.name))
		}else{
			load(file.path(temp.dir,model.name))
		}
	}


	if(is.null(zlim)){
		zlim=range(z,na.rm=TRUE)
	}


	if(plot.type[1]=="image.plot"){
      contourlevels = seq(zlim[1], zlim[2], by=contourstepsize)
	                                        
	  # Determine color.
          if(color=="heat"){
            pal=heat.colors(50)
            con.col=3
          }else if(color=="topo"){
            pal=topo.colors(50)
            con.col=2
          }else if(color=="cm"){
            pal=cm.colors(50)
            con.col=1
          }else if(color=="terrain"){
            pal=terrain.colors(50)
            con.col=2
          }else if(color=="gray"||color=="bw"||color=="grey"){
            pal=gray(seq(0.1,0.9,length=50))
            con.col=1
          }else{
	    stop("color scheme not recognised")
	  }

          x<-as.numeric(rownames(z))
          y<-as.numeric(colnames(z))

          jpeg(filename=file.path(tempdir(),"tmp.jpeg"))
	  err<-try(image.plot(x,y,z,col=pal,main=main,legend.args=legend.args,
                xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,zlim=zlim),
                silent=TRUE)
	  dev.off()

	  if(length(grep("Error",err))>0){
            if(length(unique(x))!=length(x)){
              x<-sort(jitter(x,factor=0.01))
            }
            if(length(unique(y))!=length(y)){
              y<-sort(jitter(y,factor=0.01))
            }

            jpeg(filename=file.path(tempdir(),"tmp.jpeg"))
	    err<-try(image.plot(x,y,z,col=pal,main=main,legend.args=legend.args,
                  xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,zlim=zlim),
                  silent=TRUE)
	    dev.off()
          }

          if(length(grep("Error",err))>0){
	    cat("\tplotting anyways, but will not use supplied x- and y-values ...\n")
            image.plot(z,col=pal,main=main,legend.args=legend.args,
              xlab=paste(xlab,"-- Random Units",sep=" "),ylab=paste(ylab,"-- Random Units",
              sep=" "),xlim=xlim,ylim=ylim,zlim=zlim,...)
	  }else{
            image.plot(x,y,z,col=pal,main=main,legend.args=legend.args,
              xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,zlim=zlim,...)
            contour(x,y,z,add=TRUE,nlevel=round(contourlevels,2),col=con.col,...)
          }

	  if(rug){
            xy<-expand.grid(as.numeric(rownames(z)),as.numeric(colnames(z)))
	    points(xy[,1],xy[,2],pch=19,cex=0.05)
          }

          return(invisible(list(z=z,col=pal)))
        }else if(plot.type[1]=="contour"){
		contourlevels = seq(zlim[1], zlim[2], by=contourstepsize)
			
		# Determine color.
        	if(color=="heat"){
            		pal=heat.colors(50)
            		con.col=3
        	}else if(color=="topo"){
            		pal=topo.colors(50)
            		con.col=2
        	}else if(color=="cm"){
            		pal=cm.colors(50)
            		con.col=1
        	}else if(color=="terrain"){
            		pal=terrain.colors(50)
            		con.col=2
        	}else if(color=="gray"||color=="bw"||color=="grey"){
            		pal=gray(seq(0.1,0.9,length=50))
            		con.col=1
        	}else{
			stop("color scheme not recognised")
		}

                x<-as.numeric(rownames(z))
                y<-as.numeric(colnames(z))

		jpeg(filename=file.path(tempdir(),"tmp.jpeg"))
		err<-try(image(x=x,y=x,z=z,col=pal,zlim=zlim,main=main,
                        cex.main=cex,cex.lab=cex,cex.axis=cex,xlab=xlab,
                        ylab=ylab,axes=TRUE),silent=TRUE)
		dev.off()

		if(length(grep("Error",err))>0){
                        if(length(unique(x))!=length(x)){
                            x<-sort(jitter(x,factor=0.01))
                        }
                        if(length(unique(y))!=length(y)){
                            y<-sort(jitter(y,factor=0.01))
                        }

		        jpeg(filename=file.path(tempdir(),"tmp.jpeg"))
		        err<-try(image(x=x,y=y,z=z,col=pal,zlim=zlim,main=main,
			        cex.main=cex,cex.lab=cex,cex.axis=cex,xlab=xlab,
                                ylab=ylab,axes=TRUE),silent=TRUE)
		        dev.off()
                 }

		if(length(grep("Error",err))>0){
			cat("\tplotting anyways, but will not use supplied x- and y-values ...\n")
			image(z=z,col=pal,zlim=zlim,main=main,cex.main=cex,cex.lab=cex,cex.axis=cex,
				axes=TRUE,xlab=paste(xlab,"-- Random Units",sep=" "),
				ylab=paste(ylab,"-- Random Units",sep=" "),...)
			contour(z=z,zlim=zlim,add=TRUE,levels=round(contourlevels,2),axes=FALSE,...)
		}else{
		        image(x=x,y=y,z=z,col=pal,zlim=zlim,main=main,cex.main=cex,cex.lab=cex,
                          cex.axis=cex,xlab=xlab,ylab=ylab,axes=TRUE,...)
			contour(x=x,y=y,z=z,zlim=zlim,add=TRUE,levels=round(contourlevels,2),
                            axes=FALSE,...)
                }

		rm(err)

		if(rug){
			xy<-expand.grid(as.numeric(rownames(z)),as.numeric(colnames(z)))
			points(xy[,1],xy[,2],pch=19,cex=0.05)
		}

		box()

		return(invisible(list(z=z,col=pal)))

	}else if (plot.type[1]=="persp"){
		# the color portion of this code is adapted from the persp() help page
		#par(bg="white")
		nrz<-nrow(z)
		ncz<-ncol(z)

		# Create a function interpolating colors in the range of specified colors
        	if(color=="heat"){
            		jet.colors<-colorRampPalette(heat.colors(50))
        	}else if(color=="topo"){
			#jet.colors <- colorRampPalette( c("purple","blue", "green","yellow","red","white") ) 
			jet.colors <- colorRampPalette(topo.colors(50)) 
        	}else if(color=="cm"){
            		jet.colors<-colorRampPalette(cm.colors(50))
        	}else if(color=="terrain"){
            		jet.colors<-colorRampPalette(terrain.colors(50))
        	}else if(color=="gray"||color=="bw"||color=="grey"){
            		jet.colors<-colorRampPalette(gray(seq(0.1,0.9,length=7)))
        	}else{
			stop("color scheme not recognised")
		}

		# Generate the desired number of colors from this palette
		nbcol<-100
		color<-jet.colors(nbcol)

		# Compute the z-value at the facet centres
		zfacet<-z[-1,-1]+z[-1,-ncz]+z[-nrz,-1]+z[-nrz,-ncz]

		# Recode facet z-values into color indices
		facetcol<-cut(zfacet,nbcol)

                x<-as.numeric(rownames(z))
                y<-as.numeric(colnames(z))

		jpeg(filename=file.path(tempdir(),"tmp.jpeg"))
		err<-try(persp(x=x,y=y,z=z,ticktype="detailed",
                        col=color[facetcol],phi=phi,theta=theta,
			zlab=zlab,zlim=zlim,xlab=xlab,ylab=ylab,
                        main=main,axes=TRUE),silent=TRUE)
		dev.off()


		if(length(grep("Error",err))>0){
                        if(length(unique(x))!=length(x)){
                            x<-sort(jitter(x,factor=0.01))
                        }
                        if(length(unique(y))!=length(y)){
                            y<-sort(jitter(y,factor=0.01))
                        }

		        jpeg(filename=file.path(tempdir(),"tmp.jpeg"))
		        err<-try(persp(x=x,y=y,z=z,ticktype="detailed",
                                col=color[facetcol],phi=phi,theta=theta,
			        zlab=zlab,zlim=zlim,xlab=xlab,ylab=ylab,
                                main=main,axes=TRUE),silent=TRUE)
		        dev.off()
                 }

		if(length(grep("Error",err))>0){
			cat("\tplotting anyways, but will not use supplied x- and y-values ...\n")
			persp(z=z,ticktype="detailed",col=color[facetcol],phi=phi,theta=theta,
				zlab=zlab,zlim=zlim,xlab=paste(xlab,"-- Random Units",sep=" "),
				ylab=paste(ylab,"-- Random Units",sep=" "),main=main,axes=TRUE,
                                ...)->res
		}else{
		        persp(x=x,y=y,z=z,ticktype="detailed",col=color[facetcol],
                            phi=phi,theta=theta,zlab=zlab,zlim=zlim,xlab=xlab,
                            ylab=ylab,main=main,axes=TRUE,...)->res
                }

		rm(err)


		if(rug){
			xy<-expand.grid(as.numeric(rownames(z)),as.numeric(colnames(z)))
			temp<-vector("numeric")
			for(i in 1:nrow(xy)){
				temp<-c(temp,z[as.character(xy$Var1[i]),as.character(xy$Var2[i])])
			}
			points(trans3d(xy[,1],xy[,2],temp,pmat=res),pch=19,cex=0.5)
		}

		return(invisible(list(z=z,col=color[facetcol])))

	}else{
		# the color portion of this code is adapted from the persp() help page
		#par(bg="white")
		nrz<-nrow(z)
		ncz<-ncol(z)

		# Create a function interpolating colors in the range of specified colors
        	if(color=="heat"){
            		jet.colors<-colorRampPalette(heat.colors(100))
        	}else if(color=="topo"){
			jet.colors <- colorRampPalette(topo.colors(100)) 
        	}else if(color=="cm"){
            		jet.colors<-colorRampPalette(cm.colors(100))
        	}else if(color=="terrain"){
            		jet.colors<-colorRampPalette(terrain.colors(100))
        	}else if(color=="gray"||color=="bw"||color=="grey"){
            		jet.colors<-colorRampPalette(gray(seq(0.1,0.9,length=7)))
        	}else{
			stop("color scheme not recognised")
		}

		# Generate the desired number of colors from this palette
		nbcol<-100
		color<-jet.colors(nbcol)

		# Compute the z-value at the facet centres
		zfacet<-z[-1,-1]+z[-1,-ncz]+z[-nrz,-1]+z[-nrz,-ncz]

		# Recode facet z-values into color indices
		facetcol<-cut(zfacet,nbcol)
		facetcol=color[facetcol]

		# this portion is from the persp3d() help page
		nx=length(rownames(z))
		ny=length(colnames(z))
		col <- rbind(1, cbind(matrix(facetcol, nx-1, ny-1), 1))


                x<-as.numeric(rownames(z))
                y<-as.numeric(colnames(z))

                op3d<-par3d()$cex
                par3d(cex=cex)

		# create persp3d plot
		jpeg(filename=file.path(tempdir(),"tmp.jpeg"))
		err<-try(persp3d(x=x,y=y,z=z,col=col,zlim=zlim,
				zlab=zlab,main=main,alpha=alpha,
                                smooth=FALSE,lit=lit,xlab=xlab,
                                ylab=ylab),silent=TRUE)
		dev.off()

		if(length(grep("Error",err))>0){
                        if(length(unique(x))!=length(x)){
                            x<-sort(jitter(x,factor=0.01))
                        }
                        if(length(unique(y))!=length(y)){
                            y<-sort(jitter(y,factor=0.01))
                        }

		        jpeg(filename=file.path(tempdir(),"tmp.jpeg"))
		        err<-try(persp3d(x=x,y=y,z=z,col=col,zlim=zlim,
				        zlab=zlab,main=main,alpha=alpha,
                                        smooth=FALSE,lit=lit,xlab=xlab,
                                        ylab=ylab),silent=TRUE)
		        dev.off()
                 }

		if(length(grep("Error",err))>0){
			cat("\tplotting anyways, but will not use supplied x- and y-values ...\n")
			persp3d(z=z,col=col,zlim=zlim,zlab=zlab,main=main,
                            alpha=alpha,smooth=FALSE,lit=lit,xlab=paste(xlab,
                            "-- Random Units",sep=" "),ylab=paste(ylab,
                            "-- Random Units",sep=" "),...)
		}else{
		        persp3d(x=x,y=y,z=z,col=col,zlim=zlim,zlab=zlab,
                            main=main,alpha=alpha,smooth=FALSE,lit=lit,
                            xlab=xlab,ylab=ylab,...)
                }

		if(add.raw){
			response<-gsub(" ","",strsplit(as.character(model@call)[2],"~")[[1]][1])
			xy<-ifelse(length(grep("Error",err))>0,FALSE,TRUE)
			plotRaw3d.fnc(data=model@frame,response=response,pred=pred,intr=intr,xy=xy,
				color=color.raw,alpha=alpha.raw,plot.type="persp3d",xlab="",ylab="",
				zlab="",main="",add=TRUE,shift=shift,scale=scale)
		}

		rm(err)

		if(ref.surf){
			persp3d(x=x,y=y,matrix(mean(z),ncol=ncol(z),
				nrow=nrow(z),byrow=TRUE),col="grey",
                                alpha=alpha.rs,zlim=zlim,add=TRUE,
                                lit=lit)
		}

		if(underneath){
			persp3d(x=x,y=y,z=(matrix(min(z),nrow(z),ncol(z))+zlim[1]),
				col=col,alpha=alpha.u,add=TRUE,smooth=FALSE,lit=lit,
                                zlim=zlim)	
			if(rug.u){
				for(i in 1:nrow(z)){
					plot3d(x=x[i],y=y,z=(matrix(min(z),
                                                nrow(z),ncol(z))+zlim[1])[i,],
                                                add=TRUE,col="black",size=3)
				}
			}
		}

		if(rug){
			for(i in 1:nrow(z)){
				plot3d(x=x[i],y=y,z=z[i,],add=TRUE,
                                    col="black",size=3)
			}
		}


		if(play3d || is.list(play3d)){
			if(is.list(play3d)){
				play3d(spin3d(axis=play3d[[1]], rpm=play3d[[2]]), duration=play3d[[3]])
			}else{
				play3d(spin3d(axis=c(0,0,1),rpm=4),duration=20)
			}
		}

                par3d(cex=op3d)

		return(invisible(list(x=x,y=y,z=z,col=col)))
	}
}
