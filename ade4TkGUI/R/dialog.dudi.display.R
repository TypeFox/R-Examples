################################
# Function to display a dudi
################################
"dialog.dudi.display" <- function(showCom, histCom, dudiname)
{
	tf <- tktoplevel()
	tkwm.title(tf, dudiname)
	done <- tclVar(0)
	
#
# Local Tk variables
#
	nfvar <- tclVar()
	chvar <- tclVar()
	chvar2 <- tclVar()
	
	xaxvar <- tclVar(1)
	yaxvar <- tclVar(2)

#
# intermediate functions to draw graphics
#
  
  "tablevalue" <- function(val)
	{
    g <- do.call("table.value",list(dftab = parse(text=paste(dudiname, val, sep=""))[[1]]))
    print(g@Call)
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}
  
  
	"scatterfunc" <- function()
	{
	  g <- do.call("scatter", list(x = parse(text=dudiname)[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
    assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}


	"scorefunc" <- function()
	{
    g <- do.call("score", list(x = parse(text=dudiname)[[1]], xax = parse(text=tclvalue(xaxvar))[[1]]))
    assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"plotfunc" <- function()
	{
    g <- do.call("plot", list(x = parse(text=dudiname)[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
  }

	"tabvalue" <- function()
	{
    g <- do.call("s.value", list(dfxy = parse(text=paste(dudiname, "$li", sep=""))[[1]], z = parse(text=paste(dudiname, "$tab", sep=""))[[1]], 
        xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]],
        psub.text = eval(parse(text=paste("names(",dudiname, "$tab)", sep=""))), psub.cex = 1.5))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelli" <- function()
	{
    g <- do.call("s.label",list(dfxy = parse(text=paste(dudiname, "$li", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labell1" <- function()
	{
	  g <- do.call("s.arrow",list(dfxy = parse(text=paste(dudiname, "$l1", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labeldlsdpcoa" <- function()
	{
    g <- do.call("s.label",list(dfxy = parse(text=paste(dudiname, "$dls", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labellidpcoa" <- function()
	{
    g <- do.call("s.label",list(dfxy = parse(text=paste(dudiname, "$li", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelco" <- function()
	{
		if (dclass[1] == "pca") {
			ncp <- names(dcall)
			nf1 <- names(formals(dudi.pca))
			matchres <- pmatch(ncp,nf1)
			for (i in 1:length(matchres)) {
				resi <- matchres[i]
				if (!is.na(resi) && resi == 5) k <- i
			}
			if (k != 0) {
				res1 <- pmatch(encodeString(dcall)[k], "FALSE")
			} else {
				res1 <- NA
			}
      
			if (!is.na(res1) && res1 == 1) {
				g <- do.call("s.label", list(dfxy = parse(text=paste(dudiname, "$co", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
			} else {
				g <- do.call("s.corcircle", list(dfxy = parse(text=paste(dudiname, "$co", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
			}
    } else {
			g <- do.call("s.label", list(dfxy = parse(text=paste(dudiname, "$co", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		}
    assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelc1" <- function()
	{
    g <- do.call("s.arrow", list(dfxy = parse(text=paste(dudiname, "$c1", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"corcirclefunc" <- function()
	{
    g <- do.call("s.corcircle",list(dfxy = parse(text=paste(dudiname, "$co", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelc1dpcoa" <- function()
	{
    g <- do.call("s.corcircle",list(dfxy = parse(text=paste(dudiname, "$c1", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

  "dwplotfunc" <- function()
  {
    ynames <- do.call("names", list(parse(text=paste(dudiname, "$dw", sep=""))[[1]]))
    g <- do.call("s1d.dotplot", list(score = parse(text=paste(dudiname, "$dw", sep=""))[[1]], scales = list(y = list(labels = ynames)), 
                                     paxes.draw = TRUE, paxes.y.draw =TRUE, porigin.include = FALSE))
    assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
    if (showCom) {
      print(g@Call)
      cat(substr(options("prompt")$prompt, 1, 2))
    }
    if (histCom) rewriteHistory(deparse(g@Call))
  }

  "lwplotfunc" <- function()
  {
    ynames <- do.call("names", list(parse(text=paste(dudiname, "$lw", sep=""))[[1]]))
    g <- do.call("s1d.dotplot", list(score = parse(text=paste(dudiname, "$lw", sep=""))[[1]], scales = list(y = list(labels = ynames)), 
                                     paxes.draw = TRUE, paxes.y.draw =TRUE, porigin.include = FALSE))
    assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
    if (showCom) {
      print(g@Call)
      cat(substr(options("prompt")$prompt, 1, 2))
    }
    if (histCom) rewriteHistory(deparse(g@Call))
  }

	"RaoDivplotfunc" <- function()
	{
    g <- do.call("s.value",list(dfxy = parse(text=paste(dudiname, "$li", sep=""))[[1]], z = parse(text=paste(dudiname, "$RaoDiv", sep=""))[[1]], 
        xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]], maxsize = 2, psub.text = "Rao Divcs", psub.position = "topright", psub.cex = 1.5))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}
	
  "cwplotfunc" <- function()
  {
    ynames <- do.call("names", list(parse(text=paste(dudiname, "$cw", sep=""))[[1]]))
    g <- do.call("s1d.dotplot", list(score = parse(text=paste(dudiname, "$cw", sep=""))[[1]], scales = list(y = list(labels = ynames)), paxes.draw = TRUE, paxes.y.draw =TRUE, porigin.include = FALSE))
    assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
    if (showCom) {
      print(g@Call)
      cat(substr(options("prompt")$prompt, 1, 2))
    }
    if (histCom) rewriteHistory(deparse(g@Call))
  }

	"eigplotfunc" <- function()
	{
    rank1 <- eval(parse(text=paste(dudiname,"$rank",sep=""))[[1]])
    eigvalue <- eval(parse(text = paste(dudiname, "$eig[1:", rank1, "]", sep="")))
    xax <- parse(text=tclvalue(xaxvar))[[1]]
    yax <- parse(text=tclvalue(yaxvar))[[1]]
    nf <- eval(parse(text=paste(dudiname,"$nf",sep=""))[[1]])
    g <- do.call("plotEig", list(eigvalue = eigvalue, nf = 1:nf, xax = xax, yax = yax, paxes.draw = TRUE, paxes.x.draw = FALSE))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom)
      rewriteHistory(deparse(g@Call))
	}

	"dpcoaeigplotfunc" <- function()
	{
    g <- do.call("plotEig", list(eigvalue = eval(parse(text=paste(dudiname, "$eig", sep=""))), nf = 1:eval(parse(text=paste(dudiname, "$nf", sep="")))))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labells" <- function()
	{
    g <- do.call("s.class",list(dfxy = parse(text=paste(dudiname, "$ls", sep=""))[[1]], fac = parse(text=dcall[3])[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelas" <- function()
	{
    g <- do.call("s.arrow",list(dfxy = parse(text=paste(dudiname, "$as", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelfa" <- function()
	{
    g <- do.call("s.arrow",list(dfxy = parse(text=paste(dudiname, "$fa", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labellid" <- function()
	{
    g <- do.call("s.class",list(dfxy = parse(text=paste(dudiname, "$li", sep=""))[[1]], fac = parse(text=dcall[3])[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelva" <- function()
	{
    g <- do.call("s.corcircle",list(dfxy = parse(text=paste(dudiname, "$va", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelcp" <- function()
	{
    g <- do.call("s.corcircle",list(dfxy = parse(text=paste(dudiname, "$cp", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelgc" <- function()
	{
    g <- do.call("s.label",list(dfxy = parse(text=paste(dudiname, "$gc", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelaX" <- function()
	{
    g <- do.call("s.corcircle",list(dfxy = parse(text=paste(dudiname, "$aX", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelaY" <- function()
	{
    g <- do.call("s.corcircle",list(dfxy = parse(text=paste(dudiname, "$aY", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labell1Y" <- function()
	{
    g <- do.call("s.arrow",list(dfxy = parse(text=paste(dudiname, "$l1", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelc1X" <- function()
	{
    g <- do.call("s.arrow",list(dfxy = parse(text=paste(dudiname, "$c1", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labellX" <- function()
	{
    g <- do.call("s.label",list(dfxy = parse(text=paste(dudiname, "$lX", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labellY" <- function()
	{
    g <- do.call("s.label",list(dfxy = parse(text=paste(dudiname, "$lY", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelmX" <- function()
	{
    g <- do.call("s.match",list(dfxy1 = parse(text=paste(dudiname, "$mY", sep=""))[[1]], dfxy2 = parse(text=paste(dudiname, "$mX", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelmY" <- function()
	{
    g <- do.call("s.match",list(dfxy1 = parse(text=paste(dudiname, "$mX", sep=""))[[1]], dfxy2 = parse(text=paste(dudiname, "$mY", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"tabvalueCA" <- function()
	{
    g <- do.call("table.value",list(dftab = parse(text=paste(dudiname, "$tab", sep=""))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"tabvalueX" <- function()
	{
    g <- do.call("table.value",list(dftab = parse(text=paste(dudiname, "$X", sep=""))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"tabvalueY" <- function()
	{
    g <- do.call("table.value",list(dftab = parse(text=paste(dudiname, "$Y", sep=""))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"tabvalueRaoDecodiv" <- function()
	{
    g <- do.call("table.value",list(dftab = parse(text=paste(dudiname, "$RaoDecodiv", sep=""))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}
	
	"tabvalueRaoDis" <- function()
	{
    g <- do.call("table.value",list(dftab = parse(text=paste(dudiname, "$RaoDis", sep=""))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}
	
	"labellicca" <- function()
	{
    g11 <- do.call("s.match", list(dfxy1 = parse(text=paste(dudiname, "$li", sep=""))[[1]], dfxy2 = parse(text=paste(dudiname, "$ls", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]], plot = FALSE))
    g12 <- do.call("s.label", list(dfxy = parse(text=paste(dudiname, "$co", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]], plabels.cex = 0, ppoints.cex = 2, plot = FALSE))
    g <- do.call("superpose", list(g11, g12, plot = TRUE))
    g@Call <- call("superpose", g11@Call, g12@Call, plot = TRUE)
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labellic" <- function()
	{
    g <- do.call("s.match",list(dfxy1 = parse(text=paste(dudiname, "$li", sep=""))[[1]], dfxy2 = parse(text=paste(dudiname, "$ls", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labellscca" <- function()
	{
    g11 <- do.call("s.match", list(dfxy1 = parse(text=paste(dudiname, "$ls", sep=""))[[1]], dfxy2 = parse(text=paste(dudiname, "$li", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]], plot = FALSE))
    g12 <- do.call("s.label", list(dfxy = parse(text=paste(dudiname, "$co", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]], plabels.cex = 0, ppoints.cex = 2, plot = FALSE))
    g <- do.call("superpose", list(g12, g11, plot = TRUE))
    g@Call <- call("superpose", g11@Call, g12@Call, plot = TRUE)
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labellsc" <- function()
	{
    g <- do.call("s.match",list(dfxy1 = parse(text=paste(dudiname, "$ls", sep=""))[[1]], dfxy2 = parse(text=paste(dudiname, "$li", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelcor" <- function()
	{
    g <- do.call("s.arrow",list(dfxy = parse(text=paste(dudiname, "$cor", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelasc" <- function()
	{
    g <- do.call("s.corcircle",list(dfxy = parse(text=paste(dudiname, "$as", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labelcoc" <- function()
	{
	  g <- do.call("s.label",list(dfxy = parse(text=paste(dudiname, "$co", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

	"labell1c" <- function()
	{
    g <- do.call("s.label",list(dfxy = parse(text=paste(dudiname, "$l1", sep=""))[[1]], xax = parse(text=tclvalue(xaxvar))[[1]], yax = parse(text=tclvalue(yaxvar))[[1]]))
		assign("cmdlist", c(get("cmdlist", envir=env_ade4tkgui), g@Call), envir=env_ade4tkgui)
		if (showCom) {
      print(g@Call)
			cat(substr(options("prompt")$prompt, 1, 2))
		}
		if (histCom) rewriteHistory(deparse(g@Call))
	}

#
# Frame 1 : Window title
#
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)	
	labh <- tklabel(frame1, bitmap="questhead")
	tkbind(labh, "<Button-1>", function() print(help("dudi")))
	tkgrid(tklabel(frame1,text="Duality diagram : summary and graphics", font="Times 18", foreground="red"), labh, columnspan=2)
	tkpack(frame1, fill = "x")
#
# Frame 2 : global dudi properties
#
	frame2 <- tkframe(tf, relief="groove", borderwidth=2)	
	
	dclass <- eval(parse(text=paste("class(",dudiname,")",sep="")))	

	if (dclass[1] == "cca") {
		tkgrid(tklabel(frame2,text="Canonical correspondence analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "pcaiv") {
		tkgrid(tklabel(frame2,text="Principal components analysis", font="Times 16", foreground="blue"), columnspan=2)
		tkgrid(tklabel(frame2,text="with respect to instrumental variables", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "pcaivortho") {
		tkgrid(tklabel(frame2,text="Principal components analysis", font="Times 16", foreground="blue"), columnspan=2)
		tkgrid(tklabel(frame2,text="with respect to orthogonal instrumental variables", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "coinertia") {
		tkgrid(tklabel(frame2,text="Coinertia analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "discrimin") {
		tkgrid(tklabel(frame2,text="Discriminant analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "between") {
		tkgrid(tklabel(frame2,text="Between-class analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "within") {
		tkgrid(tklabel(frame2,text="Within-class analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "pca") {
		tkgrid(tklabel(frame2,text="Principal components analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "coa") {
		tkgrid(tklabel(frame2,text="Correspondence analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "acm") {
		tkgrid(tklabel(frame2,text="Multiple correspondence analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "pco") {
		tkgrid(tklabel(frame2,text="Principal coordinate analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "fca") {
		tkgrid(tklabel(frame2,text="Fuzzy correspondence analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "fpca") {
		tkgrid(tklabel(frame2,text="Fuzzy principal components analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "mix") {
		tkgrid(tklabel(frame2,text="Mixed variables analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "nsc") {
		tkgrid(tklabel(frame2,text="Non symmetric correspondence analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "dec") {
		tkgrid(tklabel(frame2,text="Decentred correspondence analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "hillsmith") {
		tkgrid(tklabel(frame2,text="Hill & Smith analysis", font="Times 16", foreground="blue"), columnspan=2)
	} else if (dclass[1] == "dpcoa") {
		tkgrid(tklabel(frame2,text="Double PCO analysis", font="Times 16", foreground="blue"), columnspan=2)
	}
	
	tclvalue(chvar) <- dclass
	tkgrid(tklabel(frame2,text="Class:"), tklabel(frame2,text=tclvalue(chvar)), sticky="w")

	dcall <- eval(parse(text=paste(dudiname,"$call",sep="")))
	narg <- length(names(dcall))
	paramlst <- encodeString(dcall)[2:narg]
	arglst <- names(dcall)[2:narg]
	call1 <- paste(encodeString(dcall)[1],"(", paste(arglst, paramlst, sep=" = ",collapse=", "), ")", sep="")
	call.label <- tklabel(frame2)
	tkconfigure(call.label, text=call1)
	tkgrid(tklabel(frame2,text="Call:"), call.label, sticky="w")

	nf <- eval(parse(text=paste(dudiname,"$nf",sep="")))
	tclvalue(chvar) <- nf
	tkgrid(tklabel(frame2,text="Axes:"), tklabel(frame2,text=tclvalue(chvar)), sticky="w")

	if (dclass[1] != "discrimin" && dclass[1] != "dpcoa") {
		rank1 <- eval(parse(text=paste(dudiname,"$rank",sep="")))
		tclvalue(chvar) <- rank1
		tkgrid(tklabel(frame2,text="Rank:"), tklabel(frame2,text=tclvalue(chvar)), sticky="w")
	}
	
	if (dclass[1] == "coinertia") {
		rv <- eval(parse(text=paste(dudiname,"$RV",sep="")))
		tclvalue(chvar) <- signif(rv, 4)
		tkgrid(tklabel(frame2,text="RV:"), tklabel(frame2,text=tclvalue(chvar)), sticky="w")
	}
	
	eig <- eval(parse(text=paste(dudiname,"$eig",sep="")))

	if ((dclass[1] == "between") || (dclass[1] == "within")) {
		dudieig <- eval(parse(text=paste(dcall[2],"$eig",sep="")))
		ratio <- signif(sum(eig)/sum(dudieig), 4)
		tclvalue(chvar) <- ratio
		tkgrid(tklabel(frame2,text="Ratio:"), tklabel(frame2,text=tclvalue(chvar)), sticky="w")
	}
	
  l0 <- length(eig)
  eigch <- ""
  eigch <- paste(eigch, signif(eig, 4)[1:(min(5, l0))], sep="")
	if (l0 > 5) eigch[[5]] <- paste(eigch[[5]], "...", sep="")
	tclvalue(chvar) <- eigch
	tkgrid(tklabel(frame2,text="Eigenvalues:"), tklabel(frame2,text=tclvalue(chvar)), sticky="w")

	tkpack(frame2, fill = "x")
#
# Frame 3 : dudi vector components
#
	if (dclass[1] == "cca") {
		frame3 <- tkframe(tf, relief="groove", borderwidth=2)
		
		tkgrid(tklabel(frame3,text="  "), tklabel(frame3,text="Vectors:", foreground="blue"), tklabel(frame3,text="Length:", foreground="blue"),
			tklabel(frame3,text="Mode:", foreground="blue"), tklabel(frame3,text="Content:", foreground="blue"), sticky="w")
	  
		ne <- eval(parse(text=paste("length(",dudiname,"$cw)",sep="")))	
		tclvalue(chvar) <- ne
		cwplot.but <- tkbutton(frame3, text=paste(dudiname, "$cw", sep=""), anchor="w", command=function() cwplotfunc())
		tkgrid(tklabel(frame3,text="1:"), cwplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="column weights (from dudi)"), sticky="w")
	  
		ne <- eval(parse(text=paste("length(",dudiname,"$lw)",sep="")))	
		tclvalue(chvar) <- ne
		lwplot.but <- tkbutton(frame3, text=paste(dudiname, "$lw", sep=""), anchor="w", command=function() lwplotfunc())
		tkgrid(tklabel(frame3,text="2:"), lwplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="row weights (from dudi)"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$eig)",sep="")))	
		tclvalue(chvar) <- ne
		eigplot.but <- tkbutton(frame3, text=paste(dudiname, "$eig", sep=""), anchor="w", command=function() eigplotfunc())
		tkgrid(tklabel(frame3,text="3:"), eigplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="eigenvalues"), sticky="w")
	
		tkpack(frame3, fill = "x")
	} else if (dclass[1] == "coinertia") {
		frame3 <- tkframe(tf, relief="groove", borderwidth=2)
		
		tkgrid(tklabel(frame3,text="  "), tklabel(frame3,text="Vectors:", foreground="blue"), tklabel(frame3,text="Length:", foreground="blue"),
			tklabel(frame3,text="Mode:", foreground="blue"), tklabel(frame3,text="Content:", foreground="blue"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$cw)",sep="")))	
		tclvalue(chvar) <- ne
		cwplot.but <- tkbutton(frame3, text=paste(dudiname, "$cw", sep=""), anchor="w", command=function() cwplotfunc())
		tkgrid(tklabel(frame3,text="1:"), cwplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="column weights (crossed array)"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$lw)",sep="")))	
		tclvalue(chvar) <- ne
		lwplot.but <- tkbutton(frame3, text=paste(dudiname, "$lw", sep=""), anchor="w", command=function() lwplotfunc())
		tkgrid(tklabel(frame3,text="2:"), lwplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="row weights (crossed array)"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$eig)",sep="")))	
		tclvalue(chvar) <- ne
		eigplot.but <- tkbutton(frame3, text=paste(dudiname, "$eig", sep=""), anchor="w", command=function() eigplotfunc())
		tkgrid(tklabel(frame3,text="3:"), eigplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="eigenvalues"), sticky="w")
	
		tkpack(frame3, fill = "x")
	} else if (dclass[1] == "dpcoa") {
		frame3 <- tkframe(tf, relief="groove", borderwidth=2)
		
		tkgrid(tklabel(frame3,text="  "), tklabel(frame3,text="Vectors:", foreground="blue"), tklabel(frame3,text="Length:", foreground="blue"),
			tklabel(frame3,text="Mode:", foreground="blue"), tklabel(frame3,text="Content:", foreground="blue"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$dw)",sep="")))	
		tclvalue(chvar) <- ne
		cwplot.but <- tkbutton(frame3, text=paste(dudiname, "$dw", sep=""), anchor="w", command=function() dwplotfunc())
		tkgrid(tklabel(frame3,text="1:"), cwplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="weights of species"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$lw)",sep="")))	
		tclvalue(chvar) <- ne
		lwplot.but <- tkbutton(frame3, text=paste(dudiname, "$lw", sep=""), anchor="w", command=function() lwplotfunc())
		tkgrid(tklabel(frame3,text="2:"), lwplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="weights of communities"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$eig)",sep="")))	
		tclvalue(chvar) <- ne
		eigplot.but <- tkbutton(frame3, text=paste(dudiname, "$eig", sep=""), anchor="w", command=function() dpcoaeigplotfunc())
		tkgrid(tklabel(frame3,text="3:"), eigplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="eigenvalues"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$RaoDiv)",sep="")))	
		tclvalue(chvar) <- ne
		eigplot.but <- tkbutton(frame3, text=paste(dudiname, "$RaoDiv", sep=""), anchor="w", command=function() RaoDivplotfunc())
		tkgrid(tklabel(frame3,text="4:"), eigplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="diversity coefficients within communities"), sticky="w")
	
		tkpack(frame3, fill = "x")
	} else if (dclass[1] != "discrimin") {
		frame3 <- tkframe(tf, relief="groove", borderwidth=2)
		
		tkgrid(tklabel(frame3,text="  "), tklabel(frame3,text="Vectors:", foreground="blue"), tklabel(frame3,text="Length:", foreground="blue"),
			tklabel(frame3,text="Mode:", foreground="blue"), tklabel(frame3,text="Content:", foreground="blue"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$cw)",sep="")))	
		tclvalue(chvar) <- ne
		cwplot.but <- tkbutton(frame3, text=paste(dudiname, "$cw", sep=""), anchor="w", command=function() cwplotfunc())
		tkgrid(tklabel(frame3,text="1:"), cwplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="column weights"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$lw)",sep="")))	
		tclvalue(chvar) <- ne
		lwplot.but <- tkbutton(frame3, text=paste(dudiname, "$lw", sep=""), anchor="w", command=function() lwplotfunc())
		tkgrid(tklabel(frame3,text="2:"), lwplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="row weights"), sticky="w")
	
		ne <- eval(parse(text=paste("length(",dudiname,"$eig)",sep="")))	
		tclvalue(chvar) <- ne
		eigplot.but <- tkbutton(frame3, text=paste(dudiname, "$eig", sep=""), anchor="w", command=function() eigplotfunc())
		tkgrid(tklabel(frame3,text="3:"), eigplot.but, tklabel(frame3,text=tclvalue(chvar)),
			tklabel(frame3,text="numeric"), tklabel(frame3,text="eigenvalues"), sticky="w")
	
		tkpack(frame3, fill = "x")
	}
#
# Frame 3b : dudi dist components
#
	if (dclass[1] == "dpcoa") {
		frame3b <- tkframe(tf, relief="groove", borderwidth=2)
		
		tkgrid(tklabel(frame3b,text="  "), tklabel(frame3b,text="dist:", foreground="blue"), tklabel(frame3b,text="Size:", foreground="blue"),
			tklabel(frame3b,text="Content:", foreground="blue"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$RaoDis)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$RaoDis)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelfa.but <- tkbutton(frame3b, text=paste(dudiname, "$RaoDis", sep=""), anchor="w", command=function() tabvalueRaoDis())
		tkgrid(tklabel(frame3b,text="1:"), labelfa.but, tklabel(frame3b,text=tclvalue(chvar)),
			tklabel(frame3b,text="dissimilarities among communities"), sticky="w")

		tkpack(frame3b, fill = "x")
	}
#
# Frame 4 : dudi dataframe components
#
	frame4 <- tkframe(tf, relief="groove", borderwidth=2)
	
	tkgrid(tklabel(frame4,text="  "), tklabel(frame4,text="Dataframes:", foreground="blue"), tklabel(frame4,text="Nrow:", foreground="blue"),
		tklabel(frame4,text="Ncol:", foreground="blue"), tklabel(frame4,text="Content:", foreground="blue"), sticky="w")

	if (dclass[1] == "discrimin") {
		nr <- eval(parse(text=paste("dim(",dudiname,"$fa)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$fa)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelfa.but <- tkbutton(frame4, text=paste(dudiname, "$fa", sep=""), anchor="w", command=function() labelfa())
		tkgrid(tklabel(frame4,text="1:"), labelfa.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="loadings / canonical weights"), sticky="w")

		nr <- eval(parse(text=paste("dim(",dudiname,"$li)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$li)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labellid.but <- tkbutton(frame4, text=paste(dudiname, "$li", sep=""), anchor="w", command=function() labellid())
		tkgrid(tklabel(frame4,text="2:"), labellid.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="canonical scores"), sticky="w")

		nr <- eval(parse(text=paste("dim(",dudiname,"$va)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$va)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelva.but <- tkbutton(frame4, text=paste(dudiname, "$va", sep=""), anchor="w", command=function() labelva())
		tkgrid(tklabel(frame4,text="3:"), labelva.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="cos(variables, canonical scores)"), sticky="w")

		nr <- eval(parse(text=paste("dim(",dudiname,"$cp)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$cp)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelcp.but <- tkbutton(frame4, text=paste(dudiname, "$cp", sep=""), anchor="w", command=function() labelcp())
		tkgrid(tklabel(frame4,text="4:"), labelcp.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="cos(components, canonical scores)"), sticky="w")

		nr <- eval(parse(text=paste("dim(",dudiname,"$gc)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$gc)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelgc.but <- tkbutton(frame4, text=paste(dudiname, "$gc", sep=""), anchor="w", command=function() labelgc())
		tkgrid(tklabel(frame4,text="5:"), labelgc.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="class scores"), sticky="w")
	} else if (dclass[1] == "cca") {
		nr <- eval(parse(text=paste("dim(",dudiname,"$tab)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$tab)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		tab.but <- tkbutton(frame4, text=paste(dudiname, "$tab", sep=""), anchor="w", command=function() tabvalueCA())
		tkgrid(tklabel(frame4,text="1a:"), tab.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="modified array (projected variables)"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$Y)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$Y)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelY.but <- tkbutton(frame4, text=paste(dudiname, "$Y", sep=""), anchor="w", command=function() tabvalueY())
		tkgrid(tklabel(frame4,text="1b:"), labelY.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="Dependent variables"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$X)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$X)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelX.but <- tkbutton(frame4, text=paste(dudiname, "$X", sep=""), anchor="w", command=function() tabvalueX())
		tkgrid(tklabel(frame4,text="1c:"), labelX.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="Explanatory variables"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$c1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$c1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelc1.but <- tkbutton(frame4, text=paste(dudiname, "$c1", sep=""), anchor="w", command=function() labelc1())
		tkgrid(tklabel(frame4,text="2a:"), labelc1.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="PPA : Pseudo Principal Axes"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$as)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$as)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelasc.but <- tkbutton(frame4, text=paste(dudiname, "$as", sep=""), anchor="w", command=function() labelasc())
		tkgrid(tklabel(frame4,text="2b:"), labelasc.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="Principal axis of dudi$tab on PPA"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$ls)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$ls)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labellscca.but <- tkbutton(frame4, text=paste(dudiname, "$ls", sep=""), anchor="w", command=function() labellscca())
		tkgrid(tklabel(frame4,text="2c:"), labellscca.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="projection of the rows of dudi$tab on PPA"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$li)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$li)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labellicca.but <- tkbutton(frame4, text=paste(dudiname, "$li", sep=""), anchor="w", command=function() labellicca())
		tkgrid(tklabel(frame4,text="2d:"), labellicca.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="$ls predicted by X"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$l1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$l1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labell1c.but <- tkbutton(frame4, text=paste(dudiname, "$l1", sep=""), anchor="w", command=function() labell1c())
		tkgrid(tklabel(frame4,text="3a:"), labell1c.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="CPC Constraint Principal Components"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$fa)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$fa)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelfa.but <- tkbutton(frame4, text=paste(dudiname, "$fa", sep=""), anchor="w", command=function() labelfa())
		tkgrid(tklabel(frame4,text="3b:"), labelfa.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="Loadings (CPC as linear combinations of X)"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$co)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$co)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelcoc.but <- tkbutton(frame4, text=paste(dudiname, "$co", sep=""), anchor="w", command=function() labelcoc())
		tkgrid(tklabel(frame4,text="3c:"), labelcoc.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="inner product CPC - Y"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$cor)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$cor)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelcor.but <- tkbutton(frame4, text=paste(dudiname, "$cor", sep=""), anchor="w", command=function() labelcor())
		tkgrid(tklabel(frame4,text="3d:"), labelcor.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="correlation CPC - X"), sticky="w")
	
	} else if (dclass[1] == "pcaiv") {
		nr <- eval(parse(text=paste("dim(",dudiname,"$tab)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$tab)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		tab.but <- tkbutton(frame4, text=paste(dudiname, "$tab", sep=""), anchor="w", command=function() tabvalueCA())
		tkgrid(tklabel(frame4,text="1a:"), tab.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="modified array (projected variables)"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$Y)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$Y)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelY.but <- tkbutton(frame4, text=paste(dudiname, "$Y", sep=""), anchor="w", command=function() tabvalueY())
		tkgrid(tklabel(frame4,text="1b:"), labelY.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="Dependent variables"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$X)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$X)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelX.but <- tkbutton(frame4, text=paste(dudiname, "$X", sep=""), anchor="w", command=function() tabvalueX())
		tkgrid(tklabel(frame4,text="1c:"), labelX.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="Explanatory variables"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$c1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$c1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelc1.but <- tkbutton(frame4, text=paste(dudiname, "$c1", sep=""), anchor="w", command=function() labelc1())
		tkgrid(tklabel(frame4,text="2a:"), labelc1.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="PPA : Pseudo Principal Axes"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$as)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$as)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelasc.but <- tkbutton(frame4, text=paste(dudiname, "$as", sep=""), anchor="w", command=function() labelasc())
		tkgrid(tklabel(frame4,text="2b:"), labelasc.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="Principal axis of dudi$tab on PPA"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$ls)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$ls)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labellsc.but <- tkbutton(frame4, text=paste(dudiname, "$ls", sep=""), anchor="w", command=function() labellsc())
		tkgrid(tklabel(frame4,text="2c:"), labellsc.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="projection of the rows of dudi$tab on PPA"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$li)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$li)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labellic.but <- tkbutton(frame4, text=paste(dudiname, "$li", sep=""), anchor="w", command=function() labellic())
		tkgrid(tklabel(frame4,text="2d:"), labellic.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="$ls predicted by X"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$l1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$l1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labell1c.but <- tkbutton(frame4, text=paste(dudiname, "$l1", sep=""), anchor="w", command=function() labell1c())
		tkgrid(tklabel(frame4,text="3a:"), labell1c.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="CPC Constraint Principal Components"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$fa)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$fa)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelfa.but <- tkbutton(frame4, text=paste(dudiname, "$fa", sep=""), anchor="w", command=function() labelfa())
		tkgrid(tklabel(frame4,text="3b:"), labelfa.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="Loadings (CPC as linear combinations of X)"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$co)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$co)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelcoc.but <- tkbutton(frame4, text=paste(dudiname, "$co", sep=""), anchor="w", command=function() labelcoc())
		tkgrid(tklabel(frame4,text="3c:"), labelcoc.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="inner product CPC - Y"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$cor)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$cor)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelcor.but <- tkbutton(frame4, text=paste(dudiname, "$cor", sep=""), anchor="w", command=function() labelcor())
		tkgrid(tklabel(frame4,text="3d:"), labelcor.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="correlation CPC - X"), sticky="w")
	
	} else if (dclass[1] == "pcaivortho") {
		nr <- eval(parse(text=paste("dim(",dudiname,"$tab)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$tab)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		tab.but <- tkbutton(frame4, text=paste(dudiname, "$tab", sep=""), anchor="w", command=function() tabvalue())
		tkgrid(tklabel(frame4,text="1:"), tab.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="modified array"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$li)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$li)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelli.but <- tkbutton(frame4, text=paste(dudiname, "$li", sep=""), anchor="w", command=function() labelli())
		tkgrid(tklabel(frame4,text="2:"), labelli.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="row coordinates"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$l1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$l1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labell1.but <- tkbutton(frame4, text=paste(dudiname, "$l1", sep=""), anchor="w", command=function() labell1())
		tkgrid(tklabel(frame4,text="3:"), labell1.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="row normed scores"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$co)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$co)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelco.but <- tkbutton(frame4, text=paste(dudiname, "$co", sep=""), anchor="w", command=function() labelco())
		tkgrid(tklabel(frame4,text="4:"), labelco.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="column coordinates"), sticky="w")
		
		nr <- eval(parse(text=paste("dim(",dudiname,"$c1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$c1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelc1.but <- tkbutton(frame4, text=paste(dudiname, "$c1", sep=""), anchor="w", command=function() labelc1())
		tkgrid(tklabel(frame4,text="5:"), labelc1.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="column normed scores"), sticky="w")

	} else if (dclass[1] == "dpcoa") {
		nr <- eval(parse(text=paste("dim(",dudiname,"$dls)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$dls)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		tab.but <- tkbutton(frame4, text=paste(dudiname, "$dls", sep=""), anchor="w", command=function() labeldlsdpcoa())
		tkgrid(tklabel(frame4,text="1:"), tab.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="coordinates of the species"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$li)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$li)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		tab.but <- tkbutton(frame4, text=paste(dudiname, "$li", sep=""), anchor="w", command=function() labellidpcoa())
		tkgrid(tklabel(frame4,text="2:"), tab.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="coordinates of the communities"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$c1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$c1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		tab.but <- tkbutton(frame4, text=paste(dudiname, "$c1", sep=""), anchor="w", command=function() labelc1dpcoa())
		tkgrid(tklabel(frame4,text="3:"), tab.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="scores of the principal axes of the species"), sticky="w")
	
	} else if (dclass[1] == "coinertia") {
		nr <- eval(parse(text=paste("dim(",dudiname,"$tab)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$tab)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		tab.but <- tkbutton(frame4, text=paste(dudiname, "$tab", sep=""), anchor="w", command=function() tabvalueCA())
		tkgrid(tklabel(frame4,text="1:"), tab.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="crossed array (CA)"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$li)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$li)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelli.but <- tkbutton(frame4, text=paste(dudiname, "$li", sep=""), anchor="w", command=function() labelli())
		tkgrid(tklabel(frame4,text="2:"), labelli.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="Y col = CA row: coordinates"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$l1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$l1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labell1.but <- tkbutton(frame4, text=paste(dudiname, "$l1", sep=""), anchor="w", command=function() labell1Y())
		tkgrid(tklabel(frame4,text="3:"), labell1.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="Y col = CA row: normed scores"), sticky="w")
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$co)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$co)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelco.but <- tkbutton(frame4, text=paste(dudiname, "$co", sep=""), anchor="w", command=function() labelco())
		tkgrid(tklabel(frame4,text="4:"), labelco.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="X col = CA column: coordinates"), sticky="w")
		
		nr <- eval(parse(text=paste("dim(",dudiname,"$c1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$c1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelc1.but <- tkbutton(frame4, text=paste(dudiname, "$c1", sep=""), anchor="w", command=function() labelc1X())
		tkgrid(tklabel(frame4,text="5:"), labelc1.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="X col = CA column: normed scores"), sticky="w")
		
		nr <- eval(parse(text=paste("dim(",dudiname,"$lX)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$lX)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labellix.but <- tkbutton(frame4, text=paste(dudiname, "$lX", sep=""), anchor="w", command=function() labellX())
		tkgrid(tklabel(frame4,text="6:"), labellix.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="row coordinates (X)"), sticky="w")

		nr <- eval(parse(text=paste("dim(",dudiname,"$mX)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$mX)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labell1x.but <- tkbutton(frame4, text=paste(dudiname, "$mX", sep=""), anchor="w", command=function() labelmX())
		tkgrid(tklabel(frame4,text="7:"), labell1x.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="normed row scores (X)"), sticky="w")

		nr <- eval(parse(text=paste("dim(",dudiname,"$lY)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$lY)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelliy.but <- tkbutton(frame4, text=paste(dudiname, "$lY", sep=""), anchor="w", command=function() labellY())
		tkgrid(tklabel(frame4,text="8:"), labelliy.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="row coordinates (Y)"), sticky="w")

		nr <- eval(parse(text=paste("dim(",dudiname,"$mY)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$mY)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labell1y.but <- tkbutton(frame4, text=paste(dudiname, "$mY", sep=""), anchor="w", command=function() labelmY())
		tkgrid(tklabel(frame4,text="9:"), labell1y.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="normed row scores (Y)"), sticky="w")

		nr <- eval(parse(text=paste("dim(",dudiname,"$aX)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$aX)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelax.but <- tkbutton(frame4, text=paste(dudiname, "$aX", sep=""), anchor="w", command=function() labelaX())
		tkgrid(tklabel(frame4,text="10:"), labelax.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="inertia axes onto coinertia axes (X)"), sticky="w")

		nr <- eval(parse(text=paste("dim(",dudiname,"$aY)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$aY)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelay.but <- tkbutton(frame4, text=paste(dudiname, "$aY", sep=""), anchor="w", command=function() labelaY())
		tkgrid(tklabel(frame4,text="11:"), labelay.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="inertia axes onto coinertia axes (Y)"), sticky="w")
	} else {
		nr <- eval(parse(text=paste("dim(",dudiname,"$tab)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$tab)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		tab.but <- tkbutton(frame4, text=paste(dudiname, "$tab", sep=""), anchor="w", command=function() tabvalue())
		if ((dclass[1] == "between") || (dclass[1] == "within")) {
			tkgrid(tklabel(frame4,text="1:"), tab.but, tklabel(frame4,text=tclvalue(chvar)),
				tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="classes x variables array"), sticky="w")
		} else {
			tkgrid(tklabel(frame4,text="1:"), tab.but, tklabel(frame4,text=tclvalue(chvar)),
				tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="modified array"), sticky="w")
		}
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$li)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$li)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelli.but <- tkbutton(frame4, text=paste(dudiname, "$li", sep=""), anchor="w", command=function() labelli())
		if (dclass[1] == "between") {
			tkgrid(tklabel(frame4,text="2:"), labelli.but, tklabel(frame4,text=tclvalue(chvar)),
				tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="class coordinates"), sticky="w")
		} else {
			tkgrid(tklabel(frame4,text="2:"), labelli.but, tklabel(frame4,text=tclvalue(chvar)),
				tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="row coordinates"), sticky="w")
		}
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$l1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$l1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labell1.but <- tkbutton(frame4, text=paste(dudiname, "$l1", sep=""), anchor="w", command=function() labell1())
		if (dclass[1] == "between") {
			tkgrid(tklabel(frame4,text="3:"), labell1.but, tklabel(frame4,text=tclvalue(chvar)),
				tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="class normed scores"), sticky="w")
		} else {
			tkgrid(tklabel(frame4,text="3:"), labell1.but, tklabel(frame4,text=tclvalue(chvar)),
				tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="row normed scores"), sticky="w")
		}
	
		nr <- eval(parse(text=paste("dim(",dudiname,"$co)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$co)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelco.but <- tkbutton(frame4, text=paste(dudiname, "$co", sep=""), anchor="w", command=function() labelco())
		tkgrid(tklabel(frame4,text="4:"), labelco.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="column coordinates"), sticky="w")
		
		nr <- eval(parse(text=paste("dim(",dudiname,"$c1)[1]",sep="")))	
		tclvalue(chvar) <- nr
		nc <- eval(parse(text=paste("dim(",dudiname,"$c1)[2]",sep="")))	
		tclvalue(chvar2) <- nc
		labelc1.but <- tkbutton(frame4, text=paste(dudiname, "$c1", sep=""), anchor="w", command=function() labelc1())
		tkgrid(tklabel(frame4,text="5:"), labelc1.but, tklabel(frame4,text=tclvalue(chvar)),
			tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="column normed scores"), sticky="w")
		
		if ((dclass[1] == "between") || ((dclass[1] == "within"))) {
			nr <- eval(parse(text=paste("dim(",dudiname,"$ls)[1]",sep="")))	
			tclvalue(chvar) <- nr
			nc <- eval(parse(text=paste("dim(",dudiname,"$ls)[2]",sep="")))	
			tclvalue(chvar2) <- nc
			labells.but <- tkbutton(frame4, text=paste(dudiname, "$ls", sep=""), anchor="w", command=function() labells())
			tkgrid(tklabel(frame4,text="6:"), labells.but, tklabel(frame4,text=tclvalue(chvar)),
				tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="supplementary row coordinates"), sticky="w")
	
			nr <- eval(parse(text=paste("dim(",dudiname,"$as)[1]",sep="")))	
			tclvalue(chvar) <- nr
			nc <- eval(parse(text=paste("dim(",dudiname,"$as)[2]",sep="")))	
			tclvalue(chvar2) <- nc
			labelas.but <- tkbutton(frame4, text=paste(dudiname, "$as", sep=""), anchor="w", command=function() labelas())
			tkgrid(tklabel(frame4,text="7:"), labelas.but, tklabel(frame4,text=tclvalue(chvar)),
				tklabel(frame4,text=tclvalue(chvar2)), tklabel(frame4,text="projection of inertia axes"), sticky="w")
		}
	}
	tkpack(frame4, fill = "x")
#
# Frame 5 : X and Y axes choice for factor map graphics
#
	frame5 <- tkframe(tf, relief="groove", borderwidth=2)
	xax.entry <- tkentry(frame5, textvariable=xaxvar, width=4)
	yax.entry <- tkentry(frame5, textvariable=yaxvar, width=4)
	tkgrid(tklabel(frame5,text="Select X and Y axis numbers for graphics - "),
		tklabel(frame5,text="X axis : "), xax.entry, tklabel(frame5,text="Y axis : "), yax.entry)
	tkpack(frame5, fill = "x")
#
# Frame 6 : pcaiv params
#
	if ((dclass[1] == "cca") || (dclass[1] == "pcaiv") || (dclass[1] == "pcaivortho")) {
		frame6 <- tkframe(tf, relief="groove", borderwidth=2)
		txt <- tktext(frame6, bg="white", font="courier", width=50, height=2+eval(parse(text=paste(dudiname,"$nf",sep=""))))
		tkpack(txt, side="left", fill="both", expand=TRUE)
		cheval <- paste("summary(", dudiname, ")", sep="")
    sink(conn <- file("Rlisting001.tmp", open="w"))
		cat("Global analysis parameters :\n")
		print(eval(parse(text=cheval)))
		sink()
		close(conn)
		chn <- tclopen(file.path("Rlisting001.tmp"))
		tkinsert(txt, "end", tclread(chn))
		tclclose(chn)
		system("rm Rlisting001.tmp")
		tkconfigure(txt, state="disabled")
		tkmark.set(txt,"insert","0.0")
		tkfocus(txt)
		tkpack(frame6, fill = "x")
	} else if (dclass[1] == "dpcoa") {
		frame6 <- tkframe(tf, relief="groove", borderwidth=2)
		txt <- tktext(frame6, bg="white", font="courier", width=50, height=5)
		tkpack(txt, side="left", fill="both", expand=TRUE)
		sink(conn <- file("Rlisting001.tmp", open="w"))
		cheval <- paste(dudiname, "$RaoDecodiv", sep="")
		cat("Global analysis parameters :\n")
		print(eval(parse(text=cheval)))
		sink()
		close(conn)
		chn <- tclopen(file.path("Rlisting001.tmp"))
		tkinsert(txt, "end", tclread(chn))
		tclclose(chn)
		system("rm Rlisting001.tmp")
		tkconfigure(txt, state="disabled")
		tkmark.set(txt,"insert","0.0")
		tkfocus(txt)
		tkpack(frame6, fill = "x")
	}
#
# Last row buttons
#
	OK.but <- tkbutton(tf, text="Dismiss", command=function() tkdestroy(tf))
	scatter.but <- tkbutton(tf, text=paste("scatter(",dudiname,")",sep=""), default="active", command=function() scatterfunc())
	score.but <- tkbutton(tf, text=paste("score(",dudiname,")",sep=""), command=function() scorefunc())
	plot.but <- tkbutton(tf, text=paste("plot(",dudiname,")",sep=""), default="active", command=function() plotfunc())
	k <- 0
	if (dclass[1] == "pca") {
		ncp <- names(dcall)
		nf1 <- names(formals(dudi.pca))
		matchres <- pmatch(ncp,nf1)
		for (i in 1:length(matchres)) {
			resi <- matchres[i]
			if (!is.na(resi) && resi == 5) k <- i
		}
		if (k != 0) {
			res1 <- pmatch(encodeString(dcall)[k], "FALSE")
		} else {
			res1 <- NA
		}
		if (!is.na(res1) && res1 == 1) {
			tkpack(OK.but, scatter.but, score.but, side="left", expand=1, fill = "x")	
		} else {
			corcircle.but <- tkbutton(tf, text=paste("s.corcircle(",dudiname,")",sep=""), command=function() corcirclefunc())
			tkpack(OK.but, scatter.but, score.but, corcircle.but, side="left", expand=1, fill = "x")	
		}
	} else if (dclass[1] == "coa") {
    # tkpack(OK.but, scatter.but, score.but, side="left", expand=1, fill = "x")
    # the score button is deleted until score.coa was implemented in adegraphics
	  tkpack(OK.but, scatter.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "acm") {
		tkpack(OK.but, plot.but, score.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "pco") {
		tkpack(OK.but, scatter.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "fca") {
		tkpack(OK.but, plot.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "fcpa") {
		tkpack(OK.but, scatter.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "mix") {
		tkpack(OK.but, scatter.but, score.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "nsc") {
		tkpack(OK.but, scatter.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "dec") {
		tkpack(OK.but, scatter.but, side="left", expand=1, fill = "x")	
	} else if ((dclass[1] == "coinertia") || ((dclass[1] == "within"))) {
		tkpack(OK.but, plot.but, scatter.but, side="left", expand=1, fill = "x")	
	} else if ((dclass[1] == "between") || ((dclass[1] == "within"))) {
		tkpack(OK.but, plot.but, scatter.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "discrimin") {
		tkpack(OK.but, plot.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "cca") {
		tkpack(OK.but, plot.but, scatter.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "pcaiv") {
		tkpack(OK.but, plot.but, scatter.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "pcaivortho") {
		tkpack(OK.but, scatter.but, side="left", expand=1, fill = "x")	
	} else if (dclass[1] == "dpcoa") {
		tkpack(OK.but, plot.but, side="left", expand=1, fill = "x")	
	}
#
# Ending the dialog
#
	tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
	tkbind(tf, "<KeyPress-Return>", function() scatterfunc())
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	if(tclvalue(done) == "2")
    return(0)
}

