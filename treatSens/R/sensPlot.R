#############
#sensPlot - plot results of sensitivity analysis
#############

sensPlot = function(x, 
                  contour.levels = NULL,
                  col.zero = "red",
                  lty.zero = 1,
                  col.insig = "blue",
                  lty.insig = 1,
                  data.line = TRUE,
                  X.pch = NULL,  #vector of length 3: non-transformed points, transformed points, points outside plotting range
                  signif.level = 0.05,
                  labcex = 0.75,
                  limit.Xplot = FALSE, #MH: limit plotting covariates to enlarge contour
                  txtlab = FALSE,  #add text label to the plots of covariates.
                  which.txtlab = NULL, #enter numeric vector to specify which label to show. e.g. c(1:3) shows first 3 covariates.
                  ...) {
  #note in help: if contours are too rough, up nsim in sens fn
  if(class(x) == "sensitivity"){  
    sensPlotMain(x, contour.levels, col.zero, lty.zero, col.insig, lty.insig, data.line, X.pch, signif.level, labcex, limit.Xplot, txtlab, which.txtlab,...)
  }else if(class(x) == "sensitivityCombo"){
    sensPlotCombo(x, contour.levels, col.zero, lty.zero, col.insig, lty.insig, data.line, X.pch, signif.level, labcex, limit.Xplot, txtlab, which.txtlab,...)
  }else{
    stop("x must be of class sensitivity or sensitivity.combo")
  }  
} #end of sensPlot


#############
#sensPlotMain
#Plot for results from single run
#############

sensPlotMain = function(x, contour.levels, col.zero, lty.zero, col.insig, lty.insig, data.line, X.pch, signif.level, labcex, limit.Xplot, txtlab, which.txtlab,...){
  ##Add row/column for zeta = 0 if not included in grid
  null.tau=x$tau0
  null.se=x$se.tau0
  part.cors = x$sensParam == "cor"
  
  Zcors = as.numeric(dimnames(x$sp.z)[[2]]) #horizontal grids of U
  Ycors = as.numeric(dimnames(x$sp.y)[[1]]) #vertical grids of U
  
  if(part.cors){
    print.text <- "Partial correlations with U"
  }else {
    print.text <- "Coefficients on U" 
  }
  
  ##############  
  #######Only inclue X on the plot for GLM-style
  ############
  Xpart = x$Xcoef[!is.na(x$Xcoef[,1]) & !is.na(x$Xcoef[,2]),] #coefficients of null model.  
  Xpart.plot = x$Xcoef.plot[!is.na(x$Xcoef.plot[,1]) & !is.na(x$Xcoef.plot[,2]),]
  Xpart.plot2 = cbind(Xpart.plot[,1],Xpart.plot[,2], ifelse(Xpart[,2]>=0,1,2)) #MH: add sign of coef of X on Y to Xpart  
    
  #note that due to correlation among Xs, some may not appear on plot
  #because observed partial cors don't map directly to coefs in this case
  #forcing inclusion can lead to difficult to read plot.  
  if(is.null(X.pch)){
    out.pch = 1
    X.pch = ifelse(Xpart[,2]>=0,3,6) # plus sign for non-transformed plots, reverse triangle for transformed plots
  }else{
    out.pch = X.pch[3]
    X.pch = ifelse(Xpart[,2]>=0,X.pch[1],X.pch[2])	
  }
  if (any(Xpart[,2]<0)) {
    cat("Note: Predictors with negative coefficients for the response surface have been transformed through multiplication by -1 and are displayed as inverted triangles.", "\n")    
  }
    
  
  nr = length(Zcors); nc = length(Ycors)
  taus.est = t(apply(x$tau, c(1,2), mean, na.rm = T))
  K = dim(x$se.tau)[3]
  W = apply(x$se.tau^2, c(1,2), mean, na.rm = T)
  B = apply(x$tau, c(1,2), sd, na.rm = T)^2
  se.est = t(sqrt(W+(1+1/K)*B))
  taus = taus.est
  se.taus = se.est
  taus[Zcors==0,] = null.tau
  taus[,Ycors==0] = null.tau
  
  dimnames(taus)[[1]] = Zcors
  dimnames(taus)[[2]] = Ycors
  taus <<- taus ## is this right? it potentially changes the value of tau in a parent environment
  dimnames(se.taus)[[1]] = Zcors
  dimnames(se.taus)[[2]] = Ycors
  se.taus <<- se.taus
  
  if(is.null(contour.levels)){
    exTau = c(taus[dim(taus)[1], dim(taus)[2]], taus[1,dim(taus)[2]]) #extreme values of tau at right end
    clevels = round(seq(exTau[2]*.8, exTau[1]*.8, length.out = 14), 2) #vals at which contours are drawn
  }else{
    clevels = contour.levels
  } 
  
  par(mgp = c(2,.5,0)) #dist of axis label, tick mark label, tick mark
  if(part.cors){
    xlab = expression(paste("Partial cor. with U in model for treatment, ", rho^zu))
    ylab = expression(paste("Partial cor. with U in model for response, ", rho^yu))
  }else{
    xlab = expression(paste("Coef. on U in model for treatment, ", zeta^z))
    ylab = expression(paste("Coef. on U in model for response, ", zeta^y))
  }
  
  if (limit.Xplot) {
    #old codes
    plot(Xpart.plot2[,1], Xpart.plot2[,2], col=c("blue","red")[Xpart.plot2[,3]], xlim = c(min(Zcors, na.rm = T),max(Zcors, na.rm = T)), 
         ylim = c(min(Ycors, na.rm = T),max(Ycors, na.rm = T)), pch = X.pch, xlab = xlab, ylab = ylab)
    outsidePts = Xpart.plot2[(Xpart.plot[,1] < min(Zcors, na.rm = T)) | (Xpart.plot[,1] > max(Zcors, na.rm = T)) | (Xpart.plot[,2] > max(Ycors, na.rm = T)),]
    if(length(outsidePts) > 0){
      outsidePts[,1] = apply(cbind(outsidePts[,1], min(Zcors, na.rm = T)), 1, max)
      outsidePts[,1] = apply(cbind(outsidePts[,1], max(Zcors, na.rm = T)), 1, min)
      outsidePts[,2] = apply(cbind(outsidePts[,2], max(Ycors, na.rm = T)), 1, min)
      points(outsidePts[,1], outsidePts[,2], col=c("blue","red")[outsidePts[,3]], pch = out.pch, cex = 1.5, lwd = 3)
      warning("Note: predictors outside plot region plotted at margin with O")
    }
  } else {
    #MH: define max, min of plots
    xplot.min = ifelse(min(Zcors, na.rm = T)<min(Xpart.plot[,1]),min(Zcors, na.rm = T),min(Xpart.plot[,1]))
    xplot.max = ifelse(max(Zcors, na.rm = T)>max(Xpart.plot[,1]),max(Zcors, na.rm = T),max(Xpart.plot[,1]))
    yplot.max = ifelse(max(Ycors, na.rm = T)>max(Xpart.plot[,2]),max(Ycors, na.rm = T),max(Xpart.plot[,2]))  
    plot(Xpart.plot2[,1], Xpart.plot2[,2], col=c("red","blue")[as.factor(Xpart.plot2[,3])], xlim = c(xplot.min,xplot.max),
         ylim = c(0,yplot.max), pch = X.pch, xlab = xlab, ylab = ylab)
  }
  
  #codes for txtlab
  if (txtlab) {
    if (is.null(which.txtlab)) { #show all text label
      text(Xpart.plot2[,1], Xpart.plot2[,2], labels=x$varnames[-c(1,2)], cex=labcex, pos=1)
    } else { #show selected label
      which.txtlab2 = which.txtlab + 2
      varnames2 = x$varnames
      varnames2[as.numeric(paste(- which.txtlab2))] = ""
      text(Xpart.plot2[,1], Xpart.plot2[,2], labels=varnames2[-c(1,2)], cex=labcex, pos=1)
    }
  }
  
  abline(h = 0)
  abline(v = 0)
  
  legend(0.8*max(Zcors), 0, legend = round(x$tau0,2), cex = labcex,
         yjust = 0.5, x.intersp = 0, y.intersp = 0,
         bg = ifelse(par("bg")== "transparent", "white", par("bg")), 
         box.lty = 0)
  
  box()
  
  contour(Zcors, Ycors, taus, levels = clevels, 
          add = T, labcex = labcex, ...)
  
  contour(Zcors, Ycors, taus, levels = 0, lwd = 2,
          add = T, col = col.zero,lty = lty.zero,labcex = labcex,...)
  
  contour(as.numeric(dimnames(x$sp.z)[[2]]), as.numeric(dimnames(x$sp.y)[[1]]), taus.est/se.est, labels = "N.S.",
          levels = c(-1, 1)*qnorm(signif.level/2), add = T, col = col.insig,
          lty = lty.insig, labcex = labcex, lwd = 2,...)
  
  if (all(sign(Xpart[,1])!=sign(x$tau0))){
     warning("Cannot add data line because XXXXX.")
  }else{
     if(data.line & length(Xpart)>1){
      proj.pts = apply(Xpart.plot^2, 1, mean)
      max.pt = Xpart.plot[proj.pts == max(proj.pts),]
      max.pt[1] = sign(x$tau0)*abs(max.pt[1])
      Zdiff = (Zcors-max.pt[1])*sign(max.pt[1])
      Zgrt = which(Zdiff < 0)
      zcor = Zgrt[which(Zdiff[Zgrt] == max(Zdiff[Zgrt]))]  
      if((Zcors[zcor] > max.pt[1] & zcor > 1)||(zcor==length(Zcors))){ 
        zpts = c(zcor-1, zcor)
      }else{
        zpts = c(zcor, zcor+1)
      }
      Ydiff = Ycors-max.pt[2]
      Ygrt = which(Ydiff < 0)
      ycor = Ygrt[which(Ydiff[Ygrt] == max(Ydiff[Ygrt]))]
      if((Ycors[ycor] > max.pt[2] & ycor > 1)||(ycor==length(Ycors))){ 
        ypts = c(ycor-1, ycor)
      }else{
        ypts = c(ycor, ycor+1)
      }
      clevel = ((Zcors[zpts[2]] - Zcors[zpts[1]])*(Ycors[ypts[2]] - Ycors[ypts[1]]))^(-1)*
          sum(taus[zpts, ypts]*
                matrix(c(-(Zcors[zpts[2]] - max.pt[1])*(Ycors[ypts[1]] - max.pt[2]), 
                         (Zcors[zpts[1]] - max.pt[1])*(Ycors[ypts[1]] - max.pt[2]),
                         (Zcors[zpts[2]] - max.pt[1])*(Ycors[ypts[2]] - max.pt[2]),
                         -(Zcors[zpts[1]] - max.pt[1])*(Ycors[ypts[2]] - max.pt[2])), 
                       nrow = 2, byrow = T))
      contour(Zcors, Ycors, taus, levels = round(clevel,2),
                add = T, col = "grey",labcex = labcex, lwd = 2,...)
      }else{
        if(data.line)
          warning("Cannot add data line because there are no non-treatment covariates.")
      }
    }
}


#############
#sensPlotCombo
#Plot for results from combined run
#############

sensPlotCombo = function(x, contour.levels, col.zero, lty.zero, col.insig, lty.insig, data.line, X.pch, signif.level, labcex, limit.Xplot, txtlab, which.txtlab,...){
  ##Add row/column for zeta = 0 if not included in grid
  null.tau=x$tau0
  null.se=x$se.tau0
  part.cors = x$sensParam == "cor"
  
  Zcors = x$sp.z #horizontal grids of U
  Ycors = x$sp.y #vertical grids of U
  
  if(part.cors){
    print.text <- "Partial correlations with U"
  }else {
    print.text <- "Coefficients on U" 
  }
  
  ##############  
  #######Only inclue X on the plot for GLM-style
  ############
  if(x$model.type == "GLM"){
    Xpart = x$Xcoef[!is.na(x$Xcoef[,1]) & !is.na(x$Xcoef[,2]),] #coefficients of null model.  
    Xpart.plot = x$Xcoef.plot[!is.na(x$Xcoef.plot[,1]) & !is.na(x$Xcoef.plot[,2]),]
    Xpart.plot2 = cbind(Xpart.plot[,1],Xpart.plot[,2], ifelse(Xpart[,2]>=0,1,0)) #MH: add sign of coef of X on Y to Xpart  
    
    
    #note that due to correlation among Xs, some may not appear on plot
    #because observed partial cors don't map directly to coefs in this case
    #forcing inclusion can lead to difficult to read plot.  
    if(is.null(X.pch)){
      X.pch = ifelse(Xpart[,2]>=0,3,6) # plus sign for non-transformed plots, reverse triangle for transformed plots
    }else{
      X.pch = ifelse(Xpart[,2]>=0,X.pch[1],X.pch[2])  
    }
    if (any(Xpart[,2]<0)) {
      cat("Note: Predictors with negative coefficients for the response surface have been transformed through multiplication by -1 and are displayed as inverted triangles.", "\n")    
    }
  }else{
    Xpart.plot2 <- Xpart.plot <- NULL
  }  
  
  taus = x$tau
  taus[Zcors==0] = null.tau
  taus[Ycors==0] = null.tau
  taus.est = tapply(taus, list(Zcors, Ycors), mean, na.rm = T)
  W = tapply(x$se.tau^2, list(Zcors, Ycors), mean, na.rm = T)
  K = length(W)/length(Ycors)
  B = tapply(x$tau, list(Zcors, Ycors), sd, na.rm = T)^2
  se.est = sqrt(W+(1+1/K)*B)
  taus = taus.est
  se.taus = se.est
  
  taus <<- taus
  se.taus <<- se.taus
  Zcors = as.numeric(dimnames(taus)[[1]])
  Ycors = as.numeric(dimnames(taus)[[2]])

  if(sum(is.na(taus)) > 0){
    empty.cells = which(is.na(taus))
    colno = empty.cells%/%dim(taus)[1]+1
    ylow = Ycors[colno-1]
    ymis = Ycors[colno]
    yhigh = Ycors[colno+1]
    rowno = empty.cells%%dim(taus)[1]
    zlow = Zcors[rowno-1]
    zmis = Zcors[rowno]
    zhigh = Zcors[rowno+1]
    for(i in 1:length(empty.cells)){
      if(rowno[i] %in% 1:length(Zcors)){
        zinterp = (zhigh[i]-zmis[i])/(zhigh[i]-zlow[i])*taus[rowno[i]-1, colno[i]]+(zmis[i]-zlow[i])/(zhigh[i]-zlow[i])*taus[rowno[i]+1, colno[i]]
        #zseinterp = (zhigh[i]-zmis[i])/(zhigh[i]-zlow[i])*se.taus[rowno[i]-1, colno[i]]+(zmis[i]-zlow[i])/(zhigh[i]-zlow[i])*se.taus[rowno[i]+1, colno[i]]
      }else{
        zinterp <- zseinterp <- NA
      }
      
      if(colno[i] %in% 1:length(Ycors)){
        yinterp = (yhigh[i]-ymis[i])/(yhigh[i]-ylow[i])*taus[rowno[i], colno[i]-1]+(ymis[i]-ylow[i])/(yhigh[i]-ylow[i])*taus[rowno[i], colno[i]+1]
        #yseinterp = (yhigh[i]-ymis[i])/(yhigh[i]-ylow[i])*se.taus[rowno[i], colno[i]-1]+(ymis[i]-ylow[i])/(yhigh[i]-ylow[i])*se.taus[rowno[i], colno[i]+1]
      }else{
        yinterp <- yseinterp <- NA
      }
      interp = mean(c(zinterp,yinterp), na.rm = T)
      taus[rowno[i],colno[i]] = interp
      
      #seinterp = mean(c(zseinterp,yseinterp), na.rm = T)
      #se.taus[rowno[i],colno[i]] = seinterp
    }
    warning("Standard errors not interpolated; N.S. line will be fragmented")
  }
  if(is.null(contour.levels)){
    exTau = c(taus[dim(taus)[1], dim(taus)[2]], taus[1,dim(taus)[2]]) #extreme values of tau at right end
    clevels = round(seq(exTau[2]*.8, exTau[1]*.8, length.out = 14), 2) #vals at which contours are drawn
  }else{
    clevels = contour.levels
  } 
  
  par(mgp = c(2,.5,0)) #dist of axis label, tick mark label, tick mark
  if(part.cors){
    xlab = expression(paste("Partial cor. with U in model for treatment, ", rho^zu))
    ylab = expression(paste("Partial cor. with U in model for response, ", rho^yu))
  }else{
    xlab = expression(paste("Coef. on U in model for treatment, ", zeta^z))
    ylab = expression(paste("Coef. on U in model for response, ", zeta^y))
  }
  
  if (limit.Xplot || x$model.type == "BART") {
    #old codes
    plot(Xpart.plot2[,1], Xpart.plot2[,2], col=c("red","blue")[as.factor(Xpart.plot2[,3])], xlim = c(min(Zcors, na.rm = T),max(Zcors, na.rm = T)), 
         ylim = c(min(Ycors, na.rm = T),max(Ycors, na.rm = T)), pch = X.pch, xlab = xlab, ylab = ylab)    
  } else {
    #MH: define max, min of plots
    xplot.min = ifelse(min(Zcors, na.rm = T)<min(Xpart.plot[,1]),min(Zcors, na.rm = T),min(Xpart.plot[,1]))
    xplot.max = ifelse(max(Zcors, na.rm = T)>max(Xpart.plot[,1]),max(Zcors, na.rm = T),max(Xpart.plot[,1]))
    yplot.max = ifelse(max(Ycors, na.rm = T)>max(Xpart.plot[,2]),max(Ycors, na.rm = T),max(Xpart.plot[,2]))  
    plot(Xpart.plot2[,1], Xpart.plot2[,2], col=c("red","blue")[as.factor(Xpart.plot2[,3])], xlim = c(xplot.min,xplot.max),
         ylim = c(0,yplot.max), pch = X.pch, xlab = xlab, ylab = ylab)
  }
  
  #codes for txtlab
  if (txtlab) {
    if (is.null(which.txtlab)) { #show all text label
      text(Xpart.plot2[,1], Xpart.plot2[,2], labels=x$varnames[-c(1,2)], cex=labcex, pos=1)
    } else { #show selected label
      which.txtlab2 = which.txtlab + 2
      varnames2 = x$varnames
      varnames2[as.numeric(paste(- which.txtlab2))] = ""
      text(Xpart.plot2[,1], Xpart.plot2[,2], labels=varnames2[-c(1,2)], cex=labcex, pos=1)
    }
  }
  
  abline(h = 0)
  abline(v = 0)
  
  legend(0.8*max(Zcors), 0, legend = round(x$tau0,2), cex = labcex,
         yjust = 0.5, x.intersp = 0, y.intersp = 0,
         bg = ifelse(par("bg")== "transparent", "white", par("bg")), 
         box.lty = 0)
  
  box()
  
  contour(Zcors, Ycors, taus, levels = clevels, 
          add = T, labcex = labcex, ...)
  
  contour(Zcors, Ycors, taus, levels = 0, lwd = 2,
          add = T, col = col.zero,lty = lty.zero,labcex = labcex,...)
  
  contour(Zcors, Ycors, taus/se.taus, labels = "N.S.",
          levels = c(-1, 1)*qnorm(signif.level/2), add = T, col = col.insig,
          lty = lty.insig, labcex = labcex, lwd = 2,...)
  
  if(x$model.type == "GLM"){
    if (all(sign(Xpart[,1])!=sign(x$tau0))){
      warning("Cannot add data line because XXXXX.")
    }else{
      if(data.line & length(Xpart)>1){
        proj.pts = apply(Xpart.plot, 1, mean)
        max.pt = Xpart.plot[proj.pts == max(proj.pts[sign(Xpart.plot[,1])==sign(x$tau0)]),]
        zcor = (1:length(Zcors))[abs(Zcors-max.pt[1]) ==  min(abs(Zcors-max.pt[1]))]
        if((Zcors[zcor] > max.pt[1] & zcor > 1)||(zcor==length(Zcors))){ 
          zpts = c(zcor-1, zcor)
        }else{
          zpts = c(zcor, zcor+1)
        }
        ycor = (1:length(Ycors))[abs(Ycors-max.pt[2]) ==  min(abs(Ycors-max.pt[2]))]
        if((Ycors[ycor] > max.pt[2] & ycor > 1)||(ycor==length(Ycors))){ 
          ypts = c(ycor-1, ycor)
        }else{
          ypts = c(ycor, ycor+1)
        }
        clevel = ((Zcors[zpts[2]] - Zcors[zpts[1]])*(Ycors[ypts[2]] - Ycors[ypts[1]]))^(-1)*
          sum(taus[zpts, ypts]*
                matrix(c(-(Zcors[zpts[2]] - max.pt[1])*(Ycors[ypts[1]] - max.pt[2]), 
                         (Zcors[zpts[1]] - max.pt[1])*(Ycors[ypts[1]] - max.pt[2]),
                         (Zcors[zpts[2]] - max.pt[1])*(Ycors[ypts[2]] - max.pt[2]),
                         -(Zcors[zpts[1]] - max.pt[1])*(Ycors[ypts[2]] - max.pt[2])), 
                       nrow = 2, byrow = T))
        contour(Zcors, Ycors, taus, levels = round(clevel,2),
                add = T, col = "grey",labcex = labcex, lwd = 2,...)
      }else{
        if(data.line)
          warning("Cannot add data line because there are no non-treatment covariates.")
      }
    }
  }else{
    if(data.line)
      warning("Cannot add data line for BART sensitivity analysis.")
  }
}

