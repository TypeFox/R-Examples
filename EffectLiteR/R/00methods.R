
setMethod("show", "effectlite", function(object) {
  
  ng <- object@input@ng
  nk <- object@input@nk
  nz <- object@input@nz
  vnames <- object@input@vnames    
  vlevels <- object@input@vlevels
  gammas <- object@parnames@gammas
  gammalabels <- object@parnames@gammalabels
  
  label.g.function <- "(K,Z)"; label.covs <- ",K,Z"
  if(nk==1 & nz==0){label.g.function <- "()"; label.covs <- ""}
  if(nk>1 & nz==0){label.g.function <- "(K)"; label.covs <- ",K"}
  if(nk==1 & nz>0){label.g.function <- "(Z)"; label.covs <- ",Z"}
  
  cat("\n\n--------------------- Variables  --------------------- \n\n")
  cat("Outcome variable Y: ", paste0(vnames$y), "\n")
  cat("Treatment variable X: ", paste0(vnames$x), "  (Reference group: ", 
      paste0(object@input@control, ")\n"))
  if(!is.null(vnames$k)){
    cat("Categorical covariates K: ", paste0(vnames$k), "\n")
  }
  if(!is.null(vnames$z)){
    cat("Continuous covariates Z: ", paste0(vnames$z), "\n")
  }
  if(!is.null(vnames$propscore)){
    v <- vnames$propscore
    if(is(v, "formula")){v <- all.vars(v[[3]])}  
    cat("Covariates for propensity score V: ", paste0(v), "\n")
  }

  
  if(nk>1){
    cat("\nLevels of Unfolded Categorical Covariate K \n")
    tmp <- vlevels$levels.k.original
    tmp <- tmp[length(tmp):1]
    tmp <- expand.grid(tmp)
    tmp$K <- vlevels$kstar
    tmp <- tmp[,ncol(tmp):1]
    print(tmp, row.names=F, print.gap=3)
  }
  
  
  cat("\n\n --------------------- Regression Model --------------------- \n")
  
  tmp <- paste0("E(Y|X",label.covs,") = ")
  tmp <- paste0(tmp, "g0",label.g.function," + ")
  tmp <- paste0(tmp, paste0("g",1:(ng-1),label.g.function,"*I_X=",1:(ng-1), 
                            collapse=" + "))
  cat("\n",tmp, "\n")
  
  gammalabels2 <- gammalabels[,,1]
  gammalabels2[1] <- ""
  
  for(i in 1:ng){
    tmp <- paste0("  g",i-1,label.g.function," = ")
    tmp <- paste0(tmp, paste(gammas[,,i], gammalabels2, sep=" * ", collapse=" + "))
    tmp <- gsub("*  ", "", tmp, fixed=TRUE)
    if(length(gammalabels2)==1){tmp <- gsub("*", "", tmp, fixed=TRUE)}
    
    if(nchar(tmp) > 80){
      ## split g function over several lines
      tmp <- unlist(strsplit(tmp, " + ", fixed=TRUE))
      tmp <- capture.output(cat(tmp, sep=" + ", fill=80))
      tmp[2:length(tmp)] <- paste0("            + ",tmp[2:length(tmp)])
      cat(tmp, sep="\n")
    } else{
      cat(tmp, "\n")
    }
  }
  
  ## print coefficients of g-Functions  
  for(i in 1:ng){
    if(i==1){
      tmp <- paste0("Intercept Function g",i-1,label.g.function)
      cat("\n",tmp, "\n\n")
    }else{
      tmp <- paste0("Effect Function g",i-1,label.g.function)
      tmp <- paste0(tmp, "   [", object@input@vnames$x, 
                    ": ", object@input@vlevels$levels.x.original[i],
                    " vs. ",object@input@vlevels$levels.x.original[1], "]")
      cat("\n",tmp, "\n\n")
    }
    tmp <- object@results@gx[[i]]
    tmp[,2:5] <- round(tmp[,2:5], digits=3)
    print(tmp, print.gap=3, row.names=FALSE)
  }
  
  
  cat("\n\n--------------------- Cell Counts  --------------------- \n\n")

  
  if(nk>1){
    cat("\nCells \n")
    tmp <- expand.grid(K=vlevels$kstar, X=vlevels$levels.x.original)[,2:1]
    tmp$Cell <- vlevels$cell
    print(tmp, print.gap=3)
    
  }
  
  if(nk==1){
    cat("Cells \n")
    tmp <- data.frame(X=vlevels$levels.x.original)
    print(tmp, row.names=F, print.gap=3)
    
  }
  
  cat("\n")
  cat("Cell Counts \n\n")
  cat("This table shows cell counts including missings. \n")
  cat("See also output under lavaan results for number of observations \n")
  cat("actually used in the analysis. \n\n")
  
  if(nk==1){
    print(ftable(object@input@data[vnames$x]), print.gap=3)
  }else{
    cellcounts <- as.formula(paste0(paste(vnames$k, collapse="+"), 
                                    "~", vnames$x))
    print(ftable(cellcounts, data=object@input@data), print.gap=3)
  }
  
  
  cat("\n\n--------------------- Main Hypotheses --------------------- \n\n")
  if(nrow(object@results@hypotheses)==0){
    cat("Wald tests for main hypotheses are currently not available for models with \n non-standard SEs and for models with (in-)equality constraints (e.g., on interactions).")
  }else{
    hypotheses <- object@results@hypotheses
    names(hypotheses) <- c("Wald Chi-Square", "df", "p-value")
    print(hypotheses, digits=3, print.gap=3)
  }
  
  cat("\n\n --------------------- Adjusted Means --------------------- \n\n")
  namesadjmeans <- paste0("Adj.Mean",0:(ng-1))
  adjmeans <- object@results@adjmeans
  row.names(adjmeans) <- namesadjmeans
  print(adjmeans, digits=3, print.gap=3)
  
  
  cat("\n\n --------------------- Average Effects --------------------- \n\n")
  namesEgx <- paste0("E[g",1:(ng-1),label.g.function,"]")
  Egx <- object@results@Egx
  row.names(Egx) <- namesEgx
  print(Egx, digits=3, print.gap=3)
  
  
  if(!(nz==0 & nk==1)){
    cat("\n\n --------------------- Effects given a Treatment Condition --------------------- \n\n")
    tmp <- expand.grid(g=1:(ng-1), x=0:(ng-1))
    namesEgxgx <- paste0("E[g",tmp$g,label.g.function,"|X=",tmp$x, "]")
    Egxgx <- object@results@Egxgx
    row.names(Egxgx) <- namesEgxgx
    print(Egxgx, digits=3, print.gap=3)
    
  }
  
  if(nk>1){
    cat("\n\n --------------------- Effects given K=k --------------------- \n\n")
    tmp <- expand.grid(g=1:(ng-1), k=0:(nk-1))
    namesEgxgk <- paste0("E[g",tmp$g,label.g.function,"|K=",tmp$k,"]")
    Egxgk <- object@results@Egxgk
    row.names(Egxgk) <- namesEgxgk
    print(Egxgk, digits=3, print.gap=3)    
  }
  
  if(nk>1 & nz>0){
    cat("\n\n --------------------- Effects given X=x, K=k --------------------- \n\n")
    Egxgxk <- paste0("Eg",tmp$g,"gx",tmp$x,"k",tmp$k)    
    tmp <- expand.grid(g=1:(ng-1), x=0:(ng-1), k=0:(nk-1))
    namesEgxgxk <- paste0("E[g",tmp$g,label.g.function,"|X=",tmp$x,", K=",tmp$k,"]")
    Egxgxk <- object@results@Egxgxk
    row.names(Egxgxk) <- namesEgxgxk
    print(Egxgxk, digits=3, print.gap=3)    
  }
  
  
  propscore <- object@input@vnames$propscore
  if(!is.null(propscore)){
    cat("\n\n --------------------- Propensity Score Model --------------------- \n\n")
    cat("Model equation: log[P(X=1|V)/P(X=0|V)] = h1(V)", "\n")
    if(ng > 2){
      for(i in 3:ng){
        tmp <- paste0("                ", 
                      "log[P(X=", (i-1), "|V)/P(X=0|V)] = h",(i-1),"(V)",
                      "\n")
        cat(tmp)
      }
    }
    cat("R formula for nnet::multinom: ", object@input@outprop$formula, "\n")
    cat("\nEstimate\n")
    print(object@input@outprop$coef, digits=3, print.gap=3)
    cat("\nStandard Error\n")
    print(object@input@outprop$se, digits=3, print.gap=3)
    cat("\nEst./SE\n")
    print(object@input@outprop$tval, digits=3, print.gap=3)
  }
  
  
})

