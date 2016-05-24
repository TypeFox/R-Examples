predictQR <-
function(object, newdata, xreg){
#xreg o newdata have to be provided. Nessun controllo (neanche sull'ordine dei coef)
#if xreg is provided, newdata is ignored
#
#Attenzione: puo' dare problemi se la formula contiene "factor" 
#model.matrix(as.formula(m1$call$formula), data=newdata)
bspline <- function(x, ndx, xlr = NULL, knots=NULL, deg = 3, deriv = 0, outer.ok=FALSE) {
    # x: vettore di dati
    # xlr: il vettore di c(xl,xr)
    # ndx: n.intervalli in cui dividere il range
    # deg: il grado della spline
    # Restituisce ndx+deg basis functions per ndx-1 inner nodi
    #ci sono "ndx+1" nodi interni + "2*deg" nodi esterni
#    require(splines)
  if(is.null(knots)) {
    if (is.null(xlr)) {
        xl <- min(x) - 0.01 * diff(range(x))
        xr <- max(x) + 0.01 * diff(range(x))
    }
    else {
        if (length(xlr) != 2)
            stop("quando fornito, xlr deve avere due componenti")
        xl <- xlr[1]
        xr <- xlr[2]
    }
    dx <- (xr - xl)/ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
      }
      #else {
      #if(length(knots)!=(ndx+1+2*deg)) stop("errore nel numero di nodi fornito")
      #}
    B <- splineDesign(knots, x, ord = deg + 1, derivs = rep(deriv, length(x)), outer.ok=outer.ok)
    B
    }
    b<-as.matrix(object$coefficients)  #corretto in 0.3-1
    if(missing(xreg)){
      if(missing(newdata)) stop("please, provide `newdata' or 'xreg'")
      nomiCoef<-rownames(b)
      info.smooth<-object$info.smooth #estrai 
      p<- info.smooth$deg + info.smooth$ndx
      nome.smooth<-names(object$BB)
      nomiCoefPen<-paste(nome.smooth, ".ps." ,1:p,sep="")
      id.coef.smooth<-match(nomiCoefPen, nomiCoef)
      b.smooth<-b[id.coef.smooth, ]
      nomiCoefUnpen<-nomiCoef[-id.coef.smooth]
      nomiVarModello<-all.vars(object$call$formula)[-1]
      id.var<-match(nomiVarModello, names(newdata))
      if(any(is.na(id.var))) stop("`newdata' does not include all the covariates in the model")
      newdata<-newdata[,id.var,drop=FALSE] 
      x.new<-newdata[,match(nome.smooth,names(newdata))]
      newdata<-newdata[,-match(nome.smooth,names(newdata)),drop=FALSE]
      #costruire la base su x.new.. prendere il min max..
      m<-min(attr(object$BB[[1]], "covariate.35"))       #as.numeric(attr(object$BB[[1]], "covariate.n"))
      M<-max(attr(object$BB[[1]], "covariate.35"))
      B.new<-bspline(c(m, x.new, M), ndx=info.smooth$ndx, deg=info.smooth$deg)
      B.new<-B.new[-c(1,nrow(B.new)),] #???
      XREG<-as.matrix(cbind(newdata,B.new))
      colnames(XREG)<-c(nomiCoefUnpen, nomiCoefPen)
      if("(Intercept)" %in% nomiCoef ) XREG<-cbind("(Intercept)"=1,XREG)
      fit<-drop(XREG[,nomiCoef]%*%b)
      } else {
      if(!missing(newdata)) warning("`newdata' ignored when 'xreg' is provided")      
      fit<-drop(xreg%*%b)
      }
      return(fit)
      }
