svyhist<-function(formula, design, breaks = "Sturges", 
                  include.lowest = TRUE, right = TRUE, xlab=NULL,
                  main=NULL, probability=TRUE,
                  freq=!probability,...){
    if (inherits(design,"DBIsvydesign") || inherits(design,"ODBCsvydesign")){
      design$variables<-getvars(formula, design$db$connection, design$db$tablename, 
        updates = design$updates)
      class(design)<-"survey.design2"
    } 
    mf<-model.frame(formula,model.frame(design), na.action=na.pass)
    if (ncol(mf)>1) stop("Only one variable allowed.")
    variable<-mf[,1]
    varname<-names(mf)
    h <- hist(variable,  plot=FALSE, breaks=breaks,right=right)
    props <- coef(svymean(~cut(variable, h$breaks,right=right, include.lowest=include.lowest),
                          design, na.rm=TRUE))
    h$density<-props/diff(h$breaks)
    h$counts <- props*sum(weights(design,"sampling"))
    if (is.null(xlab)) xlab<-varname
    if (is.null(main)) main<-paste("Histogram of",varname)
    plot(h, ..., freq=freq,xlab=xlab,main=main)
}
