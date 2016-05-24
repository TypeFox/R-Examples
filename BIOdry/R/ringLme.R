ringLme <- structure(function# LME modeling
### LME modeling of multilevel data frames. 
##details<<This function implements \code{\link{lme}} to fit linear mixed-effects models on data frames with factor-level columns. Two kind of model formulas are available: 'lmeForm' and 'tdForm', these characters implement functions with same names (\code{\link{tdForm}} and \code{\link{lmeForm}}). Users can develop their own LME formulas by modifying arguments in any of these methods. After LME models are fitted, they  can be extended by modeling heteroscedasticity and autocorrelation of the residuals. Nevertheless, such residual modeling would take long time depending on the complexity  of the modeled levels.
##references<< Pinheiro J. C., D. M. Bates. 2000. Mixed-effects models in S and S-PLUS. Springer, New York.                                          
(
    rd, ##<<\code{data.frame} object containing factor-level columns. 
    form = 'lmeForm', ##<<\code{character}. Method name of LME
                      ##formula. Available methods are 'lmeForm' and
                      ##'tdForm' (see details).
    var.mod = FALSE, ##<< \code{logical}. If TRUE then the fitted
                     ##model is extended by modeling residual
                     ##heteroscedaticity with 'varConstPower' class in
                     ##\code{\link{nlme}}.
    arma.mod = FALSE, ##<< \code{logical}. If TRUE then the fitted
                      ##model is extended by modeling residual
                      ##autocorrelation with class 'corARMA' in
                      ##\code{\link{nlme}}, with p = 1 and q = 1,
                      ##which models residual autocorrelation for lags
                      ##> 2.
    res.data =TRUE, ##<< \code{logical}. Save residuals as a
                    ##multilevel data frame. If TRUE then a data frame
                    ##of name 'resid' is added to output list
    ... ##<< Further arguments to be passed to function in \code{form}.
) {

    pr.cov <- function(form){
        chf <- Reduce(paste,deparse(form))
        fnc <- gsub('.*~|\\|.*','',chf)
        fnc <- paste('~',fnc,sep = '')
        return(formula(fnc))}
    
    ## Implementation of lme form:
    if(grepl('~',form))
        formu <- formula(form)
    if(!grepl('~',form)){
        arf <- arguSelect(rd,fun = form,...)
        formu <- do.call(form,arf)}
        prc <- pr.cov(formu)
        environment(prc) <- .GlobalEnv
    gd <- groupedData(
        formu,data = na.omit(rd))  #<<
    mem <- do.call(
        lme,list(gd,random = pdDiag(prc),
            control = list(msMaxIter = 200)))
    if(var.mod)
        mem <- update(
        mem,weights = varConstPower(form = prc))
    if(arma.mod)
        mem <- update(mem,corr = corARMA(p = 1, q = 1))
    
    rset <- function(r.model){
        md <- r.model[['data']]
        tim <- colclass(md)[['tmp']]
        lev <- colclass(md)[['fac']]
        lg <- ncol(data.frame(getGroups(md)))
        dcum <- residuals(r.model,level = lg:1,type = 'p')
        dcum <- as.data.frame(dcum)
        nam. <- names(dcum)
        names(dcum) <- paste(nam.,'.res',sep = '')
        dres <- cbind(dcum,md[,c(tim,lev)])
        return(dres)}
    
    if(res.data)
        mem <- list(resid = rset(mem),model = mem,call = sys.call())
    
    return(mem)
### Depending on \code{res.data}, either an LME model, or a
### \code{list} with both: the LME model, and a residual data frame
### containing initial factor-level columns.
} , ex=function() {
    
    ##Multilevel data frame of tree-ring widths:
    data(Prings05,envir = environment())
    ## Radial increments measured on 2003:
    data(Pradii03,envir = environment())    
    ## Monthly precipitation sums and average temperatures:
    data(PTclim05,envir = environment())

    ##Modeling tree growth:
    mpin <- modelFrame(Prings05, y = Pradii03,
            form = NULL,# on.time = TRUE,
            MoreArgs = list(only.dup = TRUE,
            mp = c(1,1),un = c('mm','cm'),z = 2003))

    ## Detrending tree growth with a td-form model:
    rlme <- ringLme(mpin,form = 'tdForm')
    summary(rlme$model)

    ##a plot of the modeled fluctuations
    d <- groupedData(lmeForm(rlme$resid,lev.rm = 1),data = rlme$resid)
    plot(d,groups = ~ sample,auto.key = TRUE)
    
    ## A model of aridity: 
    cf <- modelFrame(PTclim05,
         lv = list('year','year'),
         fn = list('moveYr','wlai'),
         form = NULL)
    summary(cf)

    ## An lme model of aridity at 'plot' level:
    rmod <- ringLme(cf,form = 'lmeForm')
    summary(rmod$model)
    
    rk <- groupedData(lmeForm(rmod$resid),data=rmod$resid)
    plot(rk,ylab = 'detrended AI')
    
    
})
