`WGassociation` <-
function (formula, data, model=c("all"), quantitative = is.quantitative(formula, 
    data), genotypingRate = 80, level=0.95, ...)
 {
  
    if(!inherits(data,"setupSNP"))
     stop("data must be an object of class 'setupSNP'")

    if (length(attr(data,"colSNPs")) > 2000 & (length(model) > 1 | any(model%in%"all"))) 
        stop("Select only one genetic model when more than 2000 SNPs are analyzed \n or use 'scanWGassociation' function")
 
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m0)]

#
# aceptar respuesta sin formula
#    
	if( length(grep("~",mf[[2]]))==0){
		formula<-as.formula(paste(mf[[2]],"~1",sep=""))
		formula.1<- list(formula)
		mode(formula.1)<-"call"
		mf[2]<-formula.1
	}

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    temp0 <- as.character(mt)
    adj <- paste(temp0[2], temp0[1], temp0[3])

    Terms <- if (missing(data)) 
        terms(formula)
    else terms(formula, data = data)
    ord <- attr(Terms, "order")
    if (any(ord > 1)) 
     stop("interaction term is not implemented")


    association.i<-function(snp.i,adj,data,model,quantitative,genotypingRate,level, ...)
     {
       association( as.formula(paste(adj,"+",snp.i)) , data=data,
          model=model, quantitative=quantitative, genotypingRate=
          genotypingRate, level=level, ...)
     }
 

    colSNPs<-attr(data,"colSNPs")
    if (! (is.vector(colSNPs) & length(colSNPs) > 0)) stop("data should have an attribute called 'colSNPs'. Try again 'setupSNP' function")

#    if (is.vector(colSNPs) & length(colSNPs) > 1) 
#        dataSNPs <- data[, colSNPs]
#    else stop("data should have an attribute called 'colSNPs'. Try again 'setupSNP' function")


    
    type<-charmatch(model,c("codominant","dominant","recessive","overdominant","log-additive","all"))
    type<-sort(type)

    if (any(type%in%6))
      type<-1:5
  
    if(any(is.na(type)))
     stop("model must be 'codominant','dominant','recessive','overdominant', 
                     'log-additive', 'all' or any combination of them")

    SNPs<-attr(data,"label.SNPs")

    tab<-mclapply(SNPs, association.i, adj=adj, data=data, model=model,
            quantitative=quantitative, genotypingRate=
            genotypingRate, level=level, ...) 

    names(tab)<-SNPs
    attr(tab,"label.SNPs")<-attr(data,"label.SNPs")
    attr(tab,"models")<-type
    attr(tab,"quantitative")<-quantitative
    out<-extractPval(tab)

    attr(out,"tables")<-tab
    attr(out,"label.SNPs")<-attr(data,"label.SNPs")
    attr(out,"models")<-type
    attr(out,"quantitative")<-quantitative
    attr(out,"pvalues")<-out
    attr(out,"gen.info")<-attr(data,"gen.info")
    attr(out,"whole")<-attr(data,"whole")
    attr(out,"colSNPs")<-attr(data,"colSNPs")

    class(out)<-c("WGassociation","data.frame")

    out
}

