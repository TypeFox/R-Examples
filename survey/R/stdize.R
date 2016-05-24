##
##  This is how NCHS does it: postStratify to a table where proportions for by= are specified and then are applied within each cell of over=
##
svystandardize<-function(design, by, over, population, excluding.missing=NULL){
    
    if (!is.null(excluding.missing)){
         mf<-model.frame(excluding.missing, model.frame(design),na.action=na.omit)
	 naa<-attr(mf,"na.action") 
         if(!is.null(naa)) design<-design[-naa,]
    } 

    if(is.data.frame(population)) population<-population$Freq
    
    freemargins<-as.data.frame(svytable(over, design))
    fixedmargins<-as.data.frame(svytable(by,design))
    fixedmargins$Freq<-as.vector(population)/sum(as.vector(population))
    combined<-make.formula(c(attr(terms(by),"term.labels"), attr(terms(over),"term.labels")))
    allmargins<-as.data.frame(svytable(combined,design))
    allmargins$Freq<-as.vector(outer(fixedmargins$Freq, freemargins$Freq))

    design<-postStratify(design, combined, allmargins,partial=TRUE)
    design$call<-sys.call()
    design
}
