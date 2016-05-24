# design matrix for subject covs
covdesign<-function(formel,covs,ENV)
{
   # number of factor levels for terms from formula
   vars<-attr(terms(formel),"term.labels")
   order<-attr(terms(formel),"order")
            #vars<-vars[order==1]
   vars<-sort(vars[order==1]) # maybe not necessary to sort names
   levs<-apply(as.matrix(covs[,vars]),2,max) # number of factor levels


   maineffects<-as.data.frame(gfac2(levs)) #
   maineffects<-data.frame(apply(maineffects,2,factor))
   names(maineffects)<-vars
   ENV$maineffects<-maineffects            # covs design matrix with factors


   form<-formula(formel)
###   form.terms<-terms(form, keep.order=TRUE, simplify = TRUE)
   form.terms<-terms(form, simplify = TRUE)
   mmat<-model.matrix(form.terms,data=maineffects)  # intercept included
   as.matrix(unique(mmat))
}
