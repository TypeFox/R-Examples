# split data according to subject covariates
splitCovs<-function(dat,covs,formel,elim,ENV)
{
   sort.model.terms<-function(newForm){
       sort(sapply(
          strsplit(colnames(attr(terms(newForm), "factors")), ":", fixed = TRUE),
          function(x) paste(sort(x), collapse = ":")))
   }

   if(is.null(covs)){            # no covariates
      formel<-~1                 # reset formulas
      elim<-~1
   }

   if(elim=="~1") elim<-formel   # if no eliminate specification use model specification

   if(elim=="~1"){               # subject covariates
      covdesmat<-matrix(1,1)
      elimdesmat<-covdesmat
      covs<-as.matrix(rep(1,nrow(dat)),ncol=1)
   } else {                      # there are terms to eliminate
      if (formel=="~1"){         #   no covariates in actual model
           covdesmat<-matrix(1,1)
           colnames(covdesmat)<-""
      } else {                   # there are terms in actual model

           # sort formel terms first before setting up covdesmat (18.1.10)
           sorted.frml.terms<-sort.model.terms(formel)
           frmlstr<-paste("~",paste(sorted.frml.terms,collapse="+"))
           formel<-formula(frmlstr)

           covdesmat<-covdesign(formel,covs,ENV)
           ENV$model.covs<-ENV$maineffects
      }
           # sort formel terms first before setting up elimdesmat (18.1.10)
           sorted.elim.terms<-sort.model.terms(elim)
           elimstr<-paste("~",paste(sorted.elim.terms,collapse="+"))
           elim<-formula(elimstr)
           elimdesmat<-covdesign(elim,covs,ENV)
   }


#   all.term.labels<-attr(terms(formula(elimstr)),"term.labels")
   all.term.labels<-attr(terms(elim),"term.labels")
   order<-attr(terms(elim),"order")

#   all.frml.term.labels<-attr(terms(formula(frmlstr)),"term.labels")
   all.frml.term.labels<-attr(terms(formel),"term.labels")
   order.frml<-attr(terms(formel),"order")

   ENV$covlevels<-NULL

   if (length(order)>0){                       ## there are terms to eliminate


       maineffect.terms<-all.term.labels[order==1]
       if( any(table(covs[,maineffect.terms]) == 0) )
            stop("Crossclassification of subject covariates yields zero cell(s).")
       maineffect.terms.frml<-all.frml.term.labels[order.frml==1]
       if(length(maineffect.terms.frml)>0){                                        # only if cov terms in formel
         ENV$covlevels<-apply(as.matrix(covs[,maineffect.terms.frml]),2,max)
         ENV$elimcovlevels<-apply(as.matrix(covs[,maineffect.terms]),2,max)        # 23.11.09
         names(ENV$elimcovlevels)<-maineffect.terms                                # subj covs have to based on elim terms
         names(ENV$covlevels)<-maineffect.terms.frml
       }

       # split data according to subj cov groups in eliminate formula
       cList<-split(dat,covs[,maineffect.terms])

       # remove terms in elim design matrix not in formel design matrix
                                                                                 # again obsolete 18.1.10
                                                                                 # now obsolete: 23.11.09
#      enam<-colnames(elimdesmat)                                                # had to be introduced
#      cnam<-colnames(covdesmat)                                                 # since terms in
#      nam.elim<-sapply(1:ncol(elimdesmat),                                      # interactions
#                        function(i){txt<-sort(unlist(strsplit(enam[i],":")));   # may have different
#                                    paste(txt,collapse=":")}                    # order in elim and formel
#                      )                                                         # correct in R 2.6.x and
#      nam.cov<-sapply(1:ncol(covdesmat),                                        # maybe in 2.7.0
#                        function(i){txt<-sort(unlist(strsplit(cnam[i],":")));   #
#                                    paste(txt,collapse=":")}                    # now obsolete
#                      )                                                         #
#                                                                                #
#       model.terms<-nam.elim %in% nam.cov                                       #

       model.terms<-colnames(elimdesmat) %in% colnames(covdesmat)                # had been replaced


       for (i in 1:length(cList))
           cList[[i]]<-c(cList[i],list(cov=elimdesmat[i,model.terms]))
       ENV$covdesmat  <-elimdesmat[,model.terms]
       if (formel != "~1")
           colnames(ENV$model.covs)<-maineffect.terms.frml   # order of elim terms for factor matrix to
                                                        # properly label worth matrix

   } else {                                    ## nothing to eliminate - basic model

       cList<-split(dat,as.factor(rep(1,nrow(dat))))
       cList[[1]]<-c(cList[1],list(cov=1))
       ENV$covdesmat<-covdesmat
   }

   ENV$ncovpar<-ncol(covdesmat)       # number of covariates
   ENV$formel<-formel
   ENV$elim<-elim
   cList
}
