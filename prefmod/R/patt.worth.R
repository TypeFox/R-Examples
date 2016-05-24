# working version  - pw4.R -> tidied version pw5.R

patt.worth<-function(fitobj, obj.names=NULL, outmat="worth"){

  ## which model
  pattPCfitmodel <- FALSE
  pattdesignmodel <- FALSE
  pattNPMLmodel <- FALSE
  nclass<-1

  if("pattMod" %in% class(fitobj)) {
       pattPCfitmodel <- TRUE
  } else if(!is.null(fitobj$call$design)) {
     pattNPMLmodel <- TRUE
     dterm <- deparse(fitobj$call$design)            ## unten nocheinmal
     design <- tryCatch(get(dterm),error = function(e) FALSE)
     if("pattNPML" %in% class(fitobj)){
        nclass<-fitobj$call$k
        pattNPMLmodel <- TRUE
        pattdesignmodel<-TRUE
     }
  } else if(!is.null(fitobj$call$data)) {
     dterm <- deparse(fitobj$call$data)            ## unten nocheinmal
     design <- tryCatch(get(dterm),error = function(e) FALSE)
     if("pattdes" %in% class(design)) pattdesignmodel <- TRUE
  }

  if(!any(pattPCfitmodel, pattdesignmodel, pattNPMLmodel))
          stop(paste("Model result must be either from patt*.fit or having used patt.design"))




###  if("pattNPML" %in% class(fitobj)){
###     nclass<-fitobj$call$k
###     pattNPMLmodel <- TRUE
###     pattdesignmodel<-TRUE
###     design<-fitobj$data
###     if("pattdes" %in% class(design)) pattdesignmodel <- TRUE
###
###  } else {
###
###     pattdesignmodel <- FALSE
###     pattPCfitmodel <- ("pattMod" %in% class(fitobj))
###     if (!pattPCfitmodel){
###        design <- fitobj$data            ## unten nocheinmal
###        if("pattdes" %in% class(design)) pattdesignmodel <- TRUE
###     }
###
###     if (pattdesignmodel == pattPCfitmodel)  # both FALSE
###          stop(paste("Model result must be either from patt*.fit or having used patt.design"))
###  }

  ### setup array that later contains estimates for all objects x other variables
  if (pattdesignmodel) {

     ## subject covariates
     #num.scovs<-attr(design,"num.scovs")
     #if(!is.null(num.scovs))
     #    warning("Numerical subject covariates not (yet) implemented, they are ignored.\n
     #    Result is a matrix of intercepts for regression on numerical subject covariates!!")

     # extract object names from design frame attribute
     objnames<-attr(design,"objnames")
     objects<-objnames

     objcovnames<-unlist(dimnames(attr(design,"objcovs"))[2])
     objcovs<-attr(design,"objcovs")
     objnames<-c(objnames,objcovnames)

     if("gnm" %in% class(fitobj)){
       fterm<-as.character(attr(fitobj$terms,"predictor"))
     } else {
       fterm<-attr(fitobj$terms,"term.labels")
     }

     if (pattNPMLmodel){
        rterm <- deparse(fitobj$call$random)
        fterm <- deparse(fitobj$call$formula)
        fterm <- paste(fterm, rterm)
     }
     fterms<-unique(unlist(strsplit(fterm,"[ ()+*:~`]")))
     objnames<-intersect(objnames,fterms)

     objcovmodel<-any(objcovnames %in% fterms)
     used.objcovs<-intersect(objcovnames, fterms)

     if (objcovmodel) {
     objnum<-seq_along(objnames)
     names(objnum)<-objnames
     } else {
     objnum<-seq_along(objects)
     names(objnum)<-objects
     }

     # check if numeric subject covariates have been fitted
     num.scovs<-attr(design,"num.scovs")
     if(!is.null(num.scovs)){
       num.fterms<-intersect(num.scovs, fterms)
       if (length(num.fterms)>0){
         txt<-paste(
"Numerical subject covariates not (yet) implemented, they are ignored.",
"\n  Result is a matrix of intercepts for regression on numerical subject covariates!!")
         warning(txt)
       }
     }

     # extract subject covariates from design frame attribute
     cat.scovs<-attr(design,"cat.scovs")
     if(!is.null(cat.scovs))
       eterms <- cat.scovs
     else
       eterms <- vector(,0)

     # subject variables in design matrix
     if(length(eterms>0)) {
         # extract subject covariates from formula
         fterm<-deparse(fitobj$call$formula)
         fterms<-unique(unlist(strsplit(fterm,"[ ()+*:~]")))
         subjcov.names<-intersect(eterms,fterms)

         # achtung laenge 0
         if(length(subjcov.names)>0) { #subject variable in design matrix and in model
           # set up subj cov design
           if(pattNPMLmodel)
              levlist<-lapply(fitobj$data[subjcov.names],levels)
           else
              levlist<-lapply(fitobj$model[subjcov.names],levels)
           maxlev<-sapply(levlist,length)
         }
         if(nclass>1) {                                 #  add "MASS" for pattNPML models
            if(exists("levlist"))
              levlist$MASS<-paste(1:nclass)
            else
              levlist<-list(MASS=paste(1:nclass))
            maxlev<-sapply(levlist,length)
         }
           if (!exists("maxlev")){                      #subject variables in design matrix but none in model
                maxlev <- 1
                length(eterms)<-0
           }

     } else if(nclass>1) {                              # only classes no subject covariates
                levlist<-list(MASS=paste(1:nclass))
                maxlev<-sapply(levlist,length)

     # no subject variables in design matrix
     } else
          maxlev<-1

     # set up array
     array.dim<-c(OBJ=length(objects),maxlev)
     array.dimnames<-lapply(array.dim, function(i) 1:i)
     a<-array(0,dim=array.dim,dimnames=array.dimnames)
     dimnames(a)$OBJ<-objects
##}

  ### extract relevant estimates
     est<-coef(fitobj)
     estnamlist<-lapply(names(est),strsplit,":")
     obj.in.est<-sapply(estnamlist,function(x) any(objnames %in% unlist(x)))
     est<-est[obj.in.est]  # only model terms with obj
     estnamlist<-estnamlist[obj.in.est]

     # extract undecided and interactions covariate names from design frame attribute
     # if user has specified interaction with objects
     undec<-attr(design,"undec")
     ia<-attr(design,"ia")
     junk<-c(undec,ia)
     if(!is.null(junk)){
         junk.in.est<-sapply(estnamlist,function(x) any(junk %in% unlist(x)))
         est<-est[!junk.in.est]  # only model terms without junk
         estnamlist<-estnamlist[!junk.in.est]
     }

     est<-ifelse(is.na(est),0,est)
     o<-objnames
     ee<-estnamlist
     which.obj.in.est<-na.omit(unlist(lapply(ee, function(x) o[match(unlist(x),o)])))


     if(objcovmodel){

        # estnamlist without subject covariate terms
        estnamlist.wosubj<-lapply(estnamlist, function(x) intersect(objnames,unlist(x)))
        # estnamlist without object/objectcov terms
        estnamlist.subj<-lapply(estnamlist, function(x) setdiff(unlist(x),objnames))

        newest<-NULL
        newestnamlist<-NULL
        for (i in 1:length(estnamlist.wosubj)){
           if (any(objcovnames %in% estnamlist.wosubj[[i]])){
               # new estnamlist elements
               newenlel<-as.list(objects)
               newenlel<-lapply(newenlel,function(x) c(x, estnamlist.subj[[i]]))
               newestnamlist<-c(newestnamlist, newenlel)

               # new est elements
               idxterms<-unlist(estnamlist.wosubj[[i]])
               newestdes<-cbind(objcovs[,idxterms],est[i])
               newestel<-apply(newestdes,1,prod)
               names(newestel)<-sapply(newenlel, paste, collapse=":")
               newest<-c(newest, newestel)

            } else {
               newestnamlist<-c(newestnamlist, estnamlist[[i]])
               newest<-c(newest,est[i])
            }
        }

        ##if (OBJCOVMODEL)
        est<-newest
        estnamlist<-lapply(newestnamlist, function(x) list(x))

        objnum<-seq_along(objects)
        names(objnum)<-objects

        ee<-estnamlist
        which.obj.in.est<-na.omit(unlist(lapply(ee, function(x) objects[match(unlist(x),objects)])))
     }


     ### matrix with subj.cov indices and estimates
     subjcov.names<-names(maxlev)
     ldima<-length(dim(a))
     if(maxlev[1]==1)ldima<-1
     estmat<-matrix(0,nrow=length(est),ncol=ldima+1)
     rownames(estmat)<-names(est)
     colnames(estmat)<-c("OBJ",subjcov.names,"Value")
     estmat[,ncol(estmat)]<-est
     estmat[,1]<-objnum[which.obj.in.est]

     for (scovnam in subjcov.names){           # loop all subj covs
       for (i in 1:length(estnamlist)){        # loop all estimate terms
          trms<-estnamlist[[i]][[1]]           # extract single terms
          idx<-grep(scovnam,trms)              # look which single term corresponds to current subj cov
          lev<-0
          if(length(idx)>0)                    # if current subj cov is in estimate term
            lev<-as.numeric(sub(scovnam,"",trms[idx]))  # extract level of subj cov
          estmat[i,scovnam]<-lev               # write level into estmat
       }
     }


     ## fill array
     # expand estmat
     estmat.exp<-estmat

     expand<-function(i){
        row<-estmat.exp[i,,drop=FALSE]
        if (row[,scov]==0) {
          ##row<-expand.dfr(row,maxlev[scov])
           row<-row[rep(seq_len(nrow(row)), maxlev[scov]), ]
          row[,scov]<-1:maxlev[scov]
        }
        row
     }
     for (scov in subjcov.names){
        estmat.exp<-do.call("rbind",lapply(1:nrow(estmat.exp), expand))
     }
     e<-estmat.exp

     # array entries
     for (i in 1:nrow(e)){
        a[e[i,1:(length(subjcov.names)+1),drop=FALSE]]<-a[e[i,1:(length(subjcov.names)+1),drop=FALSE]]+e[i,"Value"]
        #cat(a[e[i,1:(length(subjcov.names)+1),drop=FALSE]],e[i,1:(length(subjcov.names)+1),drop=FALSE],"\n")
     }

     ## array into lambda/worth matrix
     estmat<-matrix(a,nrow=dim(a)[1])

     # labelling
     subjdes<-gfac2(maxlev)
     colnames(subjdes)<-names(maxlev)
     if(length(subjcov.names>1)) {
        x<-subjdes
        ## labels for cov groups
        xx<-mapply(function(x,y)paste(x,y,sep=""), colnames(x),data.frame(x))
        gr.labels <-apply(xx,1,paste,collapse=":")
        colnames(estmat) <- gr.labels
     } else
        colnames(estmat) <- "estimate"

     rownames(estmat)<-objects

      ## in case of objcovs the rows are collapsed according to obj covs
      if (!is.null(objcovs)) {
          # collapse rownames of estmat
          u<-1:length(objects)
          names(u)<- as.numeric(factor(apply(estmat,1,paste,collapse="")))
          nam<-aggregate(objects,list(u[names(u)]),paste)$x

          o<-as.data.frame(objcovs[,used.objcovs])
          colnames(o)<-used.objcovs
          leg<-aggregate(rownames(estmat),o,paste, simplify=FALSE)
          ngr<-ncol(leg)- length(used.objcovs)+1

          ##namleg<-paste("objgr",1:ngr,sep="")
          namleg<-"objects in groups"
          names(leg)[(length(used.objcovs)+1):ncol(leg)]<-namleg  # probably not neccessary (is only of length 1)
          estmat<-unique(estmat)

          if(is.matrix(nam)){
            rnam<-apply(nam,1,paste,sep="", collapse=",")
            if (nrow(estmat)==length(rnam))
                rownames(estmat)<-rnam
          } else {
            rnam<-lapply(nam,paste,sep="", collapse=",")
            if (nrow(estmat)==length(rnam))
                rownames(estmat)<-rnam
          }

      }

      if (exists("leg")) attr(estmat, which="objtable")<-leg


     worthmat<-apply(estmat,2,function(x) exp(2*x)/sum(exp(2*x)))
     attr(worthmat, which="objtable")<- attr(estmat, which="objtable")
     if(ncol(worthmat)==1)   colnames(worthmat) <- "worth"

     switch(outmat,
         "lambda" = {class(estmat)<-c("wmat",class(estmat));return(estmat)},
         "worth" = {class(worthmat)<-c("wmat",class(worthmat));return(worthmat)},
        # "est" = return(lambda.mat),
        stop("     outmat must be either 'worth' or 'lambda'\n")
     )

#     switch(outmat,
#        "lambda" = return(estmat),
#        "worth" = return(worthmat),
#        stop("     outmat must be either 'worth' or 'lambda'\n")
#     )

  ## end pattdesignmodel

  } else {

  ## fitobj from patt*.fit

     obj<-fitobj
     envList<-obj$envList
     ncovpar<-envList$ncovpar

     Tmod<- regexpr("T",envList$resptype)>0               # check if time model

     if(!Tmod){                                           # if not a time model
          tpoints<-1
          nobj<-obj$envList$nobj
          npar<-(nobj - 1) * ncovpar
          lambda<-obj$coefficients[1:npar]
          lmat<-matrix(lambda,nrow=nobj-1)
          lmat<-rbind(lmat,rep(0,ncol(lmat)))
     } else {                                             # time model
          tpoints<-envList$tpoints
          nobj<-(envList$nitems-1)
          npar<-nobj*ncovpar*tpoints
          lambda<-obj$coefficients[1:npar]
          lmat<-matrix(lambda,nrow=nobj)
          lmat<-rbind(lmat,rep(0,ncol(lmat)))
          lmat<-matrix(lmat,ncol=tpoints)
          nobj<-envList$nitems*tpoints
     }

     covlevels<-obj$envList$covlevels
     if (is.null(covlevels)) covlevels<-1

     ## preference parameters (summed lambdas) for cov groups
     struct <- unique(obj$envList$covdesmat)

     if (ncol(struct)==0) struct<-as.matrix(1)
     # summation matrix for objects
     dd<-diag(nobj)
     sum.mat <- struct %x% dd
     # sum up
     group.est <- sum.mat %*% as.vector(lmat)

     ## labels for cov groups
##     if(!Tmod){                                           # if not a time model - rubbish rh 2011-10-25
          x<-obj$envList$model.covs
          if (is.null(x)) {
               gr.labels <- "estimate"
          } else{
               xx<-mapply(function(x,y)paste(x,y,sep=""), colnames(x),data.frame(x))
               gr.labels <-apply(xx,1,paste,collapse=":")
          }
##     } else
##          gr.labels <- ""

     mltp<-2
     worthmatrix<-NULL
     est<-matrix(group.est,nrow=nobj/tpoints)

     ## worth matrix
     for (i in 1:ncol(est)) {

        # worth parameters

        worth<-rep(0,nobj/tpoints)
        coeff<-est[,i]
        worthdenominator<-0
        for (j in 1:(nobj/tpoints)) {
          worthdenominator<-worthdenominator+exp(mltp*coeff[j])
        }
        for (j in 1:(nobj/tpoints)) {
          worth[j]<-exp(mltp*coeff[j])/worthdenominator
        }
        worthmatrix<-cbind(worthmatrix,worth)
     }

     if (is.null(obj.names)){
        obj.names<-obj$envList$obj.names[1:(nobj/tpoints)] # default: only names for first time point are used
     } else {
        obj.names<-obj.names[1:(nobj/tpoints)]
     }
     ## label worth matrix
     if(Tmod) {
        if(gr.labels[1]=="") Tlabel<-"T" else Tlabel=":T"
        worth.names<-paste(rep(gr.labels,rep(tpoints,length(gr.labels))),paste(Tlabel,1:tpoints,sep=""),sep="")
        colnames(worthmatrix)<-worth.names
        rownames(worthmatrix)<-obj.names
     } else {
        colnames(worthmatrix)<-gr.labels
        rownames(worthmatrix)<-obj.names
     }

     colnames(est) <- colnames(worthmatrix)
     rownames(est) <- rownames(worthmatrix)

     #class(worthmatrix) <- c("pattW")                         #class: pattern worth
     #worthmatrix

     switch(outmat,
         "lambda" = {class(est)<-c("wmat",class(est));return(est)},
         "worth" = {class(worthmatrix)<-c("wmat",class(worthmatrix));return(worthmatrix)},
        # "est" = return(lambda.mat),
        stop("     outmat must be either 'worth' or 'lambda'\n")
     )

   } # end output from patt*fit
}
