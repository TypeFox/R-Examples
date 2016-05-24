llbt.worth <- function(fitobj, outmat = "worth")
{
#### updated on 2010-12-21
#### can now handle fitobjects from llbtPCfit and
#### from gnm, provided llbt.design (also new version )was used
#### object specific covariates are possible

# check if either model results from gnm (using llbt.design)
# or from llbtPC.fit

llbtdesignmodel <- FALSE
llbtPCfitmodel <- ("llbtMod" %in% class(fitobj))
if (!llbtPCfitmodel){
   design <- fitobj$data                        ## unten nocheinmal
   if("llbtdes" %in% class(design)) llbtdesignmodel <- TRUE
}
if (llbtdesignmodel == llbtPCfitmodel)
     stop("Model result must be either from llbtPC.fit or having used llbt.design")


## fitobj from llbtdesign, gnm
if (llbtdesignmodel){



     ## subject covariates

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
     if(length(eterms)>0) {
         # extract subject covariates from formula
         subjcov.names<-intersect(eterms,fterms)

         # achtung laenge 0
         if(length(subjcov.names)>0) { #subject variable in design matrix and in model
           # set up subj cov design
              levlist<-lapply(fitobj$model[subjcov.names],levels)
           maxlev<-sapply(levlist,length)
         }
           if (!exists("maxlev")){                      #subject variables in design matrix but none in model
                maxlev <- 1
                length(eterms)<-0
           }
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

     # extract category covariate names from design
     # if user has specified interaction with objects
     categ<-attr(design,"categories")
     junk<-c(categ)
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

#    switch(outmat,
#       "lambda" = return(estmat),
#       "worth" = return(worthmat),
#       stop("     outmat must be either 'worth' or 'lambda'\n")
#     )


     switch(outmat, lambda = {
         class(estmat) <- c("wmat", class(estmat))
         return(estmat)
     }, worth = {
         class(worthmat) <- c("wmat", class(worthmat))
         return(worthmat)
     }, est = {
         class(estmat) <- c("wmat", class(estmat))
         return(estmat)
     }, stop("     outmat must be either 'worth' or 'lambda'\n")
     )
  ## end llbtdesignmodel
} else {

      nobj <- fitobj$envList$nobj

      # remove category and undecided parameters from lambda
      lambda <- fitobj$coefficients[fitobj$ofInterest]
###      if(any(grep("^g[0-9]|g[0-9]$|^u$",names(lambda))))            # g1, g2 etc. do not occur in llbtPC.fit
###          lambda <- lambda[-(grep("^g[0-9]|^u$", names(lambda)))]   # undecided changed to U in llbtPC.fit
      if(any(grep("^U$",names(lambda))))
          lambda <- lambda[-(grep("^U$", names(lambda)))]
      lambda <- ifelse(is.na(lambda),0,lambda)

      # initialise lambda matrix
      npar<-length(lambda)
      lambda.mat <- matrix(, nrow=nobj, ncol=npar/nobj)

      # row and colnames for output matrix
      nam.lambda <- names(lambda)
      rownames(lambda.mat) <- nam.lambda[1:nobj]
      nam.new<-nam.lambda
      for (i  in nam.lambda[1:nobj]) # remove obj names from nam.new
          # nam.new<-gsub(paste("^",i,"$|^",i,"[^0-9]:+|[:alnum:]*:?(",i,"[^0-9]{0,2})$",sep=""),"",nam.new)
          nam.new<-gsub(paste("^",i,"$|^",i,":",sep=""),"",nam.new)
      nn<-nam.new
      nl<-nam.lambda
      colnames(lambda.mat) <- unique(nam.new)


      # fill lambda matrix lambda.mat with correct entries

      # case with covs
      if(npar>nobj){

          # sort names for model terms according to columns of model matrix
          covmat<-unique(model.matrix(fitobj$envList$formel,data=fitobj$data))
          cnam<-colnames(covmat)
          nam.cov<-sapply(1:ncol(covmat), # sort variable names within interactions terms
              function(i){txt<-sort(unlist(strsplit(cnam[i],":"))); #
              paste(txt,collapse=":")}
          )

          A<-nam.lambda[1:nobj]
          B<-nam.cov
          B[1]<-""
          lo<-as.vector(t(outer(A,B,paste,sep=":")))
          lo<-sub(":$","",lo)
          # estimates in same order as in model matrix
          lambda<-lambda[match(lo,names(lambda))]

          lambda.mat<-matrix(lambda,nrow=nobj,byrow=T)

          colnames(lambda.mat) <- unique(nam.new)
          rownames(lambda.mat) <- nam.lambda[1:nobj]

          ## calculate sums of estimates accordig to covariate groups
          lambda.vec <- as.vector(lambda.mat)
          lambda.groups.mat <- matrix((covmat %x% diag(nobj)) %*% lambda.vec, nrow=nobj)


          ## labels for cov groups
          x<-fitobj$envList$model.covs
          xx<-mapply(function(x,y)paste(x,y,sep=""), colnames(x),data.frame(x))
          gr.labels <-apply(xx,1,paste,collapse=":")
          colnames(lambda.groups.mat) <- gr.labels


          ## obsolete since option obj.names removed 2011-01-03
          #if (is.null(obj.names))
          #   rownames(lambda.groups.mat) <- nam.lambda[1:nobj]
          #else
          #   rownames(lambda.groups.mat) <- obj.names
          rownames(lambda.groups.mat) <- nam.lambda[1:nobj]

          ## worth matrix
          worth.groups.mat <- apply(lambda.groups.mat, 2, function(x) exp(2*x)/sum(exp(2*x)))

      # case w/o covs
      } else {
      lambda.mat <- matrix(lambda, nrow = nobj, dimnames = list(names(lambda)[1:nobj], "estimate"))
      lambda.groups.mat <- lambda.mat
      worth.groups.mat <- apply(lambda.groups.mat, 2, function(x) exp(2*x)/sum(exp(2*x)))
      }
}

## return
switch(outmat,
   "lambda" = {class(lambda.groups.mat)<-c("wmat",class(lambda.groups.mat));return(lambda.groups.mat)},
   "worth" = {class(worth.groups.mat)<-c("wmat",class(worth.groups.mat));return(worth.groups.mat)},
   "est" = {class(lambda.groups.mat)<-c("wmat",class(lambda.groups.mat));return(lambda.groups.mat)},
   stop("     outmat must be either 'worth' or 'lambda'\n")
)

}
