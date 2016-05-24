`patt.design` <- local({

ENV<-new.env()  # environment local to patt.design - no objects in .GlobalEnv

###function(obj,dfr=NULL)
function(obj, nitems=NULL, objnames="", objcovs=NULL,
         cat.scovs=NULL, num.scovs=NULL,
         resptype="paircomp", reverse=FALSE,
         ia=FALSE, casewise=FALSE,
         ...)

###       cov.sel="", blnIntcovs=FALSE,
###       blnGLIMcmds=FALSE,glimCmdFile="",outFile="", intFile="")
{
  env2<-new.env()

  ### 2011-03-25 deprecated options, for backwards compatibility
  dots<-as.list(substitute(list(...)))[-1]
  nn<-names(dots)
  for (i in seq(along=dots)) assign(nn[i],dots[[i]])

  if(!exists("blnIntcovs")) blnIntcovs=ia
  if(!exists("blnRevert")) blnRevert=reverse
  if(!exists("cov.sel")) cov.sel=""

  if(!exists("blnGLIMcmds")) blnGLIMcmds=FALSE
  if(!exists("glimCmdFile")) glimCmdFile=""
  if(!exists("outFile")) outFile=""
  if(!exists("intFile")) intFile=""

  ###cov.sel<-eval(cov.sel)
  blnGLIMcmds <- ifelse(blnGLIMcmds=="T",TRUE,FALSE)

  if (blnGLIMcmds){
      stop("\nGLIMcmds deprecated. Please use preftest < 0.8-24!\n\n", call.=FALSE)
  }

  # for numeric subject covariates
  if (!is.null(num.scovs)) casewise <- TRUE

  # for backwards compatibility:
  if (!is.null(cat.scovs))
      cov.sel<-cat.scovs


  #' datafile     = "",       # dataframe used
  #' nitems       = 4,
  #' blnRevert    = FALSE,
  #' blnReducecat = TRUE,
  #' blnIntcovs   = TRUE,
  #' resptype    = "ranking",
  #'
  #' cov.sel      = "",
  #'
  #' blnGLIMcmds  = FALSE,    # no GLIM output
  #' glimCmdFile  = "",
  #' outFile      = "",
  #' intFile      = ""

  # default ctrl object
  ctrl<-list(datafile="", nitems=NULL, objnames="", resptype="paircomp", blnRevert=FALSE, cov.sel="",
       blnIntcovs=FALSE, blnGLIMcmds=FALSE,glimCmdFile="",outFile="", intFile="")

  dfr<-NULL
  if(is.character(obj)){                    ## datafilename supplied
          ctrl$datafile    <-  obj
          ctrl$nitems      <-  nitems
          ctrl$objnames    <-  objnames
          ctrl$resptype    <-  resptype
          ctrl$blnReducecat <- TRUE         # always TRUE until clarification
          ctrl$blnRevert   <-  blnRevert
          ctrl$cov.sel     <-  cov.sel
          ctrl$blnIntcovs  <-  blnIntcovs
          ctrl$blnGLIMcmds <-  blnGLIMcmds
          ctrl$glimCmdFile <-  glimCmdFile
          ctrl$outFile     <-  outFile
          ctrl$intFile     <-  intFile
  } else if(is.data.frame(obj)){            ## data frame supplied
          dfr<-obj
          ctrl$datafile    <-  ""
          ctrl$nitems      <-  nitems
          ctrl$objnames    <-  objnames
          ctrl$resptype    <-  resptype
          ctrl$blnReducecat <- TRUE         # always TRUE until clarification
          ctrl$blnRevert   <-  blnRevert
          ctrl$cov.sel     <-  cov.sel
          ctrl$blnIntcovs  <-  blnIntcovs
          ctrl$blnGLIMcmds <-  blnGLIMcmds
          ctrl$glimCmdFile <-  glimCmdFile
          ctrl$outFile     <-  outFile
          ctrl$intFile     <-  intFile
  } else if (is.list(obj)) {                ## ctrl object
          for (i in names(obj))             #       replaces default values
                  ctrl[[i]]<-obj[[i]]
          ctrl$blnReducecat <- TRUE         # always TRUE until clarification
  } else {
          stop("first argument must be either ctrlobject, datafilename or dataframe")
  }

     ####from pcpatt0 - obsolete
     ####if (!is.list(ctrl))
     ####   stop(paste(deparse(substitute(ctrl)),"is not a list object - see help for pcpatt0"))

     ## initialising
     ##
     # removed 2011-01-20 cat("initialising...\n")
     flush.console()

     for (i in 1:length(ctrl))
        do.call("assign",list(names(ctrl)[i],ctrl[[i]],ENV))

     init(dfr)    # initialisation function
     dat<-ENV$dat

     ## all possible response patterns and/or difference patterns
     ##
     # removed 2011-01-20 cat("generating response patterns...\n")
     flush.console()

     resptype <- get("resptype",ENV)

     if (resptype == "rating") {
        generateLpatterns(env2)
     } else if(resptype == "ranking") {
        generateRpatterns(env2)
     } else if(resptype == "paircomp") {
        generatePCpatterns(env2)
     }
     datStr   <- env2$datStr
     dpattStr <- env2$dpattStr
     npatt    <- env2$npatt
     diffs    <- env2$diffs
     blnUndec <- env2$blnUndec

     ## designmatrix-kernel for objects, undecided/categories, interactions
     ##
     # removed 2011-01-20 cat("setting up the design matrix...\n")
     blnIntcovs <-get("blnIntcovs",ENV)
     onedesign<-designkernel(diffs,blnUndec,blnIntcovs,env2)

     # tidy up
     rm(diffs)
     rm(dat,envir=ENV)
     gc(FALSE)


     nsubj <-get("nsubj",ENV)

     ## complete design for categorical subject covariates
     ##
     blnSubjcov <-get("blnSubjcov",ENV)
     if (!blnSubjcov) {
         ## no subject covariates
         ##
         totlev <- 1
         ones.totlev <- 1
     } else {
         ## subject covariates
         ##
         covlevels <-get("covlevels",ENV)
         covnames  <-get("covnames",ENV)
         cov.case  <-get("cov.case",ENV)
         ncov      <-get("ncov",ENV)

         if (!casewise) {
              totlev <- prod(covlevels)
              # vector for kronecker products (to stack design matrix)
              ones.totlev<-rep(1,totlev)
              indx <- ncov:1
              if (ncov == 1) {                          # only 1 covariate
                   baslev <- npatt
              } else {                                  # >1 covariates
                   baslev <- c(covlevels[2:ncov],npatt)
              }
              levmult <- rev(cumprod(baslev[indx]))
         } else {
              totlev <- nsubj
         }

         # transform subject covariates data into covariate vectors
         scov<-matrix(0,nrow=totlev*npatt,ncol=ncov)
         colnames(scov)<-covnames
         scov<-data.frame(scov)

         if (!casewise) {
            for (j in 1:ncov)
              scov[,j]<-gl(covlevels[j],levmult[j],totlev*npatt)
         } else {
            for (j in 1:ncov)
              scov[,j]<-factor(rep(cov.case[,j], each=npatt))
         }
     }



     ## extension of design-kernel in case of subject covariates
     ##
     if(casewise) ones.totlev<-rep(1,nsubj)
     design<-ones.totlev %x% onedesign
     colnames(design)<-colnames(onedesign)
     # tidy up
     rm(onedesign)
     gc(FALSE)


     ## calculate response pattern frequencies based on
     ## comparison of string representation of
     ## possible patterns and observed patterns
     ##
     # removed 2011-01-20 cat("calculating response pattern frequencies...\n")
     #                    flush.console()

     # count occurrence of patterns into y
     # (according to covariates if blnSubjcov==TRUE or casewise==TRUE)
     idx<-match(datStr,dpattStr)

     if (casewise){
       idx<-idx+(0:(length(datStr)-1)*length(dpattStr))
     } else {
     if (blnSubjcov)
       #idx<-idx+t(colSums(apply(as.matrix(cov.case)-1,1,"*",levmult)))
       idx<-idx+apply(as.matrix(cov.case)-1,1,function(x) sum(x*levmult))
     }
     y<-tabulate(idx, nrow(design))

     ## old versions count response patterns
     #y<-rep(0,totlev * npatt)
     #cov.addr <- 0       # for case blnSubjcov==FALSE
     #for (i in 1:nsubj) {
     #    if (blnSubjcov)
     #       cov.addr <- sum((cov.case[i,]-1)*levmult)
     #    j <- match(datStr[i],dpattStr)
     #    y[j+cov.addr]<-y[j+cov.addr]+1
     #}
     ##alternatively:
     #j<-sapply(datStr,function(x)match(x,dpattStr),USE.NAMES=FALSE)
     #if(blnSubjcov){
     #    cov.addr<-apply(cov.case,1,function(x)crossprod((x-1),levmult))
     #    tb<-table(j+cov.addr)
     #} else {
     #    tb<-table(j)
     #}
     #y[as.numeric(names(tb))]<-tb
     #rm(j,cov.addr,tb)

     #if (casewise)
        ###dpattStr



     # tidy up
     rm(datStr,dpattStr)
     rm(cov.case, envir=ENV)
     gc(FALSE)


     ## prepare dataframe
     ##
     dm<-data.frame(cbind(y,design))
     if (blnSubjcov){
           dm<-data.frame(dm,scov)
           rm(scov,cov.case)
     }

     ## numeric subject covariates: casewise==TRUE
     if (!is.null(num.scovs)){
       num.scovs<-unique(num.scovs) # in case somebody uses same variable twice
       ncov.case<-as.matrix(dat[,c(num.scovs)], drop=FALSE)
       ncov<-matrix(,nrow=nsubj*npatt,ncol=length(num.scovs))
       for (j in seq(along=num.scovs))
              ncov[,j]<-rep(ncov.case[,j], each=npatt)
       colnames(ncov)<-num.scovs
       CASE<-gl(nsubj, npatt)
       dm<-data.frame(cbind(dm,ncov,CASE=CASE))
       rm(ncov.case, ncov)
     }


     objnames <- colnames(design[,1:nitems])

     ### new: object specific covariates
     if(!is.null(objcovs)){
        if(!is.data.frame(objcovs)) stop("object specific covariates must be specified as a data frame")
        if(any(sapply(objcovs,is.factor))) stop("object specific covariates must not be factors (use model.matrix first)")
        dm<-data.frame(dm, design[,1:nitems] %*% as.matrix(objcovs)) # add obj covs to design data frame
        objcovnames <- colnames(objcovs)

        #varnames <- c(varnames, objcovnames)
        rownames(objcovs) <- objnames
        attr(dm, which="objcovs")<-as.matrix(objcovs, drop=FALSE)
     }

     ### new: attributes in design data frame
     attr(dm, which="objnames")<-objnames
     #if(!is.null(cat.scovs)) attr(dm, which="cat.scovs")<-cat.scovs
     if(length(cov.sel)>0 && cov.sel[1]!="") attr(dm, which="cat.scovs")<-covnames
     if(!is.null(num.scovs)) attr(dm, which="num.scovs")<-num.scovs
     if(blnUndec) attr(dm, which="undec")<-compnames(nitems)
     if(blnIntcovs) attr(dm, which="ia")<-env2$labels.intpars

     class(dm)<-c("data.frame","pattdes")

     # tidy up
     rm(y,design)
     gc(FALSE)

     ### generate files for GLIM (not fixed from 0.8-24 onwards)
     ###
     #if (get("blnGLIMcmds",ENV)){
     #    cat("\nGLIMcmds deprecated. Please use preftest < 0.8-24!\n\n")
     #    # nintcovs.out<-40 # max number of values/line in interaction output file
     #    # writeGLIMcmds(dm,blnUndec,
     #    #     blnIntcovs,outFile,ncov,env2$nintpars,intFile,nintcovs.out,glimCmdFile,
     #    #     covnames,covlevels,objnames,ENV$undecnames)
     #    ##writeGLIMcmds(dm,blnUndec,ENV)
     #} else {
     #    return(dm) # R output only if blnGLIMcmds==FALSE
     #}

     # removed 2011-01-20 cat("Done\n\n")
     dm
}
}) # end local
