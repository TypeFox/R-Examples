#' Compute normalised prediction distribution errors
#' 
#' These functions compute normalised prediction distribution errors (npde) and
#' optionally prediction discrepancies (pd). \code{npde} asks the user the name
#' and structure of the files containing the data, using \code{pdemenu}, while
#' \code{autonpde} takes these variables and others as arguments.
#' 
#' Both functions compute the normalised prediction distribution errors (and/or
#' prediction discrepancies) in the same way. \code{npde} is an interactive
#' function whereas \code{autonpde} takes all required input as arguments.
#' 
#' When the computation of npde fails because of numerical problems, error
#' messages are printed out, then pd are computed instead and graphs of pd are
#' plotted so that the user may evaluate why the computation failed.
#' 
#' @aliases autonpde npde
#' @usage autonpde(namobs, namsim, iid, ix, iy, imdv = 0, icens = 0,
#' icov = 0, iipred = 0, boolsave = TRUE, namsav = "output", type.graph = "eps",
#' verbose = FALSE, calc.npde=TRUE, calc.pd=TRUE, decorr.method = "cholesky",
#'  cens.method = "cdf", units = list(x="",y=""), detect=FALSE, ties=TRUE)
#' @usage npde()
#' @param namobs name of the file containing the observed data, or a dataframe
#' containing the observed data (in both cases, the column containing the
#' various data required for the computation of the pde can be set using the
#' arguments \code{iid},\code{ix} and \code{iy} below)
#' @param namsim name of the file containing the simulated data, or a dataframe
#' containing the simulated data (the program will assume that subject ID are
#' in column 1 and simulated Y in column 3, see User Guide)
#' @param iid name/number of the column in the observed data containing the patient
#' ID; if missing, the program will attempt to detect a column named id
#' @param ix name/number of the column in the observed data containing the
#' independent variable (X); ; if missing, the program will attempt to detect a column named X
#' @param iy name/number of the column in the observed data containing the dependent
#' variable (Y); if missing, the program will attempt to detect a column with the response
#' @param imdv name/number of the column containing information about missing data
#' (MDV), defaults to 0 (column not present)
#' @param icens name/number of the column containing information about censored data
#' (cens), defaults to 0 (column not present)
#' @param icov name/number of the column(s) containing covariate information
#' defaults to 0 (no covariates)
#' @param iipred name/number of the column(s) with individual predictions
#' (ipred), defaults to 0 (individual predictions not available)
#' @param units a list with components x, y and cov (optional), specifying the
#' units respectively for the predictor (x), the response (y), and the covariates 
#' (a vector of length equal to the number of covariates). Units will default to (-) if not given.	
#' @param detect a boolean controlling whether automatic recognition of columns in
#' the dataset is on, defaults to FALSE
#' @param boolsave a boolean (TRUE if graphs and results are to be saved to a
#' file, FALSE otherwise), defaults to TRUE	
#' @param namsav name of the files to which results are to be saved (defaults
#' to "output", which will produce a file called output.eps (if the default
#' format of postscript is kept, see type.graph) for the graphs and a file
#' called output.npde for the numerical results (see value)
#' @param type.graph type of graph (one of "eps","jpeg","png","pdf"), defaults
#' to postscript ("eps")
#' @param calc.npde a boolean (TRUE if npde are to be computed, FALSE otherwise),
#' defaults to TRUE
#' @param calc.pd a boolean (TRUE if pd are to be computed, FALSE otherwise), defaults
#' to TRUE
#' @param cens.method a character string indicating the method used to handle 
#' censored data (see \code{\link{npde.cens.method}})
#' defaults to cdf
#' @param decorr.method a character string indicating the method used to decorrelate
#' observed and simulated data in the computation of npde (see \code{\link{npde.decorr.method}})
#' defaults to cholesky
#' @param ties a boolean (if FALSE, the distributions of pd and npde are smoothed by jittering the values so that there are no ties), defaults to TRUE	
#' @param verbose a boolean (TRUE if messages are to be printed as each subject is
#' processed, FALSE otherwise), defaults to FALSE
#' @return An object of class \code{\link{NpdeObject}}
#' 
#' @details The function also prints out the characteristics of the distribution of the
#' npde (mean, variance, skewness and kurtosis) as well as the results of the
#' statistical tests applied to npde. In addition, if boolsave is TRUE, two files
#' are created: 
#' \describe{
#' \item{results file}{the numerical results are saved in a file
#' with extension .npde (the name of which is given by the user). The file
#' contains the components id, xobs, ypred, npde, pd stored in columns}
#' \item{graph file}{the graphs are saved to a file with the same name as the
#' results file, and with extension depending on the format.}
#' }
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde.graphs}}, \code{\link{gof.test}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.
#' Mentre. Metrics for external model evaluation with an application to the
#' population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research},
#' 23:2036--49, 2006.
#' @keywords models
#' @export
#' @examples
#' 
#' data(theopp)
#' data(simtheopp)
#' 
#' # Calling autonpde with dataframes
#' 
#' x<-autonpde(theopp,simtheopp,1,3,4,boolsave=FALSE)
#' x
#' 
#' # Calling autonpde with names of files to be read from disk
#' 
#' write.table(theopp,"theopp.tab",quote=FALSE,row.names=FALSE)
#' write.table(simtheopp,"simtheopp.tab",quote=FALSE,row.names=FALSE)
#' x<-autonpde(namobs="theopp.tab", namsim="simtheopp.tab", iid = 1,
#' ix = 3, iy = 4, imdv=0, boolsave = FALSE)
#' 
#' head(x["results"]["res"])
#' 

npde<-function() {
    xinput<-pdemenu()    

    xdat<-npdeData(name.data=xinput$namobs,header=TRUE,name.group=xinput$iid, name.predictor=xinput$ix,name.response=xinput$iy,name.miss=xinput$imdv, name.cens=xinput$icens,name.ipred=xinput$iipred,name.covariates=xinput$icov, detect=xinput$detect)
    cat("Simulated data:",xinput$namsim,"\n")
    xsim<-npdeSimData(npde.data=xdat,name.simdata=xinput$namsim)
    	
    opt<-list(boolsave=xinput$boolsave,namsav=xinput$namfile, type.graph=xinput$type.graph,verbose=xinput$verbose,calc.npde=xinput$calc.npde, calc.pd=xinput$calc.pd,decorr.method=xinput$decorr.method,cens.method=xinput$cens.method,ties=xinput$ties)
    npde.obj<-new(Class="NpdeObject",data=xdat,sim.data=xsim,options=opt)
    npde.obj["prefs"]<-set.plotoptions(npde.obj)

    xret<-npde.main(npde.obj)
# Saving results
    if(npde.obj["options"]$boolsave) {
      npde.save(xret)
      npde.graphs(xret)
    }
    invisible(xret)
#    if (xinput$output==TRUE) invisible(xret) else return(NULL)
}

#' @export
autonpde<-function(namobs,namsim,iid,ix,iy,imdv=0,icens=0,icov=0, iipred=0,boolsave=TRUE,namsav="output",type.graph="eps",verbose=FALSE, calc.npde=TRUE,calc.pd=TRUE,decorr.method="cholesky",cens.method="cdf", units=list(x="",y=""), detect=FALSE, ties=TRUE) {
# output is deprecated, now using invisible
    if(is.data.frame(namobs)) namobs<-deparse(substitute(namobs))
    if(is.data.frame(namsim)) namsim<-deparse(substitute(namsim))
    if(missing(iid)) iid<-""
    if(missing(ix)) ix<-""
    if(missing(iy)) iy<-""
    xdat<-npdeData(name.data=namobs,header=TRUE,name.group=iid, name.predictor=ix,name.response=iy,name.covariates=icov,name.miss=imdv, name.cens=icens,name.ipred=iipred,units=units,verbose=verbose)
    xsim<-npdeSimData(npde.data=xdat,name.simdata=namsim,verbose=verbose)
#    cat("Data and sim.data object created\n")
#    if(cens.method=="cdf" & !calc.pd) {
    if(cens.method!="omit" & !calc.pd) {
      calc.pd<-TRUE
      cat("To compute npde with the",cens.method," method, pd need to be computed first, changing to calc.pd.\n")
    }
    opt<-list(boolsave=boolsave,namsav=namsav,type.graph=type.graph, verbose=verbose,calc.npde=calc.npde, calc.pd=calc.pd,decorr.method=decorr.method, cens.method=cens.method, ties=ties)
    npde.obj<-new(Class="NpdeObject",data=xdat,sim.data=xsim,options=opt)
#    cat("Npde object created\n")
    npde.obj["prefs"]<-set.plotoptions(npde.obj)
#    cat("Options added\n")
#    showall(npde.obj)
    xret<-npde.main(npde.obj)
    if(npde.obj["options"]$boolsave) {
      npde.save(xret)
      npde.graphs(xret)
    }
    invisible(xret)
}

#' Interactive menu to set the options for the npde() function
#' 
#' Interactive menu to set the options for the npde() function
#' 
#' @return A list with the information needed to compute the pd/npde
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde.graphs}}, \code{\link{gof.test}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.
#' Mentre. Metrics for external model evaluation with an application to the
#' population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research},
#' 23:2036--49, 2006.
#' @keywords internal

pdemenu<-function() {
    ick<-0
    while(ick==0) {
      namobs<-readline(prompt="Name of the file containing the observed data : ")
      datobs<-try(read.table(namobs,na.strings=c(".","NA"),nrows=10))
      if(class(datobs)!="try-error") ick<-1 else 
         cat("\n      File",namobs,"does not exist.\n")
    }

    cok<-readline(prompt="Automatic recognition of columns in the dataset (y/Y) [default=no] ? ")
    if(tolower(cok)=="y"|tolower(cok)=="yes") detect<-FALSE else detect<-TRUE
    if(!detect) {
      iid<-1;ix<-2;iy<-3;imdv<-0
      cat("I'm assuming file",namobs,"has the following structure :\n")
      cat("        ID X Y ...\n")
      cat("and does not contain a column signaling missing data.\n")
      cok<-readline(prompt="To keep, press ENTER, to change, type any letter : ")
      if(cok=="") 
        cat("Keeping this structure\n") 
      else {
      cat("In the following, you may enter the name of the column or the column number\n")
      iid<-as.integer(readline(prompt="         Column with ID information ? "))
      ix<-as.integer(readline(prompt="         Column with X (eg time) information ? "))
      iy<-as.integer(readline(prompt="         Column with Y (eg DV) information ? "))
      imdv<-as.integer(readline(prompt="         Column signaling missing data (eg MDV, press ENTER if none) ? "))
      if(is.na(imdv)) imdv<-0
      icens<-as.integer(readline(prompt="         Column signaling censoring (eg CENS, press ENTER if none) ? "))
      if(is.na(icens)) icens<-0
      iipred<-as.integer(readline(prompt="         Column with individual predictions (eg ipred, press ENTER if none) ? "))
      if(is.na(iipred)) iipred<-0
      icov<-c()
      ic<-as.integer(readline(prompt="         Columns with covariates (eg WT; enter one at a time, press ENTER if none or when finished) ? "))
      while(!is.na(ic)) {
        icov<-c(icov,ic)
        ic<-as.integer(readline(prompt="next covariate (press ENTER when finished)"))
      }
      if(is.null(icov)) icov<-0
    }
    }
    ick<-0
    while(ick==0) {
      namsim<-readline(prompt="Name of the file containing the simulated data : ")
      datobs<-try(read.table(namsim,na.strings=c(".","NA"),nrows=10))
      if(class(datobs)!="try-error") ick<-1 else 
         cat("\n      File",namsim,"does not exist.\n")
    }
    cok<-readline(prompt="Do you want results and graphs to be saved to files (y/Y) [default=yes] ? ")
    if(cok=="y"|cok=="Y"|cok=="yes"|cok=="") boolsave<-TRUE else boolsave<-FALSE
    type.graph<-"eps";namegr<-nameres<-"";decorr.method<-"cholesky";cens.method<-"cdf"
    if(boolsave) {
      cat("Different formats of graphs are possible :\n")
      cat("         1. Postscript (extension eps, default)\n")
      cat("         2. JPEG (extension jpeg)\n")
      cat("         3. PNG (extension png)\n")
      cat("         4. Acrobat PDF (extension pdf)\n")
      cok<-as.integer(readline(prompt="Which format would you like for the graph (1-4) [default 1] ? "))
      if(is.na(cok) || cok>4 || cok<0) {
        cok<-1
        cat("       postscript graph selected\n")
      }
      if(cok==2) type.graph<-"jpeg"
      if(cok==3) type.graph<-"png"
      if(cok==4) type.graph<-"pdf"
      namfile<-readline(prompt="Name of the file (extension will be added, default=output):")
      namfile<-ifelse(namfile=="","output",namfile)
    }
    ick<-0
    while(ick==0) {
      cok<-readline(prompt="Do you want to compute npde (y/Y) [default=yes] ? ")
      if(cok=="y"|cok=="Y"|cok=="yes"|cok=="") {calc.npde<-TRUE;ick<-1} else 
         calc.npde<-FALSE
      cok<-readline(prompt="Do you want to compute pd (y/Y) [default=yes] ? ")
      if(cok=="y"|cok=="Y"|cok=="yes"|cok=="") {calc.pd<-TRUE;ick<-1} else 
         calc.pd<-FALSE
      if(ick==0) cat("\n Please choose to compute at least one of npde or pd.\n") 
    }
    cat("Different decorrelation methods are available:\n")
    cat("         1. Cholesky decomposition (default)\n")
    cat("         2. Inverse using diagonalisation (as in Monolix and Nonmem)\n")
    cat("         3. Cholesky followed by polar decomposition\n")
    cok<-as.integer(readline(prompt="Which method should be used for the decorrelation (1-3)  [default 1] ? "))
    if(is.na(cok) || cok>3 || cok<0) {
      cok<-1
      cat("       Cholesky method selected\n")
    }
    if(cok==2) decorr.method<-"inverse"
    if(cok==3) decorr.method<-"polar"    
    cat("Method used to handle censored observations:\n")
    cat("         1. omit: pd will be set to NaN for missing data\n")
    cat("         2. cdf: pd will be imputed using a random sample from U(0,p_LOQ) where p_LOQ is the probability, according to the model, that a given observation is less than LOQ (default)\n")
    cat("         3. loq: an observation below the LOQ will be imputed to the LOQ\n")
    cat("         4. ppred: an observation below the LOQ will be imputed to the population model prediction\n")
    cat("         5. ipred: an observation below the LOQ will be imputed to the individual model prediction\n")
    cok<-as.integer(readline(prompt="Which method should be used (1-5)  [default 2] ? "))
    if(is.na(cok) || cok>5 || cok<0) {
      cok<-2
      cat("       cdf method selected\n")
    }
    if(cok==1) cens.method<-"omit"
    if(cok==3) cens.method<-"loq"
    if(cok==4) cens.method<-"ppred"
    if(cok==5) cens.method<-"ipred"
    if(cens.method!="omit" & !calc.pd) {
      calc.pd<-TRUE
      cat("To compute npde with the",cens.method," method, pd need to be computed first, changing calc.pd to TRUE\n")
    }
    verbose<-FALSE
    if(calc.npde) {
    cok<-readline(prompt="Do you want a message printed as the computation of npde begins in a new subject (y/Y) [default=no] ? ")
    if(cok=="y"|cok=="Y"|cok=="yes") verbose<-TRUE 
    }
    ties<-TRUE
    cok<-readline(prompt="Do you want to allow different observations to have the same value of pd/npde (y/Y) [default=yes, if no, pd and npde will be jittered] ? ")
    if(tolower(cok)=="n"|tolower(cok)=="no") ties<-TRUE 
    return(list(namobs=namobs,namsim=namsim,iid=iid,ix=ix,iy=iy,imdv=imdv,icens=icens, iipred=iipred,icov=icov,boolsave=boolsave,type.graph=type.graph,namfile=namfile, calc.pd=calc.pd,calc.npde=calc.npde,verbose=verbose,decorr.method=decorr.method, cens.method=cens.method,detect=detect,ties=ties))
}
