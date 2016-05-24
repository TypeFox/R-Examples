####################################################################################
####			SaemixData class - definition				####
####################################################################################

# ECO TODO: check on name validity

###############################
# Definition with initialise
setClass(
  Class="SaemixData",
  representation=representation(
    name.data="character",	# name of dataset
    header="logical",		# for file, whether has header
    sep="character",		# if file, separator
    na="character",		# if file, NA symbol(s)
    name.group="character",	# name of column with ID
    name.predictors="character",# name of column(s) with predictors 
    name.response="character",	# name of column with response
    name.covariates="character",# name of column(s) with covariates
    name.X="character",		# name of predictor used on X axis for graphs
    units="list",		# units (list with components for x, y, and cov)
    data="data.frame",		# the data (data frame with columns name.group (subject id), index (id renamed to 1:N), name.predictors (predictors), name.response (possibly transformed during fit), name.covariates); binary covariates are modified to 0/1
    ocov="data.frame",		# original scale for the covariates
    N="numeric",		# number of subjects
#    index="numeric",		# subject index (id renamed to 1:N)
#    id="character",		# subject id
#    xind="data.frame",		# matrix of predictors
#    y="numeric",		# vector of responses (possibly transformed during fit)
    yorig="numeric",		# vector of responses in original dataset
#    cov="data.frame",		# matrix of covariates
    ntot.obs="numeric",		# total number of observations (=dim(tab)[1])
    nind.obs="numeric"		# number of observations for each subject
  ),
  validity=function(object){
#    cat ("--- Checking SaemixData object ---\n")
    if (length(object@name.data)==0) {
      stop ("[ SaemixData : validation ] Please provide a name for the data (dataset or datafile on disk).")
    }
# Ici ou a la creation, detection automatique ?
    if (nchar(object@name.group)==0) {cat("Missing Id column\n")}
    if (nchar(object@name.predictors[1])==0) {cat("No predictors found\n")}
    if (nchar(object@name.response)==0) {cat("No response found\n")}
    N<-object@N
    if(length(object@data)>0) {
     if(N<2) {
       cat("Warning: There is only",N,"subjects in the dataset, the SAEM algorithm is a population algorithm designed to analyse longitudinal data from non-linear mixed effect models and may not work with too few subjects.\n")
     }
     if(length(unique(object@data[,object@name.response]))<3) {
       cat("Warning: The SAEM algorithm currently handles only continuous responses. It seems that the response --",object@name.response,"-- has too few modalities and the statistical model may not be appropriate.\n")
     }
    }
    return(TRUE)
  }
)

setClass(
  Class="SaemixRepData", # Saemix data, replicated for different chains
  representation=representation(
    N="numeric",		# number of subjects
    NM="numeric",		# number of subjects, replicated
    dataM="data.frame",		# replicated data with columns IdM, xM, yM
#    IdM="numeric",		# subject id
#    XM="data.frame",		# matrix of predictors
#    yM="numeric",		# vector of responses 
    nrep="numeric"		# number of replicates
  ),
  validity=function(object){
#    cat ("--- Checking SaemixData object ---\n")
    return(TRUE)
  }
)

setClass(
  Class="SaemixSimData", # Saemix predicted and simulated data
  representation=representation(
    N="numeric",		# number of subjects
    name.response="character",	# name of column with response
    name.X="character",		# name of predictor used on X axis for graphs
    units="list",		# units (list with components for x, y, and cov)
    data="data.frame",		# ECO TODO: do we need to keep it here ?
    nsim="numeric",		# number of simulations
    datasim="data.frame",	# simulated data with columns idsim (id in replications), irep (replication number), ypred (simulated predictions, without error), ysim (simulated data, with error)
    sim.psi="data.frame"	# simulated parameters
  ),
  validity=function(object){
#    cat ("--- Checking saemixSimData object ---\n")
    return(TRUE)
  }
)

###############################
# ECO validity ne semble pas etre appele automatiquement quand on cree un objet => il faut l'appeler dans initialize

setMethod(
  f="initialize",
  signature="SaemixData",
  definition= function (.Object,name.data,header,sep,na,name.group, name.predictors,name.response,name.covariates,name.X,units){
#    cat ("--- initialising SaemixData Object --- \n")
    if(missing(name.data)) stop ("Please provide a name for the data (dataset or datafile on disk).")
    .Object@name.data<-name.data
    if(missing(header)) header<-TRUE
    .Object@header<-header
    if(missing(sep)) sep<-""
    .Object@sep<-sep
    if(missing(na)) na<-"NA"
    .Object@na<-na
    if(missing(name.group)) {
      cat("   Missing ID identifier, assuming the ID is in column 1 of the dataset.\n")
      name.group<-"1"
    }
# ECO TODO: reconnaissance automatique (avant affectation a la valeur 2) ?
    if(missing(name.predictors)) {
      name.predictors<-"2"      
      cat("   Missing predictors identifier, assuming there is one predictor in column 2 of the dataset.\n")
    }
    if(missing(name.response)) {
      cat("   Missing response identifier, assuming the response is in column 3 of the dataset.\n")
      name.response<-"3"
    }
    if(missing(name.covariates)) name.covariates<-character()
    .Object@name.group<-name.group
    .Object@name.predictors<-name.predictors
    .Object@name.response<-name.response
    .Object@name.covariates<-name.covariates
    if(missing(units)) units<-list(x="-",y="-")
    if(is.null(units$x)) units$x<-"-"
    if(is.null(units$y)) units$y<-"-"
    ncov<-length(name.covariates)
    if(ncov>0) {
      nunit<-length(units$covariates)
      if(nunit==0) units$covariates<-rep("-",ncov)
      if(nunit>ncov) units$covariates<-units$covariates[1:ncov]
      if(nunit<ncov) {
        length(units$covariates)<-ncov
        units$covariates[(nunit+1):ncov]<-"-"
      }
    }
    .Object@units<-units
    .Object@name.X<-name.X
# Object validation
    validObject(.Object)
    return (.Object )
  }
)

# Initialize method for saemixRepData and saemixSimData
setMethod(
  f="initialize",
  signature="SaemixRepData",
  definition= function (.Object,data=NULL,nb.chains=1){
#    cat ("--- initialising SaemixData Object --- \n")
    if(is.null(data)) {
      .Object@N<-.Object@NM<-numeric(0)
      .Object@dataM<-data.frame()
    } else {
    N<-data@N
    .Object@N<-N
    .Object@nrep<-nb.chains
    .Object@NM<-N*nb.chains
    IdM<-kronecker(c(0:(nb.chains-1)),rep(N,data@ntot.obs))+rep(data@data[,"index"], nb.chains)
    yM<-rep(data@data[,data@name.response],nb.chains)
    XM<-do.call(rbind,rep(list(data@data[,data@name.predictors,drop=FALSE]), nb.chains))
    .Object@dataM<-data.frame(IdM=c(IdM),XM,yM=yM)
   }
# Object validation
#    validObject(.Object)
    return (.Object )
  }
)

setMethod(
  f="initialize",
  signature="SaemixSimData",
  definition= function (.Object,data=NULL,datasim=NULL) {
#    cat ("--- initialising SaemixData Object --- \n")
    if(!is.null(data)) {
      .Object@N<-data@N
      .Object@name.response<-data@name.response
      .Object@name.X<-data@name.X
      .Object@units<-data@units
      .Object@data<-data@data
    }
    if(is.null(data) || is.null(datasim) || dim(datasim)[1]==0) {
      .Object@datasim<-data.frame()
      .Object@nsim<-0
    } else {
      .Object@datasim<-datasim
      .Object@nsim<-dim(datasim)[1]/dim(data@data)[1]
    }
# Object validation
#    validObject(.Object)
    return (.Object )
  }
)

####################################################################################
####			SaemixData class - accesseur				####
####################################################################################

# Getteur
setMethod(
  f ="[",
  signature = "SaemixData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "name.data"={return(x@name.data)},
    "header"={return(x@header)},
    "sep"={return(x@sep)},
    "na"={return(x@na)},
    "name.group"={return(x@name.group)},
    "name.predictors"={return(x@name.predictors)},
    "name.response"={return(x@name.response)},
    "name.covariates"={return(x@name.covariates)},
    "name.X"={return(x@name.X)},
    "units"={return(x@units)},
    "data"={return(x@data)},
    "ocov"={return(x@ocov)},
    "N"={return(x@N)},
    "yorig"={return(x@yorig)},
    "ntot.obs"={return(x@ntot.obs)},
    "nind.obs"={return(x@nind.obs)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "name.data"={x@name.data<-value},
    "header"={x@header<-value},
    "sep"={x@sep<-value},
    "na"={x@na<-value},
    "name.group"={x@name.group<-value},
    "name.predictors"={x@name.predictors<-value},
    "name.response"={x@name.response<-value},
    "name.covariates"={x@name.covariates<-value},
    "name.X"={x@name.X<-value},
    "units"={x@units<-value},
    "data"={x@data<-value},
    "ocov"={x@ocov<-value},
    "N"={x@N<-value},
    "yorig"={x@yorig<-value},
    "ntot.obs"={x@ntot.obs<-value},
    "nind.obs"={x@nind.obs<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)


# For saemixRepData
setMethod(
  f ="[",
  signature = "SaemixRepData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "N"={return(x@N)},
    "NM"={return(x@NM)},
    "dataM"={return(x@dataM)},
    "nrep"={return(x@nrep)},
    stop("No such attribute\n")
   )
  }
)

setReplaceMethod(
  f ="[",
  signature = "SaemixRepData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "N"={x@N<-value},
    "NM"={x@NM<-value},
    "dataM"={x@dataM<-value},
    "nrep"={x@nrep<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

# For saemixSimData
setMethod(
  f ="[",
  signature = "SaemixSimData" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "N"={return(x@N)},
    "name.response"={return(x@name.response)},
    "name.X"={return(x@name.X)},
    "units"={return(x@units)},
    "data"={return(x@data)},
    "nsim"={return(x@nsim)},
    "sim.psi"={return(x@sim.psi)},
    "datasim"={return(x@datasim)},
    stop("No such attribute\n")
   )
  }
)

setReplaceMethod(
  f ="[",
  signature = "SaemixSimData" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "N"={x@N<-value},
    "name.response"={x@name.response<-value},
    "name.X"={x@name.X<-value},
    "units"={x@units<-value},
    "data"={x@data<-value},
    "sim.psi"={x@sim.psi<-value},
    "datasim"={x@datasim<-value},
    "nsim"={x@nsim<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

####################################################################################
####			SaemixData class - method to read data			####
####################################################################################

setMethod("read.saemixData","SaemixData",
  function(object) {
    ow <- options("warn")
    options("warn"=-1)
# ce test devrait aller dans la definition de la classe
    if(class(object@name.data)!="character") {
    cat("Please provide the name of the data (data.frame or path to file on disk) as a character string.\n")
    return("Creation of saemixData failed")
  }
    if(exists(object@name.data)) {
      cat("Using the object called",object@name.data,"in this R session as the data.\n")
      dat<-get(object@name.data)
    } else {
      cat("Reading data from file",object@name.data,"\n")
      header<-object@header
      if(is.null(header)) header<-TRUE
      sep<-object@sep
      if(is.null(sep)) sep<-""
      na.strings<-object@na
      if(is.null(na.strings)) na.strings<-"NA"
      dat<-try(read.table(object@name.data,header=header,sep=sep,na.strings=na.strings))
      if(class(dat)=="try-error") stop("The file ",object@name.data," does not exist. Please check the name and path.\n")      
      cat("These are the first lines of the dataset as read into R. Please check the format of the data is appropriate, if not, modify the na and/or sep items and retry:\n")
      print(head(dat))
    }
    if(dim(dat)[2]<3) {
      cat("The dataset does not contain enough data. The non-linear mixed effect model requires at least 3 columns, with subject ID, predictor (at least one) and response. \nPlease check the field separator, currently given as:", paste("sep=\"",object@sep,"\"",sep=""),"\n")
      return("Creation of saemixData failed")
    }
# Automatic recognition of columns 
#    ID (one of id, subject or sujet regardless of case)
#    response (one of Y, conc, concentration, resp, response regardless of case)
#    predictors (time and/or dose, regardless of case)
# ECO TODO: improve automatic recognition ?
# check that we have at least a column id, response and X

    if(!is.na(as.integer(object@name.group))) {
# group given as a column number
      object@name.group<-colnames(dat)[as.integer(object@name.group)]
    }
    if(is.na(object@name.group)) object@name.group<-""
    if(object@name.group=="") {
      i1<-match("id",tolower(colnames(dat)))
      if(length(i1)==0 | is.na(i1)) {
        i1<-c(grep("subject",tolower(colnames(dat)),fixed=T), grep("sujet",tolower(colnames(dat)),fixed=T))
      }
      if(length(i1)>0) {
        object@name.group<-colnames(dat)[i1[1]]
        cat("    no name for the group variable (ID) given, will use column --",object@name.group,"-- in the dataset.\n")
      }
    }
    if(is.na(match(object@name.group,colnames(dat)))) {
      cat("Can't find a column named",object@name.group,"in the data.\n")
      return("Creation of saemixData failed")
    }
    if(object@name.group=="") {
      cat("Please provide a name for the ID column.\n")
      return("Creation of saemixData failed")
    }
   i1<-as.integer(object@name.predictors[!is.na(as.integer(object@name.predictors))])
    if(length(i1)>0) { 
      object@name.predictors[!is.na(as.integer(object@name.predictors))]<- colnames(dat)[i1]
    }
    if(is.na(object@name.predictors)) object@name.predictors<-""
    if(length(object@name.predictors)==0 | (length(object@name.predictors)==1 & object@name.predictors[1]=="")) {
      i1<-c(match(c("time","temps","tps","tim","x","dose"),tolower(colnames(dat))))
      i1<-i1[!is.na(i1)]
      if(length(i1)>0) {
        object@name.predictors<-colnames(dat)[i1]
        cat("    no name for the predictor variable given, will use column(s) --",object@name.predictors,"-- in the dataset.\n")
      }
    }
    id1<-match(object@name.predictors,colnames(dat),nomatch=0)
    if(length(id1[id1==0])>0) {
      cat("    cannot find column(s) --",object@name.predictors[id1==0],"-- dropping them from the data.\n")
    }
    xnam<-object@name.predictors[id1>0]
    if(length(xnam)==0) object@name.predictors<-"" else object@name.predictors<-xnam
    if(length(xnam)==0) {
      cat("Please provide at least one predictor.\n")
      return("Creation of saemixData failed")
    }
    if(!is.na(as.integer(object@name.response))) { 
# response given as a column number
      object@name.response<-colnames(dat)[as.integer(object@name.response)]
    }
    if(object@name.response=="") {
      i1<-match("y",tolower(colnames(dat)))
      if(length(i1)==0 | is.na(i1)) { 
        i1<-c(grep("response",tolower(colnames(dat)),fixed=TRUE), match(c("resp","conc"),tolower(colnames(dat))),grep("concentration", tolower(colnames(dat)),fixed=TRUE))
        i1<-i1[!is.na(i1)]
      }
      if(length(i1)>0) {
        object@name.response<-colnames(dat)[i1[1]]
        cat("    no name for the response variable given, will use column --",object@name.response,"-- in the dataset.\n")
      }
    }
    if(is.na(object@name.response)) object@name.response<-""
    if(is.na(match(object@name.response,colnames(dat)))) {
      cat("Can't find a column named",object@name.response,"for the response column.\n")
      return("Creation of saemixData failed")
    }
    if(object@name.response=="") {
      cat("Please provide a name for the response column.\n")
      return("Creation of saemixData failed")
    }
    if(length(object@name.covariates)>0 & object@name.covariates[1]!="") {
   i1<-as.integer(object@name.covariates[!is.na(as.integer(object@name.covariates))])
      object@name.covariates[!is.na(as.integer(object@name.covariates))]<- colnames(dat)[i1]
    }
    if(nchar(object@name.group)*length(object@name.predictors)* nchar(object@name.response)<=0) {
      stop("Please check the structure of the data file and provide information concerning which columns specify the group structure (ID), the predictors (eg dose, time) and the response (eg Y, conc). See documentation for automatic recognition of column names for these elements.\n")
    }
    if(nchar(object@name.X)==0)
      object@name.X<-object@name.predictors[1]
    if(!is.na(as.integer(object@name.X))) {
      if(dim(dat)[2]<as.integer(object@name.X)) {
        cat("Attribute name.X",object@name.X,"does not correspond to a valid column in the dataset, setting the X axis for graphs to",object@name.predictors[1],".\n")
	object@name.X<-object@name.predictors[1]
      } else object@name.X<-colnames(dat)[as.integer(object@name.X)]
    } 
    if(match(object@name.X,object@name.predictors,nomatch=0)==0) {
      cat("Attribute name.X",object@name.X,"does not correspond to a valid column in the dataset, setting the X axis for graphs to",object@name.predictors[1],".\n")
      object@name.X<-object@name.predictors[1]
    }
    all.names<-c(object@name.group,object@name.predictors, object@name.response,object@name.covariates)

    dat<-dat[,all.names,drop=FALSE]
    if(class(dat)!="data.frame") dat<-as.data.frame(dat)
# Removing missing values in response & predictor columns
    if(sum(is.na(dat[,object@name.response]))>0) {
    	cat("Removing missing observations.\n")
    }
    dat<-dat[!is.na(dat[,object@name.response]),]
    for(i in object@name.predictors) {
    	if(sum(is.na(dat[,i]))>0) cat("Removing missing values for predictor",i,"\n")
    	dat<-dat[!is.na(dat[,i]),]
    }
# Saving covariates in the original format in ocov, transforming binary covariates in dat to factors
    object@ocov<-dat[,object@name.covariates,drop=FALSE]
    for(icov in object@name.covariates) {
      if(length(unique(dat[,icov]))==2) dat[,icov]<-as.integer(factor(dat[,icov]))-1
    }   
# ECO TODO: missing data in covariates kept for the moment, only excluded depending on the model
#    for(i in object@name.covariates) dat<-dat[!is.na(dat[,i]),]
    object@ntot.obs<-dim(dat)[1] # total number of observations
    id<-dat[,object@name.group]
    object@N<-length(unique(id))
    nind.obs<-tapply(id,id,length) # individual numbers of observations (1xN)
    nind.obs<-nind.obs[match(unique(id),names(nind.obs))]
    object@nind.obs<-c(nind.obs)
    dat<-cbind(index=rep(1:object@N,times=nind.obs),dat)
    object@data<-dat
    
    options(ow) # reset
    validObject(object)
    return(object)
  }
)

####################################################################################
####			SaemixData class - method to print/show data		####
####################################################################################

setMethod("print","SaemixData",
  function(x,nlines=10,...) {
    digits<-2;nsmall<-2
    cat("Object of class SaemixData\n")
    cat("    longitudinal data for use with the SAEM algorithm\n")
    cat("Dataset",x@name.data,"\n")
    st1<-paste(x@name.response," ~ ",paste(x@name.predictors,collapse=" + ")," | ", x@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    if(length(x@name.predictors)>1) {
      cat("    X variable for graphs:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
    } else  cat("    Predictor:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
    ncov<-length(x@name.covariates)
    if(ncov>0) {
      cat("    covariates:",paste(paste(x@name.covariates," (",x@units$covariates,")",sep=""),collapse=", "),"\n")
      if(length(x@ocov)>0) {
      for(icov in 1:ncov) {
      if(is.factor(x@ocov[,icov]) | length(unique(x@ocov[,icov]))==2) cat("      reference class for covariate",x@name.covariates[icov],": ",levels(as.factor(x@ocov[,icov]))[1],"\n")
      }
      }
    }
    if(FALSE) {
      cat("    Group column:",x@name.group,"\n")
      cat("    Predictors:",x@name.predictors,"\n")
      cat("    X variable for graphs:",x@name.X,paste("(",x@units$x,")",sep=""),"\n")
      cat("    Response column:",x@name.response, paste("(",x@units$y,")",sep=""),"\n")
      cat("    Covariates:",x@name.covariates,"\n")
    }
    if(length(x@data)>0) {
      if(nlines==0) return()
      cat("Dataset characteristics:\n")
      cat("    number of subjects:    ",x@N,"\n")
      if(x@N>0) {
        cat("    number of observations:",x@ntot.obs,"\n")
        cat("    average/min/max nb obs:",format(mean(x@nind.obs),digits=digits, nsmall=nsmall), " / ", min(x@nind.obs)," / ",max(x@nind.obs),"\n")
#    if(length(x@data)>0) print(x@data)
      }
      xdat<-x@data
      if(length(x@ocov)>0) xdat[,x@name.covariates]<-x@ocov
      if(nlines==(-1)) {
        cat("Data:\n")
        print(xdat)
      } else {
        cat("First",nlines,"lines of data:\n")
        nrowShow <- min (nlines , nrow(xdat ))
        print(xdat[1:nrowShow,-c(1)])
      }
    } else cat("No data.\n")
  }
)

setMethod("show","SaemixData",
  function(object) {
    cat("Object of class SaemixData\n")
    cat("    longitudinal data for use with the SAEM algorithm\n")
    cat("Dataset",object@name.data,"\n")
    st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    if(length(object@name.predictors)>1) {
      cat("    X variable for graphs:",object@name.X, paste("(",object@units$x,")",sep=""),"\n")
    }
    ncov<-length(object@name.covariates)
    if(ncov>0) {
      cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
      if(length(object@ocov)>0) {
      for(icov in 1:ncov) {
      if(is.factor(object@ocov[,icov])) cat("      reference class for covariate",object@name.covariates[icov],": ",levels(object@ocov[,icov])[1],"\n")
      }
      if(length(object@data)>0) object@data[,object@name.covariates]<-object@ocov
      }
    }
    if(length(object@data)>0) {
      cat("First lines of data:\n")
      nrowShow <- min (10 , nrow(object@data ))
      print(object@data[1:nrowShow,-c(1)])
    } else cat("No data.\n")
  }
)

# Could be print, with only head of data
setMethod("showall","SaemixData",
  function(object) {
    digits<-2;nsmall<-2
    cat("Object of class SaemixData\n")
    cat("    longitudinal data for use with the SAEM algorithm\n")
    cat("Dataset",object@name.data,"\n")
    cat("    header:",object@header,"\n")
    cat("    sep:",object@sep,"\n")
    cat("    na:",object@na,"\n")
    st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
    cat("    Structured data:",st1,"\n")
    cat("    subject identifier:    ",object@name.group,"\n")
    cat("    predictors:       ",object@name.predictors,"\n")
    cat("    response:         ",object@name.response,paste("(",object@units$y,")",sep=""),"\n")
    cat("    X variable for graphs:",object@name.X,paste("(",object@units$x,")",sep=""),"\n")
    ncov<-length(object@name.covariates)
    if(ncov>0) {
      cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
      if(length(object@ocov)>0) {
      for(icov in 1:ncov) {
      if(is.factor(object@ocov[,icov])) cat("      reference class for covariate",object@name.covariates[icov],": ",levels(object@ocov[,icov])[1],"\n")
      }
      if(length(object@data)>0) object@data[,object@name.covariates]<-object@ocov
      }
    }
    cat("Dataset characteristics:\n")
    cat("    number of subjects:    ",object@N,"\n")
    if(object@N>0) {
      cat("    number of observations:",object@ntot.obs,"\n")
      cat("    average/min/max nb obs:",format(mean(object@nind.obs),digits=digits, nsmall=nsmall), " / ", min(object@nind.obs)," / ",max(object@nind.obs),"\n")
#    if(length(object@data)>0) print(object@data)
    }
    if(length(object@data)>0) {
      cat("First lines of data:\n")
      nrowShow <- min (10 , nrow(object@data ))
      ncolShow <- min (10 , ncol(object@data))
      print(object@data[1:nrowShow,-c(1)])
    } else cat("No data.\n")
  }
)

# SaemixRepData
setMethod("show","SaemixRepData",
  function(object) {
    cat("Object of class saemixRepData\n")
    if(length(object@N)>0) {
	    cat("    replicated data used in the SAEM algorithm\n")
	    cat("    number of subjects in initial dataset",object@N,"\n")
	    cat("    number of replications",object@nrep,"\n")
	    cat("    number of subjects in replicated dataset",object@NM,"\n")
    } else cat("Empty object \n")
    } 
)

# SaemixSimData
setMethod("show","SaemixSimData",
  function(object) {
    cat("Object of class SaemixSimData\n")
    cat("    data simulated according to a non-linear mixed effect model\n")
    if(length(object@N)>0) {
    cat("Characteristics of original data\n")
    cat("    number of subjects:",object@N,"\n")
    cat("    summary of response:\n")
    print(summary(object@data[,object@name.response]))
    cat("Characteristics of simulated data\n")
    if(dim(object@datasim)[1]>0) {
      cat("    number of simulated datasets:",object@nsim,"\n")
      cat("    summary of simulated response\n")
      print(summary(object@datasim$ysim))
    } else cat("    no simulations performed yet\n")
  }}
)

####################################################################################
####				Summary method for SaemixData			####
####################################################################################

setMethod("summary","SaemixData",
  function(object, print=TRUE, ...) {
    digits<-2;nsmall<-2
    if(print) {
    	cat("Object of class SaemixData\n")
      cat("    longitudinal data for use with the SAEM algorithm\n")
      cat("Dataset",object@name.data,"\n")
      st1<-paste(object@name.response," ~ ",paste(object@name.predictors,collapse=" + ")," | ", object@name.group,sep="")
    	cat("    Structured data:",st1,"\n")
      if(length(object@name.predictors)>1) cat("    X variable for graphs:",object@name.X,paste("(",object@units$x,")",sep=""),"\n")
      if(length(object@name.covariates)>0) {
        cat("    covariates:",paste(paste(object@name.covariates," (",object@units$covariates,")",sep=""),collapse=", "),"\n")
      }
      cat("Dataset characteristics:\n")
      cat("    number of subjects:    ",object@N,"\n")
      if(object@N>0) {
        cat("    number of observations:",object@ntot.obs,"\n")
        cat("    average/min/max nb obs:",format(mean(object@nind.obs),digits=digits, nsmall=nsmall), " / ", min(object@nind.obs)," / ",max(object@nind.obs),"\n")
#    if(length(object@data)>0) print(object@data)
      }
    }
    res<-list(N=object@N,nobs=list(ntot=object@ntot.obs,nind=object@nind.obs), id=object@data[,object@name.group],x=object@data[,object@name.predictors,drop=FALSE], y=object@data[,object@name.response])
    if(length(object@name.covariates)>0) {
      res$covariates<-object@ocov
      ucov<-cbind(object@data[,object@name.group],object@ocov)
      colnames(ucov)<-object@name.group
			ucov<-ucov[match(unique(object@data$index),object@data$index),]
      res$ind.covariates<-ucov
    }
    invisible(res)
 }
)

####################################################################################
####			SaemixData class - method to plot			####
####################################################################################

saemix.data.setoptions<-function(saemix.data) {
# setting default plot options
  plot.opt<-list(
# General graphical options
    new=TRUE,				# whether a new page should be called
    ask=FALSE,				# whether the program should ask before creating a new page
    ilist=c(1:saemix.data["N"]),
    separate=FALSE,	# if TRUE, plots individual subjects (à la nlme), if FALSE plots a single plot with all the subjects
# Options for individual plots
    nmax=12,					  # maximum number of subjects
    limit=TRUE,					# limit to nmax plots
    sample=FALSE,				# if FALSE=use the (nmax) first subjects; TRUE=randomly sample (nmax) subjects from the dataset
    interactive=FALSE, 	# whether the program should prompt the user for the number of subjects to plot in the individual plots if this number exceeds nmax
# Layout and plots options
    mfrow=c(),				# page layout (if empty, defaults to the default layout for each graph type)
    main=" ",				# title
    xlab=" ",
    ylab=" ",
    col="black",
    pch=20,
    lty=1,
    lwd=1,
    xlim=c(),
    ylim=c(),
    xlog=FALSE,
    ylog=FALSE,
    type="b",
    cex=1,
    cex.axis=1,
    cex.lab=1,
    cex.main=1)
        
     if(is.null(plot.opt$name.X))
        plot.opt$name.X<-saemix.data["name.predictors"][1]
    plot.opt$xlab<-paste(plot.opt$name.X," (",plot.opt$units$x,")", sep="")
     if(length(saemix.data["name.response"])>0)
    plot.opt$ylab<-paste(saemix.data["name.response"]," (",saemix.data["units"]$y,")", sep="")
   return(plot.opt)
}

replace.data.options<-function(plot.opt,...) {
  args1<-match.call(expand.dots=TRUE)
  if(length(args1)>2) {
# Other arguments
    for(i in 3:length(args1)) {
      if(match(names(args1)[i],names(plot.opt),nomatch=0)>0)    
    plot.opt[[names(args1)[i]]]<-eval(args1[[i]]) else {
      if(!(names(args1)[i] %in% c("plot.type","individual"))) cat("Argument",names(args1)[i],"not available, check spelling.\n")
    }
   }
  }
  return(plot.opt)
}

# Plot the data, either as points or as lines grouped by x@name.group
setMethod("plot","SaemixData",
  function(x,y,...) {
    if(length(x@data)==0) {
    	cat("No data to plot.\n")
    	return()
    }
    args1<-match.call(expand.dots=TRUE)
    i1<-match("individual",names(args1))
    if(!is.na(i1)) {
    	individual<-as.logical(eval(args1[[i1]]))
    } else individual<-FALSE
    i1<-match("type",names(args1))
    if(!is.na(i1)) {
      plot.type<-as.character(args1[[i1]])
      plot.type<-plot.type[plot.type!="c"]
    } else plot.type<-c()
    if(length(plot.type)==0) plot.type<-ifelse(individual,"b","l")
    plot.opt<-saemix.data.setoptions(x)
    mainkeep<-plot.opt$main
    plot.opt$new<-TRUE
    plot.opt$xlab<-paste(x@name.X," (",x@units$x,")",sep="")
    plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
    plot.opt$type<-ifelse(individual,"b","l")
    plot.opt<-replace.data.options(plot.opt,...)
    change.main<-FALSE
    if(plot.opt$main!=mainkeep) change.main<-TRUE
    logtyp<-paste(ifelse(plot.opt$xlog,"x",""),ifelse(plot.opt$ylog,"y",""),sep="")
    if(individual) { # separate plots subject per subject
    	if(length(plot.opt$ilist)>plot.opt$nmax & plot.opt$limit) {
    		if(plot.opt$interactive) {
    			x1<-readline(prompt=paste("The number of subjects may be too large to be plotted. Should I plot only",plot.opt$nmax,"subjects ? (Y/n) \n"))
    			if(tolower(x1)=="y") {
    				plot.opt$limit<-TRUE
    				plot.opt$ilist<-plot.opt$ilist[1:plot.opt$nmax]
    				if(plot.opt$sample) plot.opt$ilist<-sort(sample(plot.opt$ilist, plot.opt$nmax)) else plot.opt$ilist<-plot.opt$ilist[1:plot.opt$nmax]
     				if(!plot.opt$ask) {
    					x1<-readline(prompt="Stop after each page of plot ? (Y/n) \n")
    			  	if(tolower(x1)=="y") plot.opt$ask<-TRUE
    				}
    			}
    		} else {
    			cat("The number of subjects is too large, I will plot only")
    			if(plot.opt$sample) cat(" the data for",plot.opt$nmax,"subjects sampled randomly;") else cat(" only the data for the first",plot.opt$nmax,"subjects;")
    			cat(" use limit=FALSE in the call to plot to force plotting all the subjects.\n")
   				if(plot.opt$sample) plot.opt$ilist<-sort(sample(plot.opt$ilist, plot.opt$nmax)) else plot.opt$ilist<-plot.opt$ilist[1:plot.opt$nmax]
    		}
    	} # end of test on length(ilist)
	    if(plot.opt$new) {
		    if(length(plot.opt$mfrow)==0) {
		    np<-length(plot.opt$ilist)
		    if(np>12) np<-12
		    n1<-round(sqrt(np))
		    n2<-ceiling(np/n1)
		    par(mfrow=c(n1,n2),ask=plot.opt$ask)
		  } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
	    }
		  xind<-x["data"][,x["name.predictors"], drop=FALSE]
		  id<-x["data"][,"index"]
		  yobs<-x["data"][,x["name.response"]]
    	for(isuj in plot.opt$ilist) {
    		if(!change.main) main<-paste("Subject",isuj) else main<-plot.opt$main
    		plot(xind[id==isuj,x@name.X],yobs[id==isuj],type=plot.type, xlab=plot.opt$xlab,ylab=plot.opt$ylab,col=plot.opt$col,pch=plot.opt$pch,log=logtyp, xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=main,cex=plot.opt$cex, cex.axis=plot.opt$cex.axis,cex.lab=plot.opt$cex.lab,lty=plot.opt$lty, lwd=plot.opt$lwd)
    	}
    } else {	# One plot for all the data
	    if(plot.opt$new) par(mfrow=c(1,1))
	      if(plot.type=="p" | plot.type=="b") {
	        plot(x@data[,x@name.X],x@data[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,pch=plot.opt$pch,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=plot.opt$main,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab) }
	      if(plot.type=="l") {
	        plot(x@data[,x@name.X],x@data[,x@name.response],xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd,type="n", log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=plot.opt$main, cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
	      }
	      if(plot.type=="l" | plot.type=="b") {
	        for(isuj in unique(x@data[,x@name.group])) {
	          lines(x@data[x@data[,x@name.group]==isuj,x@name.X], x@data[x@data[,x@name.group]==isuj,x@name.response],col=plot.opt$col, lty=plot.opt$lty,lwd=plot.opt$lwd)
	      }
			}
    }
  }
)

# Check for mirror plots
setMethod("plot","SaemixSimData",
  function(x,y,irep=-1,...) {
    args1<-match.call(expand.dots=TRUE)
    i1<-match("type",names(args1))
    if(!is.na(i1)) {
      plot.type<-as.character(args1[[i1]])
      plot.type<-plot.type[plot.type!="c"]
    } else plot.type<-"l"
    plot.opt<-saemix.data.setoptions(x)
    plot.opt$new<-TRUE
    plot.opt$plot.type<-"b"
    plot.opt$xlab<-paste(x@name.X," (",x@units$x,")",sep="")
    plot.opt$ylab<-paste(x@name.response," (",x@units$y,")",sep="")
    plot.opt<-replace.data.options(plot.opt,...)
    logtyp<-paste(ifelse(plot.opt$xlog,"x",""),ifelse(plot.opt$ylog,"y",""),sep="")
    if(length(x@sim.y)==0) cat("No simulated data.\n") else {
      if(irep<0) irep<-sample(unique(x@sim.rep),1)
      tit<-paste("Mirror plot (replication ",irep,")",sep="")
      tab<-data.frame(id=x@data[,name.group],x=x@data[,x@name.X], y=x@datasim$ysim[x@datasim$irep==irep])
      if(plot.type=="p" | plot.type=="b") {
        plot(tab[,"x"],tab[,"y"],xlab=plot.opt$xlab, ylab=plot.opt$ylab, col=plot.opt$col,pch=plot.opt$pch,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=tit,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab) }
      if(plot.type=="l") {
        plot(tab[,"x"],tab[,"y"],type="n",xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd,type="n", log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=tit, cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
      }
      if(plot.type=="l" | plot.type=="b") {
        for(isuj in unique(tab[,"id"])) {
          lines(tab[tab[,"id"]==isuj,"x"],tab[tab[,"id"]==isuj,"y"])
      }
      }
    }
  }
)

####################################################################################
####		Creating an object of SaemixData class - User-level function	####
####################################################################################

saemixData<-function(name.data,header,sep,na,name.group,name.predictors, name.response,name.X,name.covariates=c(), units=list(x="",y="",covariates=c())) {
# setting proper types for the SaemixData class
  if(missing(name.data)) {
    cat("Error in saemixData: please provide the name of the datafile or dataframe (between quotes)\n")
    return("Creation of SaemixData failed")
  }
  if(is.data.frame(name.data)) name.data<-deparse(substitute(name.data))
  if(missing(header)) header<-TRUE
  if(missing(sep)) sep<-""
  if(missing(na)) na<-"NA" else {na<-as.character(na);na[is.na(na)]<-"NA"}
  if(missing(name.group)) name.group<-"" else name.group<-as.character(name.group)
  if(missing(name.predictors)) name.predictors<-"" else name.predictors<-as.character(name.predictors)
  if(missing(name.response)) name.response<-"" else  name.response<-as.character(name.response)
  if(missing(name.X)) name.X<-"" else name.X<-as.character(name.X)
  name.covariates<-as.character(name.covariates)
  x<-new(Class="SaemixData",name.data=name.data,header=header,sep=sep,na=na, name.group=name.group,name.predictors=name.predictors,name.X=name.X, name.response=name.response,name.covariates=name.covariates,units=units)
#  showall(x)
  x1<-read.saemixData(x)
  if(class(x1)!="character") cat("\n\nThe following SaemixData object was successfully created:\n\n")
  print(x1,nlines=0)
  return(x1)
}

####################################################################################
####				saemixObject class - S3 methods			####
####################################################################################

subset.SaemixData<-function (x, subset, ...) {
    if (missing(subset)) 
        return(x)
    else {
        e <- substitute(subset)
	      xdat<-x["data"]
        r <- eval(e, xdat, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    x1<-x
    x1["data"]<-x["data"][r,,drop=FALSE]
    if(length(x["yorig"])>0) x1["yorig"]<-x["yorig"][r]
    if(length(x["ocov"])>0) x1["ocov"]<-x["ocov"][r,,drop=FALSE]
    id<-x1["data"][,x1["name.group"]]
    x1["N"]<-length(unique(id))
    nind.obs<-tapply(id,id,length) # individual numbers of observations (1xN)
    nind.obs<-c(nind.obs[match(unique(id),names(nind.obs))])
    x1["nind.obs"]<-nind.obs
    x1["ntot.obs"]<-length(id)
    x1["data"]$index<-rep(1:x1["N"],times=nind.obs)
    return(x1)
}

####################################################################################
