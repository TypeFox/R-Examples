#' @importFrom XML xmlParse
#' @importFrom XML getNodeSet
#' @importFrom XML xmlAttrs
NULL

processing_monolix  <- function(project,model=NULL,treatment=NULL,parameter=NULL,
                                output=NULL,group=NULL,r.data=TRUE,fim=NULL)
{
  ### processing_monolix
  #     takes a monolix project and extract information from
  #     mlxtran file such as model, admin, param, output, if
  #     they missed in the input parameters. 
  #     Theses informations are read in files.
  
  ##************************************************************************
  #       XML FILENUL
  #*************************************************************************
  infoProject <- getInfoXml(project)
  n.output <- length(infoProject$output)
  param <- parameter
  
  #   infoProject$resultFolder <- "project_simul"
  ##************************************************************************
  #       DATA FILE
  #**************************************************************************
  
  if (r.data==TRUE){
    datas <- readDatamlx(infoProject=infoProject)
    
    y.attr <- sapply(datas,attr,"type")
    j.long <- which(y.attr=="longitudinal")
    # dobs <- datas$observation
    # if (!is.null(names(dobs)))
    #   dobs <- list(dobs)
    datas$observation=list()
    # y <- list()
    for (iy in (1:length(j.long))){
      yi <- datas[[j.long[iy]]]
      niy <- names(yi)
      yk <- list(ylabel="observation", colNames=niy, name=niy[length(niy)], value=yi )
      datas$observation[[iy]] <- yk
    }    
    if (!is.null(datas$treatment)){
      ntr <- names(datas$treatment)
      datas$sources <- list(ylabel="sources", colNames=ntr, name="doseRegimen", value=datas$treatment )
      datas$treatment <- NULL
    }
    
    #     if (is.character(param))  {
    if (length(param)==1  && any(sapply(param,is.character))) {
      file = file.path(infoProject$resultFolder,'indiv_parameters.txt') 
      datas$parameter = readIndEstimate(file,param[which(sapply(param,is.character))])
      #        datas$parameter = readIndEstimate(file,param)
      iop_indiv=1
    }else{
      iop_indiv=0
    }
    #   data$id <- data.frame(NewId=seq(1:N),OriId=new.id)
    
    datas$id <- data.frame(newId=seq(1:datas$N),oriId=datas$id)
    if  (!is.null(param)) 
    {
      for (k in (1:length(param)))
      {
        if (isfield(param[[k]],"id"))
        {
          did <- unique(param[[k]]$id)
          datas$id <- data.frame(newId=did,oriId=did)
        }
      }
    }
    ##*********************************************************************
    #       treatment (TREATMENT)
    #**********************************************************************
    if (is.null(treatment)){
      if (is.null(datas$sources$value)){ 
        treatment = datas$sources
      } else{
        treatment = data.frame(datas$sources$value)
        names(treatment) <- datas$sources$colNames 
      }
    }
  }else{
    datas <- NULL
    iop_indiv <- 0
  }
  
  if (identical(fim,"needed")){
    if  (file.exists(file.path(infoProject$resultFolder,'correlationEstimates_sa.txt')))
      fim <- 'sa'
    else if (file.exists(file.path(infoProject$resultFolder,'correlationEstimates_lin.txt')))
      fim <- "lin"
    else
      fim <- NULL
  }
  
  ##************************************************************************
  #       PARAMETERS
  #**************************************************************************
  r = readPopEstimate(file.path(infoProject$resultFolder,'estimates.txt'),fim);
  pop_param <- r[[1]]
  paramp <- list(pop_param,datas$covariate,datas$parameter)
  
  if  (!is.null(param)) 
    paramp <- mergeDataFrame(paramp, param)
  
  ##************************************************************************
  #       FIM
  #**************************************************************************
  fim <- r[[3]]
  if (!is.null(fim)){
    pop_se <- r[[2]]
    if (fim=="sa")
      f.mat = readFIM(file.path(infoProject$resultFolder,'correlationEstimates_sa.txt'))
    else
      f.mat = readFIM(file.path(infoProject$resultFolder,'correlationEstimates_lin.txt'))
    if0 <- which(names(f.mat) %in% names(unlist(param)))
    f.mat[if0,] <- f.mat[,if0] <- NaN
    pop_se[if0] <- 0
    fim <- list(mat=f.mat,se=pop_se)
  }  
  
  #   M1 = readFIM(file.path(infoProject$resultFolder,'correlationEstimates_sa.txt'))
  #   M2 = readFIM(file.path(infoProject$resultFolder,'fimTransPop_sa.txt'))
  #   M3 = readFIM(file.path(infoProject$resultFolder,'jacobian.txt'))
  #   
  #   M1 <- as.matrix(M1)
  #   M2 <- as.matrix(M2)
  #   M3 <- as.matrix(M3)
  #   se=r[[2]]
  #   C <- diag(se)%*%M1%*%diag(se)
  #   M3%*%M2%*%M3
  
  ##************************************************************************
  #       OUTPUT 
  #**************************************************************************
  
  outputp = datas$observation
  if (is.null(output)){
    output = outputp
  }else{
    output <- formato(output)
    output <- mergeArg(outputp,output)
  }
  #   for (k in (1:length(output))){
  #     outk <- output[[k]]
  #     if (!is.null(outk$colNames)){
  #       n.outk <- outk$colNames
  #       output[[k]] <- list(name=n.outk[!(n.outk %in% c("id","time"))], 
  #                           time = outk$value[,c("id","time")])
  #     }
  #   }
  
  ##************************************************************************
  #       MODEL
  #**************************************************************************
  if (is.null(model))
  {
    # generate model from mlxtran file  
    mlxtranfile = file_path_sans_ext(basename(project))
    mlxtranpath <- dirname(project)
    model = file.path(mlxtranpath,paste0(mlxtranfile,"_model.txt"))
    session<-Sys.getenv("session.simulx")
    zz=file.path(session,'lib','lixoftLanguageTranslator')
    str=paste0('"',zz,'" --from=mlxproject --to=mlxtran')  
    str=paste0(str,' --output-file=',model,' --input-file=',project,' --option=with-observation-model') 
    system(str, wait=T)
    pathLib = Sys.getenv("LIXOFT_HOME")
    if (iop_indiv==1) {      
      # create a submodel file of model_file corresponding to the specified sections specified 
      sections       = c("LONGITUDINAL")    
      myparseModel(model, sections, model )
    }
  }
  #**************************************************************************
  #   test.colNames <- testC(list(treatment,param,output))
  #   if ((test.colNames==TRUE) | (is.null(group[[1]]$size))) {
  #     gr=group
  #   }else{
  #     gr <- list(size=c(group[[1]]$size, 1) , level=c("individual","longitudinal"))
  #   }
  gr <- group
  
  #  if (is.null(datas$regressor))
  #    ans = list(model=model, treatment=treatment, param=paramp, output=output, group=gr, id=datas$id)
  #   else
  
  
  #change regressor names and use these defined in the model in the same order 
  if (!is.null(datas$regressor))
  {
    namesReg<-names(datas$regressor)
    nbModelreg<-0
    lines <- readLines(model)
    regressorLine <-  grep('regressor', lines, fixed=TRUE, value=TRUE)
    if(length(regressorLine))
    {
      regModelNames<-c()
      for(line in seq(1:length(regressorLine)))
      {
        regModelNamesTable<-strsplit(regressorLine[line],"[\\{ \\} , ]")[[1]]
        for( i in seq(1:length(regModelNamesTable))){
          regi <- regModelNamesTable[i]
          if(!identical(regi,"")&&!length(grep("=",regi,fixed=TRUE,value=TRUE))
             &&!length(grep("regressor",regi,fixed=TRUE,value=TRUE))
             &&!length(grep("use",regi,fixed=TRUE,value=TRUE))
          ){
            regModelNames<-c(regModelNames,regi)
            nbModelreg = nbModelreg +1
          }
        }
      }
      nbregOrig<-0
      iregModel <-1
      for( i in seq(1:length(namesReg))){
        if(!identical(tolower(namesReg[i]),"id") &&
           !identical(tolower(namesReg[i]),"time")){
          namesReg[i] <- regModelNames[iregModel]
          
          iregModel <-iregModel +1 
          nbregOrig <- nbregOrig +1
        }
      }
      if(nbregOrig +1 !=  iregModel)
      {
        stop("inconsistent number of regressor between model and dregressor Field")
      }
      
      names(datas$regressor)<-namesReg
    }
    
    #---------------------------------------------------------------
  }
  
  ##set correct name of error model in parameter,  it can change  in the V2 model
  paramp[[1]]<-setErrorModelName(paramp[[1]],model)
  ##initialize latent covariates defined in the model but not used,  in parameter
  paramp[[1]]<-initLatentCov(paramp[[1]],model)
  gr    <- mklist(gr)
  #   parameter <- mklist(paramp)
  parameter <- paramp
  treatment <- mklist(treatment)
  regressor <- mklist(datas$regressor)
  occ  <- mklist(datas$occ)
  output    <- mklist(output)
  
  ans = list(model=model, 
             treatment=treatment, 
             param=parameter, 
             output=output, 
             group=gr,
             regressor=regressor, 
             id=datas$id,
             occasion=occ,
             fim=fim,
             infoParam=infoProject$parameter)
  
  return(ans)
}

getInfoXml  <- function (project)
{
  ### getInfoXml
  #
  #   getInfoXml(project)
  #     return infoProject which contains informations about the given mlxtran project
  #
  
  infoProject = list(datafile=NULL, dataformat=NULL, dataheader=NULL, output=NULL, resultFolder=NULL, mlxtranpath=NULL );
  # get path and name of monolix project
  mlxtranpath      = dirname(project);
  mlxtranpathfile = file_path_sans_ext(project)
  mlxtranfile = file_path_sans_ext(basename(project))
  infoProject$mlxtranpath = mlxtranpath
  if(file_ext(project) == "mlxtran") {
    #  project<-mlxProject2xml(project)
    session<-Sys.getenv("session.simulx")
    xmlfile <- file.path(mlxtranpath,paste0(mlxtranfile,"_tr.xmlx"))
    zz=file.path(session,'lib','lixoftLanguageTranslator')
    str=paste0('"',zz,'" --from=mlxproject --to=xmlx')  
    str=paste0(str,' --output-file=',xmlfile,' --input-file=',project,' --option=with-observation-model') 
    system(str, wait=T)
  } else {
    xmlfile <- project
  }
  
  infoResultFolder         = myparseXML(xmlfile, mlxtranpath, "resultFolder")
  infoProject$resultFolder = infoResultFolder[[1]]$uri
  ##************************************************************************
  #       GET DATA INFO
  #*************************************************************************
  #  get data format and data header used in the current project
  #   Exemple : 
  #
  #   infoProject = 
  #           datafile         : './warfarin_data.txt'
  #           dataformat       : '\t'
  #           dataheader       : {'ID'  'TIME'  'AMT'  'Y'  'YTYPE'  'COV'  'IGNORE'  'IGNORE'}
  #
  #
  infoData                = myparseXML(xmlfile, mlxtranpath, "data")
  infoProject$datafile    = infoData[[1]]$uri
  infoProject$dataformat  = infoData[[1]]$columnDelimiter
  infoProject$dataheader  = infoData[[1]]$headers
  infoOutput              = myparseXML(xmlfile, mlxtranpath, 'observationModel')
  
  for (k in 1:length(infoOutput))
    infoProject$output[[k]] = infoOutput[[k]]$name;
  
  infoParam <- myparseXML(xmlfile, mlxtranpath, "parameter")
  #   info.length <- unlist(lapply(infoParam,length))
  #   infoParam <- infoParam[info.length==2]
  i.trans <- which(!unlist(lapply(lapply(infoParam, "[[", "transformation"),is.null)))
  infoParam <- infoParam[i.trans]
  p.names <- do.call("rbind", lapply(infoParam, "[[", "name"))[,1]
  p.trans <- do.call("rbind", lapply(infoParam, "[[", "transformation"))[,1]
  infoProject$parameter <- list(name=p.names, trans=p.trans)
  
  if(file_ext(project) == "mlxtran")
    unlink(xmlfile, recursive=T)
  
  return(infoProject)
}

myparseXML  <- function (filename, mlxtranpath, node)
{
  ### myparseXML
  #
  #  myparseXML(filename, node)
  #     return information of mentioned node contained in the xmlx file (filename)
  #
  
  tree    = xmlParse(filename)
  set     = getNodeSet(tree,paste0("//", node))
  tmp=list(name=NULL)
  ans=list()
  for (i in 1 : length(set)) {
    attributs      = xmlAttrs(set[[i]])
    namesAttributs = names(attributs)
    tmp['name']= node
    for (j in 1 : length(namesAttributs)) {      
      tmp[namesAttributs[[j]]]=attributs[[j]]
      # replace '%MLXPROJECT%' by the symbol of current folder "."
      if (namesAttributs[[j]] == "uri")
        tmp[namesAttributs[[j]]]=sub("%MLXPROJECT%", mlxtranpath, tmp[namesAttributs[[j]]]) 
      # repalce '\\t' by "tab"
      if (namesAttributs[[j]] == "columnDelimiter") {
        if (tmp[namesAttributs[[j]]] == '\\t' )
          tmp[namesAttributs[[j]]]= "tab"
      }
    } 
    ans= c(ans, list(tmp))
  }
  return(ans)
}

##
readPopEstimate  <-  function(filename, fim=NULL) {
  if (file.exists(filename)) {
    data        = read.table(filename, header = TRUE, sep=";")
    if (ncol(data)==1)
      data        = read.table(filename, header = TRUE, sep="\t")
    
    name        = as.character(data[[1]])
    name        = sub(" +", "", name)
    name        = sub(" +$", "", name)
    
    #ic <- grep("corr_",name)
    #name[ic] <- sub("corr_","r_",name[ic])
    param <- as.numeric(as.character(data[['parameter']]))
    names(param) <- name
    
    if (!is.null(fim)){
      if (fim=='lin'){
        se <- as.numeric(as.character(data[['s.e._lin']]))
        if (length(se)==0)
          stop("Fisher Information matrix estimated by linearization is not available")
        names(se) <- name
      } else if (fim=='sa'){
        se <- as.numeric(as.character(data[['s.e._sa']]))
        if (length(se)==0)
          stop("Fisher Information matrix estimated by stochastic approximation is not available")
        names(se) <- name
      }
    } else {
      se <- NULL
    }
    
    return(list(param,se,fim))
    
  } else
    stop(paste("file : ",filename, " does not exist" ))
}

##
readIndEstimate  <-  function(filename, estim=NULL) {
  if (file.exists(filename)) {
    data         = read.table(filename,  header = TRUE)
    # data[[1]]    = c(1: length(data[[1]]))
    header       = names(data)
    idx          = grep(paste0("_", estim), header)
    name         = header[idx]
    name         = gsub(paste0("_", estim),"", name)
    header       = c( 'id', name)
    header       = gsub("\\.","",header)
    value        = as.matrix(data[, c(1, idx)])
    param <- data.frame(value)
    names(param) <- header
    return(param)
  } else
    stop(paste("file : ",filename, " does not exist" ))
}

##
readFIM  <-  function(filename){
  if (file.exists(filename)) {
    data        = read.table(filename, header = FALSE, sep=";")
    name        = as.character(data[[1]])
    name        = sub(" +", "", name)
    name        = sub(" +$", "", name)
    ic <- grep("corr_",name)
    name[ic] <- sub("corr_","r_",name[ic])
    data[[1]] <- NULL
    row.names(data) <- name   
    names(data) <- name   
    return(data)
  } else
    stop(paste("file : ",filename, " does not exist" ))
}

mergeDataFrame  <- function(p1,p2) {
  if  (!is.null(names(p2))) 
    p2 <- list(p2)
  for (k in (1:length(p2))){
    paramk <- p2[[k]]  
    if (is.list(paramk)){
      if (is.null(paramk$id)) {
        if (is.vector(paramk$value)){
          p.temp <- as.vector(paramk$value)
          names(p.temp) <- paramk$name
        } else {
          p.temp <- data.frame(paramk$value)
          names(p.temp) <- paramk$colNames
        }
      }else{
        #         p.temp <- paramk$value
        #         names(p.temp) <- paramk$name
        p.temp <- paramk
      }
      p2[[k]] <- p.temp
    }
  } 
  n1 = length(p1)
  i1 <- which(!unlist(lapply(p1, is.null)))
  n2 = length(p2)
  #   p  = p1
  #   np = length(p)
  for (i in 1:n2) {
    p2i=p2[[i]]
    testi  = 0
    namei2 = names(p2i)
    if ("id" %in% namei2){
      for (j in i1) {
        p1j = p1[[j]]
        i12 <- which(namei2 %in% names(p1j))
        if (length(i12)>1){
          p1j=p2i
          p1[[j]] = p1j
        }
      }      
    }else{
      for (j in i1) {
        p1j = p1[[j]]
        i12 <- which(namei2 %in% names(p1j))
        namei2=namei2[namei2!="id"]
        if (!is.null(namei2[i12]))
          p1j[namei2[i12]]=p2i[namei2[i12]]
        p1[[j]] = p1j
      }
    }
  }
  return(p1)
}

mergeArg  <- function(p1,p2)
{
  #mergeArg  
  #
  #    mergeArg(p1,p2)
  
  if (!(is.list(p1[[1]])))
    p1 = list(p1)
  if (!(is.list(p2[[1]])))
    p2 = list(p2)
  
  n1 = length(p1)
  n2 = length(p2)
  p  = p1
  np = length(p)
  for (i in 1:n2) {
    p2i=p2[[i]]
    #     if (is.null(p2i$time) || length(p2i$time)>0 ) {
    #     if (isfield(p2i,'colNames')) {
    if (!is.null(p2i$colNames)) {
      testi  = 0
      namei2 = p2i$name
      for (j in 1:n1) {
        p1j = p1[[j]];
        #         if (isfield(p1j,'colNames')) {
        if (!is.null(p1j$colNames)) {
          if (namei2==p1j$name) {
            p[[j]] = p2i
            testi  = 1
          }
        }else{
          ifs = match(namei2,p1j$name)
          if(!is.na(ifs)) {
            p[[j]]$name  = p[[j]]$name[-ifs]
            p[[j]]$value = p[[j]]$value[-ifs]
            np                = np+1
            p[[np]]           = p2i
            testi             = 1
          }
        }
      }
      if (testi==0){
        np      = np+1
        p[[np]] = p2i
      }
    }else{ 
      if (length(p2i$name)>0) {
        for (k in 1:length(p2i$name)) {
          namek2 = p2i$name[k]
          testk  = 0
          for (j in 1:n1) {
            p1i = p1[[j]]
            #             if (isfield(p1i,'colNames')) {
            if (!is.null(p1i$colNames)) {
              if (namek2==p1i$name) {
                p[[j]] =list(name= list(namek2))
                #                 if (isfield(p2i,'value'))
                if (!is.null(p2i$value))
                  p[[j]]$value=p2i$value[k]
                if ("time" %in% names(p2i))
                  #   if( isfield(p2i,'time'))
                  p[[j]]$time=p2i$time
                testk = 1
              }
            }else{
              ifs=match(namek2,p1i$name)
              if (length(ifs)>0) {
                p[[j]]$value[ifs] = p2i$value[k]
                testk             = 1
              }
            }
          }
          if (testk==0) {
            np = np+1;
            p[[np]] =list(name= list(namek2))
            #             if (isfield(p2i,'value'))
            if (!is.null(p2i$value))
              p[[np]]$value=p2i$value[k]
            #             if (isfield(p2i,'time'))
            if (!is.null(p2i$time))
              p[[np]]$time=p2i$time
          }
        }
      }
    }     
    #     }
  }
  ik <- NULL
  for (k in (1:length(p))){
    pk <- p[[k]]
    if (!is.null(pk$time) && pk$time=="none")
      ik <- c(ik,k)
  }
  p[ik] <- NULL
  return(p)
}

myparseModel  <-  function(model_file, sections, submodel_file)
{
  #   myparseModel create a submodel_file corresponding to the specified sections of model_file
  #
  #    myparseModel(model_file, sections, submodel_file)
  #       myparseModel create a submodel file corresponding to the specified sections of model_file
  #
  #       The specified sections of a model (model_file) are written in a new file (submodel_file) 
  #       In case of multiple sections, each section are concatenated in a single model. 
  #
  #       sections :  a list of string containing the name of the sections we want to write into a new file. 
  #                   could be "POPULATION",  "OBSERVATION", "INDIVIDUAL", "COVARIATE"
  #
  #   Examples
  #   --------
  #       model_file    = "home/model.txt"
  #       sections      =  c("COVARIATE", "OBSERVATION")
  #       submodel_file = "home/submodel.txt"
  #       myparseModel(model_file, sections, submodel_file)
  #
  
  #splitModel a model file_model into multiple terms corresponding to the specified sections 
  terms = splitModel( model_file,sections) 
  
  str= ""
  for (i in 1 :  length( terms))
    str= c(str, terms[[i]]$model)
  
  write(str,submodel_file)
}

splitModel  <-  function(file_model, sections)
{
  #   splitModel split a model file_model into multiple terms corresponding to the specified sections 
  #
  #   terms = splitModel(file_model, sections)
  #       splitModel split a model file_model into multiple terms corresponding to the specified sections 
  #
  #       The extracted terms are returned in form of a list of strings. 
  #       Each element of terms have two fileds :  name and model  
  #       name corresponds to the name of the section contained in model. 
  #
  #       sections :  a list of string containing the name of the sections we want to use to split the model. 
  #                   could be "POPULATION", "COVARIATE","INDIVIDUAL", "LONGITUDINAL"
  #   Examples
  #   --------
  #       file_model  = "home/model.txt"
  #       sections    =  c("COVARIATE", "LONGITUDINAL")
  #       terms       = splitModel(file_model, sections)
  #
  #      > terms[[1]]$name
  #          "COVARIATE"
  #      > terms[[1]]$model
  #          chr [1:9]
  #      > terms[[2]]$name
  #          "LONGITUDINAL"
  #      > terms[[2]]$model
  #          chr [1:20]
  #
  if (file.exists(file_model))
  {
    terms         = list()
    length(terms) = length(sections)
    for (i in 1 : length(sections))
    {
      sections_i = sections[[i]]
      con        = file(file_model, open = "r")
      lines      = readLines(con)
      close(con)
      
      idx_sections   = grep("[",lines, fixed=TRUE)
      idx_sections   = c(idx_sections,length(lines))
      idx            = grep(sections_i,lines, fixed=TRUE)
      fin_sections   = idx_sections[idx_sections>idx]
      model_temp     = ""
      while (idx < fin_sections[[1]]) {
        model_temp = c(model_temp, lines[idx])
        idx        = idx +1
      }
      terms[[i]]$name = sections_i
      terms[[i]]$model= model_temp
    }
    return(terms)
    
  }else
  {
    stop(paste("file : ",file_model, " does not exist" ))
  }
  
}

getInputSection  <-  function(model_file, section)
{
  #   getInputSection extract the input list corresponding to the specified section of model_file
  #
  #    inputList = getInputSection(model_file, section)
  #       getInputSection extract the input list corresponding to the specified section of model_file
  #
  #
  #       section :  a string containing the name of the section we want to get the input list. 
  #                   could be "POPULATION",  "LONGITUDINAL", "INDIVIDUAL", "COVARIATE"
  #
  #   Examples
  #   --------
  #       model_file    = "home/model.txt"
  #       section      =  c("INDIVIDUAL")
  #       getInputSection(model_file, section)
  #       inputList
  
  #split a model file_model into multiple terms corresponding to the specified sections 
  subsection  = splitModel( model_file,section) 
  temp        = subsection[[1]]
  # extract the input list of the subsection 
  idx         = grep("input",temp$model, fixed=TRUE)
  inputList   = temp$model[[idx]] 
  
  #   Exemple : 
  #   "input = {V_pop, Cl_pop, omega_V, omega_Cl, beta_V, weight}"
  #
  
  chaine1      = strsplit(inputList,"\\{")
  lc1 <- length(chaine1[[1]])
  chaine1      = chaine1[[1]][lc1]
  chaine2      = strsplit(chaine1,"\\}")
  chaine2      = chaine2[[1]][1]
  #   split en fonction de ","  dans chaine2 = "V_pop, Cl_pop, omega_V, omega_Cl, beta_V, weight"
  chaine3       = strsplit(chaine2,"\\,")
  chaine3       = chaine3[[1]]
  chaine3        = sub(" +", "", chaine3)
  chaine3        = sub(" +$", "", chaine3)
  inputList     = chaine3
  
  
  return(inputList)
}

sectionsModel  <-  function(file_model)
{
  sections <- c("[POPULATION]", "[COVARIATE]", "[INDIVIDUAL]",  "[LONGITUDINAL]")
  if (file.exists(file_model)) {
    terms   <- NULL
    con     <- file(file_model, open = "r")
    lines   <- readLines(con)
    close(con)
    for (i in 1 : length(sections)) {
      sections_i <- sections[[i]]      
      idx <- grep(sections_i,lines, fixed=TRUE)
      if (length(idx)>0)
        terms <- c(terms,sections_i)
    }
    return(terms)
  }else
    stop(paste("file : ",file_model, " does not exist" ))
  
}


#-------------------------------------------
formato <- function(out)
{
  if (!is.null(names(out))){  
    out=list(out) 
  }
  output <- vector("list",length(out))
  for (k in seq(1,length(out))){
    outk <- out[[k]]
    if (!isfield(outk,"name"))
      outk <- list(name=outk)
    output[[k]] <- outk
  }
  return(output)
}

#----------------------------------
testC  <- function(x)
{
  testC <- FALSE
  d <- length(x)
  for (k in seq(1,d)) {
    xk <- x[[k]]
    if (length(xk)>0){
      if (!is.null(names(xk))){
        if (any("id" %in% names(xk)))
          testC <- TRUE
      }else{
        dk <- length(xk)
        for (j in seq(1,dk)) {
          if (any( "colNames" %in% names(xk[[j]]) ))
            testC <- TRUE
        }
      }
    }
  }
  return(testC)
}


setErrorModelName<- function(param,model)
{
  ## set correct the name of error model in param, it can change  in the V2 model
  ## if user's parameter is the same as error model name  
  
  modelread <- readLines(model)
  errorline<-grep("errorModel",modelread)
  errorused<-NULL
  for(i in seq(1:length(errorline)))
  {    
    comment<-";"
    line<-strsplit(modelread[errorline[i]],comment)[[1]][1]    
    if(length(line))
    {       
      testerr<-strsplit(line,"errorModel")
      if(length(testerr[[1]])==2)
      {          
        errorsub<-strsplit(testerr[[1]][2],'[/(/)]',perl=TRUE)
        errorargs<-strsplit(errorsub[[1]][2],',')
        for(ee in seq(1:length(errorargs[[1]]))){
          errorused<-c(errorused,(errorargs[[1]][ee]))
        }        
      }       
    }
  }
  
  ##replace names without _ in param  
  replaced=FALSE
  endFlag<-"_"
  paramRead <- names(param)
  if(length(errorused))
  {    
    for(i in seq(1:length(errorused)))
    {
      erri<-errorused[i]
      erriChars <- strsplit(erri,"")[[1]]
      lastChar <-  erriChars[length(erriChars)]
      if(identical(lastChar,endFlag))
      {
        errifirstchars<-strsplit(erri,'(.$)',perl=TRUE)[[1]]
        erriline<-paste0("^",erri,"$")
        errparamf<-grep(erriline,paramRead,perl=TRUE)
        if(!length(errparamf))
        {
          for(i in seq(1:length(paramRead)))
          { 
            paramname<- NULL
            linesplitted<-strsplit(paramRead[i],'\\s',perl=TRUE)
            is<-1
            if(length(linesplitted[[1]]))
            {
              for(j in seq(1:length(linesplitted[[1]])))
              {
                if(!identical(linesplitted[[1]][j],""))
                {
                  paramname<-linesplitted[[1]][j]
                  is<-j
                  break
                }
              }
            }            
            if(identical(paramname,errifirstchars))
            {            
              paramRead[i]<-sub(linesplitted[[1]][is],erri,paramRead[i]) 
              replaced=TRUE 
              break
            }          
          }
        }
      }
    }
  }
  names(param)<-paramRead
  return(param)
}

initLatentCov<- function(param,model)
{
  ## initialize latent covariates defined in the model but nit used, 
  ## thus not present in estimates.txt
  
  modelRead <- readLines(model)
  inputLinesNum  <- plcatLine<-grep("input",modelRead)  
  inputRead <- modelRead[inputLinesNum]
  latentCovPrefix <-"plcat"
  plcatLine<-grep(latentCovPrefix,inputRead)
  plcatUsed <-NULL
  if(length(plcatLine))
  {
    for(i in seq(1:length(plcatLine)))
    {    
      comment<-";"
      line<-strsplit(inputRead[plcatLine[i]],comment)[[1]][1]    
      if(length(line)){
        lineSplit<-strsplit(line,'[/{/}," "]',perl=TRUE)[[1]]
        iPlcat<-grep(paste0("^",latentCovPrefix),lineSplit)
        for(ee in seq(1:length(iPlcat)))
        {
          plcatUsed<-c(plcatUsed,(lineSplit[iPlcat[ee]]))            
        }
      }
    }
    
    ## add plcatused in param with 0 as value    
    if(length(plcatUsed))
    {   
      namesParam<-names(param)
      plcatInParam<-NULL
      for(iused in seq(1:length(plcatUsed)))
      {
        for(iparam in seq(1:length(namesParam)))
        {
          if(identical(plcatUsed[iused],namesParam[iparam]))
          {
            plcatInParam <-c(plcatInParam,-iused)
            break
          }          
        }
      }
      if(!is.null(plcatInParam))
      {
        plcatUsed <-plcatUsed[plcatInParam]
      }
      
      if(length(plcatUsed))
      {
        for(i in seq(1:length(plcatUsed)))
        {
          namesParam <-c(namesParam,plcatUsed[i])
          param<-c(param,0)
        }
        names(param) <- namesParam
      }
    }    
  }
  return(param)
}
