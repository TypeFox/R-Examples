###################################################
### resolveVariableArguments

.resolveVariableArguments <- function(args) {
  # if the first argument is an unnamed list then just use this list
  if ( (length(args) == 1L) &&
         is.list(args[[1L]])  &&
         (is.null(names(args)) || (names(args)[[1L]] == "")) )  {
    return (args[[1L]])
  } else {
    return (args)
  }
}

###################################################
###  ManipulateClass

.Manipulate <- setRefClass("Manipulate", fields=list(.code="ANY",.controls="list"))


###################################################
###  Manipulate: validate_controls

.Manipulate$methods(
  validate_controls = function() {
  "Validate that controls are specified properly"
  ## validate that all controls have unique names
  controlNames <- names(.controls)
  duplicatedIndex <- anyDuplicated(controlNames)
  if (duplicatedIndex > 0)
    stop(paste("duplicated control name:", controlNames[[duplicatedIndex]]))
  
  ## iterate over the names and controls, adding the default values to the env
  for (name in names(.controls)) {
    ## check the name
    if (name == "")
      stop("all controls passed to manipulate must be named")
    ## confirm that this is in fact a control
    if(!is(.controls[[name]], "ManipulateControls"))
      stop(paste("argument", name, "is not a control"))
    ## default label is control name
    if(length(.controls[[name]]$label) == 0) 
      .controls[[name]]$label <<- name
  }
})


###################################################
### Manipulate_change_handler

.Manipulate$methods(
  get_values = function() {
    "Get widget values as list"
    sapply(.controls, function(control) control$get_value(), simplify=FALSE)
  },
  change_handler = function(...) {
    "Evaluate code with current values"
    values <- get_values()
    result <- withVisible(eval(.code, envir = values))
    if (result$visible) {
      eval(print(result$value))
    }
  }
)


###################################################
.Manipulate$methods(  
  execute=function(parent=NULL) {
    "Make the GUI"
    if(!is.null(names(mrv$mwindow))){
      mrv$mwindow$destroy() 
    }
    mrv$mwindow <- gtkWindow(show=FALSE)
    mrv$mwindow$setTitle("Choose Prior Density")
    ## Set up graphic device
    vpaned <- gtkVPaned()
    mrv$mwindow$add(vpaned)
    device <- gtkDrawingArea()
    device$setSizeRequest(480, 480)
    asCairoDevice(device)
    vpaned$add(device)
    ## Controls frame
    frame <- gtkFrame("Controls")
    control_table <- gtkTableNew()
    control_table$setHomogeneous(FALSE)
    control_table['column-spacing'] <- 10
    ## insert horizontal strut
    control_table$attach(strut <- gtkHBox(), 1,2,0,1, xoptions="", yoptions="shrink")
    strut$setSizeRequest(75, -1)
    frame$add(control_table)
    
    
    button2 <- gtkButton("OK")
    button2['image'] <- gtkImage(stock.id = "gtk-apply", size = "button")
    gSignalConnect(button2, "clicked", function(button) {
      para_value <- get_values()
      if(!is.null(parent)){
        if(length(parent) > length(para_value)){
          print("error: one entry should has only one value from slider")
        }
        for(i in 1:length(parent)){
          gtkEntrySetText(parent[[i]],para_value[i])
        }
      }
      mrv$mwindow$destroy()
      mrv$mwindow = NULL
    })
    vbox <- gtkVBox()
    vbox$packStart(frame, FALSE)
    vbox$packStart(button2, FALSE)
    vpaned$add(vbox)
    ## add each control
    sapply(.controls, function(control) { 
      control$make_gui(cont=control_table, handler=.self$change_handler)
    })
    mrv$mwindow$show()
    change_handler()
  })


.Manipulate$methods(  
  initialize = function(code, ...) {
    controls <- .resolveVariableArguments(list(...))
    initFields(.code = code, .controls = controls)
    validate_controls()
    callSuper()
  })


.manipulate <- function(parent,`_expr`,...) {
  manip <- .Manipulate$new(substitute(`_expr`),...)
  manip$execute(parent)
}

.ManipulateControls <- setRefClass("ManipulateControls", fields=list(l="list",widget = "ANY",label="ANY",initial="ANY"))


.ManipulateControls$methods(
  validate_inputs = function(...) {"Validate input code"},
  get_value = function(...) {"Get value of widget"}
)



.ManipulateControls$methods(
  make_gui = function(cont) {
  "Create widget, then add to table"
  ## cont a GtkTable instance
  nrows <- cont['n-rows']
  label_widget <- gtkLabel(label)
  label_widget['xalign'] <- 1
  cont$attach(label_widget, 0, 1, nrows, nrows + 1,
              xoptions = "shrink", yoptions = "shrink"
  )
  cont$attach(widget, 1, 2, nrows, nrows + 1,
              xoptions = c("expand", "fill"),
              yoptions = "")
})

# define the slider
.slider <- function(parent, min, max, initial = parent$text, label = NULL, step = -1, ticks = TRUE){
  .Slider$new(parent=parent, min, max, initial = initial, label = label, step = step, ticks = ticks)
}



.Slider <- setRefClass("Slider", contains = "ManipulateControls")

.Slider$methods(
  initialize = function(parent, min, max, initial=parent$text,
                        label = NULL, step = -1, ticks = TRUE) {
    validate_inputs(min, max, initial, step, ticks)
    ## create slider and return it
    slider <- list(type = 0, min = min, max = max, step = step, ticks = ticks)
    initFields(l = slider, label = label, initial = initial)
    .self
  })

.Slider$methods(
  make_gui = function(cont, handler, ...) {
    widget <<- gtkHScale(min = l$min, max = l$max, step = l$step)
    widget$setValue(initial)
    gSignalConnect(widget, "value-changed", handler)
    callSuper(cont)
  },
  get_value = function() {
    as.numeric(widget$getValue())
  })




####################################################################
###############
##############
#############
#########
####################################################################
####################################################################
.writeResult <- function(x){
  folder_name = paste(x$filepath,x$filename,sep="")
  if (!file.exists(folder_name)){
    dir.create(file.path(x$filepath, x$filename))
  }
  result_name = paste(folder_name,"Estimated_Result.txt",sep="/")
  a = summary(x$est)
  sink(result_name)
  print(a)
  sink()
  readLines(result_name, n = -1L, warn = FALSE)
}

.writeFitted <- function(x,accuracy.type){
  folder_name = paste(x$filepath,x$filename,sep="")
  if (!file.exists(folder_name)){
    dir.create(file.path(x$filepath, x$filename))
  }
  result_name = paste(folder_name,"Fitted_Result.txt",sep="/")
  a = fitted(x$est, accuracy.type = accuracy.type)
  sink(result_name)
  accuracy.type = tolower(a$accuracy.type)
  cat('Diagnostic accuracies ')
  if(accuracy.type=="sens" || accuracy.type=="tpr"){cat('true positive rate (sensitivity): ')}
  if(accuracy.type=="spec" || accuracy.type=="tnr"){cat('true negative rate (specificity): ')}
  if(accuracy.type=="fpr"){cat('false positive rate (1-specificity): ')}
  if(accuracy.type=="fnr"){cat('false negative rate (1-sensitivity): ')}
  if(accuracy.type=="lrpos"){cat('positive likelihood ratio (LR+): ')}
  if(accuracy.type=="lrneg"){cat('negative likelihood ratio (LR-): ')}
  if(accuracy.type=="dor"){cat('diagnostic odds ratio (DOR): \n')}
  if(accuracy.type=="ldor"){cat('log diagnostic odds ratio (LDOR): ')}
  if(accuracy.type=="rd"){cat('risk difference (RD): \n')}
  if(accuracy.type=="llrpos"){cat('log positive likelihood ratio (LLR+): ')}
  if(accuracy.type=="llrneg"){cat('log negative likelihood ratio (LLR-): ')}
  cat(paste('Based on \"nsample\" = ', x$est$misc$nsample,". \n",sep=""))
  print(round(a$fitted.value,3))
  cat("\n")
  sink()
  readLines(result_name, n = -1L, warn = FALSE)
}

.writeRcode <- function(x){
  folder_name = paste(x$filepath,x$filename,sep="")
  if (!file.exists(folder_name)){
    dir.create(file.path(x$filepath, x$filename))
  }
  rcode_name = paste(folder_name,"R_code.txt",sep="/")
  sink(rcode_name)
  cat("library(INLA)")
  cat("\n")
  cat("library(meta4diag)")
  cat("\n")
  if(x$fileform=="txt"){
    cat(paste("data = readData(\"",paste(x$filepath,x$filename,sep=""),"\")",sep=""))
  }
  if(x$fileform=="csv"){
    cat(paste("read.csv(\"",paste(x$filepath,x$filename,sep=""),"\")",sep=""))
  }
  if(x$fileform == "rdata" || x$fileform == "rda"){
    cat(paste("load(\"",paste(x$filepath,x$filename,sep=""),"\")",sep=""))
  }
  cat("\n")
  cat("\n")
  if(any(c(x$var.prior,x$var2.prior,x$cor.prior)=="Invwishart")){
    cat(paste("outprior = makePriors(var.prior = \"",x$var.prior,"\", 
                 wishart.par = c(",paste(x$wishart.par,collapse = ", "),"),  
                 init = c(0.01, 0.01, -0.1))",sep=""))
  }else{
    cat(paste("outprior = makePriors(var.prior = \"",x$var.prior,"\", 
                 var2.prior = \"",x$var2.prior,"\",
                 cor.prior = \"",x$cor.prior,"\",
                 var.par = c(",paste(x$var.par,collapse = ", "),"),  
                 var2.par = c(",paste(x$var2.par,collapse = ", "),"), 
                 cor.par = c(",paste(x$cor.par,collapse = ", "),"),
                 init = c(0.01, 0.01, -0.1))",sep=""))
  }
  cat("\n")
  cat("\n")
  if(!is.null(x$modality)){
    if(!all(is.na(x$partialdata))){
      modality_ind = x$partialdata[which(!is.na(x$partialdata))]
      if(is.null(x$covariates)){
        cat(paste("outdata = makeData(data = data[,\"",x$modality,"\" == c(\"",paste(modality_ind,collapse="\",\""),"\")], model.type = ",x$model.type,", modality = \"",x$modality,"\", covariates = NULL)",sep=""))
      }else{
        cat(paste("outdata = makeData(data = data[,\"",x$modality,"\" == c(\"",paste(modality_ind,collapse="\",\""),"\")], model.type = ",x$model.type,", modality = \"",x$modality,"\", covariates = ",x$covariates,")",sep=""))
      }
    }else{
      if(is.null(x$covariates)){
        cat(paste("outdata = makeData(data = data, model.type = ",x$model.type,", modality = NULL, covariates = NULL)",sep=""))
      }else{
        cat(paste("outdata = makeData(data = data, model.type = ",x$model.type,", modality = NULL, covariates = ",x$covariates,")",sep=""))
      }
    }
  }else{
    if(is.null(x$covariates)){
      cat(paste("outdata = makeData(data = data, model.type = ",x$model.type,", modality = NULL, covariates = NULL)",sep=""))
    }else{
      cat(paste("outdata = makeData(data = data, model.type = ",x$model.type,", modality = NULL, covariates = ",x$covariates,")",sep=""))
    }
  }
  cat("\n")
  cat("\n")
  if(all(is.na(x$flevel))){
    if(x$verbose){
      cat(paste("model = runModel(outdata = outdata, outprior = outpriors, model.type = ",x$model.type, ", link = \"",x$link,"\", verbose =  TRUE)",sep=""))
    }else{
      cat(paste("model = runModel(outdata = outdata, outprior = outpriors, model.type = ",x$model.type, ", link = \"",x$link,"\", verbose =  FALSE)",sep=""))
    }
  }else{
    if(x$verbose){
      cat(paste("model = runModel(outdata = outdata, outprior = outpriors, model.type = ",x$model.type, ", link = \"",x$link,"\", level = ",paste(x$level,collapse = ", "),", verbose =  TRUE)",sep=""))
    }else{
      cat(paste("model = runModel(outdata = outdata, outprior = outpriors, model.type = ",x$model.type, ", link = \"",x$link,"\", level = ",paste(x$level,collapse = ", "),", verbose =  FALSE)",sep=""))
    }
  }
  cat("\n")
  cat("\n")
  cat(paste("res = makeObject(outdata = outdata, outprior = outpriors, model = model, nsample = ",x$nsample,")"))
  sink()
  readLines(rcode_name, n = -1L, warn = FALSE)
}

####################################################################
###############
##############
#############            Open prior
#########
####################################################################
####################################################################
.open_prior_Sigma <- function(widget, window){
  dialog = gtkFileChooserDialog("Choose file", window, "open",
                                "gtk-cancel",GtkResponseType["cancel"],"gtk-open",
                                GtkResponseType["accept"])
  if (dialog$run() == GtkResponseType["accept"]){
    fullname = dialog$getFilename()
    filename = unlist(strsplit(basename(fullname), "[.]"))
    priorname = filename[1]
    priorform = tolower(filename[2])
    filepath = paste(dirname(fullname),"/",sep="")
    
    if(priorform == "txt"){
      df = read.table(fullname, header=TRUE)
      if(!is.data.frame(df)){
        txtwrong <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                     "Data file should be a data.frame!")
        if (txtwrong$run()==GtkResponseType["ok"]){}
        txtwrong$destroy()
      }
    }
    if(priorform == "rdata" || priorform == "rda"){
      df = load(fullname)
      df = get(df)
      if(!is.data.frame(df)){
        rdawrong <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                     "Data file should be a data.frame!")
        if (rdawrong$run()==GtkResponseType["ok"]){}
        rdawrong$destroy()
      }
    }
    if(priorform == "csv"){
      df = read.csv(fullname)
      if(!is.data.frame(df)){
        csvwrong <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                     "Data file should be a data.frame!")
        if (csvwrong$run()==GtkResponseType["ok"]){}
        csvwrong$destroy()
      }
    }
    if(is.data.frame(df)){
      hide_parent_box <- gtkVBox()
      hide_parent <- gtkEntry()
      hide_parent["text"]="8"
      hide_parent_box$packStart(hide_parent)
      hide_parent_box$hide()
      if(is.null(colnames(df))){
        df_no_col_name_dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                                 "Column names of data are not given! 
                                                 Please give the names to indicate \"x\" (variance) and \"y\" (pi(variance))!")
        if (df_no_col_name_dialog$run()==GtkResponseType["ok"]){}
        df_no_col_name_dialog$destroy()
      }
      if(dim(df)[2]!=2){
        df_no_2col_dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                                 "Data frame should have 2 columns indicating 
                                                  the value of variance points and 
                                                  the corresponding prior densities of these points!")
        if (df_no_2col_dialog$run()==GtkResponseType["ok"]){}
        df_no_2col_dialog$destroy()
      }
      
      if(!all(c("x","y") %in% colnames(df))){
        df_colname_wrong_dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                             "Data frame should have 2 columns indicating 
                                                  \"x\", the value of variance points and 
                                                  \"y\", the corresponding prior densities of these points!")
        if (df_colname_wrong_dialog$run()==GtkResponseType["ok"]){}
        df_colname_wrong_dialog$destroy()
      }
    }
    mrv$priorfile = df
    xlim = 8
    .manipulate(parent=list(hide_parent),
                .priorTableSigma(xmax = xlim,mrv=mrv),
                xlim = .slider(parent=hide_parent,min=1,max=20,step=1))
  }
  dialog$destroy()
}


.open_prior_Rho <- function(widget, window){
  dialog = gtkFileChooserDialog("Choose file", window, "open",
                                "gtk-cancel",GtkResponseType["cancel"],"gtk-open",
                                GtkResponseType["accept"])
  if (dialog$run() == GtkResponseType["accept"]){
    fullname = dialog$getFilename()
    filename = unlist(strsplit(basename(fullname), "[.]"))
    priorname = filename[1]
    priorform = tolower(filename[2])
    filepath = paste(dirname(fullname),"/",sep="")
    
    if(priorform == "txt"){
      df = read.table(fullname, header=TRUE)
      if(!is.data.frame(df)){
        txtwrong <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                     "Data file should be a data.frame!")
        if (txtwrong$run()==GtkResponseType["ok"]){}
        txtwrong$destroy()
      }
    }
    if(priorform == "rdata" || priorform == "rda"){
      df = load(fullname)
      df = get(df)
      if(!is.data.frame(df)){
        rdawrong <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                     "Data file should be a data.frame!")
        if (rdawrong$run()==GtkResponseType["ok"]){}
        rdawrong$destroy()
      }
    }
    if(priorform == "csv"){
      df = read.csv(fullname)
      if(!is.data.frame(df)){
        csvwrong <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                     "Data file should be a data.frame!")
        if (csvwrong$run()==GtkResponseType["ok"]){}
        csvwrong$destroy()
      }
    }
    if(is.data.frame(df)){
      hide_parent_box <- gtkVBox()
      hide_parent <- gtkEntry()
      hide_parent["text"]="8"
      hide_parent_box$packStart(hide_parent)
      hide_parent_box$hide()
      if(is.null(colnames(df))){
        df_no_col_name_dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                                 "Column names of data are not given! 
                                                 Please give the names to indicate \"x\" (variance) and \"y\" (pi(variance))!")
        if (df_no_col_name_dialog$run()==GtkResponseType["ok"]){}
        df_no_col_name_dialog$destroy()
      }
      if(dim(df)[2]!=2){
        df_no_2col_dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                             "Data frame should have 2 columns indicating 
                                                  the value of variance points and 
                                                  the corresponding prior densities of these points!")
        if (df_no_2col_dialog$run()==GtkResponseType["ok"]){}
        df_no_2col_dialog$destroy()
      }
      
      if(!all(c("x","y") %in% colnames(df))){
        df_colname_wrong_dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                                   "Data frame should have 2 columns indicating 
                                                  \"x\", the value of variance points and 
                                                  \"y\", the corresponding prior densities of these points!")
        if (df_colname_wrong_dialog$run()==GtkResponseType["ok"]){}
        df_colname_wrong_dialog$destroy()
      }
    }
    mrv$priorfile = df
    xlim = 1
    .manipulate(parent=list(hide_parent),
                .priorTableRho(xmax = xlim, mrv=mrv),
                xlim = .slider(parent=hide_parent,min=-1,max=1,step=1))
  }
  dialog$destroy()
}

####################################################################
###############
##############
#############             Open_cb
#########
####################################################################
####################################################################
  

.open_cb <- function(widget, window){
  dialog = gtkFileChooserDialog("Choose file", window, "open",
                                "gtk-cancel",GtkResponseType["cancel"],"gtk-open",
                                GtkResponseType["accept"])
  if (dialog$run() == GtkResponseType["accept"]){
    total_page_num = mrv$notebook$getNPages()
    if(total_page_num > 1){
      mrv$notebook$setCurrentPage(0)
      sapply(c((total_page_num-1):1), function(page_num){mrv$notebook$removePage(page_num)})
      mrv$datafile = NULL
    }
    
    fullname = dialog$getFilename()
    filename = unlist(strsplit(basename(fullname), "[.]"))
    mrv$fullname = fullname
    mrv$filename = filename[1]
    mrv$fileform = tolower(filename[2])
    mrv$filepath = paste(dirname(fullname),"/",sep="")
    
    if(mrv$fileform == "txt"){
      df = read.table(fullname, header=TRUE)
      if(!is.data.frame(df)){
        txtwrong <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                     "Data file should be a data.frame!")
        if (txtwrong$run()==GtkResponseType["ok"]){}
        txtwrong$destroy()
      }
    }
    if(mrv$fileform == "rdata" || mrv$fileform == "rda"){
      df = load(fullname)
      df = get(df)
      if(!is.data.frame(df)){
        rdawrong <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                     "Data file should be a data.frame!")
        if (rdawrong$run()==GtkResponseType["ok"]){}
        rdawrong$destroy()
      }
    }
    if(mrv$fileform == "csv"){
      df = read.csv(fullname)
      if(!is.data.frame(df)){
        csvwrong <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                     "Data file should be a data.frame!")
        if (csvwrong$run()==GtkResponseType["ok"]){}
        csvwrong$destroy()
      }
    }
    if(is.data.frame(df)){
      ######## check data
      if(is.null(colnames(df))){
        df_no_col_name_dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                                 "Column names of data are not given! 
                                                 Please give the names to indicate \"TP\", \"FN\", \"FP\", \"TN\"!")
        if (df_no_col_name_dialog$run()==GtkResponseType["ok"]){}
        df_no_col_name_dialog$destroy()
      }
      
      datanames = tolower(colnames(df))
      colnames(df) = datanames
      I = dim(df)[1]
      fic = c("tp","tn","fp","fn") # four important components
      is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
      ###### check main data
      if(!("studynames" %in% datanames)){
        #if("references" %in% datanames){
          #studynames = df$references
          #df = data.frame(studynames=studynames, df[,-which(datanames=="references")])
        #}else{
          studynames = paste("study[",c(1:I),"]",sep="")
          df = data.frame(studynames=studynames, df)
        #}
      }
      datanames = colnames(df)
      if(all(fic %in% datanames)){
        tf.integer = apply(cbind(df$tp,df$tn,df$fp,df$fn),2,function(x) all(is.wholenumber(x)))
        tf.pos = apply(cbind(df$tp,df$tn,df$fp,df$fn)>=0,2,function(x) all(x))
        if(all(tf.integer)){
          if(all(tf.pos)){
            
          }else{
            fvp = paste("data$",fic[tf.pos==FALSE],sep="",collapse=" and ")
            df_neg_dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                             paste(fvp," has some negative value!!!",sep=""))
            if (df_neg_dialog$run()==GtkResponseType["ok"]){}
            df_neg_dialog$destroy()
          }
        }else{
          fvi = paste("data$",fic[tf.integer==FALSE],sep="",collapse=" and ")
          df_ni_dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                          paste(fvi," has some non-integer value!!!",sep=""))
          if (df_ni_dialog$run()==GtkResponseType["ok"]){}
          df_ni_dialog$destroy()
        }
      }else{
        ffic = toupper(paste(fic[!(fic %in% datanames)],collapse=" "))
        df_ffic_dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                          paste("Data is not ok!",ffic,"are missing!!!",sep=" "))
        if (df_ffic_dialog$run()==GtkResponseType["ok"]){}
        df_ffic_dialog$destroy()
      }
      if(dim(df)[2]>5){
        mrv$pt_covariates = datanames[-which(datanames %in% c("studynames","tp","tn","fp","fn"))]
      } else{
        mrv$pt_covariates = NULL
      }
      if(length(mrv$pt_covariates)==0){
        mrv$pt_covariates = NULL
      }
      .load_spreadsheet(df, basename(fullname))
    }
    mrv$datafile = df
  }
  dialog$destroy()
}


####################################################################
###############
##############
#############            Some Actions
#########
####################################################################
####################################################################

.someAction <- function(action,env) 
  env$statusbar$push(env$statusbar$getContextId("message"), 
                 action$getName())

.level_box = function(level){
  box <- gtkVBox(homogeneous=FALSE,spacing=10)
  level_length = length(level)
  
  level_table <- gtkTable(rows=ceiling(0.5*level_length),columns=2,homogeneous=FALSE)
  level_table$setColSpacings(5)
  level_table$setRowSpacings(3)
  
  level_checkgp = list()
  for(i in 1:level_length){
    level_checkgp[[i]] <- gtkCheckButton(level[i])
    level_checkgp[[i]]$setActive(TRUE)
  }
  if((level_length%%2)==0){ # even number of covariates
    for(i in 1:(0.5*level_length)){
      check_ind1 = 2*i-1
      check_ind2 = 2*i
      level_table$attach(level_checkgp[[check_ind1]],left.attach=0,1, top.attach=(i-1),i,
                         xoptions=c("expand","fill"),yoptions="")
      level_table$attach(level_checkgp[[check_ind2]],left.attach=1,2, top.attach=(i-1),i,
                         xoptions=c("expand","fill"),yoptions="")
    }
  }else{ 
    if(level_length==1){
      level_table$attach(level_checkgp[[1]],left.attach=0,1, top.attach=0,1,
                         xoptions=c("expand","fill"),yoptions="")
    }else{ 
      for(i in 1:floor(0.5*level_length)){
        check_ind1 = 2*i-1
        check_ind2 = 2*i
        level_table$attach(level_checkgp[[check_ind1]],left.attach=0,1, top.attach=(i-1),i,
                           xoptions=c("expand","fill"),yoptions="")
        level_table$attach(level_checkgp[[check_ind2]],left.attach=1,2, top.attach=(i-1),i,
                           xoptions=c("expand","fill"),yoptions="")
      }
      level_table$attach(level_checkgp[[level_length]],left.attach=0,1, top.attach=floor(0.5*level_length),ceiling(0.5*level_length),
                         xoptions=c("expand","fill"),yoptions="")
    }
  }
  box$packStart(level_table, expand=FALSE,fill=FALSE)
  return(list(box = box, level_checkgp = level_checkgp))
}

.load_spreadsheet <- function(df, name){
  model <- .create_tree_model(df)
  tree_view <- .create_tree_view(model)
  entry <- .create_entry(model$getModel()$getModel())
  hbox <- gtkHBox(FALSE, 5)
  hbox$packStart(gtkLabel("Filter expression:"), FALSE, FALSE, 0)
  hbox$packStart(entry, TRUE, TRUE, 0)
  
  mrv$data_page <- gtkVBox(FALSE, 5)
  data_swindow <- gtkScrolledWindow()
  data_swindow$setPolicy("automatic", "automatic")
  data_swindow$add(tree_view)
  mrv$data_page$packStart(data_swindow, TRUE, TRUE, 0)
  mrv$data_page$packStart(hbox, FALSE, FALSE, 0)
  if (missing(name)) name <- deparse(substitute(df))
  
  data_page_event = .pageEventsLabel(name)
  mrv$notebook$appendPage(mrv$data_page, data_page_event)
  
  data_page_num = mrv$notebook$pageNum(mrv$data_page)
  mrv$notebook$setCurrentPage(data_page_num)
  
  ##### sidebar data control
  mrv$sidedata <- gtkVBox(homogeneous=FALSE)
  
  mrv$sidedata_window <- gtkScrolledWindow()
  mrv$sidedata_window["shadow-type"] = "none"
  mrv$sidedata_window$addWithViewport(mrv$sidedata)
  mrv$sidedata_window$setPolicy("never", "automatic")
  
  mrv$sidenote$removePage(mrv$sidedata_page_num)
  data_page_label = gtkLabel("Data Control Panel")
  data_page_label$setAngle(270)
  mrv$sidenote$InsertPage(mrv$sidedata_window, data_page_label, mrv$sidedata_page_num)
  mrv$sidedata_page_num = mrv$sidenote$pageNum(mrv$sidedata_window)
  mrv$sidenote$setCurrentPage(mrv$sidedata_page_num)
  
  gSignalConnect(data_page_event, "button-press-event", function(...) {
    mrv$notebook$setCurrentPage(data_page_num)
    mrv$sidenote$setCurrentPage(mrv$sidedata_page_num)
    return(TRUE)
  })
  
  mrv$statusbar$push(mrv$info, paste("Dataset", name, "loaded.")) 
  
  ### modality control
  mrv$sidedcont1 <- gtkFrame("Modality")
  mrv$sidedcont1["border-width"]=10
  
  mrv$sidedcont1_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidedcont1_main["border-width"]=10
  
#   if(!is.null(df$modality)){
#     level_modality = as.character(unique(df$modality))
#     length_modality = length(level_modality)
#     
#     mrv$sided_radio <- gtkVBox()
#     mrv$sided_radiogp <- list()
#     if(length_modality==1){
#       mrv$sided_radiogp$all <- gtkRadioButton(label = "All")
#     } else{
#       mrv$sided_radiogp[["All"]] <- gtkRadioButton(label = "All")
#       for(i in 1:length_modality){
#         mrv$sided_radiogp[[level_modality[i]]] <- gtkRadioButton(mrv$sided_radiogp, label = paste("Subdata - ", level_modality[i], sep=""))
#       }
#     }
#     sapply(mrv$sided_radiogp, mrv$sided_radio$packStart)
#     mrv$sided_radio[[1]]$setActive(TRUE) 
#     sapply(mrv$sided_radiogp,'[',"active")
#     mrv$sided_radio_align <- gtkAlignment(xalign = 0)
#     mrv$sided_radio_align$add(mrv$sided_radio)
#     
#     mrv$sidedcont1_main$packStart(mrv$sided_radio_align, expand=FALSE,fill=FALSE)
#     
#     sapply(mrv$sided_radiogp, gSignalConnect, "toggled",
#            f = function(button, ...){
#              if(button['active']){
#                mrv$choosedata = button$getLabel()
#              } 
#            })
#   }else{
#     modality_null_label = gtkLabel("NULL")
#     modality_null_label["width-request"]=100
#     modality_null_label["height-request"]=100
#     mrv$sidedcont1_main$packStart(modality_null_label, expand=FALSE,fill=FALSE)
#   }
  ######
  mrv$pt_modality = mrv$pt_covariates
  if(!is.null(mrv$pt_modality)){ 
    
    mrv$sided_radiogp <- list()
    mrv$sided_radiogp[[1]] <- gtkRadioButton(label = "None")
    for(i in 1:length(mrv$pt_modality)){
      mrv$sided_radiogp[[i+1]] <- gtkRadioButton(mrv$sided_radiogp, label = mrv$pt_modality[i])
    }
    
    length_modality = length(mrv$sided_radiogp)
    
    mrv$sidedcont1_table <- gtkTable(rows=ceiling(0.5*length_modality),columns=2,homogeneous=FALSE)
    mrv$sidedcont1_table$setColSpacings(5)
    mrv$sidedcont1_table$setRowSpacings(3)
    
    
    if((length_modality%%2)==0){ # even number of covariates
      for(i in 1:(0.5*length_modality)){
        check_ind1 = 2*i-1
        check_ind2 = 2*i
        mrv$sidedcont1_table$attach(mrv$sided_radiogp[[check_ind1]],left.attach=0,1, top.attach=(i-1),i,
                                    xoptions=c("expand","fill"),yoptions="")
        mrv$sidedcont1_table$attach(mrv$sided_radiogp[[check_ind2]],left.attach=1,2, top.attach=(i-1),i,
                                    xoptions=c("expand","fill"),yoptions="")
      }# end for loop
    }else{ # odd number of covairates
      if(length_modality==1){
        mrv$sidedcont1_table$attach(mrv$sided_radiogp[[1]],left.attach=0,1, top.attach=0,1,
                                    xoptions=c("expand","fill"),yoptions="")
      }else{ # if number of covariates is 3,5,7,9....
        for(i in 1:floor(0.5*length_modality)){
          check_ind1 = 2*i-1
          check_ind2 = 2*i
          mrv$sidedcont1_table$attach(mrv$sided_radiogp[[check_ind1]],left.attach=0,1, top.attach=(i-1),i,
                                      xoptions=c("expand","fill"),yoptions="")
          mrv$sidedcont1_table$attach(mrv$sided_radiogp[[check_ind2]],left.attach=1,2, top.attach=(i-1),i,
                                      xoptions=c("expand","fill"),yoptions="")
        }
        mrv$sidedcont1_table$attach(mrv$sided_radiogp[[length_modality]],left.attach=0,1, 
                                    top.attach=floor(0.5*length_modality),ceiling(0.5*length_modality),
                                    xoptions=c("expand","fill"),yoptions="")
      }
    }
    
    mrv$sidedcont1_main$packStart(mrv$sidedcont1_table, expand=FALSE,fill=FALSE)

    mrv$mod_level = lapply(1:length(mrv$pt_modality), function(ind) as.character(unique(df[,mrv$pt_modality[ind]])))
    bac = list()
    mrv$choosedata = list()
    for(k in 1:length(mrv$pt_modality)){
      bac[[k]] = .level_box(mrv$mod_level[[k]])
      bac[[k]]$box$hide()
      mrv$sidedcont1_main$packStart(bac[[k]]$box, expand=FALSE,fill=FALSE)
      
      mrv$choosedata[[k]] = mrv$mod_level[[k]]
    }
    
    sapply(c(1:length(mrv$pt_modality)), function(ind){
      sapply(c(1:length(mrv$mod_level[[ind]])), function(kk){
        gSignalConnect(bac[[ind]]$level_checkgp[[kk]], "toggled", f = function(button, ...){
          if(button['active']){
            mrv$choosedata[[ind]][kk] = button$getLabel()
          }else{
            mrv$choosedata[[ind]][kk] = NA
          }
        })
      })
    })
    
    
    
    
    mrv$modality = NULL
    sapply(mrv$sided_radiogp, gSignalConnect, "toggled",
           f = function(button, ...){
             if(button['active']){
               mrv$modality = button$getLabel()
               if(mrv$modality=="None"){
                 mrv$modality = NULL
                 for(ll in 1:length(mrv$pt_modality)){
                   bac[[ll]]$box$hide()
                 }
               }else{
                 which_mod = which(mrv$pt_modality %in% mrv$modality)
                 bac[[which_mod]]$box$show()
                 for(ll in c(1:length(mrv$pt_modality))[-which_mod]){
                   bac[[ll]]$box$hide()
                 }
               }
             }
             if(any(mrv$fcovariates==mrv$modality)){
               mcsamedialog <- gtkMessageDialog(mrv$main_window,"destroy-with-parent","warning","ok",
                                                "Covariates can not be same as modality!!!!")
               if (mcsamedialog$run() == GtkResponseType["ok"]){
                 mrv$sided_radiogp[[1]]$setActive(TRUE)
               }
               mcsamedialog$destroy()
             }
           })
    
    
  }else{  # no modality
    modality_null_label = gtkLabel("NULL")
    modality_null_label["width-request"]=100
    modality_null_label["height-request"]=100
    mrv$sidedcont1_main$packStart(modality_null_label, expand=FALSE,fill=FALSE)
    mrv$modality = NULL
  }
  
  mrv$sidedcont1$add(mrv$sidedcont1_main)
  mrv$sidedata$packStart(mrv$sidedcont1,expand=FALSE,fill=FALSE)
  ####### covariates
  mrv$sidedcont2 <- gtkFrame("Covariates")
  mrv$sidedcont2["border-width"]=10
  
  mrv$sidedcont2_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidedcont2_main["border-width"]=10
  if(!is.null(mrv$pt_covariates)){ 
    length_covariates = length(mrv$pt_covariates)
    
    mrv$sidedcont2_table <- gtkTable(rows=ceiling(0.5*length_covariates),columns=2,homogeneous=FALSE)
    mrv$sidedcont2_table$setColSpacings(5)
    mrv$sidedcont2_table$setRowSpacings(3)
    
    mrv$sided_checkgp = list()
    for(i in 1:length_covariates){
      mrv$sided_checkgp[[i]] <- gtkCheckButton(mrv$pt_covariates[i])
      mrv$sided_checkgp[[i]]$setActive(FALSE)
    }
    if((length_covariates%%2)==0){ # even number of covariates
      for(i in 1:(0.5*length_covariates)){
        check_ind1 = 2*i-1
        check_ind2 = 2*i
        mrv$sidedcont2_table$attach(mrv$sided_checkgp[[check_ind1]],left.attach=0,1, top.attach=(i-1),i,
                                    xoptions=c("expand","fill"),yoptions="")
        mrv$sidedcont2_table$attach(mrv$sided_checkgp[[check_ind2]],left.attach=1,2, top.attach=(i-1),i,
                                    xoptions="",yoptions="")
      }# end for loop
    }else{ # odd number of covairates
      if(length_covariates==1){
        mrv$sidedcont2_table$attach(mrv$sided_checkgp[[1]],left.attach=0,1, top.attach=0,1,
                                    xoptions=c("expand","fill"),yoptions="")
      }else{ # if number of covariates is 3,5,7,9....
        for(i in 1:floor(0.5*length_covariates)){
          check_ind1 = 2*i-1
          check_ind2 = 2*i
          mrv$sidedcont2_table$attach(mrv$sided_checkgp[[check_ind1]],left.attach=0,1, top.attach=(i-1),i,
                                      xoptions=c("expand","fill"),yoptions="")
          mrv$sidedcont2_table$attach(mrv$sided_checkgp[[check_ind2]],left.attach=1,2, top.attach=(i-1),i,
                                      xoptions="",yoptions="")
        }
        mrv$sidedcont2_table$attach(mrv$sided_checkgp[[length_covariates]],left.attach=0,1, top.attach=floor(0.5*length_covariates),ceiling(0.5*length_covariates),
                                    xoptions=c("expand","fill"),yoptions="")
      }
    }
    mrv$sidedcont2_main$packStart(mrv$sidedcont2_table, expand=FALSE,fill=FALSE)
    
    mrv$fcovariates = rep(list(NA), length_covariates)
    sapply(c(1:length_covariates), function(ind){
      gSignalConnect (mrv$sided_checkgp[[ind]], "toggled", f = function(button, ...){
        if(button['active']){
          mrv$fcovariates[[ind]] = button$getLabel()
        }else{
          mrv$fcovariates[[ind]] = NA
        }
        if(any(mrv$fcovariates==mrv$modality)){
          cmsamedialog <- gtkMessageDialog(mrv$main_window,"destroy-with-parent","warning","ok",
                                     "Covariates can not be same as modality!!!!")
          if (cmsamedialog$run() == GtkResponseType["ok"]){
            button$setActive(FALSE)
          }
          cmsamedialog$destroy()
        }
      })
    })
  }else{  # no covariates
    covariates_null_label = gtkLabel("NULL")
    covariates_null_label["width-request"]=100
    covariates_null_label["height-request"]=100
    mrv$sidedcont2_main$packStart(covariates_null_label, expand=FALSE,fill=FALSE)
    mrv$covariates = NULL
  }
  mrv$sidedcont2$add(mrv$sidedcont2_main)
  mrv$sidedata$packStart(mrv$sidedcont2,expand=FALSE,fill=FALSE)
}


.stop_cb <- function(widget, window){
  dialog <- gtkMessageDialog(window,"destroy-with-parent","question","yes-no",
                                 "Do you really want to stop the running?")
  if (dialog$run() == GtkResponseType["yes"]){
    total_page_num = mrv$notebook$getNPages()
    if(total_page_num >= 2){
      mrv$notebook$setCurrentPage(1)
      sapply(c((total_page_num-1):1), function(page_num){mrv$notebook$removePage(page_num)})
      mrv$statusbar$push(mrv$info, paste("Ready."))
      mrv$datafile = NULL
      mrv$outdata = NULL
      mrv$outpriors = NULL
    }
  }
  dialog$destroy()
}

.aboutgui <- function(widget,window){
  dialog <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                             "About meta4diag: \n
Version: 1.0.17 \n
Update date: 21 Aug 2015 \n")
  if (dialog$run() == GtkResponseType["ok"]){
  }
  dialog$destroy()
}



.quit_cb <- function(widget, window) window$destroy()

.create_tree_model <- function(df) {
  df <- cbind(rownames = rownames(df), df)
  filter_df <- cbind(filter = TRUE, df)
  model <- rGtkDataFrame(filter_df)
  filter_model <- gtkTreeModelFilterNew(model)
  filter_model$setVisibleColumn(0)
  sort_model <- gtkTreeModelSort(filter_model)
  sort_model 
}

.create_tree_view <- function(model){
  tree_view <- gtkTreeView(model)
  rdf <- model$getModel()$getModel()
  sapply(tail(seq_len(ncol(rdf)), -1), 
         function(j){
           renderer <- gtkCellRendererText()
           column <- gtkTreeViewColumn(colnames(rdf)[j], renderer, text = j - 1)
           column$setSortColumnId(j - 1)
           column$setCellDataFunc(renderer,
                                  function(column, renderer, model, iter){
                                    iter <- model$convertIterToChildIter(iter)$child.iter
                                    child <- model$getModel()
                                    iter <- child$convertIterToChildIter(iter)$child.iter
                                    i <- rdf$getPath(iter)$getIndices()[[1]] + 1
                                    renderer["text"] <- as.character(rdf[i, j])
                                  })
           tree_view$appendColumn(column)
         })
  tree_view$setHeadersClickable(TRUE)
  if (is.null(gtkCheckVersion(2, 10, 0))) tree_view$setGridLines("both")
  tree_view
}

.create_entry <- function(model){
  entry <- gtkEntry()
  gSignalConnect(entry, "activate", function(entry) {
    model[, "filter"] <<- eval(parse(text = entry$text),as.data.frame(model))
  })
  entry
}


# save_figure <- function(widget, window) {
#   dialog <- gtkFileChooserDialog("Enter a name for the file", window,
#                                  "save", "gtk-cancel", GtkResponseType["cancel"], "gtk-save",
#                                  GtkResponseType["accept"])
#   if (dialog$run() == GtkResponseType["accept"]) 
#     save_file(dialog$getFilename())
#   dialog$destroy()
# }


# 
# open_cb <- function(widget, window) {
#   dialog <- gtkFileChooserDialog("Choose a CSV file", window, "open",
#                                  "gtk-cancel", GtkResponseType["cancel"], "gtk-open",
#                                  GtkResponseType["accept"])
#   if (dialog$run() == GtkResponseType["accept"]) {
#     df <- read.table(dialog$getFilename(), header=TRUE)
#     load_spreadsheet(df, basename(dialog$getFilename()))
#     mrv$datafile = df
#   }
#   dialog$destroy()
# }

.pageEventsLabel <- function(string){
  event <- gtkEventBox()
  event["border-width"] <- 0
  label <- gtkLabel(string)
  label["height-request"]=30
  numberchar <- nchar(string)
  label["width-chars"] <- numberchar+2
  label["xalign"] <- 0.5
  label["yalign"] <- 0.5
  event$add(label)
  return(event)
}

.scrolledWindowWithViewport <- function(widget){
  viewport <- gtkViewport()
  viewport$modifyBg(viewport$getState(),as.GdkColor("white"))
  viewport$add(widget)
  window <- gtkScrolledWindow()
  window["shadow-type"] = "none"
  window$add(viewport)
  window$setPolicy("automatic", "automatic")
  return(window)
}
###########################################################
##########    main execute function
###########################################################
.executeFile <- function(widget, window){
  if (!(sum(search()=="package:INLA"))==1){
    INLA_dialog <- gtkMessageDialog(NULL,"destroy-with-parent","warning","ok",paste("R Package \"INLA\" is not loaded.","\n",
                                                                                    "You have to load it to make sure meta4diag works.","\n",
                                                                                    "Thank you!",sep=""))
    INLA_dialog["title"] <- "Please Load INLA!!!"
    if (INLA_dialog$run()==GtkResponseType["ok"]){
      INLA_dialog$destroy()
    }
  }else{
  
  total_page_num = mrv$notebook$getNPages()
  if(total_page_num > 2){
    mrv$notebook$setCurrentPage(1)
    sapply(c((total_page_num-1):2), function(page_num){mrv$notebook$removePage(page_num)})
  }
  
  transp_frame <- function(){
    frame <- gtkFrame()
    frame$setShadowType("none")
    return(frame)
  }
  
  plot_button <- function(){
    Button <- gtkButton("Plot")
    Button['image'] <- gtkImage(stock.id="gtk-refresh",size="button")
    Button$setAlignment(xalign=0.5,yalign=0.5)
    Button$setImagePosition("left")
    return(Button)
  }
  
  save_button <- function(){          
    Button <- gtkButton("Save")
    Button['image'] <- gtkImage(stock.id="gtk-save",size="button")
    Button$setAlignment(xalign=0.5,yalign=0.5)
    Button$setImagePosition("left")
    return(Button)
  }
  
  bottomLabel <- function(string, xalign){
    label <- gtkLabel(string)
    label["xalign"] <- xalign
    return(label)
  }
  
  bottomEntry <- function(string, wchars, mlength){
    entry <- gtkEntry()
    entry["text"] <- string
    entry["width-chars"] <- wchars
    entry["max-length"] <- mlength
    return(entry)
  }

  # check if data input and all prior name are given
  if(is.null(mrv$datafile)){
    nodatadialog <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                     "Please select a data file!")
    if (nodatadialog$run()==GtkResponseType["ok"]){
      mrv$sidenote$setCurrentPage(mrv$sidedata_page_num)
    }
    nodatadialog$destroy()
  }else{
    if(all(c(mrv$var.prior,mrv$var2.prior,mrv$cor.prior)!="Invwishart")){
      if(is.null(mrv$var.prior)){
        notau1dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                        "Please set prior for the first variance!")
        if (notau1dialog$run()==GtkResponseType["ok"]){
          mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
        }
        notau1dialog$destroy()
      }else{
        if(is.null(mrv$var2.prior)){
          notau2dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                          "Please set prior for the second variance!")
          if (notau2dialog$run()==GtkResponseType["ok"]){
            mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
          }
          notau2dialog$destroy()
        } else{
          if(is.null(mrv$cor.prior)){
            norhodialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                           "Please set prior for the correlation!")
            if (norhodialog$run()==GtkResponseType["ok"]){
              mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
            }
            norhodialog$destroy()
          }
        }
      }
    }
  }
  
  if(!is.null(mrv$datafile) && !is.null(mrv$var.prior) && !is.null(mrv$var2.prior) && !is.null(mrv$cor.prior)){
    #### assume all prior flag are TRUE 
    var.flag = TRUE
    var2.flag = TRUE
    cor.flag = TRUE
    if(any(c(mrv$var.prior,mrv$var2.prior,mrv$cor.prior)=="Invwishart")){
      if(mrv$sidecont4_hide2$visible){
        var.flag = FALSE
        var2.flag = FALSE
        cor.flag = FALSE
        iwdialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                    "Please accept the prior setting for the covariance matrix!")
        if (iwdialog$run()==GtkResponseType["ok"]){
          mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
        }
        iwdialog$destroy()
      }
    }else{
      if(mrv$sidecont1_hide12$visible || mrv$sidecont1_hide22$visible || mrv$sidecont1_hide32$visible || mrv$sidecont1_hide42$visible || mrv$sidecont1_hide5$visible || mrv$sidecont1_hide6$visible){
        var.flag = FALSE
        vis1dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                      "Please accept the prior setting for the first variance!")
        if (vis1dialog$run()==GtkResponseType["ok"]){
          mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
        }
        vis1dialog$destroy()
      }else{
        if(mrv$sidecont2_hide12$visible || mrv$sidecont2_hide22$visible|| mrv$sidecont2_hide32$visible || mrv$sidecont2_hide42$visible || mrv$sidecont2_hide5$visible || mrv$sidecont2_hide6$visible){
          var2.flag = FALSE
          vis2dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                        "Please accept the prior setting for the second variance!")
          if (vis2dialog$run()==GtkResponseType["ok"]){
            mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
          }
          vis2dialog$destroy()
        }else{
          if(mrv$cor.prior=="normal"){
            if(mrv$sidecont3_hide12$visible){
              cor.flag = FALSE
              vis3dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                            "Please accept the prior setting for the correlation!")
              if (vis3dialog$run()==GtkResponseType["ok"]){
                mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
              }
              vis3dialog$destroy()
            }
          }else if(mrv$cor.prior=="PC"){
            if(any(mrv$sidecont3_input2_subhide12$visible,mrv$sidecont3_input2_subhide22$visible,mrv$sidecont3_input2_subhide32$visible)){
              cor.flag = FALSE
              vis3dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                            "Please accept the prior setting for the correlation!")
              if (vis3dialog$run()==GtkResponseType["ok"]){
                mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
              }
              vis3dialog$destroy()
            } 
          }else if(mrv$cor.prior %in% c("beta","Table")){
            if(any(mrv$sidecont3_hide32$visible, mrv$sidecont3_hide4$visible)){
              cor.flag = FALSE
              vis3dialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                            "Please accept the prior setting for the correlation!")
              if (vis3dialog$run()==GtkResponseType["ok"]){
                mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
              }
              vis3dialog$destroy()
            } 
          }
        }
      }
    }
    
    ######### if all priors are proper given
    if(all(c(var.flag,var2.flag,cor.flag))){
      if(!is.null(mrv$pt_covariates)){
        if(all(unlist(lapply(mrv$fcovariates,is.na)))){
          mrv$covariates=NULL
        }else{
          covariates_temp = unlist(mrv$fcovariates)
          if(length(which(is.na(covariates_temp)))==0){
            mrv$covariates = covariates_temp
          }else{
            mrv$covariates = covariates_temp[-which(is.na(covariates_temp))]
          } 
        } 
      }
      
      ######### make priors
      mrv$statusbar$push(mrv$info, paste("Prior making.....Please wait......")) 
      mrv$outpriors = makePriors(var.prior = mrv$var.prior, var2.prior=mrv$var2.prior, cor.prior=mrv$cor.prior,
                                 var.par = mrv$var.par, var2.par=mrv$var2.par, cor.par=mrv$cor.par,wishart.par = mrv$wishart.par,
                                 init = c(0.01, 0.01, -0.1))
      
      if(is.null(mrv$outpriors)){
        mrv$statusbar$push(mrv$info, paste("Prior making problems.....code stop running....Please contact developer......"))
        priordialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                       "Prior making problem, please check prior settings!")
        if (priordialog$run()==GtkResponseType["ok"]){
          mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
        }
        priordialog$destroy()
      }
      
      ########## make data
      mrv$statusbar$push(mrv$info, paste("Prior made. Data making......Please wait......")) 
      if(!is.null(mrv$modality)){
        ind = which(mrv$modality %in% mrv$pt_modality)
        mrv$partialdata = mrv$choosedata[[ind]]
        if(all(is.na(mrv$partialdata))){
          mrv$impdata = mrv$datafile
          mrv$modality = NULL
        }else{
          partialName = mrv$partialdata[!is.na(mrv$partialdata)]
          partialInd = which(mrv$datafile[,mrv$modality] %in% partialName)
          mrv$impdata = mrv$datafile[partialInd,]
        }
      }else{
        mrv$impdata = mrv$datafile
      }
      mrv$outdata = makeData(data = mrv$impdata, model.type = mrv$model.type, covariates = mrv$covariates, modality = mrv$modality)
      if(is.null(mrv$outdata)){
        mrv$statusbar$push(mrv$info, paste("Data making problems.....code stop running....Please contact developer......"))
        datamakedialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                          "Data making problem, please check data settings!")
        if (datamakedialog$run()==GtkResponseType["ok"]){
          mrv$sidenote$setCurrentPage(mrv$sidedata_page_num)
        }
        datamakedialog$destroy()
      }
      ######### run modal
      mrv$statusbar$push(mrv$info, paste("Prior made. Data made. Model running......Please wait......")) 
      
      if(all(is.na(mrv$flevel))){
        mrv$level = c(0.5)
      }else{
        level_temp = as.numeric(unlist(mrv$flevel))
        mrv$level = level_temp[-which(is.na(level_temp))]
      }
      mrv$model = runModel(outdata = mrv$outdata, outpriors = mrv$outpriors, link = mrv$link, quantiles = mrv$level, verbose = mrv$verbose)
      if(!mrv$model$ok){
        mrv$statusbar$push(mrv$info, paste("Model running problems.....code stop running....Please contact developer......"))
        modelmakedialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                           "Model making problem, please contact developers!")
        if (modelmakedialog$run()==GtkResponseType["ok"]){}
        modelmakedialog$destroy()
      }
      
      ######### make object
      mrv$nsample = mrv$sidemcont5_spinbutton$value
      mrv$est = makeObject(mrv$outdata, mrv$outpriors, mrv$model, mrv$nsample)
      if(is.null(mrv$est)){
        mrv$statusbar$push(mrv$info, paste("Model running problems.....code stop running....Please contact developer......"))
        oedialog <-gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                    "Object making problem, please contact developers!")
        if (oedialog$run()==GtkResponseType["ok"]){}
        oedialog$destroy()
      }
      mrv$statusbar$push(mrv$info, paste("Prior made. Data made. Model made......"))
      
      
      ###################################################
      ## construct Result page
      ###################################################
      mrv$result_page = gtkFrame()
      mrv$result_page["border-width"] = 0
      
      mrv$result_view <- gtkTextView()
      mrv$result_view["left-margin"] = 10
      mrv$result_swindow <- gtkScrolledWindow()
      mrv$result_swindow["shadow-type"] = "none"
      mrv$result_swindow$add(mrv$result_view)
      mrv$result_swindow$setPolicy("automatic", "automatic")
      
      mrv$result_buffer <- mrv$result_view$getBuffer()
      mrv$result_buffer$setText(paste("\n",.writeResult(mrv),
                                      sep="",collapse = "\n"))
      result_tag_large <- mrv$result_buffer$createTag(tag.name="large", font="Monaco 10")
      
      result_iter <- mrv$result_buffer$getBounds()
      mrv$result_buffer$applyTag(result_tag_large, result_iter$start, result_iter$end)
      
      mrv$result_page$add(mrv$result_swindow)
      result_page_event = .pageEventsLabel("Estimated Result")
      mrv$notebook$insertPage(mrv$result_page,result_page_event)
      
      result_page_num = mrv$notebook$pageNum(mrv$result_page)
      
      gSignalConnect(result_page_event, "button-press-event", function(...) {
        mrv$notebook$setCurrentPage(result_page_num)
        mrv$sidenote$setCurrentPage(mrv$sidemodel_page_num)
        return(TRUE)
      })
      
      result_page_num = mrv$notebook$pageNum(mrv$result_page)
      
      ###################################################
      ## construct marginal page and active this page
      ###################################################
      length_fixed = length(mrv$est$marginals.fixed)
      names_fixed = rownames(mrv$est$summary.fixed)
      row_fixed = ceiling(length_fixed/3)
      col_fixed = 3
      
      length_hyper = length(mrv$est$marginals.hyperpar)
      names_hyper = rownames(mrv$est$summary.hyperpar)
      row_hyper = ceiling(length_fixed/3)
      col_hyper = 3
      
      marginals_fixed_plot <- gtkDrawingAreaNew()
      marginals_fixed_plot$setSizeRequest(400, 200)
      marginals_hyper_plot <- gtkDrawingAreaNew()
      marginals_hyper_plot$setSizeRequest(400, 200)
      
      marginals_fixed_box <- gtkHBox()
      marginals_fixed_box["border-width"]=10
      marginals_fixed_box$packStart(marginals_fixed_plot, expand=T, fill=T, padding=0)
      marginals_hyper_box <- gtkHBox()
      marginals_hyper_box["border-width"]=10
      marginals_hyper_box$packStart(marginals_hyper_plot, expand=T, fill=T, padding=0)
      fixed_frame <- gtkFrame("Fixed effect")
      fixed_frame["border-width"]=10
      fixed_frame$add(marginals_fixed_box)
      
      hyper_frame <- gtkFrame("Hyperpar")
      hyper_frame["border-width"]=10
      hyper_frame$add(marginals_hyper_box)
      
      marginals_plot <- gtkVBox()
      marginals_plot$packStart(fixed_frame, expand=T, fill=T, padding=0)
      marginals_plot$packStart(hyper_frame, expand=T, fill=T, padding=0)
      
      marginals_window <- gtkScrolledWindow()
      marginals_window["shadow-type"] = "none"
      marginals_window$addWithViewport(marginals_plot)
      marginals_window$setPolicy("automatic", "automatic")
      
      
      marginals_page <- gtkVBox()
      marginals_page$packStart(marginals_window, expand=T, fill=T, padding=0)
      marginals_page$packStart(gtkLabel("Bottom Button") , expand=F, fill=F, padding=0)
      
      marginals_page_event = .pageEventsLabel("Marginals Plot")
      mrv$notebook$appendPage(marginals_page, marginals_page_event)
      
      Sys.sleep(.1) 
      asCairoDevice(marginals_fixed_plot)
      par(mfrow=c(row_fixed,col_fixed))
      for(i in 1:length_fixed){
        plot(INLA::inla.smarginal(mrv$est$marginals.fixed[[i]]),type="l",xlab=names_fixed[i],ylab="")
      }
      Sys.sleep(.1) 
      asCairoDevice(marginals_hyper_plot)
      par(mfrow=c(row_hyper,col_hyper))
      for(i in 1:length_hyper){
        plot(INLA::inla.smarginal(mrv$est$marginals.hyperpar[[i]]),type="l",xlab=names_hyper[i],ylab="")
      }
      
      marginals_page_num = mrv$notebook$pageNum(marginals_page)
      
      gSignalConnect(marginals_page_event, "button-press-event", function(...) {
        mrv$notebook$setCurrentPage(marginals_page_num)
        mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
        return(TRUE)
      })
      
      mrv$statusbar$push(mrv$info, paste("Marginals ploted."))
      
      ###################################################
      ## construct Fitted page
      ###################################################
      #mrv$fitted_page = gtkFrame()
      #mrv$fitted_page["border-width"] = 0
      
      mrv$fitted_view <- gtkTextView()
      mrv$fitted_view["left-margin"] = 10
      mrv$fitted_swindow <- gtkScrolledWindow()
      mrv$fitted_swindow["shadow-type"] = "none"
      mrv$fitted_swindow$add(mrv$fitted_view)
      mrv$fitted_swindow$setPolicy("automatic", "automatic")
      
      mrv$fitted_buffer <- mrv$fitted_view$getBuffer()
      mrv$fitted_buffer$setText(paste("\n",.writeFitted(mrv,"sens"),
                                      sep="",collapse = "\n"))
      fitted_tag_large <- mrv$fitted_buffer$createTag(tag.name="large", font="Monaco 10")
      
      fitted_iter <- mrv$fitted_buffer$getBounds()
      mrv$fitted_buffer$applyTag(fitted_tag_large, fitted_iter$start, fitted_iter$end)
      
      ##### bottom box to choose accuracy type
      fitted_bottom_frame <- gtkFrame("Accuracy type: ")
      
      fitted_bottom_radio <- gtkVBox()
      fitted_bottom_radiogp <- list()
      fitted_bottom_radiogp$sens <- gtkRadioButton(label = "sens")
      fitted_bottom_radiogp$spec <- gtkRadioButton(fitted_bottom_radiogp, label = "spec")
      fitted_bottom_radiogp$FPR <- gtkRadioButton(fitted_bottom_radiogp, label = "FPR")
      fitted_bottom_radiogp$FNR <- gtkRadioButton(fitted_bottom_radiogp, label = "FNR")
      fitted_bottom_radiogp$LRpos <- gtkRadioButton(fitted_bottom_radiogp, label = "LRpos")
      fitted_bottom_radiogp$LRneg <- gtkRadioButton(fitted_bottom_radiogp, label = "LRneg")
      fitted_bottom_radiogp$RD <- gtkRadioButton(fitted_bottom_radiogp, label = "RD")
      fitted_bottom_radiogp$LLRpos <- gtkRadioButton(fitted_bottom_radiogp, label = "LLRpos")
      fitted_bottom_radiogp$LLRneg <- gtkRadioButton(fitted_bottom_radiogp, label = "LLRneg")
      fitted_bottom_radiogp$LDOR <- gtkRadioButton(fitted_bottom_radiogp, label = "LDOR")
      fitted_bottom_radiogp$DOR <- gtkRadioButton(fitted_bottom_radiogp, label = "DOR")
      #sapply(fitted_bottom_radiogp, fitted_bottom_radio$packStart)
      #fitted_bottom_radio[[1]]$setActive(TRUE) 
      sapply(fitted_bottom_radiogp,'[',"active")
      #fitted_bottom_radio_align <- gtkAlignment(xalign = 0)
      #fitted_bottom_radio_align$add(fitted_bottom_radio)
      
      
      fitted_bottom_table <- gtkTable(rows=2,columns=6,homogeneous=FALSE)
      fitted_bottom_table$setColSpacings(5)
      fitted_bottom_table$setRowSpacings(3)
      
      fitted_bottom_table$attach(fitted_bottom_radiogp$sens,left.attach=0,1, top.attach=0,1,
                                 xoptions=c("expand","fill"),yoptions="")
      fitted_bottom_table$attach(fitted_bottom_radiogp$spec,left.attach=1,2, top.attach=0,1,
                                 xoptions=c("expand","fill"),yoptions="")
      fitted_bottom_table$attach(fitted_bottom_radiogp$FPR,left.attach=2,3, top.attach=0,1,
                                 xoptions=c("expand","fill"),yoptions="")
      fitted_bottom_table$attach(fitted_bottom_radiogp$FNR,left.attach=3,4, top.attach=0,1,
                                 xoptions=c("expand","fill"),yoptions="")
      fitted_bottom_table$attach(fitted_bottom_radiogp$LRpos,left.attach=4,5, top.attach=0,1,
                                 xoptions=c("expand","fill"),yoptions="")
      fitted_bottom_table$attach(fitted_bottom_radiogp$LRneg,left.attach=5,6, top.attach=0,1,
                                 xoptions=c("expand","fill"),yoptions="")
      fitted_bottom_table$attach(fitted_bottom_radiogp$RD,left.attach=0,1, top.attach=1,2,
                                 xoptions=c("expand","fill"),yoptions="")
      fitted_bottom_table$attach(fitted_bottom_radiogp$DOR,left.attach=1,2, top.attach=1,2,
                                 xoptions=c("expand","fill"),yoptions="")
      fitted_bottom_table$attach(fitted_bottom_radiogp$LLRpos,left.attach=2,3, top.attach=1,2,
                                 xoptions=c("expand","fill"),yoptions="")
      fitted_bottom_table$attach(fitted_bottom_radiogp$LLRneg,left.attach=3,4, top.attach=1,2,
                                 xoptions=c("expand","fill"),yoptions="")
      fitted_bottom_table$attach(fitted_bottom_radiogp$LDOR,left.attach=4,5, top.attach=1,2,
                                 xoptions=c("expand","fill"),yoptions="")
      
      fitted_bottom_frame$add(fitted_bottom_table)
      
      fitted_bottom_Button <- gtkButton("Change")
      fitted_bottom_Button['image'] <- gtkImage(stock.id="gtk-refresh",size="button")
      fitted_bottom_Button$setAlignment(xalign=0.5,yalign=0.5)
      fitted_bottom_Button$setImagePosition("left")
      
      fitted_ps_box <- gtkVBox(spacing=10)
      fitted_ps_box$packStart(transp_frame(),expand=T,fill=T)
      fitted_ps_box$packStart(fitted_bottom_Button,expand=F,fill=F)
      
      fitted_bottom_box <- gtkHBox(spacing=10)
      fitted_bottom_box["border-width"]=10
      fitted_bottom_box$packStart(fitted_bottom_frame,expand=F,fill=F)
      fitted_bottom_box$packStart(transp_frame(),expand=T,fill=T)
      fitted_bottom_box$packStart(fitted_ps_box,expand=F,fill=F)
      
      mrv$fitted_page <- gtkVBox()
      mrv$fitted_page$packStart(mrv$fitted_swindow, expand=T, fill=T, padding=0)
      mrv$fitted_page$packStart(fitted_bottom_box , expand=F, fill=F, padding=0)
      
      fitted_page_event = .pageEventsLabel("Fitted Result")
      mrv$notebook$insertPage(mrv$fitted_page,fitted_page_event)
      
      fitted_page_num = mrv$notebook$pageNum(mrv$fitted_page)
      
      gSignalConnect(fitted_page_event, "button-press-event", function(...) {
        mrv$notebook$setCurrentPage(fitted_page_num)
        return(TRUE)
      })
      
      mrv$fitted.type = "sens"
      sapply(fitted_bottom_radiogp, gSignalConnect, "toggled",
             f = function(button, ...){
               if(button['active']){
                 mrv$fitted.type = button$getLabel()
               } 
             })
      
      gSignalConnect(fitted_bottom_Button, "pressed", function(...) {
        mrv$fitted_buffer <- mrv$fitted_view$getBuffer()
        mrv$fitted_buffer$setText(paste("\n",.writeFitted(mrv,mrv$fitted.type),
                                        sep="",collapse = "\n"))
        fitted_iter <- mrv$fitted_buffer$getBounds()
        mrv$fitted_buffer$applyTag(fitted_tag_large, fitted_iter$start, fitted_iter$end)
      })
      
      fitted_page_num = mrv$notebook$pageNum(mrv$fitted_page)
      
      
      ###################################################
      ## construct sroc page
      ###################################################
      mrv$sroc_plot <- gtkDrawingAreaNew()
      mrv$sroc_plot$setSizeRequest(400,400)
      
      sroc_align <- gtkAlignment(xalign = 0.5, yalign = 0.5, yscale = 0, xscale = 0)
      sroc_align$add(mrv$sroc_plot)
      
      sroc_box <- gtkHBox()
      sroc_box$packStart(sroc_align, expand=T, fill=T, padding=0)
      
      sroc_outbox <- gtkVBox()
      sroc_outbox$packStart(sroc_box, expand=T, fill=T, padding=0)
      
      mrv$sroc_window <- .scrolledWindowWithViewport(sroc_outbox)
      ########### main
      sroc_main_label <- bottomLabel("main: ",1)
      sroc_main_entry <- bottomEntry("SROC Plot", 15, 100)
      
      sorc_main_size_label <- bottomLabel("cex.main: ",1)
      sroc_main_size <- gtkSpinButtonNewWithRange(min=1,max=5,step=0.1)
      sroc_main_size$setValue(1.5)
      
      sroc_axis_label <- bottomLabel("cex.axis: ",1)
      sroc_axis_size <- gtkSpinButtonNewWithRange(min=1,max=5,step=0.1)
      sroc_axis_size$setValue(1)
      
      sroc_main_table <- gtkTable(rows=2,columns=2,homogeneous=FALSE)
      sroc_main_table$setColSpacings(5)
      sroc_main_table$setRowSpacings(3)
      
      sroc_main_table$attach(sroc_main_label,left.attach=0,1, top.attach=0,1,
                             xoptions=c("expand","fill"),yoptions="")
      sroc_main_table$attach(sroc_main_entry,left.attach=1,2, top.attach=0,1,
                             xoptions="",yoptions="")
      sroc_main_table$attach(sorc_main_size_label,left.attach=0,1, top.attach=1,2,
                             xoptions=c("expand","fill"),yoptions="")
      sroc_main_table$attach(sroc_main_size,left.attach=1,2, top.attach=1,2,
                             xoptions=c("expand","fill"),yoptions="")
      sroc_main_table$attach(sroc_axis_label,left.attach=0,1, top.attach=2,3,
                             xoptions=c("expand","fill"),yoptions="")
      sroc_main_table$attach(sroc_axis_size,left.attach=1,2, top.attach=2,3,
                             xoptions=c("expand","fill"),yoptions="")
      
      sroc_main_box <- gtkVBox()
      sroc_main_box["border-width"]=10
      sroc_main_box$packStart(sroc_main_table,expand=T,fill=F)
      
      sroc_main_frame <- gtkFrame("Title")
      sroc_main_frame$add(sroc_main_box)
      
      ########### label
      sorc_xlab_label <- bottomLabel("xlab: ",1)
      sroc_xlab_entry <- bottomEntry(mrv$est$names.fitted[2], 15, 20)
      
      sorc_ylab_label <- bottomLabel("ylab: ",1)
      sroc_ylab_entry <- bottomEntry(mrv$est$names.fitted[1], 15, 20)
      
      sorc_lab_size_label <- bottomLabel("cex: ",1)
      sroc_lab_size <- gtkSpinButtonNewWithRange(min=0.1,max=5,step=0.1)
      sroc_lab_size$setValue(1)
      
      sroc_lab_table <- gtkTable(rows=3,columns=2,homogeneous=FALSE)
      sroc_lab_table$setColSpacings(5)
      sroc_lab_table$setRowSpacings(3)
      
      sroc_lab_table$attach(sorc_xlab_label,left.attach=0,1, top.attach=0,1,
                            xoptions=c("expand","fill"),yoptions="")
      sroc_lab_table$attach(sroc_xlab_entry,left.attach=1,2, top.attach=0,1,
                            xoptions="",yoptions="")
      sroc_lab_table$attach(sorc_ylab_label,left.attach=0,1, top.attach=1,2,
                            xoptions=c("expand","fill"),yoptions="")
      sroc_lab_table$attach(sroc_ylab_entry,left.attach=1,2, top.attach=1,2,
                            xoptions="",yoptions="")
      sroc_lab_table$attach(sorc_lab_size_label,left.attach=0,1, top.attach=2,3,
                            xoptions=c("expand","fill"),yoptions="")
      sroc_lab_table$attach(sroc_lab_size,left.attach=1,2, top.attach=2,3,
                            xoptions=c("expand","fill"),yoptions="")
      
      sroc_lab_box <- gtkVBox()
      sroc_lab_box["border-width"]=10
      sroc_lab_box$packStart(sroc_lab_table,expand=T,fill=F)
      
      sroc_label_frame <- gtkFrame("Label")
      sroc_label_frame["border-width"]=0
      sroc_label_frame$add(sroc_lab_box)
      
      ########### save option
      sroc_form_label <- bottomLabel("form: ",1)
      sroc_form_combo <- gtkComboBoxNewText()
      sroc_form_combo["width-request"] <- 100
      sapply(c(".pdf",".eps",".png",".jpg"), sroc_form_combo$appendText)
      sroc_form_combo$setActive(0)
      
      mrv$saveform = ".pdf"
      gSignalConnect(sroc_form_combo, "changed", function(button, ...) {
        mrv$saveform = gtkComboBoxGetActiveText(button)
      })
      
      sroc_width_label <- bottomLabel("width: ",1)
      sroc_width_size <- gtkSpinButtonNewWithRange(min=2,max=20,step=1)
      sroc_width_size$setValue(5)
      
      sroc_height_label <- bottomLabel("height: ",1)
      sroc_height_size <- gtkSpinButtonNewWithRange(min=2,max=20,step=1)
      sroc_height_size$setValue(5)
      
      sroc_so_table <- gtkTable(rows=3,columns=2,homogeneous=FALSE)
      sroc_so_table$setColSpacings(5)
      sroc_so_table$setRowSpacings(3)
      
      sroc_so_table$attach(sroc_form_label,left.attach=0,1, top.attach=0,1,
                           xoptions=c("expand","fill"),yoptions="")
      sroc_so_table$attach(sroc_form_combo,left.attach=1,2, top.attach=0,1,
                           xoptions=c("expand","fill"),yoptions="")
      sroc_so_table$attach(sroc_width_label,left.attach=0,1, top.attach=1,2,
                           xoptions=c("expand","fill"),yoptions="")
      sroc_so_table$attach(sroc_width_size,left.attach=1,2, top.attach=1,2,
                           xoptions=c("expand","fill"),yoptions="")
      sroc_so_table$attach(sroc_height_label,left.attach=0,1, top.attach=2,3,
                           xoptions=c("expand","fill"),yoptions="")
      sroc_so_table$attach(sroc_height_size,left.attach=1,2, top.attach=2,3,
                           xoptions=c("expand","fill"),yoptions="")
      
      sroc_so_box <- gtkVBox()
      sroc_so_box["border-width"]=10
      sroc_so_box$packStart(sroc_so_table,expand=T,fill=F)
      
      sroc_so_frame <- gtkFrame("Save Options")
      sroc_so_frame["border-width"]=0
      sroc_so_frame$add(sroc_so_box)
      
      ##### plot and save button
      sroc_plotbt <- plot_button()
      sroc_savebt <- save_button()
      
      sroc_ps_box <- gtkVBox(spacing=10)
      sroc_ps_box$packStart(transp_frame(),expand=T,fill=T)
      sroc_ps_box$packStart(sroc_plotbt,expand=F,fill=F)
      sroc_ps_box$packStart(sroc_savebt,expand=F,fill=F)
      
      sroc_bottom_box <- gtkHBox(spacing=10)
      sroc_bottom_box["border-width"]=10
      sroc_bottom_box$packStart(sroc_main_frame,expand=F,fill=F)
      sroc_bottom_box$packStart(sroc_label_frame,expand=F,fill=F)
      sroc_bottom_box$packStart(sroc_so_frame,expand=F,fill=F)
      sroc_bottom_box$packStart(transp_frame(),expand=T,fill=T)
      sroc_bottom_box$packStart(sroc_ps_box,expand=F,fill=F)
      
      mrv$sroc_page <- gtkVBox()
      mrv$sroc_page$packStart(mrv$sroc_window, expand=T, fill=T, padding=0)
      mrv$sroc_page$packStart(sroc_bottom_box , expand=F, fill=F, padding=0)
      
      sroc_page_event = .pageEventsLabel("SROC Plot")      
      mrv$notebook$appendPage(mrv$sroc_page, sroc_page_event)
      
      col.collect = c("aquamarine","blueviolet","brown","red","green","darkblue",
                      "cyan","chocolate","chartreuse","darkgoldenrod","deepskyblue",
                      "darkorange","darkorchid","deeppink","hotpink","magenta")
      
      if(!is.null(mrv$modality)){
        nl = length(mrv$partialdata)
        if(nl!=1){
          ad.col = col.collect[1:(nl-1)]
          col = c(mrv$cr.col, ad.col)
          mrv$cr.col = col
          mrv$sp.col = col
          mrv$pr.col = col
          mrv$line.col = col
        }
      }
      
      Sys.sleep(.1) 
      asCairoDevice(mrv$sroc_plot)
      SROC(x=mrv$est, est.type=mrv$est.type, sp.cex=mrv$sp.cex, sp.pch=mrv$sp.pch, sp.col=mrv$sp.col,
           dataShow=mrv$dataShow, data.col=mrv$data.col, data.cex=mrv$data.cex, data.pch=mrv$data.pch, 
           lineShow=mrv$lineShow, sroc.type=mrv$sroc.type, line.lty=mrv$line.lty, line.lwd=mrv$line.lwd, line.col=mrv$line.col,
           crShow=mrv$crShow, cr.lty=mrv$cr.lty, cr.lwd=mrv$cr.lwd, cr.col=mrv$cr.col,
           prShow=mrv$prShow, pr.lty=mrv$pr.lty, pr.lwd=mrv$pr.lwd,  pr.col=mrv$pr.col,
           dataFit = T, add=FALSE)
      
      sroc_page_num = mrv$notebook$pageNum(mrv$sroc_page)
      mrv$notebook$setCurrentPage(sroc_page_num)
      
      mrv$sidenote$setCurrentPage(mrv$sidesroc_page_num)
      
      gSignalConnect(sroc_page_event, "button-press-event", function(...) {
        mrv$notebook$setCurrentPage(sroc_page_num)
        mrv$sidenote$setCurrentPage(mrv$sidesroc_page_num)
        return(TRUE)
      })
      
      gSignalConnect(sroc_plotbt, "pressed", function(...) {
        mrv$sroc_main_text = gtkEntryGetText(sroc_main_entry)
        mrv$sroc_main_cex = sroc_main_size$getValue()
        mrv$sroc_axis_cex = sroc_axis_size$getValue()
        mrv$sroc_xlabel = gtkEntryGetText(sroc_xlab_entry)
        mrv$sroc_ylabel = gtkEntryGetText(sroc_ylab_entry)
        mrv$sroc_label_cex = sroc_lab_size$getValue()
        Sys.sleep(.1) 
        asCairoDevice(mrv$sroc_plot)
        SROC(x=mrv$est, est.type=mrv$est.type, sp.cex=mrv$sp.cex, sp.pch=mrv$sp.pch, sp.col=mrv$sp.col,
             dataShow=mrv$dataShow, data.col=mrv$data.col, data.cex=mrv$data.cex, data.pch=mrv$data.pch, 
             lineShow=mrv$lineShow, sroc.type=mrv$sroc.type, line.lty=mrv$line.lty, line.lwd=mrv$line.lwd, line.col=mrv$line.col,
             crShow=mrv$crShow, cr.lty=mrv$cr.lty, cr.lwd=mrv$cr.lwd, cr.col=mrv$cr.col,
             prShow=mrv$prShow, pr.lty=mrv$pr.lty, pr.lwd=mrv$pr.lwd,  pr.col=mrv$pr.col,
             dataFit = T, add=FALSE,
             main=mrv$sroc_main_text, cex.main=mrv$sroc_main_cex, cex.axis=mrv$sroc_axis_cex,
             xlab=mrv$sroc_xlabel, ylab=mrv$sroc_ylabel, cex.lab=mrv$sroc_label_cex)
      })
      
      gSignalConnect(sroc_savebt, "pressed", function(...) {
        dialog <- gtkFileChooserDialog("Enter a name for the file", window,
                                       "save", "gtk-cancel", GtkResponseType["cancel"], "gtk-save",
                                       GtkResponseType["accept"])
        if (dialog$run() == GtkResponseType["accept"]){
          mrv$sroc_main_text = gtkEntryGetText(sroc_main_entry)
          mrv$sroc_main_cex = sroc_main_size$getValue()
          mrv$sroc_axis_cex = sroc_axis_size$getValue()
          mrv$sroc_xlabel = gtkEntryGetText(sroc_xlab_entry)
          mrv$sroc_ylabel = gtkEntryGetText(sroc_ylab_entry)
          mrv$sroc_label_cex = sroc_lab_size$getValue()
          filepath <- dialog$getFilename()
          fileform <- mrv$saveform
          filename <- paste(filepath, fileform, sep="")
          width <- sroc_width_size$getValue()
          height <- sroc_height_size$getValue()
          if(fileform==".eps"){
            cairo_ps(filename = filename, width=width, height=height)
          }else if(fileform==".pdf"){
            pdf(filename, width=width, height=height,...)
          }else if(fileform==".jpg"){
            jpeg(filename = filename,
                 width = width, height = height, units = "in", res=300)
          }else{
            png(filename = filename,
                width = width, height = height, units = "in", res=300)
          }
          SROC(x=mrv$est, est.type=mrv$est.type, sp.cex=mrv$sp.cex, sp.pch=mrv$sp.pch, sp.col=mrv$sp.col,
               dataShow=mrv$dataShow, data.col=mrv$data.col, data.cex=mrv$data.cex, data.pch=mrv$data.pch, 
               lineShow=mrv$lineShow, sroc.type=mrv$sroc.type, line.lty=mrv$line.lty, line.lwd=mrv$line.lwd, line.col=mrv$line.col,
               crShow=mrv$crShow, cr.lty=mrv$cr.lty, cr.lwd=mrv$cr.lwd, cr.col=mrv$cr.col,
               prShow=mrv$prShow, pr.lty=mrv$pr.lty, pr.lwd=mrv$pr.lwd,  pr.col=mrv$pr.col,
               dataFit = T, add=FALSE,
               main=mrv$sroc_main_text, cex.main=mrv$sroc_main_cex, cex.axis=mrv$sroc_axis_cex,
               xlab=mrv$sroc_xlabel, ylab=mrv$sroc_ylabel, cex.lab=mrv$sroc_label_cex)
          dev.off()
        }
        dialog$destroy()
      })
      
      mrv$statusbar$push(mrv$info, paste("SROC plot is done."))
      
      ###################################################
      ## construct forest page
      ###################################################
      mrv$forest_plot <- gtkDrawingAreaNew()
      mrv$forest_plot$setSizeRequest(400, 400)
      
      forest_box <- gtkHBox()
      forest_box["border-width"]=10
      forest_box$packStart(mrv$forest_plot, expand=T, fill=T, padding=0)
      
      mrv$forest_window <- .scrolledWindowWithViewport(forest_box)
      
      ########### main
      forest_main_label <- bottomLabel("main: ",1)
      forest_main_entry <- bottomEntry("Forest Plot", 15, 100)
      
      forest_main_size_label <- bottomLabel("cex.main: ",1)
      forest_main_size <- gtkSpinButtonNewWithRange(min=1,max=5,step=0.1)
      forest_main_size$setValue(1.5)
      
      forest_axis_label <- bottomLabel("cex.axis: ",1)
      forest_axis_size <- gtkSpinButtonNewWithRange(min=1,max=5,step=0.1)
      forest_axis_size$setValue(1)
      
      forest_main_table <- gtkTable(rows=2,columns=2,homogeneous=FALSE)
      forest_main_table$setColSpacings(5)
      forest_main_table$setRowSpacings(3)
      
      forest_main_table$attach(forest_main_label,left.attach=0,1, top.attach=0,1,
                               xoptions=c("expand","fill"),yoptions="")
      forest_main_table$attach(forest_main_entry,left.attach=1,2, top.attach=0,1,
                               xoptions="",yoptions="")
      forest_main_table$attach(forest_main_size_label,left.attach=0,1, top.attach=1,2,
                               xoptions=c("expand","fill"),yoptions="")
      forest_main_table$attach(forest_main_size,left.attach=1,2, top.attach=1,2,
                               xoptions=c("expand","fill"),yoptions="")
      forest_main_table$attach(forest_axis_label,left.attach=0,1, top.attach=2,3,
                               xoptions=c("expand","fill"),yoptions="")
      forest_main_table$attach(forest_axis_size,left.attach=1,2, top.attach=2,3,
                               xoptions=c("expand","fill"),yoptions="")
      
      forest_main_box <- gtkVBox()
      forest_main_box["border-width"]=10
      forest_main_box$packStart(forest_main_table,expand=T,fill=F)
      
      forest_main_frame <- gtkFrame("Title")
      forest_main_frame$add(forest_main_box)
      
      ########### save option
      forest_form_label <- bottomLabel("form: ",1)
      forest_form_combo <- gtkComboBoxNewText()
      forest_form_combo["width-request"] <- 100
      sapply(c(".pdf",".eps",".png",".jpg"), forest_form_combo$appendText)
      forest_form_combo$setActive(0)
      
      mrv$saveform_forest = ".pdf"
      gSignalConnect(forest_form_combo, "changed", function(button, ...) {
        mrv$saveform_forest = gtkComboBoxGetActiveText(button)
      })
      
      forest_width_label <- bottomLabel("width: ",1)
      forest_width_size <- gtkSpinButtonNewWithRange(min=2,max=20,step=1)
      forest_width_size$setValue(5)
      
      forest_height_label <- bottomLabel("height: ",1)
      forest_height_size <- gtkSpinButtonNewWithRange(min=2,max=20,step=1)
      forest_height_size$setValue(5)
      
      forest_so_table <- gtkTable(rows=3,columns=2,homogeneous=FALSE)
      forest_so_table$setColSpacings(5)
      forest_so_table$setRowSpacings(3)
      
      forest_so_table$attach(forest_form_label,left.attach=0,1, top.attach=0,1,
                             xoptions=c("expand","fill"),yoptions="")
      forest_so_table$attach(forest_form_combo,left.attach=1,2, top.attach=0,1,
                             xoptions=c("expand","fill"),yoptions="")
      forest_so_table$attach(forest_width_label,left.attach=0,1, top.attach=1,2,
                             xoptions=c("expand","fill"),yoptions="")
      forest_so_table$attach(forest_width_size,left.attach=1,2, top.attach=1,2,
                             xoptions=c("expand","fill"),yoptions="")
      forest_so_table$attach(forest_height_label,left.attach=0,1, top.attach=2,3,
                             xoptions=c("expand","fill"),yoptions="")
      forest_so_table$attach(forest_height_size,left.attach=1,2, top.attach=2,3,
                             xoptions=c("expand","fill"),yoptions="")
      
      forest_so_box <- gtkVBox()
      forest_so_box["border-width"]=10
      forest_so_box$packStart(forest_so_table,expand=T,fill=F)
      
      forest_so_frame <- gtkFrame("Save Options")
      forest_so_frame["border-width"]=0
      forest_so_frame$add(forest_so_box)
      
      ##### plot and save button
      forest_plotbt <- plot_button()
      forest_savebt <- save_button()
      
      forest_ps_box <- gtkVBox(spacing=10)
      forest_ps_box$packStart(transp_frame(),expand=T,fill=T)
      forest_ps_box$packStart(forest_plotbt,expand=F,fill=F)
      forest_ps_box$packStart(forest_savebt,expand=F,fill=F)
      
      forest_bottom_box <- gtkHBox(spacing=10)
      forest_bottom_box["border-width"]=10
      forest_bottom_box$packStart(forest_main_frame,expand=F,fill=F)
      forest_bottom_box$packStart(forest_so_frame,expand=F,fill=F)
      forest_bottom_box$packStart(transp_frame(),expand=T,fill=T)
      forest_bottom_box$packStart(forest_ps_box,expand=F,fill=F)
      
      mrv$forest_page <- gtkVBox()
      mrv$forest_page$packStart(mrv$forest_window, expand=T, fill=T, padding=0)
      mrv$forest_page$packStart(forest_bottom_box , expand=F, fill=F, padding=0)
      
      forest_page_event = .pageEventsLabel("Forest Plot")     
      mrv$notebook$appendPage(mrv$forest_page, forest_page_event)
      
      Sys.sleep(.1) 
      asCairoDevice(mrv$forest_plot)
      forest(x = mrv$est, accuracy.type=mrv$accuracy, est.type=mrv$forest.est.type, 
             p.cex=mrv$forest.pcex, p.pch=mrv$forest.pch, p.col=mrv$forest.pcol,
             nameShow=mrv$forest.nameShow, dataShow=mrv$forest.dataShow, ciShow=mrv$forest.ciShow,
             text.cex=mrv$forest.text.size, shade.col=mrv$forest.scol, arrow.col="black",
             main="Forest Plot", main.cex=1.5, cut=TRUE, intervals=c(0.025,0.975))
      
      forest_page_num = mrv$notebook$pageNum(mrv$forest_page)
      
      gSignalConnect(forest_page_event, "button-press-event", function(...) {
        mrv$notebook$setCurrentPage(forest_page_num)
        mrv$sidenote$setCurrentPage(mrv$sideforest_page_num)
        return(TRUE)
      })
      
      gSignalConnect(forest_plotbt, "pressed", function(...) {
        mrv$forest_main_text = gtkEntryGetText(forest_main_entry)
        mrv$forest_main_cex = forest_main_size$getValue()
        mrv$forest_axis_cex = forest_axis_size$getValue()
        Sys.sleep(.1) 
        asCairoDevice(mrv$forest_plot)
        forest(x = mrv$est, accuracy.type=mrv$accuracy, est.type=mrv$forest.est.type, 
               p.cex=mrv$forest.pcex, p.pch=mrv$forest.pch, p.col=mrv$forest.pcol,
               nameShow=mrv$forest.nameShow, dataShow=mrv$forest.dataShow, ciShow=mrv$forest.ciShow,
               text.cex=mrv$forest.text.size, shade.col=mrv$forest.scol, arrow.col="black",
               cut=TRUE, intervals=c(0.025,0.975),
               main=mrv$forest_main_text, main.cex=mrv$forest_main_cex, cex.axis=mrv$forest_axis_cex)
      })
      
      gSignalConnect(forest_savebt, "pressed", function(...) {
        dialog <- gtkFileChooserDialog("Enter a name for the file", window,
                                       "save", "gtk-cancel", GtkResponseType["cancel"], "gtk-save",
                                       GtkResponseType["accept"])
        if (dialog$run() == GtkResponseType["accept"]){
          mrv$forest_main_text = gtkEntryGetText(forest_main_entry)
          mrv$forest_main_cex = forest_main_size$getValue()
          mrv$forest_axis_cex = forest_axis_size$getValue()
          filepath <- dialog$getFilename()
          fileform <- mrv$saveform
          filename <- paste(filepath, fileform, sep="")
          width <- forest_width_size$getValue()
          height <- forest_height_size$getValue()
          if(fileform==".eps"){
            cairo_ps(filename = filename, width=width, height=height)
          }else if(fileform==".pdf"){
            pdf(filename, width=width, height=height,...)
          }else if(fileform==".jpg"){
            jpeg(filename = filename,
                 width = width, height = height, units = "in", res=300)
          }else{
            png(filename = filename,
                width = width, height = height, units = "in", res=300)
          }
          forest(x = mrv$est, accuracy.type=mrv$accuracy, est.type=mrv$forest.est.type, 
                 p.cex=mrv$forest.pcex, p.pch=mrv$forest.pch, p.col=mrv$forest.pcol,
                 nameShow=mrv$forest.nameShow, dataShow=mrv$forest.dataShow, ciShow=mrv$forest.ciShow,
                 text.cex=mrv$forest.text.size, shade.col=mrv$forest.scol, arrow.col="black",
                 cut=TRUE, intervals=c(0.025,0.975),
                 main=mrv$forest_main_text, main.cex=mrv$forest_main_cex, cex.axis=mrv$forest_axis_cex)
          dev.off()
        }
        dialog$destroy()
      })
      
      mrv$statusbar$push(mrv$info, paste("Forest plot is done."))
      
      ###################################################
      ## construct R-code page
      ###################################################
      mrv$rcode_page = gtkFrame()
      mrv$rcode_page["border-width"] = 0
      mrv$rcode_view <- gtkTextView()
      mrv$rcode_view["left-margin"] = 10
      mrv$rcode_sw <- gtkScrolledWindow()
      mrv$rcode_sw["shadow-type"] = "none"
      mrv$rcode_sw$add(mrv$rcode_view)
      mrv$rcode_sw$setPolicy("automatic", "automatic")
      
      mrv$rcode_buffer <- mrv$rcode_view$getBuffer()
      mrv$rcode_buffer$setText(paste("\n",.writeRcode(mrv),
                                     sep="",collapse = "\n"))
      
      rcode_tag_large <- mrv$rcode_buffer$createTag(tag.name="large", font="12")
      
      rcode_iter <- mrv$rcode_buffer$getBounds()
      mrv$rcode_buffer$applyTag(rcode_tag_large, rcode_iter$start, rcode_iter$end)
      
      mrv$rcode_page$add(mrv$rcode_sw)
      rcode_page_event = .pageEventsLabel("R-code")
      mrv$notebook$insertPage(mrv$rcode_page,rcode_page_event)
      
      rcode_page_num = mrv$notebook$pageNum(mrv$rcode_page)
      
      gSignalConnect(rcode_page_event, "button-press-event", function(...) {
        mrv$notebook$setCurrentPage(rcode_page_num)
        return(TRUE)
      })
      
      mrv$statusbar$push(mrv$info, paste("R-code is written."))
      
      ###################################################
      ## construct Logfile page
      ###################################################
      mrv$logfile_page = gtkFrame()
      mrv$logfile_page["border-width"] = 0
      mrv$logfile_view <- gtkTextView()
      mrv$logfile_view["left-margin"] = 10
      mrv$logfile_sw <- gtkScrolledWindow()
      mrv$logfile_sw["shadow-type"] = "none"
      mrv$logfile_sw$add(mrv$logfile_view)
      mrv$logfile_sw$setPolicy("automatic", "automatic")
      
      mrv$logfile_buffer <- mrv$logfile_view$getBuffer()
      mrv$logfile_buffer$setText(paste(mrv$est$inla.result$logfile,sep="",collapse = "\n"))
      
      mrv$logfile_page$add(mrv$logfile_sw)
      logfile_page_event = .pageEventsLabel("Logfile")
      mrv$notebook$insertPage(mrv$logfile_page,logfile_page_event)
      
      logfile_page_num = mrv$notebook$pageNum(mrv$logfile_page)
      
      gSignalConnect(logfile_page_event, "button-press-event", function(...) {
        mrv$notebook$setCurrentPage(logfile_page_num)
        return(TRUE)
      })
      
      mrv$statusbar$push(mrv$info, paste("Logfile is written."))
      
      ###################################################
      ## resize window signal
      ###################################################
      gSignalConnect(mrv$main_window, "check-resize", function(window, ...) {
        theSize <- mrv$sroc_window$getAllocation()$allocation
        size <- min(theSize$width,theSize$height)
        mrv$sroc_plot$setSizeRequest(size,size)
        
        forestSize <- mrv$forest_window$getAllocation()$allocation
        mrv$forest_plot$setSizeRequest(1.5*forestSize$width,forestSize$height)
      })
    }
  }
  }
}