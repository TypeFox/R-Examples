meta4diagGUI <- function(){
  # construct new enviroment **** very important ***
  if(!exists("mrv")){
    mrv <<- new.env() # meta4diag R-GUI environment
  }
  if(!is.null(names(mrv$main_window))){
    mrv$main_window$destroy()
  }
  mrv <<- new.env()
  
  #if (.Platform$OS.type == "windows"){
  #}else{}
  
  mrv$VERSION <- "1.0.20"
  mrv$DATE <- "2015-08-24"
  mrv$COPYRIGHT <- "Copyright (C) 2015 INLA Group."


  ########## define some functions
  mrv$search_icon = system.file("image/search.png", package="meta4diag")
  mrv$refresh_icon = system.file("image/refresh.png",package="meta4diag")
  welcome_text = readLines(system.file("etc/welcometext.txt",package="meta4diag"), n = -1L, warn = FALSE)
  
  sidebar_vector = c("prior", "model", "data", "sroc", "forest", "funnel")
  mrv$sidebar_icon =lapply(sidebar_vector, function(x){
    system.file(paste("image/sidebar_",x,"_title.gif",sep=""),package="meta4diag")
  })
  
  ####### trick the globle variable
  #shape = rate = xlim = u = variance = nu = R11 = R22 = R12 = U_min = a1 = U_max = a2 = a = b = NULL
  ####### Line Width Combo
  makePixbufForLineWidth <- function(Width){
    filename <- tempfile()
    png(filename=filename, width = 60, height = 20, units = "px", pointsize = 12,
        bg = "transparent",  res = NA,
        type = c("cairo", "cairo-png", "Xlib", "quartz"))
    grid.newpage()
    grid.draw(grid.lines(x = unit(c(0:20),"npc"),
                         y = unit(c(0.5),"npc"),
                         gp=gpar(lty=1,lwd=Width)))
    dev.off()
    image <- gdkPixbufNewFromFile(filename)
    unlink(filename)
    return(image$retval)
  }
  
  gtkLineWidthComboBox <- function(init.active){
    lwd_vector <- as.character(1:6)
    store <- gtkListStore(c("GObject","gchararray"))
    for(i in 1:length(lwd_vector)) {
      iter <- store$Append()
      store$SetValue(iter$iter, 0,  makePixbufForLineWidth(i))
      store$SetValue(iter$iter, 1, lwd_vector[i])
    }
    combo <- gtkComboBox(model=store)
    crp = gtkCellRendererPixbuf()
    crp['xalign'] <- 0
    combo$packStart(crp, expand=FALSE)
    combo$addAttribute(crp, "pixbuf", 0)
    crt = gtkCellRendererText()
    crt['xalign'] <- 1
    combo$packStart(crt)
    combo$addAttribute(crt, "text", 1)
    combo$setActive(init.active)
    return(combo)
  } 
  ####### Line Type Combo
  makePixbufForLineType <- function(Type){
    filename <- tempfile()
    png(filename=filename, width = 60, height = 20, units = "px", pointsize = 12,
        bg = "transparent",  res = NA,
        type = c("cairo", "cairo-png", "Xlib", "quartz"))
    grid.newpage()
    grid.draw(grid.lines(x = unit(c(0:20),"npc"),
                         y = unit(c(0.5),"npc"),
                         gp=gpar(lty=Type,lwd=2)))
    dev.off()
    image <- gdkPixbufNewFromFile(filename)
    unlink(filename)
    return(image$retval)
  }
  
  gtkLineTypeComboBox <- function(init.active){
    lty_vector <- as.character(1:6)
    store <- gtkListStore(c("GObject","gchararray"))
    for(i in 1:length(lty_vector)) {
      iter <- store$Append()
      store$SetValue(iter$iter, 0,  makePixbufForLineType(i))
      store$SetValue(iter$iter, 1, lty_vector[i])
    }
    combo <- gtkComboBox(model=store)
    crp = gtkCellRendererPixbuf()
    crp['xalign'] <- 0
    combo$packStart(crp, expand=FALSE)
    combo$addAttribute(crp, "pixbuf", 0)
    crt = gtkCellRendererText()
    crt['xalign'] <- 1
    combo$packStart(crt)
    combo$addAttribute(crt, "text", 1)
    combo$setActive(init.active)
    return(combo)
  } 
  ####### Points symbol Combo
  makePixbufForPointsType <- function(Type){
    filename <- tempfile()
    png(filename=filename, width = 60, height = 20, units = "px", 
        bg = "transparent",  res = NA,
        type = c("cairo", "cairo-png", "Xlib", "quartz"))
    grid.newpage()
    if(Type %in% c("@","+","%","#","*","o","O")){
      grid.draw(grid.points(x = unit(c(0.2,0.5,0.8),"npc"), y = unit(rep(0.5,3),"npc"), pch = Type,gp = gpar(ylim=c(0,1))))
    }else{
      grid.draw(grid.points(x = unit(c(0.2,0.5,0.8),"npc"), y = unit(rep(0.5,3),"npc"), pch = as.numeric(Type),gp = gpar(ylim=c(0,1))))
    }  
    dev.off()
    image <- gdkPixbufNewFromFile(filename)
    unlink(filename)
    return(image$retval)
  }
  
  gtkPointsTypeComboBox <- function(init.active){
    pch_vector <- c(as.character(1:25),"@","+","%","#","*","o","O")
    store <- gtkListStore(c("GObject","gchararray"))
    for(i in 1:length(pch_vector)){
      iter <- store$Append()
      store$SetValue(iter$iter, 0,  makePixbufForPointsType(pch_vector[i]))
      store$SetValue(iter$iter, 1, pch_vector[i])
    }
    combo <- gtkComboBox(model=store)
    crp = gtkCellRendererPixbuf()
    crp['xalign'] <- 0
    combo$packStart(crp, expand=FALSE)
    combo$addAttribute(crp, "pixbuf", 0)
    crt = gtkCellRendererText()
    crt['xalign'] <- 1
    combo$packStart(crt)
    combo$addAttribute(crt, "text", 1)
    combo$setActive(init.active)
    return(combo)
  } 
  
  ####### Points Size Combo
  makePixbufForPointsSize <- function(Size){
    filename <- tempfile()
    png(filename=filename, width = 60, height = 20, units = "px", 
        bg = "transparent",  res = NA,
        type = c("cairo", "cairo-png", "Xlib", "quartz"))
    grid.newpage()
    if(Size=="scaled"){
      grid.draw(grid.points(x = unit(c(0.2,0.5,0.8),"npc"), y = unit(rep(0.5,3),"npc"), pch = 1,gp = gpar(ylim=c(0,1),cex=c(0.3,0.8,1.5))))
    }else if(Size=="bubble"){
      grid.draw(grid.points(x = unit(c(0.2,0.5,0.8),"npc"), y = unit(rep(0.5,3),"npc"), pch = 16,gp = gpar(ylim=c(0,1),cex=c(0.3,0.8,1.5),col=rep("green",3))))
    }else{
      grid.draw(grid.points(x = unit(c(0.2,0.5,0.8),"npc"), y = unit(rep(0.5,3),"npc"), pch = 1,gp = gpar(ylim=c(0,1),cex=rep(as.numeric(Size),3))))
    }
    dev.off()
    image <- gdkPixbufNewFromFile(filename)
    unlink(filename)
    return(image$retval)
  }
  
  gtkPointsSizeComboBox <- function(init.active){
    pcex_vector <- c("scaled", as.character(seq(0.5,2.5,by=0.5)))
    store <- gtkListStore(c("GObject","gchararray"))
    for(i in 1:length(pcex_vector)){
      iter <- store$Append()
      store$SetValue(iter$iter, 0,  makePixbufForPointsSize(pcex_vector[i]))
      store$SetValue(iter$iter, 1, pcex_vector[i])
    }
    combo <- gtkComboBox(model=store)
    crp = gtkCellRendererPixbuf()
    crp['xalign'] <- 0
    combo$packStart(crp, expand=FALSE)
    combo$addAttribute(crp, "pixbuf", 0)
    crt = gtkCellRendererText()
    crt['xalign'] <- 1
    combo$packStart(crt)
    combo$addAttribute(crt, "text", 1)
    combo$setActive(init.active)
    return(combo)
  } 
  
  gtkPointsSizeComboBox_bubble <- function(init.active){
    pcex_vector <- c("bubble","scaled", as.character(seq(0.5,2.5,by=0.5)))
    store <- gtkListStore(c("GObject","gchararray"))
    for(i in 1:length(pcex_vector)){
      iter <- store$Append()
      store$SetValue(iter$iter, 0,  makePixbufForPointsSize(pcex_vector[i]))
      store$SetValue(iter$iter, 1, pcex_vector[i])
    }
    combo <- gtkComboBox(model=store)
    crp = gtkCellRendererPixbuf()
    crp['xalign'] <- 0
    combo$packStart(crp, expand=FALSE)
    combo$addAttribute(crp, "pixbuf", 0)
    crt = gtkCellRendererText()
    crt['xalign'] <- 1
    combo$packStart(crt)
    combo$addAttribute(crt, "text", 1)
    combo$setActive(init.active)
    return(combo)
  } 
  
  gtkPointsSizeComboBox_simple <- function(init.active){
    pcex_vector <- as.character(seq(0.5,2.5,by=0.5))
    store <- gtkListStore(c("GObject","gchararray"))
    for(i in 1:length(pcex_vector)){
      iter <- store$Append()
      store$SetValue(iter$iter, 0,  makePixbufForPointsSize(pcex_vector[i]))
      store$SetValue(iter$iter, 1, pcex_vector[i])
    }
    combo <- gtkComboBox(model=store)
    crp = gtkCellRendererPixbuf()
    crp['xalign'] <- 0
    combo$packStart(crp, expand=FALSE)
    combo$addAttribute(crp, "pixbuf", 0)
    crt = gtkCellRendererText()
    crt['xalign'] <- 1
    combo$packStart(crt)
    combo$addAttribute(crt, "text", 1)
    combo$setActive(init.active)
    return(combo)
  } 
  
  getValueFromFigureCombo <- function(cb){
    store <- cb$GetModel()
    iter <- cb$GetActiveIter()
    if(iter$retval){
      value = store$GetValue(iter$iter,1)$value
      if(value %in% c("bubble","scaled","@","+","%","#","*","o","O")){
        data <- value
      }else{
        data <- as.numeric(value) 
      }
    }
    return(data)
  }
  
  ########## gtk functions
  .save_cb <- function(widget, window) {
    dialog <- gtkFileChooserDialog("Enter a name for the file", window,
                                   "save", "gtk-cancel", GtkResponseType["cancel"], "gtk-save",
                                   GtkResponseType["accept"])
    if (dialog$run() == GtkResponseType["accept"]){
      filename = dialog$getFilename()
      a = mrv$est
      if(!is.null(a)){
        save(a, file=paste(filename,".Rdata",sep=""))
      }else{
        errordialog <- gtkMessageDialog(window,"destroy-with-parent","warning","ok",
                                        "Please save after loading data and running model!")
        if (errordialog$run() == GtkResponseType["ok"]){}
        errordialog$destroy()
      }
    }
    dialog$destroy()
  }
  
  ##################################################################################
  #############           START
  ##################################################################################
  mrv$main_window <- gtkWindow("toplevel",show = FALSE)
  mrv$main_window["title"] <- "meta4diag: A GUI for diagnostic meta-analysis with INLA"
  mrv$main_window["border-width"] <- 0
  computer_screen_size <- gdkScreenGetDefault()
  width <- 0.6*computer_screen_size$getWidth()
  height <- 0.9*computer_screen_size$getHeight()
  mrv$main_window$setDefaultSize(width, height)
  mrv$main_window["window-position"]="center"
  
  mrv$datafile = NULL
  mrv$prec1.prior = NULL
  mrv$prec2.prior = NULL
  mrv$rho.prior = NULL
  ### define actions
  actions <- list(
    file = list("File", NULL, "_File", NULL, NULL, NULL),
    new = list("New", "gtk-new", "_New", "<control>N", 
               "New document", .someAction),
    sub = list("Submenu", NULL, "S_ub", NULL, NULL, NULL),
    open = list("Open", "gtk-open", "_Open", "<ctrl>0", 
                "Select a TXT or Rdata or CSV file to load data", .open_cb),
    save = list("Save", "gtk-save", "_Save", "<ctrl>S", 
                "Save document", .save_cb),
    quit = list("Quit", "gtk-quit", "_Quit", "<ctrl>Q", 
                "Quit the application", .quit_cb),
    edit = list("Edit", NULL, "_Edit", NULL, NULL, NULL),
    undo = list("Undo", "gtk-undo", "_Undo", "<ctrl>Z", 
                "Undo change", .someAction),
    redo = list("Redo", "gtk-redo", "_Redo", "<ctrl>U", 
                "Redo change", .someAction),
    exe = list("Execute", "gtk-execute", "_Execute", "<ctrl>R",
               "Build model with R-INLA",.executeFile),
    stp = list("Stop", "gtk-stop", "_Stop", "<ctrl>C",
               "Stop running",.stop_cb),
    about = list("About", NULL, "_About", NULL, "About meta4diag", .aboutgui),
    help = list("Help",NULL,"_Help",NULL,NULL,NULL),
    tips = list("Tips",NULL,"_UseTooltips",NULL,NULL,NULL)
  )
  
  mrv$action_group = gtkActionGroup("meta4diagActions")
  mrv$action_group$addActions(actions, mrv$main_window)
  mrv$ui_manager <- gtkUIManager()
  mrv$ui_manager$insertActionGroup(mrv$action_group, 0)
  mrv$id <- mrv$ui_manager$addUiFromFile(system.file("etc/ui.xml", package="meta4diag"))  
  mrv$menubar <- mrv$ui_manager$getWidget("/menubar")
  mrv$toolbar <- mrv$ui_manager$getWidget("/toolbar")
  mrv$main_window$addAccelGroup(mrv$ui_manager$getAccelGroup())
  ## 
  mrv$statusbar <- gtkStatusbar()
  mrv$info <- mrv$statusbar$getContextId("info")
  mrv$statusbar$push(mrv$info, "Ready")
  
  ###########################################################
  #################       start som hiden
  ###########################################################
  hide_parent_box <- gtkVBox()
  hide_parent <- gtkEntry()
  hide_parent["text"]="8"
  hide_parent_box$packStart(hide_parent)
  hide_parent_box$hide()
  ###########################################################
  ###########################################################
  #################       start construct side bar
  ###########################################################
  ###########################################################
  mrv$sidebar <- gtkVBox(homogeneous=FALSE)
  
  mrv$sidenote <- gtkNotebook()
  mrv$sidenote["homogeneous"] = FALSE
  mrv$sidenote["show-border"] = FALSE
  mrv$sidenote["tab-border"] = 10
  mrv$sidenote$setTabPos("left")
  mrv$sidenote$setScrollable(TRUE)
  
  sidebarScrolledWindow <- function(widget){
    window <- gtkScrolledWindow()
    window["shadow-type"] = "none"
    window$addWithViewport(widget)
    window$setPolicy("never", "automatic")
    return(window)
  }
  
  sidebarLabel <- function(string, xalign){
    label <- gtkLabel(string)
    label["xalign"] <- xalign
    return(label)
  }
  
  sidebarEntry <- function(text, wchars, xalign){
    entry <- gtkEntry()
    entry["text"] <- text
    entry["width-chars"] <- wchars
    entry["xalign"] <- xalign
    return(entry)
  }
  
  sidebarTable <- function(rowSpacing, colSpacing, widgetList){
    lenw <- length(widgetList)
    table <- gtkTable(rows=ceiling(0.5*lenw),columns=2,homogeneous=FALSE)
    table$setColSpacings(colSpacing)
    table$setRowSpacings(rowSpacing)
    for(i in 1:lenw){
      rowPos <- ceiling(0.5*i)
      colPos <- i %% 2
      if(colPos==1){
        table$attach(widgetList[[i]],left.attach=0,1, top.attach=(rowPos-1),rowPos,
                     xoptions=c("expand","fill"),yoptions="")
      }
      if(colPos==0){
        table$attach(widgetList[[i]],left.attach=1,2, top.attach=(rowPos-1),rowPos,
                     xoptions="",yoptions="")
      }
    }
    return(table)
  }
  
  sidebarPlotButton <- function(){
    bt <- gtkButton("Plot")
    bt['image'] <- gtkImage(filename=mrv$search_icon,size="button")
    return(bt)
  }
  
  sidebarLoadPriorButton <- function(){
    bt <- gtkButton("Load Prior")
    bt['image'] <- gtkImage(filename=mrv$search_icon,size="button")
    return(bt)
  }
  
  sidebarAcceptButton <- function(){
    bt <- gtkButton("Accept")
    bt['image'] <- gtkImage(stock.id="gtk-apply",size="button")
    return(bt)
  }
  ##################################################################
  ##################################################################
  ######
  ######                   sidebar priors
  ######
  ##################################################################
  ##################################################################
  mrv$sideprior <- gtkVBox(homogeneous=FALSE)
  mrv$sideprior_title <- gtkImage(filename=mrv$sidebar_icon[[1]])    #### maybe change move file to folder
  mrv$sideprior$packStart(mrv$sideprior_title,expand=FALSE,fill=FALSE)
  
  mrv$sideprior_window <- sidebarScrolledWindow(mrv$sideprior)
  ##################################################################
  ###                 component 1 - first variance
  ##################################################################
  mrv$sidecont1 <- gtkFrame("First Variance")
  mrv$sidecont1["border-width"]=10
  
  mrv$sidecont1_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_main["border-width"]=10
  mrv$sidecont1_show <- gtkHBox()
  mrv$sidecont1_combo <- gtkComboBoxNewText()
  sapply(c("Inv-gamma(var)","PC(std)","Half-Cauchy(std)","Truncated-normal(std)","Uniform(std)","Table(var)","Inv-Wishart(cov.mat)"), mrv$sidecont1_combo$appendText)
  mrv$sidecont1_show$packStart(gtkLabel("Prior:  "),expand=FALSE,fill=FALSE)
  mrv$sidecont1_show$packStart(mrv$sidecont1_combo,expand=TRUE,fill=TRUE)
  mrv$sidecont1_combo$setActive(-1)
  mrv$sidecont1_main$packStart(mrv$sidecont1_show,expand=FALSE,fill=FALSE)
  
  # Inv-gamma
  mrv$sidecont1_a_label <- sidebarLabel("shape: ",1)
  mrv$sidecont1_b_label <- sidebarLabel("rate: ",1)
  mrv$sidecont1_a_entry <- sidebarEntry("0.25", 10, 0)
  mrv$sidecont1_b_entry <- sidebarEntry("0.025", 10, 0)
  sidecont1_invgamma_list <- list(mrv$sidecont1_a_label, mrv$sidecont1_a_entry, mrv$sidecont1_b_label, mrv$sidecont1_b_entry)
  # PC
  mrv$sidecont1_u_label <- sidebarLabel("u: ", 1)
  mrv$sidecont1_alpha_label <- sidebarLabel("a: ", 1)
  mrv$sidecont1_u_entry <- sidebarEntry("3", 10, 0)
  mrv$sidecont1_alpha_entry <- sidebarEntry("0.05", 10, 0)
  sidecont1_pc_list <- list(mrv$sidecont1_u_label, mrv$sidecont1_u_entry, mrv$sidecont1_alpha_label, mrv$sidecont1_alpha_entry)
  # half-cauchy
  mrv$sidecont1_gamma_label <- sidebarLabel("gamma: ", 1)
  mrv$sidecont1_gamma_entry <- sidebarEntry("3", 10, 0)
  sidecont1_hcauchy_list <- list(mrv$sidecont1_gamma_label, mrv$sidecont1_gamma_entry)
  # truncated-normal
  mrv$sidecont1_m_label <- sidebarLabel("mean: ", 1)
  mrv$sidecont1_m_entry <- sidebarEntry("3", 10, 0)
  mrv$sidecont1_v_label <- sidebarLabel("variance: ", 1)
  mrv$sidecont1_v_entry <- sidebarEntry("1", 10, 0)
  sidecont1_tnorm_list <- list(mrv$sidecont1_m_label, mrv$sidecont1_m_entry, mrv$sidecont1_v_label, mrv$sidecont1_v_entry)
  
  mrv$sidecont1_input1 <- sidebarTable(3, 3, sidecont1_invgamma_list)
  mrv$sidecont1_input2 <- sidebarTable(3, 3, sidecont1_pc_list)
  mrv$sidecont1_input3 <- sidebarTable(3, 3, sidecont1_hcauchy_list)
  mrv$sidecont1_input4 <- sidebarTable(3, 3, sidecont1_tnorm_list)
  
  # Inv-gamma
  mrv$sidecont1_input1_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont1_input1_figurebt <- sidebarPlotButton()    
  mrv$sidecont1_input1_savebt <- sidebarAcceptButton()
  mrv$sidecont1_input1_buttonbox$packStart(mrv$sidecont1_input1_figurebt,expand=T,fill=T)
  mrv$sidecont1_input1_buttonbox$packStart(mrv$sidecont1_input1_savebt,expand=T,fill=T)
  # PC
  mrv$sidecont1_input2_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont1_input2_figurebt <- sidebarPlotButton()
  mrv$sidecont1_input2_savebt <- sidebarAcceptButton()
  mrv$sidecont1_input2_buttonbox$packStart(mrv$sidecont1_input2_figurebt,expand=T,fill=T)
  mrv$sidecont1_input2_buttonbox$packStart(mrv$sidecont1_input2_savebt,expand=T,fill=T)
  # half-cauchy
  mrv$sidecont1_input3_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont1_input3_figurebt <- sidebarPlotButton()
  mrv$sidecont1_input3_savebt <- sidebarAcceptButton()
  mrv$sidecont1_input3_buttonbox$packStart(mrv$sidecont1_input3_figurebt,expand=T,fill=T)
  mrv$sidecont1_input3_buttonbox$packStart(mrv$sidecont1_input3_savebt,expand=T,fill=T)
  # t-normal
  mrv$sidecont1_input4_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont1_input4_figurebt <- sidebarPlotButton()
  mrv$sidecont1_input4_savebt <- sidebarAcceptButton()
  mrv$sidecont1_input4_buttonbox$packStart(mrv$sidecont1_input4_figurebt,expand=T,fill=T)
  mrv$sidecont1_input4_buttonbox$packStart(mrv$sidecont1_input4_savebt,expand=T,fill=T)
  # Unif
  mrv$sidecont1_input5_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont1_input5_figurebt <- sidebarPlotButton()
  mrv$sidecont1_input5_savebt <- sidebarAcceptButton()
  mrv$sidecont1_input5_buttonbox$packStart(mrv$sidecont1_input5_figurebt,expand=T,fill=T)
  mrv$sidecont1_input5_buttonbox$packStart(mrv$sidecont1_input5_savebt,expand=T,fill=T)
  # table
  mrv$sidecont1_input6_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont1_input6_figurebt <- sidebarLoadPriorButton()
  mrv$sidecont1_input6_savebt <- sidebarAcceptButton()
  mrv$sidecont1_input6_buttonbox$packStart(mrv$sidecont1_input6_figurebt,expand=T,fill=T)
  mrv$sidecont1_input6_buttonbox$packStart(mrv$sidecont1_input6_savebt,expand=T,fill=T)
  
  # hide boxes
  mrv$sidecont1_hide11 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_hide11$hide()
  mrv$sidecont1_hide12 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_hide12$hide()
  mrv$sidecont1_hide21 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_hide21$hide()
  mrv$sidecont1_hide22 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_hide22$hide()
  mrv$sidecont1_hide31 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_hide31$hide()
  mrv$sidecont1_hide32 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_hide32$hide()
  mrv$sidecont1_hide41 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_hide41$hide()
  mrv$sidecont1_hide42 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_hide42$hide()
  mrv$sidecont1_hide5 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_hide5$hide()
  mrv$sidecont1_hide6 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont1_hide6$hide()
  
  mrv$sidecont1_hide11$packStart(mrv$sidecont1_input1,expand=FALSE,fill=FALSE)
  mrv$sidecont1_hide12$packStart(mrv$sidecont1_input1_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont1_hide21$packStart(mrv$sidecont1_input2,expand=FALSE,fill=FALSE)
  mrv$sidecont1_hide22$packStart(mrv$sidecont1_input2_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont1_hide31$packStart(mrv$sidecont1_input3,expand=FALSE,fill=FALSE)
  mrv$sidecont1_hide32$packStart(mrv$sidecont1_input3_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont1_hide41$packStart(mrv$sidecont1_input4,expand=FALSE,fill=FALSE)
  mrv$sidecont1_hide42$packStart(mrv$sidecont1_input4_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont1_hide5$packStart(mrv$sidecont1_input5_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont1_hide6$packStart(mrv$sidecont1_input6_buttonbox,expand=FALSE,fill=FALSE)
  
  gSignalConnect(mrv$sidecont1_a_entry, "changed", function(entry){
    mrv$sidecont1_hide12$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont1_b_entry, "changed", function(entry){
    mrv$sidecont1_hide12$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont1_u_entry, "changed", function(entry){
    mrv$sidecont1_hide22$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont1_alpha_entry, "changed", function(entry){
    mrv$sidecont1_hide22$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont1_gamma_entry, "changed", function(entry){
    mrv$sidecont1_hide32$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont1_m_entry, "changed", function(entry){
    mrv$sidecont1_hide42$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont1_v_entry, "changed", function(entry){
    mrv$sidecont1_hide42$show()
    .checkNumEntry(entry)
  })
  
  # Inv-gamma
  gSignalConnect(mrv$sidecont1_input1_figurebt, "clicked", function(button) {
    shape = as.numeric(gtkEntryGetText(mrv$sidecont1_a_entry))
    rate = as.numeric(gtkEntryGetText(mrv$sidecont1_b_entry))
    xlim = 8
    .manipulate(parent=list(mrv$sidecont1_a_entry,mrv$sidecont1_b_entry),
                .priorInvgamma(a = shape, b = rate, xmax = xlim),
                shape = .slider(parent=mrv$sidecont1_a_entry,min=0.000001,max=4,step=0.001),
                rate = .slider(parent=mrv$sidecont1_b_entry,min=0.000001,max=4,step=0.001),
                xlim = .slider(parent=hide_parent,min=1,max=20,step=1))
  })
  
  gSignalConnect(mrv$sidecont1_input1_savebt, "clicked", function(button) {
    tau1.a = as.numeric(gtkEntryGetText(mrv$sidecont1_a_entry))
    tau1.b = as.numeric(gtkEntryGetText(mrv$sidecont1_b_entry))
    mrv$var.prior = "Invgamma"
    mrv$var.par = c(tau1.a,tau1.b)
    mrv$sidecont1_hide11$show()
    mrv$sidecont1_hide12$hide()
  })

  # PC
  gSignalConnect(mrv$sidecont1_input2_figurebt, "clicked", function(button,...) {
    tau1density = gtkComboBoxGetActiveText(mrv$sidecont1_combo)
    if(tau1density=="PC(std)"){
      u = as.numeric(gtkEntryGetText(mrv$sidecont1_u_entry))
      a = as.numeric(gtkEntryGetText(mrv$sidecont1_alpha_entry))
      xlim = 8
      .manipulate(parent=list(mrv$sidecont1_u_entry,mrv$sidecont1_alpha_entry),
                  .priorSigmaPC(u = u, alpha = a, xmax = xlim),
                  u = .slider(parent=mrv$sidecont1_u_entry,min=0.000001,max=10,step=0.01),
                  a = .slider(parent=mrv$sidecont1_alpha_entry,min=0.000001,max=0.9999999,step=0.001),
                  xlim = .slider(parent=hide_parent,min=1,max=20,step=1))
    }
  })
  
  gSignalConnect(mrv$sidecont1_input2_savebt, "clicked", function(button) {
    tau1.u = as.numeric(gtkEntryGetText(mrv$sidecont1_u_entry))
    tau1.alpha = as.numeric(gtkEntryGetText(mrv$sidecont1_alpha_entry))
    mrv$var.prior = "PC"
    mrv$var.par = c(tau1.u,tau1.alpha)
    mrv$sidecont1_hide21$show()
    mrv$sidecont1_hide22$hide()
  })
  
  # Half-cauchy
  gSignalConnect(mrv$sidecont1_input3_figurebt, "clicked", function(button,...) {
    tau1density = gtkComboBoxGetActiveText(mrv$sidecont1_combo)
    if(tau1density=="Half-Cauchy(std)"){
      gamma = as.numeric(gtkEntryGetText(mrv$sidecont1_gamma_entry))
      xlim = 8
      .manipulate(parent=list(mrv$sidecont1_gamma_entry),
                  .priorHalfCauchy(gamma=gamma, xmax = xlim),
                  gamma = .slider(parent=mrv$sidecont1_gamma_entry,min=0.000001,max=10,step=0.01),
                  xlim = .slider(parent=hide_parent,min=1,max=20,step=1))
    }
  })
  
  gSignalConnect(mrv$sidecont1_input3_savebt, "clicked", function(button) {
    tau1.gamma = as.numeric(gtkEntryGetText(mrv$sidecont1_gamma_entry))
    mrv$var.prior = "HCauchy"
    mrv$var.par = c(tau1.gamma)
    mrv$sidecont1_hide31$show()
    mrv$sidecont1_hide32$hide()
  })
  
  # Truncated-normal
  gSignalConnect(mrv$sidecont1_input4_figurebt, "clicked", function(button,...) {
    tau1density = gtkComboBoxGetActiveText(mrv$sidecont1_combo)
    if(tau1density=="Truncated-normal(std)"){
      mean = as.numeric(gtkEntryGetText(mrv$sidecont1_m_entry))
      variance = as.numeric(gtkEntryGetText(mrv$sidecont1_v_entry))
      xlim = 8
      .manipulate(parent=list(mrv$sidecont1_m_entry, mrv$sidecont1_v_entry),
                  .priorSigmaTnorm(m=mean, v=variance, xmax = xlim),
                  mean = .slider(parent=mrv$sidecont1_m_entry,min=0.000001,max=10,step=0.01),
                  variance = .slider(parent=mrv$sidecont1_v_entry,min=0.000001,max=10,step=0.01),
                  xlim = .slider(parent=hide_parent,min=1,max=20,step=1))
    }
  })
  
  gSignalConnect(mrv$sidecont1_input4_savebt, "clicked", function(button) {
    tau1.m = as.numeric(gtkEntryGetText(mrv$sidecont1_m_entry))
    tau1.v = as.numeric(gtkEntryGetText(mrv$sidecont1_v_entry))
    mrv$var.prior = "Tnormal"
    mrv$var.par = c(tau1.m, tau1.v)
    mrv$sidecont1_hide41$show()
    mrv$sidecont1_hide42$hide()
  })

  # Uniform
  gSignalConnect(mrv$sidecont1_input5_figurebt, "clicked", function(button,...) {
    tau1density = gtkComboBoxGetActiveText(mrv$sidecont1_combo)
    if(tau1density=="Uniform(std)"){
      xlim = 8
      .manipulate(parent=list(hide_parent),
                  .priorUniformSig(xmax = xlim),
                  xlim = .slider(parent=hide_parent,min=1,max=20,step=1))
    }
  })
  
  gSignalConnect(mrv$sidecont1_input5_savebt, "clicked", function(button) {
    mrv$var.prior = "Unif"
    mrv$sidecont1_hide5$hide()
  })
  
  # Table
  gSignalConnect(mrv$sidecont1_input6_figurebt, "clicked", function(button,...) {
    tau1density = gtkComboBoxGetActiveText(mrv$sidecont1_combo)
    if(tau1density=="Table(var)"){
      .open_prior_Sigma(window=mrv$main_window)
    }
  })
  
  gSignalConnect(mrv$sidecont1_input6_savebt, "clicked", function(button) {
    mrv$var.prior = "Table"
    mrv$var.par = mrv$priorfile
    mrv$sidecont1_hide6$hide()
  })

  
  mrv$sidecont1_main$packStart(mrv$sidecont1_hide11,expand=FALSE,fill=FALSE)
  mrv$sidecont1_main$packStart(mrv$sidecont1_hide12,expand=FALSE,fill=FALSE)
  mrv$sidecont1_main$packStart(mrv$sidecont1_hide21,expand=FALSE,fill=FALSE)
  mrv$sidecont1_main$packStart(mrv$sidecont1_hide22,expand=FALSE,fill=FALSE)
  mrv$sidecont1_main$packStart(mrv$sidecont1_hide31,expand=FALSE,fill=FALSE)
  mrv$sidecont1_main$packStart(mrv$sidecont1_hide32,expand=FALSE,fill=FALSE)
  mrv$sidecont1_main$packStart(mrv$sidecont1_hide41,expand=FALSE,fill=FALSE)
  mrv$sidecont1_main$packStart(mrv$sidecont1_hide42,expand=FALSE,fill=FALSE)
  mrv$sidecont1_main$packStart(mrv$sidecont1_hide5,expand=FALSE,fill=FALSE)
  mrv$sidecont1_main$packStart(mrv$sidecont1_hide6,expand=FALSE,fill=FALSE)
  mrv$sidecont1$add(mrv$sidecont1_main)
  
  
  ##################################################################
  ###                 component 2 - second variance
  ##################################################################
  mrv$sidecont2 <- gtkFrame("Second Variance")
  mrv$sidecont2["border-width"]=10
  
  mrv$sidecont2_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_main["border-width"]=10
  mrv$sidecont2_show <- gtkHBox()
  mrv$sidecont2_combo <- gtkComboBoxNewText()
  sapply(c("Inv-gamma(var)","PC(std)","Half-Cauchy(std)","Truncated-normal(std)","Uniform(std)","Table(var)","Inv-Wishart(cov.mat)"), mrv$sidecont2_combo$appendText)
  mrv$sidecont2_show$packStart(gtkLabel("Prior:  "),expand=FALSE,fill=FALSE)
  mrv$sidecont2_show$packStart(mrv$sidecont2_combo,expand=TRUE,fill=TRUE)
  mrv$sidecont2_combo$setActive(-1)
  mrv$sidecont2_main$packStart(mrv$sidecont2_show,expand=FALSE,fill=FALSE)
  # Inv-gamma
  mrv$sidecont2_a_label <- sidebarLabel("shape: ",1)
  mrv$sidecont2_b_label <- sidebarLabel("rate: ",1)
  mrv$sidecont2_a_entry <- sidebarEntry("0.25", 10, 0)
  mrv$sidecont2_b_entry <- sidebarEntry("0.025", 10, 0)
  sidecont2_invgamma_list <- list(mrv$sidecont2_a_label, mrv$sidecont2_a_entry, mrv$sidecont2_b_label, mrv$sidecont2_b_entry)
  # PC
  mrv$sidecont2_u_label <- sidebarLabel("u: ",1)
  mrv$sidecont2_alpha_label <- sidebarLabel("a: ",1)
  mrv$sidecont2_u_entry <- sidebarEntry("3", 10, 0)
  mrv$sidecont2_alpha_entry <- sidebarEntry("0.05", 10, 0)
  sidecont2_pc_list <- list(mrv$sidecont2_u_label, mrv$sidecont2_u_entry, mrv$sidecont2_alpha_label, mrv$sidecont2_alpha_entry)
  # half-cauchy
  mrv$sidecont2_gamma_label <- sidebarLabel("gamma: ", 1)
  mrv$sidecont2_gamma_entry <- sidebarEntry("3", 10, 0)
  sidecont2_hcauchy_list <- list(mrv$sidecont2_gamma_label, mrv$sidecont2_gamma_entry)
  # truncated-normal
  mrv$sidecont2_m_label <- sidebarLabel("mean: ", 1)
  mrv$sidecont2_m_entry <- sidebarEntry("3", 10, 0)
  mrv$sidecont2_v_label <- sidebarLabel("variance: ", 1)
  mrv$sidecont2_v_entry <- sidebarEntry("1", 10, 0)
  sidecont2_tnorm_list <- list(mrv$sidecont2_m_label, mrv$sidecont2_m_entry, mrv$sidecont2_v_label, mrv$sidecont2_v_entry)
  
  mrv$sidecont2_input1 <- sidebarTable(3, 3, sidecont2_invgamma_list)
  mrv$sidecont2_input2 <- sidebarTable(3, 3, sidecont2_pc_list)
  mrv$sidecont2_input3 <- sidebarTable(3, 3, sidecont2_hcauchy_list)
  mrv$sidecont2_input4 <- sidebarTable(3, 3, sidecont2_tnorm_list)
  
  # Inv-gamma
  mrv$sidecont2_input1_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont2_input1_figurebt <- sidebarPlotButton()
  mrv$sidecont2_input1_savebt <- sidebarAcceptButton()
  mrv$sidecont2_input1_buttonbox$packStart(mrv$sidecont2_input1_figurebt,expand=T,fill=T)
  mrv$sidecont2_input1_buttonbox$packStart(mrv$sidecont2_input1_savebt,expand=T,fill=T)
  # PC
  mrv$sidecont2_input2_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont2_input2_figurebt <- sidebarPlotButton()
  mrv$sidecont2_input2_savebt <- sidebarAcceptButton()
  mrv$sidecont2_input2_buttonbox$packStart(mrv$sidecont2_input2_figurebt,expand=T,fill=T)
  mrv$sidecont2_input2_buttonbox$packStart(mrv$sidecont2_input2_savebt,expand=T,fill=T)
  # half-cauchy
  mrv$sidecont2_input3_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont2_input3_figurebt <- sidebarPlotButton()
  mrv$sidecont2_input3_savebt <- sidebarAcceptButton()
  mrv$sidecont2_input3_buttonbox$packStart(mrv$sidecont2_input3_figurebt,expand=T,fill=T)
  mrv$sidecont2_input3_buttonbox$packStart(mrv$sidecont2_input3_savebt,expand=T,fill=T)
  # t-normal
  mrv$sidecont2_input4_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont2_input4_figurebt <- sidebarPlotButton()
  mrv$sidecont2_input4_savebt <- sidebarAcceptButton()
  mrv$sidecont2_input4_buttonbox$packStart(mrv$sidecont2_input4_figurebt,expand=T,fill=T)
  mrv$sidecont2_input4_buttonbox$packStart(mrv$sidecont2_input4_savebt,expand=T,fill=T)
  # Unif
  mrv$sidecont2_input5_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont2_input5_figurebt <- sidebarPlotButton()
  mrv$sidecont2_input5_savebt <- sidebarAcceptButton()
  mrv$sidecont2_input5_buttonbox$packStart(mrv$sidecont2_input5_figurebt,expand=T,fill=T)
  mrv$sidecont2_input5_buttonbox$packStart(mrv$sidecont2_input5_savebt,expand=T,fill=T)
  # table
  mrv$sidecont2_input6_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont2_input6_figurebt <- sidebarLoadPriorButton()
  mrv$sidecont2_input6_savebt <- sidebarAcceptButton()
  mrv$sidecont2_input6_buttonbox$packStart(mrv$sidecont2_input6_figurebt,expand=T,fill=T)
  mrv$sidecont2_input6_buttonbox$packStart(mrv$sidecont2_input6_savebt,expand=T,fill=T)
  
  # hide boxed
  mrv$sidecont2_hide11 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_hide11$hide()
  mrv$sidecont2_hide21 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_hide21$hide()
  mrv$sidecont2_hide12 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_hide12$hide()
  mrv$sidecont2_hide22 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_hide22$hide()
  mrv$sidecont2_hide31 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_hide31$hide()
  mrv$sidecont2_hide32 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_hide32$hide()
  mrv$sidecont2_hide41 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_hide41$hide()
  mrv$sidecont2_hide42 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_hide42$hide()
  mrv$sidecont2_hide5 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_hide5$hide()
  mrv$sidecont2_hide6 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont2_hide6$hide()
  
  mrv$sidecont2_hide11$packStart(mrv$sidecont2_input1,expand=FALSE,fill=FALSE)
  mrv$sidecont2_hide12$packStart(mrv$sidecont2_input1_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont2_hide21$packStart(mrv$sidecont2_input2,expand=FALSE,fill=FALSE)
  mrv$sidecont2_hide22$packStart(mrv$sidecont2_input2_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont2_hide31$packStart(mrv$sidecont2_input3,expand=FALSE,fill=FALSE)
  mrv$sidecont2_hide32$packStart(mrv$sidecont2_input3_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont2_hide41$packStart(mrv$sidecont2_input4,expand=FALSE,fill=FALSE)
  mrv$sidecont2_hide42$packStart(mrv$sidecont2_input4_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont2_hide5$packStart(mrv$sidecont2_input5_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont2_hide6$packStart(mrv$sidecont2_input6_buttonbox,expand=FALSE,fill=FALSE)
  
  
  
  gSignalConnect(mrv$sidecont2_a_entry, "changed", function(entry){
    mrv$sidecont2_hide12$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont2_b_entry, "changed", function(entry){
    mrv$sidecont2_hide12$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont2_u_entry, "changed", function(entry){
    mrv$sidecont2_hide22$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont2_alpha_entry, "changed", function(entry){
    mrv$sidecont2_hide22$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont2_gamma_entry, "changed", function(entry){
    mrv$sidecont2_hide32$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont2_m_entry, "changed", function(entry){
    mrv$sidecont2_hide42$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont2_v_entry, "changed", function(entry){
    mrv$sidecont2_hide42$show()
    .checkNumEntry(entry)
  })
  
  # Inv-gamma
  gSignalConnect(mrv$sidecont2_input1_figurebt, "clicked", function(button) {
    shape = as.numeric(gtkEntryGetText(mrv$sidecont2_a_entry))
    rate = as.numeric(gtkEntryGetText(mrv$sidecont2_b_entry))
    xlim = 8
    .manipulate(parent=list(mrv$sidecont2_a_entry,mrv$sidecont2_b_entry),
                .priorInvgamma(a = shape,b = rate,xmax = xlim),
                shape = .slider(parent=mrv$sidecont2_a_entry,min=0.000001,max=4,step=0.001),
                rate = .slider(parent=mrv$sidecont2_b_entry,min=0.000001,max=4,step=0.001),
                xlim = .slider(parent=hide_parent,label="xlim",min=1,max=20,step=1))
  })
  
  gSignalConnect(mrv$sidecont2_input1_savebt, "clicked", function(button) {
    tau2.a = as.numeric(gtkEntryGetText(mrv$sidecont2_a_entry))
    tau2.b = as.numeric(gtkEntryGetText(mrv$sidecont2_b_entry))
    mrv$var2.prior = "Invgamma"
    mrv$var2.par = c(tau2.a,tau2.b)
    mrv$sidecont2_hide11$show()
    mrv$sidecont2_hide12$hide()
  })
  
  # PC
  gSignalConnect(mrv$sidecont2_input2_figurebt, "clicked", function(button,...) {
    tau2density = gtkComboBoxGetActiveText(mrv$sidecont2_combo)
    if(tau2density=="PC(std)"){
      u = as.numeric(gtkEntryGetText(mrv$sidecont2_u_entry))
      a = as.numeric(gtkEntryGetText(mrv$sidecont2_alpha_entry))
      xlim = 8
      .manipulate(parent=list(mrv$sidecont2_u_entry,mrv$sidecont2_alpha_entry),
                  .priorSigmaPC(u, a, xlim),
                  u = .slider(parent=mrv$sidecont2_u_entry,min=0.000001,max=10,step=0.01),
                  a = .slider(parent=mrv$sidecont2_alpha_entry,min=0.000001,max=0.9999999,step=0.001),
                  xlim = .slider(parent=hide_parent,min=1,label="xlim",max=20,step=1))
    }
  })
  
  gSignalConnect(mrv$sidecont2_input2_savebt, "clicked", function(button) {
    tau2.u = as.numeric(gtkEntryGetText(mrv$sidecont2_u_entry))
    tau2.alpha = as.numeric(gtkEntryGetText(mrv$sidecont2_alpha_entry))
    mrv$var2.prior = "PC"
    mrv$var2.par = c(tau2.u,tau2.alpha)
    mrv$sidecont2_hide21$show()
    mrv$sidecont2_hide22$hide()
  })
  
  # Half-cauchy
  gSignalConnect(mrv$sidecont2_input3_figurebt, "clicked", function(button,...) {
    tau2density = gtkComboBoxGetActiveText(mrv$sidecont2_combo)
    if(tau2density=="Half-Cauchy(std)"){
      gamma = as.numeric(gtkEntryGetText(mrv$sidecont2_gamma_entry))
      xlim = 8
      .manipulate(parent=list(mrv$sidecont2_gamma_entry),
                  .priorHalfCauchy(gamma=gamma, xmax = xlim),
                  gamma = .slider(parent=mrv$sidecont2_gamma_entry,min=0.000001,max=10,step=0.01),
                  xlim = .slider(parent=hide_parent,min=1,max=20,step=1))
    }
  })
  
  gSignalConnect(mrv$sidecont2_input3_savebt, "clicked", function(button) {
    tau2.gamma = as.numeric(gtkEntryGetText(mrv$sidecont2_gamma_entry))
    mrv$var2.prior = "HCauchy"
    mrv$var2.par = c(tau2.gamma)
    mrv$sidecont2_hide31$show()
    mrv$sidecont2_hide32$hide()
  })
  
  # Truncated-normal
  gSignalConnect(mrv$sidecont2_input4_figurebt, "clicked", function(button,...) {
    tau2density = gtkComboBoxGetActiveText(mrv$sidecont2_combo)
    if(tau2density=="Truncated-normal(std)"){
      mean = as.numeric(gtkEntryGetText(mrv$sidecont2_m_entry))
      variance = as.numeric(gtkEntryGetText(mrv$sidecont2_v_entry))
      xlim = 8
      .manipulate(parent=list(mrv$sidecont2_m_entry, mrv$sidecont2_v_entry),
                  .priorSigmaTnorm(m=mean, v=variance, xmax = xlim),
                  mean = .slider(parent=mrv$sidecont2_m_entry,min=0.000001,max=10,step=0.01),
                  variance = .slider(parent=mrv$sidecont2_v_entry,min=0.000001,max=10,step=0.01),
                  xlim = .slider(parent=hide_parent,min=1,max=20,step=1))
    }
  })
  
  gSignalConnect(mrv$sidecont2_input4_savebt, "clicked", function(button) {
    tau2.m = as.numeric(gtkEntryGetText(mrv$sidecont2_m_entry))
    tau2.v = as.numeric(gtkEntryGetText(mrv$sidecont2_v_entry))
    mrv$var2.prior = "Tnormal"
    mrv$var2.par = c(tau2.m, tau2.v)
    mrv$sidecont2_hide41$show()
    mrv$sidecont2_hide42$hide()
  })
  
  # Uniform
  gSignalConnect(mrv$sidecont2_input5_figurebt, "clicked", function(button,...) {
    tau2density = gtkComboBoxGetActiveText(mrv$sidecont2_combo)
    if(tau2density=="Uniform(std)"){
      xlim = 8
      .manipulate(parent=list(hide_parent),
                  .priorUniformSig(xmax = xlim),
                  xlim = .slider(parent=hide_parent,min=1,max=20,step=1))
    }
  })
  
  gSignalConnect(mrv$sidecont2_input5_savebt, "clicked", function(button) {
    mrv$var2.prior = "Unif"
    mrv$sidecont2_hide5$hide()
  })
  
  # Table
  gSignalConnect(mrv$sidecont2_input6_figurebt, "clicked", function(button,...) {
    tau2density = gtkComboBoxGetActiveText(mrv$sidecont2_combo)
    if(tau2density=="Table(var)"){
      .open_prior_Sigma(window=mrv$main_window)
    }
  })
  
  gSignalConnect(mrv$sidecont2_input6_savebt, "clicked", function(button) {
    mrv$var2.prior = "Table"
    mrv$var2.par = mrv$priorfile
    mrv$sidecont2_hide6$hide()
  })
  

  mrv$sidecont2_main$packStart(mrv$sidecont2_hide11,expand=FALSE,fill=FALSE)
  mrv$sidecont2_main$packStart(mrv$sidecont2_hide12,expand=FALSE,fill=FALSE)
  mrv$sidecont2_main$packStart(mrv$sidecont2_hide21,expand=FALSE,fill=FALSE)
  mrv$sidecont2_main$packStart(mrv$sidecont2_hide22,expand=FALSE,fill=FALSE)
  mrv$sidecont2_main$packStart(mrv$sidecont2_hide31,expand=FALSE,fill=FALSE)
  mrv$sidecont2_main$packStart(mrv$sidecont2_hide32,expand=FALSE,fill=FALSE)
  mrv$sidecont2_main$packStart(mrv$sidecont2_hide41,expand=FALSE,fill=FALSE)
  mrv$sidecont2_main$packStart(mrv$sidecont2_hide42,expand=FALSE,fill=FALSE)
  mrv$sidecont2_main$packStart(mrv$sidecont2_hide5,expand=FALSE,fill=FALSE)
  mrv$sidecont2_main$packStart(mrv$sidecont2_hide6,expand=FALSE,fill=FALSE)
  mrv$sidecont2$add(mrv$sidecont2_main)
  
  ##################################################################
  ####                   component3 - correlation
  ##################################################################
  mrv$sidecont3 <- gtkFrame("Correlation")
  mrv$sidecont3["border-width"]=10
  
  mrv$sidecont3_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_main["border-width"]=10
  mrv$sidecont3_show <- gtkHBox()
  mrv$sidecont3_combo <- gtkComboBoxNewText()
  sapply(c("Gaussian(Fisher-z)","PC(cor)","Beta(transf.cor)","Table(cor)","Inv-Wishart(cov.mat)"), mrv$sidecont3_combo$appendText)
  mrv$sidecont3_show$packStart(gtkLabel("Prior:  "),expand=FALSE,fill=FALSE)
  mrv$sidecont3_show$packStart(mrv$sidecont3_combo,expand=TRUE,fill=TRUE)
  mrv$sidecont3_combo$setActive(-1)
  mrv$sidecont3_main$packStart(mrv$sidecont3_show,expand=FALSE,fill=FALSE)
  
  mrv$sidecont3_hide11 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_hide11$hide()
  mrv$sidecont3_hide12 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_hide12$hide()
  mrv$sidecont3_hide2 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_hide2$hide()
  mrv$sidecont3_hide31 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_hide31$hide()
  mrv$sidecont3_hide32 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_hide32$hide()
  mrv$sidecont3_hide4 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_hide4$hide()
  ####################
  # Normal
  ####################
  mrv$sidecont3_m_label <- sidebarLabel("mu: ",1)
  mrv$sidecont3_v_label <- sidebarLabel("variance: ",1)
  mrv$sidecont3_m_entry <- sidebarEntry("0", 10, 0)
  mrv$sidecont3_v_entry <- sidebarEntry("5", 10, 0)
  sidecont3_Gaussian_list <- list(mrv$sidecont3_m_label, mrv$sidecont3_m_entry, mrv$sidecont3_v_label, mrv$sidecont3_v_entry)
  
  mrv$sidecont3_input1 <- sidebarTable(3, 3, sidecont3_Gaussian_list)
  
  mrv$sidecont3_input1_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont3_input1_figurebt <- sidebarPlotButton()
  mrv$sidecont3_input1_savebt <- sidebarAcceptButton()
  mrv$sidecont3_input1_buttonbox$packStart(mrv$sidecont3_input1_figurebt,expand=T,fill=T)
  mrv$sidecont3_input1_buttonbox$packStart(mrv$sidecont3_input1_savebt,expand=T,fill=T)
  ####################
  # PC
  ####################
  mrv$sidecont3_input2 <- gtkTable(rows=4,columns=2,homogeneous=FALSE)
  mrv$sidecont3_input2$setColSpacings(5)
  mrv$sidecont3_input2$setRowSpacings(3)
  
  mrv$sidecont3_ref_label <- sidebarLabel("reference: ",1)
  mrv$sidecont3_per_label <- sidebarLabel("percentile: ",1)
  mrv$sidecont3_ref_entry <- sidebarEntry("0.2", 10, 0)
  mrv$sidecont3_per_entry <- sidebarEntry("0.5", 10, 0)

  
  mrv$sidecont3_radio_label <- gtkLabel("Strategy: ")
  mrv$sidecont3_radio_label["xalign"]=1
  mrv$sidecont3_radio_label["yalign"]=0
  mrv$sidecont3_radio <- gtkVBox()
  mrv$sidecont3_radiogp <- list()
  mrv$sidecont3_radiogp$le <- gtkRadioButton(label = "Left-extrem")
  mrv$sidecont3_radiogp$re <- gtkRadioButton(mrv$sidecont3_radiogp, label = "Right-extrem")
  mrv$sidecont3_radiogp$ds <- gtkRadioButton(mrv$sidecont3_radiogp, label = "Double-side")
  sapply(mrv$sidecont3_radiogp, mrv$sidecont3_radio$packStart)
  mrv$sidecont3_radio[[1]]$setActive(TRUE) 
  sapply(mrv$sidecont3_radiogp,'[',"active")
  mrv$sidecont3_radio_align <- gtkAlignment(xalign = 0)
  mrv$sidecont3_radio_align$add(mrv$sidecont3_radio)
  
  mrv$sidecont3_input2$attach(mrv$sidecont3_ref_label,left.attach=0,1, top.attach=0,1,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidecont3_input2$attach(mrv$sidecont3_ref_entry,left.attach=1,2, top.attach=0,1,
                              xoptions="",yoptions="")
  mrv$sidecont3_input2$attach(mrv$sidecont3_per_label,left.attach=0,1, top.attach=1,2,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidecont3_input2$attach(mrv$sidecont3_per_entry,left.attach=1,2, top.attach=1,2,
                              xoptions="",yoptions="")
  mrv$sidecont3_input2$attach(mrv$sidecont3_radio_label,left.attach=0,1, top.attach=2,3,
                              xoptions=c("expand","fill"),yoptions=c("expand","fill"))
  mrv$sidecont3_input2$attach(mrv$sidecont3_radio_align,left.attach=1,2, top.attach=2,3,
                              xoptions=c("expand","fill"),yoptions="")
  
  ######################## sub-hiden active by radio in correlation part
  mrv$sidecont3_input2_subhide11 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_input2_subhide11$show() # left-extrem
  mrv$sidecont3_input2_subhide12 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_input2_subhide12$show() 
  mrv$sidecont3_input2_subhide21 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_input2_subhide21$hide() # right-extrem
  mrv$sidecont3_input2_subhide22 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_input2_subhide22$hide() 
  mrv$sidecont3_input2_subhide31 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_input2_subhide31$hide() # double-side
  mrv$sidecont3_input2_subhide32 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont3_input2_subhide32$hide()
  
  ## here starts sub-hide 1: left-extrem
  mrv$sidecont3_sub1_U1_label <- sidebarLabel("u-min: ",1)
  mrv$sidecont3_sub1_A1_label <- sidebarLabel("a1: ",1)
  mrv$sidecont3_sub1_U1_entry <- sidebarEntry("-0.95", 10, 0)
  mrv$sidecont3_sub1_A1_entry <- sidebarEntry("0.05", 10, 0)
  sidecont3_sub1_list <- list(mrv$sidecont3_sub1_U1_label, mrv$sidecont3_sub1_U1_entry, mrv$sidecont3_sub1_A1_label, mrv$sidecont3_sub1_A1_entry)
  mrv$sidecont3_sub1_input <- sidebarTable(3, 3, sidecont3_sub1_list)
  
  mrv$sidecont3_sub1_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont3_sub1_figurebt <- sidebarPlotButton()
  mrv$sidecont3_sub1_savebt <- sidebarAcceptButton()
  mrv$sidecont3_sub1_buttonbox$packStart(mrv$sidecont3_sub1_figurebt,expand=T,fill=T)
  mrv$sidecont3_sub1_buttonbox$packStart(mrv$sidecont3_sub1_savebt,expand=T,fill=T)
  
  
  mrv$sidecont3_input2_subhide11$packStart(mrv$sidecont3_sub1_input,expand=FALSE,fill=FALSE)
  mrv$sidecont3_input2_subhide12$packStart(mrv$sidecont3_sub1_buttonbox,expand=FALSE,fill=FALSE)
  
  ## here starts sub-hide 2: right extreme
  mrv$sidecont3_sub2_U2_label <- sidebarLabel("u-max: ",1)
  mrv$sidecont3_sub2_A2_label <- sidebarLabel("a2: ",1)
  mrv$sidecont3_sub2_U2_entry <- sidebarEntry("0.95", 10, 0)
  mrv$sidecont3_sub2_A2_entry <- sidebarEntry("0.05", 10, 0)
  sidecont3_sub2_list <- list(mrv$sidecont3_sub2_U2_label, mrv$sidecont3_sub2_U2_entry, mrv$sidecont3_sub2_A2_label, mrv$sidecont3_sub2_A2_entry)
  mrv$sidecont3_sub2_input <- sidebarTable(3, 3, sidecont3_sub2_list)
  
  
  mrv$sidecont3_sub2_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont3_sub2_figurebt <- sidebarPlotButton()
  mrv$sidecont3_sub2_savebt <- sidebarAcceptButton()
  mrv$sidecont3_sub2_buttonbox$packStart(mrv$sidecont3_sub2_figurebt,expand=T,fill=T)
  mrv$sidecont3_sub2_buttonbox$packStart(mrv$sidecont3_sub2_savebt,expand=T,fill=T)
  
  
  mrv$sidecont3_input2_subhide21$packStart(mrv$sidecont3_sub2_input,expand=FALSE,fill=FALSE)
  mrv$sidecont3_input2_subhide22$packStart(mrv$sidecont3_sub2_buttonbox,expand=FALSE,fill=FALSE)
  
  
  # here starts sub-hide 3: double side  
  mrv$sidecont3_sub3_U1_label <- sidebarLabel("u-min: ",1)
  mrv$sidecont3_sub3_A1_label <- sidebarLabel("a1: ",1)
  mrv$sidecont3_sub3_U1_entry <- sidebarEntry("-0.95", 10, 0)
  mrv$sidecont3_sub3_A1_entry <- sidebarEntry("0.05", 10, 0)
  mrv$sidecont3_sub3_U2_label <- sidebarLabel("u-max: ",1)
  mrv$sidecont3_sub3_A2_label <- sidebarLabel("a2: ",1)
  mrv$sidecont3_sub3_U2_entry <- sidebarEntry("0.95", 10, 0)
  mrv$sidecont3_sub3_A2_entry <- sidebarEntry("0.05", 10, 0)
  sidecont3_sub3_list <- list(mrv$sidecont3_sub3_U1_label, mrv$sidecont3_sub3_U1_entry, mrv$sidecont3_sub3_A1_label, mrv$sidecont3_sub3_A1_entry,
                              mrv$sidecont3_sub3_U2_label, mrv$sidecont3_sub3_U2_entry, mrv$sidecont3_sub3_A2_label, mrv$sidecont3_sub3_A2_entry)
  mrv$sidecont3_sub3_input <- sidebarTable(3, 3, sidecont3_sub3_list)
  
  mrv$sidecont3_sub3_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont3_sub3_figurebt <- sidebarPlotButton()
  mrv$sidecont3_sub3_savebt <- sidebarAcceptButton()
  mrv$sidecont3_sub3_buttonbox$packStart(mrv$sidecont3_sub3_figurebt,expand=T,fill=T)
  mrv$sidecont3_sub3_buttonbox$packStart(mrv$sidecont3_sub3_savebt,expand=T,fill=T)
  
  
  mrv$sidecont3_input2_subhide31$packStart(mrv$sidecont3_sub3_input,expand=FALSE,fill=FALSE)
  mrv$sidecont3_input2_subhide32$packStart(mrv$sidecont3_sub3_buttonbox,expand=FALSE,fill=FALSE)
  
  ####################
  # Beta
  ####################
  mrv$sidecont3_a_label <- sidebarLabel("a: ",1)
  mrv$sidecont3_b_label <- sidebarLabel("b: ",1)
  mrv$sidecont3_a_entry <- sidebarEntry("1", 10, 0)
  mrv$sidecont3_b_entry <- sidebarEntry("5", 10, 0)
  sidecont3_Beta_list <- list(mrv$sidecont3_a_label, mrv$sidecont3_a_entry, mrv$sidecont3_b_label, mrv$sidecont3_b_entry)
  
  mrv$sidecont3_input3 <- sidebarTable(3, 3, sidecont3_Beta_list)
  
  mrv$sidecont3_input3_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont3_input3_figurebt <- sidebarPlotButton()
  mrv$sidecont3_input3_savebt <- sidebarAcceptButton()
  mrv$sidecont3_input3_buttonbox$packStart(mrv$sidecont3_input3_figurebt,expand=T,fill=T)
  mrv$sidecont3_input3_buttonbox$packStart(mrv$sidecont3_input3_savebt,expand=T,fill=T)
  
  ####################
  # Table
  ####################
  mrv$sidecont3_input4_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont3_input4_figurebt <- sidebarLoadPriorButton()
  mrv$sidecont3_input4_savebt <- sidebarAcceptButton()
  mrv$sidecont3_input4_buttonbox$packStart(mrv$sidecont3_input4_figurebt,expand=T,fill=T)
  mrv$sidecont3_input4_buttonbox$packStart(mrv$sidecont3_input4_savebt,expand=T,fill=T)
  
  ########### combine all the functional part in part 3
  mrv$sidecont3_hide11$packStart(mrv$sidecont3_input1,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide12$packStart(mrv$sidecont3_input1_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide2$packStart(mrv$sidecont3_input2,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide2$packStart(mrv$sidecont3_input2_subhide11,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide2$packStart(mrv$sidecont3_input2_subhide12,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide2$packStart(mrv$sidecont3_input2_subhide21,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide2$packStart(mrv$sidecont3_input2_subhide22,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide2$packStart(mrv$sidecont3_input2_subhide31,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide2$packStart(mrv$sidecont3_input2_subhide32,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide31$packStart(mrv$sidecont3_input3,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide32$packStart(mrv$sidecont3_input3_buttonbox,expand=FALSE,fill=FALSE)
  mrv$sidecont3_hide4$packStart(mrv$sidecont3_input4_buttonbox,expand=FALSE,fill=FALSE)
  
  mrv$sidecont3_main$packStart(mrv$sidecont3_hide11,expand=FALSE,fill=FALSE)
  mrv$sidecont3_main$packStart(mrv$sidecont3_hide12,expand=FALSE,fill=FALSE)
  mrv$sidecont3_main$packStart(mrv$sidecont3_hide2,expand=FALSE,fill=FALSE)
  mrv$sidecont3_main$packStart(mrv$sidecont3_hide31,expand=FALSE,fill=FALSE)
  mrv$sidecont3_main$packStart(mrv$sidecont3_hide32,expand=FALSE,fill=FALSE)
  mrv$sidecont3_main$packStart(mrv$sidecont3_hide4,expand=FALSE,fill=FALSE)
  mrv$sidecont3$add(mrv$sidecont3_main)
  ######### all signals for sidebar part 3
  # normal
  gSignalConnect(mrv$sidecont3_m_entry, "changed", function(entry){
    mrv$sidecont3_hide12$show()
    mrv$sidecont3_hide2$hide()
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_v_entry, "changed", function(entry){
    mrv$sidecont3_hide12$show()
    mrv$sidecont3_hide2$hide()
    .checkNumEntry(entry)
  })
  # PC
  gSignalConnect(mrv$sidecont3_ref_entry, "changed", function(entry){
    if(mrv$rho.strategy=="Double-side"){
      mrv$sidecont3_input2_subhide32$show()
    }else if(mrv$rho.strategy=="Left-extrem"){
      mrv$sidecont3_input2_subhide12$show()
    }else{
      mrv$sidecont3_input2_subhide22$show()
    }
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_per_entry, "changed", function(entry){
    if(mrv$rho.strategy=="Double-side"){
      mrv$sidecont3_input2_subhide32$show()
    }else if(mrv$rho.strategy=="Left-extrem"){
      mrv$sidecont3_input2_subhide12$show()
    }else{
      mrv$sidecont3_input2_subhide22$show()
    }
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_sub1_U1_entry, "changed", function(entry){
    mrv$sidecont3_input2_subhide12$show()
    mrv$sidecont3_input2_subhide22$hide()
    mrv$sidecont3_input2_subhide32$hide()
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_sub1_A1_entry, "changed", function(entry){
    mrv$sidecont3_input2_subhide12$show()
    mrv$sidecont3_input2_subhide22$hide()
    mrv$sidecont3_input2_subhide32$hide()
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_sub2_U2_entry, "changed", function(entry){
    mrv$sidecont3_input2_subhide22$show()
    mrv$sidecont3_input2_subhide12$hide()
    mrv$sidecont3_input2_subhide32$hide()
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_sub2_A2_entry, "changed", function(entry){
    mrv$sidecont3_input2_subhide22$show()
    mrv$sidecont3_input2_subhide12$hide()
    mrv$sidecont3_input2_subhide32$hide()
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_sub3_U1_entry, "changed", function(entry){
    mrv$sidecont3_input2_subhide32$show()
    mrv$sidecont3_input2_subhide22$hide()
    mrv$sidecont3_input2_subhide12$hide()
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_sub3_A1_entry, "changed", function(entry){
    mrv$sidecont3_input2_subhide32$show()
    mrv$sidecont3_input2_subhide22$hide()
    mrv$sidecont3_input2_subhide12$hide()
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_sub3_U2_entry, "changed", function(entry){
    mrv$sidecont3_input2_subhide32$show()
    mrv$sidecont3_input2_subhide22$hide()
    mrv$sidecont3_input2_subhide12$hide()
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_sub3_A2_entry, "changed", function(entry){
    mrv$sidecont3_input2_subhide32$show()
    mrv$sidecont3_input2_subhide22$hide()
    mrv$sidecont3_input2_subhide12$hide()
    .checkNumEntry(entry)
  })
  mrv$rho.strategy = "Left-extrem"
  sapply(mrv$sidecont3_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             mrv$rho.strategy = button$getLabel()
             if(mrv$rho.strategy=="Double-side"){
               mrv$sidecont3_input2_subhide11$hide()
               mrv$sidecont3_input2_subhide12$hide()
               mrv$sidecont3_input2_subhide21$hide()
               mrv$sidecont3_input2_subhide22$hide()
               mrv$sidecont3_input2_subhide31$show()
               mrv$sidecont3_input2_subhide32$show()
             } else if (mrv$rho.strategy=="Left-extrem"){
               mrv$sidecont3_input2_subhide21$hide()
               mrv$sidecont3_input2_subhide22$hide()
               mrv$sidecont3_input2_subhide31$hide()
               mrv$sidecont3_input2_subhide32$hide()
               mrv$sidecont3_input2_subhide11$show()
               mrv$sidecont3_input2_subhide12$show()
             } else {
               mrv$sidecont3_input2_subhide11$hide()
               mrv$sidecont3_input2_subhide12$hide()
               mrv$sidecont3_input2_subhide31$hide()
               mrv$sidecont3_input2_subhide32$hide()
               mrv$sidecont3_input2_subhide21$show()
               mrv$sidecont3_input2_subhide22$show()
             }
           } 
         })
  # Beta
  gSignalConnect(mrv$sidecont3_a_entry, "changed", function(entry){
    mrv$sidecont3_hide32$show()
    mrv$sidecont3_hide2$hide()
    .checkNumEntry(entry)
  })
  gSignalConnect(mrv$sidecont3_b_entry, "changed", function(entry){
    mrv$sidecont3_hide32$show()
    mrv$sidecont3_hide2$hide()
    .checkNumEntry(entry)
  })
  
  # normal
  gSignalConnect(mrv$sidecont3_input1_figurebt, "clicked", function(button){
    mean = as.numeric(gtkEntryGetText(mrv$sidecont3_m_entry))
    variance = as.numeric(gtkEntryGetText(mrv$sidecont3_v_entry))
    .manipulate(parent=list(mrv$sidecont3_m_entry,mrv$sidecont3_v_entry),
                .priorRhoNormalPlot(mean,variance),
                mean = .slider(parent=mrv$sidecont3_m_entry,min=-4,max=4,step=0.5),
                variance = .slider(parent=mrv$sidecont3_v_entry,min=0.000001,max=10,step=0.1))
  })
  
  gSignalConnect(mrv$sidecont3_input1_savebt, "clicked", function(button) {
    rho.mean = as.numeric(gtkEntryGetText(mrv$sidecont3_m_entry))
    rho.variance = as.numeric(gtkEntryGetText(mrv$sidecont3_v_entry))
    mrv$cor.prior = "normal"
    mrv$cor.par = c(rho.mean,rho.variance)
    mrv$sidecont3_hide11$show()
    mrv$sidecont3_hide12$hide()
  })
  # PC
  gSignalConnect(mrv$sidecont3_sub1_figurebt, "clicked", function(button){
    U_min = as.numeric(gtkEntryGetText(mrv$sidecont3_sub1_U1_entry))
    a1 = as.numeric(gtkEntryGetText(mrv$sidecont3_sub1_A1_entry))

    .manipulate(parent=list(mrv$sidecont3_sub1_U1_entry,
                            mrv$sidecont3_sub1_A1_entry),
                .priorExpRhoS1(rho.ref=as.numeric(mrv$sidecont3_ref_entry$getText()),
                               left.portion = as.numeric(gtkEntryGetText(mrv$sidecont3_per_entry)),
                               U_min, a1),
                U_min = .slider(parent=mrv$sidecont3_sub1_U1_entry,
                        min=-0.99999,
                        max=as.numeric(mrv$sidecont3_ref_entry$getText())-0.05,
                        step=0.001, label = "U_min"), 
                a1 = .slider(parent=mrv$sidecont3_sub1_A1_entry,
                        min=0.001,max=as.numeric(gtkEntryGetText(mrv$sidecont3_per_entry))-0.000001,
                        step=0.001, label = "a1"))
  })
 
  gSignalConnect(mrv$sidecont3_sub2_figurebt, "clicked", function(button){
    U_max = as.numeric(gtkEntryGetText(mrv$sidecont3_sub2_U2_entry))
    a2 = as.numeric(gtkEntryGetText(mrv$sidecont3_sub2_A2_entry))
    .manipulate(parent=list(mrv$sidecont3_sub2_U2_entry,
                            mrv$sidecont3_sub2_A2_entry),
                .priorExpRhoS2(rho.ref=as.numeric(mrv$sidecont3_ref_entry$getText()),
                               left.portion = as.numeric(gtkEntryGetText(mrv$sidecont3_per_entry)),
                               U_max,a2),
                U_max = .slider(parent=mrv$sidecont3_sub2_U2_entry,
                                min=as.numeric(gtkEntryGetText(mrv$sidecont3_ref_entry))+0.05,
                                max=0.999999,
                                step=0.001),
                a2 = .slider(parent=mrv$sidecont3_sub2_A2_entry,
                                 min=0.001,max=1-as.numeric(gtkEntryGetText(mrv$sidecont3_per_entry))-0.000001,
                                 step=0.001))
  })
  
  gSignalConnect(mrv$sidecont3_sub3_figurebt, "clicked", function(button) {
    rho.ref = as.numeric(gtkEntryGetText(mrv$sidecont3_ref_entry))
    U_min = as.numeric(gtkEntryGetText(mrv$sidecont3_sub3_U1_entry))
    a1 = as.numeric(gtkEntryGetText(mrv$sidecont3_sub3_A1_entry))
    U_max = as.numeric(gtkEntryGetText(mrv$sidecont3_sub3_U2_entry))
    a2 = as.numeric(gtkEntryGetText(mrv$sidecont3_sub3_A2_entry))
    .manipulate(parent=list(mrv$sidecont3_sub3_U1_entry,mrv$sidecont3_sub3_A1_entry,
                            mrv$sidecont3_sub3_U2_entry,mrv$sidecont3_sub3_A2_entry),
                .priorExpRhoS3(rho.ref=as.numeric(gtkEntryGetText(mrv$sidecont3_ref_entry)),
                               U_min,a1,U_max,a2),
                U_min = .slider(parent=mrv$sidecont3_sub3_U1_entry,
                                min=-0.99999,
                                max=as.numeric(gtkEntryGetText(mrv$sidecont3_ref_entry))-0.1,
                                step=0.001),
                a1 = .slider(parent=mrv$sidecont3_sub3_A1_entry,min=0.005,max=0.3,step=0.005),
                U_max = .slider(parent=mrv$sidecont3_sub3_U2_entry,
                                min=as.numeric(gtkEntryGetText(mrv$sidecont3_ref_entry))+0.1,
                                max=0.999999,
                                step=0.001),
                a2 = .slider(parent=mrv$sidecont3_sub3_A2_entry,min=0.005,max=0.3,step=0.005))
  })
  
  gSignalConnect(mrv$sidecont3_sub1_savebt, "clicked", function(button) {
    rho.ref = as.numeric(gtkEntryGetText(mrv$sidecont3_ref_entry))
    rho.per = as.numeric(gtkEntryGetText(mrv$sidecont3_per_entry))
    Umin = as.numeric(gtkEntryGetText(mrv$sidecont3_sub1_U1_entry))
    alpha1 = as.numeric(gtkEntryGetText(mrv$sidecont3_sub1_A1_entry))
    mrv$cor.prior = "PC"
    mrv$cor.par = c(1L,rho.ref,rho.per,Umin,alpha1,NA,NA)
    mrv$sidecont3_input2_subhide11$show()
    mrv$sidecont3_input2_subhide12$hide()
  })
  
  gSignalConnect(mrv$sidecont3_sub2_savebt, "clicked", function(button) {
    rho.ref = as.numeric(gtkEntryGetText(mrv$sidecont3_ref_entry))
    rho.per = as.numeric(gtkEntryGetText(mrv$sidecont3_per_entry))
    Umax = as.numeric(gtkEntryGetText(mrv$sidecont3_sub2_U2_entry))
    alpha2 = as.numeric(gtkEntryGetText(mrv$sidecont3_sub2_A2_entry))
    mrv$cor.prior = "PC"
    mrv$cor.par = c(2L,rho.ref,rho.per,NA,NA,Umax,alpha2)
    mrv$sidecont3_input2_subhide21$show()
    mrv$sidecont3_input2_subhide22$hide()
  })
  
  gSignalConnect(mrv$sidecont3_sub3_savebt, "clicked", function(button) {
    rho.ref = as.numeric(gtkEntryGetText(mrv$sidecont3_ref_entry))
    Umin = as.numeric(gtkEntryGetText(mrv$sidecont3_sub3_U1_entry))
    alpha1 = as.numeric(gtkEntryGetText(mrv$sidecont3_sub3_A1_entry))
    Umax = as.numeric(gtkEntryGetText(mrv$sidecont3_sub3_U2_entry))
    alpha2 = as.numeric(gtkEntryGetText(mrv$sidecont3_sub3_A2_entry))
    mrv$cor.prior = "PC"
    mrv$cor.par = c(3L,rho.ref,NA,Umin,alpha1,Umax,alpha2)
    mrv$sidecont3_input2_subhide31$show()
    mrv$sidecont3_input2_subhide32$hide()
  })
  # Beta
  gSignalConnect(mrv$sidecont3_input3_figurebt, "clicked", function(button){
    a = as.numeric(gtkEntryGetText(mrv$sidecont3_a_entry))
    b = as.numeric(gtkEntryGetText(mrv$sidecont3_b_entry))
    .manipulate(parent=list(mrv$sidecont3_a_entry,mrv$sidecont3_b_entry),
                .priorRhoBetaPlot(a,b),
                a = .slider(parent=mrv$sidecont3_a_entry,min=0,max=20,step=0.5),
                b = .slider(parent=mrv$sidecont3_b_entry,min=0,max=20,step=0.5))
  })
  
  gSignalConnect(mrv$sidecont3_input3_savebt, "clicked", function(button) {
    rho.a = as.numeric(gtkEntryGetText(mrv$sidecont3_a_entry))
    rho.b = as.numeric(gtkEntryGetText(mrv$sidecont3_b_entry))
    mrv$cor.prior = "beta"
    mrv$cor.par = c(rho.a,rho.b)
    mrv$sidecont3_hide31$show()
    mrv$sidecont3_hide32$hide()
  })
  
  # Table
  gSignalConnect(mrv$sidecont3_input4_figurebt, "clicked", function(button,...) {
    cordensity = gtkComboBoxGetActiveText(mrv$sidecont3_combo)
    if(cordensity=="Table(cor)"){
      .open_prior_Rho(window = mrv$main_window)
    }
  })
  
  gSignalConnect(mrv$sidecont3_input4_savebt, "clicked", function(button) {
    mrv$cor.prior = "Table"
    mrv$cor.par = mrv$priorfile
    mrv$sidecont3_hide4$hide()
  })
  
  
  ##################################################################
  ###                 component 4 - wishart parameters
  ##################################################################
  mrv$sidecont4 <- gtkFrame("Inv-Wishart Parameter")
  mrv$sidecont4["border-width"]=10
  mrv$sidecont4$hide()
  
  mrv$sidecont4_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont4_main["border-width"]=10

  mrv$sidecont4_nu_label <- sidebarLabel("d.o.f: ",1)
  mrv$sidecont4_R11_label <- sidebarLabel("R11: ",1)
  mrv$sidecont4_R22_label <- sidebarLabel("R22: ",1)
  mrv$sidecont4_R12_label <- sidebarLabel("R12: ",1)
  mrv$sidecont4_nu_entry <- sidebarEntry("4", 10, 0)
  mrv$sidecont4_R11_entry <- sidebarEntry("1", 10, 0)
  mrv$sidecont4_R22_entry <- sidebarEntry("1", 10, 0)
  mrv$sidecont4_R12_entry <- sidebarEntry("0", 10, 0)
  sidecont4_list <- list(mrv$sidecont4_nu_label, mrv$sidecont4_nu_entry, mrv$sidecont4_R11_label, mrv$sidecont4_R11_entry,
                         mrv$sidecont4_R22_label, mrv$sidecont4_R22_entry, mrv$sidecont4_R12_label, mrv$sidecont4_R12_entry)
  
  mrv$sidecont4_input <- sidebarTable(3, 3, sidecont4_list)
  
  mrv$sidecont4_input_buttonbox <- gtkHBox(homogeneous=TRUE,spacing=10)
  mrv$sidecont4_input_figurebt <- sidebarPlotButton()
  mrv$sidecont4_input_savebt <- sidebarAcceptButton()
  mrv$sidecont4_input_buttonbox$packStart(mrv$sidecont4_input_figurebt,expand=T,fill=T)
  mrv$sidecont4_input_buttonbox$packStart(mrv$sidecont4_input_savebt,expand=T,fill=T)
  
  mrv$sidecont4_hide1 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont4_hide1$show()
  mrv$sidecont4_hide2 <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidecont4_hide2$show()
  
  mrv$sidecont4_hide1$packStart(mrv$sidecont4_input,expand=FALSE,fill=FALSE)
  mrv$sidecont4_hide2$packStart(mrv$sidecont4_input_buttonbox,expand=FALSE,fill=FALSE)
  
  mrv$sidecont4_main$packStart(mrv$sidecont4_hide1,expand=FALSE,fill=FALSE)
  mrv$sidecont4_main$packStart(mrv$sidecont4_hide2,expand=FALSE,fill=FALSE)
  mrv$sidecont4$add(mrv$sidecont4_main)
  
  gSignalConnect(mrv$sidecont4_nu_entry, "changed", function(entry){
    mrv$sidecont4_hide2$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont4_R11_entry, "changed", function(entry){
    mrv$sidecont4_hide2$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont4_R22_entry, "changed", function(entry){
    mrv$sidecont4_hide2$show()
    .checkNumEntry(entry)
  })
  
  gSignalConnect(mrv$sidecont4_R12_entry, "changed", function(entry){
    mrv$sidecont4_hide2$show()
    .checkNumEntry(entry)
  })
  
  # Inv-gamma
  gSignalConnect(mrv$sidecont4_input_figurebt, "clicked", function(button) {
    nu = as.numeric(gtkEntryGetText(mrv$sidecont4_nu_entry))
    R11 = as.numeric(gtkEntryGetText(mrv$sidecont4_R11_entry))
    R22 = as.numeric(gtkEntryGetText(mrv$sidecont4_R22_entry))
    R12 = as.numeric(gtkEntryGetText(mrv$sidecont4_R12_entry))
    xlim = 8
    .manipulate(parent=list(mrv$sidecont4_nu_entry,mrv$sidecont4_R11_entry,mrv$sidecont4_R22_entry,mrv$sidecont4_R12_entry),
                .priorInvWishart(nu,R11,R22,R12,xlim),
                nu = .slider(parent=mrv$sidecont4_nu_entry,label="nu",min=4,max=20,step=1),
                R11 = .slider(parent=mrv$sidecont4_R11_entry,label="R11",min=0,max=10,step=0.1),
                R22 = .slider(parent=mrv$sidecont4_R22_entry,label="R22",min=0,max=10,step=0.1),
                R12 = .slider(parent=mrv$sidecont4_R12_entry,label="R12",min=0,max=10,step=0.1),
                xlim = .slider(parent=hide_parent,label="xlim",min=1,max=20,step=1))
  })
  
  piw.nu = as.numeric(gtkEntryGetText(mrv$sidecont4_nu_entry))
  piw.R11 = as.numeric(gtkEntryGetText(mrv$sidecont4_R11_entry))
  piw.R22 = as.numeric(gtkEntryGetText(mrv$sidecont4_R22_entry))
  piw.R12 = as.numeric(gtkEntryGetText(mrv$sidecont4_R12_entry))
  mrv$wishart.par = c(piw.nu, piw.R11, piw.R22, piw.R12)
  
  gSignalConnect(mrv$sidecont4_input_savebt, "clicked", function(button) {
    piw.nu = as.numeric(gtkEntryGetText(mrv$sidecont4_nu_entry))
    piw.R11 = as.numeric(gtkEntryGetText(mrv$sidecont4_R11_entry))
    piw.R22 = as.numeric(gtkEntryGetText(mrv$sidecont4_R22_entry))
    piw.R12 = as.numeric(gtkEntryGetText(mrv$sidecont4_R12_entry))
    mrv$wishart.par = c(piw.nu, piw.R11, piw.R22, piw.R12)
    mrv$var.prior = "Invwishart"
    mrv$var2.prior = "Invwishart"
    mrv$cor.prior = "Invwishart"
    mrv$sidecont4_hide1$show()
    mrv$sidecont4_hide2$hide()
  })
  
  #####################################################
  ##### pack all parts of priors together to side prior
  #####################################################
  mrv$sideprior$packStart(mrv$sidecont1,expand=FALSE,fill=FALSE)
  mrv$sideprior$packStart(mrv$sidecont2,expand=FALSE,fill=FALSE)
  mrv$sideprior$packStart(mrv$sidecont3,expand=FALSE,fill=FALSE)
  mrv$sideprior$packStart(mrv$sidecont4,expand=FALSE,fill=FALSE)
  
  
  #####################################################
  #####        All hiden signal
  #####################################################
  #### first component -  variance 1
  gSignalConnect(mrv$sidecont1_combo, "changed", function(button, ...) {
    tau2density = gtkComboBoxGetActiveText(mrv$sidecont2_combo)
    if(is.null(tau2density)){
      mrv$var2.prior = NULL
    }else{
      if(tau2density=="Inv-gamma(var)"){mrv$var2.prior = "Invgamma"}
      if(tau2density=="PC(std)"){mrv$var2.prior = "PC"}
      if(tau2density=="Half-Cauchy(std)"){mrv$var2.prior = "HCauchy"}
      if(tau2density=="Truncated-normal(std)"){mrv$var2.prior = "Tnormal"}
      if(tau2density=="Uniform(std)"){mrv$var2.prior = "Unif"}
      if(tau2density=="Table(var)"){mrv$var2.prior = "Table"}
    }
    
    
    cordensity = gtkComboBoxGetActiveText(mrv$sidecont3_combo)
    if(is.null(cordensity)){
      mrv$cor.prior = NULL
    }else{
      if(cordensity=="Gaussian(Fisher-z)"){mrv$cor.prior = "normal"}
      if(cordensity=="PC(cor)"){mrv$cor.prior = "PC"}
      if(cordensity=="Beta(transf.cor)"){mrv$cor.prior = "beta"}
      if(cordensity=="Table(cor)"){mrv$cor.prior = "Table"}
    }
    
    tau1density = gtkComboBoxGetActiveText(button)
    if(tau1density=="Inv-gamma(var)"){
      mrv$var.prior = "Invgamma"
      mrv$sidecont1_hide11$show()
      mrv$sidecont1_hide12$show() 
      mrv$sidecont1_hide21$hide()
      mrv$sidecont1_hide22$hide() 
      mrv$sidecont1_hide31$hide()
      mrv$sidecont1_hide32$hide() 
      mrv$sidecont1_hide41$hide()
      mrv$sidecont1_hide42$hide()
      mrv$sidecont1_hide5$hide()
      mrv$sidecont1_hide6$hide()
      mrv$sidecont2$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(tau1density=="PC(std)"){
      mrv$var.prior = "PC"
      mrv$sidecont1_hide11$hide()
      mrv$sidecont1_hide12$hide()
      mrv$sidecont1_hide21$show()
      mrv$sidecont1_hide22$show()
      mrv$sidecont1_hide31$hide()
      mrv$sidecont1_hide32$hide() 
      mrv$sidecont1_hide41$hide()
      mrv$sidecont1_hide42$hide()
      mrv$sidecont1_hide5$hide()
      mrv$sidecont1_hide6$hide()
      mrv$sidecont2$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(tau1density=="Half-Cauchy(std)"){
      mrv$var.prior = "HCauchy"
      mrv$sidecont1_hide11$hide()
      mrv$sidecont1_hide12$hide()
      mrv$sidecont1_hide21$hide()
      mrv$sidecont1_hide22$hide()
      mrv$sidecont1_hide31$show()
      mrv$sidecont1_hide32$show() 
      mrv$sidecont1_hide41$hide()
      mrv$sidecont1_hide42$hide()
      mrv$sidecont1_hide5$hide()
      mrv$sidecont1_hide6$hide()
      mrv$sidecont2$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(tau1density=="Truncated-normal(std)"){
      mrv$var.prior = "Tnormal"
      mrv$sidecont1_hide11$hide()
      mrv$sidecont1_hide12$hide()
      mrv$sidecont1_hide21$hide()
      mrv$sidecont1_hide22$hide()
      mrv$sidecont1_hide31$hide()
      mrv$sidecont1_hide32$hide() 
      mrv$sidecont1_hide41$show()
      mrv$sidecont1_hide42$show()
      mrv$sidecont1_hide5$hide()
      mrv$sidecont1_hide6$hide()
      mrv$sidecont2$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(tau1density=="Uniform(std)"){
      mrv$var.prior = "Unif"
      mrv$sidecont1_hide11$hide()
      mrv$sidecont1_hide12$hide()
      mrv$sidecont1_hide21$hide()
      mrv$sidecont1_hide22$hide()
      mrv$sidecont1_hide31$hide()
      mrv$sidecont1_hide32$hide() 
      mrv$sidecont1_hide41$hide()
      mrv$sidecont1_hide42$hide()
      mrv$sidecont1_hide5$show()
      mrv$sidecont1_hide6$hide()
      mrv$sidecont2$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(tau1density=="Table(var)"){
      mrv$var.prior = "Table"
      mrv$sidecont1_hide11$hide()
      mrv$sidecont1_hide12$hide()
      mrv$sidecont1_hide21$hide()
      mrv$sidecont1_hide22$hide()
      mrv$sidecont1_hide31$hide()
      mrv$sidecont1_hide32$hide() 
      mrv$sidecont1_hide41$hide()
      mrv$sidecont1_hide42$hide()
      mrv$sidecont1_hide5$hide()
      mrv$sidecont1_hide6$show()
      mrv$sidecont2$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else{
      mrv$cor.prior = "Invwishart"
      mrv$var.prior = "Invwishart"
      mrv$var2.prior = "Invwishart"
      mrv$sidecont1_hide11$hide()
      mrv$sidecont1_hide12$hide()
      mrv$sidecont1_hide21$hide()
      mrv$sidecont1_hide22$hide()
      mrv$sidecont1_hide31$hide()
      mrv$sidecont1_hide32$hide() 
      mrv$sidecont1_hide41$hide()
      mrv$sidecont1_hide42$hide()
      mrv$sidecont1_hide5$hide()
      mrv$sidecont1_hide6$hide()
      mrv$sidecont2$hide()
      mrv$sidecont3$hide()
      mrv$sidecont4$show()
      mrv$sidecont4_hide1$show()
      mrv$sidecont4_hide2$show()
    }
  })
  #### second component -  variance 2
  gSignalConnect(mrv$sidecont2_combo, "changed", function(button, ...) {
    tau1density = gtkComboBoxGetActiveText(mrv$sidecont1_combo)
    if(is.null(tau1density)){
      mrv$var.prior = NULL
    }else{
      if(tau1density=="Inv-gamma(var)"){mrv$var.prior = "Invgamma"}
      if(tau1density=="PC(std)"){mrv$var.prior = "PC"}
      if(tau1density=="Half-Cauchy(std)"){mrv$var.prior = "HCauchy"}
      if(tau1density=="Truncated-normal(std)"){mrv$var.prior = "Tnormal"}
      if(tau1density=="Uniform(std)"){mrv$var.prior = "Unif"}
      if(tau1density=="Table(var)"){mrv$var.prior = "Table"}
    }
    
    
    cordensity = gtkComboBoxGetActiveText(mrv$sidecont3_combo)
    if(is.null(cordensity)){
      mrv$cor.prior = NULL
    }else{
      if(cordensity=="Gaussian(Fisher-z)"){mrv$cor.prior = "normal"}
      if(cordensity=="PC(cor)"){mrv$cor.prior = "PC"}
      if(cordensity=="Beta(transf.cor)"){mrv$cor.prior = "beta"}
      if(cordensity=="Table(cor)"){mrv$cor.prior = "Table"}
    }
    
    tau2density = gtkComboBoxGetActiveText(button)
    if(tau2density=="Inv-gamma(var)"){
      mrv$var2.prior = "Invgamma"
      mrv$sidecont2_hide11$show()
      mrv$sidecont2_hide12$show() 
      mrv$sidecont2_hide21$hide()
      mrv$sidecont2_hide22$hide() 
      mrv$sidecont2_hide31$hide()
      mrv$sidecont2_hide32$hide() 
      mrv$sidecont2_hide41$hide()
      mrv$sidecont2_hide42$hide()
      mrv$sidecont2_hide5$hide()
      mrv$sidecont2_hide6$hide()
      mrv$sidecont1$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(tau2density=="PC(std)"){
      mrv$var2.prior = "PC"
      mrv$sidecont2_hide11$hide()
      mrv$sidecont2_hide12$hide()
      mrv$sidecont2_hide21$show()
      mrv$sidecont2_hide22$show()
      mrv$sidecont2_hide31$hide()
      mrv$sidecont2_hide32$hide() 
      mrv$sidecont2_hide41$hide()
      mrv$sidecont2_hide42$hide()
      mrv$sidecont2_hide5$hide()
      mrv$sidecont2_hide6$hide()
      mrv$sidecont1$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(tau2density=="Half-Cauchy(std)"){
      mrv$var2.prior = "HCauchy"
      mrv$sidecont2_hide11$hide()
      mrv$sidecont2_hide12$hide()
      mrv$sidecont2_hide21$hide()
      mrv$sidecont2_hide22$hide()
      mrv$sidecont2_hide31$show()
      mrv$sidecont2_hide32$show() 
      mrv$sidecont2_hide41$hide()
      mrv$sidecont2_hide42$hide()
      mrv$sidecont2_hide5$hide()
      mrv$sidecont2_hide6$hide()
      mrv$sidecont1$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(tau2density=="Truncated-normal(std)"){
      mrv$var2.prior = "Tnormal"
      mrv$sidecont2_hide11$hide()
      mrv$sidecont2_hide12$hide()
      mrv$sidecont2_hide21$hide()
      mrv$sidecont2_hide22$hide()
      mrv$sidecont2_hide31$hide()
      mrv$sidecont2_hide32$hide() 
      mrv$sidecont2_hide41$show()
      mrv$sidecont2_hide42$show()
      mrv$sidecont2_hide5$hide()
      mrv$sidecont2_hide6$hide()
      mrv$sidecont1$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(tau2density=="Uniform(std)"){
      mrv$var2.prior = "Unif"
      mrv$sidecont2_hide11$hide()
      mrv$sidecont2_hide12$hide()
      mrv$sidecont2_hide21$hide()
      mrv$sidecont2_hide22$hide()
      mrv$sidecont2_hide31$hide()
      mrv$sidecont2_hide32$hide() 
      mrv$sidecont2_hide41$hide()
      mrv$sidecont2_hide42$hide()
      mrv$sidecont2_hide5$show()
      mrv$sidecont2_hide6$hide()
      mrv$sidecont1$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(tau2density=="Table(var)"){
      mrv$var2.prior = "Table"
      mrv$sidecont2_hide11$hide()
      mrv$sidecont2_hide12$hide()
      mrv$sidecont2_hide21$hide()
      mrv$sidecont2_hide22$hide()
      mrv$sidecont2_hide31$hide()
      mrv$sidecont2_hide32$hide() 
      mrv$sidecont2_hide41$hide()
      mrv$sidecont2_hide42$hide()
      mrv$sidecont2_hide5$hide()
      mrv$sidecont2_hide6$show()
      mrv$sidecont1$show()
      mrv$sidecont3$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else{
      mrv$cor.prior = "Invwishart"
      mrv$var.prior = "Invwishart"
      mrv$var2.prior = "Invwishart"
      mrv$sidecont2_hide11$hide()
      mrv$sidecont2_hide12$hide() 
      mrv$sidecont2_hide21$hide()
      mrv$sidecont2_hide22$hide() 
      mrv$sidecont2_hide31$hide()
      mrv$sidecont2_hide32$hide() 
      mrv$sidecont2_hide41$hide()
      mrv$sidecont2_hide42$hide()
      mrv$sidecont2_hide5$hide()
      mrv$sidecont2_hide6$hide()
      mrv$sidecont1$hide()
      mrv$sidecont3$hide()
      mrv$sidecont4$show()
      mrv$sidecont4_hide1$show()
      mrv$sidecont4_hide2$show()
    }
  })
  #### third component -  correlation
  gSignalConnect(mrv$sidecont3_combo, "changed", function(button, ...) {
    tau1density = gtkComboBoxGetActiveText(mrv$sidecont1_combo)
    if(is.null(tau1density)){
      mrv$var.prior = NULL
    }else{
      if(tau1density=="Inv-gamma(var)"){mrv$var.prior = "Invgamma"}
      if(tau1density=="PC(std)"){mrv$var.prior = "PC"}
      if(tau1density=="Half-Cauchy(std)"){mrv$var.prior = "HCauchy"}
      if(tau1density=="Truncated-normal(std)"){mrv$var.prior = "Tnormal"}
      if(tau1density=="Uniform(std)"){mrv$var.prior = "Unif"}
      if(tau1density=="Table(var)"){mrv$var.prior = "Table"}
    }
    
    
    tau2density = gtkComboBoxGetActiveText(mrv$sidecont2_combo)
    if(is.null(tau2density)){
      mrv$var2.prior = NULL
    }else{
      if(tau2density=="Inv-gamma(var)"){mrv$var2.prior = "Invgamma"}
      if(tau2density=="PC(std)"){mrv$var2.prior = "PC"}
      if(tau2density=="Half-Cauchy(std)"){mrv$var2.prior = "HCauchy"}
      if(tau2density=="Truncated-normal(std)"){mrv$var2.prior = "Tnormal"}
      if(tau2density=="Uniform(std)"){mrv$var2.prior = "Unif"}
      if(tau2density=="Table(var)"){mrv$var2.prior = "Table"}
    }
    
    
    rhodensity = gtkComboBoxGetActiveText(button)
    if(rhodensity=="Gaussian(Fisher-z)"){
      mrv$cor.prior = "normal"
      mrv$sidecont3_hide11$show()
      mrv$sidecont3_hide12$show()  
      mrv$sidecont3_hide2$hide()
      mrv$sidecont3_hide31$hide()
      mrv$sidecont3_hide32$hide()
      mrv$sidecont3_hide4$hide()
      mrv$sidecont1$show()
      mrv$sidecont2$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(rhodensity=="PC(cor)"){
      mrv$cor.prior = "PC"
      mrv$sidecont3_hide11$hide()
      mrv$sidecont3_hide12$hide()
      mrv$sidecont3_hide2$show()
      mrv$sidecont3_hide31$hide()
      mrv$sidecont3_hide32$hide()
      mrv$sidecont3_hide4$hide()
      mrv$sidecont1$show()
      mrv$sidecont2$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(rhodensity=="Beta(transf.cor)"){
      mrv$cor.prior = "beta"
      mrv$sidecont3_hide11$hide()
      mrv$sidecont3_hide12$hide()
      mrv$sidecont3_hide2$hide()
      mrv$sidecont3_hide31$show()
      mrv$sidecont3_hide32$show()
      mrv$sidecont3_hide4$hide()
      mrv$sidecont1$show()
      mrv$sidecont2$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else if(rhodensity=="Table(cor)"){
      mrv$cor.prior = "Table"
      mrv$sidecont3_hide11$hide()
      mrv$sidecont3_hide12$hide()
      mrv$sidecont3_hide2$hide()
      mrv$sidecont3_hide31$hide()
      mrv$sidecont3_hide32$hide()
      mrv$sidecont3_hide4$show()
      mrv$sidecont1$show()
      mrv$sidecont2$show()
      mrv$sidecont4$hide()
      mrv$sidecont4_hide1$hide()
      mrv$sidecont4_hide2$hide()
    }else{
      mrv$cor.prior = "Invwishart"
      mrv$var.prior = "Invwishart"
      mrv$var2.prior = "Invwishart"
      mrv$sidecont3_hide11$hide()
      mrv$sidecont3_hide12$hide()  
      mrv$sidecont3_hide2$hide()
      mrv$sidecont3_hide31$hide()
      mrv$sidecont3_hide32$hide()
      mrv$sidecont3_hide4$hide()
      mrv$sidecont1$show()
      mrv$sidecont2$show()
      mrv$sidecont1$hide()
      mrv$sidecont2$hide()
      mrv$sidecont4$show()
      mrv$sidecont4_hide1$show()
      mrv$sidecont4_hide2$show()
    }
  })
  ##################################################################
  ##################################################################
  ######
  ######                   sidebar model
  ######
  ##################################################################
  ################################################################## 
  mrv$sidemodel <- gtkVBox(homogeneous=FALSE)
  mrv$sidemodel_title <- gtkImage(filename=mrv$sidebar_icon[[2]])    #### maybe change move file to folder
  mrv$sidemodel$packStart(mrv$sidemodel_title,expand=FALSE,fill=FALSE)
  
  mrv$sidemodel_window <- sidebarScrolledWindow(mrv$sidemodel)
  
  ############# sidebar model type
  mrv$sidemcont1 <- gtkFrame("Model Type")
  mrv$sidemcont1["border-width"]=10
  
  mrv$sidemcont1_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidemcont1_main["border-width"]=10
  
  mrv$sidem_mf_radio <- gtkVBox()
  mrv$sidem_mf_radiogp <- list()
  mrv$sidem_mf_radiogp$one <- gtkRadioButton(label = "Se & Sp")
  mrv$sidem_mf_radiogp$two <- gtkRadioButton(mrv$sidem_mf_radiogp, label = "Se & (1-Sp)")
  mrv$sidem_mf_radiogp$three <- gtkRadioButton(mrv$sidem_mf_radiogp, label = "(1-Se) & Sp")
  mrv$sidem_mf_radiogp$four <- gtkRadioButton(mrv$sidem_mf_radiogp, label = "(1-Se) & (1-Sp)")
  sapply(mrv$sidem_mf_radiogp, mrv$sidem_mf_radio$packStart)
  mrv$sidem_mf_radio[[1]]$setActive(TRUE) 
  sapply(mrv$sidem_mf_radiogp,'[',"active")
  mrv$sidem_mf_radio_align <- gtkAlignment(xalign = 0)
  mrv$sidem_mf_radio_align$add(mrv$sidem_mf_radio)
  mrv$model.type = 1
  ## pack model type
  mrv$sidemcont1_main$packStart(mrv$sidem_mf_radio_align, expand=FALSE,fill=FALSE)
  mrv$sidemcont1$add(mrv$sidemcont1_main)
  
  sapply(mrv$sidem_mf_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             if(button$getLabel()=="Se & Sp"){mrv$model.type = 1}
             if(button$getLabel()=="Se & (1-Sp)"){mrv$model.type = 2}
             if(button$getLabel()=="(1-Se) & Sp"){mrv$model.type = 3}
             if(button$getLabel()=="(1-Se) & (1-Sp)"){mrv$model.type = 4}
           } 
         })
  
  ############ sidebar link function
  mrv$sidemcont2 <- gtkFrame("Link Function")
  mrv$sidemcont2["border-width"]=10
  
  mrv$sidemcont2_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidemcont2_main["border-width"]=10
  
  mrv$sidem_lf_radio <- gtkVBox()
  mrv$sidem_lf_radiogp <- list()
  mrv$sidem_lf_radiogp$l <- gtkRadioButton(label = "logit")
  mrv$sidem_lf_radiogp$p <- gtkRadioButton(mrv$sidem_lf_radiogp, label = "probit")
  mrv$sidem_lf_radiogp$c <- gtkRadioButton(mrv$sidem_lf_radiogp, label = "cloglog")
  sapply(mrv$sidem_lf_radiogp, mrv$sidem_lf_radio$packStart)
  mrv$sidem_lf_radio[[1]]$setActive(TRUE) 
  sapply(mrv$sidem_lf_radiogp,'[',"active")
  mrv$sidem_lf_radio_align <- gtkAlignment(xalign = 0)
  mrv$sidem_lf_radio_align$add(mrv$sidem_lf_radio)
  mrv$link = "logit"
  ## pack link function
  mrv$sidemcont2_main$packStart(mrv$sidem_lf_radio_align, expand=FALSE,fill=FALSE)
  mrv$sidemcont2$add(mrv$sidemcont2_main)
  
  sapply(mrv$sidem_lf_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             mrv$link = button$getLabel()
           } 
         })
  
  ############ sidebar level
  mrv$sidemcont3 <- gtkFrame("Quantiles")
  mrv$sidemcont3["border-width"]=10
  
  mrv$sidemcont3_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidemcont3_main["border-width"]=10
  
  mrv$sidemcont3_table <- gtkTable(rows=6,columns=2,homogeneous=FALSE)
  mrv$sidemcont3_table$setColSpacings(5)
  mrv$sidemcont3_table$setRowSpacings(3)
  
  mrv$sidemcont3_checkgp <- list()
  checkgp_label = c("0.05","0.1 ","0.15","0.2 ","0.25","0.75","0.8 ","0.85","0.9 ","0.95")
  for(i in 1:10){
    mrv$sidemcont3_checkgp[[i]] <- gtkCheckButton(checkgp_label[i])
    mrv$sidemcont3_checkgp[[i]]$setActive(FALSE)
  }
  mrv$sidemcont3_label <- gtkLabel("User Input: ")
  mrv$sidemcont3_label["xalign"]=1
  mrv$sidemcont3_entry <- gtkEntry()
  mrv$sidemcont3_entry["text"]=""
  mrv$sidemcont3_entry["width-chars"]=10
  mrv$sidemcont3_entry["xalign"]=0
  
  for(i in 1:5){
    check_ind1 = 2*i-1
    check_ind2 = 2*i
    mrv$sidemcont3_table$attach(mrv$sidemcont3_checkgp[[check_ind1]],left.attach=0,1, top.attach=(i-1),i,
                                xoptions=c("expand","fill"),yoptions="")
    mrv$sidemcont3_table$attach(mrv$sidemcont3_checkgp[[check_ind2]],left.attach=1,2, top.attach=(i-1),i,
                                xoptions="",yoptions="")
  }
  mrv$sidemcont3_table$attach(mrv$sidemcont3_label,left.attach=0,1, top.attach=5,6,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidemcont3_table$attach(mrv$sidemcont3_entry,left.attach=1,2, top.attach=5,6,
                              xoptions="",yoptions="")
  ## pack level
  mrv$sidemcont3_main$packStart(mrv$sidemcont3_table, expand=FALSE,fill=FALSE)
  mrv$sidemcont3$add(mrv$sidemcont3_main)
  
  mrv$flevel = rep(list(NA),11)
  sapply(c(1:10), function(ind){
    gSignalConnect (mrv$sidemcont3_checkgp[[ind]], "toggled", f = function(button, ...){
      if(button['active']){
        mrv$flevel[[ind]] = button$getLabel()
      }else{
        mrv$flevel[[ind]] = NA
      }
    })
  })
  
  gSignalConnect(mrv$sidemcont3_entry, "changed", function(entry){
    mrv$flevel[[11]] = gtkEntryGetText(mrv$sidemcont3_entry)
    .checkNumEntry(entry)
  })
  
  
  ############ sidebar verbose
  mrv$sidemcont4 <- gtkFrame("Verbose")
  mrv$sidemcont4["border-width"]=10
  
  mrv$sidemcont4_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidemcont4_main["border-width"]=10
  
  mrv$sidem_v_radio <- gtkHBox()
  mrv$sidem_v_radiogp <- list()
  mrv$sidem_v_radiogp$t <- gtkRadioButton(label = "TRUE")
  mrv$sidem_v_radiogp$f <- gtkRadioButton(mrv$sidem_v_radiogp, label = "FALSE")
  sapply(mrv$sidem_v_radiogp, mrv$sidem_v_radio$packStart)
  mrv$sidem_v_radio[[2]]$setActive(TRUE) 
  sapply(mrv$sidem_v_radiogp,'[',"active")
  mrv$sidem_v_radio_align <- gtkAlignment(xalign = 0)
  mrv$sidem_v_radio_align$add(mrv$sidem_v_radio)
  
  ## pack verbose
  mrv$sidemcont4_main$packStart(mrv$sidem_v_radio_align, expand=FALSE,fill=FALSE)
  mrv$sidemcont4$add(mrv$sidemcont4_main)
  
  mrv$verbose = FALSE
  ### signals
  sapply(mrv$sidem_v_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             if(mrv$verbose)
               mrv$verbose = as.logical(button$getLabel())
           } 
         })
  ############ sidebar nsample
  mrv$sidemcont5 <- gtkFrame("Number of Samples")
  mrv$sidemcont5["border-width"]=10
  
  mrv$sidemcont5_main <- gtkHBox(homogeneous=FALSE,spacing=10)
  mrv$sidemcont5_main["border-width"]=10
  
  mrv$sidemcont5_button <- gtkHScale(min=1000, max=1000000, step=1000)
  mrv$sidemcont5_button['draw-value'] <- FALSE
  mrv$sidemcont5_button$setValuePos("right")
  adjustment <- mrv$sidemcont5_button$getAdjustment( )
  mrv$sidemcont5_spinbutton <- gtkSpinButton(adjustment = adjustment)
  
  ## pack verbose
  mrv$sidemcont5_main$packStart(mrv$sidemcont5_button, expand=T,fill=T,padding=5)
  mrv$sidemcont5_main$packStart(mrv$sidemcont5_spinbutton, expand=FALSE,fill=FALSE,padding=5)
  mrv$sidemcont5$add(mrv$sidemcont5_main)
  
  #gSignalConnect(mrv$sidemcont5_spinbutton, "value-changed", function(button){
  #  mrv$nsample = button$value
  #})
  ##### pack all parts of priors together to side prior
  mrv$sidemodel$packStart(mrv$sidemcont1,expand=FALSE,fill=FALSE)
  mrv$sidemodel$packStart(mrv$sidemcont2,expand=FALSE,fill=FALSE)
  mrv$sidemodel$packStart(mrv$sidemcont3,expand=FALSE,fill=FALSE)
  #mrv$sidemodel$packStart(mrv$sidemcont4,expand=FALSE,fill=FALSE)
  mrv$sidemodel$packStart(mrv$sidemcont5,expand=FALSE,fill=FALSE)
  
  ##################################################################
  ##################################################################
  ######
  ######                   sidebar data
  ######
  ##################################################################
  ##################################################################
  mrv$sidedata <- gtkVBox(homogeneous=FALSE)
  
  mrv$sidedata_window <- sidebarScrolledWindow(mrv$sidedata)
  
  mrv$sidedata_title <- gtkImage(filename=mrv$sidebar_icon[[3]])    #### maybe change move file to folder
  mrv$sidedata$packStart(mrv$sidedata_title,expand=FALSE,fill=FALSE)
  
  if(is.null(mrv$datafile)){
    mrv$sidedcont <- gtkButton("Load Data")
    mrv$sidedcont["border-width"]=10
    
    gSignalConnect(mrv$sidedcont, "clicked", function(button) {
      .open_cb(window=mrv$main_window)
    })
    mrv$sidedata$packStart(mrv$sidedcont,expand=FALSE,fill=FALSE)
  }
  
  ##################################################################
  ##################################################################
  ######
  ######                   sidebar SROC
  ######
  ##################################################################
  ##################################################################
  mrv$sideplot <- gtkVBox(homogeneous=FALSE)
  mrv$sideplot_title <- gtkImage(filename=mrv$sidebar_icon[[4]])    #### maybe change move file to folder
  mrv$sideplot$packStart(mrv$sideplot_title,expand=FALSE,fill=FALSE)
  
  mrv$sideplot_window <- sidebarScrolledWindow(mrv$sideplot)
  
  ################ summary point ###########################
  mrv$sidepcont1 <- gtkFrame("Summary Point")
  mrv$sidepcont1["border-width"]=10
  
  mrv$sidepcont1_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidepcont1_main["border-width"]=10
  
  mrv$sidep_sroc_radio <- gtkHBox()
  mrv$sidep_sroc_radiogp <- list()
  mrv$sidep_sroc_radiogp$mean <- gtkRadioButton(label = "Mean")
  mrv$sidep_sroc_radiogp$median <- gtkRadioButton(mrv$sidep_sroc_radiogp, label = "Median")
  sapply(mrv$sidep_sroc_radiogp, mrv$sidep_sroc_radio$packStart)
  mrv$sidep_sroc_radio[[1]]$setActive(TRUE) 
  sapply(mrv$sidep_sroc_radiogp,'[',"active")
  mrv$sidep_sroc_radio_align <- gtkAlignment(xalign = 0)
  mrv$sidep_sroc_radio_align$add(mrv$sidep_sroc_radio)
  
  mrv$est.type = "mean"
  # pack radio to main box
  mrv$sidepcont1_main$packStart(mrv$sidep_sroc_radio_align, expand=FALSE,fill=FALSE)
  
  sapply(mrv$sidep_sroc_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             mrv$est.type=button$getLabel()
           } 
         })
  ####### other things in a table
  mrv$sidepcont1_table <- gtkTable(rows=3,columns=2,homogeneous=FALSE)
  mrv$sidepcont1_table$setColSpacings(5)
  mrv$sidepcont1_table$setRowSpacings(7)
  
  ### summary points symbol
  mrv$sidepcont1_type_label <- sidebarLabel("Summary Point Symbol:",1)
  
  mrv$sidepcont1_type_combo <- gtkPointsTypeComboBox(0)
  
  iter = mrv$sidepcont1_type_combo$GetActiveIter()
  value = mrv$sidepcont1_type_combo$GetModel()$GetValue(iter$iter,1)$value
  if(value %in% c("@","+","%","#","*","o","O")){
    mrv$sp.pch <- value
  }else{
    mrv$sp.pch <- as.numeric(value) 
  }
  gSignalConnect(mrv$sidepcont1_type_combo, "changed",function(cb){
    mrv$sp.pch <- getValueFromFigureCombo(cb) 
  })
  
  ### summary points size
  mrv$sidepcont1_size_label <- sidebarLabel("Summary Point Size:",1)
  
  mrv$sidepcont1_size_combo <- gtkPointsSizeComboBox_simple(0)
  
  iter = mrv$sidepcont1_size_combo$GetActiveIter()
  mrv$sp.cex = as.numeric(mrv$sidepcont1_size_combo$GetModel()$GetValue(iter$iter,1)$value)
  gSignalConnect(mrv$sidepcont1_size_combo, "changed",function(cb){
    mrv$sp.cex <- getValueFromFigureCombo(cb)
  })
  ### summary points color
  mrv$sidepcont1_color_label <- sidebarLabel("Summary Point Color:",1)
  mrv$sidepcont1_color_button <- gtkColorButton(gdkColorParse(palette()[2])$color)
  mrv$sp.col = palette()[2]
  gSignalConnect(mrv$sidepcont1_color_button, "color-set", function(button, ...) {
    mrv$sp.col = rgb(red=button$color$red, green=button$color$green, blue=button$color$blue, alpha=button$alpha, maxColorValue=65535)
  })
  ## pack table
  mrv$sidepcont1_table$attach(mrv$sidepcont1_type_label,left.attach=0,1, top.attach=0,1,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont1_table$attach(mrv$sidepcont1_type_combo,left.attach=1,2, top.attach=0,1,
                              xoptions="fill",yoptions="")
  mrv$sidepcont1_table$attach(mrv$sidepcont1_size_label,left.attach=0,1, top.attach=1,2,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont1_table$attach(mrv$sidepcont1_size_combo,left.attach=1,2, top.attach=1,2,
                              xoptions="fill",yoptions="")
  mrv$sidepcont1_table$attach(mrv$sidepcont1_color_label,left.attach=0,1, top.attach=2,3,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont1_table$attach(mrv$sidepcont1_color_button,left.attach=1,2, top.attach=2,3,
                              xoptions="fill",yoptions="")
  
  ## pack table to main box
  mrv$sidepcont1_main$packStart(mrv$sidepcont1_table, expand=FALSE,fill=FALSE)
  
  mrv$sidepcont1$add(mrv$sidepcont1_main)
  
  ################ data points #######################
  mrv$sidepcont2 <- gtkFrame("Data Points")
  mrv$sidepcont2["border-width"]=10
  
  mrv$sidepcont2_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidepcont2_main["border-width"]=10
  
  mrv$sidepcont2_radio_table <- gtkTable(rows=1,columns=3,homogeneous=FALSE)
  mrv$sidepcont2_radio_table$setColSpacings(5)
  mrv$sidepcont2_radio_table$setRowSpacings(7)
  
  mrv$sidepcont2_radiogp <- list()
  mrv$sidepcont2_radiogp$o <- gtkRadioButton(label = "Original")
  mrv$sidepcont2_radiogp$f <- gtkRadioButton(mrv$sidepcont2_radiogp,label = "Fitted")
  mrv$sidepcont2_radiogp$nshow <- gtkRadioButton(mrv$sidepcont2_radiogp, label = "Not Show")
  sapply(mrv$sidepcont2_radiogp,'[',"active")
  
  mrv$sidepcont2_radio_table$attach(mrv$sidepcont2_radiogp$o,left.attach=0,1, top.attach=0,1,
                                    xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont2_radio_table$attach(mrv$sidepcont2_radiogp$f,left.attach=1,2, top.attach=0,1,
                                    xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont2_radio_table$attach(mrv$sidepcont2_radiogp$nshow,left.attach=2,3, top.attach=0,1,
                                    xoptions=c("expand","fill"),yoptions="")

  # pack data radio table in main box
  mrv$sidepcont2_main$packStart(mrv$sidepcont2_radio_table, expand=FALSE,fill=FALSE)
  mrv$dataShow = "o"
  sapply(mrv$sidepcont2_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             if(button$getLabel()=="Original"){
               mrv$dataShow = "o"
             }else if(button$getLabel()=="Fitted"){
               mrv$dataShow = "f"
             }else{
               mrv$dataShow = F
             }
           } 
         })
  ####### other things in a table
  mrv$sidepcont2_table <- gtkTable(rows=3,columns=2,homogeneous=FALSE)
  mrv$sidepcont2_table$setColSpacings(5)
  mrv$sidepcont2_table$setRowSpacings(7)
  ### data points symbol
  mrv$sidepcont2_type_label <- sidebarLabel("Data Point Symbol:",1)
  
  mrv$sidepcont2_type_combo <- gtkPointsTypeComboBox(0)
  
  iter = mrv$sidepcont2_type_combo$GetActiveIter()
  value = mrv$sidepcont2_type_combo$GetModel()$GetValue(iter$iter,1)$value
  if(value %in% c("@","+","%","#","*","o","O")){
    mrv$data.pch <- value
  }else{
    mrv$data.pch <- as.numeric(value) 
  }
  gSignalConnect(mrv$sidepcont2_type_combo, "changed",function(cb){
    mrv$data.pch <- getValueFromFigureCombo(cb) 
  })
  
  ### data points size
  mrv$sidepcont2_size_label <- sidebarLabel("Data Point Size:",1)
  
  mrv$sidepcont2_size_combo <- gtkPointsSizeComboBox_bubble(0)
  
  iter = mrv$sidepcont2_size_combo$GetActiveIter()
  value = mrv$sidepcont2_size_combo$GetModel()$GetValue(iter$iter,1)$value
  if(value %in% c("bubble","scaled")){
    mrv$data.cex <- value
  }else{
    mrv$data.cex <- as.numeric(value) 
  }
  gSignalConnect(mrv$sidepcont2_size_combo, "changed",function(cb){
    mrv$data.cex <- getValueFromFigureCombo(cb)
  })
  
  ### data bubble color
  mrv$sidepcont2_color_label <- sidebarLabel("Data Point Color:",1)
  mrv$sidepcont2_color_button <- gtkColorButton(gdkColorParse("lawngreen")$color)
  mrv$sidepcont2_color_button[["width-request"]] = 70
  tempcol = mrv$sidepcont2_color_button$color
  tempalpha = mrv$sidepcont2_color_button$alpha
  mrv$data.col = rgb(red=tempcol$red, green=tempcol$green, blue=tempcol$blue, alpha=tempalpha, maxColorValue=65535)
  gSignalConnect(mrv$sidepcont2_color_button, "color-set", function(button, ...) {
    mrv$data.col = rgb(red=button$color$red, green=button$color$green, blue=button$color$blue, alpha=button$alpha, maxColorValue=65535)
  })
  ## pack table
  mrv$sidepcont2_table$attach(mrv$sidepcont2_type_label,left.attach=0,1, top.attach=0,1,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont2_table$attach(mrv$sidepcont2_type_combo,left.attach=1,2, top.attach=0,1,
                              xoptions="fill",yoptions="")
  mrv$sidepcont2_table$attach(mrv$sidepcont2_size_label,left.attach=0,1, top.attach=1,2,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont2_table$attach(mrv$sidepcont2_size_combo,left.attach=1,2, top.attach=1,2,
                              xoptions="fill",yoptions="")
  mrv$sidepcont2_table$attach(mrv$sidepcont2_color_label,left.attach=0,1, top.attach=2,3,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont2_table$attach(mrv$sidepcont2_color_button,left.attach=1,2, top.attach=2,3,
                              xoptions="fill",yoptions="")
  
  ## pack table to main box
  mrv$sidepcont2_main$packStart(mrv$sidepcont2_table, expand=FALSE,fill=FALSE)
  
  mrv$sidepcont2$add(mrv$sidepcont2_main)
  ################ sroc lines #######################
  mrv$sidepsrocline <- gtkFrame("SROC Line")
  mrv$sidepsrocline["border-width"]=10
  
  mrv$sidepsrocline_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidepsrocline_main["border-width"]=10
  
  mrv$sidepsrocline_radio <- gtkHBox()
  mrv$sidepsrocline_radiogp <- list()
  mrv$sidepsrocline_radiogp$mean <- gtkRadioButton(label = "Show")
  mrv$sidepsrocline_radiogp$median <- gtkRadioButton(mrv$sidepsrocline_radiogp, label = "Not Show")
  sapply(mrv$sidepsrocline_radiogp, mrv$sidepsrocline_radio$packStart)
  mrv$sidepsrocline_radio[[1]]$setActive(TRUE) 
  sapply(mrv$sidepsrocline_radiogp,'[',"active")
  mrv$sidepsrocline_radio_align <- gtkAlignment(xalign = 0)
  mrv$sidepsrocline_radio_align$add(mrv$sidepsrocline_radio)
  mrv$lineShow = T
  # pack Show or not in main box
  mrv$sidepsrocline_main$packStart(mrv$sidepsrocline_radio_align, expand=FALSE,fill=FALSE)
  
  sapply(mrv$sidepsrocline_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             if(button$getLabel()=="Show"){
               mrv$lineShow = T
             }else{
               mrv$lineShow = F
             }
           } 
         })
  
  # other things in a table
  mrv$sidepsrocline_table <- gtkTable(rows=4,columns=2,homogeneous=FALSE)
  mrv$sidepsrocline_table$setColSpacings(5)
  mrv$sidepsrocline_table$setRowSpacings(7)
  
  ### sroc function
  mrv$sidepsrocline_func_label <- sidebarLabel("SROC Function:",1)
  
  mrv$sidepsrocline_func_combo <- gtkComboBoxNewText()
  sapply(c("1","2","3","4","5"), mrv$sidepsrocline_func_combo$appendText)
  mrv$sidepsrocline_func_combo$setActive(0)
  
  mrv$sroc.type=1
  gSignalConnect(mrv$sidepsrocline_func_combo, "changed", function(button, ...) {
    mrv$sroc.type = as.numeric(gtkComboBoxGetActiveText(button))
  })
  
  ### sroc line type
  mrv$sidepsrocline_type_label <- sidebarLabel("Line Type:",1)
  
  mrv$sidepsrocline_type_combo <- gtkLineTypeComboBox(0)
  
  iter = mrv$sidepsrocline_type_combo$GetActiveIter()
  mrv$line.lty = as.numeric(mrv$sidepsrocline_type_combo$GetModel()$GetValue(iter$iter,1)$value)
  gSignalConnect(mrv$sidepsrocline_type_combo, "changed",function(cb){
    mrv$line.lty <- getValueFromFigureCombo(cb) 
  })

  
  ### sroc lines width
  mrv$sidepsrocline_width_label <- sidebarLabel("Line Width:",1)
  
  mrv$sidepsrocline_width_combo <- gtkLineWidthComboBox(0)
  
  iter = mrv$sidepsrocline_width_combo$GetActiveIter()
  mrv$line.lwd = as.numeric(mrv$sidepsrocline_width_combo$GetModel()$GetValue(iter$iter,1)$value)
  gSignalConnect(mrv$sidepsrocline_width_combo, "changed",function(cb){
    mrv$line.lwd <- getValueFromFigureCombo(cb) 
  })
  
  ############## sroc lines color
  mrv$sidepsrocline_color_label <- gtkLabel("Line Color:")
  mrv$sidepsrocline_color_label["xalign"]=1
  
  gdk_color <- gdkColorParse(palette()[1])$color
  mrv$sidepsrocline_color_button <- gtkColorButton(gdk_color)
  
  tempcol = mrv$sidepsrocline_color_button$color
  tempalpha = mrv$sidepsrocline_color_button$alpha
  mrv$line.col = rgb(red=tempcol$red, green=tempcol$green, blue=tempcol$blue, alpha=tempalpha, maxColorValue=65535)
  gSignalConnect(mrv$sidepsrocline_color_button, "color-set", function(button, ...) {
    mrv$line.col = rgb(red=button$color$red, green=button$color$green, blue=button$color$blue, alpha=button$alpha, maxColorValue=65535)
  })
  ## pack table
  mrv$sidepsrocline_table$attach(mrv$sidepsrocline_func_label,left.attach=0,1, top.attach=0,1,
                                 xoptions=c("expand","fill"),yoptions="")
  mrv$sidepsrocline_table$attach(mrv$sidepsrocline_func_combo,left.attach=1,2, top.attach=0,1,
                                 xoptions="fill",yoptions="")
  mrv$sidepsrocline_table$attach(mrv$sidepsrocline_type_label,left.attach=0,1, top.attach=1,2,
                                 xoptions=c("expand","fill"),yoptions="")
  mrv$sidepsrocline_table$attach(mrv$sidepsrocline_type_combo,left.attach=1,2, top.attach=1,2,
                                 xoptions="fill",yoptions="")
  mrv$sidepsrocline_table$attach(mrv$sidepsrocline_width_label,left.attach=0,1, top.attach=2,3,
                                 xoptions=c("expand","fill"),yoptions="")
  mrv$sidepsrocline_table$attach(mrv$sidepsrocline_width_combo,left.attach=1,2, top.attach=2,3,
                                 xoptions="fill",yoptions="")
  mrv$sidepsrocline_table$attach(mrv$sidepsrocline_color_label,left.attach=0,1, top.attach=3,4,
                                 xoptions=c("expand","fill"),yoptions="")
  mrv$sidepsrocline_table$attach(mrv$sidepsrocline_color_button,left.attach=1,2, top.attach=3,4,
                                 xoptions="fill",yoptions="")
  
  ## pack line
  mrv$sidepsrocline_main$packStart(mrv$sidepsrocline_table, expand=FALSE,fill=FALSE)
  mrv$sidepsrocline$add(mrv$sidepsrocline_main)
  
  ################ confidence region #######################
  mrv$sidepcont4 <- gtkFrame("Condidence Region")
  mrv$sidepcont4["border-width"]=10
  
  mrv$sidepcont4_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidepcont4_main["border-width"]=10
  
  mrv$sidepcont4_radio <- gtkHBox()
  mrv$sidepcont4_radiogp <- list()
  mrv$sidepcont4_radiogp$bubbul <- gtkRadioButton(label = "Show")
  mrv$sidepcont4_radiogp$nshow <- gtkRadioButton(mrv$sidepcont4_radiogp, label = "Not Show")
  sapply(mrv$sidepcont4_radiogp, mrv$sidepcont4_radio$packStart)
  mrv$sidepcont4_radio[[1]]$setActive(TRUE) 
  sapply(mrv$sidepcont4_radiogp,'[',"active")
  mrv$sidepcont4_radio_align <- gtkAlignment(xalign = 0)
  mrv$sidepcont4_radio_align$add(mrv$sidepcont4_radio)
  
  mrv$crShow = T
  sapply(mrv$sidepcont4_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             if(button$getLabel()=="Show"){
               mrv$crShow = T
             }else{
               mrv$crShow = F
             }
           } 
         })
  
  # pack Show or not in main box
  mrv$sidepcont4_main$packStart(mrv$sidepcont4_radio_align, expand=FALSE,fill=FALSE)
  
  ####### other things in a table
  mrv$sidepcont4_table <- gtkTable(rows=3,columns=2,homogeneous=FALSE)
  mrv$sidepcont4_table$setColSpacings(5)
  mrv$sidepcont4_table$setRowSpacings(7)
  
  ### line type
  mrv$sidepcont4_type_label <- sidebarLabel("Line Type:",1)
  
  mrv$sidepcont4_type_combo <- gtkLineTypeComboBox(0)
  
  iter = mrv$sidepcont4_type_combo$GetActiveIter()
  mrv$cr.lty = as.numeric(mrv$sidepcont4_type_combo$GetModel()$GetValue(iter$iter,1)$value)
  gSignalConnect(mrv$sidepcont4_type_combo, "changed",function(cb){
    mrv$cr.lty <- getValueFromFigureCombo(cb)
  })
  ### lines width
  mrv$sidepcont4_width_label <- sidebarLabel("Line Width:",1)
  
  mrv$sidepcont4_width_combo <- gtkLineWidthComboBox(0)
  
  iter = mrv$sidepcont4_width_combo$GetActiveIter()
  mrv$cr.lwd = as.numeric(mrv$sidepcont4_width_combo$GetModel()$GetValue(iter$iter,1)$value)
  gSignalConnect(mrv$sidepcont4_width_combo, "changed",function(cb){
    mrv$cr.lwd <- getValueFromFigureCombo(cb) 
  })
  ### lines color
  mrv$sidepcont4_color_label <- sidebarLabel("Line Color:",1)
  
  gdk_color <- gdkColorParse("blue")$color
  mrv$sidepcont4_color_button <- gtkColorButton(gdk_color)
  
  tempcol = mrv$sidepcont4_color_button$color
  tempalpha = mrv$sidepcont4_color_button$alpha
  mrv$cr.col = rgb(red=tempcol$red, green=tempcol$green, blue=tempcol$blue, alpha=tempalpha, maxColorValue=65535)
  gSignalConnect(mrv$sidepcont4_color_button, "color-set", function(button, ...) {
    mrv$cr.col = rgb(red=button$color$red, green=button$color$green, blue=button$color$blue, alpha=button$alpha, maxColorValue=65535)
  })
  
  ## pack table
  mrv$sidepcont4_table$attach(mrv$sidepcont4_type_label,left.attach=0,1, top.attach=0,1,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont4_table$attach(mrv$sidepcont4_type_combo,left.attach=1,2, top.attach=0,1,
                              xoptions="fill",yoptions="")
  mrv$sidepcont4_table$attach(mrv$sidepcont4_width_label,left.attach=0,1, top.attach=1,2,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont4_table$attach(mrv$sidepcont4_width_combo,left.attach=1,2, top.attach=1,2,
                              xoptions="fill",yoptions="")
  mrv$sidepcont4_table$attach(mrv$sidepcont4_color_label,left.attach=0,1, top.attach=2,3,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont4_table$attach(mrv$sidepcont4_color_button,left.attach=1,2, top.attach=2,3,
                              xoptions=c("fill"),yoptions="")
  
  ## pack table to main box
  mrv$sidepcont4_main$packStart(mrv$sidepcont4_table, expand=FALSE,fill=FALSE)
  
  mrv$sidepcont4$add(mrv$sidepcont4_main)
  
  ################ prediction region #######################
  mrv$sidepcont5 <- gtkFrame("Prediction Region")
  mrv$sidepcont5["border-width"]=10
  
  mrv$sidepcont5_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidepcont5_main["border-width"]=10
  
  mrv$sidepcont5_radio <- gtkHBox()
  mrv$sidepcont5_radiogp <- list()
  mrv$sidepcont5_radiogp$show <- gtkRadioButton(label = "Show")
  mrv$sidepcont5_radiogp$nshow <- gtkRadioButton(mrv$sidepcont5_radiogp, label = "Not Show")
  sapply(mrv$sidepcont5_radiogp, mrv$sidepcont5_radio$packStart)
  mrv$sidepcont5_radio[[1]]$setActive(TRUE) 
  sapply(mrv$sidepcont5_radiogp,'[',"active")
  mrv$sidepcont5_radio_align <- gtkAlignment(xalign = 0)
  mrv$sidepcont5_radio_align$add(mrv$sidepcont5_radio)
  
  mrv$prShow = T
  sapply(mrv$sidepcont5_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             if(button$getLabel()=="Show"){
               mrv$prShow = T
             }else{
               mrv$prShow = F
             }
           } 
         })
  
  # pack Show or not in main box
  mrv$sidepcont5_main$packStart(mrv$sidepcont5_radio_align, expand=FALSE,fill=FALSE)
  
  ####### other things in a table
  mrv$sidepcont5_table <- gtkTable(rows=3,columns=2,homogeneous=FALSE)
  mrv$sidepcont5_table$setColSpacings(5)
  mrv$sidepcont5_table$setRowSpacings(7)
  
  ### line type
  mrv$sidepcont5_type_label <- sidebarLabel("Line Type:",1)
  
  mrv$sidepcont5_type_combo <- gtkLineTypeComboBox(0)
  
  iter = mrv$sidepcont5_type_combo$GetActiveIter()
  mrv$pr.lty = as.numeric(mrv$sidepcont5_type_combo$GetModel()$GetValue(iter$iter,1)$value)
  gSignalConnect(mrv$sidepcont5_type_combo, "changed",function(cb){
    mrv$pr.lty <- getValueFromFigureCombo(cb)
  })
  
  ### lines width
  mrv$sidepcont5_width_label <- sidebarLabel("Line Width:",1)
  
  mrv$sidepcont5_width_combo <- gtkLineWidthComboBox(0)
  
  iter = mrv$sidepcont5_width_combo$GetActiveIter()
  mrv$pr.lwd = as.numeric(mrv$sidepcont5_width_combo$GetModel()$GetValue(iter$iter,1)$value)
  gSignalConnect(mrv$sidepcont5_width_combo, "changed",function(cb){
    mrv$pr.lwd <- getValueFromFigureCombo(cb)
  })
  ### lines color
  mrv$sidepcont5_color_label <- sidebarLabel("Line Color:",1)
  
  gdk_color <- gdkColorParse("darkgray")$color
  mrv$sidepcont5_color_button <- gtkColorButton(gdk_color)
  
  tempcol = mrv$sidepcont5_color_button$color
  tempalpha = mrv$sidepcont5_color_button$alpha
  mrv$pr.col = rgb(red=tempcol$red, green=tempcol$green, blue=tempcol$blue, alpha=tempalpha, maxColorValue=65535)
  gSignalConnect(mrv$sidepcont5_color_button, "color-set", function(button, ...) {
    mrv$pr.col = rgb(red=button$color$red, green=button$color$green, blue=button$color$blue, alpha=button$alpha, maxColorValue=65535)
  })
  ## pack table
  mrv$sidepcont5_table$attach(mrv$sidepcont5_type_label,left.attach=0,1, top.attach=0,1,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont5_table$attach(mrv$sidepcont5_type_combo,left.attach=1,2, top.attach=0,1,
                              xoptions="fill",yoptions="")
  mrv$sidepcont5_table$attach(mrv$sidepcont5_width_label,left.attach=0,1, top.attach=1,2,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont5_table$attach(mrv$sidepcont5_width_combo,left.attach=1,2, top.attach=1,2,
                              xoptions="fill",yoptions="")
  mrv$sidepcont5_table$attach(mrv$sidepcont5_color_label,left.attach=0,1, top.attach=2,3,
                              xoptions=c("expand","fill"),yoptions="")
  mrv$sidepcont5_table$attach(mrv$sidepcont5_color_button,left.attach=1,2, top.attach=2,3,
                              xoptions=c("fill"),yoptions="")
  
  ## pack table to main box
  mrv$sidepcont5_main$packStart(mrv$sidepcont5_table, expand=FALSE,fill=FALSE)
  
  mrv$sidepcont5$add(mrv$sidepcont5_main)
  ################# 
  mrv$sideplot$packStart(mrv$sidepcont1,expand=FALSE,fill=FALSE)
  mrv$sideplot$packStart(mrv$sidepcont2,expand=FALSE,fill=FALSE)
  mrv$sideplot$packStart(mrv$sidepsrocline,expand=FALSE,fill=FALSE)
  mrv$sideplot$packStart(mrv$sidepcont4,expand=FALSE,fill=FALSE)
  mrv$sideplot$packStart(mrv$sidepcont5,expand=FALSE,fill=FALSE)
  
  ################################################
  ###### sidebar  forest 
  ################################################
  mrv$sideforest <- gtkVBox(homogeneous=FALSE)
  mrv$sideforest_title <- gtkImage(filename=mrv$sidebar_icon[[5]])    #### maybe change move file to folder
  mrv$sideforest$packStart(mrv$sideforest_title,expand=FALSE,fill=FALSE)
  
  mrv$sideforest_window <- sidebarScrolledWindow(mrv$sideforest)
  
  ################ choose type ###########################
  mrv$sidefocont1 <- gtkFrame("Choose Accuracy Type")
  mrv$sidefocont1["border-width"]=10
  
  mrv$sidefocont1_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidefocont1_main["border-width"]=10
  
  mrv$sidefocont1_l_radio <- gtkVBox(homogeneous = TRUE)
  mrv$sidefocont1_r_radio <- gtkVBox(homogeneous = TRUE)
  mrv$sidefocont1_radiogp <- list()
  mrv$sidefocont1_radiogp$sens <- gtkRadioButton(label = "Sens")
  mrv$sidefocont1_radiogp$spec <- gtkRadioButton(mrv$sidefocont1_radiogp, label = "Spec")
  mrv$sidefocont1_radiogp$fpr <- gtkRadioButton(mrv$sidefocont1_radiogp, label = "FPR")
  mrv$sidefocont1_radiogp$fnr <- gtkRadioButton(mrv$sidefocont1_radiogp, label = "FNR")
  mrv$sidefocont1_radiogp$lrp <- gtkRadioButton(mrv$sidefocont1_radiogp, label = "LR+")
  mrv$sidefocont1_radiogp$lrn <- gtkRadioButton(mrv$sidefocont1_radiogp, label = "LR-")
  mrv$sidefocont1_radiogp$dor <- gtkRadioButton(mrv$sidefocont1_radiogp, label = "DOR")
  mrv$sidefocont1_radiogp$llrp <- gtkRadioButton(mrv$sidefocont1_radiogp, label = "LLR+")
  mrv$sidefocont1_radiogp$llrn <- gtkRadioButton(mrv$sidefocont1_radiogp, label = "LLR-")
  mrv$sidefocont1_radiogp$ldor <- gtkRadioButton(mrv$sidefocont1_radiogp, label = "LDOR")
  mrv$sidefocont1_l_radio$packStart(mrv$sidefocont1_radiogp$sens)
  mrv$sidefocont1_l_radio$packStart(mrv$sidefocont1_radiogp$fpr)
  mrv$sidefocont1_l_radio$packStart(mrv$sidefocont1_radiogp$lrp)
  mrv$sidefocont1_l_radio$packStart(mrv$sidefocont1_radiogp$llrp)
  mrv$sidefocont1_l_radio$packStart(mrv$sidefocont1_radiogp$dor)
  mrv$sidefocont1_r_radio$packStart(mrv$sidefocont1_radiogp$spec)
  mrv$sidefocont1_r_radio$packStart(mrv$sidefocont1_radiogp$fnr)
  mrv$sidefocont1_r_radio$packStart(mrv$sidefocont1_radiogp$lrn)
  mrv$sidefocont1_r_radio$packStart(mrv$sidefocont1_radiogp$llrn)
  mrv$sidefocont1_r_radio$packStart(mrv$sidefocont1_radiogp$ldor)
  mrv$sidefocont1_r_radio$packStart(gtkLabel(""))
  mrv$sidefocont1_l_radio[[1]]$setActive(TRUE) 
  sapply(mrv$sidefocont1_radiogp,'[',"active")
  
  mrv$accuracy = "sens"
  sapply(mrv$sidefocont1_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             mrv$accuracy=button$getLabel()
             if(mrv$accuracy=="LR+"){mrv$accuracy="LRpos"}
             if(mrv$accuracy=="LR-"){mrv$accuracy="LRneg"}
             if(mrv$accuracy=="LLR+"){mrv$accuracy="LLRpos"}
             if(mrv$accuracy=="LLR-"){mrv$accuracy="LLRneg"}
           } 
         })
  
  mrv$sidefocont1_table <- gtkTable(rows=1,columns=2,homogeneous=FALSE)
  mrv$sidefocont1_table$setColSpacings(5)
  mrv$sidefocont1_table$setRowSpacings(7)
  
  mrv$sidefocont1_table$attach(mrv$sidefocont1_l_radio,left.attach=0,1, top.attach=0,1,
                               xoptions=c("expand","fill"),yoptions="")
  mrv$sidefocont1_table$attach(mrv$sidefocont1_r_radio,left.attach=1,2, top.attach=0,1,
                               xoptions=c("expand","fill"),yoptions="")
  
  
  mrv$sidefocont1_main$packStart(mrv$sidefocont1_table, expand=FALSE,fill=FALSE)
  
  mrv$sidefocont1$add(mrv$sidefocont1_main)
  
  ################ Estimates type ###########################
  mrv$sidefocont2 <- gtkFrame("Estimates Type")
  mrv$sidefocont2["border-width"]=10
  
  mrv$sidefocont2_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidefocont2_main["border-width"]=10
  
  mrv$sidefocont2_radio <- gtkHBox(homogeneous = TRUE)
  mrv$sidefocont2_radiogp <- list()
  mrv$sidefocont2_radiogp$mean <- gtkRadioButton(label = "Mean")
  mrv$sidefocont2_radiogp$median <- gtkRadioButton(mrv$sidefocont1_radiogp, label = "Median")
  sapply(mrv$sidefocont2_radiogp, mrv$sidefocont2_radio$packStart)
  mrv$sidefocont2_radio[[1]]$setActive(TRUE) 
  sapply(mrv$sidefocont2_radiogp,'[',"active")
  mrv$sidefocont2_radio_align <- gtkAlignment(xalign = 0)
  mrv$sidefocont2_radio_align$add(mrv$sidefocont2_radio)
  
  mrv$forest.est.type = "mean"
  
  sapply(mrv$sidefocont2_radiogp, gSignalConnect, "toggled",
         f = function(button, ...){
           if(button['active']){
             mrv$forest.est.type=button$getLabel()
           } 
         })
  
  mrv$sidefocont2_main$packStart(mrv$sidefocont2_radio_align, expand=FALSE,fill=FALSE)
  mrv$sidefocont2$add(mrv$sidefocont2_main)
  
  ################ Show name ###########################
  mrv$sidefocont3 <- gtkFrame("Show Study Names")
  mrv$sidefocont3["border-width"]=10
  
  mrv$sidefocont3_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidefocont3_main["border-width"]=10
  
  mrv$sidefocont3_combo <- gtkComboBoxNewText()
  sapply(c("No", "Left", "Center", "Right"), mrv$sidefocont3_combo$appendText)
  mrv$sidefocont3_combo$setActive(3)
  
  mrv$forest.nameShow="Right"
  
  gSignalConnect(mrv$sidefocont3_combo, "changed", function(button, ...) {
    if(gtkComboBoxGetActiveText(button)=="No"){
      mrv$forest.nameShow=F
    }else{
      mrv$forest.nameShow=gtkComboBoxGetActiveText(button)
    }
  })
  
  mrv$sidefocont3_main$packStart(mrv$sidefocont3_combo, expand=FALSE,fill=FALSE)
  mrv$sidefocont3$add(mrv$sidefocont3_main)
  
  ################ Show data ###########################
  mrv$sidefocont4 <- gtkFrame("Show Original Data")
  mrv$sidefocont4["border-width"]=10
  
  mrv$sidefocont4_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidefocont4_main["border-width"]=10
  
  mrv$sidefocont4_combo <- gtkComboBoxNewText()
  sapply(c("No", "Left", "Center", "Right"), mrv$sidefocont4_combo$appendText)
  mrv$sidefocont4_combo$setActive(0)
  
  mrv$forest.dataShow=F
  
  gSignalConnect(mrv$sidefocont4_combo, "changed", function(button, ...) {
    if(gtkComboBoxGetActiveText(button)=="No"){
      mrv$forest.dataShow=F
    }else{
      mrv$forest.dataShow=gtkComboBoxGetActiveText(button)
    }
  })
  
  mrv$sidefocont4_main$packStart(mrv$sidefocont4_combo, expand=T,fill=T)
  mrv$sidefocont4$add(mrv$sidefocont4_main)
  
  ################ Show CI ###########################
  mrv$sidefocont5 <- gtkFrame("Show Confidence Interval")
  mrv$sidefocont5["border-width"]=10
  
  mrv$sidefocont5_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidefocont5_main["border-width"]=10
  
  mrv$sidefocont5_combo <- gtkComboBoxNewText()
  sapply(c("No", "Left", "Center", "Right"), mrv$sidefocont5_combo$appendText)
  mrv$sidefocont5_combo$setActive(1)
  
  mrv$forest.ciShow="Left"
  
  gSignalConnect(mrv$sidefocont5_combo, "changed", function(button, ...) {
    if(gtkComboBoxGetActiveText(button)=="No"){
      mrv$forest.ciShow=F
      if(mrv$notebook$getNPages()>2){
        Sys.sleep(.1) 
        asCairoDevice(mrv$forest_plot)
        forest(mrv$est, accuracy=mrv$accuracy, est.type=mrv$forest.est.type, pcex=mrv$forest.pcex, pch=mrv$forest.pch, pcol=mrv$forest.pcol,
               nameShow=mrv$forest.nameShow, dataShow=mrv$forest.dataShow, ciShow=mrv$forest.ciShow,
               shadecol="gray", arrowcol="black")
      }
    }else{
      mrv$forest.ciShow=gtkComboBoxGetActiveText(button)
      if(mrv$notebook$getNPages()>2){
        Sys.sleep(.1) 
        asCairoDevice(mrv$forest_plot)
        forest(mrv$est, accuracy=mrv$accuracy, est.type=mrv$forest.est.type, pcex=mrv$forest.pcex, pch=mrv$forest.pch, pcol=mrv$forest.pcol,
               nameShow=mrv$forest.nameShow, dataShow=mrv$forest.dataShow, ciShow=mrv$forest.ciShow,
               shadecol="gray", arrowcol="black")
      }
    }
  })
  
  mrv$sidefocont5_main$packStart(mrv$sidefocont5_combo, expand=FALSE,fill=FALSE)
  mrv$sidefocont5$add(mrv$sidefocont5_main)
  
  ################ Main plot ###########################
  mrv$sidefocont6 <- gtkFrame("Main Plot")
  mrv$sidefocont6["border-width"]=10
  
  mrv$sidefocont6_main <- gtkVBox(homogeneous=FALSE,spacing=10)
  mrv$sidefocont6_main["border-width"]=10
  
  mrv$sidefocont6_table <- gtkTable(rows=4,columns=2,homogeneous=FALSE)
  mrv$sidefocont6_table$setColSpacings(5)
  mrv$sidefocont6_table$setRowSpacings(7)
  
  mrv$sidefocont6_pch_label <- gtkLabel("Points Symbol:")
  mrv$sidefocont6_pch_label["xalign"]=1
  
  mrv$sidefocont6_pch_combo <- gtkPointsTypeComboBox(0)
  
  iter = mrv$sidefocont6_pch_combo$GetActiveIter()
  value = mrv$sidefocont6_pch_combo$GetModel()$GetValue(iter$iter,1)$value
  if(value %in% c("scaled","@","+","%","#","*","o","O")){
    mrv$forest.pch <- value
  }else{
    mrv$forest.pch <- as.numeric(value) 
  }
  gSignalConnect(mrv$sidefocont6_pch_combo, "changed",function(cb){
    mrv$forest.pch <- getValueFromFigureCombo(cb)
  })
  
  ### points size
  mrv$sidefocont6_pcex_label <- gtkLabel("Points Size:")
  mrv$sidefocont6_pcex_label["xalign"]=1
  
  mrv$sidefocont6_pcex_combo <- gtkPointsSizeComboBox(0)
  
  iter = mrv$sidefocont6_pcex_combo$GetActiveIter()
  mrv$forest.pcex = mrv$sidefocont6_pcex_combo$GetModel()$GetValue(iter$iter,1)$value
  gSignalConnect(mrv$sidefocont6_pcex_combo, "changed", function(cb){
    mrv$forest.pcex = getValueFromFigureCombo(cb)
    })
  ### points color
  mrv$sidefocont6_color_label <- gtkLabel("Points Color:")
  mrv$sidefocont6_color_label["xalign"]=1
  
  gdk_color <- gdkColorParse("black")$color
  mrv$sidefocont6_color_button <- gtkColorButton(gdk_color)
  
  
  mrv$forest.pcol = "black"
  gSignalConnect(mrv$sidefocont6_color_button, "color-set", function(button, ...) {
    mrv$forest.pcol = rgb(red=button$color$red, green=button$color$green, blue=button$color$blue, alpha=button$alpha, maxColorValue=65535)
  })
  
  ### points color
  mrv$sidefocont6_scolor_label <- gtkLabel("Shade Color:")
  mrv$sidefocont6_scolor_label["xalign"]=1
  
  gdk_scolor <- gdkColorParse("darkgray")$color
  mrv$sidefocont6_scolor_button <- gtkColorButton(gdk_scolor)
  
  
  mrv$forest.scol = "darkgray"
  gSignalConnect(mrv$sidefocont6_scolor_button, "color-set", function(button, ...) {
    mrv$forest.scol = rgb(red=button$color$red, green=button$color$green, blue=button$color$blue, alpha=button$alpha, maxColorValue=65535)
  })
  
  ## pack table
  mrv$sidefocont6_table$attach(mrv$sidefocont6_pch_label,left.attach=0,1, top.attach=0,1,
                               xoptions=c("expand","fill"),yoptions="")
  mrv$sidefocont6_table$attach(mrv$sidefocont6_pch_combo,left.attach=1,2, top.attach=0,1,
                               xoptions="fill",yoptions="")
  mrv$sidefocont6_table$attach(mrv$sidefocont6_pcex_label,left.attach=0,1, top.attach=1,2,
                               xoptions=c("expand","fill"),yoptions="")
  mrv$sidefocont6_table$attach(mrv$sidefocont6_pcex_combo,left.attach=1,2, top.attach=1,2,
                               xoptions="fill",yoptions="")
  mrv$sidefocont6_table$attach(mrv$sidefocont6_color_label,left.attach=0,1, top.attach=2,3,
                               xoptions=c("expand","fill"),yoptions="")
  mrv$sidefocont6_table$attach(mrv$sidefocont6_color_button,left.attach=1,2, top.attach=2,3,
                               xoptions=c("fill"),yoptions="")
  mrv$sidefocont6_table$attach(mrv$sidefocont6_scolor_label,left.attach=0,1, top.attach=3,4,
                               xoptions=c("expand","fill"),yoptions="")
  mrv$sidefocont6_table$attach(mrv$sidefocont6_scolor_button,left.attach=1,2, top.attach=3,4,
                               xoptions=c("fill"),yoptions="")
  
  ## pack table to main box
  mrv$sidefocont6_main$packStart(mrv$sidefocont6_table, expand=FALSE,fill=FALSE)
  
  mrv$sidefocont6$add(mrv$sidefocont6_main)
  
  ### graph size
  mrv$sidefocont7 <- gtkFrame("Main Text Size")
  mrv$sidefocont7["border-width"]=10
  
  mrv$sidefocont7_main <- gtkHBox(homogeneous=FALSE,spacing=10)
  mrv$sidefocont7_main["border-width"]=10
  
  mrv$sidefocont7_button <- gtkHScale(min=1, max=3, step=0.1)
  mrv$sidefocont7_button['draw-value'] <- FALSE
  #mrv$sidefocont7_button$setValuePos("right")
  adjustment <- mrv$sidefocont7_button$getAdjustment( )
  mrv$sidefocont7_spinbutton <- gtkSpinButton(adjustment = adjustment)
  
  mrv$sidefocont7_main$packStart(mrv$sidefocont7_button, expand=T, fill=T, padding=5)
  mrv$sidefocont7_main$packStart(mrv$sidefocont7_spinbutton, expand=F, fill=F, padding=5)
  
  mrv$sidefocont7$add(mrv$sidefocont7_main)
  
  mrv$forest.text.size = mrv$sidefocont7_spinbutton$value
  gSignalConnect(mrv$sidefocont7_spinbutton, "value-changed", function(button, ...) {
    mrv$forest.text.size = mrv$sidefocont7_spinbutton$value
  })
  
  ######    

  
  
  #########
  mrv$sideforest$packStart(mrv$sidefocont1,expand=FALSE,fill=FALSE)
  mrv$sideforest$packStart(mrv$sidefocont7,expand=FALSE,fill=FALSE)
  mrv$sideforest$packStart(mrv$sidefocont6,expand=FALSE,fill=FALSE)
  mrv$sideforest$packStart(mrv$sidefocont2,expand=FALSE,fill=FALSE)
  mrv$sideforest$packStart(mrv$sidefocont3,expand=FALSE,fill=FALSE)
  mrv$sideforest$packStart(mrv$sidefocont4,expand=FALSE,fill=FALSE)
  mrv$sideforest$packStart(mrv$sidefocont5,expand=FALSE,fill=FALSE)
  #mrv$sideforest$packStart(mrv$sidefocont8,expand=FALSE,fill=FALSE)
  ########################################## 
  #####    pack all pages to side notebook
  ##########################################
  prior_page_label = gtkLabel("Prior Control Panel")
  prior_page_label$setAngle(270)
  mrv$sidenote$insertPage(mrv$sideprior_window, prior_page_label)
  model_page_label = gtkLabel("Model Control Panel")
  model_page_label$setAngle(270)
  mrv$sidenote$insertPage(mrv$sidemodel_window, model_page_label)
  data_page_label = gtkLabel("Data Control Panel")
  data_page_label$setAngle(270)
  mrv$sidenote$insertPage(mrv$sidedata_window, data_page_label)
  sroc_page_label = gtkLabel("SROC Plot Control Panel")
  sroc_page_label$setAngle(270)
  mrv$sidenote$insertPage(mrv$sideplot_window, sroc_page_label)
  forest_page_label = gtkLabel("Forest Plot Control Panel")
  forest_page_label$setAngle(270)
  mrv$sidenote$insertPage(mrv$sideforest_window, forest_page_label)
  #funnel_page_label = gtkLabel("Funnel Plot Control Panel")
  #funnel_page_label$setAngle(270)
  #mrv$sidenote$insertPage(gtkLabel("funnel"), funnel_page_label)
  
  
  mrv$sideprior_page_num = mrv$sidenote$pageNum(mrv$sideprior_window)
  mrv$sidemodel_page_num = mrv$sidenote$pageNum(mrv$sidemodel_window)
  mrv$sidedata_page_num = mrv$sidenote$pageNum(mrv$sidedata_window)
  mrv$sideforest_page_num = mrv$sidenote$pageNum(mrv$sideforest_window)
  mrv$sidesroc_page_num = mrv$sidenote$pageNum(mrv$sideplot_window)
  
  #### pack all to sidebar
  mrv$sidebar$packStart(mrv$sidenote,expand=TRUE,fill=TRUE)
  
  ###########################################################
  ###########################################################
  #################       start construct notebook
  ###########################################################
  ###########################################################
  ### construct main notebook
  mrv$notebook <- gtkNotebook()
  mrv$notebook["show-border"] = FALSE
  mrv$notebook["homogeneous"] = FALSE
  mrv$notebook["tab-border"] = 0
  mrv$notebook$setTabPos("top")
  mrv$notebook$setScrollable(TRUE)
  
  ####################### Welcome page: 1st page
  mrv$welcome_page = gtkFrame()
  mrv$welcome_page["border-width"] = 0
  # Frame["shadow"] = "none"
  mrv$welcome_view <- gtkTextView()
  mrv$welcome_view["left-margin"] = 10
  mrv$welcome_swindow <- gtkScrolledWindow()
  mrv$welcome_swindow["shadow-type"] = "none"
  mrv$welcome_swindow$add(mrv$welcome_view)
  mrv$welcome_swindow$setPolicy("never", "automatic")
  
  mrv$welcome_buffer <- mrv$welcome_view$getBuffer()
  mrv$welcome_buffer$setText(paste("\n","\n",welcome_text, sep="",collapse = ""))
  ### buffer tags
  tag_bold <- mrv$welcome_buffer$createTag(tag.name="bold", weight=PangoWeight["bold"])
  tag_emph <- mrv$welcome_buffer$createTag(tag.name="emph", style=PangoStyle["italic"])
  tag_Large <- mrv$welcome_buffer$createTag(tag.name="Large", font=" 18")
  tag_large <- mrv$welcome_buffer$createTag(tag.name="large", font=" 15")
  tag_color_red <- mrv$welcome_buffer$createTag(tag.name="red", foreground="#FF0000")
  
  
  first_line_iter <- mrv$welcome_buffer$getIterAtLine(1)
  mrv$welcome_buffer$insertWithTags(first_line_iter$iter,"Welcome to the meta4diag World! \n",tag_bold,tag_Large,tag_color_red)
  
  fifth_line_iter <- mrv$welcome_buffer$getIterAtLine(12)
  mrv$welcome_buffer$insertWithTags(fifth_line_iter$iter,"User's Guide: \n",tag_bold,tag_large)
  
  
  mrv$welcome_page$add(mrv$welcome_swindow)
  welcome_page_event <- .pageEventsLabel("Welcome!")
  mrv$notebook$insertPage(mrv$welcome_page,welcome_page_event)
  
  welcome_page_num = mrv$notebook$pageNum(mrv$welcome_page)
  
  gSignalConnect(welcome_page_event, "button-press-event", function(...) {
    mrv$notebook$setCurrentPage(welcome_page_num)
    mrv$sidenote$setCurrentPage(mrv$sideprior_page_num)
    return(TRUE)
  })
  
  mrv$midpart <- gtkHBox()
  mrv$midpart$packStart(mrv$sidebar,expand=FALSE,fill=FALSE)
  mrv$midpart$packStart(mrv$notebook,expand=TRUE,fill=TRUE)
  
  ###########################################################
  #################       ALL TOGETHER
  ###########################################################
  mrv$vbox <- gtkVBox(homogeneous = FALSE, spacing = 0)
  mrv$vbox$packStart(mrv$menubar, expand = FALSE, fill = FALSE, padding = 0)
  mrv$vbox$packStart(mrv$toolbar, FALSE, FALSE, 0)
  mrv$vbox$packStart(mrv$midpart, TRUE,TRUE,0)
  mrv$vbox$packStart(mrv$statusbar, FALSE, FALSE, 0)
  mrv$main_window$add(mrv$vbox)
  mrv$main_window$show()
  
  #########################################################
  ########### Install package
  #########################################################
  if(requireNamespace("INLA", quietly = TRUE)){
    if("INLA" %in% rownames(installed.packages()) == FALSE){
      INLA_dialog <- gtkMessageDialog(NULL,"destroy-with-parent","question","yes-no",paste("R Package \"INLA\" is not installed.","\n", "We suggest to install it.","\n",
                                                                                           "Do you want to install INLA?","\n",
                                                                                           "After installation, please load INLA!",
                                                                                           "\n","Thank you!",sep=""))
      INLA_dialog["title"] <- "INLA installation..."
      INLA_choices <- c("Stable version", "Testing version")
      INLA_radio_buttons <- NULL
      INLA_vbox <- gtkVBox(FALSE,0)
      for(choice in INLA_choices){
        INLA_choices_button <- gtkRadioButton(INLA_radio_buttons, choice)
        INLA_vbox$add(INLA_choices_button)
        INLA_radio_buttons <- c(INLA_radio_buttons, INLA_choices_button)
      }
      INLA_frame <- gtkFrame("Install amazing INLA package")
      INLA_frame$add(INLA_vbox)
      INLA_dialog[["vbox"]]$add(INLA_frame)
      INLA_radio_buttons[[1]]$setActive(TRUE) 
      # sapply(INLA_radio_buttons,'[',"active")
      sapply(INLA_radio_buttons, gSignalConnect, "toggled",
             f = function(button, ...){
               if(button['active']){
                 INLA_version = button$getLabel()
                 # print(INLA_version)
                 if(INLA_version=="Stable version"){
                   mrv$repos="http://www.math.ntnu.no/inla/R/stable"
                 } else{
                   mrv$repos="http://www.math.ntnu.no/inla/R/testing"
                 }
               } 
             })
      gSignalConnect(INLA_dialog,"response",f=function(dialog,response,user.data){
        if(response == GtkResponseType["no"]){
          INLA_no_dialog <- gtkMessageDialog(INLA_dialog,"destroy-with-parent","error","close","Wrong Choice! INLA must be installed!!!")
          INLA_no_dialog$run()
          INLA_no_dialog$destroy()
        }else{
          install.packages("INLA", repos=mrv$repos)
          INLA_dialog$destroy()
        }
      })
    }
    if (!(sum(search()=="package:INLA"))==1){
      INLA_dialog <- gtkMessageDialog(NULL,"destroy-with-parent","warning","ok",paste("R Package \"INLA\" is not loaded.","\n",
                                                                                           "You have to load it to make sure meta4diag works.","\n",
                                                                                           "Thank you!",sep=""))
      INLA_dialog["title"] <- "INLA loading..."
      if (INLA_dialog$run()==GtkResponseType["ok"]){
        INLA_dialog$destroy()
      }
    }
  }else{
    INLA_dialog <- gtkMessageDialog(NULL,"destroy-with-parent","destroy-with-parent","warning","ok",paste("R Package \"INLA\" is not installed.","\n", "We suggest to install it.","\n",
                                                                                         "After installation, please load INLA","\n",sep=""))
    if (INLA_dialog$run()==GtkResponseType["ok"]){
      INLA_dialog$destroy()
    }
  }
}