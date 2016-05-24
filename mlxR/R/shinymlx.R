#' 
#' Automatic code generation for Shiny applications
#' 
#' Creates a Shiny application for longitudinal data model  
#' 
#' shinymlx automatically generates files ui.R and server.R required 
#' for a Shiny application.
#' 
#' Elements of \code{parameters} and \code{treatment} can be either scalars or lists.
#' A scalar automatically generates a slider with default minimum and maximum values and default step. 
#' A list may contain the type of widget ("slider", "select", "numeric") and the settings defining the widget: 
#' (value, min, max, step) for \code{slider}, (selected, choices) for \code{select} and
#' value for \code{numeric}.
#'
#' See http://simulx.webpopix.org/mlxr/shinymlx/ for more details.      
#' @param model a \code{Mlxtran} or \code{PharmML} model used for the simulation
#' @param output a list - or a list of lists - with fields: 
#' \itemize{
#'   \item \code{name}: a vector of output names
#'   \item \code{time}: a vector of times, or a vector (min, max, step)
#' }
#' @param parameter a vector, or a list of shiny widgets
#' @param treatment a list with fields 
#' \itemize{
#'   \item \code{tfd} : first time of dose,
#'   \item \code{amount} : amount,
#'   \item \code{nd} : number of doses,
#'   \item \code{ii} : interdose interval,
#'   \item \code{type} : the type of input,
#' }
#' Input argument of Simulx can also be used,  i.e. a list with fields \code{time}, 
#' \code{amount}, \code{rate}, \code{tinf}, \code{type}, \code{target}.
#' 
#' @param data a datafile to display with the plot
#' @param title the title of the application
#' @param appname the name of the application (and possibly its path)
#' @param style the style of the Shiny app
#' \itemize{
#'   \item \code{"basic"}: basic Shiny app with a single side bar (default)
#'   \item \code{"navbar1"}: navigation bar and tabPanels (including outputs)
#'   \item \code{"navbar2"}: navigation bar and tabPanels (outputs separated)
#'   \item \code{"dashboard1"} : Shiny dashboard,
#' }
#' @param settings a list of settings
#' \itemize{
#'   \item \code{"tabstyle"} : look of the tabs c("tabs","pills"),
#'   \item \code{"select.x"} : display the list of variables available for the x-axis c(TRUE,FALSE),
#'   \item \code{"select.y"} : display the list of variables available for the y-axis c(TRUE,FALSE),
#'   \item \code{"select.log"} : log scale option c(TRUE,FALSE),
#'   \item \code{"select.ref"} : reference curves option c(TRUE,FALSE)
#' }   
#' @return A Shiny app with files ui.R, server.R and model.txt
#' @examples
#' \dontrun{
#' library(mlxR)
#' PKPDmodel <- inlineModel("
#' [LONGITUDINAL] 
#' input={ka,V,Cl,Imax,IC50,S0,kout}
#' EQUATION:
#' C     = pkmodel(ka, V, Cl)
#' E_0   = S0 
#' ddt_E = kout*((1-Imax*C/(C+IC50))*S0- E)  
#' ")
#' 
#' p1  <- c(ka=0.5, V=10, Cl=1)
#' p2  <- c(Imax=0.5, IC50=0.03, S0=100, kout=0.1)
#' adm <- list(tfd=5, nd=15, ii=12, amount=1)
#' f1  <- list(name = 'C', time = seq(0, 250, by=1))
#' f2  <- list(name = 'E', time = seq(0, 250, by=1))
#' f   <- list(f1, f2)
#' 
#' shinymlx(model=PKPDmodel, treatment=adm, parameter=list(p1,p2), output=f,
#'          style="dashboard1")
#' 
#' #------------------------------------------------------------------------
#' p1 <- list(
#'   ka    = list(widget="slider", value=0.5, min=0.1, max=2,  step=0.1),
#'   V     = list(widget="slider", value=10,  min=2,   max=20, step=2),
#'   Cl    = list(widget="slider", value=1,   min=0.1, max=2,  step=0.1)
#' )  
#' adm <- list(
#'   tfd    = list(widget="slider", value=5,  min=0, max=100, step=5),
#'   nd     = list(widget="numeric", value=15),
#'   ii     = list(widget="select", selected=12, choices=c(3,6,12,18,24)),
#'   amount = list(widget="slider", value=40, min=0, max=50, step=5)
#' )
#' s <- list(select.x=FALSE, select.y=FALSE)
#' shinymlx(model=PKPDmodel, treatment=adm, parameter=list(p1,p2), output=f, 
#'          style="navbar1", settings=s)
#' }
#' @importFrom utils read.csv
#' @export         
shinymlx <- function(model,parameter=NULL,output=NULL,treatment=NULL,
                     data=NULL,appname="shinymlxApp",style="basic",
                     settings=NULL,title=" ")
{
  
  select=list()
  select$x <- set.settings(settings$select.x,style)
  select$y <- set.settings(settings$select.y,style)
  select$log <- set.settings(settings$select.log,style)
  select$ref <- set.settings(settings$select.ref,style)
  
  if (is.null(settings$tabstyle)) 
    tabstyle <- "tabs"
  else  
    tabstyle <- settings$tabstyle
  
  s2r <- ""
  s2f <- ""
  unlink(file.path(mainDir=appname), recursive = TRUE, force = TRUE)
  dir.create(file.path(mainDir=appname), showWarnings = FALSE)
  file.copy(model,file.path(appname,"model.txt"),overwrite=TRUE)
  
  output <- foutput(output)  
  s1 <- out2str(output)
  if (select$x==TRUE)
    select$x <- testout(output)
  
  i.data=0
  if (!is.null(data)){
    i.data=1
    file.copy(data,file.path(appname,basename(data)),overwrite=TRUE)
    data.txt <- paste0("datax <- read.csv(file='",basename(data),"', header=TRUE, sep='\\t', quote='\')")
    s1 <- paste0(s1,'\n',data.txt)
    data <- read.csv(file=data, header=TRUE, sep="\t", quote="\"")
  }
  parameter <- fparameter(parameter)
  ptxt <- param2str(parameter)
  s2r <- paste0(s2r, ptxt[[1]])
  s2f <- paste0(s2f, ptxt[[2]])
  if (!is.null(treatment)){
    treatment <- ftreatment(treatment)
    ttxt <- adm2str(treatment)
    s2r <- paste0(s2r, ttxt[[1]])
    s2f <- paste0(s2f, ttxt[[2]])
    s3 <- (
      "    r <- simulx( model     = 'model.txt',
                 treatment = adm,
                 parameter = p,
                 output    = f)
")    
  }else{
    s3 <- (
      "    r <- simulx( model     = 'model.txt',
                 parameter = p,
                 output    = f)
")    
    treatment <- NULL
  }
  
  
  if (i.data==1)
    s4 <- paste0("    pl <- pl + geom_point(data=datax,aes(x=",names(data)[1],",y=",names(data)[2],"),size=3.5,color='#6666CC')")
  else
    s4 <- ""
  
  sui <- list2ui(parameter,treatment,output,select,style,tabstyle)
  if ((sui$i.output == FALSE) & (identical(style,"navbar2")))
    style <- "navbar1"
  ui.txt <- uiTemplate(title,sui$content,sui$output,style,tabstyle)
  write(ui.txt,file.path(appname,"ui.R"))
  server.txt <- serverTemplate(c(s1,s2r,s2f,s3,s4),select,sui$i.output,select$y)
  write(server.txt,file.path(appname,"server.R"))
  
  stools.txt <- stoolsmlx(select$x)
  write(stools.txt,file.path(appname,"shinymlxTools.R"))
  return(appname)
}

testout <- function(out){
  n <- length(out)
  select.x <- TRUE
  if (n>1){
    t <- out[[1]]$time
    for (k in (2:n)){
      if (!identical(out[[k]]$time, t))
        select.x <- FALSE
    }
  }else if (length(out[[1]]$name)==1)
    select.x <- FALSE
  
  return(select.x)
}

foutput <- function(out){
  if (!is.null(names(out)))
    out <- list(out)
  for (k in seq(1,length(out))){
    outk <- out[[k]]
    out.time <- outk$time
    if (length(out.time)>3){
      outk$time <- c(min(out.time), max(out.time), diff(range(out.time))/(length(out.time)-1))
    }else if (length(out.time)==2){
      outk$time <- c(min(out.time), max(out.time), 1)
    }
    if (is.null(names(outk)))
      names(outk <- c("min", "max", "by"))  
    out[[k]] <- outk
  } 
  return(out)
}


ftreatment <- function(trt){
  if (!is.null(names(trt)))
    trt <- list(trt)
  for (k in seq(1,length(trt))){
    trtk <- trt[[k]]
    nk <- names(trtk)
    for (j in (1: length(trtk))){
      pj <- trtk[[j]]
      if (!is.list(pj)){
        if (length(pj)==1){
          if (nk[j]=="type")
            trtk[[j]] <- list(value=pj, widget="numeric")
          else if (nk[j]=="target")
            trtk[[j]] <- list(value=pj, widget="none")
          else
            trtk[[j]] <- list(value=pj, min=pj/2, max=pj*2, step=pj*0.1, widget="slider")
        }else if (("ii" %in% nk) & (length(pj)==4)){
          trtk[[j]] <- list(value=pj[1], min=pj[2], max=pj[3], step=pj[4], widget="slider")        
        }else{
          trtk[[j]] <- list(value=pj, widget="none")        
        }
      }
    }
    trt[[k]] <- trtk
  } 
  return(trt)
}

fparameter <- function(param) 
{
  if (!is.null(names(param))){  
    param=list(param) 
  } 
  for (k in seq(1,length(param))){
    paramk <- param[[k]]
    if (is.list(paramk)){
      if (!is.null(paramk$value)){
        pk <- paramk$value
        names(pk) <- paramk$name
        paramk <- pk
      }
    }
    if (!is.list(paramk)){
      paramk <- as.list(paramk)
      for (j in (1:length(paramk)))
        paramk[[j]] <- c(1,0.5,2,0.1)*paramk[[j]]
    }
    for (j in (1: length(paramk))){
      pj <- paramk[[j]]
      if (!is.list(pj)){
        if (length(pj)==1){
          paramk[[j]] <- list(value=pj, widget="none")
        }else if (length(pj)==4){
          paramk[[j]] <- list(value=pj[1], min=pj[2], max=pj[3], step=pj[4], widget="slider")
        }else{
          stop("Error:  length of a parameter defined as a vector should be 4 (slider), or 1 (constant)")
        }        
      }
    }
    param[[k]] <- paramk
  } 
  return(param)
}


#-------------------------------------------------------

adm2str <- function(y, sp.adm="    "){
  K <- length(y)
  for (k in (1:K)){
    x <- y[[k]]
    length.x <- length(x)
    names.x  <- names(x)
    
    if ("time" %in% names.x){
      time.x <- vect2str(x[["time"]]$value)
      txt.tdose1 <- paste0(sp.adm,"t.dose <- ",time.x)
      txt.tdose2 <- txt.tdose1
    }else{
      if ("tfd" %in% names.x){
        if (x[["tfd"]]$widget == "none"){
          stfd <- x[["tfd"]]$value
        }else{
          stfd <- paste0("input$tfd",k)
        }
      }else{
        stfd <- "0"
      }
      txt.tdose1 <- paste0(sp.adm,"t1 <- ",stfd)
      txt.tdose2 <- paste0(sp.adm,"t1 <- isolate(",stfd,")")
      if ("nd" %in% names.x){
        sii <- widgtxt(x,"ii",k)
        snd <- widgtxt(x,"nd",k)
        txt.tdose1 <- paste0(txt.tdose1,"\n",sp.adm,"t2 <- ",sii,"*(",snd,"-1)+t1")
        txt.tdose1 <- paste0(txt.tdose1,"\n",sp.adm,"t.dose <- seq(t1,t2,by=",sii,")")
        sii <- paste0("isolate(",sii,")")
        snd <- paste0("isolate(",snd,")")
        txt.tdose2 <- paste0(txt.tdose2,"\n",sp.adm,"t2 <- ",sii,"*(",snd,"-1)+t1")
        txt.tdose2 <- paste0(txt.tdose2,"\n",sp.adm,"t.dose <- seq(t1,t2,by=",sii,")")
      }else{
        txt.tdose1 <- paste0(txt.tdose1,"\n",sp.adm,"t.dose <- t1")
        txt.tdose2 <- paste0(txt.tdose2,"\n",sp.adm,"t.dose <- t1")
      }
      
    }
    if (K==1){
      sadmk <- "adm"
    }else{    
      sadmk <- paste0("adm",as.character(k))
    }
    txt1 <- paste0(txt.tdose1,"\n",sp.adm,sadmk," <- list(time=t.dose, ")
    txt2 <- paste0(txt.tdose2,"\n",sp.adm,sadmk," <- list(time=t.dose, ")
    
    for (m in (1:length.x)){
      xk <- x[[m]]
      if (!(names(x[m]) %in% c("tfd","nd","ii","time"))){
        if (xk$widget=="none"){
          txtk1 <- vect2str(xk$value)
          txtk2 <- txtk1
        }else{
          txtk1 <- paste0("input$",names.x[m],k)
          txtk2 <- paste0("isolate(input$",names.x[m],k,")")
        }
        txt1 <- paste0(txt1,names.x[m],"=",txtk1,",")
        txt2 <- paste0(txt2,names.x[m],"=",txtk2,",")
      }
    }
    substr(txt1,nchar(txt1),nchar(txt1))<- ')'
    substr(txt2,nchar(txt2),nchar(txt2))<- ')'
    if (k==1){
      atxt1 <- txt1
      atxt2 <- txt2
    }else{
      atxt1 <- paste0(atxt1,"\n",txt1)
      atxt2 <- paste0(atxt2,"\n",txt2)
    }
  }
  if (K>1){
    txt <- "adm <- list(adm1"
    for (k in (2:K))
      txt <- paste0(txt,", adm",as.character(k))
    txt <- paste0(txt,")")
    atxt1 <- paste0(atxt1,"\n",sp.adm,txt)
    atxt2 <- paste0(atxt2,"\n",sp.adm,txt)
  }
  adm.txt <- list(atxt1, atxt2)
  return(adm.txt)
}

#-------------------------------
widgtxt <- function(x,txt,k){
  if (x[[txt]]$widget != "none"){
    if (x[[txt]]$widget != "select")
      s <- paste0("input$",txt,k)
    else
      s <- paste0("as.numeric(input$",txt,k,")")
  }else{
    s <- vect2str(x[[txt]]$value)
  }
  return(s)
}

#------------------------------------------------------
out2str <- function(out){
  K <- length(out)
  for (k in (1:K)){
    outk <- out[[k]]
    txt.nameo <- vect2str(outk$name)
    t.out <- outk$time
    if (K==1){
      soutk <- "f"
    }else{    
      soutk <- paste0("out",as.character(k))
    }
    txt <- paste0(soutk," <- list(time=t.dose, ")
    if (k==1){
      out.txt <- paste0(soutk," <- list(name=",txt.nameo,
                        ", time=seq(",t.out[1],",",t.out[2],",by=",t.out[3],"))")
    }else{    
      out.txt <- paste0(out.txt,"\n",soutk," <- list(name=",txt.nameo,
                        ", time=seq(",t.out[1],",",t.out[2],",by=",t.out[3],"))")
    }
  }
  if (K>1){
    txt <- "f  <- list(out1"
    for (k in (2:K))
      txt <- paste0(txt,", out",as.character(k))
    txt <- paste0(txt,")")
    out.txt <- paste0(out.txt,"\n",txt)
  }else{
    out.txt <- paste0(out.txt,"\nf <- list(f)")
  }
  return(out.txt)
}

#-------------------------------------------------------
vect2str <- function(x){
  if (is.character(x[1])){
    txt <- paste0("'",x[1],"'")
  }else{
    txt <- as.character(x[1])
  }
  if (length(x)>1){
    for (j in (2:length(x)))
      if (is.character(x[j])){
        txt <- paste0(txt,",'",x[j],"'")
      }else{
        txt <- paste0(txt,",",as.character(x[j]))
      }
    txt <- paste0("c(",txt,")")
  }
  return(txt)
}


#-------------------------------------------------------

param2str <- function(p, sp.param="    "){
  txt.input <- NULL
  txt.namep <- NULL
  for (j in (1:length(p))){
    pj <- p[[j]]
    names.p <- names(pj)
    np <- length(names.p) 
    for (k in 1:np){
      pk <- pj[[k]]
      nk <- names.p[k]
      txt.namep <- paste0(txt.namep,"'",nk,"', ")
      if (pk$widget != "none"){
        if (pk$widget == "select")
          txt.input <- paste0(txt.input,'as.numeric(input$',nk,'), ')
        else
          txt.input <- paste0(txt.input,'input$',nk,', ')
      }else{
        txt.input <- paste0(txt.input,pk$value,", ")      
      }
    }
  }
  txt.input <- substr(txt.input,1,(nchar(txt.input)-2))
  txt.namep <- substr(txt.namep,1,(nchar(txt.namep)-2))
  txtr <- paste0(sp.param,"p <- list(name  = c(",txt.namep,"),\n",
                 sp.param,"          value = c(",txt.input,"))\n")
  txtf <- paste0(sp.param,"p <- list(name  = c(",txt.namep,"),\n",
                 sp.param,"          value = isolate(c(",txt.input,")))\n")
  txt <- list(txtr, txtf)
  return(txt) 
}

#-------------------------------------------------------

uiTemplate <- function(title,s.content,s.output,style,tabstyle)
{
  if (identical(style,"basic")){
    ui.out <- paste0('
ui <- shinyUI(fluidPage(
fluidRow(
column(3,br(),
',s.content,'
),
column(9,br(),            
plotOutput("plot")
)
)
))
')    
    
  }else if (identical(style,"navbar1")){
    ui.out <- paste0('
ui <- shinyUI(fluidPage(
navbarPage("',title,'",
tabPanel("Plot",
fluidRow(
column(3,br(),
',s.content,'
),
column(9,br(),            
plotOutput("plot")
)
)
),
tabPanel("Mlxtran", pre(includeText("model.txt"))),
tabPanel("ui.R", pre(includeText("ui.R"))),
tabPanel("server.R", pre(includeText("server.R")))
)
))
')    
    
  }else if (identical(style,"navbar2")){
    ui.out <- paste0('
ui <- shinyUI(fluidPage(
navbarPage("',title,'",
tabPanel("Plot",
fluidRow(
column(2,
',s.content,'
),
column(2,
br(),
',s.output,'
),
column(8,br(),            
plotOutput("plot")
)
)
),
tabPanel("Mlxtran", pre(includeText("model.txt"))),
tabPanel("ui.R", pre(includeText("ui.R"))),
tabPanel("server.R", pre(includeText("server.R")))
)
))
')
    
  }else{
    if (!is.null(s.output)){
      sout <- paste0('
hr(),
fluidRow(               
column(1),
column(9,\n',
                     s.output,'
)
),   
')   
    }else{
      sout <- NULL
    }
    ui.out <- paste0('
library(shinydashboard)

sidebar <- dashboardSidebar(
hr(),
sidebarMenu(id="tabs",
menuItem("Plot", tabName = "plot", icon = icon("line-chart"), selected=TRUE),
menuItem("Codes",  icon = icon("file-text-o"),
menuSubItem("Mlxtran", tabName = "mlxtran", icon = icon("angle-right")),
menuSubItem("ui.R", tabName = "ui", icon = icon("angle-right")),
menuSubItem("server.R", tabName = "server", icon = icon("angle-right")),
menuSubItem("shinymlxTools.R", tabName = "tools", icon = icon("angle-right"))
)
),
',sout,'
hr()
)

body <- dashboardBody(
tabItems(
tabItem(tabName = "plot",
fluidRow(               
box(width = 3,   status = "primary",
',s.content,'          
),
box(width = 9,   height=530, status = "primary",
plotOutput("plot")
)
)
),
tabItem(tabName = "mlxtran",
box( width = NULL, status = "primary", solidHeader = TRUE, title="model.txt",                
pre(includeText("model.txt"))
)
),
tabItem(tabName = "ui",
box( width = NULL, status = "primary", solidHeader = TRUE, title="ui.R",
pre(includeText("ui.R"))
)
),
tabItem(tabName = "server",
box( width = NULL, status = "primary", solidHeader = TRUE, title="server.R",
pre(includeText("server.R"))
)
),
tabItem(tabName = "tools",
box( width = NULL, status = "primary", solidHeader = TRUE, title="server.R",
pre(includeText("shinymlxTools.R"))
)
)
#       tabItem(tabName = "about", includeMarkdown("../../about/about.Rmd"))
)
)

dashboardPage(
dashboardHeader(title = "',title,'"),
sidebar,
body
)

')
  }
  
  return(ui.out)
}

#-------------------------------------------------------

serverTemplate <- function(s, select, i.output, select.y)
{
  if (select$x==TRUE){
    spl <- paste0(
      "   pj <- paste0('pl <- pl + geom_path(data=res, aes(x=',xj,',y=',name.fj[k],',colour=",'"',"',info[[j]]$colour[k],'",'"',"),size=0.75)')
    eval(parse(text=pj))") 
    if (select$ref==TRUE){
      srf <- paste0(
        "   if (input$boxref==TRUE){
      pj <- paste0('pl <- pl + geom_path(data=ref, aes(x=',xj,',y=',name.fj[k],'),colour=",'"grey",',"size=0.75)')
      eval(parse(text=pj))
    }") 
    }else{
      srf=""
    }
  }else{
    spl <- paste0(
      "       pj <- paste0('pl <- pl + geom_path(data=res[[j]], aes(x=time,y=',name.fj[k],',colour=",'"',"',info[[j]]$colour[k],'",'"',"),size=0.75)')
       eval(parse(text=pj))") 
    if (select$ref==TRUE){
      srf <- paste0(
        "   if (input$boxref==TRUE){
      pj <- paste0('pl <- pl + geom_path(data=ref[[j]], aes(x=time,y=',name.fj[k],'),colour=",'"grey",',"size=0.75)')
      eval(parse(text=pj))
    }") 
    }else{
      srf=""
    }
  }
  srfpl <- paste0(srf,'\n',spl) 
  
  if (select$log==TRUE){
    slog <-         
      "    if (input$ilog==TRUE)
       pl=pl + scale_y_log10()"
    
  }else{
    slog=""
  }
  
  if (select$ref==TRUE){
    sref1 <- paste0(
      "  ref <- reactive({
    input$butref
",s[3],"
",s[4],"
    ref <- merge_res(r,f)
    return(ref)
  })\n")
    sref2 <- "   ref=ref()"
  }else{
    sref1 <- ""
    sref2 <- ""
  }
  
  if (i.output==TRUE){
    if (select.y)
      sply <- 'eval(parse(text=paste0("inputyj=input$out",j)))'
    else
      sply <-  'inputyj <- name.fj'
    splot <- paste0('\n   ',sply,'  
  i.plot=FALSE
  if (!is.null(inputyj)){
    ij <- which(name.fj %in% inputyj)
    if (length(ij>0)){
      eval(parse(text=paste0("inputxj=input$x",j)))
      if (!is.null(inputxj))
        xj <- inputxj
    }
    i.plot=TRUE
  }
  else if (is.null(inputyj) & length(f)==1){
    ij=1
    i.plot=TRUE
  }
  if (i.plot){
    pl <- ggplotmlx()
    nfj <- length(name.fj)
    for (k in (1:nfj)){
      if (k %in% ij){
        ',srfpl,'
      }  
    }
    pl <- pl + scale_colour_manual(values=info[[j]]$values, labels=info[[j]]$labels)
    print(pl)
    if (length(ij)>1){
      if (!is.null(input$legend) && input$legend==FALSE)
        pl <- pl + theme(legend.position="none")
      else
        pl <- pl + guides(colour=guide_legend(title=NULL)) + theme(legend.position=c(.9, .8))
      pl <- pl + ylab("")
    }else{
      pl <- pl + theme(legend.position="none")
    }
  ',s[5],'
  ',slog,'
    eval(parse(text=paste0("pl",j," <- pl")))
    gr.txt <- paste0(gr.txt,"pl",j,",")
  }')
  }else{
    splot <- paste0('
      pl <- ggplotmlx()
      nfj <- length(name.fj)
      for (k in (1:nfj)){',
                    srfpl,'
      }
      pl <- pl + scale_colour_manual(values=info[[j]]$values, labels=info[[j]]$labels)
      print(pl)
      if (length(name.fj)>1)
        pl <- pl + guides(colour=guide_legend(title=NULL)) + theme(legend.position=c(.9, .8))
      else
        pl <- pl + theme(legend.position="none")
      pl <- pl + ylab("")
      ',s[5],'
      ',slog,'
      eval(parse(text=paste0("pl",j," <- pl")))
      gr.txt <- paste0(gr.txt,"pl",j,",")')  
  }
  
  server.out <- paste0(
    'library("mlxR")
library("reshape")
library("gridExtra")
source("shinymlxTools.R")

',s[1],'
nf <- length(f)
info <- info_res(f)

server <- function(input, output) {
',sref1,'  
  res <- reactive({
',s[2],'                     
',s[4],'                     
    res <- merge_res(r,f)
    return(res)
  })  
  
  output$plot <- renderPlot({
    res=res()
',sref2,'    
    gr.txt <- "grid.arrange("
    for (j in (1:length(f))){
      xj <- "time"
      fj <- f[[j]]
      name.fj <- fj$name
',splot,'
    }
    gr.txt <- paste0(gr.txt,"ncol=1)")
    eval(parse(text=gr.txt))
  }, height = 500)
}
')
  return(server.out)
  
}

#-------------------------------------------------------
list2ui <- function(y.p=NULL,y.a=NULL, y.o=NULL, select, style, tabstyle){
  x.p <- list()
  x.a <- list()
  
  if (!is.null(y.p)){
    for (k in (1:length(y.p))){
      yk <- y.p[[k]]
      ik <- 0
      for (j in (1:length(yk))){
        if (yk[[j]]$widget != "none")
          ik <- 1
      }
      if (ik == 1)
        x.p <- c(x.p, list(yk))   
    }
  }
  if (!is.null(y.a)){
    for (k in (1:length(y.a))){
      yk <- y.a[[k]]
      ik <- 0
      for (j in (1:length(yk))){
        if (yk[[j]]$widget != "none")
          ik <- 1
      }
      if (ik == 1)
        x.a <- c(x.a, list(yk))    
    }
  }
  
  name.out <- NULL
  i.leg <- FALSE
  for (k in (1:length(y.o))){
    name.out <- c(name.out, y.o[[k]]$name)
    if (length(y.o[[k]]$name)>1 & select$y==TRUE)
      i.leg <- TRUE
  }
  
  n.p <- length(x.p)
  n.a <- length(x.a)
  n.o <- length(name.out)
  K <- n.p + n.a
  x.type <- c(rep("param",n.p),rep("adm",n.a))
  
  k.a <- k.p <- 0
  tab.content <- list()
  tab.label <- vector(length=K)
  for (k in (1:K)){
    txt <- NULL
    if (x.type[k]=="param"){
      k.p <- k.p+1
      xk <- x.p[[k.p]]
      label <- "param"
      if (n.p>1)
        label <- paste0(label,k.p)
    }else{
      k.a <- k.a+1
      xk <- x.a[[k.a]]
      label <- "adm"
      if (n.a>1)
        label <- paste0(label,k.a)
    }
    for (j in (1:length(xk))){
      xkj <- xk[[j]] 
      nxkj <- names(xk[j])
      if (x.type[k]=="param")
        id <- nxkj
      else
        id <- paste0(nxkj,k.a)
      if (xkj$widget!="none"){
        uikj <- widgetui(xkj,id, nxkj)
        txt <- paste0(txt,uikj,",\n")
      }
    }
    if (!identical(style,"basic"))
      txt <- paste0(txt,'br()')
    tab.content[[k]] <- txt
    tab.label[k] <- label
  }
  out.txt <- NULL
  if (n.o>1){
    name.x <- c("time",name.out)
    names(name.x <- name.x)
    
    if (length(y.o)==1){
      if (select$y==TRUE){
        nout <- vect2str(name.out)
        out.txt <- paste0(out.txt, 'checkboxGroupInput("out1", label="plot", choices=',
                          nout,', selected=',nout,'),')
      }
      if (select$x==TRUE)
        out.txt <- paste0(out.txt,'\nselectInput("x1", label="v.s.",',vect2str(name.x),'),')       
      if ((select$y==TRUE) | (select$x==TRUE))
        out.txt <- paste0(out.txt, 'hr(),\n')
    }else{
      for (k in (1:length(y.o))){
        if (select$y==TRUE){
          noutk <- vect2str(y.o[[k]]$name)
          out.txt <- paste0(out.txt,'checkboxGroupInput("out',k,'", label="plot ',k,'", choices=',
                            noutk,', selected=',noutk,'),\n') 
        }
        if (select$x==TRUE)
          out.txt <- paste0(out.txt,'\nselectInput("x',k,'", label="v.s.",',vect2str(name.x),'),\n')       
        if ((select$y==TRUE) | (select$x==TRUE))
          out.txt <- paste0(out.txt, 'hr(),\n')
      }
    }
    if (select$ref | select$y)
      out.txt <- paste0(out.txt, 'strong("add"),\n')
    if (i.leg)
      out.txt <- paste0(out.txt,'checkboxInput("legend", label="legend", value=TRUE),\n') 
  }
  
  if (select$ref){
    sref1 <- ('
fluidRow(
column(5,checkboxInput("boxref", label="%s")),
column(4,actionButton("butref", label = "Reset"))
),
')  
    if (identical(style,'navbar1'))
      sref1 <- sprintf(sref1,"reference")
    else
      sref1 <- sprintf(sref1,"ref.")  
    out.txt <- paste0(out.txt, sref1,'hr(),\n')
  }
  if (select$log)
    out.txt <- paste0(out.txt,'radioButtons("ilog", "scale", c("linear" = FALSE,"log" = TRUE), inline=TRUE),')       
  
  if ((!is.null(out.txt)) & (!identical(style,"basic")) ){
    out.txt <- paste0(out.txt,'\nbr()')
    i.out <- TRUE
  } else {
    i.out <- FALSE
  }
  
  if (identical(style,"basic")){
    if (is.null(out.txt))
      tc <- tab.content[[1]]
    else
      tc <- paste0(out.txt,'\nhr(),\n',tab.content[[1]])
    if (length(tab.content)>1){
      for (k in (2:length(tab.content)))
        tc <- list(paste0(tc,'\nhr(),\n',tab.content[[k]]))
    }
    tc <- paste0(tc,'br()')
    tb <- NULL
    to <- NULL
  } else if (identical(style,"navbar1")) {
    tc <- c(tab.content,out.txt)
    tb <- c(tab.label,"output")
    to <- NULL
  } else {
    tc <- tab.content
    tb <- tab.label
    to <- out.txt
  }
  if (length(tc)>1){
    txt <- paste0("tabsetPanel(type='",tabstyle,"',\n")
    for (k in (1:length(tc)))
      txt <- paste0(txt,"tabPanel('",tb[k],"',br(),\n",tc[[k]],"\n),\n")
    txt <- paste0(txt,'br()\n)\n')
    tc <- txt
  }
  tab <- list(content=tc, output=to, i.output=i.out)
  
  return(tab)
}
#--------------------------------------------------------------

widgetui <- function(xkj,id, nxkj){
  if (xkj$widget=="slider")
    uikj <-paste0('sliderInput("',id,'", label="',nxkj,'"')
  else if (xkj$widget=="numeric")                 
    uikj <-paste0('numericInput("',id,'", label="',nxkj,'"')
  else if (xkj$widget=="select")                 
    uikj <-paste0('selectInput("',id,'", label="',nxkj,'", choices=',vect2str(xkj$choices))
  
  for (k in (1:length(xkj))){
    if (names(xkj[k]) %in% c("value", "min","max","step","selected"))
      uikj <- paste0(uikj, ",", names(xkj[k]),"=", xkj[[k]] )
  }
  uikj <- paste0(uikj, ")")  
  return(uikj)
}

#-------------------------------------------------------

set.settings <- function(s,style){
  if (is.null(s)){
    if (identical(style,"basic"))
      r <- FALSE
    else
      r <- TRUE}
  else  
    r <- s
  return(r)
}
