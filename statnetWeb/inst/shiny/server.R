
library(statnetWeb)
library(RColorBrewer)
library(lattice)
library(latticeExtra)

data(faux.mesa.high)
#data(faux.magnolia.high)
data(florentine)
data(sampson)
data(samplk)
data(ecoli)
data(molecule)
data(kapferer)

BRGcol <- "darkred"
CUGcol <- "darkorange"
obsblue <- "#076EC3"
histblue <- "#83B6E1"
tgray3 <- adjustcolor("gray", alpha.f = 0.3)
tgray7 <- adjustcolor("gray", alpha.f = 0.7)

allterms <- splitargs(searchterm = "")

shinyServer(
  function(input, output, session){

    oldoptions <- options()
    on.exit(options(oldoptions))
    options(digits=3)


# Reactive Expressions ----------------------------------------------------
# These expressions contain most of the code from the ergm package that we will
# be using. Objects created with a reactive expression can be accessed from any
# other reactive expression or render functions and they only get re-run when
# their values are outdated. Since many of our render functions will be calling
# the same ergm objects, using reactive expressions will help the app run much
# faster.


values <- reactiveValues()

# when two options are available to the user, or when we need to know if one
# variable is outdated this reactive value will keep track of the state
state <- reactiveValues(symmdir = FALSE, plotperc_dd = FALSE,
                        plotperc_gd = FALSE, allterms = FALSE, gof = 0)

# To keep a list of all attributes uploaded by the user:
values$v_attrNamesToAdd <- list(1)
values$v_attrValsToAdd <- list()
values$e_attrNamesToAdd <- list(1)
values$e_attrValsToAdd <- list()
values$ev_attrNamesToAdd <- list(1)
values$ev_attrValsToAdd <- list()
values$input_termslist <- list()

#move to Help page when user clicks Help link button
observe({
  if(input$helpLink == 0) {return()}
  isolate({
    updateTabsetPanel(session, 'navbar', selected='tab8')
  })
})

#move to Data panel when user clicks Get Started button
observe({
  if(input$startButton == 0) {return()}
  isolate({
    updateTabsetPanel(session, 'navbar', selected='tab2')
  })
})

#update active tab in navbar when arrows are clicked
leftarrowclicks <- reactive({
  input$dataleft+input$plotleft+input$fitleft+input$mcmcleft+input$gofleft+input$simleft
})
rightarrowclicks <- reactive({
  input$dataright+input$plotright+input$fitright+input$mcmcright+input$gofright+input$simright
})
observe({
  if(leftarrowclicks() == 0) {return()}
  tabOptions <- c('tab1', 'tab2', 'tab3', 'tab4', 'tab5', 'tab6', 'tab7', 'tab8')
  current <- isolate(which(input$navbar==tabOptions))
  updateTabsetPanel(session, 'navbar', selected=tabOptions[current-1])
})
observe({
  if(rightarrowclicks() == 0) {return()}
  tabOptions <- c('tab1', 'tab2', 'tab3', 'tab4', 'tab5', 'tab6', 'tab7', 'tab8')
  current <- isolate(which(input$navbar==tabOptions))
  updateTabsetPanel(session, 'navbar', selected=tabOptions[current+1])
})


#nwinit is used to get the initial values of the network
nwinit <- reactive({
  #input$rawdatafile comes as a dataframe with name, size, type and datapath
  #datapath is stored in 4th column of dataframe
  #network creates a network object from the input file
  if(is.null(input$rawdatafile)){
    nw_var <- NULL
  } else {
    filepath <- input$rawdatafile[1,4]
    filename <- input$rawdatafile[1,1]
    fileext <- substr(filename,nchar(filename)-3,nchar(filename))

    if(input$filetype == 1){
      validate(
        need(fileext %in% c(".rds", ".Rds", ".RDs", ".RDS"),
             "Upload an .rds file"))
      nw_var <- readRDS(paste(filepath))
    } else if(input$filetype == 2){
      validate(
        need(fileext %in% c(".net", ".NET"),
             "Upload a .net file"))
      nw_var <- read.paj(paste(filepath))
    } else if(input$filetype == 3){
      validate(
        need(fileext %in% c(".paj",".PAJ"),
             "Upload a .paj file"))
      nws <- read.paj(paste(filepath))
      if(!is.null(pajnws())){
        nw_var <- nws$networks[[as.numeric(input$choosepajnw)]]
      }
    } else if(input$filetype == 4){
      validate(
        need(fileext %in% c(".csv",".CSV") |
               fileext %in% c(".rds", ".Rds", ".RDs", ".RDS"),
             "Upload the specified type of matrix"))
      if(fileext %in% c(".csv",".CSV")){
        header <- TRUE
        row_names<-1
        if(input$matrixtype == "edgelist"){
          header <- FALSE
          row_names<-NULL
        }
        try({nw_var <- network(read.csv(paste(filepath), sep=",", header=header,
                                        row.names=row_names),
                          directed=input$dir, loops=input$loops,
                          multiple=input$multiple, bipartite=input$bipartite,
                          matrix.type=input$matrixtype,
                          ignore.eval=FALSE, names.eval='edgevalue')
             })

      } else if(fileext %in% c(".rds", ".Rds", ".RDs", ".RDS")){
        newmx <- readRDS(paste(filepath))
        nw_var <- network(newmx,
                        directed=input$dir, loops=input$loops,
                        multiple=input$multiple, bipartite=input$bipartite,
                        matrix.type=input$matrixtype,
                        ignore.eval=FALSE, names.eval='edgevalue')

      }
    }
  }
  if(input$filetype == 5){
    if(input$samplenet == ""){
      nw_var <- NULL
    } else {
      nw_var <- eval(parse(text = input$samplenet))
      if(!is.element('bipartite',names(nw_var$gal))){
        set.network.attribute(nw_var,'bipartite',FALSE)
      }
    }
  }
  return(nw_var)
})

#list of everything in an uploaded Pajek project
pajnws <- reactive({
  nws <- NULL
  if((input$filetype == 3) & (!is.null(input$rawdatafile))){
    filename <- input$rawdatafile[1,1]
    if(substr(filename,nchar(filename)-3,nchar(filename))==".paj"){
      nws <- read.paj(paste(input$rawdatafile[1,4]))
    }
  }
  nws
})

nwname <- reactive({
  name <- input$rawdatafile[1,1]
  if(input$filetype == 5){
    name <- input$samplenet
  }
  name
})

#number of nodes in nw
nodes <- reactive({
  if(!is.network(nwinit())){return()}
  nwinit()$gal$n
})

#number of edges in initial nw
nedgesinit <- reactive({
  if(!is.network(nwinit())) return()
  network.edgecount(nwinit())
})

#initial vertex attributes
vattrinit <- reactive({
  vattrinit <- c()
  if(is.network(nwinit())){
    vattrinit<-list.vertex.attributes(nwinit())
  }
  vattrinit
})

#matrix of vertex attribute values
vattrinit.vals <- reactive({
  v <- list()
  for (j in seq(length(list.vertex.attributes(nwinit())))) {
    v[[j]] <- get.vertex.attribute(nwinit(), vattrinit()[j])
  }
  v
})

#set correct number of rows for the attribute value lists,
#so that we can add columns later
observe({
  nwinit()
  #reset lists when uploaded network changes
  vdf <- list()
  edf <- list()
  evdf <- list()
  if (is.network(nwinit())){
    n <- nodes()
    e <- nedgesinit()
    for (i in 1:n){
      vdf <- rbind(vdf,i)
    }
    for (i in 1:e){
      edf <- rbind(edf,i)
      evdf <- rbind(evdf,i)
    }
    values$v_attrValsToAdd <- vdf
    values$e_attrValsToAdd <- edf
    values$ev_attrValsToAdd <- evdf
    values$v_attrNamesToAdd <- list(1)
    values$e_attrNamesToAdd <- list(1)
    values$ev_attrNamesToAdd <- list(1)

    values$vertexnames <- network.vertex.names(nwinit())
  }
})

#names of uploaded attributes
#or helpful message that upload is incorrect
newattrnamereac <- reactive({
  newname <- ''
  try({
    path <- input$newattrvalue[1,4]
    filename <- input$newattrvalue[1,1]
    fileext <- substr(filename,nchar(filename)-3,nchar(filename))

    if(fileext %in% c(".csv", ".CSV") ){
      newattrs <- read.csv(paste(path), sep=",", header=TRUE,
                           stringsAsFactors=FALSE)
      newname <- names(newattrs)
      if(input$newattrtype == "edgeattr" & input$edgeform == "matrix"){
        newname <- substr(filename, 1, nchar(filename)-4)
      }
    } else if(fileext %in% c(".rds",".Rds",".RDs",".RDS") ){
      newattrs <- readRDS(paste(path))
      newname <- names(newattrs)
      if(class(newattrs) != "list"){
        newname <- "Attribute is not compatible, see help buttons and try again"
      }
    } else {
      newname <- "Attribute is not compatible, see help buttons and try again"
    }

  })
  if(is.null(newname)){
    newname <- "Attribute is not named,  please fix and re-upload"
  }
  newname
})

#save new vertex names
observeEvent(input$newattrButton, {
    if(input$newattrtype == "vertexnames"){
      path <- input$newattrvalue[1,4]
      filename <- input$newattrvalue[1,1]
      fileext <- substr(filename,nchar(filename)-3,nchar(filename))
      if(fileext %in% c(".csv", ".CSV")){
        newnames <- read.csv(paste(path), sep=",", header=TRUE,
                             stringsAsFactors=FALSE)
        newnames <- newnames[[1]]
      } else if(fileext %in% c(".rds",".Rds",".RDs",".RDS")){
        newnames <- readRDS(paste(path))
      }

      values$vertexnames <- newnames
    }
})

#add vertex attributes to list
observeEvent(input$newattrButton, {
      if(input$newattrtype == "vertexattr"){
          path <- input$newattrvalue[1,4]
          filename <- input$newattrvalue[1,1]
          fileext <- substr(filename,nchar(filename)-3,nchar(filename))
          if(fileext %in% c(".csv", ".CSV")){
            newattrs <- read.csv(paste(path), sep=",", header=TRUE,
                                 stringsAsFactors=FALSE)
            newname <- names(newattrs)
          } else if(fileext %in% c(".rds",".Rds",".RDs",".RDS")){
            newattrs <- readRDS(paste(path))
            newname <- names(newattrs)
          }

          namesofar <- values$v_attrNamesToAdd
          valsofar <- values$v_attrValsToAdd
          for(k in 1:length(newname)){
            namesofar <- cbind(namesofar, newname[[k]])
            valsofar <- cbind(valsofar, newattrs[[k]])
          }

          values$v_attrNamesToAdd <- namesofar
          values$v_attrValsToAdd <- valsofar
        }
})

#add edge attributes to list
observeEvent(input$newattrButton, {
    if(input$newattrtype == "edgeattr" & input$edgeform == "vector"){
      path <- input$newattrvalue[1,4]
      filename <- input$newattrvalue[1,1]
      fileext <- substr(filename,nchar(filename)-3,nchar(filename))
      if(fileext %in% c(".csv", ".CSV")){
        newattrs <- read.csv(paste(path), sep=",", header=TRUE,
                             stringsAsFactors=FALSE)
        newname <- names(newattrs)
      } else if(fileext %in% c(".rds",".Rds",".RDs",".RDS")){
        newattrs <- readRDS(paste(path))
        newname <- names(newattrs)
      }

      namesofar <- values$e_attrNamesToAdd
      valsofar <- values$e_attrValsToAdd
      for(k in 1:length(newname)){
        namesofar <- cbind(namesofar, newname[[k]])
        valsofar <- cbind(valsofar, newattrs[[k]])
      }
        values$e_attrNamesToAdd <- namesofar
        values$e_attrValsToAdd <- valsofar
      }
})

#add edge values to list
observeEvent(input$newattrButton, {
    if(input$newattrtype == "edgeattr" & input$edgeform == "matrix"){
      path <- input$newattrvalue[1,4]
      filename <- input$newattrvalue[1,1]
      fileext <- substr(filename,nchar(filename)-3, nchar(filename))
      if(fileext %in% c(".csv", ".CSV")){
        newattrs <- read.csv(paste(path), sep=",", header=TRUE,
                             row.names = 1,
                             stringsAsFactors=FALSE)
        newname <- substr(filename, 1, nchar(filename)-4)
        newattrs <- data.matrix(newattrs, rownames.force=FALSE)
      } else if(fileext %in% c(".rds",".Rds",".RDs",".RDS")){
        newattrs <- readRDS(paste(path))
        newname <- names(newattrs)
      }
      namesofar <- values$ev_attrNamesToAdd
      valsofar <- values$ev_attrValsToAdd
      j <- length(valsofar)
      for(k in 1:length(newname)){
        namesofar <- cbind(namesofar, newname[[k]])
        valsofar[[j+k]] <- newattrs[[k]]
      }
      values$ev_attrNamesToAdd <- namesofar
      values$ev_attrValsToAdd <- valsofar
    }
})

observeEvent(input$symmdir,{
  state$symmdir <- TRUE
})
observeEvent(input$symmundir,{
  state$symmdir <- FALSE
})

#attributes will be added to this network
nwmid <- reactive({
    nw_var <- nwinit()

    if (class(nw_var)=="network"){
      #preserve initial network attributes and let user choose if directed
      #after symmetrizing
      if(input$symmetrize != "Do not symmetrize"){
        symnw <- sna::symmetrize(nw_var, rule=input$symmetrize)
        nw_var <- network(symnw, matrix.type="adjacency", directed=state$symmdir,
                          hyper=nwattrinit()[2], loops=nwattrinit()[3],
                          multiple=nwattrinit()[4], bipartite=nwattrinit()[5])
        #add initial vertex attributes back after symmetrizing
        #can't add edge attributes back because number of edges has changed
        for(k in 1:length(vattrinit())){
          attr_names <- vattrinit()
          attr_list <- vattrinit.vals()
          set.vertex.attribute(nw_var, attr_names[k], attr_list[[k]])
        }
      }

      if (is.bipartite(nw_var)){
        set.vertex.attribute(nw_var, "mode", c(rep(1, nw_var$gal$bipartite),
                                               rep(2, nw_var$gal$n - nw_var$gal$bipartite)))
      }

      v_attrNamesToAdd <- values$v_attrNamesToAdd
      v_attrValsToAdd <- values$v_attrValsToAdd
      e_attrNamesToAdd <- values$e_attrNamesToAdd
      e_attrValsToAdd <- values$e_attrValsToAdd
      ev_attrNamesToAdd <- values$ev_attrNamesToAdd
      ev_attrValsToAdd <- values$ev_attrValsToAdd


      if(input$newattrButton > 0){
        try({network.vertex.names(nw_var) <- values$vertexnames})
      }
      v_numnew <- length(v_attrNamesToAdd)
      if(v_numnew > 1){
        for (j in 2:v_numnew){
          try({v_newname <- as.character(v_attrNamesToAdd[1,j])
          v_newval <- v_attrValsToAdd[,j]
          set.vertex.attribute(nw_var,v_newname,v_newval)})
        }
      }

      e_numnew <- length(e_attrNamesToAdd)
      if(e_numnew > 1){
        for (k in 2:e_numnew){
          try({e_newname <- as.character(e_attrNamesToAdd[1,k])
          e_newval <- e_attrValsToAdd[,k]
          set.edge.attribute(nw_var,e_newname,e_newval)})
        }
      }

      ev_numnew <- length(ev_attrNamesToAdd)
      if(ev_numnew > 1){
        for (l in 2:ev_numnew){
          try({ev_newname <- as.character(ev_attrNamesToAdd[1,l])
          ev_newval <- ev_attrValsToAdd[[l]]
          set.edge.value(nw_var,ev_newname,ev_newval)})
        }
      }
    }

    nw_var
	})

#use this network for future calculations
nw <- reactive({
  nw_var <- nwmid()

#deleting attributes is no longer available

  values$input_termslist <- list()
  updateTextInput(session, inputId='terms', value='edges')

  nw_var
})

elist <- reactive({
  if(!is.network(nwinit())) return()
  as.edgelist(nw())
})

#get coordinates to plot network with
coords <- reactive({
  input$refreshplot
  plot.network(nw())
  })

#initial network attributes
#returns vector of true/falses
nwattrinit <- reactive({
  if(!is.network(nwinit())){return()}
  nwattributes <- c('directed','hyper','loops','multiple','bipartite')
  unlist(lapply(nwattributes,get.network.attribute,x=nwinit()))
})

#list of all vertex attributes in nw (after adding new)
attrib <- reactive({
	  attr <- c()
    if(is.network(nw())){
		    attr<-list.vertex.attributes(nw())
    }
      attr
  })

#don't allow "na" or "vertex.names" as vertex attributes in menus on fit tab
menuattr <- reactive({
  menuattr <- attrib()
  if(is.element("na",menuattr)){
    menuattr <- menuattr[-which("na"==menuattr)]
  }
  if(is.element("vertex.names",menuattr)){
    menuattr <- menuattr[-which("vertex.names"==menuattr)]
  }
  menuattr
})

#numeric attributes only (for size menu, etc.)
numattr <- reactive({
    numattr <- c()
    if(is.network(nw())){
      for(i in 1:length(attrib())){
        if(is.numeric(get.vertex.attribute(nw(),attrib()[i]))){
          numattr <- append(numattr,attrib()[i])
        }
      }}
    numattr
})

#dataframe of nodes, their attributes, and their coordinates in nwplot
nwdf <- reactive({
  attrs <- menuattr()
  if(is.na(as.numeric(network.vertex.names(nw()))[1])){
    df <- data.frame(Names = network.vertex.names(nw()))
  } else {
    df <- data.frame(Names = as.numeric(network.vertex.names(nw())))
  }
  for(i in seq(length(attrs))){
    df[[attrs[i]]] <- get.vertex.attribute(nw(), attrs[i])
  }
  df[["Missing"]] <- get.vertex.attribute(nw(), "na")
  df[["cx"]] <- coords()[,1]
  df[["cy"]] <- coords()[,2]
  df
})

# betweenness centrality of all nodes (for sizing menu)
nodebetw <- reactive({
  if(!is.network(nw())){return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
    cmode <- 'directed'
  } else {
    gmode <- 'graph'
    cmode <- 'undirected'
  }
  sna::betweenness(nw(), gmode=gmode, diag=has.loops(nw()),

              cmode=cmode)
})

nodesize <- reactive({
  if(!is.network(nw())){return()}
  nw_var <- nw()
  #scale size of nodes onto range between .7 and 3.5
  if (input$sizeby == '1'){
    size = 1
  } else if (input$sizeby == 'Betweenness'){
    minsize <- min(nodebetw())
    maxsize <- max(nodebetw())
    size = (nodebetw()-minsize)/(maxsize-minsize)*(3.5-.7)+.7
  } else {
    minsize <- min(get.vertex.attribute(nw_var,input$sizeby))
    maxsize <- max(get.vertex.attribute(nw_var,input$sizeby))
    size = (get.vertex.attribute(nw_var,input$sizeby)-minsize)/(maxsize-minsize)*(3.5-.7)+.7
  }
  size})

#vertex color
vcol <- reactive({
  if(!is.network(nw())){return()}
  nw_var <- nw()
  if(input$colorby == 2){
    vcol <- rep(2, nodes())
  } else {
    full_list <- get.vertex.attribute(nw_var,input$colorby)
    short_list <- sort(unique(full_list))
    ncolors <- length(short_list)
    if(is.element("Other", short_list)){ #to be consistent with order of legend
      short_list <- short_list[-which(short_list=="Other")]
      short_list <- c(short_list, "Other")
    }
    full_list <- match(full_list, short_list)
    #each elt corresponds to integer position in short_list
    pal <- c('red', 'blue', 'green3', 'cyan', 'magenta3',
             'yellow', 'orange', 'black', 'grey')
    if(ncolors>9){
      pal <- colorRampPalette(brewer.pal(11,"RdYlBu"))(ncolors)
    }
    vcol <- pal[full_list]
  }
  vcol
})

legendlabels <- reactive({
  if(!is.network(nw())){return()}
  nw_var <- nw()
    if(input$colorby == 2){
      legendlabels <- NULL
    }else{
      legendlabels <- sort(unique(get.vertex.attribute(nw_var, input$colorby)))
      if(is.element("Other", legendlabels)){
        legendlabels <- legendlabels[-which(legendlabels=="Other")]
        legendlabels <- c(legendlabels, "Other")
      }
    }
    legendlabels
  })

legendfill <- reactive({
  if(input$colorby == 2){
    legendfill <- NULL
  } else {
    n <- length(legendlabels())
    pal <- c('red', 'blue', 'green3', 'cyan', 'magenta3',
             'yellow', 'orange', 'black', 'grey')
    if(n>9){
      pal <- colorRampPalette(brewer.pal(11,"RdYlBu"))(n)
    }
    legendfill <- adjustcolor(pal, alpha.f = input$transp)
  }
  legendfill
})

#simulated graphs for cug tests
observeEvent(c(nw(), input$ncugsims),{
  if(!is.null(nw())){
    s <- network.edgecount(nw())
    if (is.directed(nw())){
      mode <- "digraph"
    } else {
      mode <- "graph"
    }

    brgsims <- sna::rgraph(n = nodes(), m = input$ncugsims, tprob = sna::gden(nw()), mode = mode,
                          diag = nw()$gal$loops)
    cugsims <- sna::rgnm(n = input$ncugsims, nv = nodes(), m = s, mode = mode,
                        diag = nw()$gal$loops)

    values$cugsims <- list(brgsims, cugsims)
  }
})

brgvals <- reactive({
  apply(values$cugsims[[1]], MARGIN = 1, FUN = cugstats, term = input$cugtestterm,
        directed = nw()$gal$directed, loops = nw()$gal$loops)
})
cugvals <- reactive({
  apply(values$cugsims[[2]], MARGIN = 1, FUN = cugstats, term = input$cugtestterm,
        directed = nw()$gal$directed, loops = nw()$gal$loops)
})

#add terms to list as user enters them
#function in alert.js will click the addtermButton when user
#presses Enter from within the terms textbox
observe({
  if(input$addtermButton==0) {return()}
  isolate({
    valsofar <- values$input_termslist
    newval <- input$terms
    values$input_termslist <- rbind(valsofar, newval)
    updateTextInput(session, inputId='terms', value='')
  })
})

observe({
  if(input$resetformulaButton==0) {return()}
  isolate({
    values$input_termslist <- list()
  })
})

ergm.terms <- reactive({
  nw()
  input$resetformulaButton
  input$addtermButton
  interms <- values$input_termslist
  if(length(interms)==0) {return('NA')}
  paste(interms, collapse = '+')
})

ergm.formula <- reactive({
  if(ergm.terms()=='NA') {return()}
  formula(paste('nw() ~ ', ergm.terms(), sep = ''))})

current.simterms <- reactive({
  mod <- input$choosemodel_sim
  if(mod=="Current"){
    terms <- ergm.terms()
  } else {
    mod <- as.numeric(substr(mod,6,6))
    terms <- values$modelformulas[[mod]]
  }
  terms
})

current.simformula <- reactive ({
  terms <- current.simterms()
  formula(paste('nw() ~ ', terms, sep=''))
})

#ergm model object
model1reac <- reactive({
  if(input$fitButton == 0){
    return()
  }
  usingdefault <- isolate(input$controldefault)
  if(usingdefault){
    mod <- isolate(ergm(ergm.formula()))
  } else {
    customcontrols <- isolate(paste(input$customMCMCcontrol, sep=","))
    if(customcontrols == ""){
      mod <- isolate(ergm(ergm.formula(),
                          control=control.ergm(MCMC.interval=isolate(input$MCMCinterval),
                                               MCMC.burnin=isolate(input$MCMCburnin),
                                               MCMC.samplesize=isolate(input$MCMCsamplesize))))
    } else {
      mod <- isolate(ergm(ergm.formula(),
                          control=control.ergm(MCMC.interval=isolate(input$MCMCinterval),
                                               MCMC.burnin=isolate(input$MCMCburnin),
                                               MCMC.samplesize=isolate(input$MCMCsamplesize),
                                               eval(parse(text=customcontrols)))))
    }

  }
  return(mod)
  })

# values$modelstate keeps track of whether model fit is oudated,
# equal to 0 if the nw changes and a new model has not been fit yet
# equal to 1 otherwise
values$modelstate <- 0

#Keep track of saved models
values$modeltotal <- 0
values$modelcoefs <- list()
values$modelformulas <- list()
values$modelfits <- list()
values$modelsumstats <- list()

observeEvent(input$savemodelButton, {
  if(input$fitButton > 0){
    m <- isolate(values$modeltotal)
    if(m < 5){
      #increment label on save model button
      values$modeltotal <- m+1
    }
  }
})
observeEvent(values$modeltotal, {
  # Add to lists that hold info for each model
  # values$modelcoefs is a list object, each element is a list of
  # coefficients and stars for a single model
  m <- values$modeltotal
  if(m > 0 & m <= 5){
    values$modelcoefs[[m]] <- ergm.info(model1reac())
    values$modelformulas[[m]] <- ergm.terms()
    values$modelfits[[m]] <- model1reac()
    values$modelsumstats[[m]] <- summary(ergm.formula())
  }
})

observeEvent(input$clearmodelButton, {
  #clear saved models after button click
  values$modeltotal<-isolate(0)
  values$modelcoefs <- list()
  values$modelformulas <- list()
  values$modelfits <- list()
  values$modelsumstats <- list()
})
observe({
  #clear saved models when network changes
  if(values$modelstate==0){
    values$modeltotal<-isolate(0)
    values$modelcoefs <- list()
    values$modelformulas <- list()
    values$modelfits <- list()
    values$modelsumstats <- list()
  }
})

model1gof <- reactive({
  input$gofButton
  mod <- input$choosemodel_gof
  if(mod == "Current"){
    mod <- model1reac()
  } else {
    mod <- values$modelfits[[1]]
  }
  if(input$gofterm == 'Default'){
    #use default gof formula
    model1gof <- gof(mod)
  } else {
    gofform <- formula(paste('mod ~ ', input$gofterm, sep = ''))
    model1gof <- gof(gofform)
  }
  return(model1gof)
})

model2gof <- reactive({
  input$gofButton
  if(values$modeltotal < 2){
    return()
  } else {
    mod <- values$modelfits[[2]]
  }
  if(input$gofterm == 'Default'){
    #use default gof formula
    model2gof <- gof(mod)
  } else {
    gofform <- formula(paste('mod ~ ', input$gofterm, sep = ''))
    model2gof <- gof(gofform)
  }
  return(model2gof)
})

model3gof <- reactive({
  input$gofButton
  if(values$modeltotal < 3){
    return()
  } else {
    mod <- values$modelfits[[3]]
  }
  if(input$gofterm == 'Default'){
    #use default gof formula
    model3gof <- gof(mod)
  } else {
    gofform <- formula(paste('mod ~ ', input$gofterm, sep = ''))
    model3gof <- gof(gofform)
  }
  return(model3gof)
})

model4gof <- reactive({
  input$gofButton
  if(values$modeltotal < 4){
    return()
  } else {
    mod <- values$modelfits[[4]]
  }
  if(input$gofterm == 'Default'){
    #use default gof formula
    model4gof <- gof(mod)
  } else {
    gofform <- formula(paste('mod ~ ', input$gofterm, sep = ''))
    model4gof <- gof(gofform)
  }
  return(model4gof)
})

model5gof <- reactive({
  input$gofButton
  if(values$modeltotal < 5){
    return()
  } else {
    mod <- values$modelfits[[5]]
  }
  if(input$gofterm == 'Default'){
    #use default gof formula
    model5gof <- gof(mod)
  } else {
    gofform <- formula(paste('mod ~ ', input$gofterm, sep = ''))
    model5gof <- gof(gofform)
  }
  return(model5gof)
})

allmodelsimreac <- reactive({
  input$simButton
  mod <- input$choosemodel_sim
  if(mod=="Current"){
    mod <- model1reac()
  } else {
    mod <- as.numeric(substr(mod,6,6))
    mod <- values$modelfits[[mod]]
  }

  if(input$simcontroldefault){
    s<-isolate(simulate(mod, nsim = input$nsims))
  } else {
    customcontrols <- isolate(paste(input$simcustomMCMCcontrol, sep=","))
    if(customcontrols == ""){
      s<-isolate(simulate(mod, nsim = input$nsims,
                       control=control.simulate.ergm(MCMC.burnin=input$simMCMCburnin,
                                                    MCMC.interval=input$simMCMCinterval)))
    } else {
      s<-isolate(simulate(mod, nsim = input$nsims,
                       control=control.simulate.ergm(MCMC.burnin=input$simMCMCburnin,
                                                    MCMC.interval=input$simMCMCinterval,
                                                    eval(parse(text=customcontrols)))))
    }
  }
  return(s)
})

# Currently, the reactive statements that control the sizing/coloring/legend in
# thesimulation plots use the attributes from the original network as a point
# of reference. If the method for simulating networks changes from applying the
# same distribution of attributes, these `get.vertex.attribute` commands for
# `minsize` and `maxsize` would also need to change.

#get coordinates to plot simulations with
sim.coords.1 <- reactive({
  input$simButton
  isolate(plot.network(allmodelsimreac()))})
sim.coords.2 <- reactive({
  plot.network(allmodelsimreac()[[input$thissim]])})

nodebetw2 <- reactive({
  if(input$nsims==1){
    if(is.directed(allmodelsimreac())){
      gmode <- "digraph"
      cmode <- "directed"
    } else {
      gmode <- "graph"
      cmode <- "undirected"
    }
    b <- sna::betweenness(allmodelsimreac(), gmode=gmode, cmode=cmode)
  } else {
    if(is.directed(allmodelsimreac()[[input$thissim]])){
      gmode <- "digraph"
      cmode <- "directed"
    } else {
      gmode <- "graph"
      cmode <- "undirected"
    }
    b <- sna::betweenness(allmodelsimreac()[[input$thissim]], gmode=gmode, cmode=cmode)
  }
  return(b)
})

nodesize2 <- reactive({
  nw_var <- nw()
  #scale size of nodes onto range between .7 and 3.5
  if (input$sizeby2 == '1'){
    size = 1
  } else if (input$sizeby2 == "Betweenness"){
    minsize <- min(nodebetw2())
    maxsize <- max(nodebetw2())
    size <- (nodebetw2()-minsize)/(maxsize-minsize)*(3.5-.7)+.7
  } else {
    minsize <- min(get.vertex.attribute(nw_var,input$sizeby2))
    maxsize <- max(get.vertex.attribute(nw_var,input$sizeby2))
    if(input$nsims==1){
    size <- (get.vertex.attribute(allmodelsimreac(),input$sizeby2)-minsize) /
            (maxsize-minsize) * (3.5-.7)+.7
  }else{
    size <- (get.vertex.attribute(allmodelsimreac()[[input$thissim]],input$sizeby2)-minsize) /
            (maxsize-minsize) * (3.5-.7)+.7
  }}
  return(size)
})

vcol2 <- reactive({
  if(!is.network(nw())){return()}
  input$simButton
  isolate(nsim <- input$nsims)
  if(!is.null(input$colorby2)){
  if(input$colorby2 ==2){
    vcol <- 2
  } else {
    if(nsim == 1){
      full_list <- get.vertex.attribute(allmodelsimreac(),input$colorby2)
    } else {
      full_list <- get.vertex.attribute(allmodelsimreac()[[input$thissim]],
                                        input$colorby2)
    }
    short_list <- sort(unique(full_list))
    ncolors <- length(short_list)
    if(is.element("Other", short_list)){ #to be consistent with order of legend
      short_list <- short_list[-which(short_list=="Other")]
      short_list <- c(short_list, "Other")
    }
    full_list <- match(full_list, short_list)
    #each elt is an integer position in short_list
    pal <- c('red', 'blue', 'green3', 'cyan', 'magenta3',
             'yellow', 'orange', 'black', 'grey')
    if(ncolors>9){
      pal <- colorRampPalette(brewer.pal(11,"RdYlBu"))(ncolors)
    }
    vcol <- pal[full_list]
  }
  }
  return(vcol)
})

legendlabels2 <- reactive({
  nw_var <- nw()
  if(!is.null(input$colorby2)){
  if(input$colorby2 == 2){
    legendlabels <- NULL
  }else{
    legendlabels <- sort(unique(get.vertex.attribute(nw_var, input$colorby2)))
    if(is.element("Other", legendlabels)){
      legendlabels <- legendlabels[-which(legendlabels=="Other")]
      legendlabels <- c(legendlabels, "Other")
    }
  }
  }
  return(legendlabels)
})

legendfill2 <- reactive({
  if(!is.null(input$colorby2)){
  if(input$colorby2 == 2){
    legendfill <- NULL
  } else {
    n <- length(legendlabels2())
    pal <- c('red', 'blue', 'green3', 'cyan', 'magenta3',
             'yellow', 'orange', 'black', 'grey')
    if(n > 9){
      pal <- colorRampPalette(brewer.pal(11,"RdYlBu"))(n)
    }
    legendfill <- adjustcolor(pal, alpha.f = input$transp2)
  }
  }
  return(legendfill)
})

# Output Expressions -------------------------------------------------------

# Every piece of content that gets displayed in the app has to be
# rendered by the appropriate `render*` function, e.g. `renderPrint` for text
# and `renderPlot` for plots. Most of the render functions here call
# reactive objects that were created above. I have divided the output objects
# into sections depending on what tab of the app they are called from.


# Data Upload -------------------------------------------------------------



output$datadesc <- renderUI({
  net <- input$samplenet
  text <- div()
  if(net == "ecoli1" | net == "ecoli2"){
    text <- div(
      p("The", code("ecoli", class = "codetxt"),
        "network data set comprises two versions of a",
        "biological network in which the nodes are operons in",
        em("Escherichia Coli"), "and a directed edge from one node to another",
        "indicates that the first encodes the transcription factor that",
        "regulates the second."),
      p("The network object", code("ecoli1", class = "codetxt"),
        "is directed, with 423 nodes", "and 519 ties. The object",
        code("ecoli2", class = "codetxt"), "is an undirected",
        "version of the same network, in which the five isolated nodes",
        "(which exhibit only self-regulation in",
        code("ecoli1", class = "codetxt"), "are removed, leaving 418 nodes."),
      p("The data set is based on the RegulonDB network (Salgado et al, 2001)",
        "and was modified by Shen-Orr et al (2002)."),
      strong("References"),
      p("Salgado et al (2001), Regulondb (version 3.2): Transcriptional",
        "Regulation and Operon Organization in Escherichia Coli K-12,",
        em("Nucleic Acids Research,"), "29(1): 72-74."),
      p("Shen-Orr et al (2002), Network Motifs in the Transcriptional",
        "Regulation Network of Escerichia Coli,", em("Nature Genetics,"),
        "31(1): 64-68.")
    )
  }
  if(net == "faux.mesa.high"){
    text <- div(
      p("This data set represents a simulation of an in-school friendship",
        "network. The network is named faux.mesa.high because the school",
        "commnunity on which it is based is in the rural western US, with a",
        "student body that is largely Hispanic and Native American."),
      p(code("faux.mesa.high", class = "codetxt"), "is a network object with",
        "205 vertices (students, in this case) and 203 undirected edges",
        "(mutual friendships)."),
      p("The vertex attributes are Grade, Sex, and Race. The Grade attribute",
        "has values 7 through 12, indicating each student's grade in school.",
        "The Race attribute is based on the answers to two questions, one on",
        "Hispanic identity and one on race, and takes six possible values:",
        "White (non-Hisp.), Black (non-Hisp.), Hispanic, Asian (non-Hisp.),",
        "Native American, and Other (non-Hisp.)"),
      p("The data set is based upon a model fit to data from one school",
        "community from the AddHealth Study, Wave I (Resnick et al., 1997).",
        "The processes for constructing the network are described in Hunter,",
        "Goodreau & Handcock (2008)"),
      strong("References"),
      p("Hunter D.R., Goodreau S.M. and Handcock M.S. (2008).",
        em("Goodness of Fit of Social Network Models, Journal of the American",
           "Statistical Association.")),
      p("Resnick M.D., Bearman, P.S., Blum R.W. et al. (1997).",
        em("Protecting adolescents from harm. Findings from the National",
           "Longitudinal Study on Adolescent Health, Journal of the American",
           "Medical Association,"), "278: 823-32.")
      )
  }
  if(net == "flobusiness" | net == "flomarriage"){
    text <- div(
      p("The two", code("florentine", class = "codetxt"), "networks are of",
        "marriage and business ties among Renaissance",
        "Florentine families. The data is originally from Padgett (1994) via",
        "UCINET and stored as", code("statnet", class = "codetxt"),
        "network objects."),
      p("Breiger & Pattison (1986), in their discussion of local role analysis,",
        "use a subset of data on the social relations among Renaissance",
        "Florentine families (person aggregates) collected by John Padgett from",
        "historical documents.", code("flobusiness", class = "codetxt"),
        "contains business ties - specifically, recorded",
        "financial ties such as loans, credits and joint partnerships.",
        code("flomarriage", class = "codetxt"), "contains marriage alliances."),
      p("As Breiger & Pattison point out, the original data are symmetrically",
        "coded. This is acceptable perhaps for marital ties, but is unfortunate",
        "for the financial ties (which are almost certainly directed). Both",
        "graphs provide vertex information on (1) each family's net wealth in",
        "1427 (in thousands of lira); (2) the number of priorates (seats on the",
        "civic council) held between 1282- 1344; and (3) the total number of",
        "business or marriage ties in the total dataset of 116 families",
        "(see Breiger & Pattison (1986), p 239)."),
      p("Substantively, the data include families who were locked in a struggle",
        "for political control of the city of Florence around 1430. Two",
        "factions were dominant in this struggle: one revolved around the",
        "infamous Medicis (9), the other around the powerful Strozzis (15)."),
      strong("References"),
      p("Wasserman, S. and Faust, K. (1994)",
        em("Social Network Analysis: Methods and Applications,"),
        "Cambridge University Press, Cambridge, England."),
      p("Breiger, R. and Pattison, P. (1986).",
        em("Cumulated social roles: The duality of persons and their algebras,"),
        "Social Networks, 8, 215-256.")
      )
  }
  if(net == "kapferer" | net == "kapferer2"){
    text <- div(
      p('This well-known social network dataset, collected by Bruce Kapferer',
        'in Zambia from June 1965 to August 1965, involves interactions among',
        'workers in a tailor shop as observed by Kapferer himself. Here, an',
        'interaction is defined by Kapferer as "continuous uninterrupted social',
        'activity involving the participation of at least two persons"; only',
        'transactions that were relatively frequent are recorded. All of the',
        'interactions in this particular dataset are "sociational", as opposed',
        'to "instrumental". Kapferer explains the difference (p. 164) as follows:'),
      p('"I have classed as transactions which were sociational in content those',
        'where the activity was markedly convivial such as general conversation,',
        'the sharing of gossip and the enjoyment of a drink together. Examples',
        'of instrumental transactions are the lending or giving of money,',
        'assistance at times of personal crisis and help at work."'),
      p("Kapferer also observed and recorded instrumental transactions, many of",
        "which are unilateral (directed) rather than reciprocal (undirected),",
        "though those transactions are not recorded here. In addition, there was",
        "a second period of data collection, from September 1965 to January 1966,",
        "but these data are also not recorded here. All data are given in",
        "Kapferer's 1972 book on pp. 176-179."),
      p("During the first time period, there were 43 individuals working in this",
        "particular tailor shop; however, the better-known dataset includes only",
        "those 39 individuals who were present during both time collection",
        "periods. (Missing are the workers named Lenard, Peter, Lazarus, and",
        "Laurent.) Thus, we give two separate networks here:",
        code("kapferer", class = "codetxt"), "is the well-known 39-individual",
        "dataset, whereas", code("kapferer2", class = "codetxt"), "is the full",
        "43-individual dataset."),
      strong("References"),
      p("Kapferer, Bruce (1972), Strategy and Transaction in an African Factory,",
        "Manchester University Press.")
    )
  }
  if(net == "molecule") {
    text <- div(
      p(code("molecule", class = "codetxt"),
        "is a synthetic network of 20 nodes that is used as an example within",
        "the", code("ergm", class = "codetxt"),
        "documentation. It has an interesting elongated shape - reminencent of",
        "a chemical molecule."))
  }
  if(net == "samplike" | net == "samplk1" | net == "samplk2" | net == "samplk3"){
    text <- div(
      p('Sampson (1969) recorded the social interactions among a group of monks',
        'while resident as an experimenter on vision, and collected numerous',
        'sociometric rankings. During his stay, a political “crisis in the',
        'cloister" resulted in the expulsion of four monks (Nos. 2, 3, 17, and',
        '18) and the voluntary departure of several others - most immediately,',
        'Nos. 1, 7, 14, 15, and 16. (In the end, only 5, 6, 9, and 11 remained).',
        'Of particular interest is the data on positive affect relations',
        '(“liking"), in which each monk was asked if they had positive',
        'relations to each of the other monks.'),
      p('The data were gathered at three times to capture changes in group',
        'sentiment over time:', code("samplk1, samplk2", class = "codetxt"), "and",
        code("samplk3.", class = "codetxt"), 'They represent three time points',
        'in the period during which a new cohort entered the monastery near the',
        'end of the study but before the major conflict began. Each member',
        'ranked only his top three choices on “liking." (Some subjects offered',
        'tied ranks for their top four choices). A tie from monk A to monk B',
        'exists if A nominated B as one of his three best friends at that that',
        'time point.'),
      p(code("samplk3", class = "codetxt"),
        "is a data set of Hoff, Raftery and Handcock (2002)."),
      p(code('samplike', class = "codetxt"),
        'is the time-aggregated graph. It is the cumulative tie for “liking"',
        'over the three periods. For this, a tie from monk A to monk B exists',
        'if A nominated B as one of his three best friends at any of the three',
        'time points.'),
      p('The graphs have three vertex attributes: ',
        tags$ul(
          tags$li('Groups of novices as classified by Sampson: "Loyal",',
                  '"Outcasts", and "Turks". There is also an interstitial',
                  'group not represented here.'),
          tags$li('An indicator of attendance the minor seminary of',
                  '“Cloisterville" before coming to the monastery.'),
          tags$li('The given names of the novices.')
          )),
      strong("References"),
      p("Sampson, S.F. (1968), A novitiate in a period of change:",
        em("An experimental and case study of relationships,"),
        "Unpublished Ph.D. dissertation, Department of Sociology,",
        "Cornell University."),
      p("White, H.C., Boorman, S.A. and Breiger, R.L. (1976).",
        em("Social structure from multiple networks. I. Blockmodels of roles",
           "and positions."), "American Journal of Sociology, 81(4), 730-780."),
      p("Wouter de Nooy, Andrej Mrvar, Vladimir Batagelj (2005)",
        em("Exploratory Social Network Analysis with Pajek,"),
        "Cambridge: Cambridge University Press")

      )
  }

  text
})

output$rawdatafile <- renderPrint({
  raw <- matrix(nrow=2,ncol=1)
  rownames(raw)<-c("name:", "size:")
  if(!is.null(input$rawdatafile)){
    raw[1,1] <- input$rawdatafile[1,1]
    raw[2,1] <- paste(input$rawdatafile[1,2], " bytes")
  }
  write.table(raw, quote=FALSE, col.names=FALSE)})

output$pajchooser <- renderUI({
  pajlist <- c(None = '')
  if(!is.null(pajnws())){
    pajlist <- 1:length(pajnws()$networks)
  names(pajlist) <- names(pajnws()$networks)
  }
  selectInput('choosepajnw',
              label = 'Choose a network from the Pajek project',
              choices = pajlist)
})
outputOptions(output, "pajchooser", suspendWhenHidden = FALSE)


output$newattrname <- renderPrint({
  if(!is.null(input$newattrvalue)){
      cat(newattrnamereac())}
})

# output$modifyattrchooser <- renderUI({
#   if(!is.network(nwmid())) {return()}
#   vattr <- list.vertex.attributes(nwmid())
#   eattr <- list.edge.attributes(nwmid())
#   attrlist <- c(vattr, eattr)
#   selectInput('modifyattrs', label=NULL, choices=attrlist)
# })


#summary of network attributes
output$nwsum <- renderPrint({
  if (is.null(nw())){
    return(cat('NA'))
  }
  nw_var <- nw()
  if (class(nw_var)!="network"){
    return(cat(nw_var))
  }
  return(nw_var)
})


# Network Descriptives ------------------------------------------------------

#NETWORK PLOT

#summary of network attributes
output$attr2 <- renderPrint({
  if (!is.network(nw())){
    return(cat('NA'))
  }
  nw_var <- nw()
  return(nw_var)
})

output$dynamiccolor <- renderUI({
  selectInput('colorby',
              label = 'Color nodes according to:',
              c('None' = 2, attrib()))
})
outputOptions(output,'dynamiccolor', suspendWhenHidden=FALSE, priority=10)

# need this to know when color palette will change
output$attrlevels <- renderText({
  return(length(legendlabels()))
})
outputOptions(output,'attrlevels', suspendWhenHidden=FALSE, priority=10)

# # reactivate colorwarning when network changes
# observeEvent(nw(),{
#   tags$script(HTML(
#     "document.getElementById('colorwarning1').style.display = 'block'"))
#   tags$script(HTML(
#     "document.getElementById('closewarning1').style.display = 'block'"))
# })

output$dynamicsize <- renderUI({
  selectInput('sizeby',
              label = 'Size nodes according to:',
              c('None' = 1, 'Betweenness', numattr()))
})
outputOptions(output,'dynamicsize',suspendWhenHidden=FALSE)

# The network plot takes display options from the sidebar of the ui. Even
# though I set the value of the 'None' option in the `sizeby` menu (above) as
# 1, it gets coerced into the string '1' by the rest of the strings in the
# vector of menu options. The variable `size` takes the value 1 if the user
# wants all the nodes to be the same size, and otherwise maps the values of the
# numeric attributes into the range between .7 and 3.5 using the formula
# y = (x-a)/(b-a) * (d-c) + c, where x is the input in some range [a,b]
# and y is the output in range [c,d].

output$nwplot <- renderPlot({
  if (!is.network(nw())){
    return()
  }
  input$plottabs
  input$rawdatafile
  input$samplenet

  nw_var <- nw()
  color <- adjustcolor(vcol(), alpha.f = input$transp)
  ecolor <- 1
  vborder <- 1
  vcex <- nodesize()
  if(is.bipartite(nw())){
    sides <- c(rep(50, nw()$gal$bipartite),
               rep(3, nodes() - nw()$gal$bipartite))
  } else{
    sides <- 50
  }

  if(!is.null(values$hoverpoints)){
    if(nrow(values$hoverpoints) > 0){
      nhov <- as.numeric(rownames(values$hoverpoints))
      vcex <- rep(1, nodes())
      vcex[nhov] <- 2
    }
  }
  if(!is.null(values$clickedpoints)){
    if(nrow(values$clickedpoints) > 0){
      nclick <- as.numeric(rownames(values$clickedpoints))
      color <- adjustcolor(vcol(), alpha.f = 0.4)
      color[nclick] <- vcol()[nclick]
      ecolor <- "lightgrey"
      vborder <- rep("lightgrey", nodes())
      vborder[nclick] <- 1
    }
  }
  if(!is.null(values$dblclickpoints)){
    if(nrow(values$dblclickpoints) > 0){
      ndbl <- as.numeric(rownames(values$dblclickpoints))
      neighb <- nw()[ndbl,] == 1
      color <- adjustcolor(vcol(), alpha.f = 0.4)
      color[ndbl] <- vcol()[ndbl]
      ecolor <- rep("lightgrey", nedgesinit())
      ecolor[apply(elist(), MARGIN = 1, FUN = function(x){any(x == ndbl)})] <- "black"
      vborder <- rep("lightgrey", nodes())
      vborder[neighb] <- "black"
      vborder[ndbl] <- "black"
    }
  }

  par(mar = c(0, 0, 0, 0))
  plot.network(nw_var, coord = coords(),
               displayisolates = input$iso,
               displaylabels = input$vnames,
               vertex.col = color,
               vertex.border = vborder,
               vertex.sides = sides,
               vertex.cex = vcex,
               edge.col = ecolor)
  if(input$colorby != 2){
    legend('bottomright', title = input$colorby, legend = legendlabels(),
           fill = legendfill(), bty='n')
  }

  if(!is.null(values$clickedpoints)){
    if(nrow(values$clickedpoints) > 0){
      #isolate(legend("topleft",
      #               legend = values$clickedpoints[, c("Names", menuattr())]))
      cx <- values$clickedpoints[, "cx"]
      cy <- values$clickedpoints[, "cy"]
      name <- values$clickedpoints[, "Names"]
      attrlabel <- paste("\n", menuattr())
      text(x = cx, y = cy,
           labels = paste0(name,
                           paste(attrlabel, values$clickedpoints[, menuattr()],
                                 collapse = "")),
           pos = 4, offset = 1)
    }
  }

})

observeEvent({c(input$plot_click, input$plot_dblclick)}, {
  values$clickedpoints <- nearPoints(nwdf(), input$plot_click,
                                     xvar = "cx", yvar = "cy",
                                     threshold = 10, maxpoints = 1)
})
observeEvent(input$plot_hover, {
  values$hoverpoints <- nearPoints(nwdf(), input$plot_hover,
                                   xvar = "cx", yvar = "cy",
                                   threshold = 10, maxpoints = 1)
})
observeEvent({c(input$plot_dblclick, input$plot_click)}, {
  values$dblclickpoints <- nearPoints(nwdf(), input$plot_dblclick,
                                      xvar = "cx", yvar = "cy",
                                      threshold = 10, maxpoints = 1)
})


output$nwplotdownload <- downloadHandler(
  filename = function(){paste(nwname(),'_plot.pdf',sep='')},
  content = function(file){
    pdf(file=file, height=10, width=10)
    nw_var <- nw()
    color <- adjustcolor(vcol(), alpha.f = input$transp)
    plot.network(nw_var, coord = coords(),
                 displayisolates = input$iso,
                 displaylabels = input$vnames,
                 vertex.col = color,
                 vertex.cex = nodesize())
    if(input$colorby != 2){
      legend('bottomright', title=input$colorby, legend = legendlabels(),
             fill = legendfill())
    }
    dev.off()
  }
  )

output$attrcheck <- renderUI({
  checkboxGroupInput("attrcols",
                     label = "Include these attributes",
                     choices = c(menuattr(), "Missing"),
                     selected = c(menuattr(), "Missing"))
})
outputOptions(output, "attrcheck", suspendWhenHidden = FALSE)

output$attrtbl_lg <- renderDataTable({
  dt <- nwdf()[, c("Names", input$attrcols)]
  dt
}, options = list(pageLength = 10))

output$attrtbl_sm <- renderPrint({
  ntbl <- length(input$attrcols)
  if(ntbl == 0){return()}
  attrname <- input$attrcols
  tbl_list <- list()
  if(ntbl == 1){
    tab <- attr.info(df = nwdf(), colname = attrname,
                     numattrs = numattr(), breaks = 10)
    tbl_list[[attrname]] <- tab
  } else {
    for(a in attrname){
      tab <- attr.info(df = nwdf(), colname = a,
                       numattrs = numattr(), breaks = 10)
      tbl_list[[a]] <- tab
    }
  }
  for(a in attrname){
    print(a, quote = FALSE)
    print(tbl_list[[a]])
  }
})

output$attrhist <- renderPlot({
  nplots <- length(input$attrcols)
  if(nplots == 0){return()}
  attrname <- input$attrcols
  if(nplots == 1){
    par(mfrow = c(1, 1))
    lvls <- length(unique(nwdf()[[attrname]]))
    if(input$attrhistaxis == "density" & attrname %in% numattr() & lvls > 9){
      plot(density(nwdf()[[attrname]]), main = attrname,
           col = "#076EC3", lwd = 2)
    } else {
      tab <- attr.info(df = nwdf(), colname = attrname,
                       numattrs = numattr(), breaks = 10)
      if(input$attrhistaxis == "percent"){
        tab <- tab/sum(tab)
      }
      barplot(tab, main = attrname, col = histblue)
    }
  } else {
    r <- ceiling(nplots/2)
    par(mfrow = c(r, 2))
    for(a in attrname){
      lvls <- length(unique(nwdf()[[a]]))
      if(input$attrhistaxis == "density" & a %in% numattr() & lvls > 9){
        plot(density(nwdf()[[a]]), main = a, col = "#076EC3", lwd = 2)
      } else {
        tab <- attr.info(df = nwdf(), colname = a,
                         numattrs = numattr(), breaks = 10)
        if(input$attrhistaxis == "percent"){
          tab <- tab/sum(tab)
        }
        barplot(tab, main = a, col = histblue)
      }
    }
  }
})

output$attrhistplotspace <- renderUI({
  nplots <- length(input$attrcols)
  r <- ceiling(nplots/2)
  h <- ifelse(r == 1, 400, r * 300)
  plotOutput("attrhist", height = h)
})

#Data to use for null hypothesis overlays in network plots
uniformsamples <- reactive({
  if(!is.network(nw())){
    return()
  }
  s <- network.edgecount(nw())
  if(is.directed(nw())){
    samples <- sna::rgnm(n=50, nv=nodes(), m=s, mode='digraph',
                        diag=has.loops(nw()))
  } else {
    samples <- sna::rgnm(n=50, nv=nodes(), m=s, mode='graph',
                        diag=has.loops(nw()))
  }
  samples
})

bernoullisamples <- reactive({
  if(!is.network(nw())){
    return()
  }
  density <- sna::gden(nw())
  if(is.directed(nw())){
    samples <- sna::rgraph(n=nodes(), m=50, mode='digraph', tprob=density,
                      diag=has.loops(nw()))
  } else {
    samples <- sna::rgraph(n=nodes(), m=50, mode='graph', tprob=density,
                      diag=has.loops(nw()))
  }
  samples
})

#DEGREE DISTRIBUTION

output$dynamiccmode_dd <- renderUI({
  menu <- c()
  if(is.network(nw())){
    menu <- c("total" = "freeman")
    if(is.directed(nw())){
      menu <- c("total" = "freeman",
                "indegree",
                "outdegree")
    }
  }
  selectInput("cmode_dd",
              label = "Type of degree",
              choices = menu)
})
outputOptions(output,'dynamiccmode_dd',suspendWhenHidden=FALSE, priority=10)

output$dynamiccolor_dd <- renderUI({
  menu <- menuattr()
  if(is.network(nw())){
    if(input$cmode_dd == "freeman" & is.directed(nw())){
      menu <- c()
    }
    selectInput('colorby_dd',
                label = 'Color bars according to:',
                c('None', menu),
                selected = 'None')
  }
})
outputOptions(output,'dynamiccolor_dd',suspendWhenHidden=FALSE, priority=10)

observe({
  if(!is.null(nw())){
  if(is.network(nw())){
    if(!is.directed(nw())){
      disableWidget('cmode', session, disabled=TRUE)
    } else {
      disableWidget('cmode', session, disabled=FALSE)
    }
  }}
})

dd_plotdata <- reactive({
  if(!is.network(nw())){
    return()
  }
  if(is.directed(nw())){
    gmode <- "digraph"
  } else {
    gmode <- "graph"
  }
  if(has.loops(nw())){
    diag <- TRUE
  } else {
    diag <- FALSE
  }
  deg <- sna::degree(nw(), gmode=gmode, cmode=input$cmode_dd, diag=diag)
  data <-tabulate(deg)
  data <- append(data,sum(deg==0),after=0)
  maxdeg <- max(deg)
  names(data) <- paste(0:maxdeg)

  #for color-coded bars
  if(!is.null(input$colorby_dd) & input$colorby_dd != "None"){
    if(is.directed(nw())){
      if(input$cmode_dd=='indegree'){
        data <- summary(nw() ~ idegree(0:maxdeg, input$colorby_dd))
      } else if(input$cmode_dd=='outdegree'){
        data <- summary(nw() ~ odegree(0:maxdeg, input$colorby_dd))
      } else {
        return('Cannot color code a directed graph using total degree.')
      }
    } else {
      data <- summary(nw() ~ degree(0:maxdeg, input$colorby_dd))
    }
    data <- t(matrix(data,nrow=maxdeg+1))
    colnames(data) <- 0:maxdeg
  }
  data
})


dd_uniformoverlay <- reactive({
  if(!is.network(nw())){
    return()
  }
  reps <- 50 #number of draws
  if(is.directed(nw())){
    deg <- sna::degree(uniformsamples(), g=1:reps, gmode='digraph', cmode=input$cmode_dd)
  } else {
    deg <- sna::degree(uniformsamples(), g=1:reps, gmode='graph', cmode=input$cmode_dd)
  }
    #now deg is a matrix where each element is a degree of a node
    #each column is a different draw

  degreedata <- apply(deg, MARGIN=2, FUN=tabulate, nbins=max(deg))
    #degreedata is matrix holding tabulation of degrees (except isolates)
    #each column is different draw
  z <- apply(deg, MARGIN=2, FUN=function(x){sum(x==0)}) #tabulation of isolates
  degreedata <- rbind(z, degreedata) #complete tabulation for each draw
  degreemeans <- apply(degreedata, MARGIN=1, FUN=mean)
  names(degreemeans) <- paste(0:max(deg))
  degreesd <- apply(degreedata, MARGIN=1, FUN=sd)
  mean_and_sd <- list(degreemeans, degreesd)
})

dd_bernoullioverlay <- reactive({
  if(!is.network(nw())){
    return()
  }
  reps = 50
  density <- sna::gden(nw())
  if(is.directed(nw())){
    deg <- sna::degree(bernoullisamples(), g=1:reps, gmode='digraph',
                  cmode=input$cmode_dd)
  } else {
    deg <- sna::degree(bernoullisamples(), g=1:reps, gmode='graph',
                  cmode=input$cmode_dd)
  }
  #now deg is a matrix where each element is a degree of a node
  #each column is a different draw

  degreedata <- apply(deg, MARGIN=2, FUN=tabulate, nbins=max(deg))
  #degreedata is matrix holding tabulation of degrees (except isolates)
  #each column is different draw
  z <- apply(deg, MARGIN=2, FUN=function(x){sum(x==0)})
  #tabulation of isolates in each draw
  degreedata <- matrix(data=c(z,t(degreedata)),nrow=max(deg)+1,ncol=reps,
                       byrow=TRUE)
  #complete tabulation of degrees for each draw
  degreemeans <- apply(degreedata, MARGIN=1, FUN=mean)
  names(degreemeans) <- paste(0:max(deg))
  degreesd <- apply(degreedata, MARGIN=1, FUN=sd)
  mean_and_sd <- list(degreemeans, degreesd)
})


observeEvent(input$percButton_dd, {
  state$plotperc_dd <- TRUE
})
observeEvent(input$countButton_dd, {
  state$plotperc_dd <- FALSE
})

output$degreedist <- renderPlot({
  if(!is.network(nw())){
    return()
  }
  input$plottabs
  input$rawdatafile
  input$samplenet

  plotme <- dd_plotdata()
  color <- histblue
  ylabel <- "Count of Nodes"
  xlabel <- "Degree"
  ltext <- c()
  lcol <- c() #color for lines
  lty <- c()
  lpch <- c()
  lfill <- c() #color for boxes
  lborder <- c()
  ltitle <- NULL

  if(input$cmode_dd == "indegree"){
    xlabel <- "In Degree"
  } else if (input$cmode_dd == "outdegree"){
    xlabel <- "Out Degree"
  }

  if(!is.null(input$colorby_dd)){
    if(input$colorby_dd != "None"){
      ncolors <- dim(dd_plotdata())[1]
      if(ncolors == 2){
        color <- c("#eff3ff", "#377FBC")
      } else if(ncolors < 10){
        color <- brewer.pal(ncolors,"Blues")
      } else if(ncolors >= 10){
        color <- colorRampPalette(brewer.pal(9,"Blues"))(ncolors)
      }
      ltext <- sort(unique(get.vertex.attribute(nw(),input$colorby_dd)))
      ltext <- append(ltext, "")
      lfill <- c(color, 0)
      lborder <- append(lborder, c(rep("black", times=ncolors), 0))
      lty <- rep(0, times=ncolors+1)
      lpch <- rep(26, times=ncolors+1)
      ltitle <- input$colorby_dd
  }}

  unif_samplemeans <- dd_uniformoverlay()[[1]]
  unif_stderr <- dd_uniformoverlay()[[2]]
  unif_upperline <- unif_samplemeans + 2*unif_stderr
  unif_lowerline <- unif_samplemeans - 2*unif_stderr
  maxdeg_u <- length(unif_samplemeans)-1

  bern_samplemeans <- dd_bernoullioverlay()[[1]]
  bern_stderr <- dd_bernoullioverlay()[[2]]
  bern_upperline <- bern_samplemeans + 2*bern_stderr
  bern_lowerline <- bern_samplemeans - 2*bern_stderr
  maxdeg_b <- length(bern_samplemeans)-1

  # get maximums for y limits of plot
  if(class(dd_plotdata())=="matrix"){
    maxfreq <- max(colSums(dd_plotdata()))
    maxdeg_obs <- dim(dd_plotdata())[2]-1
  } else {
    maxfreq <- max(dd_plotdata())
    maxdeg_obs <- length(dd_plotdata())-1
  }
  maxfreq_samples <- max(max(bern_upperline), max(unif_upperline))
  ylimit <- max(maxfreq, maxfreq_samples)

  if(state$plotperc_dd) {
    plotme <- dd_plotdata()/sum(dd_plotdata())
    unif_samplemeans <- unif_samplemeans/sum(dd_plotdata())
    unif_upperline <- unif_upperline/sum(dd_plotdata())
    unif_lowerline <- unif_lowerline/sum(dd_plotdata())
    bern_samplemeans <- bern_samplemeans/sum(dd_plotdata())
    bern_upperline <- bern_upperline/sum(dd_plotdata())
    bern_lowerline <- bern_lowerline/sum(dd_plotdata())
    ylimit <- max(maxfreq/sum(dd_plotdata()), max(unif_upperline),
                  max(bern_upperline))
    ylabel <- 'Percent of Nodes'
  }

  # make sure that barplot and lines have the same length
  maxdeg_total <- max(maxdeg_obs, maxdeg_u, maxdeg_b)
  if(maxdeg_u < maxdeg_total){
    unif_samplemeans <- append(unif_samplemeans, rep(0, times=maxdeg_total-maxdeg_u))
    unif_upperline <- append(unif_upperline, rep(0, times=maxdeg_total-maxdeg_u))
    unif_lowerline <- append(unif_lowerline, rep(0, times=maxdeg_total-maxdeg_u))
  }
  if(maxdeg_b < maxdeg_total){
    bern_samplemeans <- append(bern_samplemeans, rep(0, times=maxdeg_total-maxdeg_b))
    bern_upperline <- append(bern_upperline, rep(0, times=maxdeg_total-maxdeg_b))
    bern_lowerline <- append(bern_lowerline, rep(0, times=maxdeg_total-maxdeg_b))
  }
  if(maxdeg_obs < maxdeg_total){
    if(class(plotme)=="matrix"){
      nrows <- dim(plotme)[1]
      plotme <- cbind(plotme, matrix(0, nrow=nrows, ncol= maxdeg_total-maxdeg_obs))
      colnames(plotme) <- paste(0:maxdeg_total)
    } else {
      plotme <- append(plotme, rep(0,times=maxdeg_total-maxdeg_obs))
      names(plotme) <- paste(0:maxdeg_total)
    }
  }

  #save x-coordinates of bars, so that points are centered on bars
  bar_axis <- barplot(plotme, xlab=xlabel, ylab=ylabel,
                      col=color, ylim=c(0,ylimit), plot=TRUE)
  if(input$uniformoverlay_dd){
    points(x=bar_axis-.15, y=unif_samplemeans,col=CUGcol,
           lwd=1, pch=18, cex=1.25)
    suppressWarnings(arrows(x0=bar_axis-.15, y0=unif_upperline,
                            x1=bar_axis-.15, y1=unif_lowerline,
           code=3, length=0.1, angle=90, col=CUGcol))
    ltext <- append(ltext, "CUG")
    lcol <- append(lcol, CUGcol)
    lty <- append(lty, 1)
    lpch <- append(lpch, 18)
    lfill <- append(lfill, 0)
    lborder <- append(lborder, 0)
  }
  if(input$bernoullioverlay_dd){
    points(x = bar_axis+.15, y = bern_samplemeans, col = BRGcol,
           lwd = 1, pch = 18, cex = 1.25)
    suppressWarnings(arrows(x0 = bar_axis+.15, y0 = bern_upperline,
                            x1 = bar_axis+.15, y1 = bern_lowerline,
           code = 3, length = 0.1, angle = 90, col = BRGcol))
    ltext <- append(ltext, "BRG")
    lcol <- append(lcol, BRGcol)
    lty <- append(lty, 1)
    lpch <- append(lpch, 18)
    lfill <- append(lfill, 0)
    lborder <- append(lborder, 0)
  }
  if(input$colorby_dd != "None" | input$uniformoverlay_dd | input$bernoullioverlay_dd){
    if(input$uniformoverlay_dd | input$bernoullioverlay_dd){
      lmerge <- TRUE
    } else {
      lmerge <-FALSE
      lpch <-NULL
    }
    legend(x="topright", legend=ltext, title=ltitle, fill=lfill, border=lborder,
           col=lcol, lty= lty, pch=lpch, pt.cex=1.25, bty="n", merge=lmerge)
  }
})

output$degreedistdownload <- downloadHandler(
  filename = function(){paste(nwname(),'_degreedist.pdf',sep='')},
  content = function(file){
    pdf(file=file, height=8, width=12)
    plotme <- dd_plotdata()
    color <- histblue
    ylabel <- "Count of Nodes"
    ltext <- c()
    lcol <- c() #color for lines
    lty <- c()
    lpch <- c()
    lfill <- c() #color for boxes
    lborder <- c()
    ltitle <- NULL
    if(!is.null(input$colorby_dd)){
      if(input$colorby_dd != "None"){
        ncolors <- dim(dd_plotdata())[1]
        color <- brewer.pal(ncolors,"Blues")[1:ncolors]
        color[is.na(color)] <- brewer.pal(9, "Blues")
        ltext <- sort(unique(get.vertex.attribute(nw(),input$colorby_dd)))
        ltext <- append(ltext, "")
        lfill <- c(color, 0)
        lborder <- append(lborder, c(rep("black", times=ncolors), 0))
        lty <- rep(0, times=ncolors+1)
        lpch <- rep(26, times=ncolors+1)
        ltitle <- input$colorby_dd
      }}

    unif_samplemeans <- dd_uniformoverlay()[[1]]
    unif_stderr <- dd_uniformoverlay()[[2]]
    unif_upperline <- unif_samplemeans + 2*unif_stderr
    unif_lowerline <- unif_samplemeans - 2*unif_stderr
    maxdeg_u <- length(unif_samplemeans)-1

    bern_samplemeans <- dd_bernoullioverlay()[[1]]
    bern_stderr <- dd_bernoullioverlay()[[2]]
    bern_upperline <- bern_samplemeans + 2*bern_stderr
    bern_lowerline <- bern_samplemeans - 2*bern_stderr
    maxdeg_b <- length(bern_samplemeans)-1

    # get maximums for y limits of plot
    if(class(dd_plotdata())=="matrix"){
      maxfreq <- max(colSums(dd_plotdata()))
      maxdeg_obs <- dim(dd_plotdata())[2]-1
    } else {
      maxfreq <- max(dd_plotdata())
      maxdeg_obs <- length(dd_plotdata())-1
    }
    maxfreq_samples <- max(max(bern_upperline), max(unif_upperline))
    ylimit <- max(maxfreq, maxfreq_samples)

    if(state$plotperc_dd) {
      plotme <- dd_plotdata()/sum(dd_plotdata())
      unif_samplemeans <- unif_samplemeans/sum(dd_plotdata())
      unif_upperline <- unif_upperline/sum(dd_plotdata())
      unif_lowerline <- unif_lowerline/sum(dd_plotdata())
      bern_samplemeans <- bern_samplemeans/sum(dd_plotdata())
      bern_upperline <- bern_upperline/sum(dd_plotdata())
      bern_lowerline <- bern_lowerline/sum(dd_plotdata())
      ylimit <- max(maxfreq/sum(dd_plotdata()), max(unif_upperline),
                    max(bern_upperline))
      ylabel <- 'Percent of Nodes'
    }

    # make sure that barplot and lines have the same length
    maxdeg_total <- max(maxdeg_obs, maxdeg_u, maxdeg_b)
    if(maxdeg_u < maxdeg_total){
      unif_samplemeans <- append(unif_samplemeans, rep(0, times=maxdeg_total-maxdeg_u))
      unif_upperline <- append(unif_upperline, rep(0, times=maxdeg_total-maxdeg_u))
      unif_lowerline <- append(unif_lowerline, rep(0, times=maxdeg_total-maxdeg_u))
    }
    if(maxdeg_b < maxdeg_total){
      bern_samplemeans <- append(bern_samplemeans, rep(0, times=maxdeg_total-maxdeg_b))
      bern_upperline <- append(bern_upperline, rep(0, times=maxdeg_total-maxdeg_b))
      bern_lowerline <- append(bern_lowerline, rep(0, times=maxdeg_total-maxdeg_b))
    }
    if(maxdeg_obs < maxdeg_total){
      if(class(plotme)=="matrix"){
        nrows <- dim(plotme)[1]
        plotme <- cbind(plotme, matrix(0, nrow=nrows, ncol= maxdeg_total-maxdeg_obs))
        colnames(plotme) <- paste(0:maxdeg_total)
      } else {
        plotme <- append(plotme, rep(0,times=maxdeg_total-maxdeg_obs))
        names(plotme) <- paste(0:maxdeg_total)
      }
    }

    #save x-coordinates of bars, so that points are centered on bars
    bar_axis <- barplot(plotme, xlab="Degree", ylab=ylabel,
                        col=color, ylim=c(0,ylimit), plot=TRUE)
    if(input$uniformoverlay_dd){
      points(x=bar_axis-.15, y=unif_samplemeans,col=CUGcol,
             lwd=1, pch=18, cex=1.25)
      suppressWarnings(arrows(x0=bar_axis-.15, y0=unif_upperline,
                              x1=bar_axis-.15, y1=unif_lowerline,
             code=3, length=0.1, angle=90, col=CUGcol))
      ltext <- append(ltext, "CUG")
      lcol <- append(lcol, CUGcol)
      lty <- append(lty, 1)
      lpch <- append(lpch, 18)
      lfill <- append(lfill, 0)
      lborder <- append(lborder, 0)
    }
    if(input$bernoullioverlay_dd){
      points(x=bar_axis+.15, y=bern_samplemeans,col=BRGcol,
             lwd=1, pch=18, cex=1.25)
      suppressWarnings(arrows(x0=bar_axis+.15, y0=bern_upperline,
                              x1=bar_axis+.15, y1=bern_lowerline,
             code=3, length=0.1, angle=90, col=BRGcol))
      ltext <- append(ltext, "BRG")
      lcol <- append(lcol, BRGcol)
      lty <- append(lty, 1)
      lpch <- append(lpch, 18)
      lfill <- append(lfill, 0)
      lborder <- append(lborder, 0)
    }
    if(input$colorby_dd != "None" | input$uniformoverlay_dd | input$bernoullioverlay_dd){
      if(input$uniformoverlay_dd | input$bernoullioverlay_dd){
        lmerge <- TRUE
      } else {
        lmerge <-FALSE
        lpch <-NULL
      }
      legend(x="topright", legend=ltext, title=ltitle, fill=lfill, border=lborder,
             col=lcol, lty= lty, pch=lpch, pt.cex=1.25, bty="n", merge=lmerge)
    }
    dev.off()
})

#GEODESIC DISTRIBUTION


observeEvent(input$percButton_gd, {
  state$plotperc_gd <- TRUE
})
observeEvent(input$countButton_gd, {
  state$plotperc_gd <- FALSE
})

gdata <- reactive({
  g <- sna::geodist(nw(), inf.replace=NA, count.paths=FALSE)
  gdata <- tabulate(g$gdist)
  g$gdist[is.na(g$gdist)] <- Inf
  gdata <- append(gdata, sum(g$gdist == Inf))
  maxgeo <- length(gdata)-1
  names(gdata) <- c(paste(1:maxgeo), "Inf")
  gdata
})

gd_uniformoverlay <- reactive({
  if(!is.network(nw())){
    return()
  }
  gd <- sna::geodist(uniformsamples(), count.paths=FALSE, inf.replace=NA)
  maxgeo <- max(unlist(gd), na.rm=TRUE)
  reps <- length(gd)
  gd_data <- lapply(gd,function(x){tabulate(x$gdist, nbins=maxgeo)})
    #list of tabulated geodesics for each draw, except those of 0 or Inf length
  gd_data_complete <- matrix(0, nrow=maxgeo+1, ncol=reps)
  for(k in 1:reps){
    temp <- append(gd_data[[k]], sum(is.na(gd[[k]]$gdist)))
    gd_data_complete[,k] <- temp
  }

  geomeans <- apply(gd_data_complete, MARGIN=1, FUN=mean)
  names(geomeans) <- c(paste(1:maxgeo), "Inf")
  geosd <- apply(gd_data_complete, MARGIN=1, FUN=sd)

  mean_and_sd <- list(geomeans, geosd)
})

gd_bernoullioverlay <- reactive({
  if(!is.network(nw())){
    return()
  }
  gd <- sna::geodist(bernoullisamples(), count.paths=FALSE, inf.replace=NA)
  maxgeo <- max(unlist(gd), na.rm=TRUE)
  reps <- length(gd)
  gd_data <- lapply(gd,function(x){tabulate(x$gdist, nbins=maxgeo)})
  #list of tabulated geodesics for each draw, except those of 0 or Inf length
  gd_data_complete <- matrix(0, nrow=maxgeo+1, ncol=reps)
  for(k in 1:reps){
    temp <- append(gd_data[[k]], sum(is.na(gd[[k]]$gdist)))
    gd_data_complete[,k] <- temp
  }

  geomeans <- apply(gd_data_complete, MARGIN=1, FUN=mean)
  names(geomeans) <- c(paste(1:maxgeo), "Inf")
  geosd <- apply(gd_data_complete, MARGIN=1, FUN=sd)

  mean_and_sd <- list(geomeans, geosd)
})

output$geodistplot <- renderPlot({
  if(!is.network(nw())){
    return()
  }
  input$plottabs
  input$rawdatafile
  input$samplenet

  gdata <- gdata()
  maxgeo <- length(gdata)-1

  unif_means <- gd_uniformoverlay()[[1]]
  unif_stderr <- gd_uniformoverlay()[[2]]
  unif_upperline <- unif_means + 2*unif_stderr
  unif_lowerline <- unif_means - 2*unif_stderr
  maxgeo_u <- length(unif_means)-1

  bern_means <- gd_bernoullioverlay()[[1]]
  bern_stderr <- gd_bernoullioverlay()[[2]]
  bern_upperline <- bern_means + 2*bern_stderr
  bern_lowerline <- bern_means - 2*bern_stderr
  maxgeo_b <- length(bern_means)-1

  ylabel <- "Count of Vertex Pairs"

  #for density plot
  if(state$plotperc_gd){
    unif_means <- unif_means/sum(gdata)
    unif_upperline <- unif_upperline/sum(gdata)
    unif_lowerline <- unif_lowerline/sum(gdata)
    bern_means <- bern_means/sum(gdata)
    bern_upperline <- bern_upperline/sum(gdata)
    bern_lowerline <- bern_lowerline/sum(gdata)
    gdata <- gdata/sum(gdata)
    ylimit <- max(gdata, bern_upperline, unif_upperline)
    ylabel <- "Percent of Vertex Pairs"
  }

  # make sure that barplot and lines have the same length
  maxgeo_total <- max(maxgeo, maxgeo_u, maxgeo_b)
  if(maxgeo_u < maxgeo_total){
    unif_means <- append(unif_means, rep(0, times=maxgeo_total-maxgeo_u),
                         after=length(unif_means)-1)
    unif_upperline <- append(unif_upperline, rep(0, times=maxgeo_total-maxgeo_u),
                             after=length(unif_upperline)-1)
    unif_lowerline <- append(unif_lowerline, rep(0, times=maxgeo_total-maxgeo_u),
                             after=length(unif_lowerline)-1)
  }
  if(maxgeo_b < maxgeo_total){
    bern_means <- append(bern_means, rep(0, times=maxgeo_total-maxgeo_b),
                         after=length(bern_means)-1)
    bern_upperline <- append(bern_upperline, rep(0, times=maxgeo_total-maxgeo_b),
                             after=length(bern_upperline)-1)
    bern_lowerline <- append(bern_lowerline, rep(0, times=maxgeo_total-maxgeo_b),
                             after=length(bern_lowerline)-1)
  }
  if(maxgeo < maxgeo_total){
    gdata <- append(gdata, rep(0,times=maxgeo_total-maxgeo),
                    after=length(gdata)-1)
    names(gdata) <- c(paste(1:maxgeo_total), "Inf")
  }

  if(input$excludeInfs){
    gdata <- gdata[1:maxgeo_total]
    bern_means <- bern_means[1:maxgeo_total]
    bern_upperline <- bern_upperline[1:maxgeo_total]
    bern_lowerline <- bern_lowerline[1:maxgeo_total]
    unif_means <- unif_means[1:maxgeo_total]
    unif_upperline <- unif_upperline[1:maxgeo_total]
    unif_lowerline <- unif_lowerline[1:maxgeo_total]
  }
  #get maximums to set y limits
  ylimit <- max(gdata, bern_upperline, unif_upperline)

  ltext <- c()
  lcol <- c()

  #save x-coordinates of bars, so that points are centered on bars
  bar_axis <- barplot(gdata,  col=histblue,
                      xlab = "Geodesic Value", ylab = ylabel,
                      ylim = c(0,ylimit), plot=TRUE)

  if(input$uniformoverlay_gd){
    points(x=bar_axis-.15, y=unif_means,col=CUGcol,
           lwd=1, pch=18, cex=1.25)
    suppressWarnings(arrows(x0=bar_axis-.15, y0=unif_upperline,
                            x1=bar_axis-.15, y1=unif_lowerline,
           code=3, length=0.1, angle=90, col=CUGcol))
    ltext <- append(ltext, "CUG")
    lcol <- append(lcol, CUGcol)
  }
  if(input$bernoullioverlay_gd){
    points(x=bar_axis+.15, y=bern_means,col=BRGcol,
           lwd=1, pch=18, cex=1.25)
    suppressWarnings(arrows(x0=bar_axis+.15, y0=bern_upperline,
                            x1=bar_axis+.15, y1=bern_lowerline,
           code=3, length=0.1, angle=90, col=BRGcol))
    ltext <- append(ltext, "BRG")
    lcol <- append(lcol, BRGcol)
  }
  if(input$uniformoverlay_gd | input$bernoullioverlay_gd){
    legend(x="topright", legend=ltext, col=lcol, lwd=1, pch=18,
           pt.cex=1.25, merge=TRUE,
           inset=c(.12,0), bty="n")
  }

})

output$geodistdownload <- downloadHandler(
  filename = function(){paste(nwname(),'_geodist.pdf',sep='')},
  content = function(file){
    pdf(file=file, height=8, width=12)

    gdata <- gdata()
    maxgeo <- length(gdata)-1

    unif_means <- gd_uniformoverlay()[[1]]
    unif_stderr <- gd_uniformoverlay()[[2]]
    unif_upperline <- unif_means + 2*unif_stderr
    unif_lowerline <- unif_means - 2*unif_stderr
    maxgeo_u <- length(unif_means)-1

    bern_means <- gd_bernoullioverlay()[[1]]
    bern_stderr <- gd_bernoullioverlay()[[2]]
    bern_upperline <- bern_means + 2*bern_stderr
    bern_lowerline <- bern_means - 2*bern_stderr
    maxgeo_b <- length(bern_means)-1

    ylabel <- "Count of Vertex Pairs"

    #for density plot
    if(state$plotperc_gd){
      unif_means <- unif_means/sum(gdata)
      unif_upperline <- unif_upperline/sum(gdata)
      unif_lowerline <- unif_lowerline/sum(gdata)
      bern_means <- bern_means/sum(gdata)
      bern_upperline <- bern_upperline/sum(gdata)
      bern_lowerline <- bern_lowerline/sum(gdata)
      gdata <- gdata/sum(gdata)
      ylimit <- max(gdata, bern_upperline, unif_upperline)
      ylabel <- "Percent of Vertex Pairs"
    }

    # make sure that barplot and lines have the same length
    maxgeo_total <- max(maxgeo, maxgeo_u, maxgeo_b)
    if(maxgeo_u < maxgeo_total){
      unif_means <- append(unif_means, rep(0, times=maxgeo_total-maxgeo_u),
                           after=length(unif_means)-1)
      unif_upperline <- append(unif_upperline, rep(0, times=maxgeo_total-maxgeo_u),
                               after=length(unif_upperline)-1)
      unif_lowerline <- append(unif_lowerline, rep(0, times=maxgeo_total-maxgeo_u),
                               after=length(unif_lowerline)-1)
    }
    if(maxgeo_b < maxgeo_total){
      bern_means <- append(bern_means, rep(0, times=maxgeo_total-maxgeo_b),
                           after=length(bern_means)-1)
      bern_upperline <- append(bern_upperline, rep(0, times=maxgeo_total-maxgeo_b),
                               after=length(bern_upperline)-1)
      bern_lowerline <- append(bern_lowerline, rep(0, times=maxgeo_total-maxgeo_b),
                               after=length(bern_lowerline)-1)
    }
    if(maxgeo < maxgeo_total){
      gdata <- append(gdata, rep(0,times=maxgeo_total-maxgeo), after=length(gdata)-1)
      names(gdata) <- c(paste(1:maxgeo_total), "INF")
    }

    if(input$excludeInfs){
      gdata <- gdata[1:maxgeo_total]
      bern_means <- bern_means[1:maxgeo_total]
      bern_upperline <- bern_upperline[1:maxgeo_total]
      bern_lowerline <- bern_lowerline[1:maxgeo_total]
      unif_means <- unif_means[1:maxgeo_total]
      unif_upperline <- unif_upperline[1:maxgeo_total]
      unif_lowerline <- unif_lowerline[1:maxgeo_total]
    }

    #get maximums to set y limits
    ylimit <- max(gdata, bern_upperline, unif_upperline)

    ltext <- c()
    lcol <- c()

    #save x-coordinates of bars, so that points are centered on bars
    bar_axis <- barplot(gdata,  col=histblue,
                        xlab = "Geodesic Value", ylab = ylabel,
                        ylim = c(0,ylimit), plot=TRUE)

    if(input$uniformoverlay_gd){
      points(x=bar_axis-.15, y=unif_means,col=CUGcol,
             lwd=1, pch=18, cex=1.25)
      suppressWarnings(arrows(x0=bar_axis-.15, y0=unif_upperline,
                              x1=bar_axis-.15, y1=unif_lowerline,
             code=3, length=0.1, angle=90, col=CUGcol))
      ltext <- append(ltext, "CUG")
      lcol <- append(lcol, CUGcol)
    }
    if(input$bernoullioverlay_gd){
      points(x=bar_axis+.15, y=bern_means,col=BRGcol, lwd=1,
             pch=18, cex=1.25)
      suppressWarnings(arrows(x0=bar_axis+.15, y0=bern_upperline,
                              x1=bar_axis+.15, y1=bern_lowerline,
             code=3, length=0.1, angle=90, col=BRGcol))
      ltext <- append(ltext, "BRG")
      lcol <- append(lcol, BRGcol)
    }
    if(input$uniformoverlay_gd | input$bernoullioverlay_gd){
      legend(x="topright", legend=ltext, col=lcol, lwd=1,
             pch=18, pt.cex=1.25, merge=TRUE,
             inset=c(.12,0), bty="n")
    }
    dev.off()
  })

  output$infsummary <- renderPrint({
    options(digits=3)
    gdata <- gdata()
    cols <- length(gdata)
    obs <- gdata[cols]
    unif_means <- gd_uniformoverlay()[[1]]
    unif_sd <- gd_uniformoverlay()[[2]]
    bern_means <- gd_bernoullioverlay()[[1]]
    bern_sd <- gd_bernoullioverlay()[[2]]

#     v<-data.frame(Obs=paste(obs), CUG=paste(unif_means[cols],"(",round(unif_sd[cols], digits=2),")",sep=""),
#                BRG=paste(bern_means[cols],"(",round(bern_sd[cols], digits=2),")", sep=""),
#                row.names="Infs:")
#     format(v, justify="centre")
    v <- c(paste0(unif_means[cols],"(",round(unif_sd[cols], digits=2),")"),
           paste0(bern_means[cols],"(",round(bern_sd[cols], digits=2),")"))
    cat(format("",width=4),format("Observed",width=8),
        format(c("CUG","BRG"),width=11,justify="centre"),"\n")
    cat(format("Infs:",width=3,justify="centre"),
        format(paste(obs),width=8,justify="centre"),
        format(v,width=10,justify="centre"))
  })

#MORE

observe({
  if(input$plottabs == "More"){
    updateTabsetPanel(session, 'displaytabs', selected="Network Summary")
  }
})

output$dynamiccugterm <- renderUI({
  if(!is.network(nw())){return()}
  if(is.directed(nw())){
    choices <- c("density", "isolates", "mean degree" = "meandeg", "mutual",
                 "transitive triads" = "transitive", "triangle", "twopath")
  } else {
    choices <- c("density", "concurrent", "isolates", "mean degree" = "meandeg",
                 "triangle")
  }
  #matchingterms <- splitargs(nw = nw())
  #choices <- matchingterms$names[matchingterms$args == "()"]
  selectizeInput("cugtestterm", label = "Model term",
              choices = choices)
})
outputOptions(output, 'dynamiccugterm', suspendWhenHidden = FALSE)

observeEvent(input$cugButton, {
  state$cugplot <- 1
})
observeEvent(nw(), {
  state$cugplot <- 0
})

output$cugtest <- renderPlot({
  if(!is.network(nw())) {return()}
  if(input$cugButton == 0) {return()}
  if(state$cugplot == 0){return()}

  isolate({
    term <- input$cugtestterm
    n <- nodes()
    obsval <- summary.formula(as.formula(paste("nw() ~", term)))

    # gets summary statistics of the already run simulations
    brgvals <- brgvals()
    cugvals <- cugvals()

    brghist <- hist(brgvals, plot = FALSE)

    getbreaks <- function(x){
      breaks <- brghist$breaks
      r <- range(brghist$breaks)
      int <- breaks[2] - breaks[1]
      if(min(x) < r[1]){
        toadd <- ceiling((r[1] - min(x))/int)
        front <- c((r[1]-toadd*int):(r[1]-int))
        breaks <- c(front, breaks)
      }
      if (max(x) > r[2]) {
        toadd <- ceiling((max(x)-r[2])/int)
        back <- c((r[2]+int):(r[2]+toadd*int))
        breaks <- c(breaks, back)
      }
      return(breaks)
    }

    xlims <- c(min(brgvals, cugvals, obsval) - diff(brghist$breaks)[1],
               max(brgvals, cugvals, obsval) + diff(brghist$breaks)[1])

    cughist <- hist(cugvals, breaks = getbreaks, plot = FALSE)
    par(lwd = 2)
    hist(brgvals, col = tgray3,
         border = BRGcol, ylab = NULL, main = NULL, xlab= NULL,
         xaxt = "n",
         xlim = xlims,
         ylim = c(0, max(brghist$counts, cughist$counts)),
         breaks = brghist$breaks)

    if (term == "density" | term == "meandeg"){
      abline(v = cugvals[1], col = CUGcol)
    } else {
      hist(cugvals, col = tgray7, density = 15, angle = -45,
           border = CUGcol, ylab = NULL, main = NULL, xlab= NULL,
           axes = FALSE,
           breaks = getbreaks, add = TRUE)

    }

    axis(side = 1, at = round(c(xlims[1], cughist$breaks, xlims[2]), digits = 3))
    points(x = obsval, y = 0, col = obsblue, pch = 17, cex = 2)

    legend(x = "topright", bty = "n",
           legend = c("Observed value", "CUG distribution", "BRG distribution"),
           pch = c(17, NA, NA), col = obsblue, fill = c(0, tgray7, tgray3),
           angle = -45, density = c(0, 15, 100),
           border = c(0, CUGcol, BRGcol))
    par(lwd = 1)
  })
})

output$cugtestdownload <- downloadHandler(
  filename = function(){paste0(nwname(), "_", input$cugtestterm, ".pdf")},
  content = function(file){

    term <- input$cugtestterm
    n <- nodes()
    obsval <- summary.formula(as.formula(paste("nw() ~", term)))

    # gets summary statistics of the already run simulations
    brgvals <- brgvals()
    cugvals <- cugvals()

    brghist <- hist(brgvals, plot = FALSE)

    getbreaks <- function(x){
      breaks <- brghist$breaks
      r <- range(brghist$breaks)
      int <- breaks[2] - breaks[1]
      if(min(x) < r[1]){
        toadd <- ceiling((r[1] - min(x))/int)
        front <- c((r[1]-toadd*int):(r[1]-int))
        breaks <- c(front, breaks)
      }
      if (max(x) > r[2]) {
        toadd <- ceiling((max(x)-r[2])/int)
        back <- c((r[2]+int):(r[2]+toadd*int))
        breaks <- c(breaks, back)
      }
      return(breaks)
    }

    xlims <- c(min(brgvals, cugvals, obsval) - diff(brghist$breaks)[1],
               max(brgvals, cugvals, obsval) + diff(brghist$breaks)[1])

    cughist <- hist(cugvals, breaks = getbreaks, plot = FALSE)

    pdf(file = file)

    par(lwd = 2)
    hist(brgvals, col = tgray3,
         border = BRGcol, ylab = NULL, main = NULL, xlab= NULL,
         xaxt = "n",
         xlim = xlims,
         ylim = c(0, max(brghist$counts, cughist$counts)),
         breaks = brghist$breaks)

    if (input$cugtestterm == "density" | input$cugtestterm == "meandeg"){
      abline(v = cugvals[1], col = CUGcol)
    } else {
      hist(cugvals, col = tgray7, density = 15, angle = -45,
           border = CUGcol, ylab = NULL, main = NULL, xlab= NULL,
           axes = FALSE,
           breaks = getbreaks, add = TRUE)

    }

    axis(side = 1, at = round(c(xlims[1], cughist$breaks, xlims[2]), digits = 3))
    points(x = obsval, y = 0, col = obsblue, pch = 17, cex = 2)

    legend(x = "topright", bty = "n",
           legend = c("Observed value", "CUG distribution", "BRG distribution"),
           pch = c(17, NA, NA), col = obsblue, fill = c(0, tgray7, tgray3),
           angle = -45, density = c(0, 15, 100),
           border = c(0, CUGcol, BRGcol))
    par(lwd = 1)
    dev.off()
  }
)

#since the visibility toggles between two states, set the options to
#not suspend the output when hidden
output$mixmxchooser <- renderUI({
  selectInput('mixmx', label='Choose attribute',
              choices = menuattr())
})
outputOptions(output,'mixmxchooser',suspendWhenHidden=FALSE)

output$mixingmatrix <- renderPrint({
  if(!is.network(nw())) {return()}
  if(!is.null(input$mixmx)){
  mixingmatrix(nw(), input$mixmx)}
})
outputOptions(output,'mixingmatrix',suspendWhenHidden=FALSE)

output$mixmxdownload <- downloadHandler(
  filename = function() {paste0(nwname(), "_mixingmatrix.csv")},
  contentType = "text/csv",
  content = function(file) {
    mx <- mixingmatrix(nw(), input$mixmx)[["matrix"]]
    write.csv(mx, file = file)
  }
)

# update all the menu selection options for descriptive indices when network changes
observeEvent(nw(), {
  if(is.directed(nw())){
    degmenu <- c('indegree', 'outdegree')
    betwmenu <- c('directed', 'endpoints', 'proximalsrc',
                  'proximaltar', 'proximalsum', 'lengthscaled', 'linearscaled')
    closemenu <- c('directed', 'suminvdir')
    stressmenu <- c('directed')
    hgmenu <- c('directed')
  } else {
    degmenu <- c('total')
    betwmenu <- c('undirected', 'endpoints', 'proximalsrc',
                  'proximaltar', 'proximalsum', 'lengthscaled', 'linearscaled')
    closemenu <- c('undirected', 'suminvundir')
    stressmenu <- c('undirected')
    hgmenu <- c('undirected')
  }
  updateNumericInput(session, "nodeind", label = NULL, value = 1,
                     min = 1, max = nodes())
  updateSelectInput(session, "gdegcmode", choices = degmenu)
  updateSelectInput(session, "gbetwcmode", choices = betwmenu)
  updateSelectInput(session, "gclosecmode", choices = closemenu)
  updateSelectInput(session, "gstresscmode", choices = stressmenu)
  updateSelectInput(session, "ggraphcentcmode", choices = hgmenu)
  updateSelectInput(session, "ndegcmode", choices = degmenu)
  updateSelectInput(session, "nbetwcmode", choices = betwmenu)
  updateSelectInput(session, "nclosecmode", choices = closemenu)
  updateSelectInput(session, "nstresscmode", choices = stressmenu)
  updateSelectInput(session, "ngraphcentcmode", choices = hgmenu)
})

output$gden <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  sna::gden(nw(), diag=has.loops(nw()), mode=gmode)
})
outputOptions(output,'gden',suspendWhenHidden=FALSE)

output$grecip <- renderText({
  if(!is.network(nw())) {return()}
  if(input$grecipmeas == ''){
    return()
  }
  try(sna::grecip(nw(), measure=input$grecipmeas))
})
outputOptions(output,'grecip',suspendWhenHidden=FALSE)

output$gtrans <- renderText({
  if(!is.network(nw())) {return()}
  if(input$gtransmeas == ''){
    return()
  }
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  try(sna::gtrans(nw(), diag=has.loops(nw()), mode=gmode,
                  measure=input$gtransmeas))
})
outputOptions(output,'gtrans',suspendWhenHidden=FALSE)

output$gdeg <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  d <- ""
  cmode <- input$gdegcmode
  if(cmode == 'total'){
    cmode <- 'freeman'
  }
  try(d <- sna::centralization(nw(), sna::degree, mode=gmode, diag=has.loops(nw()),
              cmode=cmode))
  d
})
outputOptions(output,'gdeg',suspendWhenHidden=FALSE)

output$gbetw <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  b <- ""
  try(b <- sna::centralization(nw(), sna::betweenness, mode=gmode, diag=has.loops(nw()),
                   cmode=input$gbetwcmode))
  b
})
outputOptions(output,'gbetw',suspendWhenHidden=FALSE)

output$gclose <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  c <- ""
  try(
    c <- sna::centralization(nw(), sna::closeness, mode=gmode, diag=has.loops(nw()),
                        cmode=input$gclosecmode))
  c
})
outputOptions(output,'gclose',suspendWhenHidden=FALSE)

output$gstress <- renderText({
  if(!is.network(nw())){ return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  s <- ""
  try(s <- sna::centralization(nw(), sna::stresscent, mode=gmode, diag=has.loops(nw()),
                          cmode=input$gstresscmode))
  s
})
outputOptions(output,'gstress',suspendWhenHidden=FALSE)

output$ggraphcent <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  g <- ""
  try(g <- sna::centralization(nw(), sna::graphcent, mode=gmode, diag=has.loops(nw()),
                          cmode=input$ggraphcentcmode))
  g
})
outputOptions(output,'ggraphcent',suspendWhenHidden=FALSE)

output$gevcent <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  e <- ""
  try(e <- sna::centralization(nw(), sna::evcent, mode=gmode, diag=has.loops(nw())))
  e
})
outputOptions(output,'gevcent',suspendWhenHidden=FALSE)

output$ginfocent <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  i<-""
  try({
    i <- sna::centralization(nw(), sna::infocent, mode=gmode, diag=has.loops(nw()),
                        cmode=input$ginfocentcmode)})
  i
})
outputOptions(output,'ginfocent',suspendWhenHidden=FALSE)

output$ndeg <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  cmode <- input$ndegcmode
  if(cmode == 'total'){
    cmode <- 'freeman'
  }
  d <- ""
  try(d <- sna::degree(nw(), nodes=input$nodeind, gmode=gmode, diag=has.loops(nw()),
                  cmode=cmode))
  d
})

output$ndegmin <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  d <- sna::degree(nw(), gmode=gmode, diag=has.loops(nw()),
              cmode=input$ndegcmode)
  min(d)
})

output$ndegmax <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  d <- sna::degree(nw(), gmode=gmode, diag=has.loops(nw()),
              cmode=input$ndegcmode)
  max(d)
})

output$nbetw <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  b <- ""
  try(b <- sna::betweenness(nw(), nodes=input$nodeind, gmode=gmode,
                       diag=has.loops(nw()),
                       cmode=input$nbetwcmode))
  b
})

output$nbetwmin <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  b <- sna::betweenness(nw(), gmode=gmode, diag=has.loops(nw()),
                   cmode=input$nbetwcmode)
  min(b)
})

output$nbetwmax <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  b <- sna::betweenness(nw(), gmode=gmode, diag=has.loops(nw()),
                   cmode=input$nbetwcmode)
  max(b)
})

output$nclose <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  c <- ""
  try(
    c <- sna::closeness(nw(), nodes=input$nodeind, gmode=gmode, diag=has.loops(nw()),
                   cmode=input$nclosecmode))
  c
})

output$nclosemin <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  c <- sna::closeness(nw(), gmode=gmode, diag=has.loops(nw()),
                 cmode=input$nclosecmode)
  min(c)
})

output$nclosemax <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  c <- sna::closeness(nw(), gmode=gmode, diag=has.loops(nw()),
                 cmode=input$nclosecmode)
  max(c)
})

output$nstress <- renderText({
  if(!is.network(nw())){ return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  s <- ""
  try(s <- sna::stresscent(nw(), nodes=input$nodeind, gmode=gmode,
                      diag=has.loops(nw()),
                      cmode=input$nstresscmode))
  s
})

output$nstressmin <- renderText({
  if(!is.network(nw())){ return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  s <- sna::stresscent(nw(), gmode=gmode, diag=has.loops(nw()),
                  cmode=input$nstresscmode)
  min(s)
})

output$nstressmax <- renderText({
  if(!is.network(nw())){ return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  s <- sna::stresscent(nw(), gmode=gmode, diag=has.loops(nw()),
                  cmode=input$nstresscmode)
  max(s)
})

output$ngraphcent <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  g <- ""
  try(g <- sna::graphcent(nw(), nodes=input$nodeind, gmode=gmode,
                          diag=has.loops(nw()),
                          cmode=input$ngraphcentcmode))
  g
})

output$ngraphcentmin <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  g <- sna::graphcent(nw(), gmode=gmode, diag=has.loops(nw()),
                      cmode=input$ngraphcentcmode)
  min(g)
})

output$ngraphcentmax <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  g <- sna::graphcent(nw(), gmode=gmode, diag=has.loops(nw()),
                      cmode=input$ngraphcentcmode)
  max(g)
})

output$nevcent <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  e <- ""
  try(e <- sna::evcent(nw(), nodes=input$nodeind, gmode=gmode, diag=has.loops(nw())))
  e
})

output$nevcentmin <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  e <- ""
  try({e <- sna::evcent(nw(), gmode=gmode, diag=has.loops(nw()))
       e <- min(e)})
  e
})

output$nevcentmax <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  e <- ""
  try({e <- sna::evcent(nw(), gmode=gmode, diag=has.loops(nw()))
       e <- max(e)})
  e
})

observeEvent(nw(), {
  state$err <- FALSE
  values$err <- c()
})

output$errstate <- renderPrint({
  print(1)
})

output$ninfocent <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  i <- ''
  i <- tryCatch({sna::infocent(nw(), nodes=input$nodeind, gmode=gmode,
                              diag=has.loops(nw()), cmode=input$ninfocentcmode)},
                error = function(e) {e})
  if("error" %in% class(i)) {
    state$err <- TRUE
    values$err[1] <- i[[1]]
    #i <- paste("Error", length(values$err))
    i <- ''
  }
  i
})
outputOptions(output, "ninfocent", priority = 50)


output$errbox <- renderPrint({
  err <- FALSE
  if(!is.null(state$err)){
    err <- state$err
    if(err){
      err <- paste0("Error ", c(1:length(values$err)), ": ", values$err, "\n")
    }
  }
  p(err)
})
outputOptions(output, "errbox", priority = 1)

output$ninfocentmin <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  i<-''
  try({
    i <- sna::infocent(nw(), gmode=gmode, diag=has.loops(nw()),
                      cmode=input$ninfocentcmode)
    i<-min(i)})
  i
})

output$ninfocentmax <- renderText({
  if(!is.network(nw())) {return()}
  if(is.directed(nw())){
    gmode <- 'digraph'
  } else {
    gmode <- 'graph'
  }
  i <- ''
  try({
    i <- sna::infocent(nw(), gmode=gmode, diag=has.loops(nw()),
                      cmode=input$ninfocentcmode)
    i<-max(i)})
  i
})


# Fit Model ---------------------------------------------------------------

observeEvent(input$matchingButton, {
  state$allterms <- FALSE
})
observeEvent(input$allButton, {
  state$allterms <- TRUE
})

output$listofterms <- renderUI({
  if(!is.network(nw())){
    return()
  }
  if(state$allterms){
    current.terms <- allterms$names
  } else {
    matchterms <- splitargs(nw = nw())
    current.terms <- matchterms$names
  }
  selectizeInput('chooseterm', label = NULL,
                 choices = c("Select a term" = "", current.terms))
})

output$termdoc <- renderUI({
  myterm <- input$chooseterm
  if(is.null(myterm)){
    return(p("Select or search for a term in the menu above."))
  } else if(myterm == ""){
    return(p("Select or search for a term in the menu above."))
  }
  chrvec <- capture.output(search.ergmTerms(name = myterm))
  desc <- strsplit(chrvec[3], split = "_")
  p(chrvec[1], br(),br(),
    strong(chrvec[2]), br(),br(),
    em(desc[[1]][2]), desc[[1]][3], br(),
    chrvec[4])
})

observe({
  if(input$controldefault){
    updateNumericInput(session, "MCMCinterval", value=1024)
    #burn-in gets updated separately, to be 16*interval
    updateNumericInput(session, "MCMCsamplesize", value=1024)
    disableWidget("MCMCinterval", session, disabled=TRUE)
    disableWidget("MCMCburnin", session, disabled=TRUE)
    disableWidget("MCMCsamplesize", session, disabled=TRUE)
    disableWidget("customMCMCcontrol", session, disable=TRUE)
  } else {
    disableWidget("MCMCinterval", session, disabled=FALSE)
    disableWidget("MCMCburnin", session, disabled=FALSE)
    disableWidget("MCMCsamplesize", session, disabled=FALSE)
    disableWidget("customMCMCcontrol", session, disable=FALSE)
  }
})

observe({
  input$MCMCinterval
  updateNumericInput(session, "burnin", value=16*input$MCMCinterval)
})

output$currentdataset1 <- renderPrint({
  if(!is.network(nw())){
    return(cat('Upload a network'))
  }
  cat(isolate(nwname()))
})

output$checkterms_fit <- renderPrint({
  if(!is.network(nw())){
    return(cat('Upload a network'))
  }
  if(ergm.terms()=='NA') return(cat('Add terms to the formula'))
  cat(ergm.terms())
})

output$prefitsum <- renderPrint({
  if(!is.network(nw()) | length(input$terms)==0){
    return(cat('NA'))
  }
  if(ergm.terms()=='NA') return(cat('Add terms to the formula'))
  options(width=140)
  summary(ergm.formula())
})

output$savemodel <- renderUI({
  m <- values$modeltotal
  actionButton('savemodelButton', label=paste0('Save Current Model (',m,'/5)'),
                 class="btn-sm")
})
outputOptions(output,'savemodel',suspendWhenHidden=FALSE)

observe({
  input$fitButton
  values$modelstate <- 1  #modelfit is up to date
})
observe({
  nw()
  input$clearmodelButton
  values$modelstate <- 0 #modelfit is outdated when nw changes or models cleared
})

output$modelfit <- renderPrint({
  if (input$fitButton == 0){
    return(cat('After adding terms to the formula, click "Fit Model" above.'))
  }
  if (values$modelstate == 0){
    return(cat('After adding terms to the formula, click "Fit Model" above.'))
  }
  model1reac()
})

output$modelfitsum <- renderPrint({
  if (input$fitButton == 0){
    return(cat('After adding terms to the formula, click "Fit Model" above.'))
  }
  if (values$modelstate == 0){
    return(cat('After adding terms to the formula, click "Fit Model" above.'))
  }
  options(width=140)
  summary(model1reac())
})

output$modelfitdownload <- downloadHandler(
  filename = function() {paste0(nwname(), "_modelfit.txt")},
  contentType = "text/csv",
  content = function(file) {
    capture.output(summary(model1reac()), file = file)
  }
)

output$modelcomparison <- renderPrint({
  options(width=140)
  x <- values$modelcoefs
  y <- values$modelsumstats
  if(length(x) == 0){return(cat(""))}
  coef.comparison(x)
  cat("\n\n","Summary Statistics", "\n")
  stat.comparison(y)
})

output$modelcompdownload <- downloadHandler(
  filename = function() {paste0(nwname(), "_modelcomparison.txt")},
  contentType = "text/csv",
  content = function(file) {
    x <- values$modelcoefs
    y <- values$modelsumstats
    capture.output(cat(nwname(),"\n"), coef.comparison(x),
                   cat("\n\n","Summary Statistics", "\n"),
                   stat.comparison(y),
                   file = file)
  }
)

# make sure that mcmc iterations output on the fitting tab by allowing
# modelfit to update even when other tab is active

# other potential methods were to use onFlush, or to set the priority of
# observers, but this is the best and least complicated (and the only
# one that works)

outputOptions(output, "modelfit", priority = 10, suspendWhenHidden = FALSE)
outputOptions(output, "modelfitsum", priority = -10)


# MCMC Diagnostics --------------------------------------------------------


# When using the `mcmc.diagnostics` function in the command line, the printed
# diagnostics and plots all output together. Instead of calling
# `mcmc.diagnositcs` in a reactive object, it gets called in both the plot
# output element and summary output element.


output$uichoosemodel_mcmc <- renderUI({
  n <- values$modeltotal
  if(n == 0){
    inlineSelectInput("choosemodel_mcmc",label=NULL,
                choices=c("Current"))
  } else {
    inlineSelectInput("choosemodel_mcmc",label=NULL,
                choices=c(paste0("Model",1:n)))
  }
})
outputOptions(output,"uichoosemodel_mcmc",suspendWhenHidden=FALSE)

output$checkterms_mcmc <- renderPrint({
  if(is.null(nw())){
    return(cat('Upload a network'))
  }
  if(input$fitButton == 0){
    return(cat('Please fit a model'))
  }
  mod <- input$choosemodel_mcmc
  if(mod=="Current"){
    cat(isolate(ergm.terms()))
  } else {
    mod <- as.numeric(substr(mod,6,6))
    cat(values$modelformulas[[mod]])
  }
})

output$currentdataset_mcmc <- renderPrint({
  if(!is.network(nw())){
    return(cat('Upload a network'))
  }
  cat(nwname())
})

output$diagnosticsplot <- renderPlot({
  if(ergm.terms()=="NA"){
    return()
  }
  mod <- input$choosemodel_mcmc
  if(mod=="Current"){
    mod <- model1reac()
  } else {
    modn <- as.numeric(substr(mod,6,6))
    mod <- values$modelfits[[modn]]
  }
  vpp <- length(mod$coef)
  tryCatch(
    mcmc.diagnostics(mod, vars.per.page = vpp),
    error = function(e) cat("MCMC was not run or MCMC sample was not stored."))
})
outputOptions(output, "diagnosticsplot", suspendWhenHidden = FALSE)

output$mcmcplotdownload <- downloadHandler(
  filename = function(){paste(nwname(),'_mcmc.pdf',sep='')},
  content = function(file){
    mod <- input$choosemodel_mcmc
    if(mod=="Current"){
      mod <- model1reac()
    } else {
      modn <- as.numeric(substr(mod,6,6))
      mod <- values$modelfits[[modn]]
    }
    vpp <- length(mod$coef)
    pdf(file=file, height=vpp*4/3, width=10)
    tryCatch(
      mcmc.diagnostics(model1reac(), vars.per.page = vpp),
      error = function(e) cat("MCMC was not run or MCMC sample was not stored."))
    dev.off()
  }
)

output$diagnosticsplotspace <- renderUI({
  if(input$fitButton == 0 | ergm.terms()=="NA"){
    return()
  }
  mod <- input$choosemodel_mcmc
  if(mod=="Current"){
    mod <- model1reac()
  } else {
    modn <- as.numeric(substr(mod,6,6))
    mod <- values$modelfits[[modn]]
  }
  vpp <- length(mod$coef)
  plotOutput('diagnosticsplot', height = vpp*400/2)
})

output$diagnostics <- renderPrint({
  if(input$fitButton == 0 | ergm.terms()=="NA"){
    return()
  }
  mod <- input$choosemodel_mcmc
  if(mod=="Current"){
    mod <- model1reac()
  } else {
    modn <- as.numeric(substr(mod,6,6))
    mod <- values$modelfits[[modn]]
  }
  isolate(tryCatch(
    mcmc.diagnostics(mod),
    error = function(e) cat("MCMC was not run or MCMC sample was not stored.")))
})
outputOptions(output, 'diagnostics', suspendWhenHidden=FALSE)


# Goodness of Fit ---------------------------------------------------------


# One drawback of the navbarPage layout option is that you can't specify
# certain elements or panels to show up on multiple pages. Furthermore,
# Shiny will not let you use the same piece of output from server.R twice
# in ui.R. Therefore, output$currentdataset2 and output$check2 are the same as
# output$currentdataset and output$check1 with different names.

# In the reactive section above the creation of model1gof depends on the term
# the user inputs. After checking that the user has already clicked the
# actionButton on the page we can output the text of the gof object and the
# plot of the gof object.

#dataset only updates after goButton on first tab has been clicked
output$currentdataset_gof <- renderPrint({
  if(!is.network(nw())){
    return(cat('Upload a network'))
  }
  cat(nwname())
})

output$uichoosemodel_gof <- renderUI({
  n <- values$modeltotal
  if(n == 0){
    inlineSelectInput("choosemodel_gof",label=NULL,
                      choices=c("Current"))
  } else {
    inlineSelectInput("choosemodel_gof",label=NULL,
                      choices=c(paste0("Model",1:n)))
  }
})
outputOptions(output,"uichoosemodel_gof",suspendWhenHidden=FALSE)
#formula only updates after fitButton has been clicked
output$checkterms_gof <- renderPrint({
  if(!is.network(nw())){
    return(cat('Upload a network'))
  }
  if(input$fitButton == 0){
    return(cat('Please fit a model'))
  }
  mod <- input$choosemodel_gof
  if(mod=="Current"){
    cat(isolate(ergm.terms()))
  } else {
    mod <- as.numeric(substr(mod,6,6))
    cat(values$modelformulas[[mod]])
  }
})

# state$gof will toggle between two states, depending on
# if gof plots are outdated compared to current ergm formula

observeEvent(input$fitButton, {
  state$gof <- 0 #gof plots are outdated
})

observeEvent(input$gofButton, {
  state$gof <- 1 #gof plots are up to date
})


output$gofsummary <- renderPrint({
  if (input$gofButton == 0 | state$gof == 0){
    return()
  }
  mod <- input$choosemodel_gof
  if(mod=="Current" | mod=="Model1"){
    gofobj <- isolate({model1gof()})
  } else {
    n <- substr(mod,6,6)
    gofobj <- isolate({switch(n, "2" = model2gof(),
                              "3" = model3gof(),
                              "4" = model4gof(),
                              "5" = model5gof())})
  }

  gofobj
  })
outputOptions(output, 'gofsummary', suspendWhenHidden=FALSE)

output$gofplot <- renderPlot({
  input$gofButton
  if(state$gof == 0){
    return()
  }
  gofterm <- isolate(input$gofterm)
  if (gofterm == 'Default'){
    if(is.directed(nw())){
      par(mfrow=c(4,1))
    } else {
      par(mfrow=c(3,1))
    }
    cex <- 1.5
  } else {
    par(mfrow=c(1,1))
    cex <- 1
  }
  mod <- input$choosemodel_gof
  if(mod=="Current" | mod=="Model1"){
    gofobj <- isolate({model1gof()})
  } else {
    n <- substr(mod,6,6)
    gofobj <- isolate({switch(n, "2" = model2gof(),
                              "3" = model3gof(),
                              "4" = model4gof(),
                              "5" = model5gof())})
  }
  par(cex.lab=cex, mar=c(5,4.2,2,2)+.1)
  isolate(plot.gofobject(gofobj, cex.axis=cex, main=NULL))
  par(mfrow=c(1,1), cex.lab=1, cex.axis=1)
})

output$gofplotdownload <- downloadHandler(
  filename = function(){paste(nwname(),'_gof.pdf',sep='')},
  content = function(file){

    pdf(file=file, height=4, width=10)
    mod <- input$choosemodel_gof
    if(mod=="Current" | mod=="Model1"){
      gofobj <- isolate({model1gof()})
    } else {
      n <- substr(mod,6,6)
      gofobj <- isolate({switch(n, "2" = model2gof(),
                                "3" = model3gof(),
                                "4" = model4gof(),
                                "5" = model5gof())})
    }
    isolate(plot.gofobject(gofobj, cex.axis=1, main=NULL))
    dev.off()
  }
)

output$gofplotspace <- renderUI({
  if(input$gofButton == 0){
    return()
  }
  if(state$gof == 0){
    return()
  }
  gofterm <- isolate(input$gofterm)
  if (gofterm == 'Default'){
    gofplotheight = 1000
  } else {
    gofplotheight = 400
  }
  plotOutput('gofplot', height=gofplotheight)
})

output$gofplotcomp <- renderPlot({
  if(state$gof == 0 | values$modeltotal == 0){
    return()
  }
  input$gofButton
  gofterm <- isolate(input$gofterm)
  n <- values$modeltotal
  if (gofterm == 'Default'){
    if(is.directed(nw())){
      cols <- isolate(4)
      bottomtext <- c("idegree","odegree","espartners","distance")
      bottommat <- c(0,(n*cols+1):(n*cols+4))
    } else {
      cols <- isolate(3)
      bottomtext <- c("degree","espartners","distance")
      bottommat <- c(0,(n*cols+1):(n*cols+3))
    }
    innermat <- matrix(1:(n*cols),ncol=cols, byrow=TRUE)

  } else {
    cols <- isolate(1)
    innermat <- matrix(1:n,ncol=1,nrow=n, byrow=TRUE)
    bottommat <- c(0,n*cols+1)
    bottomtext <- isolate(input$gofterm)
  }
  sidemat <- (bottommat[length(bottommat)]+1):(bottommat[length(bottommat)]+n)
  sidetext <- paste("Model",1:n)
  layoutmat <- cbind(sidemat,innermat)
  layoutmat <- rbind(layoutmat, bottommat)
  widths <- c(1,rep(3,cols))
  heights <- c(rep(3,n),1)
  layout(layoutmat, widths=widths, heights=heights)
  par(mar=c(1,2,2.2,1), omi=c(0,0,.3,0))

  isolate({
  for(j in 1:n){
    gofobj <- switch(j, "1" = model1gof(),
                              "2" = model2gof(),
                              "3" = model3gof(),
                              "4" = model4gof(),
                              "5" = model5gof())
    plot.gofobject(gofobj, main="", cex.axis=1)
  }
  par(mar=c(1,1,1,1))
  for(j in 1:cols){
    plot.new()
    text(x=.5,y=.7,labels=bottomtext[j], cex=1.5)
  }
  for(j in 1:n){
    plot.new()
    text(x=.7,y=.5,labels=sidetext[j], srt=90, cex=1.5)
  }
  mtext("Goodness-of-fit diagnostics",side=3,outer=TRUE,cex=1.5,padj=0)
  })
})

output$gofplotcompspace <- renderUI({
  if(input$gofButton==0 | state$gof==0){
    return()
  }
  gofterm <- isolate(input$gofterm)
  n <- isolate(values$modeltotal)
  plotheight = n*250
  if (gofterm == 'Default'){
    plotwidth = "100%"
  } else {
    plotwidth = "50%"
  }

  plotOutput('gofplotcomp', height=plotheight+100, width=plotwidth)
})

output$gofplotcompdownload <- downloadHandler(
  filename = function(){paste(nwname(),'_gofcomp.pdf',sep='')},
  content = function(file){
    gofterm <- input$gofterm
    n <- values$modeltotal
    if (gofterm == 'Default'){
      cols <- 3
      lastelt <- 15
      bottommat <- c(0,(n*cols+1):(n*cols+3))
      bottomtext <- c("degree","espartners","distance")
      if(n<=4){
        #landscape page
        pdf(file=file, height=8.5, width=11)
      } else {
        #portrait page
        pdf(file=file, height=11, width=8.5)
      }
    } else {
      cols <- 1
      lastelt <- 7
      bottommat <- c(0,n*cols+1)
      bottomtext <- input$gofterm
      pdf(file=file, height=11, width=8.5)
    }

    sidemat <- (bottommat[length(bottommat)]+1):(bottommat[length(bottommat)]+n)
    sidetext <- paste("Model",1:n)
    innermat <- matrix(1:(n*cols), ncol=cols, nrow=n, byrow=TRUE)
    layoutmat <- cbind(sidemat,innermat)
    if(n<3){
      extra <- 3-n
      last <- sidemat[length(sidemat)]
      filler <- matrix((last+1):lastelt, nrow=extra, byrow=TRUE)
      layoutmat <- rbind(layoutmat, bottommat, filler)
      heights <- c(rep(3,n),1,rep(3,3-n))
    } else {
      layoutmat <- rbind(layoutmat, bottommat)
      heights <- c(rep(3,n),1)
    }

    widths <- c(1,rep((9/cols),cols))
    layout(layoutmat, widths=widths, heights=heights)
    par(mar=c(1,2,2.2,1), omi=c(0,0,.3,0))

    for(j in 1:n){
      gofobj <- switch(j, "1" = model1gof(),
                       "2" = model2gof(),
                       "3" = model3gof(),
                       "4" = model4gof(),
                       "5" = model5gof())
      plot.gofobject(gofobj, cex.axis=1, main="")
    }
    par(mar=c(1,1,1,1))
    for(j in 1:cols){
      plot.new()
      text(x=.5,y=.7,labels=bottomtext[j], cex=1.5)
    }
    for(j in 1:n){
      plot.new()
      text(x=.7,y=.5,labels=sidetext[j], srt=90, cex=1.5)
    }
    mtext("Goodness-of-fit diagnostics",side=3,outer=TRUE,cex=1.5,padj=0)

    dev.off()
  }
)


# Simulations -------------------------------------------------------------

# On this page the user can choose how many simulations of the model to run.
# The reactive object model1simreac contains all the simulations, which we
# can output a summary of and choose one simulation at a time to plot. *Note:*
# when the user chooses to simulate one network, allmodelsimreac() is a
# reactive object of class network. When the user chooses to simulate multiple
# networks, allmodelsimreac() contains a list of the generated networks. This
# is why we have to split up the plot command in an if statement. The rest of
# the display options should look familiar from the 'Plot Network' tab.

output$uichoosemodel_sim <- renderUI({
  n <- values$modeltotal
  if(n == 0){
    inlineSelectInput("choosemodel_sim",label=NULL,
                      choices=c("Current"))
  } else {
    inlineSelectInput("choosemodel_sim",label=NULL,
                      choices=c(paste0("Model",1:n)))
  }
})
outputOptions(output,"uichoosemodel_sim",suspendWhenHidden=FALSE)

output$checkterms_sim <- renderPrint({
  if(!is.network(nw())){
    return(cat('Upload a network'))
  }
  if(input$fitButton == 0){
    return(cat('Please fit a model'))
  }
  cat(current.simterms())

})
output$currentdataset_sim <- renderPrint({
  if(!is.network(nw())){
    return(cat('Upload a network'))
  }
  cat(nwname())
})

# state$sim will toggle between two states, depending on
# if simulations are outdated compared to current ergm network/formula
state$sim <- 0

observe({
  nwname()
  input$choosemodel_sim
  input$fitButton
  state$sim <- 0 #simulations are outdated
  updateNumericInput(session, "thissim",
                     label = "Choose a simulation to plot:",
                     value = 1, min = 1, max = input$nsims)
})

observe({
  input$simButton
  state$sim <- 1 #simulations plots are up to date
})

observe({
  if(input$simcontroldefault){
    updateNumericInput(session, "simMCMCinterval", value=1024)
    #burn-in gets updated separately, to be 16*interval
    disableWidget("simMCMCinterval", session, disabled=TRUE)
    disableWidget("simMCMCburnin", session, disabled=TRUE)
    disableWidget("simcustomMCMCcontrol", session, disable=TRUE)
  } else {
    disableWidget("simMCMCinterval", session, disabled=FALSE)
    disableWidget("simMCMCburnin", session, disabled=FALSE)
    disableWidget("simcustomMCMCcontrol", session, disable=FALSE)
  }
})

observe({
  input$simMCMCinterval
  updateNumericInput(session, "simMCMCburnin", value=16*input$simMCMCinterval)
})

output$simnum <- renderText({
  input$simButton
  n <- isolate(input$nsims)
  n
})
outputOptions(output, 'simnum', suspendWhenHidden=FALSE)

output$simsummary <- renderPrint({
  if (input$simButton == 0 | state$sim == 0){
    return(cat(''))
  }
  sim <- isolate(allmodelsimreac())
  n <- isolate(input$nsims)
  terms <- isolate(current.simterms())
  sum <- list()
  sum[1] <- paste(" Number of Networks: ",n, "\n")
  sum[2] <- paste("Model: ", nwname(),' ~ ',terms, "\n")
  sum[3] <- paste("Reference: ", format(attr(sim,'reference')), "\n")
  sum[4] <- paste("Constraints: ", format(attr(sim, 'constraints')), "\n")
  sum[5] <- paste("Parameters: \n")
  cat(unlist(sum))
})

output$simsummary2 <- renderPrint({
  if (input$simButton == 0 | state$sim == 0){
    return(cat(''))
  }
  sim <- isolate(allmodelsimreac())
  sim
})

output$simcoef <- renderPrint({
  if (input$simButton == 0 | state$sim == 0){
    return(cat(''))
  }
  sim <- isolate(allmodelsimreac())
  c <- attr(sim, 'coef')
  c <- cbind(format(names(c)),format(c, digits=3))
  write.table(c, quote=FALSE, row.names=FALSE, col.names=FALSE)
})

output$simstatslabel <- renderPrint({
  if (input$simButton == 0 | state$sim == 0){
    return(cat(''))
  }
  cat(' Stored network statistics:')
})

output$simstats <- renderPrint({
  if (input$simButton == 0 | state$sim == 0){
    return(cat(''))
  }
  sim <- isolate(allmodelsimreac())
  m <- format(t(attr(sim,'stats')))
  m <- cbind(format(rownames(m)),m)
  write.table(m, quote=F, row.names=F, col.names=F)
})

output$simstats2 <- renderPrint({
  if(input$simButton == 0 | state$sim == 0){
    return()
  }
  sim <- isolate(allmodelsimreac())
  m <- attr(sim,'stats')
  mat <- cbind(apply(m,2,min),apply(m,2,max),apply(m,2,IQR),
               apply(m,2,mean),apply(m,2,sd))
  colnames(mat) <- c("min","max","IQR","mean","SD")
  mat
})

output$simstatsdownload <- downloadHandler(
  filename = function(){paste(nwname(),'_simstats',
                              input$simstatsfiletype,sep='')},
  contentType = "text/csv",
  content = function(file) {
    x<-attr(allmodelsimreac(),"stats")
    if(input$simstatsfiletype == ".csv"){
      write.csv(x, file)
    } else {
      if(input$nsims > 1){
        capture.output({
            sim <- allmodelsimreac()
            terms <- current.simterms()
            n <- input$nsims
            sum <- list()
            sum[1] <- paste(" Number of Networks: ",n, "\n")
            sum[2] <- paste("Model: ", nwname(),' ~ ',terms, "\n")
            sum[3] <- paste("Reference: ",
                            format(attr(sim,'reference')), "\n")
            sum[4] <- paste("Constraints: ",
                            format(attr(sim, 'constraints')), "\n")
            sum[5] <- paste("Parameters: \n")
            cat(unlist(sum))

            c <- attr(sim, 'coef')
            c <- cbind(format(names(c)),format(c, digits=3))
            write.table(c, quote=FALSE, row.names=FALSE, col.names=FALSE)
            cat("\n\n")
            cat(' Stored network statistics:\n')

            mat <- cbind(apply(x,2,min),apply(x,2,max),apply(x,2,IQR),
                         apply(x,2,mean),apply(x,2,sd))
            colnames(mat) <- c("min","max","IQR","mean","SD")
            cat("\n")
            print(mat)
            rownames(x) <- rep("",n)
            print(t(x))
        },file=file)
      } else {
        capture.output(allmodelsimreac(),file=file)
      }
    }

  }
  )

output$simstatsplot <- renderPlot({
  if(input$simButton==0 | state$sim == 0){
    return()
  }
  sim <- isolate(allmodelsimreac())
  simstats <- attr(sim,'stats')
  targetstats <- summary(current.simformula())
  matplot(1:nrow(simstats), simstats, pch=c(1:8),
          col=c('red', 'blue', 'green3', 'cyan', 'magenta3',
                'orange', 'black', 'grey', 'yellow'),
          xlab="Simulations",ylab="")
  abline(h=c(targetstats),
         col=c('red', 'blue', 'green3', 'cyan', 'magenta3',
               'orange', 'black', 'grey', 'yellow'))
})

output$simstatslegend <- renderPlot({
  if(input$simButton==0 | state$sim == 0){
    return()
  }
  sim <- isolate(allmodelsimreac())
  termnames <- colnames(attr(sim,'stats'))
  color <- c('red', 'blue', 'green3', 'cyan', 'magenta3',
             'orange', 'black', 'grey', 'yellow')
  plot.new()
  par(xpd=TRUE)
  legend('topleft', inset=c(-.15,0),
         legend=c(termnames, paste("target stat:",termnames)),
         pch=c(1:length(termnames), rep(NA,length(termnames))),
         lty=c(rep(NA,length(termnames)),rep(1,length(termnames))),
         col=rep(color[1:length(termnames)],2))
})

output$simstatsplotdownload <- downloadHandler(
  filename = function(){paste(nwname(),'_simstatsplot.pdf',sep='')},
  content = function(file){
    pdf(file=file, height=8, width=15)
    sim <- isolate(allmodelsimreac())
    simstats <- attr(sim,'stats')
    termnames <- colnames(attr(sim,'stats'))
    targetstats <- summary(current.simformula())
    color <- c('red', 'blue', 'green3', 'cyan', 'magenta3',
               'orange', 'black', 'grey', 'yellow')
    par(mar=c(5,4,4,12)) #increase margin on right to accommodate legend
    matplot(1:nrow(simstats), simstats, pch=c(1:8),
            col=color,
            xlab="Simulations",ylab="")
    abline(h=c(targetstats),
           col=color)
    par(xpd=TRUE)
    legend('topright', inset=c(-.2,0),
           legend=c(termnames, paste("target stat:",termnames)),
           pch=c(1:length(termnames), rep(NA,length(termnames))),
           lty=c(rep(NA,length(termnames)),rep(1,length(termnames))),
           col=rep(color[1:length(termnames)],2))
    dev.off()
  }
)

output$dynamiccolor2 <- renderUI({
  selectInput('colorby2',
              label = 'Color vertices according to:',
              c('None' = 2, attrib()),
              selected = 2)
})
outputOptions(output,'dynamiccolor2',suspendWhenHidden=FALSE, priority=10)

output$dynamicsize2 <- renderUI({
  selectInput('sizeby2',
              label = 'Size vertices according to:',
              c('None' = 1, 'Betweenness',numattr()))
})


output$simplot <- renderPlot({
  if(input$simButton == 0 | state$sim == 0){
    return()
  }
  nw_var <- nw()
  nsims <- isolate(input$nsims)
  model1sim <- isolate(allmodelsimreac())

  #can't plot simulation number greater than total sims
  if(input$thissim > nsims){
    return()
  }

  color <- adjustcolor(vcol2(), alpha.f = input$transp2)
  par(mar = c(0, 0, 0, 0))
  if (nsims == 1){
    plot(model1sim, coord = sim.coords.1(),
         displayisolates = input$iso2,
         displaylabels = input$vnames2,
         vertex.col = color,
         vertex.cex = nodesize2())
  } else {
    plot(model1sim[[input$thissim]],
         coord = sim.coords.2(),
         displayisolates = input$iso2,
         displaylabels = input$vnames2,
         vertex.col = color,
         vertex.cex = nodesize2())
  }
  if(input$colorby2 != 2){
      legend('bottomright', title=input$colorby2,
             legend = legendlabels2(), fill = legendfill2(),
             bty='n')
    }

})

output$simplotdownload <- downloadHandler(
  filename = function(){paste(nwname(),'_simplot.pdf',sep='')},
  content = function(file){
    pdf(file=file, height=10, width=10)
    color <- adjustcolor(vcol2(), alpha.f = input$transp2)
    if(input$nsims == 1){
    plot(allmodelsimreac(),
         coord = sim.coords.1(),
         displayisolates = input$iso2,
         displaylabels = input$vnames2,
         vertex.col = color,
         vertex.cex = nodesize2())
    }else{
      plot(allmodelsimreac()[[input$thissim]],
           coord = sim.coords.2(),
           displayisolates = input$iso2,
           displaylabels = input$vnames2,
           vertex.col = color,
           vertex.cex = nodesize2())
    }
    if(input$colorby2 != 2){
      legend('bottomright', title=input$colorby2,
             legend = legendlabels2(), fill = legendfill2())
    }
    dev.off()
  }
)


})
