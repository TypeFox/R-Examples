# Load packages
if(!require("shiny")){
  stop("The package 'shiny' is required to use this application.")
}
if(!require("diveRsity")){
  stop("The package 'diveRsity' is required to use this application.")
}
if(!require(qgraph)){
  stop("The package 'qgraph' is required to use this application.")
}
if(!require("shinyIncubator")){
  stop("The package 'shinyIncubator' is required to use this application.")
}
# Increase file upload limit
options(shiny.maxRequestSize = 50*(1024^2))

shinyServer(function(input, output, session){
  
  # calculate migration (only re calculate id number of bootstraps is changed)
  out <- reactive({
    if(is.null(input$file)){
      return(NULL)
    }
    withProgress(session, min=1, max=15, {
      setProgress(message = 'Calculation in progress',
                  detail = 'This may take a while...')
      for (i in 1:15) {
        setProgress(value = i)
        Sys.sleep(0.1)
      }
      #isolate({
      infile <- input$file$datapath
      diveRsity:::divMigrateOnline(infile = infile,
                                   nbs = input$nbs,
                                   stat = "all",
                                   para = FALSE)
    })
  })
  
  # population exclusion dialogue
  output$popnames <- renderUI({
    if(is.null(input$file)){
      return(NULL)
    }
    op <- out()
    if(input$nbs != 0L){
      if(input$stat == "D"){
        dat <- list(op$D, op$D_sig)
      } else if(input$stat == "Gst"){
        dat <- list(op$Gst, op$G_sig)
      } else if(input$stat == "Nm"){
        dat <- list(op$Nm, op$Nm_sig)
      }
    } else {
      if(input$stat == "D"){
        dat <- list(op$D)
      } else if(input$stat == "Gst"){
        dat <- list(op$Gst)
      } else if(input$stat == "Nm"){
        dat <- list(op$Nm)
      } 
    }
    pops <- colnames(dat[[1]])
    #pops <- as.character(sapply(pops, function(x) gsub(",", "", x)))
    checkboxGroupInput("pops", HTML("<br><h5>",
                                    "5. Exclude populations from network?",
                                    "</h5>"),
                       choices = pops,
                       selected = NULL,
                       inline = TRUE)
  })
  
  # plot download button
  output$pltDL <- renderUI({
    if(is.null(input$file)){
      return(NULL)
    }
    downloadButton("dlPlt", "Download Network")
  })
  
  # re-standardise plots following sample exclusion
  output$stdplt <- renderUI({
    if(is.null(input$file)){
      return(NULL)
    }
    radioButtons("restand", h5("6. Re-standardize network connections?"),
                 choices = c("Y", "N"),
                 selected = "N",
                 inline = TRUE) 
  })
  
  # plot download format dialogue
  output$pltFormat <- renderUI({
    if(is.null(input$file)){
      return(NULL)
    }
    radioButtons("format", h5("7. Network download format."),
                 c("pdf", "png", "eps"), inline = TRUE, 
                 selected = "pdf")
  })
  
  # network plots
  output$plt <- renderPlot({
    if(is.null(input$file)){
      return(NULL)
    }
    op <- out()
    #if(input$nbs == op$nbs){
    if(input$nbs != 0L){
      if(input$stat == "D"){
        dat <- list(op[[1]], op[[2]])
      } else if(input$stat == "Gst"){
        dat <- list(op[[3]], op[[4]])
      } else if(input$stat == "Nm"){
        dat <- list(op[[5]], op[[6]])
      }
      diag(dat[[2]]) <- FALSE
      dat[[1]][!dat[[2]]] <- 0
      dat <- dat[[1]]
      if(!is.null(input$pops)){
        idx <- sapply(input$pops, function(x){
          which(colnames(dat) == x)
        })
        dat <- dat[-idx, -idx]
        # re-standardize dat
        if(input$restand == "Y"){
          dat <- dat/max(dat, na.rm = TRUE) 
        }
      }
    } else {
      if(input$stat == "D"){
        dat <- list(op[[1]])
      } else if(input$stat == "Gst"){
        dat <- list(op[[2]])
      } else if(input$stat == "Nm"){
        dat <- list(op[[3]])
      } 
      dat <- dat[[1]]
      if(!is.null(input$pops)){
        idx <- sapply(input$pops, function(x){
          which(colnames(dat) == x)
        })
        dat <- dat[-idx, -idx]
        # re-standardize dat
        if(input$restand == "Y"){
          dat <- dat/max(dat, na.rm = TRUE) 
        }
      }
    }
    dat[is.na(dat)] <- 0
    if(input$filter_threshold != 0L){
      dat[dat <= input$filter_threshold] <- 0
    }
    #if(input$goButton == 0L) return(NULL)
    #isolate({
    qgraph::qgraph(dat, nodeNames = colnames(dat), posCol = "darkblue",
                   legend = TRUE, edge.labels = TRUE, 
                   mar = c(2,2,5,5), curve = 2.5)
    if(input$nbs != 0L){
      title(paste("Relative migration network (Filter threshold = ", 
                  input$filter_threshold, "; ", input$nbs, 
                  " bootstraps; ", input$stat, " method)", sep = ""))
    } else {
      title(paste("Relative migration network (Filter threshold = ", 
                  input$filter_threshold, "; ", input$stat, ")", 
                  sep = "")) 
    }
    #})
    
    #}
  })
  
  # Matrix table
  output$mat <- downloadHandler(
    filenames <- function(file){
      if(is.null(input$file)){
        return(NULL)
      }
      if(input$stat == "D"){
        paste("divMigrate-online_", gsub(" ", "_", date()), "-[D].txt", 
              sep = "") 
      } else if(input$stat == "Gst"){
        paste("divMigrateonline_", gsub(" ", "_", date()), "-[Gst].txt", 
              sep = "")
      } else if(input$stat == "Nm"){
        paste("divMigrateonline_", gsub(" ", "_", date()), "-[Nm].txt", 
              sep = "")
      }
    },
    content <- function(file){
      if(is.null(input$file)){
        return(NULL)
      }
      op <- out()
      if(input$nbs != 0L){
        if(input$stat == "D"){
          dat <- list(op[[1]], op[[2]])
        } else if(input$stat == "Gst"){
          dat <- list(op[[3]], op[[4]])
        } else if(input$stat == "Nm"){
          dat <- list(op[[5]], op[[6]])
        }
        diag(dat[[2]]) <- FALSE
        dat[[1]][dat[[2]]] <- paste(dat[[1]][dat[[2]]], "*", sep = "")
        dat <- dat[[1]]
      } else {
        if(input$stat == "D"){
          dat <- list(op$D)
        } else if(input$stat == "Gst"){
          dat <- list(op$Gst)
        } else if(input$stat == "Nm"){
          dat <- list(op$Nm)
        } 
        dat <- dat[[1]]
      }
      if(input$stat == "D"){
        hdr <- paste(c(" Pops", colnames(dat)), collapse = "\t")
        hdr <- c("divMigrate-online",
                 "",
                 "Relative directional Migration matrix Calculated using",
                 "D (Jost, 2008)",
                 "",
                 "If bootstrapping has been carried out, statistically",
                 "significant values (alpha = 0.5) are indicated with '*'",
                 "",
                 hdr)
        dat <- cbind(colnames(dat), dat)
        dat <- apply(dat, 1, paste, collapse = " \t")
        dat <- paste(c(hdr, dat), collapse = "\n")
        writeLines(text = dat, con = file)
      } else if(input$stat == "Gst"){
        hdr <- paste(c(" Pops", colnames(dat)), collapse = "\t")
        hdr <- c("divMigrate-online",
                 "",
                 "Relative directional Migration matrix Calculated using",
                 "Gst (Nei & Chesser, 1983)",
                 "",
                 "If bootstrapping has been carried out, statistically",
                 "significant values (alpha = 0.5) are indicated with '*'",
                 "",
                 hdr)
        dat <- cbind(colnames(dat), dat)
        dat <- apply(dat, 1, paste, collapse = " \t")
        dat <- paste(c(hdr, dat), collapse = "\n")
        writeLines(text = dat, con = file)
      } else if(input$stat == "Nm"){
        hdr <- paste(c(" Pops", colnames(dat)), collapse = "\t")
        hdr <- c("divMigrate-online",
                 "",
                 "Relative directional Migration matrix Calculated using",
                 "Nm (eqn. 12 from Alcala et al, 2014)",
                 "",
                 "If bootstrapping has been carried out, statistically",
                 "significant values (alpha = 0.5) are indicated with '*'",
                 "",
                 hdr)
        dat <- cbind(colnames(dat), dat)
        dat <- apply(dat, 1, paste, collapse = " \t")
        dat <- paste(c(hdr, dat), collapse = "\n")
        writeLines(text = dat, con = file)
      }
    }
  )
  
  # Downloadable plots
  output$dlPlt <- downloadHandler(
    #if(input$goButton==0) return(NULL)
    filename = function() {
      if(is.null(input$file)){
        return(NULL)
      }
      if(input$stat == "D"){
        paste("divMigrate-online_", gsub(" ", "_", date()), "-[D-Network].",
              input$format, sep = "") 
      } else if(input$stat == "Gst"){
        paste("divMigrate_", gsub(" ", "_", date()), "-[Gst-Network].",
              input$format, sep = "")
      } else if(input$stat == "Nm"){
        paste("divMigrate_", gsub(" ", "_", date()), "-[Nm-Network].",
              input$format, sep = "")
      }
    },
    content = function(file) {
      if(is.null(input$file)){
        return(NULL)
      }
      op <- out()
      tmp <- tempfile()
      on.exit(unlink(tmp))
      if(input$nbs != 0L){
        if(input$stat == "D"){
          dat <- list(op[[1]], op[[2]])
        } else if(input$stat == "Gst"){
          dat <- list(op[[3]], op[[4]])
        } else if(input$stat == "Nm"){
          dat <- list(op[[5]], op[[6]])
        }
        diag(dat[[2]]) <- FALSE
        dat[[1]][!dat[[2]]] <- 0
        dat <- dat[[1]]
        if(!is.null(input$pops)){
          idx <- sapply(input$pops, function(x){
            which(colnames(dat) == x)
          })
          dat <- dat[-idx, -idx]
        }
      } else{
        if(input$stat == "D"){
          dat <- list(op$D)
        } else if(input$stat == "Gst"){
          dat <- list(op$Gst)
        } else if(input$stat == "Nm"){
          dat <- list(op$Nm)
        } 
        dat <- dat[[1]]
        if(!is.null(input$pops)){
          idx <- sapply(input$pops, function(x){
            which(colnames(dat) == x)
          })
          dat <- dat[-idx, -idx]
        }
      }
      dat[is.na(dat)] <- 0
      if(input$filter_threshold != 0L){
        dat[dat <= input$filter_threshold] <- 0
      }
      if(input$format == "eps"){
        setEPS()
        postscript(file = tmp)
      } else if(input$format == "pdf"){
        pdf(file = tmp, paper = "a4r")
      } else if(input$format == "png"){
        png(file = tmp, width = 800, height = 750)
      }
      qgraph::qgraph(dat, nodeNames = colnames(dat), posCol = "darkblue",
                     legend = TRUE, edge.labels = TRUE, 
                     mar = c(2,2,5,5), curve = 2.5)
      if(input$nbs != 0L){
        title(paste("Relative migration network (Filter threshold = ", 
                    input$filter_threshold, "; ", input$nbs, 
                    " bootstraps; ", input$stat, " method)", sep = ""))
      } else {
        title(paste("Relative migration network (Filter threshold = ", 
                    input$filter_threshold, "; ", input$stat, ")", sep = "")) 
      }
      dev.off()
      bytes <- readBin(tmp, "raw", file.info(tmp)$size)
      writeBin(bytes, file)
    }
  )
  
  
})