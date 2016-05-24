RRprofReport <- function (file.name = "RRprof.out", notepad.path = "C:/Program Files/Notepad++/notepad++.exe", reportname = "my_report") {
  profdata <- readLines(file.name)                                        # Read the profile file
  interval <- as.numeric(strsplit(profdata[1L], "=")[[1L]][2L])/1e+06     # Get the interval time
  filelines <- grep("#File", profdata)                                    # Find the files that correspond to each function
  files <- profdata[filelines]
  files <- strsplit(files, ": ")
  files <- sapply(files, "[", 2)
  files <- gsub("\\\\", "/", files) # Changing the convention to slash for windows
  filesnames <- unlist(lapply(files, function(files) {
    x <- strsplit(files[[1]], "/|\\\\")
    x <- x[[1]][length(x[[1]])]
    x <- strsplit(x, ".R")[[1]][1]
  }))


  filesread <- list()
  # Some of the files correspond to functions from packages in tmp directories that dissapear after the excution.
  # Get rid of them.
  for (i in 1:length(files)) filesread[[i]] <- suppressWarnings(try(readLines(files[i],warn = FALSE), silent = TRUE))
  # filesread is a list with the content of each of the files

  fileserror <- which(unlist(lapply(filesread, function(x) !is(x,"try-error"))))

  if (length(filelines) > 1)
    profdata <- profdata[-filelines[2:length(filelines)]]
  if (length(filelines) > 0)
    profdata <- profdata[-1:-filelines[1]]
  profdata <- profdata[grep("#", profdata)] # Find the lines in the files
  total.time2 <- interval * length(profdata)
  profdata <- gsub("\\\"| $", "", profdata) # Remove inverted commas
  total.time <- list()
  self.time <- list()

  # Get the time expent on each line of code
  for (i in fileserror) {
    total.time[[i]] <- rep(0, length(filesread[[i]]))
    self.time[[i]] <- rep(0, length(filesread[[i]]))
    for (j in 1:length(filesread[[i]])) {
      tmp <- grep(paste(i, "#", j, " ", sep = ""), profdata)
      total.time[[i]][j] <- length(tmp)
      if (length(tmp) > 0) {
        for (k in 1:length(tmp)) {
          tmp2 <- strsplit(strsplit(profdata[tmp[k]],
                                    "#")[[1]][1], " ")[[1]][length(strsplit(strsplit(profdata[tmp[k]],
                                                                                     "#")[[1]][1], " ")[[1]])] == i & strsplit(strsplit(profdata[tmp[k]],
                                                                                                                                        "#")[[1]][2], " ")[[1]][1] == j
          if (tmp2) {
            self.time[[i]][j] <- self.time[[i]][j] +
              1
          }
        }
      }
    }
  }
  total.time <- lapply(total.time, function(x, y) return(x * y), interval)
  self.time <- lapply(self.time, function(x, y) return(x * y), interval)
  total.time.sum <- unlist(lapply(total.time, function(i) sum(i)))
  self.time.sum <- unlist(lapply(self.time, function(i) sum(i)))
  total.time.ptg <- as.numeric(total.time.sum/total.time2 * 100)
  self.time.ptg <- as.numeric(self.time.sum/total.time2 * 100)
  maxtime <- suppressWarnings(unlist(lapply(total.time, function(x) max(x))))
  color <- list()

  # Get the color (reddish for the lines where time is spent)
  for (i in fileserror) color[[i]] <- unlist(lapply(total.time[[i]], function(x, y) return(floor(x/y * 255)), max(maxtime)))

  # Nozzle.R1 stuff to generate the report
  r <- newCustomReport("Profile Report")

  # Scripts for the notepad++
  script1 <- "<script type=\"text/javascript\" language=\"javascript\">\n  \n  function RunFile"
  script2 <- paste("(x) {\n                 \n                 WshShell = new ActiveXObject(\"WScript.Shell\");\n                 \n                 WshShell.Run('\"",
                   notepad.path, "\"", sep = "")
  script3 <- " -n'+x);\n  \n}\n  \n  </script>"

  list2 <- fileserror
  wd <- getwd()

  setwd(Sys.getenv("HOME"))


  list1 <- grep("~", files)
  if (length(list1 > 0)) {
    list2 <- list2[-list1]
  }
  script <- list()

  # Stuff4Notepad++
  for (i in list1) {
    script[[i]] <- newHtml(paste(script1, i, script2, " \"",
                                 getwd(), strsplit(files[i], "~")[[1]][2], "\"", script3,
                                 sep = ""))
    r <- addTo(r, script[[i]])
  }
  for (i in list2) {
    script[[i]] <- newHtml(paste(script1, i, script2, " \"",
                                 files[i], "\"", script3, sep = ""))
    r <- addTo(r, script[[i]])
  }
  setwd(wd)
  # End of stuff4Notepad++

  # SummarySection
  summary.table <- data.frame(filesnames[fileserror], total.time.sum[fileserror],
                              self.time.sum[fileserror], total.time.ptg[fileserror],
                              self.time.ptg[fileserror])
  colnames(summary.table) <- c("File", "Total time", "Self time",
                               "Total time (%)", "Self time (%)")
  s1 <- newSection("Summary")
  s2 <- list()
  for (i in fileserror) s2[[i]] <- newSection(filesnames[i])
  s3 <- newSection("Call graph")
  sum <- newTable(summary.table)
  code <- list()
  for (i in fileserror) code[[i]] <- rep(" ", length(filesread[[i]]))
  for (i in fileserror) {
    for (j in 1:length(filesread[[i]])) {
      code[[i]][j] <- asCode(filesread[[i]][j])
    }
  }
  files.tables <- list()
  f.tables <- list()
  table0 <- "<div class=\"table\" id="
  table1 <- ">\n  \n  <table class=\"resulttable tablesorter sortabletable\">"
  table2 <- "                 <thead>\n  <tr>\n  <th class=\"header\">\n  time\n  </th>\n  <th class=\"header headerSortDown\">\n  line\n  </th>\n  <th class=\"header\">\n  </th>\n  </tr>\n  </thead>\n  <tbody>"
  table3 <- "</tbody>\n  </table>\n  </div>"
  # End of summary section

  # Section for each function
  for (i in fileserror) {
    files.tables[[i]] <- paste(table0, "\"File", i, "\"",
                               table1, table2, sep = "")
  }
  for (i in fileserror) {
    for (j in 1:length(filesread[[i]])) {
      #browser()
      if (total.time[[i]][j] > 0)
        # Complex part to place the colors
        files.tables[[i]] <- paste(files.tables[[i]],
                                   "<tr><td>", total.time[[i]][j], "</td><td><a onclick=\"RunFile",
                                   i, "(", j, ")\" href=\"#\">", j, "</a></td><td>",
                                   "<span style=\"background-color: #", sprintf("FF%02X%02X",
                                                                                255 - color[[i]][j], 255 - color[[i]][j]),
                                   "\">", asCode(filesread[[i]][j]), "</span>", "</td></tr>",
                                   sep = "")
      else files.tables[[i]] <- paste(files.tables[[i]],
                                      "<tr><td></td><td><a>", j, "</a></td><td>", asCode(filesread[[i]][j]),
                                      "</td></tr>", sep = "")
    }
    files.tables[[i]] <- paste(files.tables[[i]], table3)
    f.tables[[i]] <- newHtml(files.tables[[i]])
  }



  #Section call graph
  setwd(wd)
  pd <- readProfileData(file.name)
  setwd(tempdir())
  png(filename = "callgraph.png", bg = "transparent", width = 960, height = 960)
  plotProfileCallGraph(pd, style = google.style, score = "total", nodeSizeScore = "none")
  dev.off()
  fig<-newFigure(file = "callgraph.png")


  # build report

  s3 <- addTo(s3,fig)
  for (i in fileserror) {
    s2[[i]] <- addTo(s2[[i]], f.tables[[i]])
  }
  s1 <- addTo(s1, sum)
  r <- addTo(r, s1)
  r <- addTo(r,s3)
  for (i in fileserror) r <- addTo(r, s2[[i]])
  setwd(tempdir())
  writeReport(r, filename = reportname)
  # Done!!!




  # New Code to open the report in the viewer of Rstudio

  setwd(wd)
  # If Rstudio is running, use its viewer
  if (isAvailable()) {
    viewer(paste(tempdir(), "/", reportname, ".html",sep = ""))
    # Prepare the source markers for Rstudio
    markers <- list()
    for (file in 1:length(files)) {
      total.time.val <- total.time[[file]]
      dummy <- sum(total.time.val >0)
      # Retain in the markers the lines in the top 25% sorted by time (first markers are most costly lines).
      dummy <- ceiling(dummy/4)
      lines <- order(-total.time.val)[1:dummy]
      total.time.val <- total.time.val[lines]
      for (n in 1:dummy) {
        marker <- list()
        if (n==1)
          marker$type <- "warning"   # Most costly line
        else
          marker$type <- "info"      # Rest of them
        marker$file <- files[file]
        marker$line <- lines[n]
        marker$column <- 1
        message <- paste(total.time.val[n], " secs.")
        if (n==1)
          message <- paste(message, "Most consuming line in the function")
        marker$message <- message;
        markers[[length(markers)+1]] <- marker
      }
    }

    rstudioapi::callFun("sourceMarkers",
                        name = "GUIProfiler",
                        markers = markers,
                        autoSelect = "first")
  } else if (Sys.info()["sysname"] == "Windows") {
    aux <- suppressWarnings(try(system2("C:/Program Files/Internet Explorer/iexplore.exe",
                                        args = paste(tempdir(), "\\", reportname, ".html",
                                                     sep = ""), wait = FALSE)))
    if (aux != 0)
      browseURL(paste(tempdir(), "/", reportname, ".html",
                      sep = ""))
  } else browseURL(paste(tempdir(), "/", reportname, ".html",
                         sep = ""))
}
