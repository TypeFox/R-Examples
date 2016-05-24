
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



########################################################################
## Hier sind Funktionen, die vermutlich kaum ein Anwender interessiert
## und wohl eher zum internen Gebrauch, wenn auch im NAMESPACE stehend;
##
## Ausnahme: rfGenerate* stehen in generatemodels.R
########################################################################

PrintModelList <-function (operators=FALSE, internal=FALSE,
                           newstyle=TRUE) {
   stopifnot(internal <=2, internal >=0, operators<=1, operators>=0)
    .C("PrintModelList", as.integer(internal), as.integer(operators),
       as.integer(newstyle),
       PACKAGE="RandomFields")
    invisible()
}



checkExamples <- function(exclude=NULL, include=1:length(.fct.list),
                          ask=FALSE, echo=TRUE, halt=FALSE, ignore.all=FALSE,
                          path="RandomFields", package="RandomFields",
                          read.rd.files=TRUE,
                          libpath = NULL, single.runs = FALSE) {
  .exclude <- exclude
  .ask <- ask
  .echo <- echo
  .halt <- halt
  .ignore.all <- ignore.all
  .path <- path
  .package <- package
  useDynLib <- importClassesFrom <- import <-
    importFrom <- exportClasses <-
      importMethodsFrom <- exportMethods <- S3method <- function(...) NULL
  .env <- new.env()

  export <- function(...) {
    ## code from 'rm'
    dots <- match.call(expand.dots = FALSE)$...
    z <-deparse(substitute(...))
    if (length(dots) && !all(sapply(dots, function(x) is.symbol(x) || 
                                    is.character(x)))) 
      stop("... must contain names or character strings")
    z <- sapply(dots, as.character)
    assign("export", c(get("export", .env), z), .env)
  }
  assign("export", NULL, .env)
  assign("all.pattern", FALSE, .env)
  assign("p", NULL, .env)
  source(paste(.path, "NAMESPACE", sep="/"), local=TRUE)
  if (is.logical(read.rd.files) && !read.rd.files) {
    .package.env <- parent.env(.GlobalEnv)
    while (attr(.package.env, "name") != paste("package:", .package, sep="")) {
      .package.env <- parent.env(.package.env)
    }
    .orig.fct.list <- ls(envir=.package.env)
    .ok <- (get("all.pattern", .env) |
            substr(.orig.fct.list, 1, 1) %in% get("p", .env) | 
            .orig.fct.list %in% get("export", .env))
    .fct.list <- .orig.fct.list[.ok]
  } else {
    .path <- read.rd.files
    if (is.logical(.path))
      .path <- "/home/schlather/svn/RandomFields/RandomFields/man"
    .files <- dir(.path, pattern="d$")
    .fct.list <- character(length(.files))
    for (i in 1:length(.files)) {
      cat(i, .path, .files[i], "\n")
      #if (i == 152) {cat("jumped\n"); next}
      .content <- scan(paste(.path, .files[i], sep="/") , what=character(),
                       quiet=TRUE)
      .content <- strsplit(.content, "alias\\{")
      .content <- .content[which(lapply(.content, length) > 1)][[1]][2]
      .fct.list[i] <-
        strsplit(strsplit(.content,"\\}")[[1]][1], ",")[[1]][1]
    }
  }
  Print(.fct.list, length(.fct.list)) #
  .include <- include
  RFoptions(graphics.close_screen = TRUE, graphics.split_screen = TRUE)
  .RFopt <- RFoptions()
  .not_working_no <- .not_working <- NULL
  .included.fl <- .fct.list[.include]
  .not.isna <- !is.na(.included.fl)
  .include <- .include[.not.isna]
  .included.fl <- .included.fl[.not.isna]
  Print(.included.fl, get("export", .env), get("p", .env)); #return(NA)
  .max.fct.list <- max(.included.fl)
  if (single.runs) {
    file.in <- "example..R"
    file.out <- "example..Rout"
    if (file.exists(file.out)) file.remove(file.out)
  }
  
  for (.idx in .include) {
    try(repeat dev.off(), silent=TRUE)
    if (.idx %in% .exclude) next
    cat("\n\n\n\n\n", .idx, " ", .package, ":", .fct.list[.idx],
        " (total=", length(.fct.list), ") \n", sep="")
    RFoptions(LIST=.RFopt)
    .C("ResetWarnings", as.integer(FALSE))
    if (.echo) cat(.idx, "")
    .tryok <- TRUE
    if (single.runs) {
      txt <- paste("library(RandomFields", libpath, "); example(",
                .fct.list[.idx],
                ", ask =", .ask,
                ", echo =", .echo,
                ")", sep="")
      write(file=file.in,  txt)
      command <- paste("R < ", file.in, ">>", file.out)
    } else {
      if (.halt) {
         do.call(utils::example, list(.fct.list[.idx], ask=.ask, echo=.echo))
      } else {
        if (is(try(do.call(utils::example, list(.fct.list[.idx], ask=.ask,
                                                echo=.echo))), "try-error")) {
          .not_working_no<- c(.not_working_no, .idx)
          .not_working <- c(.not_working, .fct.list[.idx])
          .tryok <- FALSE
        }
      }
      RFoptions(storing = FALSE);
      cat("****** '", .fct.list[.idx], "' (", .idx, ") done. ******\n")
      if (.tryok && !is.na(RFoptions()$general$seed)) {
        Print(.not_working, paste(.not_working_no, collapse=", ")) #
        stop("seed not NA: ", .fct.list[.idx])
      }
    }
  }
  Print(.not_working, paste(.not_working_no, collapse=", ")) #
  .ret <- list(.not_working, .not_working_no)
  names(.ret) <- c(.package, "")
  return(.ret)
}

StartExample <- function(reduced = TRUE) {
  if (!interactive()) {
### do not touch next lines
    ## REDUCED <- reduced, reduced
    REDUCED <- reduced
    RFoptions(examples_reduced = REDUCED)
  }
}

FinalizeExample <- function() {
  if (!interactive()) {
    close.screen(all.screens = TRUE)  
    n <- length(dev.list()) - 2 ## otherwise R CMD check complains
    ##                             for missing graphic device
    if (n > 0) {
      for (i in 1:n) {      
        dev.off() ## OK
      }
    }
  }

  
  RFoptions(seed = NA, coord_system="auto", examples_reduced=FALSE)
  
  #Print(RFoptions()$gen)
}



showManpages <- function(path="/home/schlather/svn/RandomFields/RandomFields/man") {
  files <- dir(path)
  result <- character(length(files))
  for (i in 1:length(files)) {
    cat("\n\n\n\n\n")
    system(paste("grep \"\\examples\" -A 10000 ", path, "/", files[i], sep=""))
    result[i] <- readline(files[i])
  }
  result <- cbind(files, result)
  result <- result[result[, 2]!= "", ]
  result[result[, 2]=="t", 2] <- "To Do"
  result[result[, 2]=="s", 2] <- "dontshow"
  return(result)
}
# showManpages()

compareVersions <- function(path = "~/svn/RandomFields/cran") {
  l <- list(list(subpath="RandomFields/R", pattern="*R"),
            list(subpath="RandomFields/src", pattern="*\\.h"),
            list(subpath = "RandomFields/src", pattern = "*\\.cc"))
  
  for (i in 1:length(l)) {
    subpath <- l[[i]]$subpath
    pattern <- l[[i]]$pattern
    files <- dir(path=subpath, pattern=pattern)
    
    for (f in files) {
      cmd <- paste("diff ",  path, "/", subpath, "/", f, " ",
                   subpath,"/", f, sep="")
      Print(cmd)  #
      system(cmd)
    }
  }
}

jpegORpdf <- function(graphics=NULL) {
  GD <- getOption("device")
  
  ispdf <- is.logical(all.equal(GD, "pdf"))
  isjpg <- is.logical(all.equal(GD, "jpeg"))
  

  if (!ispdf && !isjpg) {
    ispdf <- is.function(GD) && is.logical(all.equal(args(pdf), args(GD)))
    isjpg <- is.function(GD) && is.logical(all.equal(args(jpeg), args(GD)))
  } 
  
  if (is.null(graphics)) graphics <- RFoptions()$graphics        
  file <- graphics$file
  if (file != "") {
    f <- strsplit(file, "\\.")[[1]]
    ext <- tolower(f[length(f)])
    is_pdf <- ext == "pdf"
    is_jpg <- ext == "jpeg" || ext == "jpg"
    if (is_pdf || is_jpg) graphics$file <- paste(f[-length(f)], sep=".")
   
    if (ispdf || isjpg) {
      if ((is_pdf && isjpg) || (is_jpg && ispdf)) stop("'extension' does not fit given device")
    } else {
      ispdf <- is_pdf
      isjpg <- is_jpg
    } 
  }

  if (ispdf) GD <- pdf
  else if (isjpg) GD <- jpeg

  graphics$ispdf <- ispdf
  graphics$isjpg <- isjpg
  graphics$ext <- c(".pdf", ".jpg")[isjpg + 1]
  graphics$GD <- GD

#  Print(graphics)

  return(graphics)
  
}

ScreenDevice <- function(height, width) {
  graphics <- jpegORpdf()    
  if (graphics$ispdf || graphics$isjpg) {
    if (graphics$filenumber == 99999)
      stop("number of pictures exceeds max. number of pictures")        

    file <- graphics$file
    if (file != "") {
      if (!graphics$onefile) {
        file <- paste(file,
                      formatC(format="d", flag="0", graphics$filenumber,
                              width=5),
                      graphics$ext, sep="")
        RFoptions(graphics.filenumber = graphics$filenumber + 1)
      } else file <- paste(file, graphics$ext, sep="")
      ##            Print(file, dev.list())
      if (graphics$ispdf)
        graphics$GD(height=height, width=width, onefile=FALSE, file=file)
      else if (graphics$isjpg) graphics$GD(height=height, width=width,
                                           file=file, units="in", res=graphics$resolution)
      else stop("device not recognized by 'RandomFields'")
      return()
    } else {
      if (graphics$isjpg) stop("'file' in 'RFoptions' not set.")
    }
  }
  
  args <- names(as.list(args(graphics$GD)))    
  if (all(c("width", "height") %in% args) &&
      ( !any(c("file", "filename") %in% args)) || graphics$ispdf) {
    graphics$GD(height=height, width=width)        
    ##    Print("OK", height, width, RFoptions()$graphics, par()$cex)
    return()
    }
  
  if (RFoptions()$internal$warn_aspect_ratio) {
    RFoptions(warn_aspect_ratio = FALSE)
    cat("The graphical device does not seem to be a standard screen device. Hence the\naspect ratio might not be correct. (This message appears only one per session.)")
  }
}


#restore_par <- function(oldpar) {
#  do.call(graphics::par, oldpar)
#  graphics::par(cex = oldpar$cex) ## muss extra aufgerufen werden. noch mehr?
#}

ArrangeDevice <- function(graphics, figs, dh=2.8, h.outer=1.2,
                          dw = 2.5, w.outer=0.7) {
#  plot(1,1); return;  Print("ARRANGE", graphics)
  
  H <- curH <- graphics$height
  curW <- graphics$width
  if (is.finite(H) && H>0) {
    if (is.na(graphics$width) || graphics$width<=0.0) {
      H <- H * pmin(1, graphics$increase_upto[1] / figs[1],
                    graphics$increase_upto[2] / figs[2])
      DH <- H * dh / (dh + h.outer)
      HO <- H - DH
      curH <- figs[1] * DH + HO
      W <- H * (dw + w.outer) / (dh + h.outer)
      DW <- W * dw / (dw + w.outer)
      WO <- W - DW
      curW <-  figs[2] * DW + WO
    }

    graphics <- jpegORpdf(graphics)
    jpg_pdf <- graphics$isjpg || graphics$ispdf
    
   if (is.na(graphics$always_open_device)) {
      open <- interactive() || jpg_pdf
      RFoptions(graphics.always_open_device = open)
    } else open <- graphics$always_open_device

    if (length(dev.list()) > 0 && open) {      
      close <- graphics$always_close_device
      if (is.na(close)) close <- jpg_pdf
      if (close) dev.off()
    }
    
    if (open) {
      ScreenDevice(height=curH, width=curW)
    }
    
  } else warning("device could not be opened as 'height' is not positive.")
 
  if (graphics$split_screen) {
    if (length(dev.list()) > 0 && any(par()$mfcol != 1)) par(mfcol=c(1,1))
  } else {
    if (any(figs != 1) && length(dev.list()) > 0) par(mfcol=figs)
  }
  return(c(curH, curW))
}


reverse_dependencies_with_maintainers <-
  function(packages, which = c("Depends", "Imports", "LinkingTo"),
           recursive = FALSE) {
    ## function taken from CRAN developer website. 
    repos <- getOption("repos")["CRAN"]
    Print(repos) #
    ## if (substr(repos, 1, 1) == "@") repos <- "http://cran.r-project.org"
    Print(repos) #
    contrib.url(repos, "source") # trigger chooseCRANmirror() if required
    description <- sprintf("%s/web/packages/packages.rds", repos)
    con <- if(substring(description, 1L, 7L) == "file://")
      file(description, "rb")
    else
      url(description, "rb")
    on.exit(close(con))
    db <- readRDS(gzcon(con))
    rownames(db) <- NULL
    
    rdepends <- tools::package_dependencies(packages, db, which,
                                            recursive = recursive,
                                            reverse = TRUE)
    rdepends <- sort(unique(unlist(rdepends)))
    pos <- match(rdepends, db[, "Package"], nomatch = 0L)
    
    db[pos, c("Package", "Version", "Maintainer")]
  }

Dependencies <- function(install = all(pkgs == all.pkgs),
                         check=TRUE, pkgs = all.pkgs, dir = "Dependencies") {
  all <- reverse_dependencies_with_maintainers("RandomFields", which="all")
  all.pkgs <- all[, 1]
  PKGS <- paste(all[,1], "_", all[,2], ".tar.gz", sep="")    
  
  ## getOption("repos")["CRAN"]
  URL <- "http://cran.r-project.org/src/contrib/"
  if (install) {
    system(paste("mkdir ", dir))
    system(paste("rm ", dir, "/*tar.gz*", sep=""))
    for (i in 1:length(pkgs)) {
      cat("PACKAGE:", PKGS[i], ":", i, "out of ", length(pkgs),"\n")
      x <- system(paste("(cd ", dir, "; wget ", URL, PKGS[i], ")", sep=""))
      if (x != 0) stop(PKGS[i], "not downloadable")
    ## extended version see RandomFields V 3.0.51 or earlier     
    }
  }
 
  if (check) {
    tools::check_packages_in_dir(dir=dir)
    return(NULL)

    ## old:
    for (i in 1:length(pkgs)) {
      command <- paste("(cd ", dir, "; R CMD check --as-cran", PKGS[i],")")
      Print(command) #
      x <- system(command)
      if (x != 0) stop(PKGS[i], "failed")
    }
  }
}
# R Under development (unstable) (2014-12-09 r67142) -- "Unsuffered Consequences"

#Dependencies()
