
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


## meta-function that generates functions like RMplus, RMwhittle, which
## the user will use to generate explicit covariance models, i.e. objects
## of class 'RMmodels'

param.text.fct <- function(catheg, names, havedistr=TRUE, Const=NULL,
                           ismath=FALSE){
  ifHasArg <- paste("  if (hasArg('", names,
                    "') && !is.null(subst <- substitute(", names,
                    "))) \n", sep="")
  if (ismath && any(idx <- names == "new" & Const == NN2)) {
    ifHasArg[idx] <- "  if (!(hasArg('new') && !is.null(subst <- substitute(new)))) new <- UNREDUCED\n"
  }
  
  x <- paste(ifHasArg, "\t", catheg, "[['", names, "']] <- ", sep="")
  for (i in 1:length(names)) {
    if (!ismath && names[i] == "proj")
      x[i] <- paste(x[i], "CheckProj(proj, subst)", sep="")
    else if (length(Const) > 0 && Const[i] >= NN1)
      x[i] <- paste(x[i], "CheckChar(", names[i], ", subst, ",
                    NAMES_OF_NAMES[Const[i] - NN1 + 1], ", ",
                    havedistr, ")", sep="")
    else x[i] <- paste(x[i], "CheckArg(", names[i], ", subst, ",
                       havedistr, ")", sep="")
  }
  x
}


rfGenerateModels <- function(assigning,
                             RFpath = "~/R/RF/svn/RandomFields",
                             RMmodels.file = paste(RFpath, "RandomFields/R/RMmodels.R",
                               sep="/")
                             ) {
  
  # if file already exists, remove it.
  if (assigning && file.exists(RMmodels.file))
    file.remove(RMmodels.file)
  
  write(file = RMmodels.file, append = TRUE,
        "\n## This file has been created automatically by 'rfGenerateModels'.\n\n")

  ## defined constants
  diminf <- 999999
  
  # define empty strings
  empty <- paste(rep(" ", MAXCHAR), collapse="")
  empty2 <- paste(rep(" ", MAXCHAR), collapse="")
  # inialized attribute parameter
  
  nr <- GetCurrentNrOfModels(TRUE)
  vn <- nr * MAXVARIANTS

  # get attribute parameter
  A <- .C("GetAttr", nr=integer(vn), type=integer(vn), operator=integer(vn),
          monotone=integer(vn), finiterange=integer(vn),
          simpleArguments=integer(vn), internal=integer(vn),
          domains=integer(vn), isos=integer(vn),
          maxdim=integer(vn), vdim=integer(vn),
          includevariants= as.integer(TRUE),
          paramtype = integer(vn * MAXPARAM),
          va = integer(1),
          PACKAGE="RandomFields")
  va <- A$va
  dim(A$paramtype) <- c(MAXPARAM, vn)

 
  i <- 1
  while (i <= va) {
    step <- 1
    ## sequential steps for each model

    if (A$internal[i] && A$internal[i] != INTERN_SHOW) {
      cat(i, "internal", .C("GetModelName",as.integer(A$nr[i]),
                            name=empty, nick=empty2,
                            PACKAGE="RandomFields")$name,"\n")      
      i <- i + 1
      next
    }

    domains <-  A$domains[i]
    if (domains == PREVMODELD) domains <- c(XONLY, KERNEL)
    
    # get model name
    ret <- .C("GetModelName",as.integer(A$nr[i]),
              name=empty, nick=empty2, PACKAGE="RandomFields")
    nickname <- nick <- ret$nick
    if (A$internal[i]) nickname <- paste("i", nickname, sep="")
    ismath <- substr(nick, 1, 2) == "R."
    
    #cat(i, nick, "\n");

    type <- A$type[i]
    iso <- A$isos[i]
    while (i + step <= va  &&
           nick ==  .C("GetModelName",as.integer(A$nr[i + step]),
               name=empty, nick=empty2, PACKAGE="RandomFields")$nick) {
      cat("...variant added\n")
      type <- c(type, A$type[i + step])
      iso <- c(iso, A$isos[i + step])
      step <- step + 1
    }

   
    finiterange <- as.logical(A$finiterange)
    finiterange[A$finiterange < 0] <- NA
 
    ## get names of submodels
    subname.info<- .Call("GetSubNames", as.integer(A$nr[i]), PACKAGE="RandomFields")
    subnames <- subname.info[[1]]
    subintern <- subname.info[[2]]
    subnames.notintern <- subnames[!subintern]
     
    # get names of  parameters
    paramnames <- .Call("GetParameterNames", as.integer(A$nr[i]),
                        PACKAGE="RandomFields")
    internal <- which(paramnames == INTERNAL_PARAM)
    if (length(internal) > 0) paramnames <- paramnames[-internal]
    elmnt <- which(paramnames == "element")
    if (length(elmnt) > 0)  {
      stopifnot(length(elmnt) == 1)
      paramnames <- c(paramnames[-elmnt], "element")
    }
        
    par.intern <- paramnames %in% subnames
    
    if (any(par.intern)) stop(nick, ": subnames (",
                              paste(subnames, collapse=", "),
                              ") and parameter names (",
                              paste(paramnames, collapse=", "),
                              ") match.")
   
    ex.anysub <- length(subnames)>0
    ex.sub <- length(subnames.notintern)>0
     
    ex.par <- length(paramnames)>0
    ex.std <- ((nick != ZF_DOLLAR[1] && any(isVariogram(type))) ||
               nick %in% c("RMball", "RMsum", "RMconstant", "RMfixcov",
                           "RMcovariate")
               || nick == ZF_PLUS[1] || nick[1] == ZF_MULT[1]) &&
                 !(nick %in% c("RMtrafo",  "RMsine")) && !ismath

    std.variables <-
      if (nick %in%  c("RMconstant")) "var"
      else c("var", "scale", "Aniso", "proj")
    
  #    cat( std.variables, ex.std)
   #stopifnot(i < 50)
    
    
    cat(i, "\t", nickname, ",\t",
        paste(std.variables, collapse=", "), "\t",
        ex.std, "\t",
        paste(DOMAIN_NAMES[domains+1], collapse="; "), "\t",
        paste(type, collapse="/"), "\n", sep="")
    
    if(nick == ZF_DOLLAR[1]){ 
      text.fct.head <-
        paste(nick, " <- function(phi, var, scale, Aniso, proj, anisoT)")
    } else {
      text.fct.head <-
        paste(nickname,
              " <- function(",
              if (ex.sub) {
                paste(paste(subnames.notintern, collapse=", "), sep="")
              },
              if (ex.sub && (ex.par || ex.std)) ", ",
              if (ex.par) {
                paste(paste(paramnames, collapse=", "), sep="")
              },
              if (ex.par && ex.std) ", ",
              if (ex.std) paste(std.variables, collapse =", "), 
              ")",
              sep="")
    }

     
    if (ex.par) {
      par.body <- param.text.fct("par.model", paramnames,
                                 any(isVariogram(type)) || any(type==ShapeType),
                                 A$paramtype[1:length(paramnames), i], ismath)
      if (any(idx <- paramnames == 'envir'))
        warning(nickname, ": envir not internal")
#      if (any(idx))
#        par.body[idx] <-
 #         "par.model[['envir']] <- if (hasArg(envir)) envir else new.env()"
    } else par.body <- NULL
 
    text.fct.body <-
      paste("{\n  ",
            "cl <- match.call()",
            "\n  ",
            "submodels <- par.general <- par.model <- list() \n  ",
            ## get submodels
            if (ex.anysub) {
              paste("if (hasArg(", subnames, ")) submodels[['", subnames,
                    "']] <- ", subnames, sep="", collapse="\n  ")
            },
            if (ex.anysub) "\n  ",
            "\n",
            ## get model specific parameter
            if (ex.par) paste(par.body, collapse="\n"),
            if (ex.par) "\n  ",
            ## get general model parameter
            if (ex.std) {
              paste(param.text.fct("par.general", std.variables),
                    collapse="\n  ")
            },
            "\n  ",
             # create RMmodel object
            "model <- new('", ZF_MODEL, "', ",
            "call = cl, ",
            "name = ", "'", nick, "'", ", \n  \t\t",
            "submodels = submodels, ",   "\n  \t\t",
            "par.model = par.model, ",
            "par.general = par.general)",

            "\n  ",
            "return(model)\n}\n",
            sep=""
            )

    text.fct <- paste(text.fct.head, text.fct.body)
    
    # assign class 'RMmodelgenerator' (ZF_MODEL_FACTORY) and attributes like stationarity
    # to the function:

     text.assign.class <-
      paste(nickname, " <- new('", ZF_MODEL_FACTORY,              "',", "\n\t",
         ".Data = ",        nickname,                                  ",", "\n\t",
         "type = ", "c('", paste(TYPENAMES[type+1], collapse="', '"), "'),",
            "\n\t",
         "isotropy = ", "c('", paste(ISONAMES[iso+1], collapse="', '"), "'),",
            "\n\t",
         "domain = ", "c('", paste(DOMAIN_NAMES[domains+1], collapse="', '"),   "'),", "\n\t",
         "operator = ",     as.logical(A$operator[i]),             ",", "\n\t",
         "monotone = ",    "'", MONOTONE_NAMES[A$monotone[i] + 1 - MISMATCH],
                                                                  "',", "\n\t",
         "finiterange = ",  finiterange[i],          ",", "\n\t",
         "simpleArguments = ",  as.logical(A$simpleArguments[i]),  ",", "\n\t",
         "maxdim = ", if(A$maxdim[i]>diminf) Inf else A$maxdim[i], ",", "\n\t",
         "vdim = ",         A$vdim[i],                                  "\n\t",
         ")",
         sep="")
 
    text <- paste(text.fct, "\n", text.assign.class, "\n\n\n", sep="")

    if (nickname == "RMwhittle") {
      cat(text.assign.class)
     # stop("KKKK")
    }
  
    if (assigning) {
      #sink(file = RMmodels.file, append = TRUE, type='output')
      write(file = RMmodels.file, append = TRUE, text)
      #cat(text)
      #sink()
      #unlink(RMmodels.file)
    }
    i <- i + step 
  }  ## matches for (i in 1:nr) {
 


  # if help page to the function does not exist, throw warning
  if (length(as.character(help(nick))) == 0) {
    if (file.exists("/home/schlather/R/RF/RandomFields/R/rf.RXX")||
        file.exists("do.not.rm.this.file")) {
      if (!any(nick == c("list of exceptions"))) {
        warn <- paste("Warning: help page for '", nick,"' does not exist.",
                      sep="")
        cat(warn, "\n")
      } 
    }
  }
  invisible()
}


kind <- function(Zeilen, i, start, cont="", ignore=" ", endofname=" ") {
  if (substr(Zeilen[i], 1, nchar(start)) == start) {
    s <- substring(Zeilen[i], nchar(start) + 1)
    j <- 2
    while (j <= nchar(s) && substr(s, j, j) != endofname) j <- j + 1
    stopifnot(j <= nchar(s)) 
    name <- substr(s, 1, j -1)
    u <- strsplit(substring(s, j + 1), "//")[[1]]
    kommentar <- paste(u[-1], collapse = " ")
    
    RC <- nchar(kommentar) > 0 &&
      nchar(strsplit(kommentar, "RC")[[1]][1]) < nchar(kommentar)
    value <- u <- u[1]
    i <- i + 1
    if (any(cont != "")) {
      repeat {
        j <- nchar(u)
        if (ignore != "") while (substr(u, j, j) %in% ignore) j <- j - 1
        if (substr(u, j, j) %in% cont) {
          u <- strsplit(Zeilen[i], "//")[[1]][1]
          value <- paste(value, u)
          i <- i + 1
        } else break
      }
    }
    res <- list(name=name, value=value, RC=RC, i=i)
    return(res)
  } else return(NULL)
}

clean <- function(x) paste(strsplit(paste(strsplit(x, "\t")[[1]], collapse=""), " ")[[1]], collapse="")

CC <- function(x) {
  if (is.numeric(x)) {    
    Real <- TRUE
    Integer <- x == as.integer(x)
    y <- x
  } else {
    stopifnot(is.character(x))
    if (substr(x, 1, 1) =='"') return(x)
    warn <- options()$warn
    options(warn = -1)
    y <- try(as.numeric(x), silent=TRUE)
    options(warn = warn)
    if (Real <- !class(y) == "try-error") {
      Integer <- nchar(strsplit(x, "\\.")[[1]][1]) == nchar(x)
    } else Real <- Integer <- FALSE
    y <- paste(strsplit(paste(strsplit(x, "\t")[[1]], collapse=""), " ")[[1]],
               collapse="")
  }
  if (Real) {
    y <- paste("as.", if (Integer) "integer" else "double", "(", y, ")", sep="")
  }
  return(y)
}

rfGenerateConstants <-
  function(RFpath = "~/R/RF/svn/RandomFields",
           RCauto.file = paste(RFpath, "RandomFields/R/RCauto.R", sep="/")
           ) {
  s <- scan(paste(RFpath, "RandomFields/src/AutoRandomFields.h", sep="/"),
            what=character(), sep="\n", blank.lines.skip=FALSE, skip=2)
  #for (i in 1:length(s)) cat(s[i], "\n")
  i <- 1
  typedef <- character(0)
  write(file = RCauto.file,
        "# This file has been created automatically by 'rfGenerateConstants'")
  while (i <= length(s)) {    
    if (!is.null(k <- kind(s, i, "#define", "\\", ""))) {
      value <- k$value
      if (length(typedef) > 0)
        for (j in 1:length(typedef)) {
          value <- paste(strsplit(value, typedef[j])[[1]], collapse=" ")
        }
      v <- clean(k$name)
      w <- if (k$RC) paste("RC_", v, " <-", sep="") else ""
      line <- paste(w, v , "\t<-", CC(value))
      write(file = RCauto.file, append = TRUE, line)
    } else if (!is.null(k <- kind(s, i, "typedef enum", c(",", "{"), " "))) {
      value <- strsplit(k$value, "\\{")[[1]]
      typedef <- c(typedef,
                   paste("\\(", strsplit(value[1], ";")[[1]][1] ,"\\)", sep=""))
      value <- strsplit(strsplit(value[2], "\\}")[[1]][1], ",")[[1]]
      for (j in 1:length(value)) {
        v <- clean(value[j])
        w <- if (k$RC) paste("RC_", v, " <-", sep="") else ""
        line <- paste(w, v , "\t<-", CC(j - 1))
        write(file = RCauto.file, append = TRUE, line)
      }
      write(file = RCauto.file, append = TRUE, "")
    } else if (!is.null(k <- kind(s, i, "typedef", "\\", ""))) {
      typedef <- c(typedef,
                   paste("\\(", strsplit(k$value, ";")[[1]][1] ,"\\)", sep=""))
    } else if (!is.null(k <- kind(s, i, "extern const", ",", " "))) {
      ## ignored
    } else write(file = RCauto.file, append = TRUE, "")
    i <- if (is.null(k)) i+1 else k$i
  }

  
  s <- scan(paste(RFpath, "RandomFields/src/AutoRandomFields.cc", sep="/"),
            what=character(), sep="\n", blank.lines.skip=FALSE, skip=1)
  i <- 1
  while (i <= length(s)) {
    if (!is.null(k <- kind(s, i, "  *", c(",", "{"), " ", endofname="["))) {
      value <- strsplit((k$value), "\\{")[[1]]
      value <- strsplit(value[2], "\\}")[[1]][1]
      v <- clean(k$name)
      w <- if (k$RC) paste("RC_", v, " <-", sep="") else ""
      line <- paste(w, v, "<-\nc(", value, ")")
      write(file = RCauto.file, append = TRUE, line)
      write(file = RCauto.file, append = TRUE, "")
   } else write(file = RCauto.file, append = TRUE, "")
    i <- if (is.null(k)) i+1 else k$i
  }

  define_char <- function(name, value) {
    write(file = RCauto.file, append = TRUE,
        paste("\n", name, " <- c('", sep="",
              paste(value, collapse="', '"),
              "')")
          )
  }

  define_char("list2RMmodel_Names",
              c(RFgetModelNames(group.by=NULL), ZF_INTERNALMIXED, ZF_TREND))
 
  define_char("rfgui_Names1",

  ###        library(RandomFields, lib="~/TMP")
            
              RFgetModelNames(type=TYPENAMES[c(TcfType, PosDefType) + 1],
                              isotropy=ISONAMES[ISOTROPIC + 1],
                              operator=FALSE,
                              group.by=NULL,
                              valid.in.dim = 1,#if (sim_only1dim) 1 else 2,
                              simpleArguments = TRUE,
                              vdim=1)
              

              )

  
  define_char("rfgui_Names2",
              RFgetModelNames(type=TYPENAMES[c(TcfType, PosDefType) + 1],
                              isotropy=ISONAMES[ISOTROPIC + 1],
                              operator=FALSE,
                              group.by=NULL,
                              valid.in.dim = 2,#if (sim_only1dim) 1 else 2,
                              simpleArguments = TRUE,
                              vdim=1))
 
  return(NULL)
}


rfGenerateTest <- function(files = NULL,
                           RFpath = "~/R/RF/svn/RandomFields/RandomFields") {
   start.after <-  "\\dontrun"
   end.before <- c("\\", "}")
   initial.text <- "if (RFoptions()$internal$do_tests){"
   final.text <- "}"
   comment <- "%"

   if (length(files) == 0) return()
   ncomment <- nchar(comment)[1]
   nendbefore <- nchar(end.before)[1]
   for (f in 1:length(files)) {
    cat("creating ", f, ".R\n", sep="")
    s <- scan(paste(RFpath, "/man/", files[f], ".Rd", sep=""),
              what=character(), sep="\n", blank.lines.skip=FALSE, skip=2)
    i <- 1
    while (i <= length(s) && substr(s[i], 1,
               nchar(start.after)) != start.after) i <- i + 1 
    i <- i + 1
    if (i <= length(s))
      for (j in i:length(s)) {
        if (substr(s[j], 1, ncomment) %in% comment) s[j] <- ""
        if (substr(s[j], 1, nendbefore) %in% end.before) {
          j <- j -1
          break;
        }
      }
    if (j >= i) {
      out <- paste(RFpath, "/tests/", files[f], ".R", sep="")
      write(file = out, append = FALSE, initial.text)
      write(file = out, append = TRUE, s[i:j])
      write(file = out, append = TRUE, final.text)
    }
  }

  return(NULL)
}





rfGenerateMaths <- function(files = "/usr/include/tgmath.h",
                            ## copy also in ../private/lit
                            Cfile = "QMath",
                            RFpath = "~/R/RF/svn/RandomFields/RandomFields") {
  prefix <- "R."
  start.after <-  "/* Unary functions"
  end.before <- c("#define carg")
  initial.text <- "if (RFoptions()$internal$do_tests){"
  final.text <- "}"
  comment <- c("/*", "#i", "#e")
  
  if (length(files) == 0) return()
  ncomment <- nchar(comment)[1]
  nendbefore <- nchar(end.before)[1]
  cfile <-  paste(RFpath, "/src/", Cfile, ".cc", sep="")
  write(file = cfile,
        c("// This file has been created automatically by 'rfGenerateMaths'",
          "#include <math.h>",
          "#include \"RF.h\"",
          "#include \"primitive.h\""
          ))
  manfile <- paste(RFpath, "/man/", Cfile, ".Rd", sep="")
  write(file = manfile, 
        scan(paste(manfile, 0, sep="."),
             what=character(), sep="\n", blank.lines.skip=FALSE, skip=0))  

  Rfile <- paste(RFpath, "/R/", Cfile, ".R", sep="")
  write(file = Rfile, 
          "# This file has been created automatically by 'rfGenerateMaths'")
  
  usage <- include <- list()
  for (f in 1:length(files)) {
    cat("scanning ", f, files[f],"\n")
    s <- scan(files[f], what=character(), sep="\n", blank.lines.skip=FALSE,
              skip=0)
    i <- 1
    while (i <= length(s) && substr(s[i], 1, nchar(start.after)) !=start.after){
      i <- i + 1
    }
    i <- i + 1
 
    usage[[f]] <- include[[f]] <- rep("", length(s))
    while (i <= length(s)) {

       if (substr(s[i], 1, nendbefore) %in% end.before) break
      if (!is.null(k <- kind(s, i, "#define", "\\", " ", ")"))) {
        x <- strsplit(k$name, "\\(")
        name <- clean(x[[1]][1])
      
       if (!(name %in% c("frexp", "ldexp", "remquo", "scalbn", "scalbln",
                          "ilogb", "fma"))) { 
          args <- length(strsplit(x[[1]][2], ",")[[1]])
                                        #cat(name, args, k$name, "\n")
          write(file = cfile, append = TRUE,
                paste("void Math", name,
                      "(double *x, cov_model *cov, double *v){",
                      "\nMATH_DEFAULT\n",
                      "*v = ", name, "(w[0]",                    
                      if (args == 2) ", w[1]",
                      "); \n}\n\n", sep=""))
          
          include[[f]][i] <- paste(
              'IncludeModel(".', name, '", MathDefinition, 0, 0, ', args,
              ', NULL, XONLY', ## auch fuer kernel schreiben??
              ',\n\t PREVMODELI,checkMath,rangeMath, PREF_TREND,\n\t',
              'false, SCALAR, PREVMODEL_DEP, false, false); \n',
              'nickname("', name, '");\n',
              'kappanames("a", REALSXP',
              if (args == 2) ', "b", REALSXP', ');\n',
              'addCov(Math', name, ', NULL, NULL);\n',
              'AddVariant(TrendType, PREVMODELI);\n', sep="")
          
          write(file = manfile, append = TRUE,
                paste("\\alias{", prefix, name, "}", sep=""))
          usage[[f]][i] <- paste(prefix, name, "(a", if (args == 2) ", b",
                                 ")", sep="")

          if (name %in% c(
"asin",
"atan",
"atan2",
"cos",
"sin",
"tan",
"acosh",
"asinh",
"atanh",
"cosh",
"sinh",
"tanh",
"exp",
"log",
"expm1",
"log1p",
"logb",
#"exp2",
"log2",
#"pow",
"sqrt",
#"hypot",
#"cbrt",
#"ceil",
"fabs", # abs
"floor",
#"fmod",
#"nearbyint",
"round",
"trunc",
#"lrint",
#"llrint",
#"lround",
#"llround",
#"copysign",
#"erf",
#"erfc",
"tgamma", # gamma
"lgamma",
#"rint",
#"nextafter",
#"nexttoward",
#"remainder",
#"fdim",
"fmax",
"fmin")) {
           Rname <- if (name == "fabs") "abs" else
           if (name == "tgamma") "gamma" else
           if (name == "fmin") "min" else
           if (name == "fmax") "max" else name
           
           write(file = manfile, append = TRUE,
                paste("\\alias{", Rname, "}", sep=""))
           args <- (if (Rname == "atan2") "y, x"
                    else if (Rname %in% c("min", "max")) "..."
                    else if (Rname %in% c("round")) "x, ..."
                    else "x")
           fstarg <- (if (Rname == "atan2") "y"
                    else if (Rname %in% c("min", "max")) "list(...)[[1]]"
                    else "x")
          usage[[f]][i] <- paste(Rname, "(", args, ")\n",
                                 usage[[f]][i], sep="")

           write(file = Rfile, append=TRUE,
                paste(Rname,
                      " <- function(", args, ") if (is(", fstarg,
                      ", ZF_MODEL)) ", prefix, name, 
                      "(", args, ") else ", "base::", Rname,"(",
                      args, ")", sep=""))
         }
        } # !name in
      } # !is.null k
      i <- if (is.null(k)) i+1 else k$i
    }
  }

  write(file = cfile, append=TRUE, "void includeStandardMath() {")
  write(file = manfile, append=TRUE,
        scan(paste(manfile, 1, sep="."),
             what=character(), sep="\n", blank.lines.skip=FALSE, skip=0))  
  for (f in 1:length(files)) {
    write(file = cfile, append=TRUE, include[[f]][include[[f]] != ""]);
    write(file = manfile, append=TRUE, usage[[f]][usage[[f]] != ""]);
  }
  write(file = cfile, append=TRUE, "}")
  write(file = manfile, append=TRUE,
        scan(paste(manfile, 2, sep="."),
             what=character(), sep="\n", blank.lines.skip=FALSE, skip=0))
  
  return(NULL)
 }
