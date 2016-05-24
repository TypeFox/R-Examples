##' Read Mplus output files
##'
##' @title Read Mplus output
##' @param infile Mplus output file
##' @param coef Coefficients only
##' @param \dots additional arguments to lower level functions
##' @author Klaus K. Holst
##' @export
##' @seealso getSAS
`getMplus` <-
function(infile="template.out",coef=TRUE,...) {
##  mycmd <- paste("grep -n \"Estimates     S.E.  Est./S.E.\" | cut -f1 -d:", outfile)

  if (coef) {
    start <- "MODEL RESULTS"
    end1 <- "R-SQUARE"
    res0 <- findfileseg(infile,start)[-c(seq(5))]
    res <- sapply(res0,function(x) { val <- strsplit(x," ")[[1]]; val[val!=""] })
    res <- res[unlist(lapply(res, length))!=0]

    coef.idx <- unlist(lapply(res, length))>3
    lab.idx <- which(!coef.idx)
    count <- 0
    mycoef <- c()
    myrownames <- c()
    for (i in seq_along(res)) {
      if (i %in% lab.idx) {
        count <- count+1
      } else {
        val <- as.numeric(res[[i]][-1])
        if (length(val)<5) val <- c(val,rep(0,5-length(val)))
        mycoef <- rbind(mycoef, val)
        myrownames <- c(myrownames,
                        paste(paste(res[[lab.idx[count]]],collapse=" "),res[[i]][1])
                        )

      }
    }
    rownames(mycoef) <- myrownames
    colnames(mycoef) <- c("Estimate","Std.Err","Z-value","Std","StdYX")
    return(mycoef)
  }

  start <- "Estimate       S.E.  Est./S.E."
  end1 <- "MODEL RESULTS"
  end2 <- "QUALITY OF NUMERICAL RESULTS"
##  start <- "Estimate       S.E.  Est./S.E."
##  end1 <- "Beginning Time:"
 ## end2 <- "TECHNICAL"
  res <- findfileseg(infile,start,end1);
##  res2 <- findfileseg(infile,start,end2);
##  if (length(res)>length(res2))
##    res <- res2
  cat(paste(res,"\n"))

  res <- findfileseg(infile, "TESTS OF MODEL FIT", "Chi-Square Test of Model Fit for the Baseline Model")
  cat(paste(res,"\n"))
}



`findfileseg` <-
function(infile, startstring, endstring,nlines) {
  con <- file(infile, blocking = FALSE)
  inp <- readLines(con)
  close(con)
  nullstring <- 0
  linestart <- 1; lineend <- length(inp)

  mycmd1 <- paste0("grep -n \"",startstring,"\" ", infile);  a1 <- system(mycmd1,intern=TRUE);
  if (length(a1)>0)
    linestart <- as.numeric(strsplit(a1,":")[[1]][1])

  nn <- length(inp)
  if (!missing(nlines)) nn <- linestart+nlines
  if (missing(endstring)) {
    for (i in seq(linestart,nn)) {
      lineend <- i-1
      if (inp[i]==inp[i-1]) break;
    }
  } else {
    mycmd2 <- paste0("grep -n \"",endstring,"\" ", infile);  a2 <- system(mycmd2,intern=TRUE);
    if (length(a2)>0)
      lineend <- as.numeric(strsplit(a2,":")[[1]][1])
  }

  res <- inp[linestart:lineend-1]
  return(res)
}


##################################################
### Generate code and run mplus...
##################################################

`mplus` <-
function(file="template.mplus",wait=TRUE,intern=TRUE,...) {
    if (!file.exists(file)) file <- paste0(file,".mplus")
    if (!file.exists(file)) stop("File does not exist")
    if (!exists("winecmd")) winecmd <- "wine"
    if (!exists("mplus.directory")) mplus.directory <- ""
    mycmd <- paste0(winecmd, " \"", mplus.directory, "mplus.exe\" ", file)
    system(mycmd, wait=wait, intern=TRUE)
    prefix <- strsplit(file, ".", fixed=TRUE)[[1]][1]
    return(getMplus(paste0(prefix,".out"),coef=TRUE))
}

`toMplus.data.frame` <-
function(x, datafile="data.tab",
         mplusfile="template.mplus",
         na.string=".", model="!f1 by x1;",
         analysis=NULL,
         categorical=NULL,
         group,
         run=FALSE, techout=FALSE,missing=TRUE,...) {
  write.table(x, file=datafile, sep="\t",
              quote=FALSE, row.names=FALSE, col.names=FALSE, na=na.string)
  varnames <- c()
  ngroups <- ceiling(ncol(x)/4)

  for (i in seq_len(ngroups)) {
    newline <- c("\t",colnames(x)[((i-1)*4+1):min(ncol(x), (i*4))],"\n")
    varnames <- c(varnames, newline)
  }

###   mplusfilesummary <- paste0("summary",mplusfile)
###   zz <- file(mplusfilesummary, "w")  # open an output file connection
###   cat(file=zz, "TITLE:  Summary-statistics\n")
###   cat(file=zz, "!-----------------------------------------------------\n")
###   cat(file=zz,"DATA:\n\tFILE=\"", datafile, "\";\n")
###   cat(file=zz,"VARIABLE:\n\tNAMES ARE\n")
###   cat(file=zz, varnames, ";\n\n")

###   cat(file=zz, "!-----------------------------------------------------\n")
###   cat(file=zz, "USEVARIABLES=\n!?;\n")
###   cat(file=zz, "!CATEGORICAL=?;\n")
###   cat(file=zz, "!MISSING=?;\n")
###   cat(file=zz, "!IDVARIABLE=?;\n")
###   cat(file=zz, "!-----------------------------------------------------\n")
###   cat(file=zz, "ANALYSIS:\n\tTYPE IS BASIC;\n")
###   cat(file=zz, "!-----------------------------------------------------\n")
###   cat(file=zz, "OUTPUT:\t\tstandardized sampstat;")
###   close(zz)

  zz <- file(mplusfile, "w")  # open an output file connection
  cat(file=zz, "TITLE: ...\n")
  cat(file=zz, "!-----------------------------------------------------\n")
  cat(file=zz,"DATA:\n\tFILE=\"", datafile, "\";\n")
  cat(file=zz,"VARIABLE:\n\tNAMES ARE\n")
  cat(file=zz, varnames, ";\n")
  if (!missing(group)) {
    groups <- unique(x[,group])
    mygroupdef <- paste("(",paste(groups,groups,sep="=",collapse=","),")")
    cat(file=zz, "GROUPING IS ", group, mygroupdef, ";\n", sep="")
  } else {
    cat(file=zz, "!GROUPING IS g (1=male, 2=female);\n")
  }
  cat(file=zz, "USEVARIABLES=\n", varnames,";\n")
  if (!is.null(categorical))
    cat(file=zz, paste("CATEGORICAL=",paste(categorical,collapse=" "),";\n"))
  cat(file=zz, "MISSING=",na.string,";\n",sep="")
  cat(file=zz, "!IDVARIABLE=?;\n")
  cat(file=zz, "!DEFINE: define new variables here;\n")
  cat(file=zz, "!SAVEDATA: save data and/or results;\n\n")
  if (is.null(analysis)) {
    cat(file=zz, "ANALYSIS: TYPE=MEANSTRUCTURE");
    if (missing) cat(file=zz, " MISSING;\n")
    else cat(file=zz,";\n")
    cat(file=zz, "ESTIMATOR=ML;\n")
    cat(file=zz, "INFORMATION=EXPECTED;\n")
    cat(file=zz, "ITERATIONS=5000;\n")
    cat(file=zz, "CONVERGENCE=0.00005;\n\n")
  } else {
    cat(file=zz,"ANALYSIS:\n")
    cat(file=zz, analysis,"\n")
  }
  cat(file=zz, "!-----------------------------------------------------\n")
  cat(file=zz, "MODEL:\n")
  cat(file=zz, model,"\n")
  cat(file=zz, "!-----------------------------------------------------\n")
  if (!techout)
    cat(file=zz, "OUTPUT: STANDARDIZED;\n")
  else
    cat(file=zz, "OUTPUT: MODINDICES(0); TECH1; TECH2; TECH5; STANDARDIZED;\n")
  cat(file=zz, "!\tSAMPSTAT;RESIDUAL;CINTERVAL;MODINDICES(0);\n")
  cat(file=zz, "!Other output options are:\n")
  cat(file=zz, "!\tSTANDARDIZED;     !Standardized coefficients\n")
  cat(file=zz, "!\tH1SE;             !Standard errors for the H1 model\n")
  cat(file=zz, "!\tH1TECH3;          !Estimated covar,corr matrices for par. estimates\n")
  cat(file=zz, "!\tPATTERNS;         !Summary of missing data patterns\n")
  cat(file=zz, "!\tFSCOEFFICIENT;    !Factor score coefficients and posterior covar matrix\n")
  cat(file=zz, "!\tFSDERTERMINACY;   !Factor score determinacy for each factor\n")
  cat(file=zz, "!\tTECH1;            !Parameter specifications and starting values\n")
  cat(file=zz, "!\tTECH2             !Parameter derivatives;\n")
  cat(file=zz, "!\tTECH3;            !Covar and Corr matrices for estimates\n")
  cat(file=zz, "!\tTECH4;            !Estimated means and covar for the latent variables\n")
  cat(file=zz, "!\tTECH5;            !Optimization matrix\n")
  cat(file=zz, "!\tTECH6;            !Optimization for categorical variables\n")
  cat(file=zz, "!\tTECH7;            !output for type Mixture\n")
  cat(file=zz, "!\tTECH8;            !Output for type mixture\n")
  cat(file=zz, "!\tTECH9;            !Error messages for MC study\n")
  cat(file=zz, "!\tMONTECARLO:              File is\n")
  close(zz)

  if (run & exists("mplus")) {
    res <- mplus(mplusfile)
    outfile <- paste0(strsplit(mplusfile,".",fixed=TRUE)[[1]][1],".out")
    getMplus(outfile)
    return(res)
  }
}

`toMplus.lvmfit` <-
function(x, model=NULL, data=model.frame(x), run=TRUE, categorical=NULL,##binary(Model(x)),
         mplusfile="template.mplus", ...) {
  mymodel <- ""
  M <- index(x)$M
  P <- index(x)$P
  nn <- vars(x)
  p <- length(nn)
  lat.var <- latent(x)
  lat.idx <- match(lat.var, vars(x))
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (M[i,j]!=0) {
        var1 <- nn[i]; var2 <- nn[j];
        if (i %in% lat.idx & !(j %in% lat.idx)) {## & !(j %in% lat.idx)) {
          key <- " on "
          mymodel <- paste0(mymodel, "\n", var1, " by ", var2, ";")
        } else {
          mymodel <- paste0(mymodel, "\n", var2, " on ", var1, ";")
        }
      }
    }
  }
  for (i in seq_len(p-1)) {
    for (j in ((i+1):p)) {
      if (P[i,j]!=0) {
        var1 <- nn[i]; var2 <- nn[j];
        mymodel <- paste0(mymodel, "\n", var1, " with ", var2, ";")
      }
    }
  }
  if (is.null(model))
    model <- mymodel
  mydata <- subset(as.data.frame(data), select=setdiff(nn,lat.var))
  toMplus.data.frame(mydata,model=mymodel,run=run, mplusfile=mplusfile, ...)
}
