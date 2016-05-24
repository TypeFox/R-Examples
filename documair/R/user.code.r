
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
copy8objects <- function(nfiles,nobjects,file)
#TITLE generates a file comprising indicated functions
#DESCRIPTION
# Explores a series of R code files supposed to 
# comprise the needed functions with \pkg{tagging}
# and copy them into a new not already existing file.
#DETAILS
# Standard tags defined by \samp{documair0$tags} are used.
#KEYWORDS 
#INPUTS
#{nfiles} << Files where to look for the functions.>>
#{nobjects} << A \samp{character} indicating
#           the objects to pick up from the
#           indicated files.>>
#{file} << Name of the file to create.>>
#[INPUTS]
#VALUE
# \samp{character(0)} if everything went right,
# if not, some warning messages about found
# difficulties. 
#EXAMPLE
#REFERENCE
#SEE ALSO 
#CALLING
#COMMENT See Also \samp{code7objects4text6tags}
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 14_09_07
#REVISED 14_09_07
#--------------------------------------------
{
  # checking and initialization
  rrrobject9(nobjects,"character");
  rrrobject9(nfiles,"file");
  rrrobject9(file,"character",1);
  if (rrrfidi9(file)=="f") {
    stop(paste(file,"seems to be an existing file"));
  }
  if (rrrfidi9(file)=="d") {
    stop(paste(file,"seems to be an existing directory"));
  }
  # initialization
  res <- character(0);
  # getting all the objects
  lobj <- vector("list",0);
  for (fo in nfiles) {
    lobj <- c(lobj,code7objects4text6tags(readLines(fo),addbra=TRUE));
  }
  # looping onto the functions
  tfi <- character(0);
  for (fu in nobjects) {
    if (rrrbelong9(fu,names(lobj))) {
      tfi <- c(tfi,"",lobj[[fu]]);
    } else {
      res <- c(res,paste("object",fu,"was not found"));
    }
  }
  # writing the file
  writeLines(tfi,file);
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
change8names <- function(nfiles,nnames)
#TITLE changes object names in R script files
#DESCRIPTION
# Explores a series of R script files and modify them
# changing some objects names for new ones.\cr
# Changes must be unambiguous; they are performed sequentially
# according to the order of \samp{nnames} rows. Also, if some functions 
# have got common names used in comments or messages, they will be
# changed as well since the syntaxic analysis made in
# \samp{rrrtext3places8word} is very limited.
#DETAILS
# It is difficult to ensure that the modifications
# are perfect, but at least they can give a hand.\cr
# Each file is dealt in turn. When old and new names of a file
# coincide, the old file is first saved into a temporary
# file suffixed with the moment
# of the saving and only destroyed at the end. So if something
# went wrong, you can recover the initial file.
#KEYWORDS 
#INPUTS
#{nfiles} <<A two column \samp{character} matrix,
#           giving by row old and new names of files
#           to consider. Old and new names can be
#           identical.>>
#{nnames} <<A two column \samp{character} matrix,
#           giving by row old and new names of objects.>>
#[INPUTS]
#VALUE
# \samp{character(0)} if everything went right,
# if not, some warning messages.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 14_09_04
#REVISED 14_09_04
#--------------------------------------------
{
  # checking
  rrrobject9(nfiles,"character");
  if (length(nfiles)==2) { nfiles <- matrix(nfiles,1);}
  rrrobject9(nfiles,"matrix",speci=matrix(c(NA,2,NA,2),2),
          mensaje="'nfiles' must be a 2-columns matrix");
  rrrobject9(nnames,"character");
  if (length(nnames)==2) { nnames <- matrix(nnames,1);}
  rrrobject9(nnames,"matrix",speci=matrix(c(NA,2,NA,2),2),
          mensaje="'nnames' must be a 2-columns matrix");
  # initialization
  res <- character(0);
  # looping onto the set of files
  for (iff in rrrbc(nrow(nfiles))) {
    ofi <- nfiles[iff,1]; nfi <- nfiles[iff,2];
    # precautionary saving
    if (ofi == nfi) {
      sff <- paste0(ofi,".",rrrnow("s","red",rep("",3)),".txt");
    }
    # reading the file
    tt <- readLines(ofi);
    # transforming the text
    bref <- c(letters,LETTERS,as.character(0:9));
    aref <- bref;
    tt <- rrrtext3replace(tt,nnames,
                       b.rejected=bref,a.rejected=aref);
    # writting the file
    writeLines(tt,nfi);
    # deleting the precautionary file
    if (ofi == nfi) {
      if (unlink(sff)!=0) {
        res <- c(res,paste("unable to delete ",sff," file!"));
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
display8tags <- function(tags=documair0$tags,
                         what="values",imp=FALSE)
#TITLE returns the contents of a 'tags' object
#DESCRIPTION
# Typically, \samp{tags} is either \samp{documair0$tags} or
# \samp{documair0$tgs} which defines the tagging of the R code
# to introduce the documentation to be tackled by /documair/.\cr
# Such a function is provided to help the user in seeing the 
# present constants, and possibly to change their value.
#DETAILS
# More details on the possibilities offered with argument \samp{what}
# are obtained into \samp{documair0$tags$d}.
#KEYWORDS 
#INPUTS
#[INPUTS]
#{tags}<< List to be displayed.>>
#{what}    <<a \samp{character(1)} indicating what must be returned. The acceptable
#           values are \samp{names}, \samp{definitions}, \samp{usage},
#           \samp{content}, \samp{preparation}, \samp{title}  and \samp{values}.
#           For user convenience, the argument is lowercased before checked,
#           and only the first character is taken into consideration.>>
#{imp}<< indicates if an adapted printing of the results must be done.>>
#VALUE
# A character or a list according to \samp{what}
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE Make an equivalent function to change the tags values
#AUTHOR J.-B. Denis
#CREATED 14_05_22
#REVISED 14_05_22
#--------------------------------------------
{
  # checking
  rrrobject9(what,"character",1,mensaje="Awaiting for a character(1)");
  what <- tolower(what);
  what <- substr(what,1,1);
  rrrobject9(what,"character",1,c("d","v","n","u","t","c","p"),
          paste("Must be 'names', 'definitions',",
                "'content','preparation,'",
                "'usage', 'title' or 'values'"));
  # filling
  if (what == "n") {
    res <- names(tags$v);
    if (imp) {
      print(res);
    }
  } else {
    res <- vapply(tags$v,function(a){as.character(a[[what]])},FUN.VALUE="ee");
    if (imp) {
      for (uu in rrrbf(res)) {
        cat(names(res)[uu],":",res[uu],"\n");
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
prepare8pkg <- function(pkg,perdir,pkgdir,
                        signature=2,
                        display=TRUE,
                        check=FALSE)
#TITLE prepares a package
#DESCRIPTION
# Prepares all necessary directories and files for the building of an R package named
# \samp{pkg} from files made by the user and placed in the \samp{perdir} directory.
# The built arborescence will placed into the preexisting \samp{pkgdir} directory.\cr
# When an error occurs, it is not always very simple to know the origin of it. To get
# some clues, the user is suggested to switch the  display options on within the \samp{pkg.which.txt}
# file and putting the \samp{check} argument to \samp{TRUE}.
#DETAILS
# The behaviour of \samp{prepare8pkg} with respect to files and objects coded in them, is also
# driven by a possible \samp{pkg.which.txt} file where can be indicated which 
# objects are exported or hidden, which files contains a series of aliased objects,
# which files and/or objects must be displayed. See the general description of the package for 
# details.\cr
# When no such file is present, the default behavior is that all object are exported,
# no series of objects are aliased, the displaying is just listing the explored files
# and the objects they contain. Notice that the default displaying as well as the 
# \samp{which} file display options can be cancelled with the argument \samp{display}.
#KEYWORDS 
#INPUTS
#{pkg}<< Name of the package, associated files like \samp{pkg.DESCRIPTION}
# must exist into \samp{perdir} directory.
#  >>
#{perdir} << Directory where the prepared files have to be found.>>
#{pkgdir} << Directory where the building of the package has to be prepared.
#            It is supposed to have already been created, at least empty.
#            The result of the preparation of the package will be placed
#            within it with a directory having the package name. In
#            case, it is not empty, it is first completely cleaned of
#            its contents.>>
#[INPUTS]
#{signature} << \samp{numeric(1)}. Must additional information be included in the
#               \samp{pkg.package.Rd} file? When \samp{0} no; when \samp{1}
#               just the name of \pkg{documair}; when \samp{2} plus the
#               list of the masked functions.>>
#{display} <<Must the programmed debbuging displays be performed?>>
#{check} << Must checking be done by the user after each debugging display?
#           Checking means that the program rrrpause after each display
#           awaiting for an answer from the user to continue or stop.>>
#VALUE
# Nothing but the preparation is made (see the description section)
# with possible displays to and checks from the user.
#EXAMPLE
#REFERENCE
#SEE ALSO compile8pkg display8tags
#CALLING
#COMMENT
#FUTURE Introduce more sections and improve the existing ones.
#AUTHOR J.-B. Denis
#CREATED 14_01_03
#REVISED 14_08_27
#--------------------------------------------
{
  # constants
  fields <- names(documair0$tags$v)[c(2:10,12,15)];
  # checking arguments
  rrrobject9(display,"logical",1,mensaje="'display' argument not accepted");
  if (display) { rrrform3title("Preparation by /prepare8pkg/",2);}
  rrrobject9(pkg,"character",1,mensaje="Bad package name");
  rrrobject9(perdir,"directory",1,mensaje="Bad 'perdir' name");
  rrrobject9(pkgdir,"directory",1,mensaje="Bad 'pkgdir' name");
  rrrobject9(signature,"integer",1,c(0,2),mensaje="Bad 'signature' argument");
  rrrobject9(check,"logical",1,mensaje="'check' argument not accepted");
  # normalizing the directories
  perdir <- rrrdipa(perdir);
  pkgdir <- rrrdipa(pkgdir);
  # detecting the presence of data files
  fdata <- Sys.glob(file.path(perdir,paste(pkg,".*.rda",sep="")));
  # detecting the presence of 'C' and 'Fortran' files
  fcfil <- Sys.glob(file.path(perdir,"*.c"));
  fffil <- Sys.glob(file.path(perdir,"*.f"));
  fsrc <- c(fcfil,fffil);
  # preparing some constants
  didi <- c("man","R","tests","inst");
  if (length(fdata) > 0) { didi <- c(didi,"data");}
  if (length(fsrc) > 0) { didi <- c(didi,"src");}
  names(didi) <- didi;
  #
  bpkgdir <- file.path(pkgdir,pkg);
  mpkgdir <- file.path(bpkgdir,didi[1]);
  rpkgdir <- file.path(bpkgdir,didi[2]);
  tpkgdir <- file.path(bpkgdir,didi[3]);
  ipkgdir <- file.path(bpkgdir,didi[4]);
  # starting the package directory
  if(unlink(bpkgdir,recursive=TRUE)!=0) {
    stop(paste0("Not able to delete existing ",bpkgdir," directory"));
  }
  if (!dir.create(bpkgdir)) {
    stop(paste0("Not able to create ",bpkgdir," directory"));
  }
  for (rep in didi) {
    ndir <- file.path(bpkgdir,rep);
    if (!dir.create(ndir)) {
      stop(paste0("Not able to create ",ndir," directory"));
    }
  }
  # dealing with the possible which file
  whi <- list(des=paste0(pkg,".DESCRIPTION"),
              pre=paste0(pkg,".package.Rd"),
              c.ext=c("*.code.r","*.CODE.r","*.Code.r",
                     "*.code.R","*.CODE.R","*.Code.R"),
              t.ext=c("*.test.r","*.TEST.r","*.Test.r",
                     "*.test.R","*.TEST.R","*.Test.R"),
              hif=character(0),exf=character(0),
              hio=character(0),exo=character(0),
              kwt=character(0),
              dif="_ALL_",dio=character(0)
             );
  if (display) { rrrform3title("Looking for the which file",2);}
  fiwhi <- file.path(perdir,paste0(pkg,".which.txt"));
  if (rrrfidi9(fiwhi)=="f") {
    if (display) { rrrform3title("A which file was found",2);}
    # reading it
    cowhi <- rrrfile2list(fiwhi);
    # decoding, checking and attributing it
    if (!is.null(cowhi$DESCRIPTION))  {
      rrrobject9(cowhi$DESCRIPTION,"character",1);
      whi$des <- cowhi$DESCRIPTION;
    }
    if (!is.null(cowhi$PRESENTATION)) {
      rrrobject9(cowhi$PRESENTATION,"character",1);
      whi$pre <- cowhi$PRESENTATION;
    }
    if (!is.null(cowhi$C.EXTE))         {
      rrrobject9(cowhi$C.EXTE,"character");
      whi$c.ext <- cowhi$C.EXTE;
    }
    if (!is.null(cowhi$T.EXTE))         {
      rrrobject9(cowhi$T.EXTE,"character");
      whi$t.ext <- cowhi$T.EXTE;
    }
    if (!is.null(cowhi$HIDDEN.F))     {
      rrrobject9(cowhi$HIDDEN.F,"character");
      whi$hif <- cowhi$HIDDEN.F;
    }
    if (!is.null(cowhi$EXPORTED.F))   {
      rrrobject9(cowhi$EXPORTED.F,"character");
      whi$exf <- cowhi$EXPORTED.F;
    }
    if (!is.null(cowhi$ALIASED.F))   {
      rrrobject9(cowhi$ALIASED.F,"character");
      whi$alf <- cowhi$ALIASED.F;
    }
    if (!is.null(cowhi$HIDDEN.O))     {
      rrrobject9(cowhi$HIDDEN.O,"character");
      whi$hio <- cowhi$HIDDEN.O;
    }
    if (!is.null(cowhi$EXPORTED.O))   {
      rrrobject9(cowhi$EXPORTED.O,"character");
      whi$exo <- cowhi$EXPORTED.O;
    }
    if (!is.null(cowhi$DISPLAY.F))     {
      rrrobject9(cowhi$DISPLAY.F,"character");
      whi$dif <- cowhi$DISPLAY.F;
    }
    if (!is.null(cowhi$DISPLAY.O))     {
      rrrobject9(cowhi$DISPLAY.O,"character");
      whi$dio <- cowhi$DISPLAY.O;
    }
    if (!is.null(cowhi$KEYWORDS))     {
      rrrobject9(cowhi$KEYWORDS,"ncharacter");
      whi$kwt <- cowhi$KEYWORDS;
    }
  } else {
    if (display) { rrrform3title("No 'which' file found",2);}
  }
  # dealing with the description file
  if (display) { rrrform3title("Getting the DESCRIPTION file",2);}
  fides <- file.path(perdir,whi$des);
  if (rrrfidi9(fides)!="f") {
    stop(paste0("No '",fides,"' was found"));
  } else {
    if(!file.copy(fides,file.path(bpkgdir,"DESCRIPTION"))) {
      stop(paste0("It was not possible to copy ",fides," into ",bpkgdir));
    }
    despac <- analyse8description(fides);
  }
  # looking for the available code files
  if (display) {
    rrrform3title("Getting the available files",2);
    cat("     (( * )) for hidden objects\n");
    cat("     [[ * ]] for exported objects with a user made 'Rd' file\n");
    cat("     << * >> for exported objects with a /documair/ made 'Rd file\n");
  }
  ofi <- character(0);
  for (cexte in whi$c.ext) {
    ofi <- c(ofi,list.files(path=perdir,pattern=glob2rx(paste0("*",cexte))));
  }
  if (length(ofi) == 0) {
    stop(paste0("No function files found into ",perdir));
  }
  tfi <- character(0);
  for (texte in whi$t.ext) {
    tfi <- c(tfi,list.files(path=perdir,pattern=glob2rx(paste0("*",texte))));
  }
  # copying the possible data files
  for (fda in fdata) {
    # getting the data name from the rda file
    nda <- strsplit(fda,paste(pkg,".",sep=""),fixed=TRUE)[[1]];
    nda <- nda[length(nda)];
    nda <- strsplit(nda,".",fixed=TRUE)[[1]][1];
    # getting the necessary documentation file
    Rd1 <- file.path(perdir,paste(pkg,nda,"Rd",sep="."));
    if (!file.exists(Rd1)) {
      stop(paste(Rd1,"the documentation file to be associated to",fda,"is not present"));
    }
    # getting the destination files
    rad2 <- file.path(bpkgdir,paste0(didi["data"],"/",nda,".rda"));
    Rd2  <- file.path(mpkgdir,paste0(nda,".Rd"));
    if(!file.copy(fda,rad2)) {
      stop(paste0("It was not possible to copy ",fda," as ",rad2));
    }
    if(!file.copy(Rd1,Rd2)) {
      stop(paste0("It was not possible to copy ",Rd1," as ",Rd2));
    }
  }
  # copying the possible 'C' and 'Fortran'  files
  for (fex in fsrc) {
    # getting the original file
    fex1 <- strsplit(fex,"/",fixed=TRUE)[[1]];
    fex1 <- fex1[length(fex1)];
    # getting the destination file
    fex2 <- file.path(bpkgdir,paste0(didi["src"],"/",fex1));
    if(!file.copy(fex,fex2)) {
      stop(paste0("It was not possible to copy ",fex," as ",fex2));
    }
  }
  # copying the 'perso' directory in the 'inst' one.
  if (display) { rrrform3title("Copying 'perso' into 'inst'",2);}
  ## looking for complete set of files
  pfi <- list.files(path=perdir,pattern=glob2rx("*"));
  ## excluding some files (like *~)
  for (ii in rev(rrrbf(pfi))) {
    if (substr(pfi[ii],nchar(pfi[ii]),nchar(pfi[ii])) == "~") {
      pfi <- pfi[-ii];
    }
  }
  ## copying the remaining files if any
  for (ii in rrrbf(pfi)) {
    fifu <- file.path(perdir,pfi[ii]);
    if(!file.copy(fifu,file.path(ipkgdir,pfi[ii]))) {
      stop(paste0("It was not possible to copy ",
                  fifu," into ",ipkgdir));
    }
  }
  # looking at the available objects
  if (display) { rrrform3title("Getting and Reading the available objects",2);}
  nobjt <- nobje <- nobjud <- 0; maobj <- exobj <- character(0);
  for (of in ofi) {
    # type of display for this file
    if (whi$dif == "_ALL_") {
      ofdif <- TRUE;
    } else {
      ofdif <- (of %in% whi$dif);
    }
    if (display) { if (ofdif) {
      rrrform3title(paste0("   Code file: ",of),1);
    }}
    # type of this file exported, hidden or aliased?
    of.exp <- TRUE;
    if (length(whi$hif) > 0) {
      if ("_ALL_" %in% whi$hif) { of.exp <- FALSE;
      } else { if (of %in% whi$hif) { of.exp <- FALSE;}}
    }
    if (!of.exp) {
      if (length(whi$exf) > 0) {
        if ("_ALL_" %in% whi$exf) { of.exp <- TRUE;
        } else { if (of %in% whi$exf) { of.exp <- TRUE;}}
      }
    }
    of.typ <- "h";
    if (of.exp) {of.typ <- "e";}
    if (length(whi$alf) > 0) {
      if (of %in% whi$alf) { of.typ <- "a";}
    }
    # extracting the objects within this file
    fifi <- file.path(perdir,of);
    bibi <- rrrfile2text(fifi,clean=FALSE,comment=character(0));
    anan <- code7objects4text6tags(bibi);
    if (of.typ %in% c("e","h")) {
      # non aliased situation
      for (ico in rrrbf(anan)) {
        # analyzing the object code
        lyse <- parse8code(anan[[ico]]);
        icon <- lyse$name;
        # type of display for this object
        oodif <- FALSE;
        if (length(whi$dio)>0) {
          if ("_ALL_" %in% whi$dio) {
            oodif <- TRUE;
          } else {
            oodif <- (icon %in% whi$dio);
          }
        }
        # is this file exported of hidden?
        ico.exp <- of.exp;
        if (length(whi$hio) > 0) {
          if ("_ALL_" %in% whi$hio) { ico.exp <- FALSE;
          } else { if (icon %in% whi$hio) { ico.exp <- FALSE;}}
        }
        if (!ico.exp) {
          if (length(whi$exo) > 0) {
            if ("_ALL_" %in% whi$exo) { ico.exp <- TRUE;
            } else { if (icon %in% whi$exo) { ico.exp <- TRUE;}}
          }
        }
        if (ico.exp) {
          # exported object
          fi1rd <- file.path(perdir,paste0(icon,".Rd"));
          fi2rd <- file.path(bpkgdir,paste0(didi[1],"/",icon,".Rd"));
          if (rrrfidi9(fi1rd)=="f") {
            # user documented
            if(!file.copy(fi1rd,fi2rd)) {
              stop(paste0("It was not possible to copy ",fi1rd," into ",fi2rd));
            }
            if (display) { if (ofdif) {
              cat("[[",icon,"]] ",sep="");
            }}
            nobjud <- nobjud+1;
          } else {
            # /documair/ documented
            # adding elementary component as keyword
            kkww <- components9(icon);
            lyse$kwd <- unique(c(lyse$kwd,kkww));
            # translating the keywords
            lyse$kwd <- rrrtext3translate(lyse$kwd,whi$kwt);
            # removing the designated '_NO_'
            wuwu <- which(lyse$kwd=="_NO_");
            if (length(wuwu)>0) {lyse$kwd <- lyse$kwd[-wuwu];}
            # removing the possible alias
            lyse$ali <- NULL;
            # making and writing the Rd file
            if (display) { if (oodif) {
              write8lyse(lyse,fields,check);
            }}
            rdf <- make8rd(lyse);
            rrrtext2file(rdf,fi2rd);
            nobje <- 1+nobje;
            if (display) { if (ofdif) {
              cat("<<",icon,">> ",sep="");
            }}
          }
          exobj <- c(exobj,icon);
        } else {
          # hidden object
          if (display) { if (ofdif) {
            cat("((",icon,")) ",sep="");
          }}
          maobj <- c(maobj,icon);
          if (display) { if (oodif) {
            write8lyse(lyse,fields,check);
          }}
        }
        nobjt <- 1+nobjt;
      }
    }
    if (of.typ =="a") {
      # aliased situation
      #===================
      # analyzing the repetition of different component
      # to know if parentherized name must be added.
      compte <- c(des=0,exa=0,val=0,det=0)
      for (ico in rrrbf(anan)) {
        yse <- parse8code(anan[[ico]]);
        if (length(yse$des)>0) { compte["des"] <- compte["des"]+1;}
        if (length(yse$exa)>0) { compte["exa"] <- compte["exa"]+1;}
        if (length(yse$val)>0) { compte["val"] <- compte["val"]+1;}
        if (length(yse$det)>0) { compte["det"] <- compte["det"]+1;}
      }
      compte <- (compte>1);
      #===================
      for (ico in rrrbf(anan)) {
        # analyzing the object code
        llyse <- parse8code(anan[[ico]]);
        icon <- llyse$name;
        if (ico==1) {
          # type of display for this object
          oodif <- FALSE;
          if (length(whi$dio)>0) {
            if ("_ALL_" %in% whi$dio) {
                oodif <- TRUE;
            } else {
              oodif <- (icon %in% whi$dio);
            }
          }
          # this is the parent alias
          #----------------------
          lyse <- llyse;
          # aliased are exported objects
          fi1rd <- file.path(perdir,paste0(icon,".Rd"));
          fi2rd <- file.path(bpkgdir,paste0(didi[1],"/",icon,".Rd"));
          if (rrrfidi9(fi1rd)=="f") {
            # user documented
            if(!file.copy(fi1rd,fi2rd)) {
              stop(paste0("It was not possible to copy ",fi1rd," into ",fi2rd));
            }
            docu <- TRUE;
            if (display) { if (ofdif) {
              cat("[[",icon,"]] ",sep="");
            }}
            fi2rali <- fi2rd
            nobjud <- nobjud+1;
          } else {
            # /documair/ documented
            docu <- FALSE;
            # adding elementary component as keyword
            kkww <- components9(icon);
            lyse$kwd <- unique(c(lyse$kwd,kkww));
            # translating the keywords
            lyse$kwd <- rrrtext3translate(lyse$kwd,whi$kwt);
            # removing the designated '_NO_'
            wuwu <- which(lyse$kwd=="_NO_");
            if (length(wuwu)>0) {lyse$kwd <- lyse$kwd[-wuwu];}
            ## completing the description by adding the name
            while (rrrvoid9(lyse$des[1])&(length(lyse$des)>0)) { lyse$des <- lyse$des[-1];}
            if (length(lyse$des)>0) {
              if (compte["des"]) {
                lyse$des <- c(paste0(" (",lyse$name,")"),lyse$des);
              }
            }
            ## completing the examples by adding the name
            while (rrrvoid9(lyse$exa[1])&(length(lyse$exa)>0)) { lyse$exa <- lyse$exa[-1];}
            if (length(lyse$exa)>0) {
              if (compte["exa"]) {
                lyse$exa <- c("##",paste("##",lyse$name),lyse$exa);
              }
            }
            ## completing the value by adding the name
            while (rrrvoid9(lyse$val[1])&(length(lyse$val)>0)) { lyse$val <- lyse$val[-1];}
            if (length(lyse$val)>0) {
              if (compte["val"]) {
                lyse$val <- c(paste0(" (",lyse$name,")"),lyse$val);
              }
            }
            ## completing the details by adding the name
            if (length(lyse$det)>0) {
              if (compte["det"]) {
                lyse$det <- c(paste0(" (",lyse$name,")"),lyse$det);
              }
            }
            exobj <- c(exobj,icon);
            ## taking count of the first proposed alias (if any)
            ## as main name
            if (length(lyse$ali)>0) { 
              fi2rali <- file.path(bpkgdir,paste0(didi[1],"/",lyse$ali[1],".Rd"));
              lyse$ali <- unique(c(lyse$ali,lyse$name));
              lyse$name <- lyse$ali[1];
            } else {
              fi2rali <- fi2rd;
              lyse$ali <- lyse$name;
            }
            if (display) { if (ofdif) {
              cat("<<",icon,">> ",sep="");
            }}
          }
          # counting
          exobj <- c(exobj,llyse$name);
          nobje <- 1+nobje;
          nobjt <- 1+nobjt;
          #-------------------
        } else {
          # this is a child alias
          #----------------------
          pali <- lyse$name;
          # aliased are exported objects
          fi1rd <- file.path(perdir,paste0(icon,".Rd"));
          fi2rd <- file.path(bpkgdir,paste0(didi[1],"/",icon,".Rd"));
          if (rrrfidi9(fi1rd)=="f") {
            # user documented
            cat("The parent child is",pali,"\n");
            if (docu) {
              cat("its documentation was provided with a specific .Rd file\n");
            } else {
              cat("its documentation was automatically prepared from the code.\n");
            }
            stop(paste0("The .Rd file (",fi1rd,") was found for the child aliased object:",
                        icon,"! Documentation must be provided at the parent level.\n"));
            nobjud <- nobjud+1;
          } else {
            if (display) { if (ofdif) {
              if (docu) {
                cat("[[",icon,"]] ",sep="");
              } else {
                cat("<<",icon,">> ",sep="");
              }
            }}
            ## Adding the see also
            lyse$sal <- c(lyse$sal,llyse$sal);
            ## Adding the keywords
            kkww <- components9(llyse$name);
            lyse$kwd <- unique(c(lyse$kwd,kkww,llyse$kwd));
            # translating the keywords
            lyse$kwd <- rrrtext3translate(lyse$kwd,whi$kwt);
            # removing the designated '_NO_'
            wuwu <- which(lyse$kwd=="_NO_");
            if (length(wuwu)>0) {lyse$kwd <- lyse$kwd[-wuwu];}
            ## Adding the description
            while (rrrvoid9(llyse$des[1])&(length(llyse$des)>0)) { llyse$des <- llyse$des[-1];}
            if (length(llyse$des)>0) {
              if (compte["des"]) {
                lyse$des <- c(lyse$des,paste0("\\cr (",llyse$name,") "),llyse$des);
              } else {
                lyse$des <- c(lyse$des,llyse$des);
              }
            }
            ## Adding the examples and possibly introducing the name
            while (rrrvoid9(llyse$exa[1])&(length(llyse$exa)>0)) { llyse$exa <- llyse$exa[-1];}
            if (length(llyse$exa)>0) {
              if (compte["exa"]) {
                lyse$exa <- c(lyse$exa,"##",paste("##",llyse$name),llyse$exa);
              } else {
                lyse$exa <- c(lyse$exa,llyse$exa);
              }
            }
            ## Adding the value
            while (rrrvoid9(llyse$val[1])&(length(llyse$val)>0)) { llyse$val <- llyse$val[-1];}
            if (length(llyse$val)>0) {
              if (compte["val"]) {
                lyse$val <- c(lyse$val,paste0("\\cr (",llyse$name,") "),llyse$val);
              } else {
                lyse$val <- c(lyse$val,llyse$val);
              }
            }
            ## Adding the arguments
            lyse$in0 <- unique(c(lyse$in0,llyse$in0));
            lyse$in1 <- unique(c(lyse$in1,llyse$in1));
            ## Adding the usage
            lyse$usage <- c(lyse$usage,
                                llyse$usage);
            ## Adding the details
            if(length(lyse$det)==0) {care <- "";} else {care <- "\\cr";}
            if (length(llyse$det)>0) {
              if (compte["des"]) {
                lyse$det <- c(lyse$det,
                              paste0(care,"(",llyse$name,")"),
                              llyse$det);
              } else {
                lyse$det <- c(lyse$det,llyse$det);
              }
            }
            ## Adding the aliases
            lyse$ali <- unique(c(lyse$ali,llyse$ali,llyse$name));
            # counting
            exobj <- c(exobj,llyse$name);
            nobje <- 1+nobje;
            nobjt <- 1+nobjt;
          } # end of a child alias
        } # end of parent or child alias
        #-------------------
      }
      if (!docu) {
        # making and writing the Rd file
        rdf <- make8rd(lyse);
        rrrtext2file(rdf,fi2rali);
      }
      if (display) { if (oodif) {
        write8lyse(lyse,fields,check);
      }}
    } # of.typ == "a"
    #===================
    if (display) {
      if (ofdif) {
        cat("\n");
        if (check) { rrrpause("file OK?");}
      }
    }
    # copying the code files
    if(!file.copy(fifi,rpkgdir)) {
      stop(paste0("It was not possible to copy ",fifi," into ",rpkgdir));
    }
  } # end of 'of' loop
  #
  # copying the present test file
  for (tf in tfi) {
    fifi <- file.path(perdir,tf);
    fofo <- file.path(tpkgdir,paste0(tf,".R"));
    if (display) { rrrform3title(paste0("   Test file: ",tf),1);}
    if (display) { cat("\n");}
    if(!file.copy(fifi,fofo)) {
      cat("fifi =",fifi,"and fofo =",fofo,"\n");
      stop(paste0("It was not possible to copy ",tf," into ",tpkgdir));
    }
  }
  #
  if (display) {
    rrrform3title(paste0(nobjt, " objects were found in total"),1);
    rrrform3title(paste0(nobje, " of them will be exported to the users"),1);
    rrrform3title(paste0(nobjud," of them have been documented by hand"),1);
  }
  if (nobjt == 0) {
    stop(paste0("No '*.code.r' files found into ",perdir));
  }
  # dealing with the presentation file
  if (display) { rrrform3title("Getting the general presentation file",2);}
  fipres <- file.path(perdir,whi$pre);
  fipr   <- file.path(mpkgdir,paste0(pkg,"-package.Rd"));
  if (!file.exists(fipres)) {
    stop(paste0("No '",fipres,"' was found"));
  } else {
    if (signature >0 ) {
      # adding some information to the presentation file
      presentation <- rrrfile2text(fipres);
      presentation <- c(presentation,"\\section{Additional Information}{");
      presentation <- c(presentation,"  \\itemize{");
      presentation <- c(presentation,paste("    \\item This package was built with /documair/ package (version ",
                                           documair0$version$v,")",sep=""));
      if (signature > 1) {
        presentation <- c(presentation,paste("       on",rrrnow("d")));
        presentation <- c(presentation,"    \\item");
        presentation <- c(presentation,paste("  There are",nobjt,"object(s) in total."));
        if (nobjt != nobje) {
          presentation <- c(presentation,paste0(nobje," are exported"));
          presentation <- c(presentation,paste0(" and there exist ",nobjt-nobje," masked object(s): "));
          presentation <- c(presentation,paste0(paste(sort(maobj),collapse=", "),"."));
        }
        if (nobjud > 0) {
          presentation <- c(presentation,paste0("    \\item ",nobjud," object(s) where documented by hand.)"));
        }
        if (length(fcfil) > 0) {
          presentation <- c(presentation,paste0("    \\item There were ",length(fcfil)," C file(s)."));
        }
        if (length(fffil) > 0) {
          presentation <- c(presentation,paste0("    \\item There were ",length(fffil)," Fortran file(s)."));
        }
        presentation <- c(presentation,paste0("   \\item They were provided through ",length(ofi)," code file(s)."));
        if (length(tfi)>0) {
          presentation <- c(presentation,paste0("    \\item There were also ",length(tfi)," test file(s)."));
        }
      }
      presentation <- c(presentation,"  }");
      presentation <- c(presentation,"}");
      rrrtext2file(presentation,fipr);
      if (!file.exists(fipr)) {
        stop(paste0("The creation of ",fipr," have failed."));
      }
    } else {
      # the presentation file is copied as it is
      if (!file.copy(fipres,fipr)) {
        stop(paste0("It was not possible to copy ",fipres," into ",bpkgdir));
      }
    }
  }
  # preparing the NAMESPACE file
  if (display) { rrrform3title("making the NAMESPACE",2);}
  finsp <- file.path(bpkgdir,"NAMESPACE");
  sink(finsp);
  cat("####### EXPORTED OBJECTS ###\n");
  cat("export(",exobj[1]);
  for (ii in rrrbf(exobj[-1])) { cat(",\n       ",exobj[1+ii]);}
  cat(");\n");
  if (length(fsrc) > 0) {
    cat("####### EXPORTED OBJECTS ###\n");
    cat("useDynLib(",pkg,sep='');
    for (fex in fsrc) {
      # getting the original file
      fex1 <- strsplit(fex,"/",fixed=TRUE)[[1]];
      fex1 <- fex1[length(fex1)];
      fess <- strsplit(fex1,".",fixed=TRUE)[[1]];
      fff <- paste(fess[rrrbc(length(fess)-1)],sep=".");
      cat(",",fff,sep="");
    }
    cat(");\n");
  }
  sink();
  #
  # returning
  invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
compile8pkg <- function(pkg,pkgdir,
                        chkdir=getwd(),resdir=getwd(),
                        what="pz",
                        display=FALSE)
#TITLE compiles a package 
#DESCRIPTION
# From the standard directory containing the package sources, launches 
# the checking and building by \pkg{R} itself. Pdf manual and tar ball
# can be placed in another directory.
#DETAILS
#KEYWORDS 
#INPUTS
#{pkg}<< Name of the package.>>
#{pkgdir} << Directory where the building of the package has been prepared.
#            Probably by the function \samp{prepare8pkg}.>>
#[INPUTS]
#{chkdir} << Directory where to put the directory where "R" makes the checking.
#            Must exist.  >>
#{resdir} << Directory where to copy the resulting files. Must exist.>>
#{what} << \samp{character(1)}. Which resulting file to copy: \samp{p} 
#           for the pdf manual and \samp{z} for the tar.gz ball.>>
#{display} << Must indications about the process progress be displayed on the screen?>>
#VALUE
# Nothing but the check, the compilation, [the copy of the manual,
#         the tar.gz] are performed.
#EXAMPLE
#REFERENCE
#SEE ALSO prepare8pkg
#CALLING
#COMMENT
#FUTURE  Make \samp{display} more efficient when \samp{FALSE}.
#AUTHOR J.-B. Denis
#CREATED 14_01_22
#REVISED 14_09_01
#--------------------------------------------
{
  # initialization
  res <- character(0);
  # checking arguments
  if (display) { rrrform3title("Preparation by /compile8pkg/",2);}
  rrrobject9(pkg,"character",1,mensaje="Bad package name");
  rrrobject9(pkgdir,"directory",1,mensaje="Bad 'pkgdir' directory");
  rrrobject9(chkdir,"directory",1,mensaje="Bad 'chkdir' directory");
  rrrobject9(resdir,"directory",1,mensaje="Bad 'resdir' directory");
  rrrobject9( what,"character",1,mensaje="Bad 'what' option");
  rrrobject9(display,"logical",1,mensaje="Bad 'display' option");
  # completing the pkg directory
  pkgdir <- paste0(rrrdipa(pkgdir),"/",pkg);
  # checking the package
  if (display) { rrrform3title(paste0("Checking ",pkg),2);}
  system(paste0("R CMD check -o ",chkdir," ",pkgdir),intern=!display);
  # looking for the version number
  if (display) { rrrform3title("Looking for the version number",2);}
  pkgdirD <- file.path(pkgdir,"DESCRIPTION");
  vv <- rrrfile2text(pkgdirD);
  versi <- "";
  for (ii in rrrbf(vv)) {
    if (substr(vv[ii],1,8)=="Version:") {
      versi <- substr(vv[ii],9,nchar(vv[ii]));
      versi <- rrrform3crop(versi," ");
    }
  }
  if (versi == "") {
    rrrerreur(pkgdirD,"No version number found in that file",w=TRUE);
  }
  # specific checking for /documair/
  if (pkg=="documair") {
    if (documair0$version$v != versi) {
      rrrform3display(versi);
      rrrform3display(documair0$version$v);
      stop("The version constant is inconsistent for /documair/!");
    }
  }
  # Checking directory
  check7dir <- file.path(rrrdipa(chkdir),paste0(pkg,paste0(".Rcheck/")));
  # copying the documentation
  if (rrrbelong9("p",what,lower=TRUE)) {
    if (display) { rrrform3title("Copying the pdf manual",2);}
    fi1doc <- file.path(check7dir,paste0(pkg,"-manual.pdf"));
    fi2doc <- file.path(rrrdipa(resdir),paste0(pkg,"-manual.pdf"));
    if (file.exists(fi2doc)) { file.remove(fi2doc);}
    if (!file.copy(fi1doc,fi2doc)) {
      res <- c(res,"The almost certain reason for this failure is that",
                    "{R} was not able to create the pdf manual due to previous errors!",
                    paste0("It was not possible to copy ",fi1doc," as ",fi2doc));
    }
  }
  # making the tar.gz
  if (rrrbelong9("z",what,lower=TRUE)) {
    if (display) { rrrform3title("Making the tar ball file",2);}
    system(paste0("R CMD build --resave-data ",pkgdir),intern=!display);
    tgzaa <- file.path(rrrdipa(getwd()),paste0(pkg,"_",versi,".tar.gz"));
    tgzoo <- file.path(rrrdipa(resdir),paste0(pkg,"_",versi,".tar.gz"));
    if (file.exists(tgzoo)) {file.remove(tgzoo);}
    if (!file.rename(tgzaa,tgzoo)) {
      res <- c(res,"The almost certain reason for this failure is that",
                   "{R} was not able to create the tar.gz due to previous errors!",
                   paste0("It was not possible to copy ",tgzaa," as ",tgzoo));
    }
  }
  # deleting the checking directory
  if (!rrrbelong9("l",what,lower=TRUE)) {
    if (unlink(check7dir,recursive=TRUE,force=TRUE)!=0) {
      res <- c(res,paste0("It was not possible to destroy the directory '",check7dir,"'."));
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
build8pkg <- function(pkg,documair7dir,destination7dir,
                      what="pz",display=TRUE,
                      signature=2,check=FALSE)
#TITLE builds a package 
#DESCRIPTION
# From the documair directory containing at least the description file (\samp{pkg.DESCRIPTION}), 
# the presentation file (\samp{pkk.pacakge.Rd}) and the code files (\samp{*.code.r}) launches 
# the preparation, the checking and the building by \pkg{R} itself.\cr
# Resulting \samp{pdf} manual and \samp{pkg.tar.gz} files are placed into the indicated 
# destination directory. Some intermediate directory (especially the \samp{pkg.Rcheck} directory)
# are placed then deleted (in case of no error) into the within directory: so the possibility
# to create and remove them there is necessary. To look at the intermediate results, \samp{prepare8pkg} 
# and \samp{compile8pkg} must be directly used, or \samp{l} for \samp{left} must be included in
# the argument \samp{what}.
#DETAILS
# Whatever is the argument \samp{what}, the preparation of the arborescence is made.\cr
# Checking with fatal errors can be issued by the functions called by \samp{build8pkg}. If this
# occurs, it means that the check made by \samp{build8pkg} are not sufficient: please indicate us it. 
#KEYWORDS 
#INPUTS
#{pkg}<< Name of the package.>>
#{documair7dir} << Directory where are gathered all necessary files for \pkg{documair}.>>
#{destination7dir} << Directory where the resulting files have to be placed. If not
#                     existing will be created.>>
#[INPUTS]
#{what} << \samp{character(1)}. Indicates the actions to perform. It is lowercased
#          before decoding. Wwhen it comprises \samp{c} the check is performed; 
#          when it comprises \samp{z} the tar.gz is generated; 
#          when it comprises \samp{p} the pdf manual is generated;
#          when it comprises \samp{l} intermediate directories are left.>>
#{display} << Must indications about the process progress be displayed onto the screen?>>
#{signature} << \samp{numeric(1)}. Must additional information be included in the
#               \samp{pkg.package.Rd} file? When \samp{0} no; when \samp{1}
#               just the name of \pkg{documair}; when \samp{2} plus the
#               list of the masked functions.>>
#{check} << Must checking be done by the user after each debugging display?
#           Checking means that the program rrrpause after each display
#           awaiting for an answer from the user to continue or stop.>>
#VALUE
# \samp{character(0)} when all went right if not a \samp{character}
# reporting the motives of difficulties.
#EXAMPLE
#REFERENCE
#SEE ALSO prepare8pkg compile8pkg
#CALLING
#COMMENT
#FUTURE  Make \samp{display} more efficient when \samp{FALSE}.
#AUTHOR J.-B. Denis
#CREATED 14_09_01
#REVISED 14_09_01
#--------------------------------------------
{
  ## STARTING
  # initializing
  res <- character(0);
  # checking arguments and preparation
  res <- c(res,rrrobject9(display,"logical",1,mensaje="'display' must be TRUE or FALSE",fatal=FALSE));
  res <- c(res,rrrobject9(pkg,"character",1,mensaje="Bad package name",fatal=FALSE));
  res <- c(res,rrrobject9(documair7dir,"directory",1,mensaje="Bad 'documair7dir' directory",fatal=FALSE));
  ddir <- rrrobject9(destination7dir,"character",1,mensaje="Bad 'destination7dir' string");
  if (length(ddir)==0) {
    if (rrrfidi9(destination7dir)!="d") {
      if (!dir.create(destination7dir)) {
        res <- c(res,paste("Impossible to create",destination7dir,"as a directory"));
      }
    }
  }
  res <- c(res,rrrobject9(what,"character",1,mensaje="Bad 'what' option",fatal=FALSE));
  res <- c(res,rrrobject9(signature,"integer",1,c(0,2),mensaje="Bad 'signature' argument",fatal=FALSE));
  res <- c(res,rrrobject9(check,"logical",1,mensaje="'check' argument not accepted",fatal=FALSE));
  if (length(res)>0) { return(res);}
  what <- tolower(what);
  ## PREPARATION
  # making the directory where to place the arborescence
  prepdir <- paste0(getwd(),"/","preparation7dir.",rrrnow("d"));
  if (rrrfidi9(prepdir)=="d") {
    if (unlink(prepdir,recursive=TRUE,force=TRUE)!=0) {
      res <- c(paste("The directory",prepdir,"was not empty and could not be removed."));
      return(res);
    }
  }
  if (!dir.create(prepdir)) {
      res <- c(paste("It was not possible to create the directory",prepdir,"."));
      return(res);
  }
  # making the arborescence
  prepare8pkg(pkg,perdir=documair7dir,pkgdir=prepdir,
              signature=signature,display=display,check=check);
  # ELABORATION OF THE PACKAGE
  # building the package
  res <- c(res,compile8pkg(pkg,pkgdir=prepdir,resdir=destination7dir,
                           what=what,display=display));
  # REMOVING THE INTERMEDIATE DIRECTORY
  if (!rrrbelong9("l",what,lower=TRUE)) {
    if (unlink(prepdir,recursive=TRUE,force=TRUE)!=0) {
      res <- c(res,paste0("It was not possible to destroy the directory '",prepdir,"'."));
    }
  }
  # RETURNING
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
