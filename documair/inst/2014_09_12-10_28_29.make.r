#
# 14_01_21 14_01_28 14_01_29 14_05_20 14_05_21
# 14_06_10 14_06_14 14_06_15 14_06_16 14_06_24
# 14_06_25 14_08_04 14_08_10 14_08_11 14_08_25
# 14_08_26 14_08_27 14_08_29 14_09_01 14_09_04
# 14_09_08 14_09_12
#
# Running this script within its directory will launch the
# construction of several versions of /documair/ packages:
#
#   /documair/ the standard version for publication on the CRAN
#
#   /documair1/ the version generated without 'documair.which.txt'
#               file.
#
#   /documair2/ the version generated with a very limited general
#               description without any signature and also with
#               a very limited 'documair2.which.txt' file, also
#               silently build.
#
#   /documair3/ Version similar to /documair/ but with
#               two examples of C and Fortran functions.
#
#  For the last building to be effective, necessary tools to
#               compile 'C' and 'Fortran' functions must be
#               available to \pkg{R}.
# 
# To do so a "../tmp" directory will be created (after
#               deletion if already existing).
#
# All resulting pdf manuals and tar.gz files will be placed in the
#    directory "../resu" (possibly replacing existing ones).
#
#
# Some other facilities for debugging can be found in the script itself.
#
#
## parameters
#
nbad <- 3;
build9 <- TRUE; # use build9pkg instead of (prepare8pkg and compile8pkg)
#
#
## Loading /rbsa/ one way or the other
#
rbsapa <- TRUE;
rbsach <- "/home/jbdenis/attente.liana/inra/paquets/rbsa/pro/perso/";
rbsach <- "/home/jbdenis/bananier/inra/p/r/paquets/rbsa/pro/perso/";
if (rbsapa) {
  library(rbsa);
} else {
  rbsaso <- Sys.glob(paste0(rbsach,"*.code.r"));
  for (rbsafi in rbsaso) {
    source(rbsafi);
  }
}
fins <- c("",bc(nbad));
#
#
## preparing and the directories
#
# test directory
tedi <- "../tmp";
if (fidi9(tedi)=="d") { unlink(tedi,recursive=TRUE);}
dir.create(tedi);
for (fs in fins) {
  dir.create(paste0(tedi,"/documair",fs));
  #
  perper <- paste0(tedi,"/perso",fs);
  dir.create(perper);
  file.copy(paste0("documair",fs,".DESCRIPTION"),perper);
  file.copy("documair.package.Rd",paste0(perper,"/documair",fs,".package.Rd"));
  if (fs!="1") { file.copy(paste0("documair",fs,".which.txt"),perper);}
  #
  file.copy(Sys.glob("*.code.r"),perper);
  file.copy(Sys.glob("*.test.r"),perper);
  if (fs=="3") {
    file.copy(Sys.glob("sum*"),perper);
  } else {
    file.remove(paste0(perper,"/exterieur.code.r"));
  }
}
#
# result directory
redi <- "../resu";
if (fidi9(redi)!="d") { dir.create(redi);}
#
#
## sourcing the coding files
#     
source8file(Sys.glob("*.code.r"));
#
#     
## building the package(s)
#
quels <- bf(fins);
quels <- 1;
#
for (pkg in quels) {
  fs <- fins[pkg];
  paquet <- paste0("documair",fs);
  perdir <- paste0(tedi,"/perso",fs);
  pkgdir <- paste0(tedi,"/documair",fs)
  tra <- (pkg!=3);
  sig <- abs(3-pkg);
  form3title(paste(paquet,"started"),3);
  if (build9) {
    resu <- build8pkg(paquet,documair7dir=perdir,
                      destination7dir=redi,what="pzl");
browser()
    if (length(resu)>0) {
      print(resu);
      stop("Something Wrong");
    }
  }else {
    prepare8pkg(paquet,perdir=perdir,pkgdir=pkgdir,signature=sig,display=tra,check=FALSE);
    compile8pkg(paquet,pkgdir,chkdir=tedi,resdir=redi,display=tra);
  }
  #
  form3title(paste(paquet,"was built"),3);
}
#
form3title("*make.r* finished its job!",8);
