if (interactive()) {
  savehistory();
} else {
  # GLAD v1.12.0 depends on the 'tcltk' package, but that
  # cannot be loaded if there is no display.  A workaround
  # is to fake that the 'tcltk' package is loaded:
  attach(list(), name="package:tcltk");
}
library(aroma.affymetrix);
#source("init.R");

# Use special file cache for testing
options("R.cache::rootPath"="~/.Rcache,scratch");
options("R.cache::touchOnLoad"=TRUE);


args <- commandArgs(asValues=TRUE, excludeReserved=TRUE, excludeEnvVars=TRUE);
print(args);

if (isTRUE(args$devel)) {
  setOption(aromaSettings, "devel/dropRootPathTags", TRUE);
}
print(getOption(aromaSettings, "devel/dropRootPathTags"));


paths <- c();
allPaths <- c("testScripts/replication/chipTypes", 
              "testScripts/system/chipTypes");
for (path in allPaths) {
  path <- Arguments$getReadablePath(path, mustExist=TRUE);
  paths0 <- list.files(path=path, full.names=TRUE);
  paths <- c(paths, paths0);
}

..pathnames <- lapply(paths, FUN=list.files, pattern="[.]R$", full.names=TRUE);
names(..pathnames) <- basename(paths);
#..pathnames <- ..pathnames[names(..pathnames)];

..chipTypes <- c("Mapping10K_Xba142",
                 "Test3",
                 "HG-U133_Plus_2",
                 "Mapping50K_Hind240",
                 "Hs_PromPR_v02",
                 "Mapping250K_Nsp",
                 "Mapping250K_Sty",
                 "HuEx-1_0-st-v2",
                 "GenomeWideSNP_5",
                 "GenomeWideSNP_6",
                 "Cytogenetics_Array",
                 "MOUSEDIVm520650");


if (!is.null(args$chipTypes)) {
  cat("User-specified chip types:\n");
  ..chipTypes <- trim(unlist(strsplit(args$chipTypes, split=";")));
  print(..chipTypes);
}

if (isTRUE(args$reverse)) {
  cat("Reversing chip types:\n");
  ..chipTypes <- rev(..chipTypes);
  print(..chipTypes);
}

if (isTRUE(args$shuffle)) {
  cat("Reshuffling chip types:\n");
  ..chipTypes <- sample(..chipTypes);
  print(..chipTypes);
}

cat("Processing chip types:\n");
print(..chipTypes);


for (..chipType in ..chipTypes) {
  keep <- (names(..pathnames) == ..chipType);
  ..pathnamesT <- unlist(..pathnames[keep], use.names=FALSE);
  cat("PATHNAMES CHIPTYPE: \n");
  print(..pathnamesT);
  for (pathname in ..pathnamesT) {
    if (regexpr("hetero", pathname) != -1)
      next;
    if (regexpr("expectile", pathname) != -1)
      next;

    cat("** PATHNAME: ", pathname, "\n", sep="");
    tryCatch({
      source(pathname, echo=TRUE);
      cat("** PATHNAME DONE: ", pathname, "\n", sep="");
    }, error = function(ex) {
      cat("************************************************\n");
      cat("** ", rep(c(" ER", "ROR "), times=6), " **\n", sep="");
      print(ex);
      print(sessionInfo());
      cat("************************************************\n");
      cat("** PATHNAME FAILED: ", pathname, "\n", sep="");
    });
    
    rm(list=setdiff(ls(), c("..chipType", "..chipTypes", 
                    "pathname", "..pathnames", "..pathnamesT")));
    gc();
  }
}
