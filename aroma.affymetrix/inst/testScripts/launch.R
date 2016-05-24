############################################################################
# EXAMPLES:
#
# // Run all system test sets for the Test3 chip type
# Rscript testScripts/launch.R --pattern=Test3
#
# // Run a random system (default) test set
# Rscript testScripts/launch.R --order=random --nbrOfSets 1
#
# // Run a random replication test set
# Rscript testScripts/launch.R --group=replication --order=random --nbrOfSets 1
#
# // Run one of the complete analyses
# Rscript testScripts/launch.R --group=complete --pattern=GSE12702
############################################################################
library("R.utils");

cat("==========================================================\n");
cat("BEGIN OF SESSION:\n");

cat("Command line arguments:\n");
args <- cmdArgs();
print(args);

printf("Hostname: %s\n", System$getHostname());
printf("Username: %s\n", System$getUsername());
printf("Current directory: %s\n", getwd());

cat("R system environment variables:\n");
sysenv <- c("R_LIBS", "R_LIBS_USER", "R_LIBS_SITE");
sysenv <- sapply(sysenv, FUN=Sys.getenv, USE.NAMES=TRUE, simplify=FALSE);
str(sysenv);

cat(".libPaths():\n");
print(.libPaths());

print(sessionInfo());

cat("Memory statistics:\n");
print(gc());
if (.Platform$OS.type == "windows") {
  printf("Current memory usage: %g MB\n", memory.size(max=FALSE));
  printf("Maximum memory usage: %g MB\n", memory.size(max=TRUE));
  printf("Maximum memory limit: %g MB\n", memory.size(max=NA));
}
cat("==========================================================\n");


# Load aroma.affymetrix (in a fault-tolerant way)
kk <- 1L;
while (kk < 10L) {
  printf("#%02d. Trying to load aroma.affymetrix...\n", kk);
  tryCatch({
    library("aroma.affymetrix");
    break;
  }, error = function(ex) {
    print(traceback());
    print(ex);
    # Sleep for a while and try again
    Sys.sleep(10);
    FALSE;
  });
  kk <- kk + 1L;
} # while()
if (kk >= 10L) throw("Failed to load aroma.affymetrix.");

path <- Arguments$getReadablePath("testScripts/R", mustExist=FALSE);
if (!isDirectory(path)) {
  path <- system.file(package="aroma.affymetrix");
  path <- file.path(path, "testScripts", "R");
}
path <- Arguments$getReadablePath(path);

pathname <- file.path(path, "launchUtils.R");
pathname <- Arguments$getReadablePathname(pathname);

source(pathname);

do.call(launchTestGroups, args);


cat("==========================================================\n");
cat("END OF SESSION:\n");
# Override default settings with command line arguments
args <- commandArgs(asValues=TRUE);
print(args);

printf("Hostname: %s\n", System$getHostname());
printf("Username: %s\n", System$getUsername());
printf("Current directory: %s\n", getwd());

cat("R system environment variables:\n");
sysenv <- c("R_LIBS", "R_LIBS_USER", "R_LIBS_SITE");
sysenv <- sapply(sysenv, FUN=Sys.getenv, USE.NAMES=TRUE, simplify=FALSE);
str(sysenv);

cat(".libPaths():\n");
print(.libPaths());

print(sessionInfo());

cat("Memory statistics:\n");
print(gc());
if (.Platform$OS.type == "windows") {
  printf("Current memory usage: %g MB\n", memory.size(max=FALSE));
  printf("Maximum memory usage: %g MB\n", memory.size(max=TRUE));
  printf("Maximum memory limit: %g MB\n", memory.size(max=NA));
}
cat("==========================================================\n");


############################################################################
# HISTORY:
# 2014-01-28
# o Now the testScripts/launch.R script looks for a "local"
#   ./testScripts/R/ directory first before turning to ditto in the
#   installed packages.
# o BUG FIX: Using R.utils::cmdArgs(), which also drops the R executable,
#   which R.utils::commandArgs() used before did not.
# 2012-11-30
# o Now outputting session information useful for debugging.
# 2012-11-21
# o Now using R.cache root path specific to the aroma tests.
# 2012-10-19
# o ROBUSTNESS: Better pathname validation.
# o Adding system details to output at the beginning.
# 2012-09-14
# o Created.
############################################################################
