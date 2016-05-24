if (interactive()) savehistory();
library("aroma.affymetrix");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Allocate UGC file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chipType <- "Mapping250K_Nsp";
cdf <- AffymetrixCdfFile$byChipType(chipType);
ugc <- AromaUnitGcContentFile$allocateFromCdf(cdf, tags="na28,h=500kb,HB20090602");
print(ugc);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import GC contents from NetAffx files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
csv <- AffymetrixNetAffxCsvFile$byChipType(chipType, tags=".na28");

colClasses <- c("^(probeSetID|%GC)$"="character");
data <- readDataFrame(csv, colClasses=colClasses);
units <- indexOf(cdf, names=data$probeSetID);
stopifnot(all(is.finite(units)));
values <- as.double(data[["%GC"]]);
ugc[units,1] <- values;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
srcFileTags <- list();
srcFiles <- c(list(cdf), list(csv));
for (kk in seq(along=srcFiles)) {
  srcFile <- srcFiles[[kk]];
  tags <- list(
    filename=getFilename(srcFile), 
    filesize=getFileSize(srcFile), 
    checksum=getChecksum(srcFile)
  );
  srcFileTags[[kk]] <- tags;
}
print(srcFileTags);

footer <- readFooter(ugc);
footer$createdOn <- format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE);
footer$createdBy = list(
  fullname = "Henrik Bengtsson", 
  email = sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org")
);
names(srcFileTags) <- sprintf("srcFile%d", seq(along=srcFileTags));
footer$srcFiles <- srcFileTags;
footer$gcBinWidth <- as.integer(500e3);
writeFooter(ugc, footer);


print(ugc);
print(summary(ugc));
print(range(ugc[,1], na.rm=TRUE));


