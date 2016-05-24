if (interactive()) savehistory();
library("aroma.affymetrix");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Allocate UGC file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full");
ugc <- AromaUnitGcContentFile$allocateFromCdf(cdf, tags="na27,h=500kb,HB20090322");
print(ugc);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import GC contents from NetAffx files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
csvList <- list(
  AffymetrixNetAffxCsvFile$byChipType("GenomeWideSNP_6", tags=".cn.na27"),
  AffymetrixNetAffxCsvFile$byChipType("GenomeWideSNP_6", tags=".na27.1")
);

colClasses <- c("^(probeSetID|%GC)$"="character");
for (csv in csvList) {
  data <- readDataFrame(csv, colClasses=colClasses);
  units <- indexOf(cdf, names=data$probeSetID);
  ugc[units,1] <- as.double(data[["%GC"]]);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
srcFileTags <- list();
srcFiles <- c(list(cdf), csvList);
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
print(range(ugc[,1]));

