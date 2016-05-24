.Workspace <- new.env()
.Workspace$deserialisers <- list()
.Workspace$pathHandlers <- list()

.Analyze <- list(
    datatypes=list(codes=c(     2,          4,          8,          16,         64),
                   rTypes=c(   "integer",  "integer",  "integer",  "double",   "double"),
                   sizes=c(     1,          2,          4,          4,          8),
                   isSigned=c(  FALSE,      TRUE,       TRUE,       TRUE,       TRUE)))

.Dicom <- list(
    nonCharTypes=list(codes=c("OF", "FL", "FD", "SL", "SS", "UL", "US", "AT"),
                      rTypes=c(rep("double",3), rep("integer",5)),
                      sizes=c(4, 4, 8, 4, 2, 4, 2, 2),
                      counts=c(1, 1, 1, 1, 1, 1, 1, 2),
                      isSigned=c(rep(TRUE,5), rep(FALSE,3))),
    longTypes=c("OB", "OW", "OF", "SQ", "UT", "UN"),
    convertibleTypes=c("OF", "FL", "FD", "SL", "SS", "UL", "US", "AT", "DS", "IS"),
    transferSyntaxes=list("1.2.840.10008.1.2"   = list(endian="little",explicitTypes=FALSE),
                          "1.2.840.10008.1.2.1" = list(endian="little",explicitTypes=TRUE),
                          "1.2.840.10008.1.2.2" = list(endian="big",explicitTypes=TRUE)))

.Nifti <- list(
    datatypes=list(codes=c(     2,          4,          8,          16,         64,         256,        512,        768),
                   rTypes=c(   "integer",  "integer",  "integer",  "double",   "double",   "integer",  "integer",  "integer"),
                   sizes=c(     1,          2,          4,          4,          8,          1,          2,          4),
                   isSigned=c(  FALSE,      TRUE,       TRUE,       TRUE,       TRUE,       TRUE,       FALSE,      FALSE)),
    units=list(unknown=0, m=1, mm=2, um=3, s=8, ms=16, us=24),
    xformCodes=list(unknown=0, scannerAnatomical=1, alignedAnatomical=2, talairach=3, mni=4),
    magicStrings=list(list(c(charToRaw("ni1"),as.raw(0)), c(charToRaw("n+1"),as.raw(0))),
                      list(c(charToRaw("ni2"),as.raw(0)), c(charToRaw("n+2"),as.raw(0)))))

.Mgh <- list(
    datatypes=list(codes=c(     0,          4,          1,          3      ),
                   rTypes=c(   "integer",  "integer",  "integer",  "double"),
                   sizes=c(     1,          2,          4,          4      ),
                   isSigned=c(  FALSE,      TRUE,       TRUE,       TRUE   )))

.FileTypes <- list(
    typeNames=c(       "ANALYZE",  "NIFTI",    "NIFTI_PAIR",   "ANALYZE_GZ",   "NIFTI_GZ", "NIFTI_PAIR_GZ", "MGH",  "MGH_GZ"),
    formatNames=c(     "Analyze",  "Nifti",    "Nifti",        "Analyze",      "Nifti",    "Nifti",         "Mgh",  "Mgh"),
    singleFile=c(       NA,         TRUE,       FALSE,          NA,             TRUE,       FALSE,           NA,     NA),
    gzipped=c(          FALSE,      FALSE,      FALSE,          TRUE,           TRUE,       TRUE,            FALSE,  TRUE),
    headerSuffixes=c(  "hdr",      "nii",      "hdr",          "hdr.gz",       "nii.gz",   "hdr.gz",        "mgh",  "mgz"),
    imageSuffixes=c(   "img",      "nii",      "img",          "img.gz",       "nii.gz",   "img.gz",        "mgh",  "mgz"))
