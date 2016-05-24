sd_section("Overview",
           "",
           c("alakazam"))

sd_section("File I/O",
           "",
           c("readChangeoDb", "writeChangeoDb"))

sd_section("Sequence cleaning",
           "",
           c("maskSeqEnds", "maskSeqGaps", "collapseDuplicates"))

sd_section("Lineage reconstruction",
           "",
           c("makeChangeoClone", "buildPhylipLineage", "ChangeoClone-class"))

sd_section("Diversity analysis",
           "",
           c("countClones", "estimateAbundance", "rarefyDiversity", "testDiversity",
             "calcCoverage", "calcDiversity", "plotAbundance", "plotDiversityCurve", 
             "DiversityCurve-class", "DiversityTest-class"))

sd_section("Ig and TCR sequence annotation",
           "",
           c("countGenes", "extractVRegion", "getSegment"))

sd_section("Sequence distance calculation",
           "",
           c("getSeqDistance", "getSeqMatrix", "testSeqEqual", 
             "getAAMatrix", "getDNAMatrix"))

sd_section("Amino acid propertes",
           "",
           c("translateDNA", "aminoAcidProperties", "countPatterns", "isValidAASeq",
             "aliphatic", "bulk", "charge", "gravy", "polar"))

sd_section("Data and constants",
           "",
           c("IUPAC_CODES", "ABBREV_AA", "DEFAULT_COLORS", "IMGT_REGIONS", 
             "ExampleTrees"))