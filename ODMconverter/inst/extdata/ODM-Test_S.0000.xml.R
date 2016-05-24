odmdata <- data.frame( 
I.1001=logical(),
I.1002=integer(),
I.1003=character(),
I.1004=integer(),
I.1005=character(),
I.1006=character(),
I.1007=numeric(),
I.1008=character(),
stringsAsFactors = F)

attr(odmdata, "StudyOID")         <- "S.0000"
attr(odmdata, "Sponsor")          <- "Testsponsor"
attr(odmdata, "Condition")        <- "Testcondition"
attr(odmdata, "StudyName")        <- "ODM Test Study"
attr(odmdata, "StudyDescription") <- "Test of ODM tools"
attr(odmdata, "Form")             <- "ODM-Test"
attr(odmdata, "FirstName")        <- "Test"
attr(odmdata, "LastName")         <- "Testname"
attr(odmdata, "Organization")     <- "Test organization"

attr(odmdata, "varlabels") <- c(
"Willingness to participate in clinicial trials",
"Age",
"Date of Birth",
"Gender",
"Diagnosis text",
"Diagnosis code",
"Creatinine",
"Time of lab value")
attr(odmdata, "varlabels_en") <- c(
"Willingness to participate in clinicial trials",
"Age",
"Date of Birth",
"Gender",
"Diagnosis text",
"Diagnosis code",
"Creatinine",
"Time of lab value")
attr(odmdata, "varlabels_de") <- c(
"Bereitschaft zur Teilnahme an klinischen Studien",
"Alter",
"Geburtsdatum",
"Geschlecht",
"Diagnosetext",
"Diagnosekode",
"Kreatinin",
"Zeitpunkt des Laborwerts")
attr(odmdata, "DataTypes") <- c(
"boolean",
"integer",
"date",
"integer",
"string",
"string",
"float",
"time")
attr(odmdata, "itemgroups") <- c("IG.1","IG.1","IG.1","IG.1","IG.1","IG.1","IG.1","IG.1")
attr(odmdata, "itemgroups_unique") <- c("IG.1")
attr(odmdata, "itemgroupnames_unique") <- c("General Info")

attr(odmdata, "itemgroupnames_unique_en") <- c("General Information")

attr(odmdata, "itemgroupnames_unique_de") <- c("Allgemeine Angaben")

attr(odmdata, "IGAlias") <- c("IG.1","UMLS CUI","C0332118","IG.1","SNOMEDCT_2012_01_31","106227002")

odmdata$I.1004 <- factor(levels=c("1","2"))
attr(odmdata$I.1004,"labels") <- c("male","female")
attr(odmdata$I.1004,"labels2") <- c("maennlich","weiblich")
attr(odmdata$I.1004,"Alias_Contexts") <- c("UMLS CUI","SNOMED CT 2010_0731")
attr(odmdata$I.1004,"UMLS CUI") <- c("C0024554","C0015780")
attr(odmdata$I.1004,"SNOMED CT 2010_0731") <- c("248153007","248152002")
attr(odmdata, "Alias_Items") <- matrix(ncol=3,byrow=T, data=c(
"I.1001","UMLS CUI","C1516879",
"I.1002","SNOMED CT 2010_0731","102518004",
"I.1003","SNOMED CT 2010_0731","152322001",
"I.1004","SNOMED CT 2010_0731","139865004",
"I.1005","SNOMED CT 2010_0731","439401001",
"I.1007","LOINC","38483-4" ))
