## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
library(alakazam)

# Load Change-O file
file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
db <- readChangeoDb(file)
db <- db[db$SAMPLE == "RL01", ]

## ---- eval=TRUE, warning=FALSE, fig.width=7.5, fig.height=6--------------
db_props <- aminoAcidProperties(db, seq="JUNCTION", nt=TRUE, trim=TRUE, 
                                label="CDR3")

# The full set of properties are calculated by default
dplyr::select(db_props[1:3, ], starts_with("CDR3"))

# Plot a subset of the properties
tmp_theme <- theme_bw() + theme(legend.position="bottom")
g1 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_LENGTH)) + tmp_theme +
    ggtitle("CDR3 length") + 
    xlab("Isotype") + ylab("Amino acids") +
    scale_fill_manual(name="Isotype", values=IG_COLORS) +
    geom_boxplot(aes(fill=ISOTYPE))
g2 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_GRAVY)) + tmp_theme + 
    ggtitle("CDR3 hydrophobicity") + 
    xlab("Isotype") + ylab("GRAVY") +
    scale_fill_manual(name="Isotype", values=IG_COLORS) +
    geom_boxplot(aes(fill=ISOTYPE))
g3 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_BASIC)) + tmp_theme +
    ggtitle("CDR3 basic residues") + 
    xlab("Isotype") + ylab("Basic residues") +
    scale_y_continuous(labels=scales::percent) +
    scale_fill_manual(name="Isotype", values=IG_COLORS) +
    geom_boxplot(aes(fill=ISOTYPE))
g4 <- ggplot(db_props, aes(x=ISOTYPE, y=CDR3_AA_ACIDIC)) + tmp_theme +
    ggtitle("CDR3 acidic residues") + 
    xlab("Isotype") + ylab("Acidic residues") +
    scale_y_continuous(labels=scales::percent) +
    scale_fill_manual(name="Isotype", values=IG_COLORS) +
    geom_boxplot(aes(fill=ISOTYPE))
multiggplot(g1, g2, g3, g4, ncol=2)

## ---- eval=TRUE, warning=FALSE-------------------------------------------
db_props <- aminoAcidProperties(db, seq="JUNCTION", property=c("gravy", "charge"),
                                nt=TRUE, trim=TRUE, label="CDR3")
dplyr::select(db_props[1:3, ], starts_with("CDR3"))

## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
# Load the relevant data objects from the seqinr package
library(seqinr)
data(aaindex)
data(pK)
h <- aaindex[["KIDA850101"]]$I
p <- setNames(pK[["Murray"]], rownames(pK))
# Rename the hydrophobicity vector to use single-letter codes
names(h) <- translateStrings(names(h), ABBREV_AA)
db_props <- aminoAcidProperties(db, seq="JUNCTION", property=c("gravy", "charge"), 
                                nt=TRUE, trim=TRUE, label="CDR3", 
                                hydropathy=h, pK=p)
dplyr::select(db_props[1:3, ], starts_with("CDR3"))

## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
# Translate junction DNA sequences to amino acids and trim first and last codons
cdr3 <- translateDNA(db$JUNCTION[1:3], trim=TRUE)

# Grand average of hydrophobicity
gravy(cdr3)

# Average bulkiness
bulk(cdr3)

# Average polarity
polar(cdr3)

# Normalized aliphatic index
aliphatic(cdr3)
# Unnormalized aliphatic index
aliphatic(cdr3, normalize=FALSE)

# Normalized net charge
charge(cdr3)
# Unnormalized net charge
charge(cdr3, normalize=FALSE)

# Count of acidic amino acids
# Takes a named list of regular expressions
countPatterns(cdr3, c(ACIDIC="[DE]"), label="CDR3")

