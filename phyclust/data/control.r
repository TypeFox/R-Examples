### This file contains all control functions for EM.

### These orders are matter for C.
.init.procedure <- c("exhaustEM", "emEM", "RndEM", "RndpEM")
.init.method <- c("randomMu", "NJ", "randomNJ", "PAM", "K-Medoids", "manualMu")
.substitution.model <- data.frame(
  model = c("JC69", "K80", "F81", "HKY85", "SNP_JC69", "SNP_F81",
            "E_F81", "E_HKY85", "E_SNP_F81"),
  code.type = c(rep("NUCLEOTIDE", 4), rep("SNP", 2),
                rep("NUCLEOTIDE", 2), rep("SNP", 1))
)
.edist.model <- c("D_JC69", "D_K80", "D_HAMMING", "D_HAMMING_WOGAP")
.identifier <- c("EE", "EV", "VE", "VV")
.code.type <- c("NUCLEOTIDE", "SNP", "CODON", "AMINO_ACID")
.em.method <- c("EM", "ECM", "AECM")
.boundary.method <- c("ADJUST", "IGNORE")
.label.method <- c("NONE", "SEMI", "GENERAL")
.se.model <- c("CONVOLUTION")
