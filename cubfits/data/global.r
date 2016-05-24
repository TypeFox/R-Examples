# global variables.

.CF.GV <- list(
  amino.acid = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                 "M", "N", "P", "Q", "R", "S", "T", "V", "W", "X", "Y"),
  amino.acid.3 = c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His",
                   "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln",
                   "Arg", "Ser", "Thr", "Val", "Trp", "Stop", "Tyr"),
  synonymous.codon = list(
    A = c("GCA", "GCC", "GCG", "GCT"),
    C = c("TGC", "TGT"),
    D = c("GAC", "GAT"),
    E = c("GAA", "GAG"),
    F = c("TTC", "TTT"),
    G = c("GGA", "GGC", "GGG", "GGT"),
    H = c("CAC", "CAT"),
    I = c("ATA", "ATC", "ATT"),
    K = c("AAA", "AAG"),
    L = c("CTA", "CTC", "CTG", "CTT", "TTA", "TTG"),
    M = c("ATG"),
    N = c("AAC", "AAT"),
    P = c("CCA", "CCC", "CCG", "CCT"),
    Q = c("CAA", "CAG"),
    R = c("AGA", "AGG", "CGA", "CGC", "CGG", "CGT"),
    S = c("AGC", "AGT", "TCA", "TCC", "TCG", "TCT"),
    T = c("ACA", "ACC", "ACG", "ACT"),
    V = c("GTA", "GTC", "GTG", "GTT"),
    W = c("TGG"),
    X = c("TAA", "TAG", "TGA"),
    Y = c("TAC", "TAT")
  ),

  amino.acid.split = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "X", "Y",
                       "Z"),
  amino.acid.split.3 = c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His",
                         "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln",
                         "Arg", "Ser4", "Thr", "Val", "Trp", "Stop", "Tyr",
                         "Ser2"),
  synonymous.codon.split = list(
    A = c("GCA", "GCC", "GCG", "GCT"),
    C = c("TGC", "TGT"),
    D = c("GAC", "GAT"),
    E = c("GAA", "GAG"),
    F = c("TTC", "TTT"),
    G = c("GGA", "GGC", "GGG", "GGT"),
    H = c("CAC", "CAT"),
    I = c("ATA", "ATC", "ATT"),
    K = c("AAA", "AAG"),
    L = c("CTA", "CTC", "CTG", "CTT", "TTA", "TTG"),
    M = c("ATG"),
    N = c("AAC", "AAT"),
    P = c("CCA", "CCC", "CCG", "CCT"),
    Q = c("CAA", "CAG"),
    R = c("AGA", "AGG", "CGA", "CGC", "CGG", "CGT"),
    S = c("TCA", "TCC", "TCG", "TCT"),  # split 2 codons to Z.
    T = c("ACA", "ACC", "ACG", "ACT"),
    V = c("GTA", "GTC", "GTG", "GTT"),
    W = c("TGG"),
    X = c("TAA", "TAG", "TGA"),  # stop codons.
    Y = c("TAC", "TAT"),
    Z = c("AGC", "AGT")  # obtain 2 codons from S.
  )
) # End of .CF.GV.
