# $Id: biostrings.R 152 2014-06-08 07:38:00Z larssn $


translate <- function( nuc.strings, Mstart=TRUE ){
# Kristian Hovde Liland and Lars Snipen
# Biostatistics, Norwegian University of Life Sciences
  
  codons <- c("GCA","GCC","GCG","GCT","GCX",
              "CGA","CGC","CGG","CGT","AGA","AGG","CGX",
              "AAC","AAT",
              "GAC","GAT",
              "TGC","TGT",
              "CAA","CAG",
              "GAA","GAG",
              "GGA","GGC","GGG","GGT","GGX",
              "CAC","CAT",
              "ATA","ATC","ATT",
              "ATG",
              "TTA","TTG","CTA","CTC","CTG","CTT","CTX",
              "AAA","AAG",
              "TTT","TTC",
              "CCA","CCC","CCG","CCT","CCX",
              "TCA","TCC","TCG","TCT","AGT","AGC","TCX",
              "ACA","ACC","ACG","ACT","ACX",
              "TGG",
              "TAT","TAC",
              "GTA","GTC","GTG","GTT","GTX",
              "TAA","TGA","TAG",
              "XAA","XAC","XAG","XAT","XCA","XCC","XCG","XCT","XGA","XGC","XGG","XGT","XTA","XTC","XTG","XTT",
              "AXA","AXC","AXG","AXT","CXA","CXC","CXG","CXT","GXA","GXC","GXG","GXT","TXA","TXC","TXG","TXT",
              "AAX","AGX","ATX","CAX","GAX","TAX","TGX","TTX",
              "AXX","CXX","GXX","TXX",
              "XXA","XXC","XXG","XXT",
              "XAX","XCX","XGX","XTX",
              "XXX")
  codons.translated <- c("A","A","A","A","A",
                         "R","R","R","R","R","R","R",
                         "N","N",
                         "D","D",
                         "C","C",
                         "Q","Q",
                         "E","E",
                         "G","G","G","G","G",
                         "H","H",
                         "I","I","I",
                         "M",
                         "L","L","L","L","L","L","L",
                         "K","K",
                         "F","F",
                         "P","P","P","P","P",
                         "S","S","S","S","S","S","S",
                         "T","T","T","T","T",
                         "W",
                         "Y","Y",
                         "V","V","V","V","V",
                         "*","*","*",
                         "X","X","X","X","X","X","X","X","X","X","X","X","X","X","X","X",
                         "X","X","X","X","X","X","X","X","X","X","X","X","X","X","X","X",
                         "X","X","X","X","X","X","X","X",
                         "X","X","X","X",
                         "X","X","X","X",
                         "X","X","X","X",
                         "X")
  nuc.strings <- toupper( nuc.strings )
  nuc.strings <- gsub( "U", "T", nuc.strings )
  nuc.strings <- gsub( "[^ACGTX]", "X", nuc.strings )
  ns <- length( nuc.strings )
  aa.strings <- character( ns )
                     
  n <- nchar( nuc.strings )
  is.ok <- (n %% 3) == 0
  if( sum( is.ok ) < ns ){
    idx.not <- which( !is.ok )
    warning( "Sequence ", paste( idx.not, collapse=" " ), " not divisible by 3, translation may be out-of-frame!\n" )
    nn <- floor( n[idx.not]/3 ) * 3
    nuc.strings[idx.not] <- substring( nuc.strings[idx.not], rep( 1, length(idx.not) ), nn )
  }
  n <- nchar( nuc.strings )
  
  aa.strings <- vapply( 1:ns, function( i ){
    s1 <- seq( from=1, to=n[i], by=3 )
    s2 <- seq( from=3, to=n[i], by=3 )
    paste( codons.translated[match( substring( nuc.strings[i], s1, s2 ), codons )], sep="", collapse="" )
  }, "", USE.NAMES=FALSE)
  if( Mstart ) aa.strings <- gsub( "^[VL]", "M", aa.strings )
  
  return( aa.strings )
}



reverseComplement <- function( nuc.strings ){
# Kristian Hovde Liland and Lars Snipen
# Biostatistics, Norwegian University of Life Sciences
  
  nuc.strings <- toupper( nuc.strings )
  nuc.strings <- gsub( "U", "T", nuc.strings )
  nuc.strings <- gsub( "[^ACGTX]", "X", nuc.strings )
  ACGT <- c("A","C","G","T","X")
  TGCA <- c("T","G","C","A","X")
  nuc.strings <- toupper( nuc.strings )
  nuc.strings <- strsplit( nuc.strings, "" )
  nuc.strings <- lapply( nuc.strings, match, ACGT )
  nuc.strings <- lapply( nuc.strings, rev )

  rc.strings <- vapply( 1:length( nuc.strings ), function( i ){
      cc <- TGCA[nuc.strings[[i]]]
      return( paste( cc, collapse="" ) )
    }, "", USE.NAMES=FALSE )
  return( rc.strings )
}
