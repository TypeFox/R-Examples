#' @title Translating Codons to Amino Acids
#' @description Translate a RNA 3-character codon to an amino acid (one character)
#' @param codon a RNA 3-character codon
#' @seealso \code{\link{codonToAAthreeRNA}}
#' @export
#' @examples
#' codonToAAoneRNA('AAA')
codonToAAoneRNA <- function(codon){
    if(grepl('UCA', codon, ignore.case = TRUE)){return('S')}         # Serine
    else if(grepl('UCC', codon, ignore.case = TRUE)){return('S')}    # Serine
    else if(grepl('UCG', codon, ignore.case = TRUE)){return('S')}    # Serine
    else if(grepl('UCU', codon, ignore.case = TRUE)){return('S')}    # Serine
    else if(grepl('UUC', codon, ignore.case = TRUE)){return('F')}    # Phenylalanine
    else if(grepl('UUU', codon, ignore.case = TRUE)){return('F')}    # Phenylalanine
    else if(grepl('UUA', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('UUG', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('UAC', codon, ignore.case = TRUE)){return('Y')}    # Tyrosine
    else if(grepl('UAU', codon, ignore.case = TRUE)){return('Y')}    # Tyrosine
    else if(grepl('UAA', codon, ignore.case = TRUE)){return('Stop')} # Stop
    else if(grepl('UAG', codon, ignore.case = TRUE)){return('Stop')} # Stop
    else if(grepl('UGC', codon, ignore.case = TRUE)){return('C')}    # Cysteine
    else if(grepl('UGU', codon, ignore.case = TRUE)){return('C')}    # Cysteine
    else if(grepl('UGA', codon, ignore.case = TRUE)){return('Stop')} # Stop
    else if(grepl('UGG', codon, ignore.case = TRUE)){return('W')}    # Tryptophan
    else if(grepl('CUA', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('CUC', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('CUG', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('CUU', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('CCA', codon, ignore.case = TRUE)){return('P')}    # Proline
    else if(grepl('CCC', codon, ignore.case = TRUE)){return('P')}    # Proline
    else if(grepl('CCG', codon, ignore.case = TRUE)){return('P')}    # Proline
    else if(grepl('CCU', codon, ignore.case = TRUE)){return('P')}    # Proline
    else if(grepl('CAC', codon, ignore.case = TRUE)){return('H')}    # Histidine
    else if(grepl('CAU', codon, ignore.case = TRUE)){return('H')}    # Histidine
    else if(grepl('CAA', codon, ignore.case = TRUE)){return('Q')}    # Glutamine
    else if(grepl('CAG', codon, ignore.case = TRUE)){return('Q')}    # Glutamine
    else if(grepl('CGA', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('CGC', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('CGG', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('CGU', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('AUA', codon, ignore.case = TRUE)){return('I')}    # Isoleucine
    else if(grepl('AUC', codon, ignore.case = TRUE)){return('I')}    # Isoleucine
    else if(grepl('AUU', codon, ignore.case = TRUE)){return('I')}    # Isoleucine
    else if(grepl('AUG', codon, ignore.case = TRUE)){return('M')}    # Methionine
    else if(grepl('ACA', codon, ignore.case = TRUE)){return('T')}    # Threonine
    else if(grepl('ACC', codon, ignore.case = TRUE)){return('T')}    # Threonine
    else if(grepl('ACG', codon, ignore.case = TRUE)){return('T')}    # Threonine
    else if(grepl('ACU', codon, ignore.case = TRUE)){return('T')}    # Threonine
    else if(grepl('AAC', codon, ignore.case = TRUE)){return('N')}    # Asparagine
    else if(grepl('AAU', codon, ignore.case = TRUE)){return('N')}    # Asparagine
    else if(grepl('AAA', codon, ignore.case = TRUE)){return('K')}    # Lysine
    else if(grepl('AAG', codon, ignore.case = TRUE)){return('K')}    # Lysine
    else if(grepl('AGC', codon, ignore.case = TRUE)){return('S')}    # Serine
    else if(grepl('AGU', codon, ignore.case = TRUE)){return('S')}    # Serine
    else if(grepl('AGA', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('AGG', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('GUA', codon, ignore.case = TRUE)){return('V')}    # Valine
    else if(grepl('GUC', codon, ignore.case = TRUE)){return('V')}    # Valine
    else if(grepl('GUG', codon, ignore.case = TRUE)){return('V')}    # Valine
    else if(grepl('GUU', codon, ignore.case = TRUE)){return('V')}    # Valine
    else if(grepl('GCA', codon, ignore.case = TRUE)){return('A')}    # Alanine
    else if(grepl('GCC', codon, ignore.case = TRUE)){return('A')}    # Alanine
    else if(grepl('GCG', codon, ignore.case = TRUE)){return('A')}    # Alanine
    else if(grepl('GCU', codon, ignore.case = TRUE)){return('A')}    # Alanine
    else if(grepl('GAC', codon, ignore.case = TRUE)){return('D')}    # Aspartic Acid
    else if(grepl('GAU', codon, ignore.case = TRUE)){return('D')}    # Aspartic Acid
    else if(grepl('GAA', codon, ignore.case = TRUE)){return('E')}    # Glutamic Acid
    else if(grepl('GAG', codon, ignore.case = TRUE)){return('E')}    # Glutamic Acid
    else if(grepl('GGA', codon, ignore.case = TRUE)){return('G')}    # Glycine
    else if(grepl('GGC', codon, ignore.case = TRUE)){return('G')}    # Glycine
    else if(grepl('GGG', codon, ignore.case = TRUE)){return('G')}    # Glycine
    else if(grepl('GGU', codon, ignore.case = TRUE)){return('G')}    # Glycine
    else {stop(paste("Bad code: ", codon, ". code should only contains three of 'A','U','G','C'",sep=""))}
}
