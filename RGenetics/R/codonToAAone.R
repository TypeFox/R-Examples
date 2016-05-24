#' @title Translating Codons to Amino Acids
#' @description Translate a DNA 3-character codon to an amino acid (one character)
#' @param codon a DNA 3-character codon
#' @seealso \code{\link{codonToAAthree}}
#' @export
#' @examples
#' codonToAAone('AAA')
codonToAAone <- function(codon){
    if(grepl('TCA', codon, ignore.case = TRUE)){return('S')}         # Serine
    else if(grepl('TCC', codon, ignore.case = TRUE)){return('S')}    # Serine
    else if(grepl('TCG', codon, ignore.case = TRUE)){return('S')}    # Serine
    else if(grepl('TCT', codon, ignore.case = TRUE)){return('S')}    # Serine
    else if(grepl('TTC', codon, ignore.case = TRUE)){return('F')}    # Phenylalanine
    else if(grepl('TTT', codon, ignore.case = TRUE)){return('F')}    # Phenylalanine
    else if(grepl('TTA', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('TTG', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('TAC', codon, ignore.case = TRUE)){return('Y')}    # Tyrosine
    else if(grepl('TAT', codon, ignore.case = TRUE)){return('Y')}    # Tyrosine
    else if(grepl('TAA', codon, ignore.case = TRUE)){return('Stop')} # Stop
    else if(grepl('TAG', codon, ignore.case = TRUE)){return('Stop')} # Stop
    else if(grepl('TGC', codon, ignore.case = TRUE)){return('C')}    # Cysteine
    else if(grepl('TGT', codon, ignore.case = TRUE)){return('C')}    # Cysteine
    else if(grepl('TGA', codon, ignore.case = TRUE)){return('Stop')} # Stop
    else if(grepl('TGG', codon, ignore.case = TRUE)){return('W')}    # Tryptophan
    else if(grepl('CTA', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('CTC', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('CTG', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('CTT', codon, ignore.case = TRUE)){return('L')}    # Leucine
    else if(grepl('CCA', codon, ignore.case = TRUE)){return('P')}    # Proline
    else if(grepl('CCC', codon, ignore.case = TRUE)){return('P')}    # Proline
    else if(grepl('CCG', codon, ignore.case = TRUE)){return('P')}    # Proline
    else if(grepl('CCT', codon, ignore.case = TRUE)){return('P')}    # Proline
    else if(grepl('CAC', codon, ignore.case = TRUE)){return('H')}    # Histidine
    else if(grepl('CAT', codon, ignore.case = TRUE)){return('H')}    # Histidine
    else if(grepl('CAA', codon, ignore.case = TRUE)){return('Q')}    # Glutamine
    else if(grepl('CAG', codon, ignore.case = TRUE)){return('Q')}    # Glutamine
    else if(grepl('CGA', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('CGC', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('CGG', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('CGT', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('ATA', codon, ignore.case = TRUE)){return('I')}    # Isoleucine
    else if(grepl('ATC', codon, ignore.case = TRUE)){return('I')}    # Isoleucine
    else if(grepl('ATT', codon, ignore.case = TRUE)){return('I')}    # Isoleucine
    else if(grepl('ATG', codon, ignore.case = TRUE)){return('M')}    # Methionine
    else if(grepl('ACA', codon, ignore.case = TRUE)){return('T')}    # Threonine
    else if(grepl('ACC', codon, ignore.case = TRUE)){return('T')}    # Threonine
    else if(grepl('ACG', codon, ignore.case = TRUE)){return('T')}    # Threonine
    else if(grepl('ACT', codon, ignore.case = TRUE)){return('T')}    # Threonine
    else if(grepl('AAC', codon, ignore.case = TRUE)){return('N')}    # Asparagine
    else if(grepl('AAT', codon, ignore.case = TRUE)){return('N')}    # Asparagine
    else if(grepl('AAA', codon, ignore.case = TRUE)){return('K')}    # Lysine
    else if(grepl('AAG', codon, ignore.case = TRUE)){return('K')}    # Lysine
    else if(grepl('AGC', codon, ignore.case = TRUE)){return('S')}    # Serine
    else if(grepl('AGT', codon, ignore.case = TRUE)){return('S')}    # Serine
    else if(grepl('AGA', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('AGG', codon, ignore.case = TRUE)){return('R')}    # Arginine
    else if(grepl('GTA', codon, ignore.case = TRUE)){return('V')}    # Valine
    else if(grepl('GTC', codon, ignore.case = TRUE)){return('V')}    # Valine
    else if(grepl('GTG', codon, ignore.case = TRUE)){return('V')}    # Valine
    else if(grepl('GTT', codon, ignore.case = TRUE)){return('V')}    # Valine
    else if(grepl('GCA', codon, ignore.case = TRUE)){return('A')}    # Alanine
    else if(grepl('GCC', codon, ignore.case = TRUE)){return('A')}    # Alanine
    else if(grepl('GCG', codon, ignore.case = TRUE)){return('A')}    # Alanine
    else if(grepl('GCT', codon, ignore.case = TRUE)){return('A')}    # Alanine
    else if(grepl('GAC', codon, ignore.case = TRUE)){return('D')}    # Aspartic Acid
    else if(grepl('GAT', codon, ignore.case = TRUE)){return('D')}    # Aspartic Acid
    else if(grepl('GAA', codon, ignore.case = TRUE)){return('E')}    # Glutamic Acid
    else if(grepl('GAG', codon, ignore.case = TRUE)){return('E')}    # Glutamic Acid
    else if(grepl('GGA', codon, ignore.case = TRUE)){return('G')}    # Glycine
    else if(grepl('GGC', codon, ignore.case = TRUE)){return('G')}    # Glycine
    else if(grepl('GGG', codon, ignore.case = TRUE)){return('G')}    # Glycine
    else if(grepl('GGT', codon, ignore.case = TRUE)){return('G')}    # Glycine
    else {stop(paste("Bad code: ", codon, ". code should only contains three of 'A','T','G','C'",sep=""))}
}
