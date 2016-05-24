# Amino acid sequence properties

#' @include Alakazam.R
NULL

#### Chemical property functions ####

#' Calculates the average bulkiness of amino acid sequences
#'
#' \code{bulk} calculates the average bulkiness score of amino acid sequences. 
#' Non-informative positions are excluded, where non-informative is defined as any 
#' character in \code{c("X", "-", ".", "*")}.
#' 
#' @param    seq         vector of strings containing amino acid sequences.
#' @param    bulkiness   named numerical vector defining bulkiness scores for 
#'                       each amino acid, where names are single-letter amino acid 
#'                       character codes. If \code{NULL}, then the Zimmerman et al, 1968
#'                       scale is used.
#' 
#' @return   A vector of bulkiness scores for the sequence(s).
#' 
#' @references
#' \enumerate{
#'   \item  Zimmerman JM, Eliezer N, Simha R. The characterization of amino acid sequences 
#'            in proteins by statistical methods. J Theor Biol 21, 170-201 (1968).
#' }
#' @seealso 
#' For additional size related indices see \code{\link[seqinr]{aaindex}}.
#'
#' @examples
#' # Default bulkiness scale
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
#' bulk(seq)
#'
#' # Use the Grantham, 1974 side chain volumn scores from the seqinr package
#' library(seqinr)
#' data(aaindex)
#' x <- aaindex[["GRAR740103"]]$I
#' # Rename the score vector to use single-letter codes
#' names(x) <- translateStrings(names(x), ABBREV_AA)
#' # Calculate average volume
#' bulk(seq, bulkiness=x)
#' 
#' @export
bulk <- function(seq, bulkiness=NULL) {
    # Get bulkiness scores
    if (is.null(bulkiness)) {
        bulkiness <- BULKINESS_ZIMJ68
    }
    # Remove non-informative positions
    seq <- gsub("[X\\.\\*-]", "", as.character(seq))
    # Create character vector from string
    aa <- strsplit(seq, "")
    # Calculate average bulkiness
    aa_bulk <- sapply(aa, function(x) sum(bulkiness[x]) / length(x))
    
    return(aa_bulk)
}


#' Calculates the average polarity of amino acid sequences
#'
#' \code{polar} calculates the average polarity score of amino acid sequences. 
#' Non-informative positions are excluded, where non-informative is defined as any 
#' character in \code{c("X", "-", ".", "*")}.
#' 
#' @param    seq         vector of strings containing amino acid sequences.
#' @param    polarity    named numerical vector defining polarity scores for 
#'                       each amino acid, where names are single-letter amino acid 
#'                       character codes. If \code{NULL}, then the Grantham, 1974
#'                       scale is used.
#' 
#' @return   A vector of bulkiness scores for the sequence(s).
#' 
#' @references
#' \enumerate{
#'   \item  Grantham R. Amino acid difference formula to help explain protein evolution. 
#'            Science 185, 862-864 (1974).
#' }
#' @seealso 
#' For additional size related indices see \code{\link[seqinr]{aaindex}}.
#'
#' @examples
#' # Default scale
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
#' polar(seq)
#'
#' # Use the Zimmerman et al, 1968 polarity scale from the seqinr package
#' library(seqinr)
#' data(aaindex)
#' x <- aaindex[["ZIMJ680103"]]$I
#' # Rename the score vector to use single-letter codes
#' names(x) <- translateStrings(names(x), ABBREV_AA)
#' # Calculate polarity
#' polar(seq, polarity=x)
#' 
#' @export
polar <- function(seq, polarity=NULL) {
    # Get bulkiness scores
    if (is.null(polarity)) {
        polarity <- POLARITY_GRAR74
    }
    # Remove non-informative positions
    seq <- gsub("[X\\.\\*-]", "", as.character(seq))
    # Create character vector from string
    aa <- strsplit(seq, "")
    # Calculate average polarity
    aa_polar <- sapply(aa, function(x) sum(polarity[x]) / length(x))
    
    return(aa_polar)
}


#' Calculates the hydrophobicity of amino acid sequences
#'
#' \code{gravy} calculates the Grand Average of Hydrophobicity (GRAVY) index 
#' of amino acid sequences using the method of Kyte & Doolittle. Non-informative
#' positions are excluded, where non-informative is defined as any character in 
#' \code{c("X", "-", ".", "*")}.
#' 
#' @param    seq         vector of strings containing amino acid sequences.
#' @param    hydropathy  named numerical vector defining hydropathy index values for 
#'                       each amino acid, where names are single-letter amino acid 
#'                       character codes. If \code{NULL}, then the Kyte & Doolittle
#'                       scale is used.
#' 
#' @return   A vector of GRAVY scores for the sequence(s).
#' 
#' @references
#' \enumerate{
#'   \item  Kyte J, Doolittle RF. A simple method for displaying the hydropathic character 
#'            of a protein. J Mol Biol. 157, 105-32 (1982).
#' }
#' @seealso 
#' For additional hydrophobicity indices see \code{\link[seqinr]{aaindex}}.
#'
#' @examples
#' # Default scale
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT")
#' gravy(seq)
#'
#' # Use the Kidera et al, 1985 scores from the seqinr package
#' library(seqinr)
#' data(aaindex)
#' x <- aaindex[["KIDA850101"]]$I
#' # Rename the score vector to use single-letter codes
#' names(x) <- translateStrings(names(x), ABBREV_AA)
#' # Calculate hydrophobicity
#' gravy(seq, hydropathy=x)
#' 
#' @export
gravy <- function(seq, hydropathy=NULL) {
    # Get hydrophobicity scores
    if (is.null(hydropathy)) {
        hydropathy <- HYDROPATHY_KYTJ82
    }

    # Remove non-informative positions
    seq <- gsub("[X\\.\\*-]", "", as.character(seq))
    # Create character vector from string
    aa <- strsplit(seq, "")
    # Calculate GRAVY
    aa_gravy <- sapply(aa, function(x) sum(hydropathy[x]) / length(x))
    
    return(aa_gravy)
}


#' Calculates the aliphatic index of amino acid sequences
#' 
#' \code{aliphatic} calculates the aliphatic index of amino acid sequences using 
#' the method of Ikai. Non-informative positions are excluded, where non-informative 
#' is defined as any character in \code{c("X", "-", ".", "*")}.
#'
#' @param    seq        vector of strings containing amino acid sequences.
#' @param    normalize  if \code{TRUE} then divide the aliphatic index of each amino acid 
#'                      sequence by the number of informative positions. Non-informative 
#'                      position are defined by the presence any character in 
#'                      \code{c("X", "-", ".", "*")}. If \code{FALSE} then return the raw
#'                      aliphatic index.
#'                      
#' @return   A vector of the aliphatic indices for the sequence(s).
#'
#' @references 
#' \enumerate{
#'   \item  Ikai AJ. Thermostability and aliphatic index of globular proteins. 
#'            J Biochem. 88, 1895-1898 (1980).
#' }
#' 
#' @examples 
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", NA, "XXTQMYVRT")
#' aliphatic(seq)
#'
#' @export
aliphatic <- function(seq, normalize=TRUE) {
    # Calculate aliphatic index for valid amino acids
    ala <- countOccurrences(seq, "[A]")
    val <- countOccurrences(seq, "[V]")
    leu_ile <- countOccurrences(seq, "[LI]")
    
    aa_aliphatic = ala + 2.9 * val + 3.9 * leu_ile
    
    if (normalize) {
        aa_aliphatic <- aa_aliphatic / stri_length(gsub("[X\\.\\*-]", "", seq))
    }
    
    return(aa_aliphatic)
}


#' Calculates the net charge of amino acid sequences.
#'
#' \code{charge} calculates the net charge of amino acid sequences using 
#' the method of Moore, 1985, with exclusion of the C-terminus and N-terminus charges.
#' 
#' @param    seq        vector strings defining of amino acid sequences.
#' @param    pH         environmental pH.
#' @param    pK         named vector defining pK values for each charged amino acid,
#'                      where names are the single-letter amino acid character codes
#'                      \code{c("R", "H", "K", "D", "E", "C", "Y")}). If \code{NULL}, 
#'                      then the EMBOSS scale is used.
#' @param    normalize  if \code{TRUE} then divide the net charge of each amino acid 
#'                      sequence by the number of informative positions. Non-informative 
#'                      position are defined by the presence any character in 
#'                      \code{c("X", "-", ".", "*")}. If \code{FALSE} then return the raw
#'                      net charge.
#' 
#' @return   A vector of net charges for the sequence(s).
#'
#' @references
#' \enumerate{
#'   \item  Moore DS. Amino acid and peptide net charges: A simple calculational procedure. 
#'            Biochem Educ. 13, 10-11 (1985).
#'   \item  \url{http://emboss.sourceforge.net/apps/cvs/emboss/apps/iep.html}
#'   }
#'   
#' @seealso 
#' For additional pK scales see \code{\link[seqinr]{pK}}.
#' 
#' @examples 
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVRT") 
#' # Normalized charge
#' charge(seq)
#' # Unnormalized charge
#' charge(seq, normalize=FALSE)
#' 
#' # Use the Murray et al, 2006 scores from the seqinr package
#' library(seqinr)
#' data(pK)
#' x <- setNames(pK[["Murray"]], rownames(pK))
#' # Calculate charge
#' charge(seq, pK=x, normalize=FALSE)
#'
#' @export
charge <- function(seq, pH=7.4, pK=NULL, normalize=TRUE) {
    
    # Get charge data
    if(is.null(pK)) {
        pK <- PK_EMBOSS
    }
    
    # Calculate charge
    arg <- countOccurrences(seq, "R") * (1/(1 + 10^(1 * (pH - pK["R"]))))
    his <- countOccurrences(seq, "H") * (1/(1 + 10^(1 * (pH - pK["H"]))))
    lys <- countOccurrences(seq, "K") * (1/(1 + 10^(1 * (pH - pK["K"]))))
    asp <- countOccurrences(seq, "D") * (-1/(1 + 10^(-1 * (pH - pK["D"]))))
    glu <- countOccurrences(seq, "E") * (-1/(1 + 10^(-1 * (pH - pK["E"]))))
    cys <- countOccurrences(seq, "C") * (-1/(1 + 10^(-1 * (pH - pK["C"]))))
    tyr <- countOccurrences(seq, "Y") * (-1/(1 + 10^(-1 * (pH - pK["Y"]))))
    aa_charge <- arg + lys + his + asp + glu + tyr + cys
    
    if (normalize) {
        aa_charge <- aa_charge / stri_length(gsub("[X\\.\\*-]", "", seq))
    }
    
    return(aa_charge)
}

#' Validate amino acid sequences
#'
#' \code{isValidAASeq} checks that a set of sequences are valid non-ambiguous 
#' amino acid sequences. A sequence is considered valid if it contains only 
#' characters in the the non-ambiguous IUPAC character set or any characters in 
#' \code{c("X", ".", "-", "*")}.
#'  
#' @param    seq  character vector of sequences to check.
#' 
#' @return   A logical vector with \code{TRUE} for each valid amino acid sequences 
#'           and \code{FALSE} for each invalid sequence.
#' @seealso 
#' See \code{\link{ABBREV_AA}} for the set of non-ambiguous amino acid characters.
#' See \code{\link{IUPAC_AA}} for the full set of ambiguous amino acid characters.
#' 
#' @examples 
#' seq <- c("CARDRSTPWRRGIASTTVRTSW", "XXTQMYVR--XX", "CARJ", "10") 
#' isValidAASeq(seq)
#' 
#' @export
isValidAASeq <- function(seq) {
    
    # Get valid amino acids from seqinr
    # for consistency with `gravy` and other
    # amino acid properties that don't consider
    # amino acid ambiguities and special encoded amino acids
    # http://pir.georgetown.edu/resid/faq.shtml#q01
    # Also include here characters for non informative positions
    valid_AA <- c(names(ABBREV_AA), "X", ".", "*", "-")
    
    .isValid <- function(aa) {
        all(aa %in% valid_AA)
    }
    
    return(sapply(strsplit(seq, ""), .isValid))
    
#    valid_AA <- paste(c(names(ABBREV_AA),"X.*-"),collapse="")
#    valid <- !grepl(paste0("[^",valid_AA,"]"), seq) & !is.na(seq)
#    valid
    
}

# Count patterns
# 
# Counts the number of times a "pattern" occurs in "x", a string
#
# @param   x        a string (usually amino acids)
# @param   pattern  regular expression to be matched in string
#  
# @return  number of times the regular expression was found
countOccurrences <- function(x, pattern) {
    return(sapply(gregexpr(pattern, x), function(y) { sum(y > 0) }))
}


#' Count sequence patterns
#' 
#' \code{countPatterns} counts the fraction of times a set of character patterns occur 
#' in a set of sequences.
#'
#' @param   seq         character vector of either DNA or amino acid sequences.
#' @param   patterns    list of sequence patterns to count in each sequence. If the 
#'                      list is named, then names will be assigned as the column names of 
#'                      output data.frame.
#' @param   nt		    if \code{TRUE} then \code{seq} are DNA sequences and and will be 
#'                      translated before performing the pattern search.
#' @param   trim        if \code{TRUE} remove the first and last codon or amino acid from 
#'                      each sequence before the pattern search. If \code{FALSE} do
#'                      not modify the input sequences.
#' @param   label       string defining a label to add as a prefix to the output 
#'                      column names.
#'  
#' @return  A data.frame containing the fraction of times each sequence pattern was 
#'          found.
#' 
#' @examples 
#' seq <- c("TGTCAACAGGCTAACAGTTTCCGGACGTTC",
#'          "TGTCAGCAATATTATATTGCTCCCTTCACTTTC",
#'          "TGTCAAAAGTATAACAGTGCCCCCTGGACGTTC")
#' patterns <- c("A", "V", "[LI]")
#' names(patterns) <- c("ARG", "VAL", "ISO_LEU")
#' countPatterns(seq, patterns, nt=TRUE, trim=TRUE, label="CDR3")
#'             
#' @export
countPatterns <- function(seq, patterns, nt=FALSE, trim=FALSE, label="REGION") {
    # Translate sequences if nucleotide
    region_aa <- if (nt) { translateDNA(seq, trim=trim) } else { seq }
    
    # TODO: Check that NA is passed through correctly.
    # TODO: What is the proper length normalization? With or without non-informative position?
    
    # Calculate region lengths
    aa_length <- stri_length(region_aa)
    # Count occurrence of each amino acid pattern for each sequence
    out_df <- data.frame(matrix(0, nrow=length(region_aa), ncol=length(patterns)))
    # If patterns are unnamed, make the names X1...Xn
    if(is.null(names(patterns))) { names(patterns) <- names(out_df) }
    # If region name, append to names of patterns
    if(label != '') { 
        names(out_df) <- paste(label, names(patterns), sep="_") 
    } else {
        names(out_df) <- names(patterns)
    }
    # Iterate over patterns
    for(i in 1:length(patterns)) {
        out_df[, i] <- countOccurrences(region_aa, patterns[i]) / aa_length
    }
    return(out_df)
}


#' Calculates amino acid chemical properties for sequence data
#'
#' \code{aminoAcidProperties} calculates amino acid sequence physicochemical properties, including
#' length, hydrophobicity, bulkiness, polarity, aliphatic index, net charge, acidic residue
#' content, basic residue content, and aromatic residue content.
#'
#' @param   data          \code{data.frame} containing sequence data.
#' @param   property      vector strings specifying the properties to be calculated. Defaults
#'                        to calculating all defined properties.
#' @param   seq           \code{character} name of the column containing input 
#'                        sequences.
#' @param   nt      	  boolean, TRUE if the sequences (or sequence) are DNA and will be translated.
#' @param   trim          if \code{TRUE} remove the first and last codon/amino acids from each
#'                        sequence before calculating properties. If \code{FALSE} do
#'                        not modify input sequences.
#' @param   label         name of sequence region to add as prefix to output column names.
#' @param   ...           additional named arguments to pass to the functions 
#'                        \link{gravy}, \link{bulk}, \link{aliphatic}, \link{polar} or \link{charge}.
#' 
#' @return  A modified \code{data} data.frame with the following columns:
#'          \itemize{
#'            \item  \code{*_AA_LENGTH}:     number of amino acids.
#'            \item  \code{*_AA_GRAVY}:      grand average of hydrophobicity (GRAVY) index.
#'            \item  \code{*_AA_BULK}:       average bulkiness of amino acids.
#'            \item  \code{*_AA_ALIPHATIC}:  aliphatic index.
#'            \item  \code{*_AA_POLARITY}:   average polarity of amino acids.
#'            \item  \code{*_AA_CHARGE}:     normalized net charge.
#'            \item  \code{*_AA_BASIC}:      fraction of informative positions that are 
#'                                           Arg, His or Lys.
#'            \item  \code{*_AA_ACIDIC}:     fraction of informative positions that are 
#'                                           Asp or Glu.
#'            \item  \code{*_AA_AROMATIC}:   fraction of informative positions that are 
#'                                           His, Phe, Trp or Tyr.
#'            
#'          }
#'          
#'          Where \code{*} is the value from \code{label} or the name specified for 
#'          \code{seq} if \code{label=NULL}.
#'          
#' @details 
#' For all properties except for length, non-informative positions are excluded, 
#' where non-informative is defined as any character in \code{c("X", "-", ".", "*")}.
#' 
#' The scores for GRAVY, bulkiness and polarity are calculated as simple averages of the 
#' scores for each informative positions. The basic, acid and aromatic indices are 
#' calculated as the fraction of informative positions falling into the given category.
#' 
#' The aliphatic index is calculated using the Ikai, 1980 method.
#' 
#' The net charge is calculated using the method of Moore, 1985, excluding the N-terminus and
#' C-terminus charges, and normalizing by the number of informative positions.  The default 
#' pH for the calculation is 7.4.
#' 
#' The following data sources were used for the default property scores:
#' \itemize{
#'   \item  hydropathy:  Kyte & Doolittle, 1982.  
#'   \item  bulkiness:   Zimmerman et al, 1968. 
#'   \item  polarity:    Grantham, 1974.
#'   \item  pK:          EMBOSS.
#' }
#' 
#' @references
#' \enumerate{
#'   \item  Zimmerman JM, Eliezer N, Simha R. The characterization of amino acid sequences 
#'            in proteins by statistical methods. J Theor Biol 21, 170-201 (1968).
#'   \item  Grantham R. Amino acid difference formula to help explain protein evolution. 
#'            Science 185, 862-864 (1974).
#'   \item  Ikai AJ. Thermostability and aliphatic index of globular proteins. 
#'            J Biochem 88, 1895-1898 (1980).
#'   \item  Kyte J, Doolittle RF. A simple method for displaying the hydropathic character 
#'            of a protein. J Mol Biol 157, 105-32 (1982).
#'   \item  Moore DS. Amino acid and peptide net charges: A simple calculational procedure. 
#'            Biochem Educ 13, 10-11 (1985).
#'   \item  Wu YC, et al. High-throughput immunoglobulin repertoire analysis distinguishes 
#'            between human IgM memory and switched memory B-cell populations. 
#'            Blood 116, 1070-8 (2010).
#'   \item  Wu YC, et al. The relationship between CD27 negative and positive B cell 
#'            populations in human peripheral blood. 
#'            Front Immunol 2, 1-12 (2011).
#'   \item  \url{http://emboss.sourceforge.net/apps/cvs/emboss/apps/iep.html}
#' }
#' 
#' @seealso 
#' See \link{countPatterns} for counting the occurance of specific amino acid subsequences.
#' See \link{gravy}, \link{bulk}, \link{aliphatic}, \link{polar} and \link{charge} for functions 
#' that calculate the included properties individually.
#' 
#' @examples
#' # Load example data
#' file <- system.file("extdata", "ExampleDb.gz", package="alakazam")
#' df <- readChangeoDb(file)
#' df <- df[c(1,10,100), c("SEQUENCE_ID", "JUNCTION")]
#' 
#' # Calculate default amino acid properties from amino acid sequences
#' # Use a custom output column prefix
#' df$JUNCTION_TRANS <- translateDNA(df$JUNCTION)
#' aminoAcidProperties(df, seq="JUNCTION_TRANS", label="JUNCTION")
#
#' # Calculate default amino acid properties from DNA sequences
#' aminoAcidProperties(df, seq="JUNCTION", nt=TRUE)
#' 
#' # Use the Grantham, 1974 side chain volume scores from the seqinr package
#' # Set pH=7.0 for the charge calculation
#' # Calculate only average volume and charge
#' # Remove the head and tail amino acids from the junction, thus making it the CDR3
#' library(seqinr)
#' data(aaindex)
#' x <- aaindex[["GRAR740103"]]$I
#' # Rename the score vector to use single-letter codes
#' names(x) <- translateStrings(names(x), ABBREV_AA)
#' # Calculate properties
#' aminoAcidProperties(df, property=c("bulk", "charge"), seq="JUNCTION", nt=TRUE, 
#'                     trim=TRUE, label="CDR3", bulkiness=x, pH=7.0)
#'
#' @export
aminoAcidProperties <- function(data, property=c("length", "gravy", "bulk",
                                                 "aliphatic","polarity","charge",
                                                 "basic","acidic", "aromatic"),
                                seq="JUNCTION", nt=FALSE, trim=FALSE, label=NULL, ...) {
    # Check arguments
    property <- match.arg(property, several.ok=TRUE)
    
    # Define the data.frame that will be returned with amino acid properties
    prop_colnames <- list(
        "length"    = "AA_LENGTH",
        "gravy"     = "AA_GRAVY",
        "bulk"      = "AA_BULK",
        "aliphatic" = "AA_ALIPHATIC",
        "polarity"  = "AA_POLARITY",
        "charge"    = "AA_CHARGE",
        "basic"     = "AA_BASIC",
        "acidic"    = "AA_ACIDIC",
        "aromatic"  = "AA_AROMATIC"
    )
    # If no label, use sequence column name
    if (is.null(label)) { label <- seq }
    prop_colnames <- lapply(prop_colnames, function(x) { paste(label,x,sep="_") })
    
    out_df <- data.frame(matrix(NA, nrow=nrow(data), ncol=length(property)))
    colnames(out_df) <- prop_colnames[property]
    
    # Check if out_df column names already existed in data
    # if yes, throw warning
    check <- checkColumns(data, colnames(out_df))
    if (any(check == TRUE)) { warning("Duplicated columns found. Overwriting previous values.")}
    # Check input
    if (length(seq) > 1) {
        stop("You may specify only one sequence column; seq must be a vector of length 1.")
    }
    check <- checkColumns(data, seq)
    if (check != TRUE) { stop(check) }
    
    # Assign ellipsis arguments to correct function
    dots <- list(...)
    args_gravy <- dots[names(dots) %in% names(formals(gravy))]
    args_bulk <- dots[names(dots) %in% names(formals(bulk))]
    args_aliphatic <- dots[names(dots) %in% names(formals(aliphatic))]
    args_polar <- dots[names(dots) %in% names(formals(polar))]
    args_charge <- dots[names(dots) %in% names(formals(charge))]
    
    # Get sequence vector and translate if required
    region <- as.character(data[[seq]])
    region_aa <- if (nt) { translateDNA(region, trim=trim) } else { region }
    
    ## Will retrieve properties for valid sequences only
    ## keep index to fill results data.frame
    valid_seq <- isValidAASeq(region_aa)
    if (any(valid_seq == F) ){
        not_valid_num <- sum(!valid_seq)
        warning(paste0("Found ", not_valid_num , " sequences with non valid amino acid symbols"))
    }
    valid_seq_idx <- which(valid_seq)
    region_aa <- region_aa[valid_seq_idx]
    
    # Calculate region lengths
    if ("length" %in% property) {
        aa_length <- stri_length(region_aa)
        out_df[valid_seq_idx , prop_colnames$length] <- aa_length
    }
    # Average hydrophobicity
    if ("gravy" %in% property) {
        #aa_gravy <- gravy(region_aa, hydropathy)
        aa_gravy <- do.call('gravy', c(list(seq=region_aa), args_gravy))
        out_df[valid_seq_idx , prop_colnames$gravy] <- aa_gravy
    }
    # Average bulkiness
    if ("bulk" %in% property) {
        #aa_bulk <- bulk(region_aa)
        aa_bulk <- do.call('bulk', c(list(seq=region_aa), args_bulk))
        out_df[valid_seq_idx , prop_colnames$bulk] <- aa_bulk
    }
    if ("aliphatic" %in% property) {
        # Normalizes aliphatic index
        aa_aliphatic <- do.call('aliphatic', c(list(seq=region_aa), args_aliphatic))
        out_df[valid_seq_idx , prop_colnames$aliphatic] <- aa_aliphatic
    }
    # Average polarity
    if ("polarity" %in% property) {
        #aa_polarity <- polar(region_aa)
        aa_polarity <- do.call('polar', c(list(seq=region_aa), args_polar))
        out_df[valid_seq_idx , prop_colnames$polarity] <- aa_polarity
    }
    # Normalized net charge
    if ("charge" %in% property) {
        #aa_charge <- charge(region_aa)
        aa_charge <- do.call('charge', c(list(seq=region_aa), args_charge))
        out_df[valid_seq_idx , prop_colnames$charge] <- aa_charge
    }
    
    # Count of informative positions
    aa_info <-  stri_length(gsub("[X\\.\\*-]", "", region_aa))
    # Fraction of amino acid that are basic
    if ("basic" %in% property) {
        aa_basic <- countOccurrences(region_aa, "[RHK]") / aa_info
        out_df[valid_seq_idx , prop_colnames$basic] <- aa_basic
    }
    # Fraction of amino acid that are acidic
    if ("acidic" %in% property) {
        aa_acidic <- countOccurrences(region_aa, "[DE]") / aa_info
        out_df[valid_seq_idx , prop_colnames$acidic] <- aa_acidic
    }
    # Count fraction of aa that are aromatic
    if ("aromatic" %in% property) {
        aa_aromatic <- countOccurrences(region_aa, "[FWHY]") / aa_info
        out_df[valid_seq_idx , prop_colnames$aromatic] <- aa_aromatic
    }
    
    data_cols <- colnames(data) %in% colnames(out_df) == FALSE
    return(cbind(data[, data_cols], out_df))
}

