#' Plasmid dilution series raw data
#' 
#' These are the raw data from the \code{pds_raw} data set as measured by the
#' BioRad QX100 Droplet Digital PCR System.
#' 
#' The results can be as calculated by the BioRad QX100 Droplet Digital PCR
#' System are to be found in \code{pds}.
#' 
#' Setup: Duplex assay with constant amount of genomic DNA and six 10-fold
#' dilutions of plasmid DNA with 4 replicates, ranging theoretically from ~
#' 10^4 to 10^-1 copies/ micro L plus 4 replicates without plasmid DNA.
#' Included are No-gDNA-control and No-template-control, 2 replicates each.
#' 
#' Annotation: FX.Y (X = dilution number, Y = replicate number). Hardware:
#' Bio-Rad QX100 Droplet digital PCR system Details: Genomic DNA isolated from
#' Pseudomonas putida KT2440. Plasmid is pCOM10-StyA::EGFP StyB [Jahn et al.,
#' 2013, Curr Opin Biotechnol, Vol. 24 (1): 79-87]. Template DNA was heat
#' treated at 95 degree Celsius for 5 min prior to PCR. Channel 1, primers for
#' genomic DNA marker ileS, Taqman probes (FAM labelled). Channel 2, primers
#' for plasmid DNA marker styA, Taqman probes (HEX labelled).
#' 
#' 
#' @name pds_raw
#' @docType data
#' @details
#' The results can be as calculated by the BioRad QX100 Droplet Digital PCR
#' System are to be found in \code{pds}.
#'
#' Setup: Duplex assay with constant amount of genomic DNA and six 10-fold
#' dilutions of plasmid DNA with 4 replicates, ranging theoretically from ~
#' 10^4 to 10^-1 copies/ micro L plus 4 replicates without plasmid DNA.
#' Included are No-gDNA-control and No-template-control, 2 replicates each.
#' 
#' Annotation: FX.Y (X = dilution number, Y = replicate number). Hardware:
#'   Bio-Rad QX100 Droplet digital PCR system Details: Genomic DNA isolated from
#' Pseudomonas putida KT2440. Plasmid is pCOM10-StyA::EGFP StyB [Jahn et al.,
#'                                                               2013, Curr Opin Biotechnol, Vol. 24 (1): 79-87]. Template DNA was heat
#' treated at 95 degree Celsius for 5 min prior to PCR. Channel 1, primers for
#' genomic DNA marker ileS, Taqman probes (FAM labelled). Channel 2, primers
#' for plasmid DNA marker styA, Taqman probes (HEX labelled).
#' @format A list of 32 data frames.
#' \describe{ 
#'   \item{Well }{| ExptType }{| Experiment }{| Sample + Dilution step }{| TypeAssay }{| Assay}
#'   \item{A01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^4 }{| Ch1Unknown }{| ileS}
#'   \item{B01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^4 }{| Ch1Unknown }{| ileS}
#'   \item{C01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^4 }{| Ch1Unknown }{| ileS}
#'   \item{D01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^4 }{| Ch1Unknown }{| ileS}
#'   \item{E01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^3 }{| Ch1Unknown }{| ileS}
#'   \item{F01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^3 }{| Ch1Unknown }{| ileS}
#'   \item{G01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^3 }{| Ch1Unknown }{| ileS}
#'   \item{H01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^3 }{| Ch1Unknown }{| ileS}
#'   \item{A02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^2 }{| Ch1Unknown }{| ileS}
#'   \item{B02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^2 }{| Ch1Unknown }{| ileS}
#'   \item{C02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^2 }{| Ch1Unknown }{| ileS}
#'   \item{D02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^2 }{| Ch1Unknown }{| ileS}
#'   \item{E02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^1 }{| Ch1Unknown }{| ileS}
#'   \item{F02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^1 }{| Ch1Unknown }{| ileS}
#'   \item{G02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^1 }{| Ch1Unknown }{| ileS}
#'   \item{H02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^1 }{| Ch1Unknown }{| ileS}
#'   \item{A03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^0 }{| Ch1Unknown }{| ileS}
#'   \item{B03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^0 }{| Ch1Unknown }{| ileS}
#'   \item{C03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^0 }{| Ch1Unknown }{| ileS}
#'   \item{D03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^0 }{| Ch1Unknown }{| ileS}
#'   \item{E03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^-1 }{| Ch1Unknown }{| ileS}
#'   \item{F03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^-1 }{| Ch1Unknown }{| ileS}
#'   \item{G03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^-1 }{| Ch1Unknown }{| ileS}
#'   \item{H03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^-1 }{| Ch1Unknown }{| ileS}
#'   \item{A04 }{| Absolute Quantification }{| ABS }{| gDNA }{| Ch1Unknown }{| ileS}
#'   \item{B04 }{| Absolute Quantification }{| ABS }{| gDNA }{| Ch1Unknown }{| ileS}
#'   \item{C04 }{| Absolute Quantification }{| ABS }{| gDNA }{| Ch1Unknown }{| ileS}
#'   \item{D04 }{| Absolute Quantification }{| ABS }{| gDNA }{| Ch1Unknown }{| ileS}
#'   \item{E04 }{| Absolute Quantification }{| ABS }{| B + P 10^2 }{| Ch1NTC }{| ileS}
#'   \item{F04 }{| Absolute Quantification }{| ABS }{| B + P 10^2 }{| Ch1NTC }{| ileS}
#'   \item{G04 }{| Absolute Quantification }{| ABS }{| B }{| Ch1NTC }{| ileS}
#'   \item{H04 }{| Absolute Quantification }{| ABS }{| B }{| Ch1NTC }{| ileS}
#'   \item{A01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^4 }{| Ch2Unknown }{| styA}
#'   \item{B01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^4 }{| Ch2Unknown }{| styA}
#'   \item{C01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^4 }{| Ch2Unknown }{| styA}
#'   \item{D01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^4 }{| Ch2Unknown }{| styA}
#'   \item{E01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^3 }{| Ch2Unknown }{| styA}
#'   \item{F01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^3 }{| Ch2Unknown }{| styA}
#'   \item{G01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^3 }{| Ch2Unknown }{| styA}
#'   \item{H01 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^3 }{| Ch2Unknown }{| styA}
#'   \item{A02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^2 }{| Ch2Unknown }{| styA}
#'   \item{B02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^2 }{| Ch2Unknown }{| styA}
#'   \item{C02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^2 }{| Ch2Unknown }{| styA}
#'   \item{D02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^2 }{| Ch2Unknown }{| styA}
#'   \item{E02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^1 }{| Ch2Unknown }{| styA}
#'   \item{F02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^1 }{| Ch2Unknown }{| styA}
#'   \item{G02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^1 }{| Ch2Unknown }{| styA}
#'   \item{H02 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^1 }{| Ch2Unknown }{| styA}
#'   \item{A03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^0 }{| Ch2Unknown }{| styA}
#'   \item{B03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^0 }{| Ch2Unknown }{| styA}
#'   \item{C03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^0 }{| Ch2Unknown }{| styA}
#'   \item{D03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^0 }{| Ch2Unknown }{| styA}
#'   \item{E03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^-1 }{| Ch2Unknown }{| styA}
#'   \item{F03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^-1 }{| Ch2Unknown }{| styA}
#'   \item{G03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^-1 }{| Ch2Unknown }{| styA}
#'   \item{H03 }{| Absolute Quantification }{| ABS }{| gDNA + P 10^-1 }{| Ch2Unknown }{| styA}
#'   \item{A04 }{| Absolute Quantification }{| ABS }{| gDNA }{| Ch2NTC }{| styA}
#'   \item{B04 }{| Absolute Quantification }{| ABS }{| gDNA }{| Ch2NTC }{| styA}
#'   \item{C04 }{| Absolute Quantification }{| ABS }{| gDNA }{| Ch2NTC }{| styA}
#'   \item{D04 }{| Absolute Quantification }{| ABS }{| gDNA }{| Ch2NTC }{| styA}
#'   \item{E04 }{| Absolute Quantification }{| ABS }{| B + P 10^2 }{| Ch2Unknown }{| styA}
#'   \item{F04 }{| Absolute Quantification }{| ABS }{| B + P 10^2 }{| Ch2Unknown }{| styA}
#'   \item{G04 }{| Absolute Quantification }{| ABS }{| B }{| Ch2NTC }{| styA}
#'   \item{H04 }{| Absolute Quantification }{| ABS }{| B }{| Ch2NTC }{| styA}
#' }
#' @author Michael Jahn, Stefan Roediger, Michal Burdukiewcz
#' @references Jahn et al., 2013, \emph{Curr Opin Biotechnol}, Vol. 24 (1):
#' 79-87
#' 
#' Jahn M, Vorpahl C, Tuerkowsky D, Lindmeyer M, Buehler B, Harms H, et al. 
#' Accurate Determination of Plasmid Copy Number of Flow-Sorted Cells using 
#' Droplet Digital PCR. \emph{Anal Chem} 2014; 86:5969--76. doi:10.1021/ac501118v.
#' @source Michael Jahn Flow cytometry group / Environmental microbiology
#' Helmholtz Centre for Environmental Research - UFZ Permoserstrasse 15 / 04318
#' Leipzig / Germany phone +49 341 235 1318 michael.jahn [at] ufz.de /
#' www.ufz.de
#' @keywords datasets
#' @examples
#' 
#' #str(pds_raw)
#' bioamp(data = pds_raw[["A01"]], main = "Well A01", pch = 19)
#' 
NULL