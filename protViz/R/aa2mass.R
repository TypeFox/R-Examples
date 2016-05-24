#R
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/aa2mass.R $
# $Id: aa2mass.R 6675 2014-09-17 11:12:25Z cpanse $
# $Date: 2014-09-17 13:12:25 +0200 (Wed, 17 Sep 2014) $


aa2mass <- function(peptideSequence, 
    mass=AA$Monoisotopic, 
    letter1=AA$letter1){

    return(.Call("aa2mass_main", peptideSequence, mass, letter1, PACKAGE = "protViz")$output)
}
