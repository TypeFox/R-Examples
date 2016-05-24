#' @title Utilities and kinship information for Behavior Genetics and Developmental research using the NLSY.
#' 
#' @description Utilities and kinship information for Behavior Genetics and Developmental research using the NLSY.
#' 
#' @docType package
#' @name NlsyLinks-package
#' @aliases NlsyLinks
#' @note This package considers both Gen1 and Gen2 subjects.  "Gen1" refers to
#' subjects in the original NLSY79 sample
#' (\url{http://www.bls.gov/nls/nlsy79.htm}).  "Gen2" subjects are the
#' biological children of the Gen1 females -ie, those in the NLSY79 Children
#' and Young Adults sample (\url{http://www.bls.gov/nls/nlsy79ch.htm}).
#' 
#' The release version is available through \href{https://cran.r-project.org/package=NlsyLinks}{CRAN} by   
#' running \code{install.packages('NlsyLinks')}.  
#' The most recent development version is available through \href{https://github.com/LiveOak/NlsyLinks}{GitHub} by   
#' running 
#' \code{devtools::install_github} \code{(repo = 'LiveOak/NlsyLinks')}   
#' (make sure \href{https://cran.r-project.org/package=devtools}{devtools} is already installed).
#' If you're having trouble with the package, please install the development version.  If this doesn't solve
#' your problem, please create a \href{https://github.com/LiveOak/NlsyLinks/issues}{new issue}, or email Will.
#' 
#' @author 
#'   \href{http://scholar.google.com/citations?user=ffsJTC0AAAAJ}{William Howard Beasley} (\href{http://howardliveoak.com/}{Howard Live Oak LLC}, Norman)
#'   
#'  \href{http://www.vanderbilt.edu/psychological_sciences/bio/joe-rodgers}{Joseph Lee Rodgers} (Vanderbilt University, Nashville)
#'  
#'  \href{http://find.ouhsc.edu/Faculty.aspx?FacultyID=1041}{David Bard} (University of Oklahoma Health Sciences Center, OKC)
#'  
#'  Kelly Meredith (Oklahoma City University, OKC)
#'  
#'  \href{http://students.ou.edu/H/Michael.D.Hunter-1/}{Michael D. Hunter} (University of Oklahoma, Norman)
#' 
#' Maintainer: Will Beasley <wibeasley@@hotmail.com>
#' 
#' @references This package's development was largely supported by the NIH
#' Grant 1R01HD65865, \href{http://taggs.hhs.gov/AwardDetail.cfm?s_Award_Num=R01HD065865&n_Prog_Office_Code=50}{``NLSY Kinship Links: Reliable and Valid Sibling Identification"} (PI: Joe Rodgers).  A more complete list of research articles
#' using NLSY Kinship Links is maintained on our \href{http://liveoak.github.io/NlsyLinks/research-publications.html}{package's website}.
#' 
#' Rodgers, Joseph Lee, & Kohler, Hans-Peter (2005).  \href{http://www.springerlink.com/content/n3x1v1q282583366/}{Reformulating and
#' simplifying the DF analysis model.}
#' \emph{Behavior Genetics, 35} (2), 211-217.
#' 
#' Rodgers, J.L., Bard, D., Johnson, A., D'Onofrio, B., & Miller, W.B. (2008).
#' \href{http://www.ncbi.nlm.nih.gov/pubmed/18825497}{The Cross-Generational Mother-Daughter-Aunt-Niece Design: Establishing
#' Validity of the MDAN Design with NLSY Fertility Variables.} \emph{Behavior
#' Genetics, 38}, 567-578.
#' 
#' D'Onofrio, B.M., Van Hulle, C.A., Waldman, I.D., Rodgers, J.L., Rathouz,
#' P.J., & Lahey, B.B. (2007). \href{http://www.ncbi.nlm.nih.gov/pubmed/17984398}{Causal inferences regarding prenatal alcohol
#' exposure and childhood externalizing problems.} \emph{Archives of General
#' Psychiatry, 64}, 1296-1304.
#' 
#' Rodgers, J.L. & Doughty, D. (2000).  \href{http://link.springer.com/chapter/10.1007\%2F978-1-4615-4467-8_6}{Genetic and environmental influences
#' on fertility expectations and outcomes using NLSY kinship data.}  In J.L.
#' Rodgers, D. Rowe, & W.B. Miller (Eds.) \emph{Genetic influences on
#' fertility and sexuality.} Boston: Kluwer Academic Press.
#' 
#' Cleveland, H.H., Wiebe, R.P., van den Oord, E.J.C.G., & Rowe, D.C. (2000).
#' \href{http://www.ncbi.nlm.nih.gov/pubmed/10953940}{Behavior problems among children from different family structures: The
#' influence of genetic self-selection.} \emph{Child Development, 71}, 733-751.
#' 
#' Rodgers, J.L., Rowe, D.C., & Buster, M. (1999). \href{http://www.ncbi.nlm.nih.gov/pubmed/10081235}{Nature, nurture, and first
#' sexual intercourse in the USA: Fitting behavioural genetic models to NLSY
#' kinship data.} \emph{Journal of Biosocial Sciences, 31}.
#' 
#' Rodgers, J.L., Rowe, D.C., & Li, C. (1994). \href{http://psycnet.apa.org/journals/dev/30/3/374/}{Beyond nature versus nurture:
#' DF analysis of nonshared influences on problem behaviors.}
#' \emph{Developmental Psychology, 30}, 374-384.
#' @keywords package
#' @examples
#' 
#' library(NlsyLinks) #Load the package into the current R session.
#' summary(Links79Pair)  #Summarize the five variables.
#' hist(Links79Pair$R)  #Display a histogram of the Relatedness values.
#' table(Links79Pair$R)  #Create a table of the Relatedness values for the whole sample.
#' 
#' \dontrun{
#' # Install/update NlsyLinks with the release version from CRAN.
#' install.packages('NlsyLinks')
#' 
#' # Install/update NlsyLinks with the development version from GitHub
#' #install.packages('devtools') #Uncomment if `devtools` isn't installed already.
#' devtools::install_github('LiveOak/NlsyLinks')
#' }
NULL
