#' Data sets from The Analysis of Biological Data
#' 
#' The \code{abd} package contains data sets and sample code for the book,
#' \emph{The Analysis of Biological Data} by Michael C. Whitlock and Dolph
#' Schluter (2009; Roberts and Company Publishers).
#' 
#' \tabular{ll}{ Package: \tab abd\cr 
#'               Type: \tab Package\cr 
#'               Version: \tab 0.2-8\cr 
#'               Date: \tab 2015-07-02\cr 
#'               License: \tab GPL\cr 
#'               LazyLoad: \tab yes\cr 
#'               LazyData: \tab yes\cr }
#' 
#' @name abd-package
#' @aliases abd-package abd
#' @docType package
#' @author Kevin M. Middleton (\url{middletonk@@missouri.edu}); Randall Pruim
#' (\url{rpruim@@calvin.edu})
#' @references Whitlock, M.C. and D. Schluter. 2009. \emph{The Analysis of
#' Biological Data}. Roberts and Company Publishers. ISBN: 0981519407.
#' \url{http://www.roberts-publishers.com/biology/the-analysis-of-biological-data.html}
#' @keywords package
#' @import nlme lattice grid mosaic
#' @importFrom grDevices colorRampPalette gray
#' @examples
#' 
#' trellis.par.set(theme=col.abd())  # set color theme
#' show.settings()
#' abdData(3)                        # look for data sets in chapter 3
#' abdData('Finch')                  # look for data sets with 'finch' in name
#' 
NULL


#' \code{abd} Data Sets
#' 
#' Information about the location of data sets in \emph{Analysis of Biological
#' Data}
#' 
#' 
#' @name dataInfo
#' @docType data
#' @format A data frame with 143 observations on the following 5 variables.
#' \describe{ \item{name}{name of data set}
#' \item{chapter}{chapter in which data set appears}
#' \item{type}{used in an \code{Example} or a \code{Problem}}
#' \item{number}{example or problem number}
#' \item{sub}{sub-problem: \code{} \code{a} \code{b} \code{c}} }
#' @seealso \code{\link{abdData}}
#' @keywords datasets
#' @examples
#' 
#' str(dataInfo)
#' 
NULL


#' Carbon Dioxide and Growth Rate in Algae
#' 
#' Growth rates of the unicellular alga \emph{Chlamydomonas} after 1,000
#' generations of selection under \code{High} and \code{Normal} levels of
#' carbon dioxide.
#' 
#' 
#' @name AlgaeCO2
#' @docType data
#' @format A data frame with 14 observations on the following 2 variables.
#' \describe{ \item{treatment}{a factor with levels \code{Normal} and
#' \code{High}} \item{growthrate}{a numeric vector} }
#' @source Collins, S. and G. Bell. 2004. Phenotypic consequences of 1,000
#' generations of selection at elevated CO\eqn{_{2}}{2} in a green alga.
#' \emph{Nature} 431: 566-569.
#' @keywords datasets
#' @examples
#' 
#' AlgaeCO2
#' xyplot(growthrate ~ treatment, AlgaeCO2, type = c('p', 'a'))
#' 
NULL





#' Antilles Bird Immigration Dates
#' 
#' Approximate dates of immigration for 37 species of birds in the Lesser
#' Antilles.
#' 
#' 
#' @name Antilles
#' @docType data
#' @format A data frame with 37 observations of one variable. \describe{
#' \item{immigration.date}{approximate immigration date (in millions of
#' years)} }
#' @references
#' \url{http://www.sciencemag.org/cgi/content/abstract/sci;294/5546/1522}
#' @source \emph{inferred from} Ricklefs, R.E. and E. Bermingham. 2001.
#' Nonequilibrium diversity dynamics of the Lesser Antillean avifauna.
#' \emph{Science} 294: 1522-1524.
#' @keywords datasets
#' @examples
#' 
#' histogram(~immigration.date, Antilles,n=15)
#' densityplot(~immigration.date, Antilles)
#' 
NULL





#' Effects of Aspirin on Cancer Rates
#' 
#' Frequency of cancer in 39,876 women taking and not taking aspirin.
#' 
#' 
#' @name Aspirin
#' @docType data
#' @format A data frame with 39876 observations on the following \describe{
#' \item{treatment}{a factor with levels \code{Aspirin} and
#' \code{Placebo}} \item{cancer}{a factor with levels \code{no} and
#' \code{yes}} }
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/15998890}
#' @source Cook, N.R., I. Lee, J.M. Gaziano, D. Gordon, P.M. Ridker, J.E.
#' Manson, C.H. Hennekens, and J.E. Buring. 2005. Low-dose aspirin in the
#' primary prevention of cancer. \emph{Journal of the American Medical
#' Association} 294: 47-55.
#' @keywords datasets
#' @examples
#' 
#' demo(sec9.2)
#' 
NULL





#' Foraging Gene Expression
#' 
#' Levels of expression of the foraging gene (\emph{for}; \code{Expression}) in
#' two worker types (\code{type}) in three bee colonies (\code{colony}). Note
#' that \code{colony} is not coded as a factor.
#' 
#' 
#' @name BeeGenes
#' @docType data
#' @format A data frame with 6 observations on the following 3 variables.
#' \describe{ \item{type}{a factor with levels \code{forager}
#' \code{nurse}} \item{colony}{a numeric identifier}
#' \item{expression}{expression level of the \emph{for} gene} }
#' @source Ben-Shahar, Y., A. Robichon, M.B. Sokolowski, and G.E. Robinson.
#' 2002. Influence of gene action across different time scales on behavior.
#' \emph{Science} 296: 741-744.
#' @keywords datasets
#' @examples
#' 
#' str(BeeGenes)
#' BeeGenes
#' xtabs( expression ~ type + colony, BeeGenes )
#' 
NULL





#' Bee Lifespans
#' 
#' Lifespan of 33 foraging honey bees.
#' 
#' 
#' @name BeeLifespans
#' @docType data
#' @format A data frame with 33 observations on the following variable.
#' \describe{ \item{hours}{a numeric vector} }
#' @source \emph{inferred from} Visscher, P.K. and R. Dukas. 1997. Survivorship
#' of foraging honey bees. \emph{Insectes Sociaux} 44: 1-5.
#' @keywords datasets
#' @examples
#' 
#' histogram(~hours, BeeLifespans, n=10)
#' densityplot(~hours, BeeLifespans)
#' 
NULL





#' Beetle Wings and Horns
#' 
#' Relative size of the horns and wings in 19 female \emph{Onthophagus
#' sagittarius} beetles.
#' 
#' 
#' @name Beetles
#' @docType data
#' @format A data frame with 19 observations on the following 2 variables.
#' \describe{ \item{horn.size}{a numeric vector}
#' \item{wing.mass}{a numeric vector} }
#' @references
#' \url{http://www.scienceonline.org/cgi/content/abstract/291/5508/1534}
#' @source Emlen, D.J. 2001. Costs and the diversification of exaggerated
#' animal structures. \emph{Science} 291: 1534-1536.
#' @keywords datasets
#' @examples
#' 
#' str(Beetles)
#' xyplot(wing.mass ~ horn.size, Beetles)
#' 
NULL





#' Sex Ratios in Birds
#' 
#' Correlation coefficient of sex ratio in bird offspring.
#' 
#' 
#' @name BirdSexRatio
#' @docType data
#' @format A data frame with 15 observations of one variable \describe{
#' \item{corr.coeff}{correlation coefficient of sex ratio in bird
#' offspring} }
#' @source West, S.A. and B.C. Sheldon. 2002. Constraints in the evolution of
#' sex ratio adjustment. \emph{Science} 295: 1695-1688.
#' @keywords datasets
#' @examples
#' 
#' histogram(~corr.coeff, BirdSexRatio, n = 10,
#'   xlab = "Correlation Coefficient")
#' 
NULL





#' Testosterone Levels in Blackbirds
#' 
#' Experimental manipulation of testosterone levels in male Red-winged
#' Blackbirds (\emph{Agelaius phoeniceus}) and resulting changes in antibody
#' levels
#' 
#' 
#' @name Blackbirds
#' @docType data
#' @format A data frame with 13 observations on the following 6 variables.
#' \describe{ \item{before}{a numeric vector} \item{after}{a
#' numeric vector} \item{log.before}{a numeric vector}
#' \item{log.after}{a numeric vector} \item{diff.in.logs}{a
#' numeric vector} \item{diff}{a numeric vector} }
#' @source Hasselquist, D., J.A. Marsh, P.W. Sherman, and J.C. Wingfield. 1999.
#' Is avian immunocompetence suppressed by testosterone? \emph{Behavioral
#' Ecology and Sociobiology} 45: 167-175.
#' @keywords datasets
#' @examples
#' 
#' Blackbirds
#' xyplot(log.after ~ log.before, data = Blackbirds,
#'   ylab = "log Antibody production after implant",
#'   xlab = "log Antibody production before implant"
#' )
#' 
NULL





#' Heat Loss and Body Fat
#' 
#' Heat loss during exercise and relative body fat in 12 boys.
#' 
#' 
#' @name BodyFatHeatLoss
#' @docType data
#' @format A data frame with 12 observations on the following 2 variables.
#' \describe{ \item{leanness}{a numeric vector}
#' \item{lossrate}{a numeric vector} }
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/4732330}
#' @source Sloan, R.E.G. and W.R. Keatinge. 1973. Cooling rates of young people
#' swimming in cold water. \emph{Journal of Applied Physiology} 35: 371-375.
#' @keywords datasets
#' @examples
#' 
#' xyplot(lossrate ~ leanness, BodyFatHeatLoss)
#' 
NULL





#' Proteolipid Protein 1 Gene Expression
#' 
#' Expression levels of the proteolipid protein 1 gene (PLP1;
#' \code{PLP1.expression}) in 45 individuals in one of three \code{group}s.
#' 
#' 
#' @name BrainExpression
#' @docType data
#' @format A data frame with 45 observations on the following 2 variables.
#' \describe{ \item{group}{a factor with levels: \code{bipolar},
#' \code{control}, and \code{schizo}} \item{PLP1.expression}{a numeric
#' vector} }
#' @source \emph{inferred from} Tkachev, D., M.L. Mimmack, M.M. Ryan, M.
#' Wayland, T. Freeman, P.B. Jones, M. Starkey, M.J. Webster, R.H. Yolken, S.
#' Bahn. 2003. Oligodendrocyte dysfunction in schizophrenia and bipolar
#' disorder. \emph{Lancet} 362(9386): 798-805.
#' @keywords datasets
#' @examples
#' 
#' bwplot(PLP1.expression ~ group, BrainExpression)
#' 
NULL





#' Salmon Survival in the Presence of Brook Trout
#' 
#' Total numbers of salmon released (\code{salmon.released}) and surviving
#' (\code{salmon.surviving}) in 12 streams, 6 with brook trout \code{present}
#' and 6 with brook trout \code{absent}. The proportion of salmon surviving
#' (\code{proportion.surviving}) is given for each stream.
#' 
#' 
#' @name BrookTrout
#' @aliases BrookTrout BrookTrout2
#' @docType data
#' @format \code{BrookTrout} is a data frame with 12 observations on the
#' following 4 variables.  \code{BrookTrout2} is a different summary of the
#' same study and gives survival rates for chinook in different years.
#' \describe{ \item{trout}{a factor with levels \code{absent} and
#' \code{present} indicating whether brook trout are absent or present in the
#' stream} \item{salmon.released}{a numeric vector of the total number
#' of salmon released} \item{salmon.surviving}{a numeric vector of the
#' number of salmon surviving} \item{proportion.surviving}{a numeric
#' vector of the proportion of salmon surviving} }
#' @source Levin, P.S., S. Achord, B.E. Fiest, and R.W. Zabel. 2002.
#' Non-indigenous brook trout and the demise of Pacific salmon: a forgotten
#' threat? \emph{Proceedings of the Royal Society of London, Series B,
#' Biological Sciences} 269: 1663-1670.
#' @keywords datasets
#' @examples
#' 
#' str(BrookTrout)
#' str(BrookTrout2)
#' 
#' bwplot(proportion.surviving ~ trout, BrookTrout)
#' 
#' aggregate(proportion.surviving ~ trout, BrookTrout, FUN = favstats)
#' summary(proportion.surviving ~ trout, BrookTrout, fun = favstats)
#' 
NULL





#' Deaths from Horse Kicks
#' 
#' Numbers of deaths resulting from horse kicks per regiment-years for the
#' Prussian army.
#' 
#' 
#' @name Cavalry
#' @docType data
#' @format A data frame with 5 observations on the following 2 variables.
#' \describe{ \item{deaths}{a numeric vector} \item{count}{a
#' numeric vector} }
#' @source Bortkiewicz, L. 1898. Das Gesetz der Kleinen Zahlen (Teubner,
#' Leipzig), \emph{as cited in} Larson, R.J. and M.L. Marx. 1981. \emph{An
#' Introduction to Mathematical Statistics and its Applications}.
#' Prentice-Hall: Englewood Cliffs, NJ.
#' @keywords datasets
#' @examples
#' 
#' Cavalry
#' xyplot(count ~ deaths, Cavalry, type='h', lwd=4)
#' barchart(count ~ deaths, Cavalry, horizontal = FALSE,
#'          box.ratio = 1000, origin=0)
#' 
NULL





#' Alarm Calls in Chickadees
#' 
#' Number of "dee" notes per call in Black-capped Chickadees (\emph{Poecile
#' atricapilla}) for 13 predator species with differing body masses.
#' 
#' 
#' @name Chickadees
#' @docType data
#' @format A data frame with 13 observations on the following 3 variables.
#' \describe{ \item{species}{a character vector} \item{mass}{a
#' numeric vector} \item{dees}{a numeric vector} }
#' @references \url{http://www.sciencemag.org/cgi/content/short/308/5730/1934}
#' @source Templeton, C.N., E. Greene, and K. Davis. 2005. Allometry of alarm
#' calls: Black-capped Chickadees encode information about predator size.
#' \emph{Science} 308: 1934-1937.
#' @keywords datasets
#' @examples
#' 
#' str(Chickadees)
#' Chickadees
#' 
#' xyplot(dees ~ mass, data = Chickadees,
#'    xlab = "Predator body mass (kg)",
#'    ylab = "'Dees' per call", type=c('p','r')
#' )
#' 
NULL





#' Brodmann's Area 44 in Chimps
#' 
#' Asymmetry of Brodmann's area 44 in 20 chimpanzees.
#' 
#' 
#' @name ChimpBrains
#' @docType data
#' @format A data frame with 20 observations on the following 3 variables.
#' \describe{ \item{name}{name of chimp} \item{sex}{a factor
#' with levels \code{F} and \code{M}} \item{asymmetry}{asymmetry score}
#' }
#' @source Cantalupo, C. and W.D. Hopkins. 2001. Asymmetric Broca's area in
#' great apes. \emph{Nature} 414: 505.
#' @keywords datasets
#' @examples
#' 
#' xyplot(asymmetry ~ sex, ChimpBrains)
#' aggregate(asymmetry ~ sex, ChimpBrains, FUN = favstats)
#' summary(asymmetry ~ sex, ChimpBrains, fun = favstats)
#' 
NULL





#' Cichlid Mating Preference
#' 
#' Preference index in F1 and F2 crosses of two species of cichlids from Lake
#' Victoria, \emph{Pundamilia pundamilia} and \emph{P. nyererei}.
#' 
#' 
#' @name Cichlids
#' @docType data
#' @format A data frame with 53 observations on the following 2 variables.
#' \describe{ \item{genotype}{a factor with levels \code{F1} and
#' \code{F2}} \item{preference}{a numeric vector} }
#' @references
#' \url{http://rspb.royalsocietypublishing.org/content/272/1560/237.full.pdf}
#' @source Haeslery, M.P. and O. Seehausen. 2005. Inheritance of female mating
#' preference in a sympatric sibling species pair of Lake Victoria cichlids:
#' implications for speciation. \emph{Proceedings of the Royal Society of
#' London, Series B, Biological Sciences} 272: 237-245.
#' @keywords datasets
#' @examples
#' 
#' str(Cichlids)
#'
#' summary(preference ~ genotype, Cichlids, fun = favstats)
#' 
#' if (require(plyr)) {
#' ddply(Cichlids, .(genotype),
#'   function(df)c(mean = mean(df$preference),
#'                 standard.deviation = sd(df$preference),
#'                 n = length(df$preference)))
#' }
#' 
NULL





#' GnRH Levels in Cichlids
#' 
#' Levels of mRNA for gonadotropin-releasing hormone in cichlids
#' (\emph{Haplochromis burtoni}) that are (\emph{n} = 5) and are not (\emph{n}
#' = 6) territorial.
#' 
#' 
#' @name CichlidsGnRH
#' @docType data
#' @format A data frame with 11 observations on the following 2 variables.
#' \describe{ \item{territorial}{a factor with levels \code{No} and
#' \code{Yes}} \item{GnRH.mRNA}{a numeric vector} }
#' @references \url{http://jeb.biologists.org/cgi/content/abstract/205/17/2567}
#' @source White, S.A., T. Nguyen, and R.D. Fernald. 2002. Social regulation of
#' gonadotropin-releasing hormone. \emph{Journal of Experimental Biology} 205:
#' 2567-2581.
#' @keywords datasets
#' @examples
#' 
#' xyplot(GnRH.mRNA ~ territorial, CichlidsGnRH, type=c('p','a'))
#' 
NULL





#' Biomass Change in Rainforests near Clearcuts
#' 
#' Biomass change in 36 Amazonian rainforests following clearcuts ranging from
#' 50 m to several kilometers.
#' 
#' 
#' @name Clearcuts
#' @docType data
#' @format A data frame with 36 observations of one variable. \describe{
#' \item{biomass.change}{} }
#' @references
#' \url{http://www.sciencemag.org/cgi/content/abstract/278/5340/1117}
#' @source Laurance, W.F., S.G. Laurance, L.V. Ferreira, J.M. Rankin-de Merona,
#' C. Gascon, T.E. Lovejoy. 1997. Biomass collapse in Amazonian forest
#' fragments. \emph{Science} 278: 1117-1118.
#' @keywords datasets
#' @examples
#' 
#' str(Clearcuts)
#' histogram(~biomass.change, Clearcuts)
#' 
NULL





#' Effects of Cocaine on Dopamine Receptors
#' 
#' Percent of dopamine receptors blocked (\code{percent.bocked}) and the
#' perceived level of high as determined by PET scans (\code{high}) in 34
#' humans.
#' 
#' 
#' @name CocaineDopamine
#' @docType data
#' @format A data frame with 34 observations on the following 2 variables.
#' \describe{ \item{percent.blocked}{a numeric vector}
#' \item{high}{a numeric vector} }
#' @references
#' \url{http://www.nature.com/nature/journal/v386/n6627/abs/386827a0.html}
#' @source Volkow, N.D., G.-J. Wang, R.W. Foltin, J.S. Fowler, N.N. Abumrad, S.
#' Vitkun, J. Logan, S.J. Gatley, N. Pappas, R. Hitzemann, and C.E. Shea. 1997.
#' Relationship between subjective effects of cocaine and dopamine transporter
#' occupancy. \emph{Nature} 386: 827-830.
#' @keywords datasets
#' @examples
#' 
#' str(CocaineDopamine)
#' xyplot(high ~ percent.blocked, CocaineDopamine)
#' 
NULL





#' Frequency of Convictions for a Cohort of English Boys
#' 
#' Data on frequency of convictions for a cohort of 395 boys.
#' 
#' 
#' @name Convictions
#' @docType data
#' @format A data frame with 15 observations on the following 2 variables.
#' \describe{ \item{convictions}{number of convictions}
#' \item{boys}{number of boys with given number of convictions} }
#' @references \url{http://www.icpsr.umich.edu/icpsrweb/NACJD/archive.jsp}
#' @source Farrington, D.P. 1994. \emph{Cambridge Study in Delinquent
#' Development} [Great Britain], 1961-1981. 2nd ICPSR ed. Inter-university
#' Consortium for Political and Social Research, Ann Arbor, MI.
#' @keywords datasets
#' @examples
#' 
#' str(Convictions)
#' barchart(boys ~ as.factor(convictions), Convictions, horizontal = FALSE, origin=0)
#' xyplot( boys ~ convictions, Convictions, type = "h", lwd = 20)
#' 
NULL





#' Convictions and Income Level in a Cohort of English Boys
#' 
#' Data reporting the number of individual with and without convictions per
#' income level.
#' 
#' 
#' @name ConvictionsAndIncome
#' @docType data
#' @format A data frame with 395 observations on the following 2 variables.
#' \describe{ \item{convicted}{a factor with levels \code{no} and
#' \code{yes}} \item{income}{a factor with levels \code{adequate},
#' \code{comfortable}, and \code{inadequate}} }
#' @references \url{http://www.icpsr.umich.edu/icpsrweb/NACJD/archive.jsp}
#' @source Farrington, D.P. 1994. \emph{Cambridge Study in Delinquent
#' Development} [Great Britain], 1961-1981. 2nd ICPSR ed. Inter-university
#' Consortium for Political and Social Research, Ann Arbor, MI.
#' @keywords datasets
#' @examples
#' 
#' str(ConvictionsAndIncome)
#' ConvictionsAndIncome
#' 
#' xtabs(~ convicted + income, data = ConvictionsAndIncome)
#' 
NULL





#' Immunity and Sperm Viability in Crickets
#' 
#' Sperm viability and immune function, measured by lysozyme activity in
#' crickets. Each observation is a mean for a single family of males.
#' 
#' 
#' @name Crickets
#' @docType data
#' @format A data frame with 41 observations on the following 2 variables.
#' \describe{ \item{sperm.viability}{a numeric vector}
#' \item{lysozyme}{a numeric vector} }
#' @source Simmons, L.W. and B. Roberts. 2005. Bacterial immunity traded for
#' sperm viability in male crickets. \emph{Science} 309: 2031.
#' @keywords datasets
#' @examples
#' 
#' Crickets
#' xyplot(lysozyme ~ sperm.viability, Crickets)
#' 
NULL





#' DEET and Mosquito Bites
#' 
#' Administered dose of DEET and number of mosquito bites for 52 women.
#' 
#' 
#' @name DEET
#' @docType data
#' @format A data frame with 52 observations on the following 2 variables.
#' \describe{ \item{dose}{a numeric vector} \item{bites}{a
#' numeric vector} }
#' @source Golenda, C.F., V.B. Solberg, R. Burge, J.M. Gambel, and R.A. Wirtz.
#' 1999. Gender-related efficacy difference to an extended duration formulation
#' of topical N,N-diethyl-\emph{m}-toluamide (DEET). \emph{American Journal of
#' Tropical Medicine and Hygiene} 60: 654-657.
#' @keywords datasets
#' @examples
#' 
#' str(DEET)
#' xyplot(bites ~ dose, DEET)
#' 
NULL





#' Daphnia Longevity
#' 
#' Number of spores and host longevity in the crustacean \emph{Daphnia magna}.
#' 
#' 
#' @name DaphniaLongevity
#' @docType data
#' @format A data frame with 32 observations on the following 2 variables.
#' \describe{ \item{longevity}{a numeric vector}
#' \item{sqrt.spores}{a numeric vector} }
#' @references
#' \url{http://www.plosbiology.org/article/info:doi/10.1371/journal.pbio.0040197}
#' @source Jensen, K.H., T.J. Little, A. Skorping, and D. Ebert. 2006.
#' Empirical support for optimal virulence in a castrating parasite. \emph{PLoS
#' Biology} 4(7): e197
#' @keywords datasets
#' @examples
#' 
#' str(DaphniaLongevity)
#' xyplot(sqrt.spores ~ longevity, DaphniaLongevity)
#' 
NULL





#' Daphnia Resistance to Cyanobacteria
#' 
#' Resistance of \emph{Daphnia} eggs to different levels of cyanobacteria
#' (\code{cyandensity}) from 1962-1997.
#' 
#' 
#' @name DaphniaResistance
#' @docType data
#' @format A data frame with 32 observations on the following 2 variables.
#' \describe{ \item{density}{a factor with levels: \code{high},
#' \code{low}, and \code{med}} \item{resistance}{a numeric vector} }
#' @source \emph{inferred from} Hairston, N.G., Jr., W. Lampert, C.E. Cáceres,
#' C.L. Holtmeier, L.J. Weider, U. Gaedke, J.M. Fischer, J.A. Fox, and D.M.
#' Post. 1999. Dormant eggs record rapid evolution. \emph{Nature} 401: 446.
#' @keywords datasets
#' @examples
#' 
#' str(DaphniaResistance)
#' 
#' 
#' bwplot(resistance ~ density, DaphniaResistance)
#' # with such a small data set, we can display all the data
#' # rather than a summary
#' xyplot(resistance ~ density, DaphniaResistance)
#' histogram( ~ resistance| density, DaphniaResistance, 
#' 	strip=FALSE, strip.left = TRUE,
#' 	layout=c(1,3)
#' 	)
#' 
NULL





#' Day of Birth
#' 
#' Day of the week for 350 U.S. births in 1999.
#' 
#' 
#' @name DayOfBirth
#' @docType data
#' @format A data frame with 7 observations on the following 2 variables.
#' \describe{ \item{day}{a character vector} \item{births}{a
#' numeric vector} }
#' @references \url{http://cdc.gov/NCHS/products/nvsr.htm}
#' @source Ventura, S.J., J.A. Martin, S.C. Curtin, F. Menacker, and B.E.
#' Hamilton. 2001. Births: final data for 1999. \emph{National Vital Statistics
#' Reports} Vol. 49, No. 1.
#' @keywords datasets
#' @examples
#' 
#' DayOfBirth
#' barchart( day ~ births, DayOfBirth, origin=0)
#' 
#' # Fix bad ordering of days
#' DayOfBirth$oday <- with(DayOfBirth, ordered(day, levels = day))
#' barchart( oday ~ births, DayOfBirth, origin=0)
#' barchart( births ~ oday, DayOfBirth, horizontal = FALSE, origin=0)
#' barchart( births ~ oday, DayOfBirth, horizontal = FALSE, origin=0, 
#'  scales = list(x=list(rot=45)))
#' 
#' barplot(DayOfBirth$births,
#'   ylim = c(0, 70),
#'   names.arg = DayOfBirth$day,
#'   las = 2,
#'   mgp = c(3, 0.75, 0))
#' 
#' 
NULL





#' Desert Bird Census Data
#' 
#' Census data for desert birds.
#' 
#' 
#' @name DesertBirds
#' @docType data
#' @format A data frame with 43 observations on the following 2 variables.
#' \describe{ \item{species}{a character vector} \item{count}{a
#' numeric vector} }
#' @references \url{http://www.mbr-pwrc.usgs.gov/bbs/}
#' @source Sauer, J.R., J.E. Hines, and J. Fallon. 2003. The North American
#' breeding bird survey, results and analysis 1966-2002. Version 2003.1. USGS
#' Patuxent Wildlife Research Center, Laurel, MD.
#' @keywords datasets
#' @examples
#' 
#' str(DesertBirds)
#' histogram(~ count, DesertBirds,
#'   xlab = "Abundance"
#'   )
#' 
NULL





#' Dioecy vs. Monomorphism in Plants
#' 
#' Number of \code{dioecious} and \code{monomorphic} taxa among pairs of
#' closely related plants.
#' 
#' 
#' @name Dioecy
#' @docType data
#' @format A data frame with 28 observations on the following 3 variables.
#' \describe{ \item{dioecious}{a numeric vector}
#' \item{monomorphic}{a numeric vector}
#' \item{taxon.pair}{identifier for pair} }
#' @source Heilbuth, J.C. 2000. Lower species richness in dioecious clades.
#' \emph{The American Naturalist} 156: 221-241.
#' @keywords datasets
#' @examples
#' 
#' xyplot(dioecious ~ monomorphic, Dioecy, alpha = 0.65, pch = 16)
#' 
NULL





#' Dolphin Swimming Behavior
#' 
#' Percentage of time 8 sleeping dolphins from the Southern Hemisphere spent
#' swimming clockwise.
#' 
#' 
#' @name Dolphins
#' @docType data
#' @format A data frame with 8 observations on one variable. \describe{
#' \item{percent.clockwise}{percent of time spent swimming clockwise
#' while sleeping.} }
#' @references \url{http://faculty.washington.edu/chudler/dolp.html}
#' @source Stafne, G.M. and P.R. Manger. 2004. Predominance of clockwise
#' swimming during rest in Southern Hemisphere dolphins. \emph{Physiology and
#' Behavior} 82: 919-926.
#' @keywords datasets
#' @examples
#' 
#' Dolphins
#' hist(Dolphins$percent.clockwise)
#' histogram(~ percent.clockwise, Dolphins)
#' 
NULL





#' Heritability of Body Condition in Dung Beetles
#' 
#' Body condition (\code{offspring.condition}) in 36 dung beetles
#' (\emph{Onthophagus taurus}) from 12 \code{male}s each mated to 3 different
#' virgin females.
#' 
#' 
#' @name DungBeetles
#' @docType data
#' @format A data frame with 36 observations on the following 2 variables.
#' \describe{ \item{id}{a numeric vector}
#' \item{offspring.condition}{a numeric vector} }
#' @references \url{http://en.wikipedia.org/wiki/Dung_beetle}
#' 
#' \url{http://www.nature.com/nature/journal/v410/n6829/abs/410684a0.html}
#' @source \emph{inferred from} Kotiaho, J.S., L.W. Simmons, and J.L. Tomkins.
#' 2001. Towards a resolution of the lek paradox. \emph{Nature} 410: 684-686.
#' @keywords datasets
#' @examples
#' 
#' str(DungBeetles)
#' xyplot(offspring.condition ~ factor(id), DungBeetles, 
#'   xlab='Dung Beetle', 
#'   ylab='offspring condition')
#' 
NULL





#' Earthworm Diversity and Soil Nitrogen Levels
#' 
#' Number of earthworm species and total nitrogen content in the soil in 39
#' hardwood forest plots.
#' 
#' 
#' @name Earthworms
#' @docType data
#' @format A data frame with 39 observations on the following 2 variables.
#' \describe{ \item{worm.species}{a numeric vector}
#' \item{nitrogen}{a numeric vector} }
#' @references
#' \url{http://www3.interscience.wiley.com/journal/118701215/abstract}
#' @source Gundale, M.J., W.M. Jolly, and T.H. Deluca. 2005. Susceptibility of
#' a northern hardwood forest to exotic earthworm invasion. \emph{Conservation
#' Biology} 19: 1075-1083.
#' @keywords datasets
#' @examples
#' 
#' str(Earthworms)
#' xyplot(nitrogen ~ worm.species, Earthworms)
#' 
#' 
NULL





#' Earwig Density and Forceps
#' 
#' Earwig (\emph{Forficula auricularia}) density and the proportion of trapped
#' earwigs with abdominal forceps (used for fighting and courtship).
#' 
#' 
#' @name Earwigs
#' @docType data
#' @format A data frame with 7 observations on the following 2 variables.
#' \describe{ \item{density}{a numeric vector}
#' \item{proportion.forceps}{a numeric vector} }
#' @references \url{http://en.wikipedia.org/wiki/Forficula_auricularia}
#' 
#' \url{http://www.arkive.org/common-european-earwig/forficula-auricularia/}
#' 
#' \url{http://eol.org/pages/473785}
#' @source Tomkins, J.L.  and G.S. Brown. 2004. Population density drives the
#' local evolution of a threshold dimorphism. \emph{Nature} 431: 1099-1103.
#' @keywords datasets
#' @examples
#' 
#' xyplot(proportion.forceps ~ density, data=Earwigs, type='h', lwd=6)
#' 
#' 
NULL





#' Eelgrass Genotypes
#' 
#' Number of shoots (\code{shoots}) surviving in each of 32 experimental plots
#' planted with 1, 3, or 6 different genotypes of eelgrass
#' (\code{treatment.genotypes}).
#' 
#' 
#' @name Eelgrass
#' @docType data
#' @format A data frame with 32 observations on the following 2 variables.
#' \describe{ \item{genotypes}{a numeric vector of the number of
#' genotypes planted in each plot} \item{shoots}{a numeric vector of
#' the total number of shoots in each plot} }
#' @references \url{http://www.pnas.org/content/102/8/2826.abstract}
#' @source \emph{inferred from} Reusch, T.B.H., A. Ehlers, A. Hämmerli, and B.
#' Worm. 2005. Ecosystem recovery after climatic extremes enhanced by genotypic
#' diversity. \emph{Proceedings of the National Academy of Sciences (USA)} 102:
#' 2826-2831.
#' @keywords datasets
#' @examples
#' 
#' Eelgrass
#' 
#' # Convert treatment.genotypes to a factor
#' Eelgrass$genotypesF <-
#'   factor(Eelgrass$genotypes)
#' str(Eelgrass)
#' xyplot(shoots ~ genotypes, Eelgrass)
#' xyplot(shoots ~ genotypesF, Eelgrass)
#' 
NULL





#' Diet Breadth in a Rainforest Community
#' 
#' Number of different species (\code{breadth}) in 127 species
#' (\code{no.species}) in the rainforest community at El Verde, Puerto Rico
#' 
#' 
#' @name ElVerde
#' @docType data
#' @format A data frame with 38 observations on the following 2 variables.
#' \describe{ \item{breadth}{a numeric vector}
#' \item{num.species}{a numeric vector} }
#' @source Waide R.B. and W.B. Reagan, eds. 1996. \emph{The Food Web of a
#' Tropical Rainforest}. University of Chicago Press, Chicago.
#' @keywords datasets
#' @examples
#' 
#' ElVerde
#' xyplot(num.species ~ breadth, ElVerde, type='h',lwd=3)
#' 
NULL





#' Electric Fish
#' 
#' Species abundance of electric fish upstream and downstream of the entrance
#' of a tributary in the Amazon basin.
#' 
#' 
#' @name ElectricFish
#' @docType data
#' @format A data frame with 12 observations on the following 3 variables.
#' \describe{ \item{tributary}{a character vector}
#' \item{species.upstream}{a numeric vector of the number of species of
#' electric fish present upstream of the tributary}
#' \item{species.downstream}{a numeric vector of the number of species
#' of electric fish present downstream of the tributary} }
#' @references
#' \url{http://www.sciencemag.org/cgi/content/abstract/305/5692/1960}
#' @source Fernandes, C.C., J. Podos, and J.G. Lundberg. 2004. Amazonian
#' ecology: tributaries enhance the diversity of electric fishes.
#' \emph{Science} 305: 1960-1962.
#' @keywords datasets
#' @examples
#' 
#' ElectricFish
#' require(grid)
#' xyplot(species.upstream ~ species.downstream, data = ElectricFish,
#'   panel=function(x, y, ...){
#'     grid.text(ElectricFish$tributary, x=x, y=y, 
#'       rot = 45,
#'       gp = gpar(cex=.6),
#'       default.units = 'native')
#'     }
#'   )
#' 
NULL





#' Endangered and Threatened Species
#' 
#' Frequency of taxon groups on the U.S. Fish and Wildlife Service list of
#' endangered and threatened species (2002).
#' 
#' 
#' @name EndangeredSpecies
#' @docType data
#' @format A data frame with 11 observations on the following 2 variables.
#' \describe{ \item{taxon}{a character vector}
#' \item{num.species}{a numeric vector} }
#' @references \url{http://www.fws.gov/endangered/}
#' @source U.S. Fish and Wildlife Service. 2001. Number of U.S. listed species
#' per calendar year.
#' @keywords datasets
#' @examples
#' 
#' str(EndangeredSpecies)
#' EndangeredSpecies
#' 
NULL





#' 2D:4D Finger Ratio
#' 
#' The ratio of the lengths of the index finger to the ring finger in 46 males
#' and the number of CAG repeats for each.
#' 
#' 
#' @name FingerRatio
#' @docType data
#' @format A data frame with 46 observations on the following 2 variables.
#' \describe{ \item{CAGrepeats}{a numeric vector of the number of CAG
#' repeats} \item{finger.ratio}{a numeric vector of the ratio of digit
#' 2 to digit 4} }
#' @references \url{http://en.wikipedia.org/wiki/Digit_ratio}
#' @source \emph{inferred from} Manning, J.T., P.E. Bundred, D.J. Newton, and
#' B.F. Flanagan. 2003. The second to fourth digit ratio and variation in the
#' androgen receptor gene. \emph{Evolution and Human Behavior} 24: 399-405.
#' @keywords datasets
#' @examples
#' 
#' str(FingerRatio)
#' xyplot(finger.ratio ~ CAGrepeats, FingerRatio,
#'   xlab = "Number of CAG Repeats",
#'   ylab = "2D:4D Ratio"
#' )
#' 
NULL





#' Spermatophore Mass in Fireflies
#' 
#' Measurements of spermatophore mass (milligrams) in 35 fireflies
#' (\emph{Photinus ignitus}).
#' 
#' 
#' @name Fireflies
#' @docType data
#' @format A data frame with 35 observations of one variable. \describe{
#' \item{sp.mass}{} }
#' @references
#' \url{http://beheco.oxfordjournals.org/cgi/content/abstract/14/1/135}
#' 
#' \url{http://en.wikipedia.org/wiki/Firefly}
#' @source \emph{inferred from} Cratsley, C.K. and S.M. Lewis. 2003. Female
#' preference for male courtship flashes in \emph{Photinus ignitus} fireflies.
#' \emph{Behavioral Ecology} 14: 135-140.
#' @keywords datasets
#' @examples
#' 
#' str(Fireflies)
#' histogram(~sp.mass, Fireflies, n=12)
#' 
NULL





#' Firefly Flash Duration
#' 
#' Flash duration (measured in milliseconds) of a sample of male fireflies
#' (\emph{Photinus ignitus}; \emph{n} = 35).
#' 
#' 
#' @name FireflyFlash
#' @docType data
#' @format A data frame with 35 observations of one variable. \describe{
#' \item{flash}{duration of flash (milliseconds)} }
#' @source \emph{inferred from} Cratsley, C.K. and S.M. Lewis. 2003. Female
#' preference for male courtship flashes in \emph{Photinus ignitus} fireflies.
#' \emph{Behavioral Ecology} 14: 135-140.
#' @keywords datasets
#' @examples
#' 
#' str(FireflyFlash)
#' histogram(~flash, FireflyFlash)
#' 
NULL





#' Testes Size in Flies
#' 
#' Testes size (square mm; \code{Testes.area}) in 8 populations of common
#' yellow dung flies (\emph{Scathophaga stercoraria}) with different mating
#' systems (\code{Mating.system}).
#' 
#' 
#' @name FlyTestes
#' @docType data
#' @format A data frame with 8 observations on the following 2 variables.
#' \describe{ \item{mating}{a factor with levels \code{Monogamous}
#' \code{Polyandrous}} \item{testes.area}{a numeric vector} }
#' @references \url{http://en.wikipedia.org/wiki/Scathophaga_stercoraria}
#' @source Hosken, D.J. and P.I. Ward. 2001.  Experimental evidence for testis
#' size evolution via sperm competition. \emph{Ecology Letters} 4: 10-13.
#' @keywords datasets
#' @examples
#' 
#' str(FlyTestes)
#' FlyTestes
#' 
NULL





#' Forehead Patch Size in Collared Flycatachers
#' 
#' Forehead patch size in 30 male Collared Flycatachers measured in two
#' consecutive years.
#' 
#' 
#' @name FlycatcherPatch
#' @docType data
#' @format A data frame with 30 observations on the following 2 variables.
#' \describe{ \item{patch98}{a numeric vector} \item{patch99}{a
#' numeric vector} }
#' @source Griffith, S.C. and B.C. Sheldon. 2001. Phenotypic plasticity in the
#' expression of a sexually selected trait: neglected components of variation.
#' \emph{Animal Behaviour} 61: 987-993.
#' @keywords datasets
#' @examples
#' 
#' str(FlycatcherPatch)
#' xyplot(patch99 ~ patch98, FlycatcherPatch)
#' 
NULL





#' Gene Regulation in Saccharomyces
#' 
#' Number of genes regulated by 109 regulatory genes of \emph{Saccharomyces
#' cerevisiae}.
#' 
#' 
#' @name GeneRegulation
#' @docType data
#' @format A data frame with 26 observations on the following 2 variables.
#' \describe{ \item{genes.regulated}{a numeric vector}
#' \item{count}{a numeric vector} }
#' @source Guelzim, N., S. Bottani, P. Bourgine and F. Képès. 2002. Topological
#' and causal structure of the yeast transcriptional regulatory network.
#' \emph{Nature Genetics} 31: 60-63.
#' @keywords datasets
#' @examples
#' 
#' str(GeneRegulation)
#' xyplot(count ~ genes.regulated, GeneRegulation, type='h', lwd=3)
#' 
NULL





#' GlidingSnakes
#' 
#' Undulation rate (\emph{Hz}) of 8 paradise tree snakes (\emph{Chrysopelea
#' paradisi}).
#' 
#' 
#' @name GlidingSnakes
#' @docType data
#' @format A data frame with eight observations of one variable. \describe{
#' \item{undulation.rate}{undulation rate} }
#' @references
#' \url{http://www.nature.com/nature/journal/v418/n6898/abs/418603a.html}
#' 
#' \url{http://www.flyingsnake.org/}
#' @source Socha, J.J. 2002. Gliding flight in the paradise tree snake.
#' \emph{Nature} 418: 603-604.
#' @keywords datasets
#' @examples
#' 
#' histogram(~undulation.rate , data=GlidingSnakes, n=7,
#'   xlab = "Undulation rate (Hz)",
#'   type='count')
#' 
#' 
NULL





#' Godwit Arrival Dates
#' 
#' Arrival dates for males and females in 10 pairs of Black-tailed Godwits
#' (\emph{Limosa limosa})
#' 
#' 
#' @name GodwitArrival
#' @docType data
#' @format A data frame with 10 observations on the following 2 variables.
#' \describe{ \item{female}{a numeric vector} \item{male}{a
#' numeric vector} }
#' @references \url{http://en.wikipedia.org/wiki/Black-tailed_godwit}
#' @source Gunnarsson, T.G., J.A. Gill, T. Sigurbjörnsson, and W.J. Sutherland.
#' 2004. Pair bonds: arrival synchrony in migratory birds. \emph{Nature} 431:
#' 646.
#' @keywords datasets
#' @examples
#' 
#' xyplot(male~female, GodwitArrival, main='Arrival of Godwit pairs')
#' 
NULL





#' Grassland Diversity
#' 
#' Species diversity in 10 experimental plots in the Park Grass Experiment at
#' Rothamsted Experimental Station to which varying numbers of nutrients have
#' been added.
#' 
#' 
#' @name Grassland
#' @docType data
#' @format A data frame with 10 observations on the following 2 variables.
#' \describe{ \item{nutrients}{a numeric vector}
#' \item{num.species}{a numeric vector} }
#' @references \url{http://www.rothamsted.ac.uk/}
#' @source Harpole, W. S. and D. Tilman. 2007. Grassland species loss due to
#' reduced niche dimension. \emph{Nature} 446: 791-793.
#' @keywords datasets
#' @examples
#' 
#' xyplot(num.species ~ jitter(nutrients, amount=0.1), Grassland, pch=16)
#' 
NULL





#' Malaria in Populations of Great Tit
#' 
#' Two-by-two contingency table of malaria (\emph{Plasmodium}) infection status
#' in control and egg-removal populations of Great Tit (\emph{Parus major}).
#' 
#' 
#' @name GreatTitMalaria
#' @docType data
#' @format A data frame with 65 observations on the following 2 variables.
#' \describe{ \item{treatment}{a factor with levels \code{Control} and
#' \code{Egg removal}} \item{response}{a ordered factor with levels
#' \code{No Malaria} and \code{Malaria}} }
#' @references
#' \url{http://www.nature.com/nature/journal/v381/n6583/abs/381565a0.html}
#' @source Oppliger, A., P. Christe, and H. Richner. 1996. Clutch size and
#' malaria resistance. \emph{Nature} 381: 565.
#' @keywords datasets
#' @examples
#' 
#' str(GreatTitMalaria)
#' 
#' table(GreatTitMalaria)
#' 
#' if(require(vcd)) {
#'   mosaic(~treatment + response, GreatTitMalaria)
#' }
#' 
#' 
NULL





#' Diversity in Urban Green Space
#' 
#' Measures of biodiversity in 15 urban green spaces in Sheffield, England.
#' 
#' 
#' @name Greenspace
#' @docType data
#' @format A data frame with 15 observations on the following 6 variables.
#' \describe{ \item{site}{a factor with levels \code{A} - \code{O}}
#' \item{attachment}{a numeric vector} \item{area}{a numeric
#' vector} \item{butterfly}{a numeric vector} \item{bird}{a
#' numeric vector} \item{ln.plant}{a numeric vector} }
#' @references
#' \url{http://rsbl.royalsocietypublishing.org/content/3/4/390.abstract}
#' @source Fuller, R.A., K.N. Irvine, P. Devine-Wright, P.H. Warren, and K.J.
#' Gaston. 2007. Psychological benefits of greenspace increase with
#' biodiversity. \emph{Biology Letters} 3: 390-394.
#' @keywords datasets
#' @examples
#' 
#' str(Greenspace)
#' splom(Greenspace[,2:6])
#' 
NULL





#' Ornamentation and Attractiveness in Guppies
#' 
#' The father's ornamentation (composite score of color and brightness) and
#' son's attractiveness (relative rates of visits by females) in male guppies
#' (\emph{Poecilia reticulata}).
#' 
#' 
#' @name Guppies
#' @docType data
#' @format A data frame with 36 observations on the following 2 variables.
#' \describe{ \item{father.ornament}{a numeric vector}
#' \item{son.attract}{a numeric vector} }
#' @references
#' \url{http://www.nature.com/nature/journal/v406/n6791/abs/406067a0.html}
#' @source \emph{inferred from} Brooks, R. 2000. Negative genetic correlation
#' between male sexual attractiveness and survival. \emph{Nature} 406: 67-70.
#' @keywords datasets
#' @examples
#' 
#' str(Guppies)
#' xyplot(son.attract ~ father.ornament,
#'   Guppies,
#'   xlab = "Father's ornamentation",
#'   ylab = "Son's attractiveness"
#'   )
#' 
NULL





#' Hemoglobin Levels in High Altitude Populations
#' 
#' Relative rates of hemoglobin concentration in four populations of humans
#' living at different altitudes.
#' 
#' 
#' @name Hemoglobin
#' @docType data
#' @format A data frame with 40 observations on the following 3 variables.
#' \describe{ \item{hemoglobin}{a numeric vector}
#' \item{group}{a factor with levels: \code{Andes}, \code{Ethiopia},
#' \code{Tibet}, and \code{USA}} \item{relative.frequency}{a numeric
#' vector} }
#' @references \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC139295/}
#' @source \emph{inferred from} Beall, C.M., M.J. Decker, G.M. Bittenham, I.
#' Kushner, A. Gebremedhin, K.P. Strohl. 2002. An Ethiopian pattern of human
#' adaptation to high-altitude hypoxia. \emph{Proceeding of the National
#' Academy of Sciences (USA)} 99(26): 17215-17218.
#' @keywords datasets
#' @examples
#' 
#' str(Hemoglobin)
#' 
#' xyplot(relative.frequency ~ hemoglobin | group, Hemoglobin,
#'   type ='h', lwd=4, layout=c(1,4))
#' 
NULL





#' Memory and the Hippocampus
#' 
#' Spatial memory score (\code{memory}) and percent lesion of the hippocampus
#' (\code{lesion}).
#' 
#' 
#' @name HippocampusLesions
#' @docType data
#' @format A data frame with 57 observations on the following 2 variables.
#' \describe{ \item{lesion}{a numeric vector} \item{memory}{a
#' numeric vector} }
#' @source Broadbent, N.J., L.R. Squire, and R.E. Clark. 2004. Spatial memory,
#' recognition memory, and the hippocampus. \emph{Proceedings of the National
#' Academy of Sciences (USA)} 101: 14515-14520.
#' @keywords datasets
#' @examples
#' 
#' HippocampusLesions
#' 
#' xyplot(memory ~ lesion, data = HippocampusLesions,
#'   pch = 16, col = "red")
#' 
#' plot(memory ~ lesion, data = HippocampusLesions,
#'   pch = 16, col = "red")
#' 
NULL





#' Horn Length and Predation Status of Horned Lizards
#' 
#' Squamosal horn length (mm; \code{horn.length}) and predation status
#' (\code{group}; living or killed) for 184 horned lizards (\emph{Phrynosoma
#' mcalli}).
#' 
#' 
#' @name HornedLizards
#' @docType data
#' @format A data frame with 184 observations on the following 2 variables.
#' \describe{ \item{horn.length}{a numeric vector}
#' \item{group}{a numeric vector} }
#' @references \url{http://www.sciencemag.org/cgi/pdf_extract/304/5667/65}
#' @source Young, K.V., E.D. Brodie, Jr., and E.D. Brodie, III. 2004. How the
#' horned lizard got its horns. \emph{Science} 304: 65.
#' @keywords datasets
#' @examples
#' 
#' str(HornedLizards)
#' 
#' histogram(~horn.length | group, HornedLizards, 
#'   layout=c(1,2),
#'   xlab="Horn Length (mm)")
#' 
NULL





#' Human Body Temperature
#' 
#' Body temperature for 25 randomly chosen health people
#' 
#' 
#' @name HumanBodyTemp
#' @docType data
#' @format A data frame with 25 observations of one variable. \describe{
#' \item{temp}{body temperature (degrees F)} }
#' @references
#' \url{http://www.amstat.org/publications/jse/v4n2/datasets.shoemaker.html}
#' 
#' Mackowiak, P.A., Wasserman, S.S., and Levine, M.M. 1992. A critical
#' appraisal of 98.6 degrees F, the upper limit of the normal body temperature,
#' and other legacies of Carl Reinhold August Wunderlich. \emph{Journal of the
#' American Medical Association} 268: 1578-1580.
#' @source Shoemaker, A. L. 1996. What's normal? -- Temperature, gender, and
#' heart rate. \emph{Journal of Statistics Education} 4(2).
#' @keywords datasets
#' @examples
#' 
#' histogram(~temp, HumanBodyTemp)
#' stem(HumanBodyTemp$temp, scale = 2)
#' favstats(HumanBodyTemp$temp)
#' 
NULL





#' Human Gene Lengths
#' 
#' Lengths in number of nucleotides (\code{gene.length}) for 20,290 human genes
#' 
#' 
#' @name HumanGeneLengths
#' @docType data
#' @format A data frame with 20,290 observations on the following variable.
#' \describe{ \item{gene.length}{a numeric vector} }
#' @references \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC540092/}
#' 
#' \url{http://www.ensembl.org/}
#' @source Hubbard, T., D. Andrews, M. Caccamo, G. Cameron, Y. Chen, M. Clamp,
#' L. Clarke, G. Coates, T. Cox, F. Cunningham, V. Curwen, T. Cutts, T. Down,
#' R. Durbin, X. M. Fernandez-Suarez, J. Gilbert, M. Hammond, J. Herrero, H.
#' Hotz, K. Howe, V. Iyer, K. Jekosch, A. Kahari, A. Kasprzyk, D. Keefe, S.
#' Keenan, F. Kokocinsci, D. London, I. Longden, G. McVicker, C. Melsopp, P.
#' Meidl, S. Potter, G. Proctor, M. Rae, D. Rios, M. Schuster, S. Searle, J.
#' Severin, G. Slater, D. Smedley, J. Smith, W. Spooner, A. Stabenau, J.
#' Stalker, R. Storey, S. Trevanion, A. Ureta-Vidal, J. Vogel, S. White, C.
#' Woodwark, and E. Birne. 2005. Ensembl 2005. \emph{Nucleic Acids Research}
#' 33: D447-D453.
#' @keywords datasets
#' @examples
#' 
#' str(HumanGeneLengths)
#' histogram(~ gene.length, HumanGeneLengths,
#'           subset = gene.length < 15000)
#' 
NULL





#' Intense Hurricanes
#' 
#' Number of hurricanes greater than or equal to category 3 on the
#' Safir-Simpson scale during the 20th century.
#' 
#' 
#' @name Hurricanes
#' @docType data
#' @format A data frame with 4 observations on the following 2 variables.
#' \describe{ \item{hurricanes}{a numeric vector}
#' \item{count}{a numeric vector} }
#' @references
#' \url{http://www.aoml.noaa.gov/hrd/Landsea/Blakeetal_noaamemoApr2007.pdf}
#' @source Blake, E.S., E.N. Rappaport, J.D. Jarrell, and C.W. Landsea. 2005.
#' The deadliest, costliest, and most intense United States tropical cyclones
#' from 1851 to 2006 (and other frequently requested hurricane facts).
#' \emph{NOAA Technical Memorandum NWS TPC-4}.
#' @keywords datasets
#' @examples
#' 
#' Hurricanes
#' 
NULL





#' Iguana Body Length Changes
#' 
#' Body size change in 64 Galápagos marine iguanas (\emph{Amblyrhynchus
#' cristatus}) that survived the 1992-1993 El Niño event.
#' 
#' 
#' @name Iguanas
#' @docType data
#' @format A data frame with 64 observations of one variable. \describe{
#' \item{change.in.length}{} }
#' @references \url{http://en.wikipedia.org/wiki/Marine_iguana}
#' @source Wikelski, M. and C. Thom. 2000. Marine iguanas shrink to survive El
#' Niño. \emph{Nature} 403: 37-38.
#' @keywords datasets
#' @examples
#' 
#' str(Iguanas)
#' histogram(~ change.in.length, Iguanas, n = 10)
#' 
NULL





#' Intertidal Algae
#' 
#' Area coverage of red algae (\emph{Mazzaella parksii}) in two herbivore
#' treatments (\code{herbivores}) at two tide levels (\code{height}).
#' 
#' 
#' @name IntertidalAlgae
#' @docType data
#' @format A data frame with 64 observations on the following 3 variables.
#' \describe{ \item{height}{a factor with levels \code{low} and
#' \code{mid}} \item{herbivores}{a factor with levels \code{minus} and
#' \code{plus}} \item{sqrt.area}{a numeric vector} }
#' @source Harley, C.D.G. 2003. Individualistic vertical responses of
#' interacting species determine range limits across a horizontal gradient.
#' \emph{Ecology} 84: 1477-1488.
#' @keywords datasets
#' @examples
#' 
#' str(IntertidalAlgae)
#' 
#' # Using * includes the main effects and the interaction
#' aov.fit <- aov(sqrt.area ~ herbivores * height, data = IntertidalAlgae)
#' summary(aov.fit)
#' lm.fit <- lm(sqrt.area ~ herbivores * height, data = IntertidalAlgae)
#' anova(lm.fit)
#' 
NULL





#' Circadian Rhythm Phase Shift
#' 
#' Shift in circadian rhythm (hours; \code{shift}) in three light treatments
#' (\code{treatment}).
#' 
#' 
#' @name JetLagKnees
#' @docType data
#' @format A data frame with 22 observations on the following 2 variables.
#' \describe{ \item{treatment}{a factor with levels \code{control},
#' \code{eyes}, and \code{knee}} \item{shift}{a numeric vector} }
#' @references \url{http://www.sciencemag.org/cgi/content/full/297/5581/571}
#' @source \emph{inferred from} Wright, K.P., Jr. and C.A. Czeisler 2002.
#' Absence of circadian phase resetting in response to bright light behind the
#' knees. \emph{Science} 297: 571.
#' @keywords datasets
#' @examples
#' 
#' demo(sec15.1)
#' 
NULL





#' Body Mass and Beak Length in Three Species of Finches in Kenya
#' 
#' Data on body mass and beak length in three species of finches:
#' Crimson-rumped waxbill (\code{CRU.WAXB}), Cutthroat finch (\code{CUTTHROA}),
#' and White-browed sparrow weaver (\code{WB.SPARW}).
#' 
#' 
#' @name KenyaFinches
#' @docType data
#' @format A data frame with 45 observations on the following 3 variables.
#' \describe{ \item{species}{a factor with levels: \code{CRU.WAXB},
#' \code{CUTTHROA}, and \code{WB.SPARW}} \item{mass}{mass (g)}
#' \item{beak.length}{beak length (mm)} }
#' @source Schluter, D. 1988. The evolution of finch communities on islands and
#' continents: Kenya vs. Galapagos. \emph{Ecological Monographs} 58: 229-249.
#' @keywords datasets
#' @examples
#' 
#' table(KenyaFinches$species)
#' xyplot(beak.length ~ species, KenyaFinches)
#' bwplot(beak.length ~ species, KenyaFinches)
#' 
NULL





#' Brain Structure in Bilingual Humans
#' 
#' Proficiency score (summary of reading, writing, and speech) in subjects'
#' second language and density of gray matter in the left inferior parietal
#' region.
#' 
#' 
#' @name LanguageBrains
#' @docType data
#' @format A data frame with 22 observations on the following 2 variables.
#' \describe{ \item{proficiency}{a numeric vector}
#' \item{greymatter}{a numeric vector} }
#' @source Mechelli, A., J.T. Crinion, U. Noppeney, J. O'Doherty, J. Ashburner,
#' R.S. Frackowiak, and C.J. Price. 2004. Structural plasticity in the
#' bilingual brain. \emph{Nature} 431: 757.
#' @keywords datasets
#' @examples
#' 
#' str(LanguageBrains)
#' xyplot(proficiency ~ greymatter, LanguageBrains)
#' 
NULL





#' Exploited Larval Fish
#' 
#' Age (\code{age}) and coefficient of variation (\code{cv}) in larval fish
#' from exploited and unexploited species (\code{exploited}).
#' 
#' 
#' @name LarvalFish
#' @docType data
#' @format A data frame with 28 observations on the following 3 variables.
#' \describe{ \item{age}{a numeric vector} \item{cv}{a numeric
#' vector} \item{exploited}{a factor with levels \code{no} and
#' \code{yes}} }
#' @source Hsieh, C.H., C.S. Reiss, J.R. Hunter, J.R. Beddington, R.M. May, and
#' G. Sugihara. 2006. Fishing elevates variability in the abundance of
#' exploited species. \emph{Nature} 443: 859-862.
#' @keywords datasets
#' @examples
#' 
#' str(LarvalFish)
#' xyplot(cv ~ age | exploited, LarvalFish)
#' xyplot(cv ~ age, groups=exploited, LarvalFish)
#' 
NULL





#' Left-handedness and Rates of Violence
#' 
#' Prevalence of left-handedness (\code{percent.left}) and homicide rates
#' (\code{murder}) for 8 societies.
#' 
#' 
#' @name Lefthanded
#' @docType data
#' @format A data frame with 8 observations on the following 2 variables.
#' \describe{ \item{percent.left}{a numeric vector}
#' \item{murder.rate}{a numeric vector} }
#' @references
#' \url{http://rspb.royalsocietypublishing.org/content/272/1558/25.abstract}
#' @source Faurie, C. and M. Raymond. 2005. Handedness, homicide and negative
#' frequency-dependent selection. \emph{Proceedings of the Royal Society of
#' London B} 272: 25-28.
#' @keywords datasets
#' @examples
#' 
#' str(Lefthanded)
#' xyplot(murder.rate ~ percent.left, Lefthanded)
#' 
#' 
NULL





#' Time to Reproduction in Female Lions
#' 
#' Time to reproduction (\code{Days}) based on whether death of previous cubs
#' was due to infanticide (\code{New}) or accidental (\code{Same}).
#' 
#' 
#' @name LionCubs
#' @docType data
#' @format A data frame with 14 observations on the following 2 variables.
#' \describe{ \item{cause.of.death}{a factor with \code{accident} and
#' \code{infanticide}} \item{days.to.next.cub}{a numeric vector} }
#' @source Packer, C. and A.E. Pusey. 1983. Adaptations of female lions to
#' infanticide by incoming males. \emph{The American Naturalist} 121: 716-728.
#' @keywords datasets
#' @examples
#' 
#' xyplot(days.to.next.cub ~ cause.of.death, LionCubs)
#' 
NULL





#' Lion Age and Nose Coloration
#' 
#' Ages (in years; \code{age}) of 32 male lions and relative coloration of
#' their noses (\code{proportion.black}).
#' 
#' 
#' @name LionNoses
#' @docType data
#' @format A data frame with 32 observations on the following 2 variables.
#' \describe{ \item{age}{a numeric vector}
#' \item{proportion.black}{a numeric vector} }
#' @references
#' \url{http://www.nature.com/nature/journal/v428/n6979/abs/nature02395.html}
#' @source Whitman, K., A.M. Starfield, H.S. Quadling and C. Packer. 2004.
#' Sustainable trophy hunting of African lions. \emph{Nature} 428: 175-178.
#' @keywords datasets
#' @examples
#' 
#' xyplot(age ~ proportion.black, LionNoses)
#' 
NULL





#' Liver Preparation
#' 
#' The unbound fraction of taurocholate for each of five concentrations of
#' administered taurocholate.
#' 
#' 
#' @name LiverPreparation
#' @docType data
#' @format A data frame with 5 observations on the following 2 variables.
#' \describe{ \item{concentration}{a numeric vector}
#' \item{unbound.fraction}{a numeric vector} }
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/3199317}
#' @source Smallwood, R.H., D.J. Morgan, G.W. Mihaly, and R.A. Smallwood. 1998.
#' Effect of a protein binding change on unbound and total plasma
#' concentrations for drugs of intermediate hepatic extraction. \emph{Journal
#' of Pharmacokinetics and Pharmacodynamics} 16:397-411.
#' @keywords datasets
#' @examples
#' 
#' str(LiverPreparation)
#' xyplot(unbound.fraction ~ concentration, LiverPreparation)
#' 
#' 
NULL





#' Bite Force in Collard Lizards
#' 
#' Bite force (N) and territory area in 11 male collared lizards
#' (\emph{Crotaphytus collaris}).
#' 
#' 
#' @name LizardBite
#' @docType data
#' @format A data frame with 11 observations on the following 2 variables.
#' \describe{ \item{bite}{force of bite (N)}
#' \item{territory}{area of territory} }
#' @source Lappin, A. K., and J. F. Husak. 2005. Weapon performance, not size,
#' determines mating success and potential reproductive output in the collared
#' lizard (\emph{Crotaphytus collaris}). \emph{The American Naturalist} 166:
#' 426-436.
#' @note In the original publication (Lappin and Husak, 
#' 2005; Figure 3A), the data are presented in log-10 units. The 
#' data in \code{LizardBite} and in chapter 17, question 9 was
#' back-transformed using \emph{e} (i.e., \code{exp()}). To recover the data
#' from the original publication, use \code{10^(log(LizardBite$territory))}
#' and \code{10^(log(LizardBite$bite))}.
#' @keywords datasets
#' @examples
#' 
#' str(LizardBite)
#' xyplot(territory ~ bite, LizardBite)
#' 
NULL





#' Sprint Speeds in Canyon Lizards
#' 
#' Sprint speeds (\code{speed}) in 34 canyon lizards 
#' (\emph{Sceloporous merriami}) measured in
#' successive years in Big Bend National Park. Note
#' that \code{lizard} is not coded as a factor.
#' 
#' 
#' @name LizardSprint
#' @docType data
#' @format A data frame with 68 observations on the following 2 variables.
#' \describe{ \item{lizard}{a numeric vector} \item{speed}{a
#' numeric vector} }
#' @references \url{http://en.wikipedia.org/wiki/Sceloporus_merriami}
#' @source \emph{inferred from} Huey, R.B. and A.E. Dunham. 1987. The
#' repeatability of locomotor performance in natural populations of the lizard
#' \emph{Sceloporus merriami}. \emph{Evolution} 42: 1116-1120.
#' @keywords datasets
#' @examples
#' 
#' histogram(~ speed, LizardSprint)
#' Lizard2 <- aggregate(speed ~ lizard, LizardSprint, mean)
#' histogram(~ speed, Lizard2)
#' 
NULL





#' Lobster Orientation
#' 
#' Orientation of 15 lobsters relative to initial position.
#' 
#' 
#' @name Lobsters
#' @docType data
#' @format A data frame with 15 observations of one variable. \describe{
#' \item{orientation}{} }
#' @references
#' \url{http://www.unc.edu/depts/geomag/PDFGeomag/BolesandLohmann2003.pdf}
#' @source Boles, L.C. and K.J. Lohmann. 2003. True navigation and magnetic
#' maps in spiny lobsters. \emph{Nature} 421: 60-63.
#' @keywords datasets
#' @examples
#' 
#' histogram(~ orientation, Lobsters)
#' dotplot(~ orientation, Lobsters)
#' 
NULL





#' Lodgepole Pine Cone Masses
#' 
#' Masses of cones of lodgepole pines (\code{conemass}) from 16 different
#' habitat types (\code{habitat}) in western North America.
#' 
#' 
#' @name LodgepolePines
#' @docType data
#' @format A data frame with 16 observations on the following 4 variables.
#' \describe{ \item{habitat}{a factor with levels: \code{island
#' absent}, \code{island present}, and \code{mainland present}}
#' \item{conemass}{mass of cone} \item{location}{\code{island}
#' or \code{mainland}} \item{squirrels}{\code{absent} or
#' \code{present}} }
#' @references \url{http://en.wikipedia.org/wiki/Lodgepole_pine}
#' 
#' \url{http://en.wikipedia.org/wiki/Red_crossbill}
#' @source Edelaar, P. and C.W. Benkman. 2006. Replicated population divergence
#' caused by localised coevolution? A test of three hypotheses in the Red
#' Crossbill-lodgepole pine system. \emph{Journal of Evolutionary Biology} 19:
#' 1651-1659.
#' @keywords datasets
#' @examples
#' 
#' LodgepolePines
#' str(LodgepolePines)
#' xyplot(conemass ~ habitat, LodgepolePines)
#' 
NULL





#' Autoimmune Reactivity in Lupus-prone Mice
#' 
#' Autoimmune reactivity (\code{dilution} at which reactivity could be
#' detected) in three \code{treatment}s of lupus-prone mice.
#' 
#' 
#' @name LupusMice
#' @docType data
#' @format A data frame with 20 observations on the following 2 variables.
#' \describe{ \item{treatment}{a factor with levels: \code{enhanced},
#' \code{sham}, and \code{untreated}} \item{dilution}{a numeric vector
#' of the dilution level at which reactivity could be detected} }
#' @source McGaha, T.L., B. Sorrentino, and J.V. Ravetch. 2005. Restoration of
#' tolerance in lupus by targeted inhibitory receptor expression.
#' \emph{Science} 307: 590-593.
#' @keywords datasets
#' @examples
#' 
#' str(LupusMice)
#' 
NULL





#' Population Cycles of Lynx in Canada 1752-1819
#' 
#' Number of lynx pelts (\code{pelts}) reported in Canada per year from 1752 to
#' 1819.
#' 
#' 
#' @name Lynx
#' @docType data
#' @format A data frame with 68 observations on the following 2 variables.
#' \describe{ \item{year}{a numeric vector} \item{pelts}{a
#' numeric vector} }
#' @source Elton, C. and M. Nicholson. 1942. The ten-year cycle in numbers of
#' the lynx in Canada. \emph{Journal of Animal Ecology} 11: 215-244.
#' @keywords datasets
#' @examples
#' 
#' xyplot(pelts ~ year, Lynx, type=c('p','l'))
#' 
NULL





#' Marine Reserve Biomass
#' 
#' Relative biomass in 32 marine reserves.
#' 
#' 
#' @name MarineReserve
#' @docType data
#' @format A data frame with 32 observations of one variable.  \describe{
#' \item{biomass.ratio}{} }
#' @source Halpern, B.S. 2003. The impact of marine reserves: do reserves work
#' and does reserve size matter? \emph{Ecological Applications} 13: S117-S137.
#' @keywords datasets
#' @examples
#' 
#' str(MarineReserve)
#' histogram(~ biomass.ratio, MarineReserve)
#' 
NULL





#' Mass Extinction Frequency
#' 
#' The frequency of mass extinctions in the fossil record.
#' 
#' 
#' @name MassExtinctions
#' @docType data
#' @format A data frame with 21 observations on the following 2 variables.
#' \describe{ \item{num.extinctions}{a numeric vector}
#' \item{count}{a numeric vector} }
#' @references
#' \url{http://www.sciencemag.org/cgi/content/abstract/sci;215/4539/1501}
#' @source Raup, D.M. and J.J. Sepkoski, Jr. 1982. Mass extinctions in the
#' marine fossil record. \emph{Science} 215: 1501-1503.
#' @keywords datasets
#' @examples
#' 
#' MassExtinctions
#' 
NULL





#' Energy Expenditure in Mole Rats
#' 
#' Energy expenditure (\code{ln.energy}) in two castes (\code{caste}) of
#' Damaraland mole rats (\emph{Cryptomys damarensis}) with body mass
#' (\code{ln.mass}) as a covariate.
#' 
#' 
#' @name MoleRats
#' @docType data
#' @format A data frame with 35 observations on the following 3 variables.
#' \describe{ \item{caste}{a factor with levels \code{lazy} and
#' \code{worker}} \item{ln.mass}{a numeric vector}
#' \item{ln.energy}{a numeric vector} }
#' @references
#' \url{http://www.nature.com/nature/journal/v440/n7085/abs/nature04578.html}
#' @source \emph{inferred from} Scantlebury, M., J.R. Speakman, M.K.
#' Oosthuizen, T.J. Roper and N.C. Bennett. 2006. Energetics reveals
#' physiologically distinct castes in a eusocial mammal. \emph{Nature} 440:
#' 795-797.
#' @keywords datasets
#' @examples
#' 
#' MoleRats
#' 
NULL





#' Body Size in Anopheles Mosquitoes
#' 
#' Weights of female and male mosquitos (\emph{Anopheles darlingi})
#' 
#' 
#' @name Mosquitoes
#' @docType data
#' @format A data frame with 20 observations on the following 2 variables.
#' \describe{ \item{weight}{a numeric vector} \item{sex}{a
#' factor with levels \code{female} and \code{male}} }
#' @references \url{http://www.bioline.org.br/request?oc95154}
#' @source Lounibos, L.P., N. Nishimura, J. Conn, and R. Lourenco-de-Oliveira.
#' 1995. Life history correlates of adult size in the malaria vector
#' \emph{Anopheles darlingi}. \emph{Memórias do Instituto Oswaldo Cruz} 90:
#' 769-774.
#' @keywords datasets
#' @examples
#' 
#' xyplot(weight ~ sex, Mosquitoes)
#' 
NULL





#' Mouse Empathy
#' 
#' Percentage of time spent stretching in three treatments of mice. Both
#' \code{condition} and \code{treatment} code for the same variable.
#' 
#' 
#' @name MouseEmpathy
#' @docType data
#' @format A data frame with 42 observations on the following 3 variables.
#' \describe{ \item{treatment}{a factor with levels \code{Both
#' Writhing}, \code{Isolated}, and \code{One Writhing}}
#' \item{percent.stretching}{a numeric vector} \item{trt}{a
#' factor with levels \code{bw}, \code{isolated}, and \code{ow}} }
#' @source Langford, D.J., S.E. Crager, Z. Shehzah, S.B. Smith, S.G. Sotocinal,
#' J.S. Levenstadt, M.L. Chande, D.J. Levitin, J.S. Mogill. 2006. Social
#' modulation of pain as evidence for empathy in mice. \emph{Science} 312:
#' 1967-1970.
#' @keywords datasets
#' @examples
#' 
#' str(MouseEmpathy)
#' 
#' aov.fit <- aov(percent.stretching ~ treatment, data = MouseEmpathy)
#' summary(aov.fit)
#' lm.fit <- lm(percent.stretching ~ treatment, data = MouseEmpathy)
#' anova(lm.fit)
#' 
NULL





#' Cranial Capacity in Neanderthals and Modern Humans
#' 
#' Brain size (\code{lnbrain}) and body mass (\code{lnmass}) in Neanderthals
#' and early modern humans (\code{species}).
#' 
#' 
#' @name NeanderthalBrains
#' @docType data
#' @format A data frame with 39 observations on the following 3 variables.
#' \describe{ \item{ln.mass}{log of body mass (kg)}
#' \item{ln.brain}{log of brain size} \item{species}{a factor
#' with levels \code{neanderthal} \code{recent}} }
#' @source Ruff, C.B., E. Trinkaus, and T.W. Holliday. 1997. Body mass and
#' encephalization in Pleistocene \emph{Homo}. \emph{Nature} 387: 173-176.
#' @keywords datasets
#' @examples
#' 
#' xyplot(ln.brain ~ ln.mass, data=NeanderthalBrains, groups=species)
#' 
NULL





#' Effects of Trimethadione on Lifespan in Nematodes
#' 
#' \code{lifespan} of the nematode \emph{Caenorhabditis elegans} in control and
#' three experimental \code{treatment}s of the anticonvulsant drug
#' trimethadione.
#' 
#' 
#' @name NematodeLifespan
#' @docType data
#' @format A data frame with 200 observations on the following 2 variables.
#' \describe{ \item{treatment}{a factor with levels: \code{adult},
#' \code{larva}, \code{larva+adult}, and \code{water}}
#' \item{lifespan}{a numeric vector of lifespan} }
#' @references
#' \url{http://www.sciencemag.org/cgi/content/abstract/307/5707/258}
#' @source \emph{inferred from} Evason, K., C. Huang, I. Yamben, D.F. Covey,
#' and K. Kornfeld. 2005. Anticonvulsant medications extend worm life-span.
#' \emph{Science} 307: 258-262.
#' @keywords datasets
#' @examples
#' 
#' str(NematodeLifespan)
#' 
NULL





#' Photosynthesis in Neotropical Trees
#' 
#' Photosynthetic capacity (\code{photosynthetic.capacity}) and number of
#' fruits produced in the previous season (\code{previous.fruits}) of 9 females
#' of \emph{Ocotea tenera}.
#' 
#' 
#' @name NeotropicalTrees
#' @docType data
#' @format A data frame with 9 observations on the following 2 variables.
#' \describe{ \item{previous.fruits}{a numeric vector}
#' \item{photosynthetic.capacity}{a numeric vector} }
#' @references \url{http://www.pnas.org/content/101/21/8051.long}
#' @source \emph{inferred from} Wheelwright, N.T. and B.A. Logan. 2004.
#' Previous-year reproduction reduces photosynthetic capacity and slows
#' lifetime growth in females of a neotropical tree. \emph{Proceedings of the
#' National Academy of Sciences (USA)} 101: 8051-8055.
#' @keywords datasets
#' @examples
#' 
#' str(NeotropicalTrees)
#' NeotropicalTrees
#' 
NULL





#' Tetrodotoxin Resistance in Garter Snakes
#' 
#' Percent reduction in crawl speed (\code{resistance}) in the garter snake
#' after injection of the neurotoxin tetrodotoxin from the rough-skinned newt
#' (\emph{Taricha granulosa}).
#' 
#' 
#' @name Newts
#' @docType data
#' @format A data frame with 12 observations on the following 2 variables.
#' \describe{ \item{locality}{a factor with levels: \code{Benton} and
#' \code{Warrenton}} \item{resistance}{a numeric vector} }
#' @references
#' \url{http://www.sciencemag.org/cgi/content/abstract/297/5585/1336}
#' @source Geffeney, S., E.D. Brodie, Jr., P.C. Ruben, and E.D. Brodie III.
#' 2002. Mechanisms of adaptation in a predator-prey arms race: TTX-resistant
#' sodium channels. \emph{Science} 297: 1336-1339.
#' @keywords datasets
#' @examples
#' 
#' Newts
#' 
NULL





#' No Smoking Day
#' 
#' Number of workplace injuries on No Smoking Day (\code{Injuries.on.NSD})
#' compared to the same Wednesday in the previous year
#' (\code{Injuries.before.NSD}) for 1987-1996.
#' 
#' 
#' @name NoSmokingDay
#' @docType data
#' @format A data frame with 10 observations on the following 3 variables.
#' \describe{ \item{year}{a numeric vector}
#' \item{injuries.before.NSD}{a numeric vector}
#' \item{injuries.on.NSD}{a numeric vector} }
#' @references \url{http://www.nosmokingday.org.uk/}
#' @source Waters, A.J., M.J. Jarvis, and S.R. Sutton. 1998. Nicotine
#' withdrawal and accident rates. \emph{Nature} 394: 137.
#' @keywords datasets
#' @examples
#' 
#' NoSmokingDay
#' 
NULL





#' Atlantic Cod Recruits
#' 
#' Number (\eqn{\log_{10}}{log10} transformed) of Atlantic cod (\emph{Gadus
#' morhua}) that recruited (grew to catchable size) in the North Sea over a 39
#' years span.
#' 
#' 
#' @name NorthSeaCod
#' @docType data
#' @format A data frame with 39 observations of one variable. \describe{
#' \item{log10.recruits}{} }
#' @references
#' \url{http://www.nature.com/nature/journal/v426/n6967/abs/nature02164.html}
#' @source \emph{inferred from} Beaugrand, G., K.M. Brander, J.A. Lindley, S.
#' Souissi, and P.C. Reid. 2003. Plankton effect on cod recruitment in the
#' North Sea. \emph{Nature} 426: 661-664.
#' @keywords datasets
#' @examples
#' 
#' favstats(NorthSeaCod$log10.recruits)
#' 
NULL





#' Ostrich Body and Brain Temperatures
#' 
#' Body and brain temperatures (\eqn{^{\circ}}{degrees }C) in free-ranging
#' ostriches (\emph{Struthio camelus}) at the the Lichtenburg Game Breeding
#' Centre, Lichtenburg, South Africa.
#' 
#' 
#' @name OstrichTemp
#' @docType data
#' @format A data frame with 6 observations on the following 3 variables.
#' \describe{ \item{ostrich}{a numeric vector identifying ostrich
#' number} \item{body.temp}{a numeric vector of body temperature in
#' \eqn{^{\circ}}{degrees }C} \item{brain.temp}{a numeric vector of
#' brain temperature in \eqn{^{\circ}}{degrees }C} }
#' @references \url{http://jeb.biologists.org/cgi/content/abstract/206/7/1171}
#' 
#' \url{http://www.sa-venues.com/game-reserves/nwp_lichtenburg.htm}
#' @source Fuller, A., P.R. Kamerman, S.K. Maloney, G. Mitchell, and D.
#' Mitchell. 2003. Variability in brain and arterial blood temperatures in
#' free-ranging ostriches in their natural habitat. \emph{Journal of
#' Experimental Biology} 206: 1171-1181.
#' @keywords datasets
#' @examples
#' 
#' xyplot(brain.temp ~ body.temp, OstrichTemp)
#' 
NULL





#' Penguin Heart Rate
#' 
#' Slope of regressions of mass-specific metabolic rate on heart rate for three
#' groups of Macaroni Penguins.
#' 
#' 
#' @name Penguins
#' @docType data
#' @format A data frame with 24 observations on the following 2 variables.
#' \describe{ \item{group}{a factor with levels \code{BF}, \code{BM},
#' and \code{MF}} \item{slope}{a numeric vector} }
#' @source Green, J. A., P. J. Butler, A. J. Woakes, I. L. Boyd and R. L.
#' Holder. 2001. Heart rate and rate of oxygen consumption of exercising
#' macaroni penguins. \emph{Journal of Experimental Biology} 204: 673-684.
#' @keywords datasets
#' @examples
#' 
#' str(Penguins)
#' dotplot(slope ~ group, Penguins)
#' 
NULL





#' Population Persistence Times
#' 
#' Persistence times (\code{generations}) in the annual plant \emph{Cardamine
#' pensylvanica} in four experimental populations (\code{treatment}).
#' 
#' 
#' @name PlantPersistence
#' @docType data
#' @format A data frame with 16 observations on the following 2 variables.
#' \describe{ \item{generations}{a numeric vector}
#' \item{treatment}{a factor with levels: \code{Isolated},
#' \code{Medium}, \code{Long}, and \code{Continuous}} }
#' @source Molofsky, J. and J.-B. Ferdy. 2005. Extinction dynamics in
#' experimental metapopulations. \emph{Proceedings of the National Academy of
#' Sciences (USA)} 102: 3726-3731.
#' @keywords datasets
#' @examples
#' 
#' xyplot(generations~treatment, PlantPersistence)
#' 
NULL





#' Sterility in Hybrid Pollens
#' 
#' Genetic distance between pairs of species of the genus \emph{Silene} and
#' proportion of their hybrid offspring that are sterile.
#' 
#' 
#' @name Pollen
#' @docType data
#' @format A data frame with 23 observations on the following 2 variables.
#' \describe{ \item{genetic.distance}{a numeric vector}
#' \item{proportion.sterile}{a numeric vector} }
#' @source Moyle, L.C., M.S. Olson, and P. Tiffin. 2004. Patterns of
#' reproductive isolation in three angiosperm genera. \emph{Evolution} 58:
#' 1195-1208.
#' @keywords datasets
#' @examples
#' 
#' str(Pollen)
#' xyplot(proportion.sterile ~ genetic.distance, Pollen)
#' 
NULL





#' Powerball Tickets Sold
#' 
#' The number of Powerball tickets sold per day of the week for three years.
#' 
#' 
#' @name Powerball
#' @docType data
#' @format A data frame with 7 observations on the following 2 variables.
#' \describe{ \item{day}{a character vector}
#' \item{millions.of.tickets.sold}{a numeric vector} }
#' @references
#' \url{http://www.dartmouth.edu/~chance/chance_news/recent_news/chance_news_13.02.html}
#' @source Oster, E. 2004. Dreaming big: Why do people play Powerball?
#' \emph{Chance News} 13.02.
#' @keywords datasets
#' @examples
#' 
#' Powerball
#' xyplot(millions.of.tickets.sold ~ day, Powerball)
#' 
NULL





#' Primate Metabolic Rates
#' 
#' Body mass (g) and metabolic rate (watts) for 17 species of primates.
#' 
#' 
#' @name PrimateMetabolism
#' @docType data
#' @format A data frame with 17 observations on the following 2 variables.
#' \describe{ \item{mass}{mass (g) } \item{bmr}{metabolic rate
#' (watts)} }
#' @references \url{http://jeb.biologists.org/cgi/content/abstract/160/1/25}
#' @source Heusner, A.A. 1991. Size and power in mammals. \emph{Journal of
#' Experimental Biology} 160: 25-54.
#' @keywords datasets
#' @examples
#' 
#' str(PrimateMetabolism)
#' xyplot(bmr ~ mass, PrimateMetabolism)
#' xyplot(bmr ~ mass, PrimateMetabolism, scales=list(log=TRUE))
#' 
NULL





#' Primate White Blood Cell Counts and Promiscuity
#' 
#' White blood cell (WBC) counts in pairs of closely related primate species
#' 
#' 
#' @name PrimateWBC
#' @docType data
#' @format A data frame with 9 observations on the following 2 variables.
#' \describe{ \item{WBC.less}{a numeric vector}
#' \item{WBC.more}{a numeric vector} }
#' @source Nunn, C.L., J.L. Gittleman, and J. Antonovics. 2000. Promiscuity and
#' the primate immune system. \emph{Science} 290: 1168-1170.
#' @keywords datasets
#' @examples
#' 
#' xyplot(WBC.more ~ WBC.less, PrimateWBC)
#' 
NULL





#' Progesterone and Exercise
#' 
#' Progesterone levels and rates of ventilation during submaximal exercise in
#' 30 women.
#' 
#' 
#' @name ProgesteroneExercise
#' @docType data
#' @format A data frame with 30 observations on the following 2 variables.
#' \describe{ \item{progesterone}{a numeric vector}
#' \item{ventilation}{a numeric vector} }
#' @references \url{http://jeb.biologists.org/cgi/content/abstract/205/2/233}
#' @source Brutsaert, T.D., H. Spielvogel, E. Caceres, M. Araoz, R.T.
#' Chatterton, V.J. Vitzthum. 2002. Effect of menstrual cycle phase on exercise
#' performance of high-altitude native women at 3600 m. \emph{Journal of
#' Experimental Biology} 205: 233-239
#' @keywords datasets
#' @examples
#' 
#' str(ProgesteroneExercise)
#' xyplot(ventilation ~ progesterone, ProgesteroneExercise)
#' 
NULL





#' Multiple Mating in Pseudoscorpions
#' 
#' Successful numbers of broods (\code{Number.of.successful.broods}) in two
#' groups of female pseudoscrpions (\emph{Cordylochernes scorpioides}), one
#' mated to the same male twice and one to two different males.
#' 
#' 
#' @name Pseudoscorpions
#' @docType data
#' @format A data frame with 36 observations on the following 2 variables.
#' \describe{ \item{treatment}{a factor with levels \code{DM}
#' \code{SM}} \item{successful.broods}{a numeric vector} }
#' @references \url{http://www.pnas.org/content/96/18/10236.long}
#' @source Newcomer, S.D., J.A. Zeh, and D.W. Zeh. 1999. Genetic benefits
#' enhance the reproductive success of polyandrous females. \emph{Proceedings
#' of the National Academy of Sciences (USA)} 96: 10236-10241.
#' @keywords datasets
#' @examples
#' 
#' str(Pseudoscorpions)
#' bwplot(successful.broods ~ treatment, Pseudoscorpions)
#' aggregate(successful.broods ~ treatment, Pseudoscorpions, favstats)
#' 
NULL





#' Pufferfish Mimicry
#' 
#' Number of \code{predators} approaching models painted to resemble pufferfish
#' (\emph{Canthigaster valentini}) across a range of similarities
#' (\code{resemblance})
#' 
#' 
#' @name Pufferfish
#' @docType data
#' @format A data frame with 20 observations on the following 2 variables.
#' \describe{ \item{resemblance}{a numeric vector}
#' \item{predators}{a numeric vector} }
#' @references
#' \url{http://rspb.royalsocietypublishing.org/content/270/1516/667.full.pdf}
#' 
#' \url{http://en.wikipedia.org/wiki/Canthigaster_valentini}
#' 
#' \url{http://www.fishbase.org/Summary/SpeciesSummary.php?id=6544}
#' @source Caley, M.J. and D. Schluter. 2003. Predators favour mimicry in a
#' tropical reef fish. \emph{Proceedings of the Royal Society of London Series
#' B, Biological Sciences} 270: 667-672.
#' @keywords datasets
#' @examples
#' 
#' str(Pufferfish)
#' xyplot(predators ~ jitter(resemblance, amount = 0.1), Pufferfish)
#' Pufferfish
#' 
NULL





#' Temperature Change and Meal Size in Rattlesnakes
#' 
#' Temperature change after a meal (% of body mass) in 17 South American
#' rattlesnakes (\emph{Crotalus durissus}).
#' 
#' 
#' @name Rattlesnakes
#' @docType data
#' @format A data frame with 17 observations on the following 2 variables.
#' \describe{ \item{meal.size}{a numeric vector}
#' \item{temp.change}{a numeric vector} }
#' @references \url{http://jeb.biologists.org/cgi/content/abstract/207/4/579}
#' @source Tattersall, G.J., W.K. Milsom, A.S. Abe, S.P. Brito, and D.V.
#' Andrade. 2004. The thermogenesis of digestion in rattlesnakes. \emph{Journal
#' of Experimental Biology} 207: 579-585.
#' @keywords datasets
#' @examples
#' 
#' str(Rattlesnakes)
#' xyplot(meal.size ~ temp.change, Rattlesnakes)
#' 
NULL





#' Rigormortis and Time of Death
#' 
#' Number of bodies reaching rigormortis in each hour after death.
#' 
#' 
#' @name Rigormortis
#' @docType data
#' @format A data frame with 12 observations on the following 2 variables.
#' \describe{ \item{hours}{a numeric vector} \item{count}{a
#' numeric vector} }
#' @source Pounder, D.J. 1995. Postmortem changes and time of death. University
#' of Dundee.
#' @keywords datasets
#' @examples
#' 
#' xyplot(count ~ hours, Rigormortis, type='h', lwd=3)
#' barchart(count ~ hours, Rigormortis, horizontal=FALSE, origin=0)
#' 
NULL





#' Indian Rope Trick
#' 
#' Perceived impressiveness (\code{impressiveness}) of a written account of the
#' Indian Rope Trick and the corresponding number of \code{year}s since it was
#' witnessed.
#' 
#' 
#' @name RopeTrick
#' @docType data
#' @format A data frame with 21 observations on the following 2 variables.
#' \describe{ \item{years}{a numeric vector}
#' \item{impressiveness}{a numeric vector} }
#' @references \url{http://www.richardwiseman.com/resources/ropeJSPR.pdf}
#' @source Wiseman, R. and P. Lamont. 1996. Unravelling the Indian rope-trick.
#' \emph{Nature} 383: 212-213.
#' @keywords datasets
#' @examples
#' 
#' xyplot(impressiveness ~ years, RopeTrick)
#' 
#' 
NULL





#' Sagebrush Cricket Mating Times
#' 
#' Time to mating (\code{time.to.mating}) in fed and unfed (\code{treatment})
#' sagebrush crickets (\emph{Cyphoderris strepitans}).
#' 
#' 
#' @name SagebrushCrickets
#' @docType data
#' @format A data frame with 24 observations on the following 2 variables.
#' \describe{ \item{treatment}{a factor with levels: \code{fed} and
#' \code{starved}} \item{time.to.mating}{a numeric vector} }
#' @source Chadwick Johnson, J., T.M. Ivy, and S.K. Sakaluk. 1999. Female
#' remating propensity contingent on sexual cannibalism in sagebrush crickets,
#' \emph{Cyphoderris strepitans}: a mechanism of cryptic female choice.
#' \emph{Behavioral Ecology} 10: 227-233.
#' @keywords datasets
#' @examples
#' 
#' SagebrushCrickets
#' str(SagebrushCrickets)
#' 
NULL





#' Pacific Salmon Color
#' 
#' Skin color sockeye and kokanee morphs of the Pacific salmon
#' (\emph{Oncorhynchus nerka}) raised in a low carotenoid environment.
#' 
#' 
#' @name SalmonColor
#' @docType data
#' @format A data frame with 35 observations on the following 2 variables.
#' \describe{ \item{species}{a factor with levels \code{kokanee} and
#' \code{sockeye}} \item{skin.color}{a numeric vector} }
#' @source Craig, J.K. and C. Foote. 2001. Countergradient variation and
#' secondary sexual color: phenotypic convergence promotes genetic divergence
#' in carotenoid use between sympatric anadromous and nonanadromous morphs of
#' sockeye salmon (\emph{Oncorhynchus nerka}). \emph{Evolution} 55: 380-391.
#' @keywords datasets
#' @examples
#' 
#' SalmonColor
#' histogram(~ skin.color | species, SalmonColor)
#' bwplot(skin.color ~ species, SalmonColor)
#' 
NULL





#' Number of Seedlings Per Quadrat
#' 
#' Data on frequency of seeding per quadrat for 80 hypothetical quadrats.
#' 
#' 
#' @name Seedlings
#' @docType data
#' @format A data frame with 8 observations on the following 2 variables.
#' \describe{ \item{seedlings}{a numeric vector} \item{count}{a
#' numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' Seedlings
#' 
NULL





#' Data for Meta-analysis
#' 
#' Data for meta-analysis on the relationship between testosterone and
#' aggression.
#' 
#' 
#' @name Selection
#' @docType data
#' @format A data frame with 814 observations on the following 8 variables.
#' 
#' \describe{ \item{species}{species investigated}
#' \item{traitname}{trait investigated}
#' \item{strength.of.selection}{strength of selection}
#' \item{sample.size}{size of sample} \item{authors}{authors of
#' publication} \item{year}{year of publication}
#' \item{journal}{journal of publication}
#' \item{volume.pages}{volume and pages} }
#' @source Kingsolver, J.G., H.E. Hoekstra, J.M. Hoekstra, D. Berrigan, S.N.
#' Vignieri, C.E. Hill, A. Hoang, P. Gibert, and P. Beerli. 2001. The strength
#' of phenotypic selection in natural populations. \emph{The American
#' Naturalist} 157: 245-261.
#' @keywords datasets
#' @examples
#' 
#' histogram(~ strength.of.selection, Selection,n=40)
#' table(Selection$species) -> s
#' table(s)
#' s[s>10] # most common species
#' table(Selection$traitname) -> t
#' table(t)
#' t[t>10] # most common traits
#' 
NULL





#' Sexual Conflict
#' 
#' Number of species in each of two taxa in closely related taxon pairings and
#' the difference between the two groups. One taxon has multiple matings
#' (\code{polyandrous.species}) and one has only single matings
#' (\code{monandrous.species}).
#' 
#' 
#' @name SexualSelection
#' @docType data
#' @format A data frame with 25 observations on the following 4 variables.
#' \describe{ \item{polyandrous.species}{a numeric vector}
#' \item{monandrous.species}{a numeric vector}
#' \item{difference}{a numeric vector}
#' \item{taxon.pair}{identifier} }
#' @references \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC27046/}
#' @source Arnqvist, G., M. Edvardsson, U. Friberg, and T. Nilsson. 2000.
#' Sexual conflict promotes speciation in insects. \emph{Proceedings of the
#' National Academy of Sciences (USA)} 97: 10460-10464.
#' @keywords datasets
#' @examples
#' 
#' SexualSelection
#' 
#' histogram(~ difference, SexualSelection, n = 20)
#' 
#' hist(SexualSelection$difference, breaks = 20)
#' 
#' # Calculate the number of tests and the number of negative tests
#' (n <- length(SexualSelection$difference))
#' (n.neg <- sum(SexualSelection$difference < 0))
#' 
#' 2 * pbinom(q = n.neg, size = n, prob = 0.5)
#' 
#' # With a binomial test
#' binom.test(n.neg, n, p = 0.5)
#' 
NULL





#' Shad Parasites
#' 
#' Frequency of the nematode \emph{Camallanus oxycephalus} per fish.
#' 
#' 
#' @name ShadParasites
#' @docType data
#' @format A data frame with 7 observations on the following 2 variables.
#' \describe{ \item{parasites}{a numeric vector} \item{count}{a
#' numeric vector} }
#' @references
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/9881385}
#' @source Shaw, D.J., B.T. Grenfell, and A.P. Dobson. 1998. Patterns of
#' macroparasite aggregation in wildlife host populations. \emph{Parasitology}
#' 117: 597-610.
#' @keywords datasets
#' @examples
#' 
#' ShadParasites
#' 
NULL





#' Seal Body Lengths and Age
#' 
#' Body length (cm) and age (days) for 9,665 female Northern fur seals
#' (\emph{Callorhinus ursinus}).
#' 
#' 
#' @name ShrinkingSeals
#' @docType data
#' @format A data frame with 9,665 observations on the following 2 variables.
#' \describe{ \item{age}{age (days)} \item{length}{body length
#' (cm) } }
#' @source Trites, A.W. and M.A. Bigg. 1996. Physical growth of northern fur
#' seals: seasonal fluctuations and migratory influences. \emph{Journal of
#' Zoology (London)} 238: 459-482.
#' @keywords datasets
#' @examples
#' 
#' str(ShrinkingSeals)
#' 
#' plot(ShrinkingSeals, pch = 16, cex = 0.5)
#' xyplot(length ~ age, ShrinkingSeals, pch=16, alpha=0.65, cex=0.6)
#' 
NULL





#' Ambient Temperature and O-Ring Failures
#' 
#' Data on \code{temperature} and number of O-ring \code{failures} for 23 space
#' shuttle launches.
#' 
#' 
#' @name ShuttleDisaster
#' @docType data
#' @format A data frame with 23 observations on the following 2 variables.
#' \describe{ \item{temperature}{a numeric vector}
#' \item{failures}{a numeric vector} }
#' @references Tufte, E.R. 1997. \emph{Visual Explanations: Images and
#' Quantities, Evidence and Narrative}. Graphics Press.
#' 
#' \url{http://www.edwardtufte.com/tufte/}
#' @source Dalal, S.R., E.B. Fowlkes, and B. Hoadley. 1989. Risk analysis of
#' the Space Shuttle: Pre-Challenger prediction of failure. \emph{Journal of
#' the American Statistical Association} 408: 945-957.
#' @keywords datasets
#' @examples
#' 
#' str(ShuttleDisaster)
#' xyplot( jitter(failures, amount=0.1) ~ temperature, ShuttleDisaster,
#'   ylab='number of failures'
#'   )
#' 
NULL





#' Rate of Speciation in Silverswords
#' 
#' Speciation "waiting times" in Hawaiian silverswords (\emph{Dubautia}).
#' 
#' 
#' @name Silversword
#' @docType data
#' @format A data frame with 21 observations on the following variable.
#' \describe{ \item{waiting.time}{a numeric vector} }
#' @source \emph{inferred from} Baldwin, B. G. and M. J. Sanderson 1998. Age
#' and rate of diversification of the Hawaiian silversword alliance
#' (Compositae). \emph{Proceedings of the National Academy of Sciences (USA)}
#' 95: 9402-9406.
#' @keywords datasets
#' @examples
#' 
#' Silversword
#' 
NULL





#' Sleep and Learning
#' 
#' The increase in "slow-wave" sleep and improvements in spatial learning tasks
#' in 10 humans.
#' 
#' 
#' @name SleepAndPerformance
#' @docType data
#' @format A data frame with 10 observations on the following 2 variables.
#' \describe{ \item{sleep}{a numeric vector}
#' \item{improvement}{a numeric vector} }
#' @references \url{http://www.ncbi.nlm.nih.gov/pubmed/15184907}
#' @source Huber, R., M.F. Ghilardi, M. Massimini, and G. Tononi. 2004. Local
#' sleep and learning. \emph{Nature} 430: 78-81.
#' @keywords datasets
#' @examples
#' 
#' str(SleepAndPerformance)
#' xyplot(improvement ~ sleep, SleepAndPerformance)
#' 
NULL





#' Body Masses of Female Sockeye Salmon
#' 
#' Body Masses of 228 female Sockeye Salmon (\emph{Oncorhynchus nerka};
#' \url{http://www.nmfs.noaa.gov/pr/species/fish/sockeyesalmon.htm})
#' 
#' 
#' @name SockeyeFemales
#' @docType data
#' @format A data frame with 228 observations of a single variable. \describe{
#' \item{mass}{body mass (kg)} }
#' @source Hendry, A.P., O.K. Berg, and T.P. Quinn. 1999. Condition dependence
#' and adaptation-by-time: Breeding date, life history, and energy allocation
#' within a population of salmon. \emph{Oikos} 85: 499-514.
#' @keywords datasets
#' @examples
#' 
#' str(SockeyeFemales)
#' summary(SockeyeFemales)
#' 
NULL





#' Lifetime Reproductive Success in House Sparrows
#' 
#' A cross table of lifetime reproductive success (\code{LifetimeRS}) in female
#' and male house sparrows \emph{Passer domesticus} in Norway.
#' 
#' @name Sparrows
#' @docType data
#' @format A data frame with 9 observations on the following 3 variables.
#' \describe{ \item{lifetimeRS}{a numeric vector}
#' \item{females}{a numeric vector} \item{males}{a numeric vector} }
#' @source Jensen, H., B.-E. Saether, T.H. Ringsby, J. Tufto, S.C.
#' Griffith, and H. Ellegren. 2004. Lifetime reproductive success in relation
#' to morphology in the House Sparrow \emph{Passer domesticus}. \emph{Journal
#' of Animal Ecology} 73: 599-611.
#' @keywords datasets
#' @examples
#' 
#' Sparrows
#' 
NULL





#' Social Spiders
#' 
#' Web height above ground (cm) and colony size for 17 colonies of the spider
#' \emph{Cryptophora citricola} in Gabon.
#' 
#' 
#' @name SpiderColonies
#' @docType data
#' @format A data frame with 17 observations on the following 3 variables.
#' \describe{ \item{colony}{identifier} \item{height}{height of
#' web above ground (cm)} \item{spiders}{number of spiders in colony} }
#' @source Rypstra, A. L. 1979. Foraging folks of spiders, a study of aggregate
#' behavior in \emph{Cryptophora citricola} Forskal (Araneae: Araneidae) in
#' West Africa. \emph{Behavioral Ecology and Sociobiology} 5: 291-300.
#' @keywords datasets
#' @examples
#' 
#' str(SpiderColonies)
#' SpiderColonies
#' 
NULL





#' Spider Running Speeds after Amputation
#' 
#' Data on speed before and after amputation of a pedipalp in the spider genus
#' \emph{Tidarren}.
#' 
#' 
#' @name SpiderSpeed
#' @docType data
#' @format A data frame with 32 observations on the following 2 variables.
#' \describe{ \item{speed.before}{speed (cm/s) before amputation }
#' \item{speed.after}{speed (cm/s) after amputation } }
#' @references \url{http://en.wikipedia.org/wiki/Pedipalp}, \url{http://en.wikipedia.org/wiki/Tidarren}, \url{http://www.pnas.org/content/101/14/4883}
#'
#' @source Ramos, M., D.J. Irschick, and T.E. Christenson. 2004. Overcoming an
#' evolutionary conflict: Removal of a reproductive organ greatly increases
#' locomotor performance. \emph{Proceedings of the National Academy of Sciences
#' (USA)} 101: 4883-4887.
#' @keywords datasets
#' @examples
#' xyplot(speed.after ~ speed.before, SpiderSpeed)
#' favstats(SpiderSpeed$speed.before)
#' favstats(SpiderSpeed$speed.after)
#' favstats(SpiderSpeed$speed.after - SpiderSpeed$speed.before)
#' 
NULL





#' Eye Widths in Stalk-Eyed Flies
#' 
#' Eye width in 9 male stalk-eyed flies (\emph{Cyrtodiopsis dalmanni}).
#' 
#' 
#' @name Stalkies1
#' @docType data
#' @format a data frame with 9 observations of 1 variable \describe{
#' \item{eye.span}{eye span (mm)} }
#' @source Data provided by Kevin Fowler, University College, London.
#' @keywords datasets
#' @examples
#' 
#' Stalkies1
#' 
NULL





#' Stalk-eyed Fly Eyespan
#' 
#' Eyespan width (mm; \code{Eye.span}) in 45 stalk-eyed flies
#' (\emph{Cyrtodiopsis dalmanni}) fed a corn or cotton diet (\code{Food}).
#' 
#' 
#' @name Stalkies2
#' @docType data
#' @format A data frame with 45 observations on the following 2 variables.
#' \describe{ \item{food}{a factor with levels \code{Corn}
#' \code{Cotton}} \item{eye.span}{a numeric vector} }
#' @source David, P., T. Bjorksten, K. Fowler, and A. Pomiankowski. 2000.
#' Condition-dependent signalling of genetic variation in stalk-eyed flies.
#' \emph{Nature} 406: 186-188.
#' @keywords datasets
#' @examples
#' 
#' str(Stalkies2)
#' xyplot(eye.span ~ food, Stalkies2)
#' aggregate(eye.span ~ food, Stalkies2, FUN = favstats)
#' 
NULL





#' Number of Lateral Plates in Sticklebacks
#' 
#' Number of lateral plates (\code{plates}) in threespine sticklebacks
#' (\emph{Gasterosteus aculeatus}) with three different \emph{Ectodysplasin}
#' genotypes (\code{mm}, \code{Mm}, and \code{MM}).
#' 
#' 
#' @name SticklebackPlates
#' @docType data
#' @format A data frame with 344 observations on the following 2 variables.
#' \describe{ \item{genotype}{a factor with levels \code{mm},
#' \code{Mm}, and \code{MM}} \item{plates}{number of plates} }
#' @references Colosimo P.F., K.E. Hosemann, S. Balabhadra, G. Villarreal, M.
#' Dickson, J. Grimwood, J Schmutz, R.M. Myers, D. Schluter, D.M. Kingsley.
#' 2005. Widespread parallel evolution in sticklebacks by repeated fixation of
#' ectodysplasin alleles. \emph{Science }307: 1928-33.
#' \url{http://www.sciencemag.org/cgi/content/full/307/5717/1928}
#' @source Colosimo, P.F., C.L. Peichel, K. Nereng, B.K. Blackman, M.D.
#' Shapiro, D. Schluter, and D.M. Kingsley. 2004. The genetic architecture of
#' parallel armor plate reduction in threespine sticklebacks. \emph{PLoS
#' Biology} 2: 635-641.
#' \url{http://www.plosbiology.org/article/info:doi/10.1371/journal.pbio.0020109}
#' @keywords datasets
#' @examples
#' 
#' aggregate(plates ~ genotype, SticklebackPlates, FUN = favstats)
#' 
#' histogram( ~ plates | genotype, SticklebackPlates, 
#'   layout = c(1,3),
#'   n = 15,
#'   xlab = "Number of Lateral Body Plates"
#'   )
#' 
#' densityplot( ~ plates | genotype, SticklebackPlates, 
#'   xlab = "Number of Lateral Body Plates",
#'   layout = c(1,3)
#'   )
#' 
#' 
NULL





#' Mating Preferences in Sticklebacks
#' 
#' Mating preference in 9 populations of three-spined sticklebacks.
#' 
#' 
#' @name SticklebackPreference
#' @docType data
#' @format A data frame with 9 observations of one variable. \describe{
#' \item{preference.index}{a numeric vector} }
#' @references
#' \url{http://www.nature.com/nature/journal/v429/n6989/abs/nature02556.html}
#' @source McKinnon, J. S., S. Mori, B.K. Blackman, L. David, D.M. Kingsley, L.
#' Jamieson, J. Chou, and D. Schluter. 2004. Evidence for ecology's role in
#' speciation. \emph{Nature} 429: 294-298.
#' @keywords datasets
#' @examples
#' 
#' SticklebackPreference
#' histogram(~ preference.index, SticklebackPreference)
#' dotplot(~ preference.index, SticklebackPreference)
#' 
NULL





#' Sumo Wrestling Wins
#' 
#' Counts of number of wins for sumo wrestlers.
#' 
#' 
#' @name Sumo
#' @docType data
#' @format A data frame with 16 observations on the following 2 variables.
#' \describe{ \item{wins}{a numeric vector} \item{count}{a
#' numeric vector} }
#' @source Duggan, M. and S.D. Leavitt. 2002. Winning isn't everything:
#' Corruption in sumo wrestling. \emph{The American Economic Review} 92:
#' 1594-1605.
#' @keywords datasets
#' @examples
#' 
#' xyplot(count ~ wins, Sumo, type='h', lwd=4)
#' 
NULL





#' Syrup Swimming
#' 
#' Relative swimming speed (speed in syrup / speed in water) for 18 swimmers.
#' 
#' 
#' @name SyrupSwimming
#' @docType data
#' @format A data frame with 18 observations of one variable. \describe{
#' \item{relative.speed}{ratio of speed in syrup to speed in water} }
#' @references \url{http://www3.interscience.wiley.com/journal/109665380/issue}
#' @source Gettelfinger, B. and E. L. Cussler. 2004. Will Humans Swim Faster or
#' Slower in Syrup? \emph{AIChE Journal} 50: 2646-2647.
#' @keywords datasets
#' @examples
#' 
#' SyrupSwimming
#' histogram(~ relative.speed, SyrupSwimming)
#' dotplot(~ relative.speed, SyrupSwimming)
#' 
NULL





#' Causes of Teenage Deaths
#' 
#' Data from Table 1 (p. 14) on causes of death for all races, both sexes, ages
#' 15-19.
#' 
#' 
#' @name TeenDeaths
#' @docType data
#' @format A data frame with 11 observations on the following 2 variables.
#' \describe{ \item{cause}{a character vector} \item{deaths}{a
#' numeric vector} }
#' @source Anderson, R.N. 2001. Deaths: Leading causes for 1999. \emph{National
#' vital statistics reports} 49(11):1-88. National Center for Health
#' Statistics; Hyattsville, Maryland.
#' @keywords datasets
#' @examples
#' 
#' str(TeenDeaths)
#' TeenDeaths
#' 
#' barchart(deaths ~ cause, TeenDeaths, 
#'   horizontal = FALSE,
#'   ylab = "Number of Deaths",
#'   xlab = "Cause of Death", origin=0,
#'   scales = list(x = list(rot=45)))
#' 
#' barchart(deaths~ordered(cause, levels=cause), TeenDeaths, 
#'   horizontal = FALSE,
#'   ylab = "Number of Deaths",
#'   xlab = "Cause of Death", origin=0,
#'   scales=list(x=list(rot=45))
#'   )
#' 
NULL





#' Telomere Shortening
#' 
#' Telomere length (ratio) and years since their child's diagnosis with chronic
#' illness.
#' 
#' 
#' @name Telomeres
#' @docType data
#' @format A data frame with 39 observations on the following 2 variables.
#' \describe{ \item{years}{a numeric vector}
#' \item{telomere.length}{a numeric vector} }
#' @references \url{http://www.pnas.org/content/101/49/17312}
#' @source Epel, E.S., E.H. Blackburn, J. Lin, F.S. Dhabhar, N.E. Adler, J.D.
#' Morrow, and R.M. Cawthon. 2004. Accelerated telomere shortening in response
#' to life stress. \emph{Proceedings of the National Academy of Sciences (USA)}
#' 101: 17312-17315.
#' @keywords datasets
#' @examples
#' 
#' xyplot(years ~ telomere.length, Telomeres,
#'   xlab = "Time since diagnosis (years)",
#'   ylab = "Telomere length (ratio)"
#' )
#' 
#' 
NULL





#' Hypoxanthine and Time Since Death
#' 
#' Hypoxanthine levels in the vitreous humour of the eye and time since death
#' (hours) for 48 subjects.
#' 
#' 
#' @name TimeOfDeath
#' @docType data
#' @format A data frame with 48 observations on the following 2 variables.
#' \describe{ \item{hours}{a numeric vector}
#' \item{hypoxanthine}{a numeric vector} }
#' @source James, R.A., P.A. Hoadley, and B.G. Sampson. 1997. Determination of
#' postmortem interval by sampling vitreous humor. \emph{American Journal of
#' Forensic Medicine and Pathology} 18: 158-162.
#' @keywords datasets
#' @examples
#' 
#' xyplot(hypoxanthine ~ hours, TimeOfDeath, type=c('p','r'))
#' 
NULL





#' Right-handed Toads
#' 
#' Hypothetical probability of a toad being right-handed
#' 
#' 
#' @name Toads
#' @docType data
#' @format A data frame with 19 observations on the following 2 variables.
#' \describe{ \item{n.toads}{a numeric vector} \item{prob}{a
#' numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' Toads
#' # generate this data manually
#' cbind(0:18, dbinom(0:18, 18, 0.5))
#' xyplot(prob ~ n.toads, Toads, type = 'h', lwd = 4)
#' barchart(prob ~ n.toads, Toads, origin=0, horizontal=FALSE)
#' plotDist('binom', params = list(18,0.5), kind = 'hist')
#' 
#' 
NULL





#' Flower Length in Tobacco Plants
#' 
#' Distribution of flow lengths in F1 and F2 populations of \emph{Nicotiana}.
#' 
#' 
#' @name Tobacco
#' @docType data
#' @format A data frame with 13 observations on the following 3 variables.
#' \describe{ \item{flower.length}{a numeric vector of flower length in
#' mm} \item{f1.count}{a numeric vector of the number of F1 plants with
#' flower lengths in this size range} \item{f2.count}{a numeric vector
#' of the number of F2 plants with flower lengths in this size range} }
#' @seealso \code{\link{Tobacco2}}
#' @references \url{http://www.genetics.org/content/vol1/issue2/}
#' 
#' \url{http://en.wikipedia.org/wiki/Nicotiana}
#' @source East, E.M. 1916. Studies on size inheritance in \emph{Nicotiana}.
#' \emph{Genetics} 1: 164-176.
#' @keywords datasets
#' @examples
#' 
#' Tobacco
#' 
NULL





#' Flower Length in Tobacco Plants
#' 
#' Distribution of flow lengths in F1 and F2 populations of \emph{Nicotiana}.
#' 
#' 
#' @name Tobacco2
#' @docType data
#' @format A data frame with 617 observations on the following 2 variables.
#' \describe{ \item{flower.length}{a numeric vector}
#' \item{generation}{a factor with levels \code{F1} \code{F2}} }
#' @seealso \code{\link{Tobacco}}
#' @references \url{http://www.genetics.org/content/vol1/issue2/}
#' 
#' \url{http://en.wikipedia.org/wiki/Nicotiana}
#' @source East, E.M. 1916. Studies on size inheritance in \emph{Nicotiana}.
#' \emph{Genetics} 1: 164-176.
#' @keywords datasets
#' @examples
#' 
#' xtabs(~ flower.length + generation, Tobacco2)
#' bwplot(flower.length ~ generation, Tobacco2)
#' 
NULL





#' Radioactive Teeth
#' 
#' Actual birth year and birth year estimated from relative radioactivity of
#' the enamel for 20 samples.
#' 
#' 
#' @name ToothAge
#' @docType data
#' @format A data frame with 20 observations on the following 2 variables.
#' \describe{ \item{actual}{a numeric vector}
#' \item{estimated}{a numeric vector} }
#' @source Spalding, K.L., B.A. Buchholz, L.-E. Bergman, H. Druid, and J.
#' Frisén. 2005. Age written in teeth by nuclear tests. \emph{Nature} 437:
#' 333-334.
#' @keywords datasets
#' @examples
#' 
#' str(ToothAge)
#' xyplot(actual ~ estimated, ToothAge)
#' 
NULL





#' Tree Seedlings and Sunflecks
#' 
#' Fleck duration (min) and relative seedling growth rate (mm/mm/week) for 21
#' seedlings of \emph{Shorea leprosula}.
#' 
#' 
#' @name TreeSeedlings
#' @docType data
#' @format A data frame with 21 observations on the following 2 variables.
#' \describe{ \item{fleck.duration}{a numeric vector}
#' \item{growth}{a numeric vector} }
#' @references \url{http://jxb.oxfordjournals.org/cgi/content/short/56/411/469}
#' @source Leakey, A.D.B., J.D. Scholes, and M.C. Press. 2005. Physiological
#' and ecological significance of sunflecks for dipterocarp seedlings.
#' \emph{Journal of Experimental Botany} 56: 469-482.
#' @keywords datasets
#' @examples
#' 
#' str(TreeSeedlings)
#' splom(TreeSeedlings)
#' 
NULL





#' Frequencies of Fish Eaten by Trematode Infection Level
#' 
#' Frequencies of killifish (\emph{Fundulus parvipinnis}) eaten by birds
#' depending on level of infection by the trematode \emph{Euhaplorchis
#' californiensis}.
#' 
#' 
#' @name Trematodes
#' @docType data
#' @format A data frame with 141 observations on the following 2 variables.
#' \describe{ \item{infection.status}{a factor with levels:
#' \code{high}, \code{light}, and \code{uninfected}} \item{eaten}{a
#' factor with levels: \code{no} and \code{yes}} }
#' @source Lafferty, K.D. and A.K. Morris. 1996. Altered behavior of
#' parasitized killifish increases susceptibility to predation by bird final
#' hosts. \emph{Ecology} 77: 1390-1397.
#' @keywords datasets
#' @examples
#' 
#' demo(sec9.3)
#' 
NULL





#' Trillium Recruitment near Clearcuts
#' 
#' Recruitment of \emph{Trillium} and distance to nearest clearcut in eight
#' populations in southwestern Oregon.
#' 
#' 
#' @name Trillium
#' @docType data
#' @format A data frame with 8 observations on the following 3 variables.
#' \describe{ \item{population}{a numeric vector}
#' \item{edge.dist}{a numeric vector} \item{recruitment}{a
#' numeric vector} }
#' @source Jules, E.S. and B.J. Rathcke. 1999.  Mechanisms of reduced trillium
#' recruitment along edges of old-growth forest fragments. \emph{Conservation
#' Biology} 13: 784-793
#' @keywords datasets
#' @examples
#' 
#' str(Trillium)
#' splom(Trillium)
#' 
NULL





#' Truffle Distribution
#' 
#' Number of truffles per plot for 288 plots in an old growth forest in
#' northeastern California.
#' 
#' 
#' @name Truffles
#' @docType data
#' @format A data frame with 5 observations on the following 2 variables.
#' \describe{ \item{truffles}{a numeric vector} \item{count}{a
#' numeric vector} }
#' @references \url{http://www.fs.fed.us/psw/publications/watersj/waters2.PDF}
#' @source Waters, J.R., K.S. McKelvey, D.L. Luoam, and C.J. Zabel. 1997.
#' Truffle production in old-growth and mature fir stands in northeastern
#' California. \emph{Forest Ecology and Management} 96: 155-166.
#' @keywords datasets
#' @examples
#' 
#' Truffles
#' xyplot(count ~ truffles, Truffles, type='h', lwd=4)
#' barchart(count ~ truffles, Truffles, origin=0, horizontal=FALSE)
#' 
NULL





#' Dietary Learning in Tsetse Flies
#' 
#' Dietary conditioning \code{treatment} and subsequent proportion of tsetse
#' flies (\emph{Glossina palpalis}) feeding on cow blood in each of 13 cohorts.
#' 
#' 
#' @name TsetseLearning
#' @docType data
#' @format A data frame with 13 observations on the following 2 variables.
#' \describe{ \item{treatment}{a factor with levels \code{cow} and
#' \code{lizard}} \item{proportion.cow}{a numeric vector} }
#' @references
#' \url{http://rsbl.royalsocietypublishing.org/content/3/2/113.full}
#' @source \emph{inferred from} Bouyer, J., M. Pruvot, Z. Bengaly, P.M. Guerin,
#' and R. Lancelot. 2007. Learning influences host choice in tsetse.
#' \emph{Biology Letters} 3: 113-116.
#' @keywords datasets
#' @examples
#' 
#' xyplot(proportion.cow ~ treatment, TsetseLearning)
#' 
NULL





#' Number of Boys in Two-Child Families
#' 
#' The number of boys in a sample of 2,444 two-child families.
#' 
#' 
#' @name TwoKids
#' @docType data
#' @format A data frame with 3 observations on the following 2 variables.
#' \describe{ \item{num.boys}{a numeric vector} \item{count}{a
#' numeric vector} }
#' @references
#' \url{http://www.dartmouth.edu/~chance/chance_news/recent_news/chance_news_10.11.html#item13}
#' @source Rodgers, J.L. and D. Doughty. 2001. Does having boys or girls run in
#' the family? \emph{Chance Magazine} Fall 2001: 8-13.
#' @keywords datasets
#' @examples
#' 
#' TwoKids
#' observed <- TwoKids$count
#' expected <- c(585.3, 1221.4, 637.3)
#' chisq.test(observed, p = expected, rescale.p = TRUE)
#' 
#' # Alternate calculation, using Pr[male] = 0.512
#' # and rbinom. See Figure 5.7-1
#' n <- sum(observed)
#' pr.m <- 0.512
#' pr.f <- 0.488
#' 
#' # Calculate the probabilities of 0, 1, and 2 males
#' (pr.0 <- pr.f^2)
#' (pr.1 <- pr.m * pr.f + pr.f * pr.m)
#' (pr.2 <- pr.m^2)
#' 
#' set.seed(1)
#' (expected2 <- c(rbinom(1, n, pr.0),
#'                 rbinom(1, n, pr.1),
#'                 rbinom(1, n, pr.2)))
#' chisq.test(observed, p = expected2, rescale.p = TRUE)
#' 
NULL





#' Vampire Bat Bites
#' 
#' Numbers of cattle bitten by the cow's estrous cycle.
#' 
#' 
#' @name VampireBites
#' @docType data
#' @format A data frame with 4 observations on the following 3 variables.
#' \describe{ \item{estrous}{a factor with levels: \code{no} and
#' \code{yes}} \item{bitten}{a factor with levels: \code{no} and
#' \code{yes}} \item{count}{a numeric vector} }
#' @source Turner, D.C. 1975. \emph{The Vampire Bat: a Field Study in Behavior
#' and Ecology}. Johns Hopkins Press: Baltimore, MD.
#' @keywords datasets
#' @examples
#' 
#' demo(sec9.4)
#' 
NULL





#' Vasopressin Manipulation in the Meadow Vole
#' 
#' Time spent with a female (\code{percent}) in control and
#' vasopressin-enhanced groups (\code{treatment}) of meadow voles
#' (\emph{Microtus pennsylvanicus}).
#' 
#' @name VasopressinVoles
#' @docType data
#' @format A data frame with 31 observations on the following 2 variables.
#' \describe{ \item{treatment}{a factor with levels \code{control} and
#' \code{enhanced}} \item{percent}{a numeric vector} }
#' @source \emph{inferred from} Lim, M.M., Z. Wang, D.E.
#' Olazabal, X. Ren, E.F. Terwilliger, and L.J. Young. 2004.
#' Enhanced partner preference in a promiscuous species by manipulating the
#' expression of a single gene. \emph{Nature} 429: 754-757.
#' @keywords datasets
#' @examples
#' xyplot(percent ~ treatment, VasopressinVoles, type=c('p','a'))
#' bwplot(percent ~ treatment, VasopressinVoles)
#' 
NULL





#' Climbing Vines
#' 
#' Number of \code{climbing} and \code{nonclimbing} species within closely
#' related general of plants.
#' 
#' 
#' @name Vines
#' @docType data
#' @format A data frame with 48 observations on the following 2 variables.
#' \describe{ \item{climbing}{a numeric vector}
#' \item{nonclimbing}{a numeric vector} }
#' @references
#' \url{http://rspb.royalsocietypublishing.org/content/271/1552/2011.full.pdf}
#' @source Gianoli, E. 2004. Evolution of a climbing habit promotes
#' diversification in flowering plants. \emph{Proceedings of the Royal Society
#' of London, Series B, Biological Sciences} 271: 2011-2015.
#' @keywords datasets
#' @examples
#' 
#' xyplot(nonclimbing ~ climbing, Vines, scales=list(log=TRUE))
#' 
NULL





#' Home Range Size in Field Voles
#' 
#' Home range size size in field voles (\emph{Microtus agrestis}).
#' 
#' 
#' @name VoleDispersal
#' @docType data
#' @format A data frame with 5 observations on the following 3 variables.
#' \describe{ \item{homeranges}{a numeric vector}
#' \item{count}{a numeric vector} \item{sex}{a factor with
#' levels \code{female} and \code{male}} }
#' @source Sandell, M., J. Agrell, S. Erlinge, and J. Nelson. 1991. Adult
#' philopatry and dispersal in the field vole \emph{Microtus agrestis}.
#' \emph{Oecologia} 86: 153-158.
#' @keywords datasets
#' @examples
#' 
#' xtabs(count~sex+homeranges,VoleDispersal)
#' barchart( xtabs(count~sex+homeranges,VoleDispersal), origin=0, auto.key=TRUE)
#' barchart(count~sex+homeranges,VoleDispersal, origin=0)
#' barchart(count~sex,groups=homeranges,VoleDispersal, origin=0)
#' barchart(count~sex,groups=homeranges,VoleDispersal, origin=0,stack=TRUE)
#' 
NULL





#' Walking Stick Femur Length
#' 
#' Two measures of femur length \code{femur.length} for each of 25 
#' walking sticks (\emph{Timema cristinae}). Note that \code{specimen}
#' is not coded as a factor.
#' 
#' 
#' @name WalkingStickFemurs
#' @docType data
#' @format A data frame with 50 observations on the following 2 variables.
#' \describe{ \item{specimen}{a integer denoting specimen number.}
#' \item{femur.length}{a numeric vector of femur length} }
#' @references
#' \url{http://www.sfu.ca/biology/faculty/crespi/pdfs/96-Nosil&Crespi2006PNAS.pdf}
#' @source Nosil, P. and B.J. Crespi. 2006. Experimental evidence that
#' predation promotes divergence in adaptive radiation. \emph{Proceedings of
#' the National Academy of Sciences (USA)} 103: 9090-9095.
#' @keywords datasets
#' @examples
#' 
#' demo(sec15.6)
#' 
NULL





#' Walking Stick Head Width
#' 
#' Two measures of head width (\code{head.width}) for each of 25 walking sticks
#' (\emph{Timema cristinae}).
#' 
#' 
#' @name WalkingStickHeads
#' @docType data
#' @format A data frame with 50 observations on the following 2 variables.
#' \describe{ \item{specimen}{a factor with levels \code{1-25}}
#' \item{head.width}{a numeric vector} }
#' @references
#' \url{http://www.sfu.ca/biology/faculty/crespi/pdfs/96-Nosil&Crespi2006PNAS.pdf}
#' @source Nosil, P. and B.J. Crespi. 2006. Experimental evidence that
#' predation promotes divergence in adaptive radiation. \emph{Proceedings of
#' the National Academy of Sciences (USA)} 103: 9090-9095.
#' @keywords datasets
#' @examples
#' 
#' aggregate(head.width ~ specimen, data=WalkingStickHeads, mean) -> WS
#' histogram(~ head.width, WS)
#' 
NULL





#' Energetic Cost of Diving
#' 
#' Comparison of oxygen consumption in feeding vs. non-feeding dives of the
#' same length in the Weddell seal (\emph{Leptonychotes weddellii}).
#' 
#' 
#' @name WeddellSeals
#' @docType data
#' @format A data frame with 10 observations on the following 3 variables.
#' \describe{ \item{individual}{a numeric vector}
#' \item{oxygen.use.nonfeeding}{a numeric vector}
#' \item{oxygen.use.feeding}{a numeric vector} }
#' @references \url{http://jeb.biologists.org/cgi/content/full/207/6/973}
#' @source Williams, T.M., L.A. Fuiman, M. Horning, and R.W. Davis. 2004. The
#' cost of foraging by a marine predator, the Weddell seal \emph{Leptonychotes
#' weddellii}: pricing by the stroke. \emph{Journal of Experimental Biology}
#' 207: 973-982.
#' @keywords datasets
#' @examples
#' 
#' xyplot(oxygen.use.nonfeeding ~ oxygen.use.feeding, WeddellSeals)
#' 
NULL





#' Presidential "Wills"
#' 
#' Number of times a presidential candidate said "will," "shall," or "going to"
#' in presidential debates from 1960-2004 (years incomplete).
#' 
#' 
#' @name WillsDebates
#' @docType data
#' @format A data frame with 8 observations on the following 6 variables.
#' \describe{ \item{year}{year of presidential debate(s)}
#' \item{winner}{winner of the popular vote (may not be winner of
#' election)} \item{loser}{loser of popular vote (may not be loser of
#' election)} \item{winner.wills}{number of times will/shall used by
#' winner during debates} \item{loser.wills}{number of times will/shall
#' used by loser during debates} \item{diff.wills}{difference between
#' number of times will/shall used by two candidates} }
#' @keywords datasets
#' @examples
#' 
#' WillsDebates
#' 
NULL





#' Presidential "Wills"
#' 
#' Number of times a presidential candidate said "will," "shall," or "going to"
#' in presidential debates from 1960-2004 (years incomplete).
#' 
#' 
#' @name WillsPresidents
#' @docType data
#' @format A data frame with 16 observations on the following 3 variables.
#' \describe{ \item{candidate}{a character vector with the candidate's
#' name} \item{winner}{a factor with levels \code{n} \code{y}
#' indicating whether the candidate won the election \code{y} or not.}
#' \item{wills}{a numeric vector} \item{loser.wills}{a numeric
#' vector} \item{difference}{a numeric vector} \item{year}{a
#' numeric vector} }
#' @seealso \code{\link{WillsDebates}}
#' @keywords datasets
#' @examples
#' 
#' WillsPresidents
#' 
NULL





#' Wolf Tooth Measurements
#' 
#' Measurement (cm) of the distance between the canine and last molar teeth in
#' 35 wolves.
#' 
#' 
#' @name WolfTeeth
#' @docType data
#' @format A data frame with 35 observations of one variable.  \describe{
#' \item{length}{distance from canine to last molar teach (cm)} }
#' @references
#' \url{http://rspb.royalsocietypublishing.org/content/263/1372/849.abstract}
#' @source Whitlock, M. 1996. The heritability of fluctuating asymmetry and the
#' genetic control of developmental stability. \emph{Proceedings of the Royal
#' Society, Series B} 263: 849-853.
#' @keywords datasets
#' @examples
#' 
#' histogram(~ length, WolfTeeth)
#' 
NULL





#' Inbreeding in Wolves
#' 
#' Inbreeding coefficient and the number of pups produced in 24 mated pairs of
#' wolves (\emph{Canis lupus}) from 1983-2002.
#' 
#' 
#' @name Wolves
#' @docType data
#' @format A data frame with 24 observations on the following 2 variables.
#' \describe{ \item{inbreeding.coefficient}{a numeric vector}
#' \item{pups}{a numeric vector} }
#' @source Liberg, O.H., H. Andrén, H.-C. Pedersen, H. Sand, D. Sejberg, P.
#' Wabakken, M. Åkesson, and S. Bensch. 2005. Severe inbreeding depression in a
#' wild wolf (\emph{Canis lupus}) population. \emph{Biology Letters} 1: 17-20.
#' @keywords datasets
#' @examples
#' 
#' Wolves
#' xyplot(inbreeding.coefficient ~ jitter(pups, amount=0.15), Wolves) 
#' 
NULL





#' World Cup Goals
#' 
#' Number of goals per team during the 2002 World Cup.
#' 
#' 
#' @name WorldCup
#' @docType data
#' @format A data frame with 7 observations on the following 2 variables.
#' \describe{ \item{score}{a numeric vector} \item{count}{a
#' numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' xyplot(count ~ score, WorldCup, type='h', lwd=4)
#' barchart(count ~ score, WorldCup, origin=0, horizontal=FALSE)
#' 
NULL





#' Distribution of Wrasses
#' 
#' Number and sex of adult wrasses in a section of the Great Barrier Reef.
#' 
#' 
#' @name WrasseSexes
#' @docType data
#' @format A data frame with 3 observations on the following 3 variables.
#' \describe{ \item{males}{a numeric vector} \item{females}{a
#' numeric vector} \item{count}{a numeric vector} }
#' @keywords datasets
#' @examples
#' 
#' xtabs(count ~ males + females, WrasseSexes)
#' 
NULL





#' Yeast Regulatory Genes
#' 
#' Number of genes regulated by 109 yeast regulatory genes.
#' 
#' 
#' @name YeastGenes
#' @docType data
#' @format A data frame with 6 observations on the following 2 variables.
#' \describe{ \item{genes.controlled}{a numeric vector}
#' \item{count}{a numeric vector} }
#' @source Guelzim, N., S. Bottani, P. Bourgine and F. Képès. 2002. Topological
#' and causal structure of the yeast transcriptional regulatory network.
#' \emph{Nature Genetics} 31: 60-63.
#' @keywords datasets
#' @examples
#' 
#' str(YeastGenes)
#' barchart(count ~ genes.controlled , origin=0, YeastGenes, horizontal=FALSE)
#' 
NULL





#' Mate Preference in Zebra Finches
#' 
#' Percentage of time that a female spent next to a carotenoid-supplemented
#' male Zebra Finch compared to his non-supplemented brother.
#' 
#' 
#' @name ZebraFinchBeaks
#' @docType data
#' @format A numeric vector with 10 observations.
#' @references
#' \url{http://www.sciencemag.org/cgi/content/abstract/300/5616/125}
#' @source Blount, J.D., N.B. Metcalfe, T.R. Birkhead, P.F. Surai. 2003.
#' Carotenoid modulation of immune function and sexual attractiveness in Zebra
#' Finches. \emph{Science} 300: 125-127.
#' @keywords datasets
#' @examples
#' 
#' ZebraFinchBeaks
#' 
NULL





#' Zebra Finch Carotenoids
#' 
#' Data on cell-mediated immunocompetence (\code{PHA}) and humoral immunity
#' (\code{SRBC}) in Zebra Finches that received supplemental carotenoids
#' (\code{CAROT}) and those that did not (\code{NO}).
#' 
#' 
#' @name ZebraFinches
#' @docType data
#' @format A data frame with 20 observations on the following 3 variables.
#' \describe{ \item{treatment}{a factor with levels: \code{CAROT} and
#' \code{NO}} \item{PHA}{a numeric vector} \item{SRBC}{a
#' numeric vector} }
#' @source McGraw, K.J. and D.R. Ardia. 2003. Carotenoids, immunocompetence,
#' and the information content of sexual colors: an experimental test.
#' \emph{The American Naturalist} 162: 704-712.
#' @keywords datasets
#' @examples
#' 
#' ZebraFinches
#' 
NULL





#' Home Range Size and Mortality
#' 
#' Home range size (\eqn{\log_{10}}{log10}) and captive infant mortality (%)
#' for 20 species of carnivores.
#' 
#' 
#' @name ZooMortality
#' @docType data
#' @format A data frame with 20 observations on the following 2 variables.
#' \describe{ \item{log.homerange}{a numeric vector}
#' \item{mortality}{a numeric vector} }
#' @source Clubb, R. and G. Mason. 2003. Captivity effects on wide ranging
#' carnivores. \emph{Nature} 425: 473-474.
#' @keywords datasets
#' @examples
#' 
#' str(ZooMortality)
#' 
NULL





#' Zooplankton Depredation
#' 
#' Diversity of zooplankton (\code{zooplankton}) prey in each of 5 replicate
#' blocks (\code{block}) of three treatment levels (\code{treatment}). By
#' default, \code{block} is not coded as a factor.
#' 
#' 
#' @name Zooplankton
#' @docType data
#' @format A data frame with 15 observations on the following 3 variables.
#' \describe{ \item{treatment}{a factor with levels \code{control},
#' \code{high}, and \code{low}} \item{zooplankton}{a numeric vector}
#' \item{block}{a numeric vector} }
#' @source \emph{inferred from} Svanbäck, R. and D.I. Bolnick. 2007.
#' Intraspecific competition drives increased resource use diversity within a
#' natural population. \emph{Proceedings of the Royal Society of London Series
#' B, Biological Sciences} 274: 839-844.
#' @keywords datasets
#' @examples
#' 
#' Zooplankton
#' 
#' Zooplankton$block <- factor(Zooplankton$block)
#' str(Zooplankton)
#' 
#' aov.fit <- aov(zooplankton ~ block + treatment,
#'   data = Zooplankton)
#' summary(aov.fit)
#' 
NULL

