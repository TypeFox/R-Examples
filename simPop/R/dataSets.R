#' Population totals Region times Gender for Austria 2006
#'
#' Population characteristics Region times Gender from Austria.
#'
#'
#' @name totalsRG
#' @aliases totalsRG totalsRGtab
#' @docType data
#' @format totalsRG: A data frame with 18 observations on the following 3
#' variables.  \describe{ \item{list("rb090")}{gender; a factor with levels
#' \code{female} \code{male}} \item{list("db040")}{region; a factor with levels
#' \code{Burgenland} \code{Carinthia} \code{Lower Austria,} \code{Salzburg}
#' \code{Styria} \code{Tyrol} \code{Upper Austria} \code{Vienna}
#' \code{Vorarlberg}} \item{list("Freq")}{totals; a numeric vector} }
#' totalsRGtab: a two-dimensional table holding the same information
#' @source StatCube - statistical data base,
#' \url{http://www.statistik.at}
#' @keywords datasets
#' @examples
#'
#' data(totalsRG)
#' totalsRG
#' data(totalsRGtab)
#' totalsRGtab
#'
NULL


#' Synthetic EU-SILC data
#'
#' This data set is synthetically generated from real Austrian EU-SILC
#' (European Union Statistics on Income and Living Conditions) data.
#'
#' The data set is used as population data in some of the examples in package
#' \code{simFrame}.  Note that it is included for illustrative purposes only.
#' It consists of 25 000 households, hence it does not represent the true
#' population sizes of Austria and its regions.
#'
#' Only a few of the large number of variables in the original survey are
#' included in this example data set.  Some variable names are different from
#' the standardized names used by the statistical agencies, as the latter are
#' rather cryptic codes.  Furthermore, the variables \code{hsize},
#' \code{eqsize}, \code{eqIncome} and \code{age} are not included in the
#' standardized format of EU-SILC data, but have been derived from other
#' variables for convenience.  Moreover, some very sparse income components
#' were not included in the the generation of this synthetic data set. Thus the
#' equivalized household income is computed from the available income
#' components.
#'
#' @name eusilcP
#' @aliases eusilcP
#' @docType data
#' @format A \code{data.frame} with 58 654 observations on the following 28
#' variables: \describe{ \item{hid}{integer; the household ID.}
#' \item{region}{factor; the federal state in which the household is
#' located (levels \code{Burgenland}, \code{Carinthia}, \code{Lower Austria},
#' \code{Salzburg}, \code{Styria}, \code{Tyrol}, \code{Upper Austria},
#' \code{Vienna} and \code{Vorarlberg}).} \item{hsize}{integer; the
#' number of persons in the household.} \item{eqsize}{numeric; the
#' equivalized household size according to the modified OECD scale.}
#' \item{eqIncome}{numeric; a simplified version of the equivalized
#' household income.} \item{pid}{integer; the personal ID.}
#' \item{id}{the household ID combined with the personal ID.  The first five
#' digits represent the household ID, the last two digits the personal ID (both
#' with leading zeros).} \item{age}{integer; the person's age.}
#' \item{gender}{factor; the person's gender (levels \code{male} and
#' \code{female}).} \item{ecoStat}{factor; the person's economic status
#' (levels \code{1} = working full time, \code{2} = working part time, \code{3}
#' = unemployed, \code{4} = pupil, student, further training or unpaid work
#' experience or in compulsory military or community service, \code{5} = in
#' retirement or early retirement or has given up business, \code{6} =
#' permanently disabled or/and unfit to work or other inactive person, \code{7}
#' = fulfilling domestic tasks and care responsibilities).}
#' \item{citizenship}{factor; the person's citizenship (levels
#' \code{AT}, \code{EU} and \code{Other}).} \item{py010n}{numeric;
#' employee cash or near cash income (net).} \item{py050n}{numeric;
#' cash benefits or losses from self-employment (net).}
#' \item{py090n}{numeric; unemployment benefits (net).}
#' \item{py100n}{numeric; old-age benefits (net).}
#' \item{py110n}{numeric; survivor's benefits (net).}
#' \item{py120n}{numeric; sickness benefits (net).}
#' \item{py130n}{numeric; disability benefits (net).}
#' \item{py140n}{numeric; education-related allowances (net).}
#' \item{hy040n}{numeric; income from rental of a property or land
#' (net).} \item{hy050n}{numeric; family/children related allowances
#' (net).} \item{hy070n}{numeric; housing allowances (net).}
#' \item{hy080n}{numeric; regular inter-household cash transfer
#' received (net).} \item{hy090n}{numeric; interest, dividends, profit
#' from capital investments in unincorporated business (net).}
#' \item{hy110n}{numeric; income received by people aged under 16
#' (net).} \item{hy130n}{numeric; regular inter-household cash transfer
#' paid (net).} \item{hy145n}{numeric; repayments/receipts for tax
#' adjustment (net).} \item{main}{logical; indicates the main income
#' holder (i.e., the person with the highest income) of each household.} }
#' @references Eurostat (2004) Description of target variables: Cross-sectional
#' and longitudinal. \emph{EU-SILC 065/04}, Eurostat.
#' @source This is a synthetic data set based on Austrian EU-SILC data from
#' 2006.  The original sample was provided by Statistics Austria.
#' @keywords datasets
#' @examples
#'
#' data(eusilcP)
#' summary(eusilcP)
#'
NULL





#' Synthetic EU-SILC survey data
#'
#' This data set is synthetically generated from real Austrian EU-SILC
#' (European Union Statistics on Income and Living Conditions) data.
#'
#' The data set consists of 4641 households and is used as sample data in some
#' of the examples in package \code{simPopulation}.  Note that it is included
#' for illustrative purposes only.  The sample weights do not reflect the true
#' population sizes of Austria and its regions.  The resulting population data
#' is about 100 times smaller than the real population size to save computation
#' time.
#'
#' Only a few of the large number of variables in the original survey are
#' included in this example data set.  The variable names are rather cryptic
#' codes, but these are the standardized names used by the statistical
#' agencies.  Furthermore, the variables \code{hsize}, \code{age} and
#' \code{netIncome} are not included in the standardized format of EU-SILC
#' data, but have been derived from other variables for convenience.
#'
#' @name eusilcS
#' @docType data
#' @format A data frame with 11725 observations on the following 18 variables.
#' \describe{ \item{db030}{integer; the household ID.}
#' \item{hsize}{integer; the number of persons in the household.}
#' \item{db040}{factor; the federal state in which the household is
#' located (levels \code{Burgenland}, \code{Carinthia}, \code{Lower Austria},
#' \code{Salzburg}, \code{Styria}, \code{Tyrol}, \code{Upper Austria},
#' \code{Vienna} and \code{Vorarlberg}).} \item{age}{integer; the
#' person's age.} \item{rb090}{factor; the person's gender (levels
#' \code{male} and \code{female}).} \item{pl030}{factor; the person's
#' economic status (levels \code{1} = working full time, \code{2} = working
#' part time, \code{3} = unemployed, \code{4} = pupil, student, further
#' training or unpaid work experience or in compulsory military or community
#' service, \code{5} = in retirement or early retirement or has given up
#' business, \code{6} = permanently disabled or/and unfit to work or other
#' inactive person, \code{7} = fulfilling domestic tasks and care
#' responsibilities).} \item{pb220a}{factor; the person's citizenship
#' (levels \code{AT}, \code{EU} and \code{Other}).}
#' \item{netIncome}{numeric; the personal net income.}
#' \item{py010n}{numeric; employee cash or near cash income (net).}
#' \item{py050n}{numeric; cash benefits or losses from self-employment
#' (net).} \item{py090n}{numeric; unemployment benefits (net).}
#' \item{py100n}{numeric; old-age benefits (net).}
#' \item{py110n}{numeric; survivor's benefits (net).}
#' \item{py120n}{numeric; sickness benefits (net).}
#' \item{py130n}{numeric; disability benefits (net).}
#' \item{py140n}{numeric; education-related allowances (net).}
#' \item{db090}{numeric; the household sample weights.}
#' \item{rb050}{numeric; the personal sample weights.} }
#' @references Eurostat (2004) Description of target variables: Cross-sectional
#' and longitudinal. \emph{EU-SILC 065/04}, Eurostat.
#' @source This is a synthetic data set based on Austrian EU-SILC data from
#' 2006.  The original sample was provided by Statistics Austria.
#' @keywords datasets
#' @examples
#'
#' data(eusilcS)
#' summary(eusilcS)
#'
NULL


#' Synthetic EU-SILC 2013 survey data
#'
#' This data set is synthetically generated from real Austrian EU-SILC
#' (European Union Statistics on Income and Living Conditions) data 2013.
#'
#' The data set consists of 5977 households and is used as sample data in some
#' of the examples in package \code{simPop}.  Note that it is included
#' for illustrative purposes only.  The sample weights do not reflect the true
#' population sizes of Austria and its regions.
#'
#' 62 variables of the original survey are
#' simulated for this example data set.  The variable names are rather cryptic
#' codes, but these are the standardized names used by the statistical
#' agencies.  Furthermore, the variables \code{hsize}, \code{age} and
#' \code{netIncome} are not included in the standardized format of EU-SILC
#' data, but have been derived from other variables for convenience.
#'
#' @name eusilc13puf
#' @docType data
#' @format A data frame with 13513 observations on the following 62 variables.
#' \describe{
#' \item{db030}{integer; the household ID.}
#' \item{hsize}{integer; the number of persons in the household.}
#' \item{db040}{factor; the federal state in which the household is
#' located (levels \code{Burgenland}, \code{Carinthia}, \code{Lower Austria},
#' \code{Salzburg}, \code{Styria}, \code{Tyrol}, \code{Upper Austria},
#' \code{Vienna} and \code{Vorarlberg}).}
#' \item{age}{integer; the
#' person's age.}
#' \item{rb090}{factor; the person's gender (levels
#' \code{male} and \code{female}).}
#' \item{pid}{personal ID}
#' \item{weight}{sampling weights}
#' \item{pl031}{factor; the person's
#' economic status (levels \code{1} = working full time, \code{2} = working
#' part time, \code{3} = unemployed, \code{4} = pupil, student, further
#' training or unpaid work experience or in compulsory military or community
#' service, \code{5} = in retirement or early retirement or has given up
#' business, \code{6} = permanently disabled or/and unfit to work or other
#' inactive person, \code{7} = fulfilling domestic tasks and care
#' responsibilities).}
#' \item{pb220a}{factor; the person's citizenship
#' (levels \code{AT}, \code{EU} and \code{Other}).}
#' \item{pb190}{for details, see Eurostat's code book}
#' \item{pe040}{for details, see Eurostat's code book}
#' \item{pl111}{for details, see Eurostat's code book}
#' \item{pgrossIncomeCat}{for details, see Eurostat's code book}
#' \item{pgrossIncome}{for details, see Eurostat's code book}
#' \item{hgrossIncomeCat}{for details, see Eurostat's code book}
#' \item{hgrossIncome}{for details, see Eurostat's code book}
#' \item{hgrossminusCat}{for details, see Eurostat's code book}
#' \item{hgrossminus}{for details, see Eurostat's code book}
#' \item{py010g}{for details, see Eurostat's code book}
#' \item{py021g}{for details, see Eurostat's code book}
#' \item{py050g}{for details, see Eurostat's code book}
#' \item{py080g}{for details, see Eurostat's code book}
#' \item{py090g}{for details, see Eurostat's code book}
#' \item{py100g}{for details, see Eurostat's code book}
#' \item{py110g}{for details, see Eurostat's code book}
#' \item{py120g}{for details, see Eurostat's code book}
#' \item{py130g}{for details, see Eurostat's code book}
#' \item{py140g}{for details, see Eurostat's code book}
#' \item{hy040g}{for details, see Eurostat's code book}
#' \item{hy050g}{for details, see Eurostat's code book}
#' \item{hy060g}{for details, see Eurostat's code book}
#' \item{hy070g}{for details, see Eurostat's code book}
#' \item{hy080g}{for details, see Eurostat's code book}
#' \item{hy090g}{for details, see Eurostat's code book}
#' \item{hy100g}{for details, see Eurostat's code book}
#' \item{hy110g}{for details, see Eurostat's code book}
#' \item{hy120g}{for details, see Eurostat's code book}
#' \item{hy130g}{for details, see Eurostat's code book}
#' \item{hy140g}{for details, see Eurostat's code book}
#' \item{rb250}{for details, see Eurostat's code book}
#' \item{p119000}{for details, see Eurostat's code book}
#' \item{p038003f}{for details, see Eurostat's code book}
#' \item{p118000i}{for details, see Eurostat's code book}
#' \item{aktivi}{for details, see Eurostat's code book}
#' \item{erwintensneu}{for details, see Eurostat's code book}
#' \item{rb050}{for details, see Eurostat's code book}
#' \item{pb040}{for details, see Eurostat's code book}
#' \item{hb030}{for details, see Eurostat's code book}
#' \item{px030}{for details, see Eurostat's code book}
#' \item{rx030}{for details, see Eurostat's code book}
#' \item{pb030}{for details, see Eurostat's code book}
#' \item{rb030}{for details, see Eurostat's code book}
#' \item{hx040}{for details, see Eurostat's code book}
#' \item{pb150}{for details, see Eurostat's code book}
#' \item{rx020}{for details, see Eurostat's code book}
#' \item{px020}{for details, see Eurostat's code book}
#' \item{hx050}{for details, see Eurostat's code book}
#' \item{eqInc}{for details, see Eurostat's code book}
#' \item{hy010}{for details, see Eurostat's code book}
#' \item{hy020}{for details, see Eurostat's code book}
#' \item{hy022}{for details, see Eurostat's code book}
#' \item{hy023}{for details, see Eurostat's code book}
#' }
#' @references Eurostat (2013) Description of target variables: Cross-sectional
#' and longitudinal.
#' @source This is a synthetic data set based on Austrian EU-SILC data from
#' 2013.  The original sample was provided by Statistics Austria.
#' @keywords datasets
#' @author Matthias Templ
#' @examples
#' data(eusilc13puf)
#' str(eusilc13puf)
NULL


#' Synthetic GLSS survey data
#'
#' This data set is synthetically generated from real GLSS (Ghana Living
#' Standards Survey) data.
#'
#' The data set consists of 8700 households and is used as sample data in some
#' of the examples in package \code{simPopulation}.  Note that it is included
#' for illustrative purposes only.  The sample weights do not reflect the true
#' population sizes of Ghana and its regions.  The resulting population data is
#' about 100 times smaller than the real population size to save computation
#' time.
#'
#' Only some of the variables in the original survey are included in this
#' example data set.  Furthermore, categories are aggregated for certain
#' variables due to the large number of possible outcomes in the original
#' survey data.
#'
#' @name ghanaS
#' @docType data
#' @format A data frame with 36970 observations on the following 14 variables.
#' \describe{ \item{hhid}{integer; the household ID.}
#' \item{hsize}{integer; the number of persons in the household.}
#' \item{region}{factor; the region in which the household is located
#' (levels \code{western}, \code{central}, \code{greater accra}, \code{volta},
#' \code{eastern}, \code{ashanti}, \code{brong ahafo}, \code{northern},
#' \code{upper east} and \code{upper west}).} \item{clust}{factor; the
#' enumeration area.} \item{age}{integer; the person's age.}
#' \item{sex}{factor; the person's sex (levels \code{male} and
#' \code{female}).} \item{relate}{factor; the relationship with the
#' household head (levels \code{head}, \code{spouse}, \code{child},
#' \code{grandchild}, \code{parent/parentlaw}, \code{son/daughterlaw},
#' \code{other relative}, \code{adopted child}, \code{househelp} and
#' \code{non_relative}).} \item{nation}{factor; the person's
#' nationality (levels \code{ghanaian birth}, \code{ghanaian naturalise},
#' \code{burkinabe}, \code{malian}, \code{nigerian}, \code{ivorian},
#' \code{togolese}, \code{liberian}, \code{other ecowas}, \code{other africa}
#' and \code{other}).} \item{ethnic}{factor; the person's ethnicity
#' (levels \code{akan}, \code{all other tribes}, \code{ewe}, \code{ga-dangbe},
#' \code{grusi}, \code{guan}, \code{gurma}, \code{mande} and
#' \code{mole-dagbani}).} \item{religion}{factor; the person's religion
#' (levels \code{catholic}, \code{anglican}, \code{presbyterian},
#' \code{methodist}, \code{pentecostal}, \code{spiritualist}, \code{other
#' christian}, \code{moslem}, \code{traditional}, \code{no religion} and
#' \code{other}).} \item{highest_degree}{factor; the person's highest
#' degree of education (levels \code{none}, \code{mlsc}, \code{bece},
#' \code{voc/comm}, \code{teacher trng a}, \code{teacher trng b}, \code{gce 'o'
#' level}, \code{ssce}, \code{gce 'a' level}, \code{tech/prof cert},
#' \code{tech/prof dip}, \code{hnd}, \code{bachelor}, \code{masters},
#' \code{doctorate} and \code{other}).} \item{occupation}{factor; the
#' person's occupation (levels \code{armed forces and other security
#' personnel}, \code{clerks}, \code{craft and related trades workers},
#' \code{elementary occupations}, \code{legislators, senior officials and
#' managers}, \code{none}, \code{plant and machine operators and assemblers},
#' \code{professionals}, \code{service workers and shop and market sales
#' workers}, \code{skilled agricultural and fishery workers}, and
#' \code{technicians and associate professionals}).}
#' \item{income}{numeric; the person's annual income.}
#' \item{weight}{numeric; the sample weights.} }
#' @references Ghana Statistical Service (2008) Ghana Living Standards Survey:
#' Report of the fifth round.
#' @source This is a synthetic data set based on GLSS data from 2006.  The
#' original sample was provided by Ghana Statistical Service.
#' @keywords datasets
#' @examples
#'
#' data(ghanaS)
#' summary(ghanaS)
#'
NULL





#' Extract and modify variables from population or sample data stored in an
#' object of class \code{\link{simPopObj-class}}.
#'
#' Using \code{\link{samp}} \code{\link{samp<-}} it is possible to extract or
#' rather modify variables of the sample data within slot \code{data} in slot
#' \code{sample} of the \code{\link{simPopObj-class}}-object. Using
#' \code{\link{pop}} \code{\link{pop<-}} it is possible to extract or rather
#' modify variables of the synthetic population within in slot \code{data} in
#' slot \code{sample} of the \code{\link{simPopObj-class}}-object.
#'
#'


#' Population totals Region times Gender for Austria 2006
#'
#' Population characteristics Region times Gender from Austria.
#'
#'
#' @name totalsRG
#' @aliases totalsRG totalsRGtab
#' @docType data
#' @format totalsRG: A data frame with 18 observations on the following 3
#' variables.  \describe{ \item{list("rb090")}{gender; a factor with levels
#' \code{female} \code{male}} \item{list("db040")}{region; a factor with levels
#' \code{Burgenland} \code{Carinthia} \code{Lower Austria,} \code{Salzburg}
#' \code{Styria} \code{Tyrol} \code{Upper Austria} \code{Vienna}
#' \code{Vorarlberg}} \item{list("Freq")}{totals; a numeric vector} }
#' totalsRGtab: a two-dimensional table holding the same information
#' @source StatCube - statistical data base,
#' \url{http://www.statistik.at/}
#' @keywords datasets
#' @examples
#'
#' data(totalsRG)
#' totalsRG
#' data(totalsRGtab)
#' totalsRGtab
#'
NULL
