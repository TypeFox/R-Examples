#' Species-by-Trait Matrix for Late Ordovician Marine Fossils.
#'
#' Sample data set of life habit codings (functional traits) for fossil taxa
#' from the Late Ordovician (Type Cincinnatian) Kope and Waynesville Formations
#' from Ohio, Indiana, and Kentucky (U.S.A.). The faunal list was compiled from
#' the Paleobiololgy Database (\url{paleobiodb.org/}).
#'
#' @name KWTraits
#' @rdname KWTraits
#'
#' @details Binary traits are coded with 0=absent and 1=present. Five ordered
#'   numeric traits (body volume, mobility, distance from seafloor
#'   [stratification]) were rescaled to range from 0 to 1 with discrete bins at
#'   equally spaced intermediate values.
#'
#'   See Novack-Gottshall (2007: especially online Supplementary Appendix A at
#'   \url{www.ben.edu/faculty/pnovack-gottshall/2007_Novack-Gottshall_Ecospace.pdf})
#'    for definition each functional trait, justifications, explanations, and
#'   examples. Novack-Gottshall (2007: Supplementary Appendix B; In press:
#'   Supplementary Appendix A) provides examples of how traits were coded using
#'   inferences derived from functional morphology, body size, ichnology,
#'   \emph{in situ} preservation, biotic associations recording direct
#'   interactions, and interpretation of geographic and depositional environment
#'   patterns.
#'
#'   Indeterminate taxa (e.g., trepostome bryozoan indet. or
#'   \emph{Platystrophia} sp.) that occurred within individual samples within
#'   these formations were excluded from the aggregate species pool unless their
#'   occurrence was the sole member of that taxon. Such indeterminate taxa and
#'   genera lacking a species identification were coded for a particular state
#'   only when all other members of that taxon within the Kope-Waynesville
#'   species pool unanimously shared that common state; otherwise, the state was
#'   listed as NA (missing).
#'
#' @format A data frame with 237 rows (taxa) and 40 columns (3 taxonomic
#'   identifiers and 37 functional traits): \describe{\item{Class}{Taxonomic
#'   class(character)} \item{Genus}{Taxonomic genus (character)}
#'   \item{sp.}{Taxonomic species (character)} \item{SEXL}{Sexual reproduction
#'   (binary)} \item{ASEX}{Asexual reproduction (binary)} \item{BVOL}{Skeletal
#'   body volume of typical adult (ordered numeric with 7 bins). Estimated using
#'   methods of Novack-Gottshall (2008).\itemize{ \item 1.000: >= 100 cm^3\item
#'   0.833: 100-10 cm^3 \item 0.667: 10-1 cm^3 \item 0.500: 1-0.1 cm^3 \item
#'   0.333: 0.1-0.01 cm^3 \item 0.167: 0.01-0.001 cm^3 \item 0: < 0.001 cm^3}}
#'   \item{BIOT}{Biotic substrate composition (binary)} \item{LITH}{Lithic
#'   substrate composition (binary)} \item{FLUD}{Fluidic medium (binary)}
#'   \item{HARD}{Hard substrate consistency (binary)} \item{SOFT}{Soft substrate
#'   consistency (binary)} \item{INSB}{Insubstantial medium consistency
#'   (binary)} \item{SPRT}{Supported on other object (binary)}
#'   \item{SSUP}{Self-supported (binary)} \item{ATTD}{Attached to substrate
#'   (binary)} \item{FRLV}{Free-living (binary)} \item{MOBL}{Mobility (ordered
#'   numeric with 5 bins): \itemize{ \item 1: habitually mobile \item 0.75:
#'   intermittently mobile \item 0.50: facultatively mobile \item 0.25:
#'   passively mobile (i.e., planktonic drifting) \item 0: sedentary
#'   (immobile)}} \item{ABST}{Primary microhabitat stratification: absolute
#'   distance from seafloor (ordered numeric with 5 bins): \itemize{ \item 1: >=
#'   100 cm \item 0.75: 100-10 cm: \item 0.50: 10-1 cm\item 0.25: 1-0.1 cm \item
#'   0: <0.1 cm}} \item{AABS}{Primary microhabitat is above seafloor (i.e.,
#'   epifaunal)} \item{IABS}{Primary microhabitat is within seafloor (i.e.,
#'   infaunal)} \item{RLST}{Immediately surrounding microhabitat stratification:
#'   relative distance from substrate (ordered numeric with 5 bins): \itemize{
#'   \item 1: >= 100 cm \item 0.75: 100-10 cm: \item 0.50: 10-1 cm \item 0.25:
#'   1-0.1 cm \item 0: <0.1 cm}} \item{AREL}{Lives above immediate substrate}
#'   \item{IREL}{Lives within immediate substrate} \item{FAAB}{Food is above
#'   seafloor} \item{FIAB}{Food is within seafloor} \item{FAST}{Primary feeding
#'   microhabitat stratification: absolute distance of food from seafloor
#'   (ordered numeric with 5 bins): \itemize{ \item 1: >= 100 cm \item 0.75:
#'   100-10 cm: \item 0.50: 10-1 cm \item 0.25: 1-0.1 cm \item 0: <0.1 cm}}
#'   \item{FARL}{Food is above immediate substrate} \item{FIRL}{Food is within
#'   immediate substrate} \item{FRST}{Immediately surrounding feeding
#'   microhabitat stratification: relative distance of food from substrate
#'   (ordered numeric with 5 bins): \itemize{ \item 1: >= 100 cm \item 0.75:
#'   100-10 cm: \item 0.50: 10-1 cm\item 0.25: 1-0.1 cm \item 0: <0.1 cm}}
#'   \item{AMBT}{Ambient foraging habit} \item{FILT}{Filter-feeding foraging
#'   habit} \item{ATTF}{Attachment-feeding foraging habit}
#'   \item{MASS}{Mass-feeding foraging habit} \item{RAPT}{Raptorial foraging
#'   habit} \item{AUTO}{Autotrophic diet} \item{MICR}{Microbivorous (bacteria,
#'   protists, algae) diet} \item{CARN}{Carnivorous diet} \item{INCP}{Food has
#'   incorporeal physical condition} \item{PART}{Food consumed as particulate
#'   matter} \item{BULK}{Food consumed as bulk matter}}
#'
#' @source Novack-Gottshall, P.M. In review at \emph{Paleobiology}, submitted
#'   Oct. 5, 2015. General models of ecological diversification. II. Simulations
#'   and empirical applications.
#'
#' @references Novack-Gottshall, P.M. 2007. Using a theoretical ecospace to
#'   quantify the ecological diversity of Paleozoic and modern marine biotas.
#'   \emph{Paleobiology} 33: 274-295.
#' @references Novack-Gottshall, P.M. 2008. Using simple body-size metrics to
#'   estimate fossil body volume: empirical validation using diverse Paleozoic
#'   invertebrates. \emph{PALAIOS} 23(3):163-173.
#' @references Novack-Gottshall, P.M. In review at \emph{Paleobiology},
#'   submitted Oct. 5, 2015. General models of ecological diversification. II.
#'   Simulations and empirical applications.
#' @references Villeger, S., P. M. Novack-Gottshall, and D. Mouillot. 2011. The
#'   multidimensionality of the niche reveals functional diversity changes in
#'   benthic marine biotas across geological time. \emph{Ecology Letters}
#'   14(6):561-568.
#'
#' @keywords datasets
"KWTraits"
