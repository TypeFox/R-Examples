#' Soc.ca a package for specific correspondence analysis
#' 
#' This package is optimized to the needs of scientists within the social 
#' sciences. The soc.ca package produces specific and class specific multiple
#' correspondence analysis on survey-like data. Soc.ca is optimized to only give
#' the most essential statistical output sorted so as to help in analysis. 
#' Seperate functions exists for near publication-ready plots and tables.
#' 
#' We are in debt to the work of others, especially Brigitte Le Roux and Henry
#' Rouanet for the mathematical definitions of the method and their examples.
#' Furthermore this package was initially based on code from the ca package
#' written by Michael Greenacre and Oleg Nenadic.
#' 
#' If you are looking for features that are absent in soc.ca, it may be 
#' available in some of these packages for correspondence analysis: \pkg{ca},
#' \pkg{anacor} and \pkg{FactoMineR}.
#' 
#' @references Le Roux, Brigitte, and Henry Rouanet. 2010. Multiple 
#'   correspondence analysis. Thousand Oaks: Sage.
#' @references Le Roux, Brigitte, and Henry Rouanet. 2004. Geometric Data 
#'   Analysis from Correspondence Analysis to Structured Data Analysis. 
#'   Dordrecht: Kluwer Academic Publishers.
#'   
#' @docType package
#' @name soc.ca
#' @import ggplot2
#' @import gridExtra
#' @import ellipse
#' @import stats
#' @import utils
#' @import shiny
#' @import reshape2
#' @import ggrepel
#' @examples
#' data(taste)
#' # Create a data frame of factors containing all the active variables 
#' taste          <- taste[which(taste$Isup == 'Active'), ]
#'
#' attach(taste)
#' active         <- data.frame(TV, Film, Art, Eat)
#' sup            <- data.frame(Gender, Age, Income)
#' detach(taste)
#' 
#' # Runs the analysis
#' result         <- soc.mca(active, sup)
NULL

#'Directors dataset
#'
#'Prosopographical data on the top 100 CEO's from the 82 largest Danish 
#'corporations.
#'@details The directors dataset is prosopographical data collected from a wide 
#'  array of sources on biographic and corporate information. Sources include 
#'  the Danish variant of Who's Who (Blaa Bog), a private business information 
#'  database (Greens Erhvervsinformation), journalistic portrait articles, 
#'  article search engines, bibliographic databases and financial reports. CEOs 
#'  from 82 corporations were selected according to their position as CEO in 
#'  December 2007. 18 executives are included on other criteria, taking into 
#'  account the magnitude of the corporations and issues regarding ownership and
#'  control, resulting in a final population of 100 CEOs. The 82 corporations 
#'  have formal ownership and management located in Denmark and were selected 
#'  through either financial capital, measured as having a turnover of over five
#'  billion DKK (650 million Eur.), or organizational capital, defined as having
#'  at least 5000 employees; 34 corporations were included on both criteria, 45 
#'  on financial capital and three on organizational capital alone. To avoid 
#'  including investors, rather than executives, a minimum of 500 employees was 
#'  also required, excluding 12 firms. Companies acting only as subsidiaries 
#'  were also excluded. Data is for public use  and no author permission is 
#'  needed, but we would love to hear from you if you find the data useful. The 
#'  following example is based on the analysis from the article: "A Very 
#'  Economic Elite: The Case of the Danish Top CEOs".
#'  
#'@name directors
#'@docType data
#'@author Christoph Ellersgaard
#'@author Anton Grau Larsen
#'@references Ellersgaard, Christoph, Anton Grau Larsen, og Martin D. Munk. 
#'  2012. "A Very Economic Elite: The Case of the Danish Top CEOs". Sociology.
#'@references Ellersgaard, Christoph Houman, og Anton Grau Larsen. 2010. 
#'  "Firmaets Maend". Master Thesis, Copenhagen: University of Copenhagen.
#'@references Ellersgaard, Christoph Houman, og Anton Grau Larsen. 2011. 
#'  "Kulturel kapital blandt topdirektoerer i Danmark - En domineret 
#'  kapitalform?" Dansk Sociologi 22(3):9-29.
#'@references Larsen, Anton Grau, og Christoph Houman Ellersgaard. 2012. "Status
#'  og integration paa magtens felt for danske topdirektoerer". Praktiske 
#'  Grunde. Nordisk tidsskrift for kultur- og samfundsvidenskab 2012(2-3).
#'@keywords data
#' @examples
#' \dontrun{
#' data(directors)
#' attach(directors)
#' 
#' 
#' active      <- data.frame(careerprofile_maclean_cat, careerfoundation_maclean_cat,
#'                           years_between_edu_dir_cat, time_in_corp_before_ceo_cat,
#'                           age_as_ceo_cat, career_changes_cat2, mba, abroad, hd, phd,
#'                           education, author, placeofbirth, familyclass_bourdieu,
#'                           partnersfamily_in_whoswho, family_in_whoswho)
#' 
#' sup       	<- data.frame(size_prestige, ownership_cat_2, sector, location)
#' 
#' id          <- navn
#' 
#' options(passive = c("MISSING", "Missing", "Irrelevant", "residence_value_cat2: Udlandet"))
#' 
#' result      <- soc.mca(active, sup, id)
#' 
#' result
#' 
#' # Contribution
#' contribution(result, 1)
#' contribution(result, 2)
#' contribution(result, 3)
#' contribution(result, 1, all = TRUE)
#' contribution(result, 1, indices = TRUE)
#' contribution(result, 1, mode = "mod")
#' contribution(result, mode = "variable")
#' 
#' # Individuals
#' contribution(result, 1, mode = "ind")
#' contribution(result, 2, mode = "ind")
#' 
#' 
#' # Table of variance
#' variance(result)
#' 
#' # Invert
#' result      <- invert(result, c(1, 2, 3))
#' 
#' # Export and assign label
#' # export.label(result)
#' 
#' # result      <- assign.label(result,
#' #  file = "https://raw.github.com/Rsoc/soc.ca/master/extra/director_labels.csv")
#' 
#' 
#' 
#' # Add.n
#' result      <- add.to.label(result)
#' contribution(result, 2)
#' 
#' 
#' # The result object or "soc.ca" object
#' str(result)
#' dim1 <- result$coord.ind[, 1]
#' qplot(dim1)
#' 
#' # Quadrant
#' quad      <- create.quadrant(result)
#' table(quad)
#' quad      <- create.quadrant(result, cut.min = 0, cut.max = 0)
#' table(quad)
#' 
#' 
#' # Map of individuals
#' map.ind(result)
#' map.ind(result, dim = c(2, 1), label = TRUE)
#' map.ind(result, dim = c(2, 1), point.size = 3, point.shape = 2)
#' map.ind(result, dim = c(2, 1), map.title = "The top 100 Danish CEO's",
#' point.color = quad)
#' # Map of the individuals colored by contribution
#' map.ind(result, point.color = result$ctr.ind[, 1],
#'  point.shape = 18) + scale_color_continuous(low = "white", high = "red")
#' 
#' 
#' # Map of contributing modalities
#' map.ctr(result, dim = c(2, 1))
#' map.ctr(result, dim = c(2, 1), ctr.dim = 2)
#' map.ctr(result, point.size = 3)
#'
#' map.active(result, dim = c(2, 1))
#' map.sup(result, dim = c(2, 1))
#' 
#' # Plot.list
#' 
#' # Selecting specific active modalities
#' select      <- c("Career start: Corporation (n:57)", "No Phd (n:92)")
#' boo.select  <- match(select, result$names.mod)
#' map.select(result, list.mod = boo.select)
#' 
#' highcor     <- which(result$cor.mod[, 1] >= 0.2) 
#' map.select(result, list.mod = highcor)
#' 
#' # Selecting specific supplementary modalities
#' 
#' highdim3    <- which(sqrt(result$coord.sup[, 3]^2) >= 0.5)
#' map.select(result, list.sup = highdim3)
#' 
#' # Selecting specific individuals based on a certain criteria
#' 
#' forfatter   <- author == "Forfatter"
#' map.select(result, list.ind = forfatter)
#' 
#' # Combining it all
#' map.select(result, list.mod = highcor, list.sup = highdim3, list.ind = forfatter)
#' 
#' # Add points to an existing plot
#' ctrplot     <- map.ctr(result, ctr.dim = 1, point.color = "red")
#' map.add(result, ctrplot, data.type = "ctr", ctr.dim = 2, point.color = "blue")
#' 
#' # Using the list option in add.points
#' forfatter    <- author == "Forfatter"
#' map.add(result, ctrplot, data.type = "select", list.ind = forfatter, colour = "purple")
#' 
#' # Using the list option in add.points to add labels to only a part of the cloud of individuals
#' forfatter     <- author == "Forfatter"
#' notforfatter  <- author != "Forfatter"
#' map.forfatter <- map.select(result, list.ind = notforfatter, label = FALSE)
#' map.forfatter
#' map.forfatter <- map.add(result, map.forfatter, data.type = "select", list.ind = forfatter)
#' map.forfatter
#' 
#' # Plotting all the modalities of one individual
#' result2       <- soc.ca(active, sup, id)
#' individual    <- which(id == "Lars Larsen")
#' ind.mat       <- indicator(active)
#' modalities    <- names(which(ind.mat[individual, ] == 1))
#' mod.ind       <- match(modalities, result2$names.mod)
#' 
#' lars          <- map.select(result2, list.mod = mod.ind)
#' map.add(result2, lars, data.type = "select", list.ind = individual, colour = "red")
#' 
#' # Adding concentration ellipses to an existing plot
#' el.forfatter  <- map.ellipse(result, map.forfatter, author)
#' el.forfatter
#' }
NULL

#' Taste dataset
#' 
#' The taste example dataset used by Le Roux & Rouanet(2010):
#' @return
#' The variables included in the dataset:
#' \item{Preferred TV program}{(8 categories): news, comedy, police, nature, sport, films, drama, soap operas}
#' \item{Preferred Film}{(8 categories): action, comedy, costume drama, documentary, horror, musical, romance, SciFi}
#' \item{Preferred type of Art}{(7 categories): performance, landscape, renaissance, still life, portrait, modern, impressionsism}
#' \item{Preferred place to Eat out}{(6 categories): fish & chips, pub, Indian restuarant, Italian restaurant, French restaurant, steak house}
#' @name taste
#' @docType data
#' @author Brigitte Le Roux
#' @references Le Roux, Brigitte, Henry Rouanet, Mike Savage, og Alan Warde. 2008. "Class and Cultural Division in the UK". Sociology 42(6):1049-1071.
#' @references Le Roux, B., og H. Rouanet. 2010. Multiple correspondence analysis. Thousand Oaks: Sage.
#' @keywords data
#' @examples
#' \dontrun{
#' # The taste example
#' data(taste)
#' data_taste           <- taste[which(taste$Isup == 'Active'), ]
#' active               <- data.frame(data_taste$TV, data_taste$Film, data_taste$Art, data_taste$Eat)
#' sup                  <- data.frame(data_taste$Gender, data_taste$Age, data_taste$Income)
#' 
#' # Multiple Correspondence Analysis
#' result.mca     <- soc.mca(active, sup)
#' str(result.mca)
#' result.mca
#' 
#' variance(result.mca) # See p.46 in Le Roux(2010)
#' 
#' contribution(result.mca, 1)
#' contribution(result.mca, 2)
#' contribution(result.mca, 1:3, mode = "variable")
#' 
#' map.active(result.mca, point.fill = result.mca$variable)
#' map.active(result.mca,
#'  map.title="Map of active modalities with size of contribution to 1. dimension",
#'  point.size=result.mca$ctr.mod[, 1])
#' map.active(result.mca, 
#'  map.title="Map of active modalities with size of contribution to 2. dimension",
#'  point.size=result.mca$ctr.mod[, 2])
#' 
#' map.ind(result.mca)
#' map.ind(result.mca, dim=c(1, 2), point.color=result.mca$ctr.ind[, 1],
#'  point.shape=18) + scale_color_continuous(low="white", high="black")
#' 
#' # Plot of all dublets
#' map.ind(result.mca, map.title="Map of all unique individuals", point.color=duplicated(active))
#' map.ind(result.mca, map.title="Map with individuals colored by the TV variable",
#'  point.color=active$TV)
#' 
#' # Ellipse 
#' map             <- map.ind(result.mca)
#' map.ellipse(result.mca, map, as.factor(data_taste$Age == '55-64'))
#' 
#' ##### Specific Multiple Correspondence Analysis
#' options(passive= c("Film: CostumeDrama", "TV: Tv-Sport"))
#' result.smca     <- soc.mca(active, sup)
#' result.smca
#' result.smca$names.passive
#' 
#' ##### Class Specific Correspondence Analysis
#' options(passive=NULL)
#' 
#' class.age     <- which(data_taste$Age == '55-64')
#' result.csca   <- soc.csa(result.mca, class.age, sup)
#' str(result.csca)
#' # Correlations
#' csa.measures(result.csca)
#' variance(result.csca)
#' contribution(result.csca, 1)
#' contribution(result.csca, 2)
#' contribution(result.csca, 1:3, mode = "variable")
#' 
#' # Plots
#' map.ind(result.csca)
#' map.csa.mca(result.csca)
#' map.csa.mca.array(result.csca)
#' }
NULL