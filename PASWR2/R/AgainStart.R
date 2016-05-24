#' The PASWR2 Package
#' 
#' Data sets and functions for \emph{Probability and Statistics with R}, Second Edition.
#' @name PASWR2-package
#' @docType package
#' @title The PASWR2 Package
#' @keywords package
NULL
#####################################################################################
#' @name AGGRESSION
#' @title TV and Behavior
#' @aliases AGGRESSION
#' @docType data
#' @description Data regarding the aggressive behavior in relation to exposure to violent television programs.
#' @format A data frame with 16 observations on the following two variables: 
#' \itemize{
#' \item \code{violence} (an integer vector)
#' \item \code{noviolence} (an integer vector)
#' }
#' @details This is data regarding aggressive behavior in relation to exposure to violent television programs from Gibbons (1997) with the following exposition: \dQuote{\ldots a group of children are matched as well as possible as regards home environment, genetic  factors, intelligence, parental attitudes, and so forth, in an effort to minimize factors other than TV that might influence a tendency for aggressive behavior.  In each of the resulting 16 pairs, one child is randomly selected to view the most violent shows on TV, while the other watches cartoons, situation comedies, and the like.  The children are then subjected to a series of tests designed to produce an ordinal measure of their aggression factors.} (pages 143-144)
#' @source Gibbons, J. D. (1977) \emph{Nonparametric Methods for Quantitavie Analysis}. American Science Press.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' AL <- reshape(AGGRESSION, varying = c("violence", "noviolence"), 
#' v.names = "aggression", direction = "long")
#' ggplot(data = AL, aes(x = factor(time), y = aggression, fill = factor(time))) + 
#' geom_boxplot() + labs(x = "") + scale_x_discrete(breaks = c(1, 2), 
#' labels = c("Violence", "No Violence")) + guides(fill = FALSE) + scale_fill_brewer()
#' rm(AL)
#' with(data = AGGRESSION, 
#' wilcox.test(violence, noviolence, paired = TRUE, alternative = "greater"))
#' @keywords datasets
"AGGRESSION"
#####################################################################################
#' @name APPLE
#' @title Apple Hardness
#' @aliases APPLE
#' @docType data
#' @description An experiment was undertaken where seventeen apples were randomly selected from an orchard (\code{fresh}) and measured for hardness.  Seventeen apples were also randomly selected from a warehouse (\code{warehouse}) where the apples had been stored for one week and measured for hardness. 
#' @format A data frame with 34 observations on the following two variables: 
#' \itemize{
#' \item \code{hardness} (hardness rating measured in \eqn{\texttt{kg}/\texttt{meter}^2} for both the \code{fresh} and \code{warehouse} apples)
#' \item \code{location} (\code{factor} with two levels \code{fresh} and \code{warehouse})  
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # ggplot2 approach
#' ggplot(data = APPLE, aes(sample = hardness)) + stat_qq() + facet_grid(. ~ location)
#' ggplot(data = APPLE, aes(sample = hardness, color = location)) + stat_qq()
#' ggplot(data = APPLE, aes(x = hardness, fill = location)) + geom_density(alpha = 0.4) + 
#' scale_fill_brewer()
#' # lattice approach
#' qqmath(~hardness|location, data = APPLE)
#' qqmath(~hardness, group = location, type = c('p', 'r'), auto.key = TRUE, data = APPLE)
#' @keywords datasets
"APPLE"
#####################################################################################
#' @name APTSIZE
#' @title Apartment Size
#' @aliases APTSIZE
#' @docType data
#' @description Size of apartments in Mendebaldea, Spain, and San Jorge, Spain
#' @format A data frame with 15 observations on the following two variables:
#' \itemize{
#' \item \code{size} (apartment size in square meters)
#' \item \code{location} (\code{factor} with two levels \code{SanJorge} and \code{Mendebaldea})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' p <- ggplot(data = APTSIZE, aes(x = location, y = size, fill = location)) + 
#' labs(x = "", y = "Apartment size (square meters)") + 
#' scale_x_discrete(breaks = c("Mendebaldea", "SanJorge"), 
#' labels =c("Mendebaldea", "San Jorge")) + scale_fill_brewer()
#' p + geom_boxplot()
#' # remove the legend
#' p + geom_boxplot() + guides(fill = FALSE)
#' # violin plot
#' p + geom_violin(scale = 'area') + guides(fill = FALSE)
#' p + geom_violin(scale = 'count') + guides(fill = FALSE)
#' p + geom_violin() + geom_boxplot(width = 0.15, fill = 'black') + guides(fill = FALSE) + 
#' stat_summary(fun.y = median, geom = "point", fill = "white", shape = 23, size = 3)
#' # dotplot
#' p + geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 3) + 
#' guides(fill = FALSE)
#' p + geom_boxplot(width = 0.4) + geom_dotplot(binaxis = "y", stackdir = "center", 
#' binwidth = 3) + guides(fill = FALSE) + scale_fill_brewer(type = "qual", palette = 1)
#' # base graphics
#' boxplot(size ~ location, data = APTSIZE, col = c("red", "yellow"), 
#' ylab = "Apartment size (square meters)")
#' @keywords datasets
"APTSIZE"
#####################################################################################
#' @name BABERUTH
#' @title George Herman Ruth
#' @aliases BABERUTH
#' @docType data
#' @description Baseball statistics for George Herman Ruth (The Bambino or the Sultan of Swat)
#' @format A data frame with 22 observation of the following 14 variables:
#' \itemize{
#' \item\code{year} (year in which the season occurred) 
#' \item \code{team} (team for which he played \code{Bos-A}, \code{Bos-N}, or \code{NY-A}) 
#' \item \code{g} (games played) 
#' \item \code{ab} (at bats) 
#' \item \code{r} (runs scored) 
#' \item \code{h} (hits) 
#' \item \code{X2b} (doubles) 
#' \item \code{X3b} (triples) 
#' \item \code{hr} (home runs) 
#' \item \code{RBI} (runs batted in) 
#' \item \code{sb} (stolen bases) 
#' \item \code{bb} (base on balls or walks) 
#' \item \code{ba} (batting average = h/ab)
#' \item \code{slg} (slugging percentage = total bases/at bats)
#' }
#' @source \url{http://www.baseball-reference.com/about/bat_glossary.shtml}
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = BABERUTH, aes(x = ba)) + geom_histogram(binwidth = 0.03) + 
#' facet_grid(team ~ .) + labs(x = "Batting average")
#' ggplot(data = BABERUTH, aes(x = g, y = ab, color = rbi)) + geom_point() + 
#' labs(x = "Number of Games Played", y = "Times at Bat", color = "Runs\n Batted In", 
#' title = "George Herman Ruth")
#' @keywords datasets
"BABERUTH"
#####################################################################################
#' @name BAC
#' @title Blood Alcohol Content
#' @aliases BAC
#' @docType data
#' @description Two volunteers weighing 180 pounds each consumed a twelve ounce beer every fifteen minutes for one hour. One hour after the fourth beer was consumed, each volunteer's blood alcohol was measured with ten different breathalyzers from the same company. The numbers recorded in data frame \code{BAC} are the sorted blood alcohol content values reported with breathalyzers from company \code{X} and company \code{Y}. 
#' @format A data frame with 10 observations of the following 2 variables: 
#' \itemize{
#' \item \code{X} (blood alcohol content measured in g/L) 
#' \item \code{Y} (blood alcohol content measured in g/L)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' with(data = BAC, 
#' var.test(X, Y, alternative = "less"))
#' # Convert data from wide to long format 
#' # library(reshape2)
#' # BACL <- melt(BAC, variable.name = "company", value.name = "bac")
#' # ggplot(data = BACL, aes(x = company, y = bac, fill = company)) + 
#' # geom_boxplot() + guides(fill = FALSE) + scale_fill_brewer() + 
#' # labs(y = "blood alcohol content measured in g/L")
#' # Convert with reshape()
#' BACL <- reshape(BAC, varying = c("X", "Y"), v.names = "bac", timevar = "company", 
#' direction = "long")
#' ggplot(data = BACL, aes(x = factor(company), y = bac, fill = factor(company))) + 
#' geom_boxplot() + guides(fill = FALSE) + scale_fill_brewer() + 
#' labs(y = "blood alcohol content measured in g/L", x = "") + 
#' scale_x_discrete(breaks = c(1, 2), labels = c("Company X", "Company Y"))
#' 
#' # Base graphics
#' boxplot(BAC$Y, BAC$X)
#' @keywords datasets
"BAC"
#####################################################################################
#' @name BATTERY
#' @title Lithium Batteries
#' @aliases BATTERY
#' @docType data
#' @description A manufacturer of lithium batteries has two production facilities, \code{A} and \code{B}. Facility \code{A} batteries have an advertised life of 180 hours.  Facility \code{B} batteries have an advertised life of 200 hours. Fifty randomly selected batteries from Facility \code{A} are selected and tested. Fifty randomly selected batteries from Facility \code{B} are selected and tested. The lifetimes for the tested batteries are stored in the variable \code{lifetime}.  
#' @format A data frame with 100 observations on the following two variables: 
#' \itemize{
#' \item \code{lifetime} (life time measured in hours)
#' \item \code{facility} (\code{factor} with two levels \code{A} and \code{B})  
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' p <- ggplot(data = BATTERY, aes(x = lifetime, color = facility))
#' p + geom_density()
#' q <- ggplot(data = BATTERY, aes(x = facility, y = lifetime))
#' q + geom_violin()
#' ggplot(data = BATTERY, aes(x = facility, y = lifetime, fill = facility)) + 
#' geom_violin() + scale_fill_brewer() + guides(fill = FALSE)
#' ggplot(data = BATTERY, aes(sample = lifetime)) + stat_qq() + facet_grid(. ~ facility)
#' ggplot(data = BATTERY, aes(sample = lifetime, color = facility)) + stat_qq()
#' # lattice approach
#' qqmath(~ lifetime|facility, data = BATTERY)
#' qqmath(~ lifetime, group = facility, type = c('p', 'r'), auto.key=TRUE, data = BATTERY)
#' @keywords datasets
"BATTERY"
#####################################################################################
#' @name BIOMASS
#' @title Beech Trees
#' @aliases BIOMASS
#' @docType data
#' @description Several measurements of 42 beech trees (\emph{Fagus Sylvatica}) taken from a forest in Navarre (Spain)
#' @format A data frame with 42 observations on the following 4 variables: 
#' \itemize{
#' \item \code{diameter} (diameter of the stem in centimeters)
#' \item \code{height} (height of the tree in meters)
#' \item \code{stemweight} (weight of the stem in kilograms)
#' \item \code{aboveweight} (aboveground weight in kilograms)  
#' }
#' @source \emph{Gobierno de Navarra} and \emph{Gestion Ambiental Viveros y Repoblaciones de Navarra}, 2006.  The data were obtained within the European Project FORSEE.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' pairs(BIOMASS, col = "red", cex = 0.75)
#' plot(log(aboveweight) ~ log(diameter), data = BIOMASS)
#' # logarithmic axes
#' ggplot(data = BIOMASS, aes(x = diameter, y = aboveweight, color = log(stemweight))) + 
#' geom_point() + scale_x_log10() + scale_y_log10() + 
#' labs(x = "diameter of the stem in centimeters", y = "above ground weight in kilograms")
#' @keywords datasets
"BIOMASS"
#####################################################################################
#' @name BODYFAT
#' @title Body Fat Composition
#' @aliases BODYFAT
#' @docType data
#' @description Values from a study reported in the \emph{American Journal of Clinical Nutrition} that investigated a new method for measuring body composition
#' @format A data frame with 18 observations on the following 3 variables: 
#' \itemize{
#' \item \code{age} (age in years)
#' \item \code{fat} (percent body fat composition)
#' \item \code{sex} (a factor with levels \code{F} for female and \code{M} for male) 
#' }
#' @source Mazess, R. B., Peppler, W. W., and Gibbons, M. (1984) \dQuote{Total Body Composition by Dual-Photon (153 Gd) Absorptiometry.} \emph{American Journal of Clinical Nutrition}, \bold{40}, \bold{4}: 834-839.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # base graphics
#' boxplot(fat ~ sex, data = BODYFAT)
#' # ggplot2 approach
#' ggplot(data=BODYFAT, aes(x = sex, y = fat, fill = sex)) + geom_boxplot() + 
#' labs(x = "",y = "Percent body fat") + scale_x_discrete(breaks=c("F", "M"), 
#' labels =c("Female", "Male")) + guides(fill = FALSE) + 
#' scale_fill_manual(values = c("red", "green"))
#' # Brewer Colors
#' ggplot(data=BODYFAT, aes(x = sex, y = fat, fill = sex)) + geom_boxplot() + 
#' labs(x = "", y = "Percent body fat") + scale_x_discrete(breaks=c("F", "M"), 
#' labels =c("Female", "Male")) + guides(fill = FALSE) + scale_fill_brewer()
#' ggplot(data=BODYFAT, aes(x = fat, fill = sex)) + geom_density(alpha = 0.4) + 
#' scale_fill_brewer() 
#' @keywords datasets
"BODYFAT"
#####################################################################################
#' @name CALCULUS
#' @title Calculus Assessment Scores
#' @aliases CALCULUS
#' @docType data
#' @description Mathematical assessment scores for 36 students enrolled in a biostatistics course according to whether or not the students had successfully completed a calculus course prior to enrolling in the biostatistics course
#' @format A data frame with 36 observations on the following 2 variables: 
#' \itemize{
#' \item \code{score} (assessment score for each student)
#' \item \code{calculus} (a factor with levels \code{NO} and \code{YES} for students who did not  and did successfully complete calculus prior to enrolling in the biostatistics course) 
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # ggplot2 approach
#' ggplot(data = CALCULUS, aes(sample = score)) + stat_qq() + facet_grid(. ~ calculus)
#' ggplot(data = CALCULUS, aes(x = calculus, y = score, fill = calculus)) + geom_boxplot() + 
#' guides(fill = FALSE) + scale_fill_brewer()
#' ggplot(data = CALCULUS, aes(sample = score, color = calculus)) + stat_qq()
#' # lattice approach
#' qqmath(~score|calculus, data = CALCULUS)
#' qqmath(~score, group = calculus, type = c('p', 'r'), auto.key=TRUE, data = CALCULUS)
#' @keywords datasets
"CALCULUS"
#####################################################################################
#' @name CARS2004
#' @title Cars in the European Union (2004)
#' @aliases CARS2004
#' @docType data
#' @description The numbers of cars per 1000 inhabitants (\code{cars}), the total number of known mortal accidents (\code{deaths}), and the country  population/1000 (\code{population}) for the 25 member countries of the European Union for the year 2004
#' @format A data frame with 25 observations on the following 4 variables: 
#' \itemize{
#' \item \code{country} (a factor with levels \code{Austria}, \code{Belgium}, \code{Cyprus}, \code{Czech Republic}, \code{Denmark}, \code{Estonia}, \code{Finland}, \code{France}, \code{Germany}, \code{Greece}, \code{Hungary}, \code{Ireland}, \code{Italy}, \code{Latvia}, \code{Lithuania}, \code{Luxembourg}, \code{Malta}, \code{Netherlands}, \code{Poland}, \code{Portugal}, \code{Slovakia}, \code{Slovenia}, \code{Spain}, \code{Sweden}, and \code{United Kingdom})
#' \item \code{cars} (number of cars per 1000 inhabitants)
#' \item \code{deaths} (total number of known mortal accidents) 
#' \item \code{population} (country population/1000)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' plot(deaths ~ cars, data = CARS2004)
#' ggplot(data = CARS2004, aes(x = population, y = deaths, color = cars)) + geom_point()
#' @keywords datasets
"CARS2004"
#####################################################################################
#' @name CHIPS
#' @title Silicon Chips
#' @aliases CHIPS
#' @docType data
#' @description Two techniques of splitting chips are randomly assigned to 28 sheets so that each technique is applied to 14 sheets. The the number of usable chips from each silicon sheet is stored in the variable \code{number}.
#' @format A data frame with 28 observations on the following 2 variables: 
#' \itemize{
#' \item \code{number} (number of usable chips from each silicon sheet)
#' \item \code{method} (a factor with levels \code{techniqueI} and \code{techniqueII}) 
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # ggplot2 approach
#' ggplot(data = CHIPS, aes(sample = number)) + stat_qq() + facet_grid(. ~ method)
#' ggplot(data = CHIPS, aes(sample = number, color = method)) + stat_qq()
#' ggplot(data=BODYFAT, aes(x = fat, fill = sex)) + geom_density(alpha = 0.4) + 
#' scale_fill_brewer() 
#' 
#' # lattice approach
#' qqmath(~ number|method, data = CHIPS)
#' qqmath(~ number, group = method, type = c('p', 'r'), auto.key = TRUE, data = CHIPS)
#' @keywords datasets
"CHIPS"
###########################################################################
#' @name CIRCUIT
#' @title Circuit Design Lifetime
#' @aliases CIRCUIT
#' @docType data
#' @description Results from an accelerated life test used to estimate the lifetime of four different circuit designs (lifetimes in thousands of hours)
#' @format A data frame with 26 observations on the following 2 variables: 
#' \itemize{
#' \item \code{lifetime} (lifetimes in thousands of hours)
#' \item \code{design} (a factor with levels \code{DesignI} and \code{DesignII}) 
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # ggplot2 approach
#' ggplot(data = CIRCUIT, aes(x = design, y = lifetime, fill = design)) + geom_boxplot() + 
#' labs(x = "", y = "Lifetime in thousands of hours") + guides(fill = FALSE) + 
#' scale_fill_brewer()
#' ggplot(data = CIRCUIT, aes(x = design, y = lifetime, fill = design)) + geom_violin() + 
#' labs(x = "", y = "Lifetime in thousands of hours") + guides(fill = FALSE) + 
#' scale_fill_brewer()
#' # Reorder the boxplots by medians
#' ggplot(data = CIRCUIT, aes(x = reorder(design, lifetime, FUN = median),  lifetime, 
#' fill = design)) + geom_boxplot() + labs(x = "", y = "Lifetime in thousands of hours") + 
#' guides(fill = FALSE) + scale_fill_brewer()
#' @keywords datasets
"CIRCUIT"
#####################################################################################
#' @name COSAMA
#' @title Cosmed Versus Amatek
#' @aliases COSAMA
#' @docType data
#' @description The Cosmed is a portable metabolic system. A study at Appalachian State University compared the metabolic values obtained from the Cosmed to those of a reference unit (Amatek) over a range of workloads from easy to maximal to test the validity and reliability of the Cosmed. A small portion of the results for maximal oxygen consumption (VO2 in ml/kg/min) measurements taken at a 150 watt workload are stored in \code{COSAMA}.
#' @format A data frame with 14 observations on the following 3 variables: 
#' \itemize{\item \code{subject} (subject number)
#' \item \code{cosmed} (measured VO2 with Cosmed)
#' \item \code{amatek} (measured VO2 with Amatek) 
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # ggplot2 approach
#' ggplot(data = COSAMA, aes(factor(1), y = cosmed - amatek)) + geom_boxplot() + 
#' labs(x = "")
#' # Line Plots: First change data format from wide to long with melt() from reshape2
#' # library(reshape2)
#' # CA <- melt(COSAMA, id.vars = "subject", variable.name = "treatment", 
#' # value.count = "VO2")
#' # ggplot(data = CA, aes(x = subject, y = value, color = treatment)) + geom_line()
#' # rm(CA)
#' # Convert to long format with reshape()
#' CA <- reshape(COSAMA, varying = c("cosmed", "amatek"), v.names = "VO2", 
#' timevar = "treatment", idvar = "subject", direction = "long")
#' ggplot(data = CA, aes(x = subject, y = VO2, color = factor(treatment))) + geom_line() + 
#' labs(color = "Treatment") + scale_color_discrete(labels = c("Cosmed", "Amatek"))
#' rm(CA)
#' # lattice approach
#' bwplot(~ (cosmed - amatek), data = COSAMA)
#' @keywords datasets
"COSAMA"
###########################################################################
#' @name COWS
#' @title Butterfat of Cows
#' @aliases COWS
#' @docType data
#' @description Random samples of ten mature (five-years-old and older) and ten two-year-old cows were taken from each of five breeds. The average butterfat percentage of these 100 cows is stored in the variable \code{butterfat} with the type of cow stored in the variable \code{breed} and the age of the cow stored in the variable \code{age}.
#' @format A data frame with 100 observations on the following 3 variables: 
#' \itemize{
#' \item \code{butterfat} (average butterfat percentage)
#' \item \code{age} (a factor with levels \code{2 years old} and \code{Mature})
#' \item \code{breed} (a factor with levels \code{Ayrshire}, \code{Canadian}, \code{Guernsey}, \code{Holstein-Friesian}, and \code{Jersey}) 
#' }
#' @source Canadian record book of purebred dairy cattle.
#' @references \itemize{ \item Sokal, R. R. and Rohlf, F. J. 1994. \emph{Biometry}.  W. H. Freeman, New York, third edition. \item Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.}
#' @examples
#' ggplot(data = COWS, aes(x = breed, y = butterfat, fill = age)) + 
#' geom_boxplot(position = position_dodge(1.0)) + 
#' labs(x = "", y = "Average butterfat percentage") + scale_fill_brewer()
#' summary(aov(butterfat ~ breed + age, data = COWS))
#' @keywords datasets
"COWS"
#####################################################################################
#' @name DEPEND
#' @title Number of Dependent Children for 50 Families
#' @aliases DEPEND
#' @docType data
#' @description Number of dependent children for 50 randomly selected families
#' @format A data frame with 50 observations on 1 variable: 
#' \itemize{
#' \item \code{number} (number of dependent children)
#' }
#' @source Kitchens, L. J. 2003. \emph{Basic Statistics and Data Analysis}. Duxbury.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' xtabs(~number, data = DEPEND)
#' ggplot(data = DEPEND, aes(x = factor(number))) + 
#' geom_bar(fill = "cornsilk", color = "orange") + labs(x = "Number of Dependent Children")
#' ggplot(data = DEPEND, aes(x = number)) + geom_density(fill = "pink", alpha = 0.3, 
#' color = "red") + labs(x = "Number of Dependent Children")
#' @keywords datasets
"DEPEND"
###########################################################################
#' @name DROSOPHILA
#' @title Drosophila Melanogaster
#' @aliases DROSOPHILA
#' @docType data
#' @description \code{DROSOPHILA} contains per diem fecundity (number of eggs laid per female per day for the first 14 days of life) for 25 females from each of three lines of \emph{Drosophila melanogaster}. The three lines are \code{Nonselected} (control), \code{Resistant}, and \code{Susceptible}.
#' @format A data frame with 75 observations on the following 2 variables: 
#' \itemize{
#' \item \code{fecundity} (number of eggs laid per female per day for the first 14 days of life)
#' \item \code{line} (a factor with levels \code{Nonselected}, \code{Resistant}, and \code{Susceptible})
#' }
#' @source The original measurements are from an experiment conducted by R. R. Sokal (\emph{Biometry} by Sokal and Rohlf, 1994, p. 237).
#' @references \itemize{ \item Sokal, R. R. and Rohlf, F. J. 1994. \emph{Biometry}.  W. H. Freeman, New York, third edition. \item Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.}
#' @examples
#' ggplot(data = DROSOPHILA, aes(x = reorder(line, fecundity, FUN = median),  
#' y = fecundity, fill = line)) + geom_boxplot() + guides(fill = FALSE) + 
#' labs(y ="number of eggs laid per female \n per day for the first 14 days of life", 
#' x = "") + scale_fill_brewer()
#' ggplot(data = DROSOPHILA, aes(x = reorder(line, fecundity, FUN = median), 
#' y = fecundity, fill = line)) + geom_violin() + guides(fill = FALSE) + 
#' labs(y ="number of eggs laid per female \n per day for the first 14 days of life", 
#' x = "") + scale_fill_brewer()
#' summary(aov(fecundity ~ line, data = DROSOPHILA))
#' @keywords datasets
"DROSOPHILA"
###########################################################################
#' @name ENGINEER
#' @title Engineers' Salaries
#' @aliases ENGINEER
#' @docType data
#' @description Salaries for engineering graduates 10 years after graduation
#' @format A data frame with 51 observations on the following 2 variables: 
#' \itemize{
#' \item \code{salary} (salary 10 years after graduation in thousands of dollars)
#' \item \code{university} (one of three different engineering universities)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = ENGINEER, aes(x = university, y = salary, fill = university)) + 
#' geom_boxplot() + guides(fill = FALSE) + scale_fill_brewer() + 
#' labs(y = "salary 10 years after graduation \n in thousands of dollars")
#' # Violin Plots
#' ggplot(data = ENGINEER, aes(x = university, y = salary, fill = university)) + 
#' geom_violin() + guides(fill = FALSE) + scale_fill_brewer() + 
#' labs(y = "salary 10 years after graduation \n in thousands of dollars")
#' @keywords datasets
"ENGINEER"
#####################################################################################
#' @name EPIDURAL
#' @title Traditional Sitting Position Versus Hamstring Stretch Position
#' @aliases EPIDURAL
#' @docType data
#' @description Initial results from a study to determine whether the traditional sitting position or the hamstring stretch position is superior for administering epidural anesthesia to pregnant women in labor as measured by the number of obstructive (needle to bone) contacts (\code{oc})
#' @format A data frame with 85 observations on the following 7 variables: 
#' \itemize{
#' \item \code{doctor} (a factor with levels \code{Dr. A}, \code{Dr. B}, \code{Dr. C}, and \code{Dr. D})
#' \item \code{kg} (weight in kg of patient)
#' \item \code{cm} (height in cm of patient)
#' \item \code{ease} (a factor with levels \code{Difficult}, \code{Easy}, and \code{Impossible} indicating the physicians' assessments of how well bone landmarks could be felt in the patient)
#' \item {treatment} (a factor with levels \code{Hamstring Stretch} and \code{Traditional Sitting})
#' \item \code{oc} (number of obstructive contacts)
#' \item \code{complications} (a factor with levels \code{Failure - person got dizzy}, \code{Failure - too many OCs}, \code{None}, \code{Paresthesia}, and \code{Wet Tap})
#' }
#' @source Fisher, K. S., Arnholt, A. T., Douglas, M. E., Vandiver, S. L., Nguyen, D. H. 2009. \dQuote{A Randomized Trial of the Traditional Sitting Position Versus the Hamstring Stretch Position for Labor Epidural Needle Placement.} \emph{Journal of Anesthesia & Analgesia}, Vol 109, No. 2: 532-534.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' xtabs(~ doctor + ease, data = EPIDURAL)
#' xtabs(~ doctor + factor(ease, levels = c("Easy", "Difficult", "Impossible")), 
#' data = EPIDURAL)
#' @keywords datasets
"EPIDURAL"
#####################################################################################
#' @name EPIDURALF
#' @title Traditional Sitting Position Versus Hamstring Stretch Position
#' @aliases EPIDURALF
#' @docType data
#' @description Intermediate results from a study to determine whether the traditional sitting position or the hamstring stretch position is superior for administering epidural anesthesia to pregnant women in labor as measured by the number of obstructive (needle to bone) contacts (\code{oc})
#' @format A data frame with 342 observations on the following 7 variables: 
#' \itemize{
#' \item \code{doctor} (a factor with levels \code{Dr. A}, \code{Dr. B}, \code{Dr. C}, and \code{Dr. D})
#' \item \code{kg} (weight in kg of patient)
#' \item \code{cm} (height in cm of patient)
#' \item \code{ease} (a factor with levels \code{Difficult}, \code{Easy}, and \code{Impossible} indicating the physicians' assessments of how well bone landmarks could be felt in the patient)
#' \item \code{treatment} (a factor with levels \code{Hamstring Stretch} and \code{Traditional Sitting})
#' \item \code{oc} (number of obstructive contacts)
#' \item \code{complications} (a factor with levels \code{Failure - person got dizzy}, \code{Failure - too many OCs}, \code{None}, \code{Paresthesia}, and \code{Wet Tap})
#' }
#' @source Fisher, K. S., Arnholt, A. T., Douglas, M. E., Vandiver, S. L., Nguyen, D. H. 2009. \dQuote{A Randomized Trial of the Traditional Sitting Position Versus the Hamstring Stretch Position for Labor Epidural Needle Placement.} \emph{Journal of Anesthesia & Analgesia}, Vol 109, No. 2: 532-534.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = EPIDURALF, aes(x = treatment, y = oc, fill = treatment)) +
#'  geom_boxplot() + guides(fill = FALSE) + scale_fill_brewer() + 
#'  labs(y = "number of obstructive contacts")
#' @keywords datasets
"EPIDURALF"
###########################################################################
#' @name EURD
#' @title European Union Research and Development
#' @aliases EURD
#' @docType data
#' @description A random sample of 15 countries' research and development investments for the years 2002 and 2003 was taken, and the results in millions of Euros are stored in \code{EURD}.
#' @format A data frame with 15 observations on the following 3 variables: 
#' \itemize{
#' \item \code{country} (a character vector with values \code{Bulgaria}, \code{Croatia}, \code{Cyprus}, \code{Czech Republic}, \code{Estonia}, \code{France}, \code{Hungary}, \code{Latvia}, \code{Lithuania}, \code{Malta}, \code{Portugal}, \code{Romania}, \code{Slovakia}, and \code{Slovenia})
#' \item \code{rd2002} (research and development investments in millions of Euros for 2002)
#' \item \code{rd2003} (research and development investments in millions of Euros for 2003)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = EURD, aes(x = rd2002, y =  rd2003)) + geom_point() +
#' geom_smooth(method = "lm")
#' ggplot(data = EURD, aes(sample = rd2003 - rd2002)) + stat_qq()
#' # lattice approach
#' qqmath(~ (rd2003 - rd2002), data = EURD, type =c("p", "r"))
#' @keywords datasets
"EURD"
###########################################################################
#' @name FAGUS
#' @title Retained Carbon in Beech Trees
#' @aliases FAGUS
#' @docType data
#' @description The carbon retained by leaves measured in kg/ha is recorded for forty-one different plots of mountainous regions of Navarre (Spain), depending on the forest classification: areas with 90\% or more beech trees (\emph{Fagus Sylvatica}) are labeled monospecific, while areas with many species of trees are labeled multispecific.
#' @format A data frame with 41 observations on the following 3 variables: 
#' \itemize{
#' \item \code{plot} (plot number)
#' \item \code{carbon} (carbon retained by leaves measured in kg/ha)
#' \item \code{type} (a factor with levels \code{monospecific} and \code{multispecific})
#' }
#' @source \emph{Gobierno de Navarra} and \emph{Gestion Ambiental Viveros y Repoblaciones de Navarra}, 2006.  The data were obtained within the European Project FORSEE.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = FAGUS, aes(x = type, y = carbon)) + geom_boxplot()
#' @keywords datasets
"FAGUS"
###########################################################################
#' @name FCD
#' @title Fat Cats
#' @aliases FCD
#' @docType data
#' @description In a weight loss study on obese cats, overweight cats were randomly assigned to one of three groups and boarded in a kennel.  In each of the three groups, the cats' total caloric intake was strictly controlled (1 cup of generic cat food) and monitored for 10 days. The difference between the groups was that group A was given 1/4 of a cup of cat food every six hours, group B was given 1/3 a cup of cat food every eight hours, and group C was given 1/2 a cup of cat food every twelve hours. The weights of the cats at the beginning and end of the study were recorded, and the difference in weights (grams) was stored in the variable \code{Weight} of the data frame \code{FCD}.
#' @format A data frame with 36 observations on the following 2 variables: 
#' \itemize{
#' \item \code{weight} (difference in weight (grams))
#' \item \code{diet} (a factor with levels \code{A}, \code{B}, and \code{C})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # checking.plots()?
#' p <- ggplot(data = FCD, aes(x = diet, y = weight))
#' p + geom_violin(fill = "blue")
#' aov(weight ~ diet, data = FCD)
#' @keywords datasets
"FCD"
###########################################################################
#' @name FERTILIZE
#' @title Cross and Auto Fertilization
#' @aliases FERTILIZE
#' @docType data
#' @description Plants' heights in inches obtained from two seeds, one obtained by cross fertilization and the other by auto fertilization, in two opposite but separate locations of a pot are recorded.
#' @format A data frame with 30 observations on the following 3 variables: 
#' \itemize{
#' \item \code{height} (height of plant in inches)
#' \item \code{fertilization} (a factor with levels \code{cross} and \code{self})
#' \item \code{pot} (a factor with fifteen levels)
#' }
#' @source Darwin, C. 1876. \emph{The Effect of Cross and Self-Fertilization in the Vegetable Kingdom.} D. Appleton and Company. 
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' p <- ggplot(data = FERTILIZE, aes(x = height, color = fertilization))
#' p + geom_density()
#' t.test(height ~ fertilization, data = FERTILIZE)
#' @keywords datasets
"FERTILIZE"
###########################################################################
#' @name FOOD
#' @title Carrot Shear
#' @aliases FOOD
#' @docType data
#' @description Shear measured in kN on frozen carrots from four randomly selected freezers
#' @format A data frame with 16 observations on the following 2 variables: 
#' \itemize{
#' \item \code{shear} (carrot shear measured in kN)
#' \item \code{freezer} (a factor with levels \code{A}, \code{B}, \code{C}, and \code{D})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' summary(aov(shear ~ freezer, data = FOOD))
#' @keywords datasets
"FOOD"
###########################################################################
#' @name FORMULA1
#' @title Pit Stop Times
#' @aliases FORMULA1
#' @docType data
#' @description Pit stop times for two teams at 10 randomly selected Formula 1 races
#' @format A data frame with 10 observations on the following 3 variables: 
#' \itemize{
#' \item \code{race} (number corresponding to a race site)
#' \item \code{team1} (pit stop times for team one)
#' \item \code{team2} (pit stop times for team two)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # Change data format from wide to long
#' # library(reshape2)
#' # F1L <- melt(data = FORMULA1, id.vars = "race", variable.name = "team", 
#' # value.name = "time")
#' # ggplot(data = F1L, aes(x = team, y = time)) + geom_boxplot()
#' # Using reshape()
#' F1L <- reshape(FORMULA1, varying = c("team1", "team2"), v.names = "time", 
#' timevar = "team", idvar = "race", direction = "long")
#' ggplot(data = F1L, aes(x = factor(team), y = time, fill = factor(team))) + 
#' geom_boxplot() + guides(fill = FALSE) + scale_x_discrete(breaks = 1:2, 
#' labels = c("Team 1", "Team 2")) + labs(x = "", y = "Pit stop times in seconds")
#' with(data = FORMULA1, 
#' boxplot(team1, team2, col = c("red", "blue")))
#' @keywords datasets
"FORMULA1"
###########################################################################
#' @name GD
#' @title Times Until Failure
#' @aliases GD
#' @docType data
#' @description Contains time until failure in hours for a particular electronic component subjected to an accelerated stress test
#' @format A data frame with 100 observations on the following variable: 
#' \itemize{
#' \item \code{attf} (times until failure in hours)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = GD, aes(x = attf, y = ..density..)) + 
#' geom_histogram(binwidth = 2, fill = "cornsilk", color = "orange") + 
#' geom_density(color = "gray", size = 1) + labs( x = "time until failure in hours")
#' @keywords datasets
"GD"
###########################################################################
#' @name GLUCOSE
#' @title Blood Glucose Levels
#' @aliases GLUCOSE
#' @docType data
#' @description Fifteen diabetic patients were randomly selected, and their blood glucose levels were measured in mg/100 ml with two different devices.
#' @format A data frame with 15 observations on the following 3 variables: 
#' \itemize{
#' \item \code{patient} (patient number)
#' \item \code{old} (blood glucose level in mg/100 ml using an old device)
#' \item \code{new} (blood glucose level in mg/100 ml using a new device)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' with(data = GLUCOSE,
#' boxplot(old, new, col = c("red", "blue")))
#' @keywords datasets
"GLUCOSE"
###########################################################################
#' @name GRADES
#' @title GPA and SAT Scores
#' @aliases GRADES
#' @docType data
#' @description The admissions committee of a comprehensive state university selected, at random, the records of 200 second semester freshmen. The results, first semester college GPA and high school SAT scores, are stored in the data frame \code{GRADES}.
#' @format A data frame with 200 observations on the following 2 variables: 
#' \itemize{
#' \item \code{sat} (SAT score)
#' \item \code{gpa} (grade point average)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # base scatterplot
#' plot(gpa ~ sat, data = GRADES)
#' # lattice scatterplot
#' xyplot(gpa ~ sat, data = GRADES, type = c("p", "smooth"))
#' # ggplot scatterplot
#' ggplot(data = GRADES, aes(x = sat, y = gpa)) + geom_point() + geom_smooth()
#' @keywords datasets
"GRADES"
###########################################################################
#' @name GROCERY
#' @title Grocery Spending
#' @aliases GROCERY
#' @docType data
#' @description The consumer expenditure survey, created by the U.S. Department of Labor, was administered to 30 households in Watauga County, North Carolina, to see how the cost of living in Watauga county with respect to total dollars spent on groceries compares with other counties. The amount of money each household spent per week on groceries is stored in the variable \code{amount}. 
#' @format A data frame with 30 observations on the following variable: 
#' \itemize{
#' \item \code{amount} (total dollars spent on groceries)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' with(data = GROCERY, 
#' z.test(amount, sigma.x = 25, mu = 100, alternative = "greater"))
#' hist(GROCERY$amount, xlab = "Weekly grocery bill", main = "")
#' ggplot(data = GROCERY, aes(x = amount, y = ..density..)) + 
#' geom_histogram(binwidth = 8, fill = "cornsilk", color = "gray80") + 
#' geom_density(color = "lightblue", size = 1, fill = "lightblue", alpha = .2) + 
#' labs(x = "Weekly grocery bill (in dollars)")
#' @keywords datasets
"GROCERY"
###########################################################################
#' @name HARDWATER
#' @title Mortality and Water Hardness
#' @aliases HARDWATER
#' @docType data
#' @description Mortality and drinking water hardness for 61 cities in England and Wales
#' @format A data frame with 61 observations on the following 4 variables: 
#' \itemize{
#' \item \code{location} (a factor with levels \code{North} and \code{South} indicating whether the town is as far north as Derby or further)
#' \item \code{town} (the name of the town)
#' \item \code{mortality} (average annual mortality per 100,000 males)
#' \item \code{hardness} (calcium concentration (in parts per million))
#' }
#' @source D. J. Hand, F. Daly, A. D. Lunn, K. J. McConway and E. Ostrowski. 1994. \emph{A Handbook of Small Datasets}. Chapman and Hall/CRC, London.
#' @details These data were collected in an investigation of environmental causes of disease.  They show the annual mortality rate per 100,000 for males, averaged over the years 1958-1964, and the calcium concentration (in parts per million) in the drinking water supply for 61 large towns in England and Wales.  (The higher the calcium concentration, the harder the water.)
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = HARDWATER, aes(x = hardness, y = mortality, color  = location)) + 
#' geom_point() + labs(y = "averaged annual mortality per 100,000 males", 
#' x = "calcium concentration (in parts per million)")
#' @keywords datasets
"HARDWATER"
###########################################################################
#' @name HOUSE
#' @title House Prices
#' @aliases HOUSE
#' @docType data
#' @description Random sample of house prices (in thousands of dollars) for three bedroom/two bath houses in Watauga County, NC
#' @format A data frame with 14 observations on the following 2 variables: 
#' \itemize{
#' \item \code{neighborhood} (a factor with levels \code{Blowing Rock}, \code{Cove Creek}, \code{Green Valley}, \code{Park Valley}, \code{Parkway}, and \code{Valley Crucis})
#' \item \code{price} (price of house in thousands of dollars)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' with(data = HOUSE,
#' t.test(price, mu = 225))
#' @keywords datasets
"HOUSE"
###########################################################################
#' @name HSWRESTLER
#' @title High School Wrestlers
#' @aliases HSWRESTLER
#' @docType data
#' 
#' @description The body fat percentage of 78 high school wrestlers was measured using three separate techniques, and the results are stored in the data frame \code{HSWRESTLER}. The techniques used were hydrostatic weighing (\code{hwfat}), skin fold measurements (\code{skfat}), and the Tanita body fat scale (\code{tanfat}).
#' 
#' @format A data frame with 78 observations on the following 9 variables: 
#' \itemize{
#' \item \code{age} (age of wrestler in years)
#' \item \code{ht} (height of wrestler in inches)
#' \item \code{wt} (weight of wrestler in pounds)
#' \item \code{abs} (abdominal fat)
#' \item \code{triceps} (tricep fat)
#' \item \code{subscap} (subscapular fat)
#' \item \code{hwfat} (hydrostatic measure of percent fat)
#' \item \code{tanfat} (Tanita measure of percent fat)
#' \item \code{skfat} (skin fold measure of percent fat)
#' }
#' @source Data provided by Dr. Alan Utter, Department of Health Leisure and Exercise Science, Appalachian State University
#' 
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' 
#' @examples
#' FAT <- c(HSWRESTLER$hwfat, HSWRESTLER$tanfat, HSWRESTLER$skfat)
#' GROUP <- factor(rep(c("hwfat", "tanfat", "skfat"), rep(78, 3)))
#' BLOCK <- factor(rep(1:78, 3))
#' friedman.test(FAT ~ GROUP | BLOCK)
#' rm(FAT, BLOCK, GROUP)
#' ggplot(data = HSWRESTLER, aes(x = tanfat, y = hwfat, color = age)) + geom_point() + 
#' geom_smooth() + labs(x = "Tanita measure of percent fat", 
#' y = "hydrostatic measure of percent fat")
#' 
#' @keywords datasets
"HSWRESTLER"
###########################################################################
#' @name HUBBLE
#' @title Hubble Telescope
#' @aliases HUBBLE
#' @docType data
#' @description The Hubble Space Telescope was put into orbit on April 25, 1990. Unfortunately, on June 25, 1990, a spherical aberration was discovered in Hubble's primary mirror. To correct this, astronauts had to work in space. To prepare for the mission, two teams of astronauts practiced  making repairs under simulated space conditions. Each team of astronauts went through 15 identical scenarios. The times to complete each scenario were recorded in days.
#' @format A data frame with 15 observations on the following 2 variables: 
#' \itemize{
#' \item \code{team1} (days to complete scenario)
#' \item \code{team2} (days to complete scenario)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' with(data = HUBBLE, 
#' qqnorm(team1 - team2))
#' with(data = HUBBLE, 
#' qqline(team1 - team2))
#' # Trellis Approach
#' qqmath(~(team1 - team2), data = HUBBLE, type=c("p", "r"))
#' # ggplot approach
#' ggplot(data = HUBBLE, aes(sample = team1 - team2)) + stat_qq(color = "blue")
#' @keywords datasets
"HUBBLE"
###########################################################################
#' @name INSURQUOTES
#' @title Insurance Quotes
#' @aliases INSURQUOTES
#' @docType data
#' @description Insurance quotes for two insurers of hazardous waste jobs
#' @format A data frame with 15 observations on the following 2 variables: 
#' \itemize{
#' \item \code{companyA} (quotes from company A in Euros)
#' \item \code{companyB} (quotes from company B in Euros)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = INSURQUOTES, aes(sample = companyA - companyB)) + 
#' stat_qq(col = "orange", size = 4)
#' with(data = INSURQUOTES,
#' t.test(companyA, companyB))
#' @keywords datasets
"INSURQUOTES"
###########################################################################
#' @name JANKA
#' @title Australian Eucalypt Hardwoods
#' @aliases JANKA
#' @docType data
#' @description The dataset consists of density and hardness measurements from 36 Australian Eucalypt hardwoods.
#' @format A data frame with 36 observations on the following 2 variables: 
#' \itemize{
#' \item \code{density} (a measure of density of the timber)
#' \item \code{hardness} (the Janka hardness of the timber)
#' }
#' @details Janka hardness is a structural property of Australian hardwood timbers. The Janka hardness test measures the force required to imbed a steel ball into a piece of wood.
#' @source Williams, E.J. 1959. \emph{Regression Analysis}.  John Wiley & Sons, New York.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = JANKA, aes(x = density, y = hardness)) + geom_point() + geom_smooth()
#' @keywords datasets
"JANKA"
###########################################################################
#' @name KINDER
#' @title Kindergarten Class
#' @aliases KINDER
#' @docType data
#' @description The data frame \code{KINDER} contains the height in inches and weight in pounds of 20 children from a kindergarten class.
#' @format A data frame with 20 observations on the following 2 variables: 
#' \itemize{
#' \item \code{ht} (height in inches of each child)
#' \item \code{wt} (weight in pounds of each child)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = KINDER, aes(x = ht, y = wt)) + geom_point(color = "blue") + 
#' geom_smooth(method = "lm", color = "red") + labs(x = "height in inches", 
#' y = "weight in pounds")
#' @keywords datasets
"KINDER"
###########################################################################
#' @name LEDDIODE
#' @title LED Diodes
#' @aliases LEDDIODE
#' @docType data
#' @description The diameter in millimeters for a random sample of 15 diodes from each of the two suppliers is stored in the data frame \code{LEDDIODE}.
#' @format A data frame with 30 observations on the following 2 variables: 
#' \itemize{
#' \item \code{diameter} (diameter of diode measured in millimeters)
#' \item \code{supplier} (factor with levels \code{supplierA} and \code{supplierB})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = LEDDIODE, aes(supplier, diameter)) + geom_boxplot()
#' @keywords datasets
"LEDDIODE"
###########################################################################
#' @name LOSTR
#' @title Lost Revenue Due to Worker Illness
#' @aliases LOSTR
#' @docType data
#' @description Data set containing the lost revenue in dollars/day and number of workers absent due to illness for a metallurgic company
#' @format A data frame with 25 observations on the following 2 variables: 
#' \itemize{
#' \item \code{numbersick} (number of absent workers due to illness)
#' \item \code{lostrevenue} (lost revenue in dollars)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = LOSTR, aes(x = numbersick, y = lostrevenue)) + geom_point(color = "red", 
#' pch = 21, fill = "pink", size = 4) + geom_smooth(method = "lm") + 
#' labs(x = "number of absent workers due to illness", y = "lost revenue in dollars")
#' @keywords datasets
"LOSTR"
###########################################################################
#' @name MILKCARTON
#' @title Milk Carton Drying Times
#' @aliases MILKCARTON
#' @docType data
#' @description A plastics manufacturer makes two sizes of milk containers: half gallon and gallon sizes. The time required for each size to dry is recorded in seconds in the data frame \code{MILKCARTON}.
#' @format A data frame with 80 observations on the following 2 variables: 
#' \itemize{
#' \item \code{seconds} (drying time in seconds)
#' \item \code{size} (factor with levels \code{halfgallon} and \code{wholegallon})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = MILKCARTON, aes(x = size, y = seconds)) + geom_boxplot()
#' ggplot(data = MILKCARTON, aes(x = size, y = seconds, fill = size)) + geom_boxplot() + 
#' guides(fill = FALSE) + scale_fill_brewer() + 
#' labs(x = "size of container", y = "drying time in seconds")
#' @keywords datasets
"MILKCARTON"
#####################################################################################
#' @name NC2010DMG
#' @title North Carolina Demographics
#' @aliases NC2010DMG
#' @docType data
#' @description North Carolina county demographics for 2010 and county voter information for the North Carolina Amendment 1 ballot initiative which took place May 8, 2012, are stored in the data frame \code{NC2010DMG}.
#' @format A data frame with 100 observations (counties) on the following 32 variables: 
#' \itemize{
#' \item \code{countyName} (Name of North Carolina county)
#' \item \code{pop2010} (Total population of the county in 2010)
#' \item \code{medage} (Median age of the county in 2010)
#' \item \code{divorced} (Number of divorced adults in 2010)
#' \item \code{pctrural} (The percent of the population that lived in a rural area of the county in 2010)
#' \item \code{edu_baorup} (The total number of people with a Bachelor's degree in 2010)
#' \item \code{medinc} (The median household income adjusted for inflation in 2010)
#' \item \code{col_enroll} (The number of people enrolled in college in 2010)
#' \item \code{age18-24} (The number of people between the ages of 18 and 24 in the county in 2010)
#' \item \code{age25-29} (The number of people between the ages of 25 and 29 in the county in 2010)
#' \item \code{age60up} (The number of people over the age of 60 in the county in 2010)
#' \item \code{white} (The number of white people in the county in 2010)
#' \item \code{black} (The number of black people in the county in 2010)
#' \item \code{MaleBachelor} (The number of males with a Bachelor's degree in 2010)
#' \item \code{MaleMaster} (The number of males with a Master's degree in 2010)
#' \item \code{MaleProfessional} (The number of males with a professional degree in 2010)
#' \item \code{MaleDoctorate} (The number of males with a Doctorate degree in 2010)
#' \item \code{FemaleBachelor} (The number of females with a Bachelor's degree in 2010)
#' \item \code{FemaleMaster} (The number of females with a Master's degree in 2010)
#' \item \code{FemaleProfessional} (The number of females with a professional degree in 2010)
#' \item \code{FemaleDoctorate} (The number of females with a Doctorate degree in 2010)
#' \item \code{Owneroccupied} (The number of homes that are owner occupied in 2010)
#' \item \code{Renteroccupied} (The number of homes that are renter occupied in 2010)
#' \item \code{popden} (The number of people per square mile in 2010)
#' \item \code{pctfor} (The percent of voters that voted for Amendment 1 on May 8, 2012)
#' \item \code{turnout} (The percent of registered voters who voted May 8, 2012)
#' \item \code{obama08} (The percent of voters who voted for Barrack Obama in the 2008 presidential election)
#' \item \code{mccain08} (The percent of voters who voted for John McCain in the 2008 presidential election)
#' \item \code{evanrate} (Evangelical rates of adherence per 1,000 population in 2010)
#' \item \code{churches} (The number of churches in the county in 2010)
#' \item \code{colleges} (The number of colleges in the county in 2010)
#' }
#' @source The original data was provided by E.L. Davison, Department of Sociology, Appalachian State University.  Variables
#' \code{countyName} through \code{popden} were obtained from \url{http://factfinder2.census.gov/faces/nav/jsf/pages/searchresults.xhtml?refresh=t#"} and further cleaned by Maureen O'Donnell and Eitan Lees. The variables \code{pctfor} through \code{mccain08} were obtained from \url{http://www.ncsbe.gov/}. The variables \code{evanrate} and \code{churches} were obtained from \url{http://thearda.com}, while the information for \code{colleges} was obtained from \url{http://collegestats.org/colleges/north-carolina}.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = MILKCARTON, aes(x = size, y = seconds)) + geom_boxplot()
#' ggplot(data = MILKCARTON, aes(x = size, y = seconds, fill = size)) + geom_boxplot() + 
#' guides(fill = FALSE) + scale_fill_brewer() + 
#' labs(x = "size of container", y = "drying time in seconds")
#' @keywords datasets
"NC2010DMG"
#####################################################################################
#' @name PAMTEMP
#' @title Pamplona Temperatures
#' @aliases PAMTEMP
#' @docType data
#' @description The data frame \code{PAMTEMP} has records of the temperature and precipitation for Pamplona, Spain from January 1, 1990 to December 31, 2010. 
#' @format A data frame with 7547 observations on the following 7 variables: 
#' \itemize{
#' \item \code{tmax} (maximum daily temperature in Celsius)
#' \item \code{tmin} (minimum daily temperature in Celsius)
#' \item \code{precip} (daily precipitation in mm)
#' \item \code{day} (day of the month)
#' \item \code{month} (month of the year)
#' \item \code{year} (year)
#' \item \code{tmean} (the average of \code{tmax} and \code{tmin})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' str(PAMTEMP)
#' levels(PAMTEMP$month)
#' PAMTEMP$month <- factor(PAMTEMP$month, levels = month.abb[1:12])
#' levels(PAMTEMP$month)
#' ggplot(data = PAMTEMP, aes(x = 1:dim(PAMTEMP)[1], y = tmean)) + 
#' geom_line() + 
#' theme_bw() + 
#' labs(x = "", y = "Average Temperature (Celcius)")
#' @keywords datasets
"PAMTEMP"
#####################################################################################
#' @name PHENYL
#' @title Phenylketonuria
#' @aliases PHENYL
#' @docType data
#' @description The data frame \code{PHENYL} records the level of Q10 at four different times for 46 patients diagnosed with phenylketonuria. The variable \code{Q10.1} contains the level of Q10 measured in micromoles for the 46 patients. \code{Q10.2}, \code{Q10.3}, and \code{Q10.4} are the values recorded at later times, respectively, for the 46 patients.
#' @format A data frame with 46 observations on the following 4 variables: 
#' \itemize{
#' \item \code{Q10.1} (level of Q10 at time 1 in micromoles)
#' \item \code{Q10.2} (level of Q10 at time 2 in micromoles)
#' \item \code{Q10.3} (level of Q10 at time 3 in micromoles)
#' \item \code{Q10.4} (level of Q10 at time 4 in micromoles)
#' }
#' @details Phenylketonuria (PKU) is a genetic disorder that is characterized by an inability of the body to utilize the essential amino acid, phenylalanine. Research suggests patients with phenylketonuria have deficiencies in coenzyme Q10.
#' @source Artuch, R., \emph{et. al.} 2004. \dQuote{Study of Antioxidant Status in Phenylketonuric Patients.} \emph{Clinical Biochemistry}, \bold{37}: 198-203.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' PL <- stack(PHENYL)
#' PL$sub <- factor(rep(1:46, 4))
#' ggplot(data = PL, aes(x= ind, y = values, group = sub, color = sub)) + geom_line() + 
#' guides(color = FALSE)
#' with(data = PHENYL,
#' t.test(Q10.1, conf.level = 0.99))
#' @keywords datasets
"PHENYL"
###########################################################################
#' @name PHONE
#' @title Telephone Call Times
#' @aliases PHONE
#' @docType data
#' @description \code{PHONE} contains times in minutes of long distance telephone calls during a one month period for a small business.
#' @format A data frame with 23 observations on the following variable: 
#' \itemize{
#' \item \code{call.time} (time spent on long distance calls in minutes)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' with(data = PHONE, 
#' SIGN.test(call.time, md = 2.1))
#' @keywords datasets
"PHONE"
###########################################################################
#' @name RAT
#' @title Rat Survival Time
#' @aliases RAT
#' @docType data
#' @description The survival time in weeks of 20 male rats exposed to high levels of radiation
#' @format A data frame with 20 observations on the following variable: 
#' \itemize{
#' \item \code{survival.time} (number of weeks survived)
#' }
#' @source Lawless, J. 1982. \emph{Statistical Models and Methods for Lifetime Data}. John Wiley, New York.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = RAT, aes(sample = survival.time)) + stat_qq()
#' ggplot(data = RAT, aes(x = survival.time)) + geom_density(alpha = 0.2, fill = "blue") +
#' labs(x = "Survival time in weeks")
#' @keywords datasets
"RAT"
###########################################################################
#' @name RATBP
#' @title Rat Blood Pressure
#' @aliases RATBP
#' @docType data
#' @description Twelve rats were chosen, and a drug was administered to six rats, the treatment group, chosen at random. The other six rats, the control group, received a placebo. The drops in blood pressure (mmHg) for the treatment group (with probability distribution F) and the control group (with probability distribution G), respectively, were recorded.
#' @format A data frame with 12 observations on the following 2 variables: 
#' \itemize{
#' \item \code{mmHg} (drops in blood pressure in mm of Hg where positive values are decreases, negative values are increases)
#' \item \code{group} (factor with levels \code{control} and \code{treatment})
#' }
#' @source The data is originally from Ott and Mendenhall (\emph{Understanding Statistics}, 1985, problem 8.17).
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # Boxplot
#' ggplot(data = RATBP, aes(x = group, y = mmHg)) + geom_boxplot()
#' ggplot(data = RATBP, aes(x = group, y = mmHg, fill = group)) + geom_boxplot() + 
#' guides(fill = FALSE) + labs(x = "", y = "drops in blood pressure in mm  of Hg") + 
#' scale_fill_brewer()
#' @keywords datasets
"RATBP"
###########################################################################
#' @name REFRIGERATOR
#' @title Refrigerator Energy Consumption
#' @aliases REFRIGERATOR
#' @docType data
#' @description Sixty 18 cubic feet refrigerators were randomly selected from a company's warehouse. The first thirty had their motors modified while the last thirty were left intact. The energy consumption (kilowatts) for a 24 hour period for each refrigerator was recorded and stored in the variable \code{kilowatts}.
#' @format A data frame with 60 observations on the following 2 variables: 
#' \itemize{
#' \item \code{kilowatts} (energy consumption in kilowatts for a 24 hour period)
#' \item \code{group} (factor with levels \code{original} and \code{modified})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # Boxplot
#' ggplot(data = REFRIGERATOR, aes(x = group, y = kilowatts)) + geom_boxplot()
#' ggplot(data = REFRIGERATOR, aes(x = group, y = kilowatts, fill = group)) + 
#' geom_boxplot() + labs(y = "energy consumption in kilowatts for a 24 hour period") + 
#' guides(fill = FALSE) + scale_fill_brewer()
#' @keywords datasets
"REFRIGERATOR"
###########################################################################
#' @name ROACHEGGS
#' @title Oriental Cockroaches
#' @aliases ROACHEGGS
#' @docType data
#' @description A laboratory is interested in testing a new child friendly pesticide on \emph{Blatta orientalis} (oriental cockroaches). Scientists apply the new pesticide to 81 randomly selected Blatta orientalis oothecae (eggs). The results from the experiment are stored in the data frame \code{ROACHEGGS} in the variable \code{eggs}.  A zero in the variable \code{eggs} indicates that nothing hatched from the egg while a 1 indicates the birth of a cockroach.
#' @format A data frame with 81 observations on the following variable: 
#' \itemize{
#' \item \code{eggs} (numeric vector where a 0 indicates nothing hatched while a 1 indicates the birth of a cockroach.)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' p <- seq(0.1, 0.9, 0.001) 
#' negloglike <- function(p){ 
#' -(sum(ROACHEGGS$eggs)*log(p) + sum(1-ROACHEGGS$eggs)*log(1-p))
#' }
#' nlm(negloglike, .2)
#' rm(p, negloglike)
#' @keywords datasets
"ROACHEGGS"
###########################################################################
#' @name SALINITY
#' @title Surface-Water Salinity
#' @aliases SALINITY
#' @docType data
#' @description Surface-water salinity measurements were taken in a bottom-sampling project in Whitewater Bay, Florida. 
#' @format A data frame with 48 observations on the following variable: 
#' \itemize{
#' \item \code{salinity} (surface-water salinity measurements)
#' }
#' @source Davis, J. 1986. \emph{Statistics and Data Analysis in Geology}. John Wiley, New York.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # Boxplot
#' ggplot(data = SALINITY, aes(x = salinity)) + geom_density(fill = "yellow", alpha = 0.3)
#' @keywords datasets
"SALINITY"
###########################################################################
#' @name SATFRUIT
#' @title Fruit Trees
#' @aliases SATFRUIT
#' @docType data
#' @description To estimate the total surface occupied by fruit trees in 3 small areas (R63, R67, and R68) of Navarre (Spain) in 2001, a sample of 47 square segments has been taken. The experimental units are square segments or quadrats of 4 hectares, obtained by random sampling after overlaying a square grid on the study domain. 
#' @format A data frame with 47 observations on the following 17 variables: 
#' \itemize{
#' \item \code{quadrat} (number of the sampled segment or quadrat)
#' \item \code{smallarea} (the small area, a factor with levels \code{R63}, \code{R67}, and \code{R68})
#' \item \code{wheat} (area classified as wheat in the sampled segment)
#' \item \code{barley} (area classified as barley in the sampled segment)
#' \item \code{nonarable} (area classified as non-arable in the sampled segment)
#' \item \code{corn} (area classified as corn in the sampled segment)
#' \item \code{sunflower} (area classified as sunflower in the sampled segment)
#' \item \code{vineyard} (area classified as vineyard in the sampled segment)
#' \item \code{grass} (area classified as grass in the sampled segment)
#' \item \code{asparagus} (area classified as asparagus in the sampled segment)
#' \item \code{alfalfa} (area classified as alfalfa in the sampled segment)
#' \item \code{rape} (area classified as rape in the sampled segment)
#' \item \code{rice} (area classified as rice in the sampled segment)
#' \item \code{almonds} (area classified as almonds in the sampled segment)
#' \item \code{olives} (area classified as olives in the sampled segment)
#' \item \code{fruit} (area classified as fruit trees in the sampled segment)
#' \item \code{observed} (the observed area of fruit trees in the sampled segment)
#' }
#' @source Militino, A. F., \emph{et. al.} 2006. \dQuote{Using Small Area Models to Estimate the Total Area Occupied by Olive Trees.} \emph{Journal of Agricultural, Biological and Environmental Statistics}, \bold{11}: 450-461.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' pairs(SATFRUIT[,15:17])
#' @keywords datasets
"SATFRUIT"
###########################################################################
#' @name SBIQ
#' @title County IQ
#' @aliases SBIQ
#' @docType data
#' @description A school psychologist administered the Stanford-Binet intelligence quotient (IQ) test in two counties. Forty randomly selected, gifted and talented students were selected from each county. The Stanford-Binet IQ test is said to follow a normal distribution with a mean of 100 and standard deviation of 16.
#' @format A data frame with 80 observations on the following 2 variables: 
#' \itemize{
#' \item \code{score} (IQ score)
#' \item \code{county} (factor with levels \code{County1} and \code{County2})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = SBIQ, aes(sample = score, color = county)) + stat_qq()
#' @keywords datasets
"SBIQ"
###########################################################################
#' @name SCHIZO
#' @title Dopamine Activity
#' @aliases SCHIZO
#' @docType data
#' @description Twenty-five patients with schizophrenia were classified as psychotic or nonpsychotic after being treated with an antipsychotic drug.  Samples of cerebral fluid were taken from each patient and assayed for dopamine \eqn{\beta}-hydroxylase (DBH) activity. The dopamine measurements for the two groups are in nmol/ml-hour per milligram of protein.
#' @format A data frame with 25 observations on the following 2 variables: 
#' \itemize{
#' \item \code{dopamine} (dopamine activity level)
#' \item \code{classification} (factor with levels \code{psychotic} and \code{nonpsychotic})
#' }
#' @source Sternberg, D. E., Van Kammen, D. P., and Bunney,W. E. 1982. \dQuote{Schizophrenia: Dopamine \eqn{\beta}-Hydroxylase Activity and Treatment Response.} \emph{Science}, \bold{216}: 1423-1425.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = SCHIZO, aes(x = classification, y = dopamine)) + geom_boxplot()
#' @keywords datasets
"SCHIZO"
###########################################################################
#' @name SCORE
#' @title Standardized Test Scores
#' @aliases SCORE
#' @docType data
#' @description Standardized test scores from a random sample of twenty college freshmen
#' @format A data frame with 20 observations on the following  variable: 
#' \itemize{
#' \item \code{scores} (standardized test score)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = SCORE, aes(sample = scores)) + stat_qq()
#' @keywords datasets
"SCORE"
###########################################################################
#' @name SDS4
#' @title M1 Motorspeedway Times
#' @aliases SDS4
#' @docType data
#' @description The times recorded are those for 41 successive vehicles travelling northwards along the M1 motorway in England when passing a fixed point near Junction 13 in Bedfordshire on Saturday, March 23, 1985. After subtracting the times, the following 40 interarrival times reported to the nearest second are stored in \code{SDS4} under the variable \code{times}.
#' @format A data frame with 40 observations on the following variable: 
#' \itemize{
#' \item \code{times} (interarrival times to the nearest second)
#' }
#' @source Hand, D. J., \emph{et. al.} 1994. \emph{A Handbook of Small Data Sets}. Chapman & Hall, London.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = SDS4, aes(x = times)) + geom_histogram(binwidth = 2)
#' ggplot(data = SDS4, aes(x = times, y = ..density..)) + 
#' geom_histogram(binwidth = 2, color = "red", fill = "pink", alpha = 0.5) + 
#' geom_density(fill = "cornsilk", alpha = 0.5) + 
#' labs(x = "interarrival times to the nearest second", y = "")
#' @keywords datasets
"SDS4"
###########################################################################
#' @name SIMDATAST
#' @title Simulated Data (Predictors)
#' @aliases SIMDATAST
#' @docType data
#' @description Simulated data for five variables
#' @format A data frame with 200 observations on the following 5 variables: 
#' \itemize{
#' \item \code{y1} (a numeric vector)
#' \item \code{y2} (a numeric vector)
#' \item \code{x1} (a numeric vector)
#' \item \code{x2} (a numeric vector)
#' \item \code{x3} (a numeric vector)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = SIMDATAST, aes(x = x1, y = y1)) + geom_point() + geom_smooth()
#' @keywords datasets
"SIMDATAST"
###########################################################################
#' @name SIMDATAXT
#' @title Simulated Data (Logarithms)
#' @aliases SIMDATAXT
#' @docType data
#' @description Simulated data for four variables
#' @format A data frame with 200 observations on the following 4 variables: 
#' \itemize{
#' \item \code{y1} (a numeric vector)
#' \item \code{y2} (a numeric vector)
#' \item \code{x1} (a numeric vector)
#' \item \code{x2} (a numeric vector)
#' \item \code{x3} (a numeric vector)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = SIMDATAXT, aes(x = x1, y = y)) + geom_point() + geom_smooth()
#' @keywords datasets
"SIMDATAXT"
###########################################################################
#' @name SOCCER
#' @title World Cup  Soccer
#' @aliases SOCCER
#' @docType data
#' @description \code{SOCCER} contains how many goals were scored in the regulation 90 minute periods of World Cup soccer matches from 1990 to 2002.
#' @format A data frame with 575 observations on the following 3 variables: 
#' \itemize{
#' \item \code{cgt} (cumulative goal time in minutes - total time accumulated when a particular goal is scored)
#' \item \code{game} (game in which goals were scored)
#' \item \code{goals} (number of goals scored in regulation period)
#' }
#' @details The World Cup is played once every four years. National teams from all over the world compete. In 2002 and in 1998, thirty-six teams were invited; whereas, in 1994 and in 1990, only 24 teams participated. The data frame \code{SOCCER} contains three columns: \code{cgt}, \code{game}, and \code{goals}. All of the information contained in \code{Soccer} is indirectly available from the FIFA World Cup website, located at \url{http://fifaworldcup.yahoo.com/}.
#' @source Chu, S. 2003. \dQuote{Using Soccer Goals to Motivate the Poisson Process.} \emph{INFORMS} Transaction on Education, \bold{3}, \bold{2}: 62-68.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' xtabs(~goals, data = SOCCER)
#' @keywords datasets
"SOCCER"
###########################################################################
#' @name STATTEMPS
#' @title Student Temperatures
#' @aliases STATTEMPS
#' @docType data
#' @description In a study conducted at Appalachian State University, students used digital oral thermometers to record their temperatures each day they came to class.  A randomly selected day of student temperatures is provided in \code{STATTEMPS}.  Information is also provided with regard to subject gender and the hour of the day when the students' temperatures were measured.
#' @format A data frame with 34 observations on the following 3 variables: 
#' \itemize{
#' \item \code{temperature} (temperature in Fahrenheit)
#' \item \code{gender} (a factor with levels \code{Female} and \code{Male})
#' \item \code{class} (a factor with levels \code{8 a.m.} and \code{9 a.m.})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' p <- ggplot(data = STATTEMPS, aes(x = gender, y =temperature, fill = class))
#' p + geom_violin()
#' @keywords datasets
"STATTEMPS"
###########################################################################
#' @name STSCHOOL
#' @title School Satisfaction
#' @aliases STSCHOOL
#' @docType data
#' @description A questionnaire is randomly administered to 11 students from State School \code{x} and to 15 students from State School \code{y}. The results have been ordered and stored in the data frame \code{STSCHOOL}.
#' @format A data frame with 26 observations on the following 4 variables: 
#' \itemize{
#' \item \code{x} (satisfaction score)
#' \item \code{y} (satisfaction score)
#' \item \code{satisfaction} (combined satisfaction scores)
#' \item \code{school} (a factor with levels \code{x} and \code{y})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' with(data = STSCHOOL, t.test(x, y, var.equal=TRUE))
#' @keywords datasets
"STSCHOOL"
############################################################################  
#' @name SUNDIG
#' @title Workstation Comparison
#' @aliases SUNDIG
#' @docType data
#' @description To compare the speed differences between two different brands of workstations (Sun and Digital), the times each brand took to complete complex simulations were recorded. Five complex simulations were selected, and the five selected simulations were run on both workstations. The resulting times in minutes for the five simulations are stored in data frame \code{SUNDIG}.
#' @format A data frame with 5 observations on the following 3 variables: 
#' \itemize{
#' \item \code{sun} (time in seconds for a Sun workstation to complete a simulation)
#' \item \code{digital} (time in seconds for a Digital workstation to complete a simulation)
#' \item \code{difference} (difference between \code{sun} and \code{digital})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' with(data = SUNDIG, t.test(sun, digital, paired=TRUE)$conf)
#' @keywords datasets
"SUNDIG"
###########################################################################
#' @name SUNFLOWER
#' @title Sunflower Defoliation
#' @aliases SUNFLOWER
#' @docType data
#' @description Seventy-two field trials were conducted by applying four defoliation treatments (non-defoliated control, 33\%, 66\%, and 100\%) at different growth stages (\code{stage}) ranging from pre-flowering (1) to physiological maturity (5) in four different locations of Navarre, Spain: Carcastillo (1), Melida (2), Murillo (3), and Unciti (4). There are two response variables: \code{yield} in kg/ha of the sunflower and \code{numseed}, the number of seeds per sunflower head. Data are stored in the data frame \code{SUNFLOWER}.
#' @format A data frame with 72 observations on the following 5 variables: 
#' \itemize{
#' \item \code{location} (a factor with levels \code{A}, \code{B}, \code{C}, and \code{D} for locations Carcastillo, Melida, Murillo, and Unciti, respectively)
#' \item \code{stage} (a factor with levels \code{stage1}, \code{stage2}, \code{stage3}, \code{stage4}, and \code{stage5})
#' \item \code{defoli} (a factor with levels \code{control}, \code{treat1}, \code{treat2}, and \code{treat3})
#' \item \code{yield} (sunflower yield in kg/ha)
#' \item \code{numseed} (number of seeds per sunflower head)
#' }
#' @source Muro, J., \emph{et. al.} 2001. \dQuote{Defoliation Effects on Sunflower Yield Reduction.} \emph{Agronomy Journal}, \bold{93}: 634-637.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' summary(aov(yield ~ stage + defoli + stage:defoli, data = SUNFLOWER))
#' ggplot(data = SUNFLOWER, aes(numseed, yield, color = defoli)) + geom_point() + 
#' geom_smooth(method = "lm", se = FALSE) + facet_grid(location ~ .)
#' @keywords datasets
"SUNFLOWER"
###########################################################################
#' @name SURFACESPAIN
#' @title Surface Area for Spanish Communities
#' @aliases SURFACESPAIN
#' @docType data
#' @description Surface area (\eqn{\texttt{km}^2}) for seventeen autonomous Spanish communities.
#' @format A data frame with 17 observations on the following 2 variables: 
#' \itemize{
#' \item \code{community} (a factor with levels \code{Andalucia}, \code{Aragon},\code{Asturias}, \code{Baleares}, \code{C.Valenciana}, \code{Canarias}, \code{Cantabria}, \code{Castilla-La Mancha}, \code{Castilla-Leon}, \code{Cataluna}, \code{Extremadura}, \code{Galicia}, \code{La Rioja}, \code{Madrid}, \code{Murcia}, \code{Navarre}, and \code{P.Vasco})
#' \item \code{surface} (surface area in \eqn{\texttt{km}^2})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' # Base Graphs
#' with(data = SURFACESPAIN, barplot(surface, names.arg = community, las = 2))
#' # ggplot2
#' ggplot(data = SURFACESPAIN, aes(x = reorder(community, surface), y = surface)) + 
#' geom_bar(stat = "identity", fill = "yellow", color = "gold") + coord_flip() + 
#' labs(x = "", y = "squared kilometers")
#' # Trellis Approach
#' barchart(community ~ surface, data = SURFACESPAIN)
#' @keywords datasets
"SURFACESPAIN"
###########################################################################
#' @name SWIMTIMES
#' @title Swim Times
#' @aliases SWIMTIMES
#' @docType data
#' @description Swimmers' improvements in seconds for two diets are stored in the data frame \code{SWIMTIMES}. The values in \code{seconds} represent the time improvement in seconds for swimmers.
#' @format A data frame with 28 observations on the following 2 variables: 
#' \itemize{
#' \item \code{seconds} (time improvement in seconds)
#' \item \code{diet} (a factor with levels \code{lowfat} and \code{highfat})
#' }
#' @details Times for the thirty-two swimmers for the 200 yard individual medley were taken right after the swimmers' conference meet. The swimmers were randomly assigned to follow one of the diets.  One group followed a low fat diet the entire year but lost two swimmers along the way. The other group followed a high fat diet the entire year and also lost two swimmers.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' wilcox.test(seconds ~ diet, data = SWIMTIMES)
#' ggplot(data = SWIMTIMES, aes(x = diet, y = seconds, fill = diet)) + geom_violin() + 
#' guides(fill = FALSE) + scale_fill_brewer()
#' @keywords datasets
"SWIMTIMES"
###########################################################################
#' @name TENNIS
#' @title Speed Detector
#' @aliases TENNIS
#' @docType data
#' @description The Yonalasee tennis club has two systems to measure the speed of a tennis ball. The local tennis pro suspects one system (\code{speed1}) consistently records faster speeds. To test her suspicions, she sets up both systems and records the speeds of 12 serves (three serves from each side of the court). The values are stored in the data frame \code{TENNIS} in the variables \code{speed1} and  \code{speed2}. The recorded speeds are in kilometers per hour.
#' @format A data frame with 12 observations on the following 2 variables: 
#' \itemize{
#' \item \code{speed1} (speed in kilometers per hour)
#' \item \code{speed2} (speed in kilometers per hour)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' with(data = TENNIS, boxplot(speed1, speed2))
#' @keywords datasets
"TENNIS"
###########################################################################
#' @name TESTSCORES 
#' @title Statistics Grades
#' @aliases TESTSCORES
#' @docType data
#' @description Test grades of 29 students taking a basic statistics course
#' @format A data frame with 29 observations on the following variable: 
#' \itemize{
#' \item \code{grade} (test score)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = TESTSCORES, aes(x = grade)) + geom_histogram(binwidth = 5, 
#' fill = "cornsilk", color = "gray60", alpha = 0.7)
#' @keywords datasets
"TESTSCORES"
###########################################################################
#' @name TIRE
#' @title Stopping Distance
#' @aliases TIRE
#' @docType data
#' @description The data frame \code{TIRE} has the stopping distances measured to the nearest foot for a standard sized car to come to a complete stop from a speed of sixty miles per hour. There are six measurements of the stopping distance for four different tread patterns labeled A, B, C, and D. The same driver and car were used for all twenty-four measurements.
#' @format A data frame with 24 observations on the following 3 variables: 
#' \itemize{
#' \item \code{stopdist} (stopping distance measured to the nearest foot)
#' \item \code{tire} (a factor with levels \code{A}, \code{B}, \code{C}, and \code{D})
#' \item \code{order} (order the experiment was conducted)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = TIRE, aes(x = reorder(tire, stopdist, FUN = median), y = stopdist, 
#' fill = tire)) + geom_boxplot() + guides(fill = FALSE) + 
#' labs(y = "Stopping distance in feet", x = "Tire Brand") + scale_fill_brewer()
#' summary(aov(stopdist ~ tire, data = TIRE))
#' p <- ggplot(data = TIRE, aes(x = reorder(tire, stopdist, FUN = mean), 
#' y = stopdist, fill = tire))
#' p + geom_boxplot(width = 0.6) + geom_dotplot(binaxis = "y", stackdir = "center",
#' binwidth = 2) + guides(fill = FALSE) + scale_fill_brewer() + 
#' stat_summary(fun.y = mean, geom = "point", fill = "black", shape = 23, size = 3) +
#' labs(x = "Tire Brand", y = "Stopping distance in feet")
#' @keywords datasets
"TIRE"
###########################################################################
#' @name TIREWEAR
#' @title Tire Wear
#' @aliases TIREWEAR
#' @docType data
#' @description The data frame \code{TIREWEAR} contains measurements for the amount of tread loss in thousandths of an inch after 10,000 miles of driving.
#' @format A data frame with 16 observations on the following 3 variables: 
#' \itemize{
#' \item \code{wear} (tread loss measured in thousandths of an inch)
#' \item \code{treat} (a factor with levels \code{A}, \code{B}, \code{C}, and \code{D})
#' \item \code{block} (a factor with levels \code{Car1}, \code{Car2}, \code{Car3}, and \code{Car4})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' par(mfrow=c(1, 2), cex = 0.8) 
#' with(data = TIREWEAR, 
#' interaction.plot(treat, block, wear, type = "b", legend = FALSE)) 
#' with(data = TIREWEAR, 
#' interaction.plot(block, treat, wear, type = "b", legend = FALSE)) 
#' par(mfrow=c(1, 1), cex = 1)
#' @keywords datasets
"TIREWEAR"
###########################################################################
#' @name TITANIC3
#' @title Titanic Survival Status
#' @aliases TITANIC3
#' @docType data
#' @description The \code{TITANIC3} data frame describes the survival status of individual passengers on the Titanic.  The  \code{TITANIC3} data frame does not contain information for the crew, but it does contain actual and estimated ages for almost 80\% of the passengers.
#' @format A data frame with 1309 observations on the following 14 variables: 
#' \itemize{
#' \item \code{pclass} (a factor with levels \code{1st}, \code{2nd}, and \code{3rd})
#' \item \code{survived} (Survival where 0 = No; 1 = Yes)
#' \item \code{name} (Name)
#' \item \code{sex} (a factor with levels \code{female} and \code{male})
#' \item \code{age} (age in years)
#' \item \code{sibsp} (Number of Siblings/Spouses Aboard)
#' \item \code{parch} (Number of Parents/Children Aboard)
#' \item \code{ticket} (Ticket Number)
#' \item \code{fare} (Passenger Fare)
#' \item \code{cabin} (Cabin)
#' \item \code{embarked} (a factor with levels \code{Cherbourg}, \code{Queenstown}, and \code{Southampton})
#' \item \code{boat} (Lifeboat Number)
#' \item \code{body} (Body Identification Number)
#' \item \code{home.dest} (Home/Destination)
#' }
#' @details Thomas Cason from the University of Virginia has greatly updated and improved the \code{titanic} data frame using the \emph{Encyclopedia Titanica} and created a new dataset called \code{TITANIC3}. This dataset reflects the state of data available as of August 2, 1999. Some duplicate passengers have been dropped; many errors have been corrected; many missing ages have been filled in; and new variables have been created.
#' @source \url{http://biostat.mc.vanderbilt.edu/twiki/pub/Main/DataSets/titanic.html}
#' @references \itemize{ 
#' \item Harrell, F. E. 2001.
#'  \emph{Regression Modeling Strategies with Applications to Linear Models, Logistic Regression, and Survival Analysis}.  Springer.
#' \item Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' }   
#' @examples
#' with(TITANIC3, table(pclass, sex))
#' @keywords datasets
"TITANIC3"
###########################################################################
#' @name TOE
#' @title Nuclear Energy
#' @aliases TOE
#' @docType data
#' @description Nuclear energy (in TOE, tons of oil equivalent) produced in 12 randomly selected European countries during 2003
#' @format A data frame with 12 observations on the following variable: 
#' \itemize{
#' \item \code{energy} (nuclear energy measured in tons of oil equivalent)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = TOE, aes(x = energy)) + geom_density(color = "red", alpha = 0.3, 
#' fill = "pink")
#' @keywords datasets
"TOE"
###########################################################################
#' @name TOP20
#' @title Tennis Income
#' @aliases TOP20
#' @docType data
#' @description \code{TOP20} contains data (in millions of dollars) corresponding to the earnings of 15 randomly selected tennis players whose earnings fall somewhere in positions 20 through 100 of ranked earnings.
#' @format A data frame with 15 observations on the following variable: 
#' \itemize{
#' \item \code{income} (yearly income in millions of dollars)
#' }
#' @source \url{http://www.atpworldtour.com/}
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = TOP20, aes(x = income)) + 
#' geom_histogram(binwidth = 1, fill = "lightblue", color = "blue") + 
#' labs(x = "yearly income in millions of dollars")
#' @keywords datasets
"TOP20"
###########################################################################
#' @name URLADDRESS
#' @title Megabytes Downloaded
#' @aliases URLADDRESS
#' @docType data
#' @description The manager of a URL commercial address is interested in predicting the number of megabytes downloaded, \code{megasd}, by clients according to the number minutes they are connected, \code{mconnected}. The manager randomly selects (megabyte, minute) pairs, and records the data.  The pairs (\code{megasd}, \code{mconnected}) are stored in the data frame \code{URLADDRESS}.
#' @format A data frame with 30 observations on the following 2 variables: 
#' \itemize{
#' \item \code{megasd} (megabytes downloaded)
#' \item \code{mconnected} (number of minutes connected)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = URLADDRESS, aes(x = mconnected, y = megasd)) + 
#' geom_point(color = "blue") + 
#' labs(x = "number of minutes connected", y = "megabytes downloaded")
#' @keywords datasets
"URLADDRESS"
###########################################################################
#' @name VIT2005
#' @title Apartments in Vitoria
#' @aliases VIT2005
#' @docType data
#' @description Descriptive information and the appraised total price (in Euros) for apartments in Vitoria, Spain
#' @format A data frame with 218 observations on the following 5 variables: 
#' \itemize{
#' \item \code{totalprice} (the market total price (in Euros) of the apartment including garage(s) and storage room(s))
#' \item \code{area} (the total living area of the apartment in square meters)
#' \item \code{zone} (a factor indicating the neighborhood where the apartment is located with levels \code{Z11}, \code{Z21}, \code{Z31}, \code{Z32}, \code{Z34}, \code{Z35}, \code{Z36}, \code{Z37}, \code{Z38}, \code{Z41}, \code{Z42}, \code{Z43}, \code{Z44}, \code{Z45}, \code{Z46}, \code{Z47}, \code{Z48}, \code{Z49}, \code{Z52}, \code{Z53}, \code{Z56}, \code{Z61}, and \code{Z62})
#' \item \code{category} (a factor indicating the condition of the apartment with levels \code{2A}, \code{2B}, \code{3A}, \code{3B}, \code{4A}, \code{4B}, and \code{5A} ordered so that \code{2A} is the best and \code{5A} is the worst)
#' \item \code{age} (age of the apartment in years)
#' \item \code{floor} (floor on which the apartment is located)
#' \item \code{rooms} (total number of rooms including bedrooms, dining room, and kitchen)
#' \item \code{out} (a factor indicating the percent of the apartment exposed to the elements: The levels \code{E100}, \code{E75}, \code{E50}, and \code{E25}, correspond to complete exposure, 75\% exposure, 50\% exposure, and 25\% exposure, respectively.)
#' \item \code{conservation} (is an ordered factor indicating the state of conservation of the apartment.  The levels \code{1A}, \code{2A}, \code{2B}, and \code{3A} are ordered from best to worst conservation.)
#' \item \code{toilets} (the number of bathrooms)
#' \item \code{garage} (the number of garages)
#' \item \code{elevator} (indicates the absence (0) or presence (1) of elevators.)
#' \item \code{streetcategory} (an ordered factor from best to worst indicating the category of the street with levels \code{S2}, \code{S3}, \code{S4}, and \code{S5})
#' \item \code{heating} (a factor indicating the type of heating with levels \code{1A}, \code{3A}, \code{3B}, and \code{4A} which correspond to: no heating, low-standard private heating, high-standard private heating, and central heating, respectively.)
#' \item \code{storage} (the number of storage rooms outside of the apartment)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = VIT2005, aes(x = area, y = totalprice, color = factor(elevator))) + 
#' geom_point()
#' modTotal <- lm(totalprice ~ area + as.factor(elevator) + area:as.factor(elevator), 
#' data = VIT2005)
#' modSimpl <- lm(totalprice ~ area, data = VIT2005)
#' anova(modSimpl, modTotal)
#' rm(modSimpl, modTotal)
#' @keywords datasets
"VIT2005"
###########################################################################
#' @name WAIT
#' @title Waiting Time
#' @aliases WAIT
#' @docType data
#' @description A statistician records how long he must wait for his bus each morning. 
#' @format A data frame with 15 observations on the following variable: 
#' \itemize{
#' \item \code{minutes} (waiting time in minutes)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' with(data= WAIT, wilcox.test(minutes, mu = 6, alternative = "less"))
#' @keywords datasets
"WAIT"
###########################################################################
#' @name WASHER
#' @title Washer Diameter
#' @aliases WASHER
#' @docType data
#' @description Diameter of circular metal disk
#' @format A data frame with 20 observations on the following variable: 
#' \itemize{
#' \item \code{diameter} (diameter of washer in cm)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = WASHER, aes(x = diameter)) + geom_density(fill = "blue", alpha = 0.2)
#' @keywords datasets
"WASHER"
###########################################################################
#' @name WATER
#' @title Sodium Content of Water
#' @aliases WATER
#' @docType data
#' @description An independent agency measures the sodium content in 20 samples from source \code{x} and in 10 samples from source \code{y} and stores them in the data frame \code{WATER}.
#' @format A data frame with 30 observations on the following 4 variables: 
#' \itemize{
#' \item \code{x} (sodium content measured in mg/L)
#' \item \code{y} (sodium content measured in mg/L)
#' \item \code{sodium} (combined sodium content measured in mg/L)
#' \item \code{source} (a factor with levels \code{x} and \code{y})
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = WATER, aes(x = sodium, y = ..density.., fill = source)) + 
#' geom_density(alpha = 0.2)
#' t.test(sodium ~ source, data = WATER, alternative = "less")
#' @keywords datasets
"WATER"
###########################################################################
#' @name WCST
#' @title Wisconsin Card Sorting Test
#' @aliases WCST
#' @docType data
#' @description The following data are the test scores from a group of 50 patients from the \emph{Virgen del Camino} Hospital (Pamplona, Spain) on the Wisconsin Card Sorting Test.
#' @format A data frame with 50 observations on the following variable: 
#' \itemize{
#' \item \code{score} (score on the Wisconsin Card Sorting Test)
#' }
#' @details The \dQuote{Wisconsin Card Sorting Test} is widely used by psychiatrists, neurologists, and neuropsychologists with patients who have a brain injury, neurodegenerative disease, or a mental illness such as schizophrenia. Patients with any sort of frontal lobe lesion generally do poorly on the test.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = WCST, aes(x = score)) + geom_density(fill = "lightblue", alpha = 0.8, 
#' color ="blue")
#' @keywords datasets
"WCST"
##########################################################################
#' @name WEIGHTGAIN
#' @title Weight Gain in Rats
#' @aliases WEIGHTGAIN
#' @docType data
#' @description The data come from an experiment to study the gain in weight of rats fed on four different diets, distinguished by amount of protein (low and high) and by source of protein (beef and cereal).
#' @format A data frame with 40 observations on the following 3 variables: 
#' \itemize{
#' \item \code{proteinsource} (a factor with levels \code{Beef} and \code{Cereal})
#' \item \code{proteinamount} (a factor with levels \code{High} and \code{Low})
#' \item \code{weightgain} (weight gained in grams)
#' }
#' @details The design of the experiment is a completely randomized design with ten rats in each of the four treatments.
#' @source Hand, D. J., F. Daly, A. D. Lunn, K. J. McConway, and E. Ostrowski. 1994. \emph{A Handbook of Small Datasets}. Chapman and Hall/CRC, London.
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. London: Chapman & Hall.
#' @examples
#' ggplot(data = WEIGHTGAIN, aes(x = proteinamount, y = weightgain, 
#' fill = proteinsource)) + geom_boxplot()
#' aov(weightgain ~ proteinsource*proteinamount, data = WEIGHTGAIN)
#' @keywords datasets
"WEIGHTGAIN"
###########################################################################
#' @name WHEATSPAIN
#' @title Wheat Surface Area in Spain
#' @aliases WHEATSPAIN
#' @docType data
#' @description Seventeen Spanish communities and their corresponding surface area (in hecatares) dedicated to growing wheat
#' @format A data frame with 17 observations on the following 3 variables: 
#' \itemize{
#' \item \code{community} (a factor with levels \code{Andalucia}, \code{Aragon}, \code{Asturias}, \code{Baleares}, \code{C.Valenciana}, \code{Canarias}, \code{Cantabria}, \code{Castilla-La Mancha}, \code{Castilla-Leon}, \code{Cataluna}, \code{Extremadura}, \code{Galicia}, \code{La Rioja}, \code{Madrid}, \code{Murcia}, \code{Navarre}, and \code{P.Vasco})
#' \item \code{hectares} (surface area measured in hectares)
#' \item \code{acres} (surface area measured in acres)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = WHEATSPAIN, aes(x = reorder(community, acres), y = acres)) + 
#' geom_bar(stat="identity", color = "orange", fill = "gold") + coord_flip() + 
#' labs(x = "") 
#' @keywords datasets
"WHEATSPAIN"
###########################################################################
#' @name WHEATUSA2004
#' @title USA Wheat Surface 2004
#' @aliases WHEATUSA2004
#' @docType data
#' @description USA's 2004 harvested wheat surface by state
#' @format A data frame with 30 observations on the following 2 variables: 
#' \itemize{
#' \item \code{states} (a factor with levels \code{AR}, \code{CA}, \code{CO}, \code{DE}, \code{GA}, \code{ID}, \code{IL}, \code{IN}, \code{KS}, \code{KY}, \code{MD}, \code{MI}, \code{MO}, \code{MS}, \code{MT}, \code{NC}, \code{NE}, \code{NY}, \code{OH}, \code{OK}, \code{OR}, \code{Other}, \code{PA}, \code{SC}, \code{SD}, \code{TN}, \code{TX}, \code{VA}, \code{WA}, and \code{WI})
#' \item \code{acres} (wheat surface area measured in thousands of acres)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = WHEATUSA2004, aes(x = reorder(states, acres), y = acres)) + 
#' geom_bar(stat = "identity", color = "gold", fill = "yellow") + coord_flip() + 
#' labs(x = "") 
#' @keywords datasets
"WHEATUSA2004"
###########################################################################
#' @name WOOL
#' @title Wool Production
#' @aliases WOOL
#' @docType data
#' @description Random sample of wool production in thousands of kilograms on 5 different days at two different locations 
#' @format A data frame with 30 observations on the following 2 variables: 
#' \itemize{
#' \item \code{production} (wool production in thousands of kilograms)
#' \item \code{location} (a factor with levels \code{textileA} and \code{textileB}.)
#' }
#' @references Ugarte, M. D., Militino, A. F., and Arnholt, A. T. 2015. \emph{Probability and Statistics with R}, Second Edition. Chapman & Hall / CRC.
#' @examples
#' ggplot(data = WOOL, aes(location, production, fill = location)) + geom_boxplot() + 
#' guides(fill = FALSE) + scale_fill_brewer()
#' t.test(production ~ location, data = WOOL)
#' @keywords datasets
"WOOL"
#####################################################################################
#####################################################################################
#' @import e1071 lattice grid
#' @importFrom ggplot2 ggplot 
#' @importFrom graphics abline axis box boxplot dotchart hist legend lines mtext par plot plot.design points polygon segments text title
#' @importFrom stats complete.cases dbinom density dnorm fitted fivenum interaction.plot ks.test median pnorm pt qchisq qnorm qqline qqnorm qt quantile rbinom rnorm rstandard sd setNames shapiro.test var
#' @importFrom utils combn
NULL