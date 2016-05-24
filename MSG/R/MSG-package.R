#' Modern Statistical Graphics
#'
#' Datasets and functions for the Chinese book ``Modern Statistical Graphics''.
#' @name MSG-package
#' @aliases MSG-package MSG
#' @docType package
#' @import graphics
#' @author Yihui Xie <\url{http://yihui.name}>
#' @keywords package
NULL


#' Random numbers containing a ``circle''
#'
#' The data was generated from two independent random varialbes (standard Normal
#' distribution) and further points on a circle were added to the data. The
#' order of the data was randomized.
#'
#' See the example section for the code to generate the data.
#' @format A data frame with 20000 observations on the following 2 variables.
#'   \describe{ \item{V1}{the first random variable with the x-axis coordinate
#'   of the circle} \item{V2}{the second random variable with the y-axis
#'   coordinate of the circle} }
#'
#' @source \url{http://yihui.name/en/2008/09/to-see-a-circle-in-a-pile-of-sand/}
#' @name BinormCircle
#' @docType data
#' @examples data(BinormCircle)
#'
#' ## original plot: cannot see anything
#' plot(BinormCircle)
#'
#' ## transparent colors (alpha = 0.1)
#' plot(BinormCircle, col = rgb(0, 0, 0, 0.1))
#'
#' ## set axes lmits
#' plot(BinormCircle, xlim = c(-1, 1), ylim = c(-1, 1))
#'
#' ## small symbols
#' plot(BinormCircle, pch = ".")
#'
#' ## subset
#' plot(BinormCircle[sample(nrow(BinormCircle), 1000), ])
#'
#' ## 2D density estimation
#' library(KernSmooth)
#' fit = bkde2D(as.matrix(BinormCircle), dpik(as.matrix(BinormCircle)))
#' # perspective plot by persp()
#' persp(fit$x1, fit$x2, fit$fhat)
#'
#' if (interactive() && require('rgl')) {
#' # perspective plot by OpenGL
#' rgl.surface(fit$x1, fit$x2, fit$fhat)
#' # animation
#' M = par3d("userMatrix")
#' play3d(par3dinterp(userMatrix = list(M, rotate3d(M,
#'    pi/2, 1, 0, 0), rotate3d(M, pi/2, 0, 1, 0), rotate3d(M, pi,
#'    0, 0, 1))), duration = 20)
#' }
#'
#' ## data generation
#' x1 = rnorm(10000); y1 = rnorm(10000)
#' x2 = rep(0.5 * cos(seq(0, 2 * pi, length = 500)), 20)
#' y2 = rep(0.5 * sin(seq(0, 2 * pi, length = 500)), 20)
#' x = cbind(c(x1, x2), c(y1, y2))
#' BinormCircle = as.data.frame(round(x[sample(20000), ], 3))
NULL


#' Life Expectancy and the Number of People with Higher Education in China
#' (2005)
#'
#' This data contains the life expectancy and number of people with higher
#' education in the 31 provinces and districts in China (2005).
#' @format A data frame with 31 observations on the following 2 variables.
#'   \describe{ \item{Life.Expectancy}{Life expectancy}
#'   \item{High.Edu.NO}{Number of people with higher education} }
#' @source China Statistical Yearbook 2005. National Bureau of Statistics.
#' @name ChinaLifeEdu
#' @docType data
#' @examples
#' data(ChinaLifeEdu)
#' x = ChinaLifeEdu
#' plot(x, type = "n", xlim = range(x[, 1]), ylim = range(x[, 2]))
#' u = par("usr")
#' rect(u[1], u[3], u[2], u[4], col = "antiquewhite",
#'     border = "red")
#' library(KernSmooth)
#' est = bkde2D(x, apply(x, 2, dpik))
#' contour(est$x1, est$x2, est$fhat, nlevels = 15, col = "darkgreen",
#'     add = TRUE, vfont = c("sans serif", "plain"))
NULL


#' Export of US and China from 1999 to 2004 in US dollars
#' @format A data frame with 13 observations on the following 3 variables.
#'   \describe{ \item{Export}{amount of export} \item{Year}{year from 1999 to
#'   2004} \item{Country}{country: US or China} }
#' @source \url{http://stat.wto.org}
#' @name Export.USCN
#' @docType data
#' @examples
#' data(Export.USCN)
#' par(mar = c(4, 4.5, 1, 4.5))
#' plot(1:13, Export.USCN$Export, xlab = "Year / Country",
#'     ylab = "US Dollars ($10^16)", axes = FALSE, type = "h",
#'     lwd = 10, col = c(rep(2, 6), NA, rep(4, 6)), lend = 1, panel.first = grid())
#' xlabel = paste(Export.USCN$Year, "\n", Export.USCN$Country)
#' xlabel[7] = ""
#' xlabel
#' abline(v = 7, lty = 2)
#' axis(1, at = 1:13, labels = xlabel, tick = FALSE, cex.axis = 0.75)
#' axis(2)
#' (ylabel = pretty(Export.USCN$Export * 8.27))
#' axis(4, at = ylabel/8.27, labels = ylabel)
#' mtext("Chinese RMB", side = 4, line = 2)
#' box()
NULL


#' Percentage data in Chinese government websites
#'
#' This data was collected from Google by searching for percentages in Chinese
#' goverment websites.
#'
#' We can specify the domain when searching in Google. For this data, we used
#' \samp{site:gov.cn}, e.g. to search for \samp{87.53\% site:gov.cn}.
#' @format A data frame with 10000 observations on the following 4 variables.
#'   \describe{ \item{percentage}{a numeric vector: the percentages}
#'   \item{count}{a numeric vector: the number of webpages corresponding to a
#'   certain percentage} \item{round0}{a logical vector: rounded to integers?}
#'   \item{round1}{a logical vector: rounded to the 1st decimal place?} }
#' @source Google (date: 2009/12/17)
#' @name gov.cn.pct
#' @docType data
#' @examples
#' data(gov.cn.pct)
#' pct.lowess = function(cond) {
#'     with(gov.cn.pct, {
#' plot(count ~ percentage, pch = ifelse(cond, 4, 20), col = rgb(0:1,
#'         0, 0, c(0.04, .5))[cond + 1], log = "y")
#'     lines(lowess(gov.cn.pct[cond, 1:2], f = 1/3), col = 2, lwd = 2)
#'     lines(lowess(gov.cn.pct[!cond, 1:2], f = 1/3), col = 1, lwd = 2)
#' })
#' }
#' par(mar = c(3.5, 3.5, 1, 0.2), mfrow = c(2, 2))
#' with(gov.cn.pct, {
#'     plot(percentage, count, type = "l", panel.first = grid())
#'     plot(percentage, count, type = "l", xlim = c(10, 11), panel.first = grid())
#'     pct.lowess(round0)
#'     pct.lowess(round1)
#' })
#' if(interactive()){
#' devAskNewPage(ask = TRUE)
#'
#' with(gov.cn.pct, {
#'     plot(count ~ percentage, type = "l")
#'     grid()
#'
#'     devAskNewPage(ask = FALSE)
#'
#'     for (i in 0:99) {
#'         plot(count ~ percentage, type = "l", xlim = i + c(0,
#'             1), panel.first = grid())
#'     }
#'
#'     devAskNewPage(ask = TRUE)
#'
#'     plot(count ~ percentage, pch = 20, col = rgb(0:1, 0, 0, c(0.07,
#'         1))[round0 + 1], log = "y")
#'     lines(lowess(gov.cn.pct[round0, 1:2], f = 1/3), col = "red",
#'         lwd = 2)
#'     lines(lowess(gov.cn.pct[!round0, 1:2], f = 1/3), col = "black",
#'         lwd = 2)
#'
#'     plot(count ~ percentage, pch = 20, col = rgb(0:1, 0, 0, c(0.07,
#'         1))[round1 + 1], log = "y")
#'     lines(lowess(gov.cn.pct[round1, 1:2], f = 1/3), col = "red",
#'         lwd = 2)
#'     lines(lowess(gov.cn.pct[!round1, 1:2], f = 1/3), col = "black",
#'         lwd = 2)
#' })
#' }
NULL


#' Number of plants corresponding to altitude
#'
#' For each altitude, the number of plants is recorded.
#' @format A data frame with 600 observations on the following 2 variables.
#'   \describe{ \item{altitude}{altitude of the area} \item{counts}{number of
#'   plants} }
#' @source
#' \url{http://cos.name/2008/11/lowess-to-explore-bivariate-correlation-by-yihui/}
#' @name PlantCounts
#' @docType data
#' @examples
#' ## different span for LOWESS
#' data(PlantCounts)
#' par(las = 1, mar = c(4, 4, 0.1, 0.1), mgp = c(2.2, 0.9, 0))
#' with(PlantCounts, {
#'     plot(altitude, counts, pch = 20, col = rgb(0, 0, 0, 0.5),
#'         panel.first = grid())
#'     for (i in seq(0.01, 1, length = 70)) {
#'         lines(lowess(altitude, counts, f = i), col = rgb(0, i,
#'             0), lwd = 1.5)
#'     }
#' })
NULL


#' The differences of P-values in t test assuming equal or unequal variances
#'
#' Given that the variances of two groups are unequal, we compute the difference
#' of P-values assuming equal or unequal variances respectively by simulation.
#'
#' See the Examples section for the generation of this data.
#' @source By simulation.
#' @format A data frame with 1000 rows and 99 columns.
#' @name t.diff
#' @docType data
#' @references Welch B (1947). ``The generalization of Student's problem when
#'   several different population variances are involved.'' Biometrika, 34(1/2),
#'   28--35.
#' @examples
#' data(t.diff)
#' boxplot(t.diff, axes = FALSE, xlab = expression(n[1]))
#' axis(1)
#' axis(2)
#' box()
#'
#' ## reproducing the data
#' if (interactive()) {
#' set.seed(123)
#' t.diff = NULL
#' for (n1 in 2:100) {
#'     t.diff = rbind(t.diff, replicate(1000, {
#'         x1 = rnorm(n1, mean = 0, sd = runif(1, 0.5, 1))
#'         x2 = rnorm(30, mean = 1, sd = runif(1, 2, 5))
#'         t.test(x1, x2, var.equal = TRUE)$p.value - t.test(x1,
#'             x2, var.equal = FALSE)$p.value
#'     }))
#' }
#' t.diff = as.data.frame(t(t.diff))
#' colnames(t.diff) = 2:100
#' }
NULL


#' Results of a Simulation to Tukey's Fast Test
#'
#' For the test of means of two samples, we calculated the P-values and recorded
#' the counts of Tukey's rule of thumb.
#'
#' See the reference for details.
#' @format A data frame with 10000 observations on the following 3 variables.
#'   \describe{ \item{pvalue.t}{P-values of t test} \item{pvalue.w}{P-values of
#'   Wilcoxon test} \item{count}{Tukey's counts} }
#' @source Simulation; see the Examples section below.
#' @name tukeyCount
#' @docType data
#' @references D. Daryl Basler and Robert B. Smawley. Tukey's Compact versus
#'   Classic Tests. \emph{The Journal of Experimental Education}, Vol. 36, No.
#'   3 (Spring, 1968), pp. 86-88
#' @examples
#' data(tukeyCount)
#'
#' ## does Tukey's rule of thumb agree with t test and Wilcoxon test?
#' with(tukeyCount, {
#'     ucount = unique(count)
#'     stripchart(pvalue.t ~ count, method = "jitter", jitter = 0.2,
#'         pch = 19, cex = 0.7, vertical = TRUE, at = ucount - 0.2,
#'         col = rgb(1, 0, 0, 0.2), xlim = c(min(count) - 1, max(count) +
#'             1), xaxt = "n", xlab = "Tukey Count", ylab = "P-values")
#'     stripchart(pvalue.w ~ count, method = "jitter", jitter = 0.2,
#'         pch = 21, cex = 0.7, vertical = TRUE, at = ucount + 0.2,
#'         add = TRUE, col = rgb(0, 0, 1, 0.2), xaxt = "n")
#'     axis(1, unique(count))
#'     lines(sort(ucount), tapply(pvalue.t, count, median), type = "o",
#'         pch = 19, cex = 1.3, col = "red")
#'     lines(sort(ucount), tapply(pvalue.w, count, median), type = "o",
#'         pch = 21, cex = 1.3, col = "blue", lty = 2)
#'     legend("topright", c("t test", "Wilcoxon test"), col = c("red",
#'         "blue"), pch = c(19, 21), lty = 1:2, bty = "n", cex = 0.8)
#' })
#'
#' if (interactive()) {
#'
#' ## this is how the data was generated
#' set.seed(402)
#' n = 30
#' tukeyCount = data.frame(t(replicate(10000, {
#'     x1 = rweibull(n, runif(1, 0.5, 4))
#'     x2 = rweibull(n, runif(1, 1, 5))
#'     c(t.test(x1, x2)$p.value, wilcox.test(x1, x2)$p.value, with(rle(rep(0:1,
#'         each = n)[order(c(x1, x2))]), ifelse(head(values, 1) ==
#'         tail(values, 1), 0, sum(lengths[c(1, length(lengths))]))))
#' })))
#' colnames(tukeyCount) = c("pvalue.t", "pvalue.w", "count")
#'
#' }
NULL



#' The scores of the game Canabalt from Twitter
#' @name canabalt
#' @docType data
#' @references
#' \samp{http://www.neilkodner.com/2011/02/visualizations-of-canabalt-scores-scraped-from-twitter/}
#' (the URL is not longer accessible)
#' @examples library(ggplot2)
#' data(canabalt)
#' print(qplot(device,score,data=canabalt))
#' print(qplot(reorder(death,score,median),score,data=canabalt,
#' geom='boxplot')+coord_flip())
NULL


#' Attributes of some music clips
#' @name music
#' @docType data
#' @references Cook D, Swayne DF (2007). Interactive and Dynamic Graphics for
#'   Data Analysis With R and GGobi. Springer. ISBN 978-0-387-71761-6.
#' @examples data(music)
NULL


#' Country power indicators of China vs America
#' @name cn_vs_us
#' @docType data
#' @references \url{http://www.guardian.co.uk/news/datablog/2011/jan/19/china-social-media}
#' @examples data(cn_vs_us)
NULL


#' Top TV earners
#'
#' The pay per episode for actors as well as other information.
#' @name tvearn
#' @docType data
#' @references \url{http://flowingdata.com/2011/02/15/visualize-this-tvs-top-earners/}
#' @examples data(tvearn)
#' plot(pay ~ rating, data=tvearn)
#' library(ggplot2)
#' qplot(pay,data=tvearn,geom='histogram',facets=gender~.,binwidth=20000)
#' qplot(rating,pay,data=tvearn,geom=c('jitter','smooth'),color=type)
NULL


#' Assists between players in CLE and LAL
#'
#' The players in the rows assisted the ones in the columns.
#' @name assists
#' @docType data
#' @references \url{http://www.basketballgeek.com/data/}
#' @examples data(assists)
#'
#' if (require('sna')) {
#' set.seed(2011)
#' gplot(assists,displaylabels=TRUE,label.cex = .7)
#' }
NULL


#' Earth quakes from 1973 to 2010
#'
#' The time, location and magnitude of all the earth quakes with magnitude being
#' greater than 6 since 1973.
#' @name quake6
#' @docType data
#' @references \url{http://cos.name/cn/topic/101510}
#' @examples data(quake6)
#' library(ggplot2)
#' qplot(year, month, data = quake6) + stat_sum(aes(size = ..n..)) +
#' scale_size(range = c(1, 10))
NULL


#' Composition of Soil from Murcia Province, Spain
#'
#' The proportions of sand, silt and clay in soil samples are given for 8
#' contiguous sites. The sites extended over the crest and flank of a low rise
#' in a valley underlain by marl near Albudeite in the province of Murcia,
#' Spain. The sites were small areas of ground surface of uniform shape
#' internally and delimited by relative discontinuities externally. Soil samples
#' were obtained for each site at 11 random points within a 10m by 10m area
#' centred on the mid-point of the site. All samples were taken from the same
#' depth. The data give the sand, silt and clay content of each sample,
#' expressed as a percentage of the total sand, silt and clay content.
#' @name murcia
#' @docType data
#' @references \url{http://www.statsci.org/data/general/murcia.html}
#' @examples data(murcia)
#' boxplot(sand~site,data=murcia)
NULL


#' Longitude and latitude of earthquakes in the Sichuan Province
#' @name eq2010
#' @docType data
#' @examples data(eq2010)
#' plot(lat~long,data=eq2010)
NULL
