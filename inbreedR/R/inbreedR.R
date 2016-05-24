#' inbreedR: Workflows for analysing variance in inbreeding and HFCs based on SNP or microsatellite markers.
#'
#' @description 
#' 
#' \code{inbreedR} contains the following functions:
#'
#' \link{g2_microsats}
#' \link{g2_snps}
#' \link{convert_raw}
#' \link{check_data}
#' \link{r2_hf}
#' \link{r2_Wf}
#' \link{HHC}
#' \link{sMLH}
#' \link{MLH}
#' \link{simulate_g2}
#' \link{simulate_r2_hf}
#' \link{plot.inbreed}
#' \link{print.inbreed}
#' 
#'
#' @details 
#' 
#' A correlation between heterozygosity (h) and fitness (W) requires a simultaneous effect of 
#' inbreeding level (f) on both of them.
#' A heterozygosity-fitness correlation (HFC) thus is the product of two correlations, 
#' which can be summarized in the following equation:
#' 
#' \deqn{r(W, h) = r(W, f)r(h, f)}
#' 
#' Estimating these parameters and their sensitivity towards the number and type of genetic markers used is the central framework
#' of the inbreedR package. At the heart of measuring inbreeding based on genetic markers is the g2 statistic, which estimates 
#' the correlation of heterozygosity across markers, called identity disequilibrium (ID). ID is a proxy for inbreeding.
#' 
#' The package has three main goals:
#' 
#' \itemize{
#' \item Assessing identity disequilibria and the potential to detect heterozygosity-fitness correlations 
#' \item Providing insights on the sensitivity of these measures based on the number/type of molecular markers used
#' \item Implementing computationally efficient functions in a flexible environment for analysing 
#' inbreeding and HFC`s with both small and large datasets.
#' }
#' 
#' For a short introduction to inbreedR start with the vignette:
#' \code{browseVignettes(package = "inbreedR")}
#'
#' @author  Martin Stoffel (martin.adam.stoffel@@gmail.com), Mareike Esser (messer@@uni-bielefeld.de)
#'
#' @references
#' Slate, J., David, P., Dodds, K. G., Veenvliet, B. A., Glass, B. C., Broad, T. E., & McEwan, J. C. (2004). 
#' Understanding the relationship between the inbreeding coefficient 
#' and multilocus heterozygosity: theoretical expectations and empirical data. Heredity, 93(3), 255-265.
#' 
#' Szulkin, M., Bierne, N., & David, P. (2010). HETEROZYGOSITY-FITNESS CORRELATIONS: A TIME FOR REAPPRAISAL. 
#' Evolution, 64(5), 1202-1217.
#' 
#' David, P., Pujol, B., Viard, F., Castella, V. and Goudet, J. (2007),
#' Reliable selfing rate estimates from imperfect population genetic data. Molecular Ecology,
#' 16: 2474
#'
#' Hoffman, J.I., Simpson, F., David, P., Rijks, J.M., Kuiken, T., Thorne, M.A.S., Lacey, R.C. & Dasmahapatra, K.K. (2014) High-throughput sequencing reveals inbreeding depression in a natural population.
#' Proceedings of the National Academy of Sciences of the United States of America, 111: 3775-3780. 
#'
#' @docType package
#' @name inbreedR
NULL
