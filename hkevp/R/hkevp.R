#' @title
#' A hierarchical model for spatial extreme values
#' 
#' @author
#' Quentin Sebille
#'
#'
#'
#' @description
#' This package focus on a particular hierarchical model for spatial extreme values: the HKEVP (Hierarchical Kernel Extreme Value Process) of \cite{Reich and Shaby (2012)}.
#' 
#' This model is designed for inference on point-referenced block maxima across a region of space. Both marginal parameters and spatial dependence structure are modelled in a Bayesian framework (see details). Estimation procedure of the HKEVP is performed via a Metropolis-within-Gibbs algorithm that returns posterior distributions of all dependence and marginal parameters.
#' 
#' Simulation and fitting procedure for this model are provided along with several other tools as marginal distribution extrapolation and spatial prediction, which is used for instance in \cite{Shaby and Reich (2012)}.
#' 
#' 
#' 
#'
#' @docType package
#' @name hkevp
#' 
#' 
#' @details
#' The functions included in \code{hkevp} package are listed below:
#' 
#' \enumerate{
#' \item \code{\link{hkevp.fit}}: fits the HKEVP to sptial block maxima data.
#' \item \code{\link{hkevp.rand}}: simulates data from the HKEVP.
#' \item \code{\link{hkevp.expmeasure}}: computes the exponent measure of the HKEVP.
#' \item \code{\link{mcmc.fun}}: applies a function to the main Markov chains obtained by \code{hkevp.fit}. Useful to compute the posterior means or quantiles for instance.
#' \item \code{\link{mcmc.plot}}: plots the resulting Markov chains, in order to visually assess convergence.
#' \item \code{\link{mcmc.accept}}: computes and plots the acceptance rates of the MCMC procedure for each parameter of the HKEVP.
#' \item \code{\link{posterior.expmeasure}}: computes the posterior distribution of the exponent measure from the output of \code{hkevp.fit}.
#' \item \code{\link{posterior.extrapol}}: computes the value of each GEV parameter at target positions (marginal extrapolation).
#' \item \code{\link{posterior.predict}}: predicts the value of the max-stable process at target positions (spatial prediction).
#' }
#' 
#' 
#' 
#' The HKEVP of \cite{Reich and Shaby (2012)} is a hierarchical spatial max-stable model for extreme values. Max-stable models arise as the limiting distribution of renormalized maxima of stochastic processes and they generalize the Extreme Value Theory (EVT) to the infinite-dimensional case. For more information about EVT and max-stable processes, see for instance \cite{Beirlant et al. (2004)}, \cite{de Haan and Ferreira (2006)} and \cite{Coles (2001)}. Statistical inference on spatial extremes are reviewed in \cite{Cooley et al. (2012)} and \cite{Davison et al. (2012)}.
#' 
#' Let \eqn{Y(\cdot)} be the process of block maxima recorded over a spatial region. Assume this process is max-stable, then all marginal distributions are necessarily GEV, i.e:
#' \deqn{Y(s) \sim GEV\{\mu(s),\sigma(s),\xi(s)\}}
#' with location, scale and shape parameters \eqn{\mu(s)}, \eqn{\sigma(s)} and \eqn{\xi(s)} respectively, indexed by position \eqn{s} in space. Without loss of generality, one may look at the \emph{simple max-stable process} \eqn{Z(\cdot)} with \eqn{GEV(1,1,1)} margins, where spatial dependence is contained unconditionally of the marginals.
#' 
#' The HKEVP is thus defined by assuming that there exist \eqn{\alpha\in(0,1]} and a set of knots \eqn{\{v_1,...,v_L\}} with associated kernel functions \eqn{\{\omega_1(\cdot),...,\omega_L(\cdot)\}} such that:
#' \deqn{Z(s) = U(s)\theta(s)}
#' where \eqn{U(s)} is a spatially-independent process with \eqn{GEV(1,\alpha,\alpha)} margins and
#' \deqn{\theta(s) = \left[ \sum_{\ell=1}^L A_{\ell} \omega_{\ell}(s)^{1/\alpha}\right]^{\alpha}}
#' is the \emph{residual dependence process}, defined with a random variable \eqn{A_{\ell} \sim PS(\alpha)}, the positive stable distribution with characteristic exponent \eqn{\alpha}. Note that the kernel functions must satisfy the condition \eqn{\sum_{\ell=1}^L \omega_{\ell}(s) = 1}, for all \eqn{s}.
#' 
#' Under those assumptions, \cite{Reich and Shaby (2012)} showed that this model lead to an explicit formula for the distribution of a joint vector of maxima, which is not the case for max-stable processes since no general parametric representation exist (cf. \cite{de Haan and Ferreira (2006)}). The joint distribution of the vector \eqn{\{Z(s_1), \ldots, Z(s_n)\}} under the HKEVP assumptions can indeed be written:
#' \deqn{P\{ Z_1(s_1)<z_1, \ldots, Z_n(s_n)<z_n \} = \sum_{\ell=1}^L \left[ \sum_{i=1}^n \left(\frac{\omega_\ell(s_i)}{z_i}\right)^{1/\alpha}\right]^{\alpha} ~.}
#' 
#' The HKEVP can be seen as an approximation of the \cite{Smith (1990)} model, but with an additional dependence parameter \eqn{\alpha} which controls the strength of spatial dependence. Low value for \eqn{\alpha} means strong dependence structure while \eqn{\alpha=1} means independence.
#' 
#' In the article of \cite{Reich and Shaby (2012)}, the kernels are chosen to be standardized Gaussian kernels, i.e:
#' \deqn{\omega_\ell(s) = \frac{K(s_\ell|v_\ell,\tau)}{\sum_{j=1}^L K(s_j|v_j,\tau)}}
#' where \eqn{K(\cdot|v,\tau)} is the Gaussian kernel centered at \eqn{v_\ell} and \eqn{\tau} is a bandwidth parameter.
#' 
#' Conditionally on the marginal and the dependence parameters, we obtain independent responses which allow the use of statistical inference. Since the model is structured via multiple layers, Bayesian inference is preferred. Details of the MCMC algorithm can be found in \cite{Reich and Shaby (2012)}.
#' 
#' Marginal parameters \eqn{\mu(s)}, \eqn{\gamma(s)=\log \sigma(s)} and \eqn{\xi(s)} are modelled through latent Gaussian processes.
#' 
#' 
#' 
#'
#' @note
#' I would like to thank Brian Reich and Benjamin Shaby for their help regarding the implementation of the inference procedure of the HKEVP, and Mathieu Ribatet for introducing me to the Bayesian world. I also acknowledge Electricite de France (EDF) for the financial support and the helpful discussions around the HKEVP with Anne Dutfoy, Marie Gallois, Thi Thu Huong Hoang and Sylvie Parey. Finally, I thank my two PhD supervisors Anne-Laure Fougeres and Cecile Mercadier for their reviews of this package.
#' 
#' 
#' 
#' 
#' 
#' @references 
#' Beirlant, J., Goegebeur, Y., Segers, J. J. J., & Teugels, J. (2004). Statistics of Extremes: Theory and Applications. <DOI:10.1002/0470012382>
#' 
#' Coles, S. (2001). An Introduction to Statistical Modeling of Extreme Values. Springer Science & Business Media. <10.1007/978-1-4471-3675-0>
#' 
#' Cooley, D., Cisewski, J., Erhardt, R. J., Jeon, S., Mannshardt, E., Omolo, B. O., & Sun, Y. (2012). A survey of spatial extremes: Measuring spatial dependence and modeling spatial effects. Revstat, 10(1), 135-165.
#' 
#' Davison, A. C., Padoan, S. A., & Ribatet, M. (2012). Statistical modeling of spatial extremes. Statistical Science, 27(2), 161-186. <DOI:10.1214/11-STS376>
#' 
#' de Haan, L., & Ferreira, A. (2006). Extreme Value Theory: An Introduction. Springer Science & Business Media. <DOI:10.1007/0-387-34471-3>
#' 
#' Reich, B. J., & Shaby, B. A. (2012). A hierarchical max-stable spatial model for extreme precipitation. The annals of applied statistics, 6(4), 1430. <DOI:10.1214/12-AOAS591>
#' 
#' Shaby, B. A., & Reich, B. J. (2012). Bayesian spatial extreme value analysis to assess the changing risk of concurrent high temperatures across large portions of European cropland. Environmetrics, 23(8), 638-648. <DOI:10.1002/env.2178>
#' 
#' Smith, R. L. (1990). Max-stable processes and spatial extremes. Unpublished manuscript.
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
NULL