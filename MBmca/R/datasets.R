#' Bimodal melting curve experiment on the surface of microbeads.
#' 
#' A melting curve experiment with six microbead populations and short
#' oligonucleotide probes (direct hybridization). Detection probes for human
#' VIM (vimentin), MLC-2v (myosin regulatory light chain 2, ventricular/cardiac
#' muscle isoform), SERCA2 (Atp2a2 - ATPase, Calcium-transporting ATPase
#' sarcoplasmic reticulum type, slow twitch skeletal muscle isoform), HRPT1
#' (hyperparathyroidism 1) and the artificial sequences Poly(dA)20 (20 nt of
#' dA) and aCS (artificial Control Sequence).
#' 
#' 
#' @name DMP
#' @docType data
#' @format A data frame with the melting curves of six different capture and
#' detection probe pairs on six microbeads populations for VIM, MLC-2v, SERCA2,
#' Poly(dA)20, aCS and HPRT1. First column contains the temperature (in degree
#' Celsius, 1 degree Celsius per step) followed by melting curves of VIM,
#' MLC-2v, SERCA2, Poly(dA)20, aCS and HPRT1 with bimodal melting pattern. The
#' dyes and quencher used were Atto 647N and BHQ2.  \describe{
#' \item{T}{a numeric vector, Temperature in degree Celsius.}
#' 
#' \item{VIM.no.Q}{a numeric vector, VIM without quencher and without
#' Poly(dT)20 region.}
#' 
#' \item{MLC2v.w.Q.w.dT20}{a numeric vector, MLC-2v with
#' quencher-labeled detection probe and fluorescent Poly(dA)20 detection
#' probe.}
#' 
#' \item{SERCA2.no.Q.w.dT20}{a numeric vector, SERCA2 without
#' quencher-labeled detection probe and Poly(dA)20 detection probe.}
#' 
#' \item{PolydA20.w.Q}{a numeric vector, Poly(dT)20 with fluorescent
#' Poly(dA)20 detection probe (quencher labeled).}
#' 
#' \item{aCS.w.Q.w.dT20}{a numeric vector, artificial Control Sequence
#' without quencher-labeled detection probe and fluorescent Poly(dA)20
#' detection probe.}
#' 
#' \item{HPRT1.no.Q.w.dT20}{a numeric vector, HPRT1 without
#' quencher-labeled detection probe and fluorescent Poly(dA)20 detection
#' probe.} }
#' @seealso \code{\link{MFIerror}}, \code{\link{mcaSmoother}},
#' \code{\link{diffQ}}, \code{\link{diffQ}}, \code{\link{MultiMelt}},
#' \code{\link{DualHyb}}
#' @source Data were measured with the VideoScan platform:
#' 
#' A Highly Versatile Microscope Imaging Technology Platform for the Multiplex
#' Real-Time Detection of Biomolecules and Autoimmune Antibodies. S. Roediger,
#' P. Schierack, A. Boehm, J. Nitschke, I. Berger, U. Froemmel, C. Schmidt, M.
#' Ruhland, I. Schimke, D. Roggenbuck, W. Lehmann and C. Schroeder.
#' \emph{Advances in Biochemical Bioengineering/Biotechnology}. 133:35--74,
#' 2013. \url{http://www.ncbi.nlm.nih.gov/pubmed/22437246}
#' 
#' Surface Melting Curve Analysis with R. S. Roediger, A. Boehm and I.
#' Schimke. \emph{The R Journal}. 5(2):37--52, 2013.
#' \url{http://journal.r-project.org/}
#' @keywords datasets
#' @examples
#' 
#' data(DMP)
#' 
NULL





#' Surface melting curve data from direct hybridization experiment of short
#' oligonucleotide probes with bimodal melting curve pattern.
#' 
#' A melting curve experiment with four microbead populations and short
#' oligonucleotide probes (direct hybridization) and longer probes
#' (dual-hybridization probes) capture probe. Detection probes for human VIM
#' (vimentin), MLC-2v (myosin regulatory light chain 2, ventricular/cardiac
#' muscle isoform) and SERCA2 (Atp2a2 - ATPase, Calcium-transporting ATPase
#' sarcoplasmic reticulum type, slow twitch skeletal muscle isoform). One
#' sequence of VIM contained a mutation at position 41.
#' 
#' The melting curve was conducted with short oligonucleotide probes (direct
#' hybridization) and longer probes (dual-hybridization probes) on the surface
#' of microbeads (sequences and materials according to Roediger et al.  (2012))
#' using the VideoScan platform by Roediger et al. (2012). The dyes and
#' quencher used were Atto 647N and BHQ2.
#' 
#' @name DualHyb
#' @docType data
#' @format A data frame with the melting curves of three different capture and
#' detection probe pairs for HRPT1 and MLC-2v. First column contains the
#' temperature (in degree Celsius, 0.5 degree Celsius per step) followed by
#' melting curves of HRPT1 on 12 microbead populations and melting curves of
#' MLC-2v on 12 microbead populations.  \describe{ \item{T}{a numeric
#' vector, Temperature in degree Celsius.}
#' 
#' \item{MLC2v}{a numeric vector, MLC-2v with quencher-labeled
#' detection probe}
#' 
#' \item{SERCA2}{a numeric vector, SERCA2 without quencher-labeled
#' detection probe}
#' 
#' \item{VIM.w.Mutation}{a numeric vector, mutated VIM with
#' quencher-labeled detection probe}
#' 
#' \item{VIM.wo.Mutation}{a numeric vector, native VIM with
#' quencher-labeled detection probe} }
#' @seealso \code{\link{MFIerror}}, \code{\link{mcaSmoother}},
#' \code{\link{diffQ}}, \code{\link{diffQ2}}, \code{\link{DMP}},
#' \code{\link{MultiMelt}}
#' @source Data were measured with the VideoScan platform:
#' 
#' A Highly Versatile Microscope Imaging Technology Platform for the Multiplex
#' Real-Time Detection of Biomolecules and Autoimmune Antibodies. S. Roediger,
#' P. Schierack, A. Boehm, J. Nitschke, I. Berger, U. Froemmel, C. Schmidt, M.
#' Ruhland, I. Schimke, D. Roggenbuck, W. Lehmann and C. Schroeder.
#' \emph{Advances in Biochemical Bioengineering/Biotechnology}. 133:35--74,
#' 2013. \url{http://www.ncbi.nlm.nih.gov/pubmed/22437246}
#' 
#' Surface Melting Curve Analysis with R. S. Roediger, A. Boehm and I.
#' Schimke. \emph{The R Journal}. 5(2):37--52, 2013.
#' \url{http://journal.r-project.org/}
#' 
#' Nucleic acid detection based on the use of microbeads: a review. S.
#' Roediger, C. Liebsch, C. Schmidt, W. Lehmann, U. Resch-Genger, U. Schedler,
#' P. Schierack. \emph{Microchim Acta} 2014:1--18. DOI:
#' 10.1007/s00604-014-1243-4
#' @keywords datasets
#' @examples
#' 
#' data(DualHyb)
#' 
NULL





#' Surface melting curve data from direct hybridization experiment of short
#' oligonucleotides.
#' 
#' A melting curve experiment with twelve microbead populations and the short
#' oligonucleotide capture probe and detection probe for human HRPT1
#' (hyperparathyroidism 1) and MLC-2v (myosin regulatory light chain 2,
#' ventricular/cardiac muscle isoform).
#' 
#' The melting curve was conducted with short oligonucleotide probes on the
#' surface of microbeads using the VideoScan platform according to Roediger et
#' al. (2012). The dyes and quencher used were Atto 647N and BHQ2.
#' 
#' @name MultiMelt
#' @docType data
#' @format A data frame with the melting curves of two different capture and
#' detection probe pairs for HRPT1 and MLC-2v. First column contains the
#' temperature (in degree Celsius, 1 degree Celsius per step) followed by
#' melting curves of HRPT1 on twelve microbead populations and melting curves
#' of MLC-2v on twelve microbead populations.  \describe{ \item{T}{a
#' numeric vector for the temperature in degree Celsius}
#' \item{HPRT1.1}{a numeric vector, as HPRT1.1 of detection/capture
#' probe HPRT1/HPRT1-cap on microbead population 1} \item{HPRT1.2}{a
#' numeric vector, as HPRT1.2 on microbead population 2}
#' \item{HPRT1.3}{a numeric vector, as HPRT1.3 on microbead population
#' 3} \item{HPRT1.4}{a numeric vector, as HPRT1.4 on microbead
#' population 4} \item{HPRT1.5}{a numeric vector, as HPRT1.5 on
#' microbead population 5} \item{HPRT1.6}{a numeric vector, as HPRT1.6
#' on microbead population 6} \item{HPRT1.7}{a numeric vector, as
#' HPRT1.7 on microbead population 7} \item{HPRT1.8}{a numeric vector,
#' as HPRT1.8 on microbead population 8} \item{HPRT1.9}{a numeric
#' vector, as HPRT1.9 on microbead population 9} \item{HPRT1.10}{a
#' numeric vector, as HPRT1.10 on microbead population 10}
#' \item{HPRT1.11}{a numeric vector, as HPRT1.11 on microbead
#' population 11} \item{HPRT1.12}{a numeric vector, as HPRT1.12 on
#' microbead population 12} \item{MLC2v1}{a numeric vector, as MLC2v1
#' of detection/capture probe MLC-2v/MLC-2v-cap on microbead population 1}
#' \item{MLC2v2}{a numeric vector, as MLC2v2 on microbead population 2}
#' \item{MLC2v3}{a numeric vector, as MLC2v3 on microbead population 3}
#' \item{MLC2v4}{a numeric vector, as MLC2v4 on microbead population 4}
#' \item{MLC2v5}{a numeric vector, as MLC2v5 on microbead population 5}
#' \item{MLC2v6}{a numeric vector, as MLC2v6 on microbead population 6}
#' \item{MLC2v7}{a numeric vector, as MLC2v7 on microbead population 7}
#' \item{MLC2v8}{a numeric vector, as MLC2v8 on microbead population 8}
#' \item{MLC2v9}{a numeric vector, as MLC2v9 on microbead population 9}
#' \item{MLC2v10}{a numeric vector, as MLC2v10 on microbead population
#' 10} \item{MLC2v11}{a numeric vector, as MLC2v11 on microbead
#' population 11} \item{MLC2v12}{a numeric vector, as MLC2v12 on
#' microbead population 12} }
#' @seealso \code{\link{MFIerror}}, \code{\link{mcaSmoother}},
#' \code{\link{diffQ}}, \code{\link{diffQ2}}, \code{\link{DMP}},
#' \code{\link{DualHyb}}
#' @source Data were measured with the VideoScan platform:
#' 
#' A Highly Versatile Microscope Imaging Technology Platform for the Multiplex
#' Real-Time Detection of Biomolecules and Autoimmune Antibodies. S. Roediger,
#' P. Schierack, A. Boehm, J. Nitschke, I. Berger, U. Froemmel, C. Schmidt, M.
#' Ruhland, I. Schimke, D. Roggenbuck, W. Lehmann and C. Schroeder.
#' \emph{Advances in Biochemical Bioengineering/Biotechnology}. 133:35--74,
#' 2013. \url{http://www.ncbi.nlm.nih.gov/pubmed/22437246}
#' 
#' Surface Melting Curve Analysis with R. S. Roediger, A. Boehm and I.
#' Schimke. \emph{The R Journal}. 5(2):37--52, 2013, 2013.
#' \url{http://journal.r-project.org/}
#' @keywords datasets
#' @examples
#' 
#' data(MultiMelt)
#' 
NULL