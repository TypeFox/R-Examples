#' An example dataset
#'
#' A large dataset of provenance data from Namibia comprised of 14
#' sand samples from the Namib Sand Sea and 2 samples from the Orange
#' River.
#'
#' \code{Namib} is a list containing the following 6 items:
#'
#' \code{DZ}: a \code{distributional} dataset containing the zircon
#' U-Pb ages for ca. 100 grains from each sample, as well as their
#' (1-sigma) analytical uncertainties.
#'
#' \code{PT}: a \code{compositional} dataset with the bulk petrography
#' of the samples, i.e. the quartz ('Q'), K-feldspar ('KF'),
#' plagioclase ('P'), and lithic fragments of metamorphic ('Lm'),
#' volcanic ('Lv') and sedimentary ('Ls') origin.
#'
#' \code{HM}: a \code{compositional} dataset containing the heavy
#' mineral composition of the samples, comprised of zircon ('zr'),
#' tourmaline ('tm'), rutile ('rt'), Ti-oxides ('TiOx'), titanite
#' ('sph'), apatite ('ap'), epidote ('ep'), garnet ('gt'), staurolite
#' ('st'), andalusite ('and'), kyanite ('ky'), sillimanite ('sil'),
#' amphibole ('amp'), clinopyroxene ('cpx') and orthopyroxene ('opx').
#'
#' \code{PTHM}: a \code{compositional} dataset combining the variables
#' contained in \code{PT} and \code{HM} plus 'mica', 'opaques',
#' 'turbids' and 'other' transparent heavy minerals ('LgM'),
#' normalised to 100.
#' 
#' \code{Major}: a \code{compositional} dataset listing the
#' concentrations (in wt%) of SiO2, Al2O3, Fe2O3, MgO, CaO, Na2O, K2O,
#' TiO2, P2O5 and MnO.
#'
#' \code{Trace}: a \code{compositional} data listing the concentrations
#' (in ppm) of Rb, Sr, Ba, Sc, Y, La, Ce, Pr, Nd, Sm, Gd, Dy, Er, Yb, Th,
#' U, Zr, Hf, V, Nb, Cr, Co, Ni, Cu, Zn, Ga and Pb.
#' 
#' @name Namib
#' @docType data
#' @examples
#' data(Namib)
#' samp <- Namib$DZ$x[['N1']]
#' dens <- KDE(samp,0,3000)
#' plot(dens)
#' @author Pieter Vermeesch and Eduardo Garzanti
#' @references Vermeesch, P. and Garzanti, E., Making geological sense
#' of 'Big Data' in sedimentary provenance analysis, Chemical Geology
#' 409 (2015) 20-27
NULL

#' A list of rock and mineral densities
#'
#' List of rock and mineral densities using the following
#' abbreviations: Q (quartz), KF (K-feldspar), P (plagioclase), F
#' (feldspar), Lvf (felsic/porfiritic volcanic rock fragments), Lvm
#' (microlithic / porfiritic / trachitic volcanic rock fragments), Lcc
#' (calcite), Lcd (dolomite), Lp (marl), Lch (chert), Lms
#' (argillaceous / micaceous rock fragments), Lmv (metavolcanics), Lmf
#' (metasediments), Lmb (metabasites), Lv (volcanic rock fragments),
#' Lc (carbonates), Ls (sedimentary rock fragments), Lm (metamorphic
#' rock fragments), Lu (serpentinite), mica, opaques, FeOx
#' (Fe-oxides), turbids, zr (zircon), tm (tourmaline), rt (rutile),
#' TiOx (Ti-oxides), sph (titanite), ap (apatite), mon (monazite), oth
#' (other minerals), ep (epidote), othLgM (prehnite + pumpellyite +
#' lawsonite + carpholite), gt (garnet), ctd (chloritoid), st
#' (staurolite), and (andalusite), ky (kyanite), sil (sillimanite),
#' amp (amphibole), px (pyroxene), cpx (clinopyroxene), opx
#' (orthopyroxene), ol (olivine), spinel and othHM (other heavy
#' minerals).
#' 
#' @name densities
#' @seealso restore, minsorting
#' @docType data
#' @examples
#' data(Namib,densities)
#' N8 <- subset(Namib$HM,select="N8")
#' distribution <- minsorting(N8,densities,phi=2,sigmaphi=1,medium="air",by=0.05)
#' plot(distribution)
#' @author Alberto Resentini and Pieter Vermeesch
#' @references Resentini, A, Malusa M G and Garzanti, E. "MinSORTING:
#' An Excel worksheet for modelling mineral grain-size distribution in
#' sediments, with application to detrital geochronology and
#' provenance studies." Computers & Geosciences 59 (2013): 90-97.
#'
#' Garzanti, E, Ando, S and Vezzoli, G. "Settling
#' equivalence of detrital minerals and grain-size dependence of
#' sediment composition." Earth and Planetary Science Letters 273.1
#' (2008): 138-151.
NULL

#' Petrographic end-member compositions
#'
#' A compositional dataset comprising the mineralogical compositions
#' of the following end-members: \code{undissected_magmatic_arc},
#' \code{dissected_magmatic_arc}, \code{ophiolite},
#' \code{recycled_clastic},
#' \code{undissected_continental_block},
#' \code{transitional_continental_block},
#' \code{dissected_continental_block},
#' \code{subcreted_axial_belt} and
#' \code{subducted_axial_belt}
#' @name endmembers
#' @seealso minsorting
#' @docType data
#' @examples
#' data(endmembers,densities)
#' ophiolite <- subset(endmembers,select="ophiolite")
#' plot(minsorting(ophiolite,densities,by=0.05))
#' @author Alberto Resentini and Pieter Vermeesch
#' @references Resentini, A, Malusa M G and Garzanti, E. "MinSORTING:
#' An Excel worksheet for modelling mineral grain-size distribution in
#' sediments, with application to detrital geochronology and
#' provenance studies." Computers & Geosciences 59 (2013): 90-97.
#'
#' Garzanti, E, Ando, S and Vezzoli, G. "Settling
#' equivalence of detrital minerals and grain-size dependence of
#' sediment composition." Earth and Planetary Science Letters 273.1
#' (2008): 138-151.
NULL
