#  ypd.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#  
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


################################################
# Function: ypd
#
# There are currently two different versions of an
# in sillico YPD medium:
#
# 1) Harrison, R.; Papp, B.; Pal, C.; Oliver, S. G. and Delnert, D.:
#    "Plasticity of genetic interactions in metabolic networks of yeast."
#    2007, PNAS, 104(7):2307--2312.
#    ver = harrison2007
#
# 2) Bilu, Y.; Shlomi, T.; Barkai, N. and Ruppin, E.:
#    "Conservation of expression and sequence of metabolic genes is
#     reflected by activity across metabolic states."
#    2006, PLoS Comput Biol, 2(8):932--938.
#    ver = bilu2006
#
# One can choose between the two by changing the "ver" option.


ypd <- function(model,
                def_bnd = SYBIL_SETTINGS("MAXIMUM"),
                ver = "harrison2007") {

  if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
  }

  if (mod_name(model) != "Sc_iND750") {
      warning("use Saccharomyces cerevisiae mod_name(model) == Sc_iND750")
  }

  medium <- matrix(c(
                "EX_nh4(e)",   def_bnd,
                "EX_pi(e)",    def_bnd,
                "EX_so4(e)",   def_bnd,
                "EX_glc(e)",   20,
                "EX_o2(e)",    2,
                "EX_ala_L(e)", 0.5,
                "EX_arg_L(e)", 0.5,
                "EX_asn_L(e)", 0.5,
                "EX_asp_L(e)", 0.5,
                "EX_cys_L(e)", 0.5,
                "EX_his_L(e)", 0.5,
                "EX_leu_L(e)", 0.5,
                "EX_lys_L(e)", 0.5,
                "EX_met_L(e)", 0.5,
                "EX_pro_L(e)", 0.5,
                "EX_ser_L(e)", 0.5,
                "EX_thr_L(e)", 0.5,
                "EX_trp_L(e)", 0.5,
                "EX_tyr_L(e)", 0.5,
                "EX_dcyt(e)",  0.5,
                "EX_gly(e)",   0.5,
                "EX_gua(e)",   0.5,
                "EX_thymd(e)", 0.5
                     ), ncol = 2, byrow = TRUE)

  colnames(medium) <- c("react_id", "lowbnd")

  medium_unique = NULL
  
  if (ver == "bilu2006") {
      medium_unique <- matrix(c(
                    "EX_h2o(e)",   def_bnd,
                    "EX_na1(e)",   def_bnd,
                    "EX_k(e)",     def_bnd,
                    "EX_co2(e)",   def_bnd,
                    "EX_ade(e)",   0.5,
                    "EX_gln_L(e)", 0.5,
                    "EX_ile_L(e)", 0.5,
                    "EX_phe_L(e)", 0.5,
                    "EX_val_L(e)", 0.5
                         ), ncol = 2, byrow = TRUE)
  }

  if (ver == "harrison2007") {
      medium_unique <- matrix(c(
                    "EX_h(e)",     def_bnd,
                    "EX_ura(e)",   0.5,
                    "EX_ttdca(e)", 0.5,
                    "EX_hdca(e)",  0.5,
                    "EX_ocdca(e)", 0.5,
                    "EX_glu_L(e)", 0.5,
                    "EX_ergst(e)", 0.5,
                    "EX_chol(e)",  0.5
                         ), ncol = 2, byrow = TRUE)
  }

  medium  <- rbind(medium, medium_unique)
  react   <- checkReactId(model, medium[,"react_id"])
  if (!is(react, "reactId")) {
      stop("model seems not compatible with ypd")
  }
  ex      <- findExchReact(model)

  num_med <- length(react_pos(react))
  num_ex  <- length(react_pos(ex))

  modelNEW <- changeBounds(model, ex, rep(0, num_ex), rep(def_bnd, num_ex))

  modelNEW <- changeBounds(modelNEW,
                           react,
                           -1 * as.numeric(medium[,"lowbnd"]),
                           rep(def_bnd, num_med))

  return(modelNEW)

}

