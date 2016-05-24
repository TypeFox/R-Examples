#-----------------------------------------------------------------------
#     Copyright (C) 2012-2015  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#

#-----------------------------------------------------------------------
#' Create a vector of diagonal Gaussian mixture model names.
#'
#' In a diagonal Gaussian mixture model, we assume that the variance
#' matrices are diagonal in each cluster. Assumptions on the proportions
#' and standard deviations give rise to 8 models:
#' \enumerate{
#'  \item {The proportions can be equal or free.}
#'  \item {The standard deviations can be equal or free for all the variables.}
#'  \item {The standard deviations can be equal or free for all the clusters.}
#' }
#'
#' The model names are summarized in the following array:
#' \tabular{llll}{
#'  Model Name      \tab Proportions \tab s.d. in variables \tab s.d. in clusters \cr
#'  gaussian_p_sjk  \tab Equal       \tab Free              \tab Free  \cr
#'  gaussian_p_sj   \tab Equal       \tab Free              \tab Equal \cr
#'  gaussian_p_sk   \tab Equal       \tab Equal             \tab Free  \cr
#'  gaussian_p_s    \tab Equal       \tab Equal             \tab Equal \cr
#'  gaussian_pk_sjk \tab Free        \tab Free              \tab Free  \cr
#'  gaussian_pk_sj  \tab Free        \tab Free              \tab Equal \cr
#'  gaussian_pk_sk  \tab Free        \tab Equal             \tab Free  \cr
#'  gaussian_pk_s   \tab Free        \tab Equal             \tab Equal \cr
#' }
#'
#' @param prop A character string equal to "equal", "free" or "all". Default is "all".
#' @param sdInCluster A character string equal to "equal", "free" or "all". Default is "all".
#' @param sdBetweenCluster A character string equal to "equal", "free" or "all". Default is "all".
#'
#' @return A vector of character with the model names.
#' @examples
#' clusterDiagGaussianNames()
#' ## same as c("gaussian_p_sk", "gaussian_pk_sk")
#' clusterDiagGaussianNames(prop="all", sdInCluster="equal", sdBetweenCluster= "free")
#'
#' @rdname clusterDiagGaussianNames
#' @export
clusterDiagGaussianNames <- function(prop = "all", sdInCluster="all", sdBetweenCluster = "all")
{
  if(sum(prop %in% c("equal","free","all")) != 1)
  { stop("prop is not valid. See ?clusterDiagGaussianNames for the list of prop.")}
  if(sum(sdInCluster %in% c("equal","free","all")) != 1)
  { stop("sdInCluster is not valid. See ?clusterDiagGaussianNames for the list of sdInCluster.")}
  if(sum(sdBetweenCluster %in% c("equal","free","all")) != 1)
  { stop("sdBetweenCluster is not valid. See ?clusterDiagGaussianNames for the list of sdBetweenCluster.")}

  all = c( "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_pk_sk", "gaussian_pk_s"
         , "gaussian_p_sjk", "gaussian_p_sj", "gaussian_p_sk", "gaussian_p_s")
  propFree  = c( "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_pk_sk", "gaussian_pk_s")
  propEqual = c( "gaussian_p_sjk", "gaussian_p_sj", "gaussian_p_sk", "gaussian_p_s")
  sdJFree  = c( "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_p_sjk", "gaussian_p_sj")
  sdJEqual = c( "gaussian_pk_sk", "gaussian_pk_s", "gaussian_p_sk", "gaussian_p_s")
  sdKFree  = c( "gaussian_pk_sjk", "gaussian_pk_sk", "gaussian_p_sjk", "gaussian_p_sk")
  sdKEqual = c( "gaussian_pk_sj", "gaussian_pk_s", "gaussian_p_sj", "gaussian_p_s")

  res = all;
  if (prop == "free")  { res = intersect(res, propFree);}
  if (prop == "equal") { res = intersect(res, propEqual);}
  if (sdInCluster =="free")  { res = intersect(res, sdJFree);}
  if (sdInCluster == "equal") { res = intersect(res, sdJEqual);}
  if (sdBetweenCluster =="free")  { res = intersect(res, sdKFree);}
  if (sdBetweenCluster =="equal") { res = intersect(res, sdKEqual);}

  res
}

#' check if a vector of diagonal Gaussian mixture model name is correct.
#' @param names a vector of character
#' @rdname clusterDiagGaussianNames
clusterValidDiagGaussianNames <- function(names)
{
  nb = length(names)
  if ( nb == 0 ) { return(FALSE);}

  all = clusterDiagGaussianNames();
  for (i in 1:nb)
  {  if ( sum(names[i] %in% all) != 1 ) { return(FALSE);}}
  return(TRUE)
}

#-----------------------------------------------------------------------
#' Create a vector of gamma mixture model names.
#'
#' In a gamma mixture model, we can assume that the shapes are equal in each/all
#' cluster(s) or not. We can also assume that the scales are equal in each/all
#' cluster(s) or not.
#'
#' Some configuration are impossibles. If the shapes are equal between all the
#' clusters, then the scales cannot be equal between all the clusters. Conversely
#' if the scales are equal between all the cluster, then the shapes cannot be equal
#' between all the clusters.
#'
#' This gives rise to 24 models:
#' \enumerate{
#'  \item {The proportions can be equal or free.}
#'  \item {The shapes can be equal or free in each clusters.}
#'  \item {The shapes can be equal or free between all clusters.}
#'  \item {The scales can be equal or free for each clusters.}
#'  \item {The scales can be equal or free between all clusters.}
#' }
#'
#' The model names are summarized in the following array:
#' \tabular{lllll}{
#'          &  ajk                &  ak            &  aj            &  a \cr
#'bjk       & gamma_*_ajk_bjk     & gamma_*_ak_bjk & gamma_*_aj_bjk & gamma_*_a_bjk \cr
#'bk        & gamma_*_ajk_bk      & gamma_*_ak_bk  & gamma_*_aj_bk  & gamma_*_a_bk \cr
#'bj        & gamma_*_ajk_bj      & gamma_*_ak_bj  & NA             & NA  \cr
#'b         & gamma_*_ajk_b       & gamma_*_ak_b   & NA             & NA \cr
#'}
#'
#' @param prop A character string equal to "equal", "free" or "all". Default is "all".
#' @param shapeInCluster A character string equal to "equal", "free" or "all". Default is "all".
#' @param shapeBetweenCluster A character string equal to "equal", "free" or "all". Default is "all".
#' @param scaleInCluster A character string equal to "equal", "free" or "all". Default is "all".
#' @param scaleBetweenCluster A character string equal to "equal", "free" or "all". Default is "all".
#'
#' @return A vector of character with the model names.
#' @examples
#' clusterGammaNames()
#' ## same as c("gamma_p_ak_bj", "gamma_pk_ak_bj")
#' clusterGammaNames("all", "equal", "free", "free", "equal")
#'
#' @rdname clusterGammaNames
#' @export
clusterGammaNames <- function( prop = "all", shapeInCluster="all", shapeBetweenCluster="all"
                             , scaleInCluster= "all", scaleBetweenCluster="all")
{
  if(sum(prop %in% c("equal","free","all")) != 1)
  { stop("prop is not valid. See ?clusterGammaNames for the list of prop.")}
  if(sum(shapeInCluster %in% c("equal","free","all")) != 1)
  { stop("shapeInCluster is not valid. See ?clusterGammaNames for the list of shapeInCluster.")}
  if(sum(shapeBetweenCluster %in% c("equal","free","all")) != 1)
  { stop("shapeBetweenCluster is not valid. See ?clusterGammaNames for the list of shapeBetweenCluster.")}
  if(sum(scaleInCluster %in% c("equal","free","all")) != 1)
  { stop("scaleInCluster is not valid. See ?clusterGammaNames for the list of scaleInCluster.")}
  if(sum(scaleBetweenCluster %in% c("equal","free","all")) != 1)
  { stop("scaleBetweenCluster is not valid. See ?clusterGammaNames for the list of scaleBetweenCluster.")}

  all = c( "gamma_p_ajk_bjk",  "gamma_p_ajk_bk",  "gamma_p_ajk_bj",  "gamma_p_ajk_b"
         , "gamma_p_ak_bjk",  "gamma_p_ak_bk",  "gamma_p_ak_bj",  "gamma_p_ak_b"
         , "gamma_p_aj_bjk", "gamma_p_aj_bk"
         , "gamma_p_a_bjk", "gamma_p_a_bk"
         , "gamma_pk_ajk_bjk", "gamma_pk_ajk_bk", "gamma_pk_ajk_bj", "gamma_pk_ajk_b"
         , "gamma_pk_ak_bjk", "gamma_pk_ak_bk", "gamma_pk_ak_bj", "gamma_pk_ak_b"
         , "gamma_pk_aj_bjk", "gamma_pk_aj_bk"
         , "gamma_pk_a_bjk", "gamma_pk_a_bk"
         )
  propFree = c( "gamma_pk_ajk_bjk", "gamma_pk_ajk_bk", "gamma_pk_ajk_bj", "gamma_pk_ajk_b"
              , "gamma_pk_ak_bjk", "gamma_pk_ak_bk", "gamma_pk_ak_bj", "gamma_pk_ak_b"
              , "gamma_pk_aj_bjk", "gamma_pk_aj_bk"
              , "gamma_pk_a_bjk", "gamma_pk_a_bk"
              )
  propEqual =c( "gamma_p_ajk_bjk",  "gamma_p_ajk_bk",  "gamma_p_ajk_bj",  "gamma_p_ajk_b"
              , "gamma_p_ak_bjk",  "gamma_p_ak_bk",  "gamma_p_ak_bj",  "gamma_p_ak_b"
              , "gamma_p_aj_bjk", "gamma_p_aj_bk"
              , "gamma_p_a_bjk", "gamma_p_a_bk"
              )
  shapeJFree = c( "gamma_p_ajk_bjk",  "gamma_p_ajk_bk",  "gamma_p_ajk_bj",  "gamma_p_ajk_b"
                , "gamma_p_aj_bjk",   "gamma_p_aj_bk"
                , "gamma_pk_ajk_bjk", "gamma_pk_ajk_bk", "gamma_pk_ajk_bj", "gamma_pk_ajk_b"
                , "gamma_pk_aj_bjk",  "gamma_pk_aj_bk"
                )
  shapeJEqual = c( "gamma_p_ak_bjk", "gamma_p_ak_bk",  "gamma_p_ak_bj",  "gamma_p_ak_b"
                 , "gamma_p_a_bjk",  "gamma_p_a_bk"
                 , "gamma_pk_ak_bjk","gamma_pk_ak_bk",  "gamma_pk_ak_bj",  "gamma_pk_ak_b"
                 , "gamma_pk_a_bjk", "gamma_pk_a_bk"
                 )
  shapeKFree = c( "gamma_p_ajk_bjk", "gamma_p_ajk_bk",  "gamma_p_ajk_bj",  "gamma_p_ajk_b"
                , "gamma_p_ak_bjk",  "gamma_p_ak_bk",   "gamma_p_ak_bj",   "gamma_p_ak_b"
                , "gamma_pk_ajk_bjk","gamma_pk_ajk_bk", "gamma_pk_ajk_bj", "gamma_pk_ajk_b"
                , "gamma_pk_ak_bjk", "gamma_pk_ak_bk",  "gamma_pk_ak_bj",  "gamma_pk_ak_b"
                )
  shapeKEqual = c( "gamma_p_aj_bjk",  "gamma_p_aj_bk"
                 , "gamma_p_a_bjk",   "gamma_p_a_bk"
                 , "gamma_pk_aj_bjk", "gamma_pk_aj_bk"
                 , "gamma_pk_a_bjk",  "gamma_pk_a_bk"
                 )
   scaleJFree = c( "gamma_p_ajk_bjk",  "gamma_p_ajk_bj"
                 , "gamma_p_ak_bjk",  "gamma_p_ak_bj"
                 , "gamma_p_aj_bjk"
                 , "gamma_p_a_bjk"
                 , "gamma_pk_ajk_bjk", "gamma_pk_ajk_bj"
                 , "gamma_pk_ak_bjk", "gamma_pk_ak_bj"
                 , "gamma_pk_aj_bjk"
                 , "gamma_pk_a_bjk"
                 )
   scaleJEqual = c( "gamma_p_ajk_bk", "gamma_p_ajk_b"
                  , "gamma_p_ak_bk", "gamma_p_ak_b"
                  , "gamma_p_aj_bk"
                  , "gamma_p_a_bk"
                  , "gamma_pk_ajk_bk", "gamma_pk_ajk_b"
                  , "gamma_pk_ak_bk", "gamma_pk_ak_b"
                  , "gamma_pk_aj_bk"
                  , "gamma_pk_a_bk"
           )
   scaleKFree = c( "gamma_p_ajk_bjk",  "gamma_p_ajk_bk"
                 , "gamma_p_ak_bjk",  "gamma_p_ak_bk"
                 , "gamma_p_aj_bjk", "gamma_p_aj_bk"
                 , "gamma_p_a_bjk", "gamma_p_a_bk"
                 , "gamma_pk_ajk_bjk", "gamma_pk_ajk_bk"
                 , "gamma_pk_ak_bjk", "gamma_pk_ak_bk"
                 , "gamma_pk_aj_bjk", "gamma_pk_aj_bk"
                 , "gamma_pk_a_bjk", "gamma_pk_a_bk"
                 )
   scaleKEqual = c( "gamma_p_ajk_bj",  "gamma_p_ajk_b"
                  , "gamma_p_ak_bj",  "gamma_p_ak_b"
                  , "gamma_pk_ajk_bj", "gamma_pk_ajk_b"
                  , "gamma_pk_ak_bj", "gamma_pk_ak_b"
                  )

  res = all;
  if (prop == "free")  { res = intersect(res, propFree);}
  if (prop == "equal") { res = intersect(res, propEqual);}
  if (shapeInCluster =="free")   { res = intersect(res, shapeJFree);}
  if (shapeInCluster == "equal") { res = intersect(res, shapeJEqual);}
  if (shapeBetweenCluster =="free")  { res = intersect(res, shapeKFree);}
  if (shapeBetweenCluster =="equal") { res = intersect(res, shapeKEqual);}
  if (scaleInCluster =="free")   { res = intersect(res, scaleJFree);}
  if (scaleInCluster == "equal") { res = intersect(res, scaleJEqual);}
  if (scaleBetweenCluster =="free")  { res = intersect(res, scaleKFree);}
  if (scaleBetweenCluster =="equal") { res = intersect(res, scaleKEqual);}

  res
}

#' check if a vector of gamma mixture model name is correct.
#' @param names a vector of character
#' @rdname clusterGammaNames
clusterValidGammaNames <- function(names)
{
  nb = length(names)
  if ( nb == 0 ) { return(FALSE);}

  all = clusterGammaNames();
  for (i in 1:nb)
  {  if ( sum(names[i] %in% all) != 1 ) { return(FALSE);}}
  return(TRUE)
}

#-----------------------------------------------------------------------
#' Create a vector of Categorical mixture model names.
#'
#' In a Categorical mixture model, we can build 4 models:
#' \enumerate{
#'  \item {The proportions can be equal or free.}
#'  \item {The probabilities can be equal or free for all the variables.}
#' }
#'
#' The model names are summarized in the following array:
#' \tabular{lll}{
#'  Model Name         \tab Proportions \tab Probabilities between variables  \cr
#'  categorical_p_pjk  \tab Equal       \tab Free                             \cr
#'  categorical_p_pk   \tab Equal       \tab Equal                            \cr
#'  categorical_pk_pjk \tab Free        \tab Free                             \cr
#'  categorical_pk_pk  \tab Free        \tab Equal                            \cr
#' }
#'
#' @param prop A character string equal to "equal", "free" or "all". Default is "all".
#' @param probabilities A character string equal to "equal", "free" or "all". Default is "all".
#'
#' @return A vector of character with the model names.
#' @examples
#' clusterCategoricalNames()
#' clusterCategoricalNames("all", "equal") # same as c( "categorical_pk_pk", "categorical_p_pk")
#'
#' @rdname clusterCategoricalNames
#' @export
clusterCategoricalNames <- function(prop = "all", probabilities="all")
{
  if(sum(prop %in% c("equal","free","all")) != 1)
  { stop("prop is not valid. See ?clusterCategoricalNames for the list of prop.")}
  if(sum(probabilities %in% c("equal","free","all")) != 1)
  { stop("probabilities is not valid. See ?clusterCategoricalNames for the list of probabilities.")}

  all = c( "categorical_pk_pjk", "categorical_pk_pk", "categorical_p_pjk", "categorical_p_pk")
  propFree  = c( "categorical_pk_pjk", "categorical_pk_pk")
  propEqual = c( "categorical_p_pjk", "categorical_p_pk")
  probFree  = c( "categorical_pk_pjk", "categorical_p_pjk")
  probEqual = c( "categorical_pk_pk", "categorical_p_pk")

  res = all;
  if (prop == "free")  { res = intersect(res, propFree);}
  if (prop == "equal") { res = intersect(res, propEqual);}
  if (probabilities =="free")  { res = intersect(res, probFree);}
  if (probabilities =="equal") { res = intersect(res, probEqual);}

  res
}


#' check if a vector of Categorical mixture model name is correct.
#' @param names a vector of character
#' @rdname clusterCategoricalNames
clusterValidCategoricalNames <- function(names)
{
  nb = length(names)
  if ( nb == 0 ) { return(FALSE);}

  all = c( "categorical_pk_pjk", "categorical_p_pjk",  "categorical_pk_pk", "categorical_p_pk")
  for (i in 1:nb)
  {  if ( sum(names[i] %in% all) != 1 ) { return(FALSE);}}
  return(TRUE)
}

#' Create a vector of Poisson mixture model names.
#'
#' In a Poisson mixture model, we can build 4 models:
#' \enumerate{
#'  \item {The proportions can be equal or free.}
#'  \item {The means can be equal, free or proportional for all the variables.}
#' }
#'
#' The model names are summarized in the following array:
#' \tabular{lll}{
#'  Model Name      \tab Proportions \tab Mean between variables \cr
#'  poisson_p_ljk   \tab Equal       \tab Free                   \cr
#'  poisson_p_lk    \tab Equal       \tab Equal                  \cr
#'  poisson_p_ljlk  \tab Equal       \tab Proportional           \cr
#'  poisson_pk_ljk  \tab Free        \tab Free                   \cr
#'  poisson_pk_lk   \tab Free        \tab Equal                  \cr
#'  poisson_pk_ljlk \tab Free        \tab Proportional
#' }
#'
#' @param prop A character string equal to "equal", "free" or "all". Default is "all".
#' @param mean A character string equal to "equal", "free", "proportional or "all". Default is "all".
#'
#' @return A vector of character with the model names.
#' @examples
#' clusterPoissonNames()
#' clusterPoissonNames("all", "proportional") # same as c( "poisson_pk_ljlk", "poisson_p_ljlk")
#'
#' @rdname clusterPoissonNames
#' @export
clusterPoissonNames <- function(prop = "all", mean="all")
{
  if(sum(prop %in% c("equal","free","all")) != 1)
  { stop("prop is not valid. See ?clusterPoissonNames for the list of prop.")}
  if(sum(mean %in% c("equal","free","proportional","all")) != 1)
  { stop("mean is not valid. See ?clusterPoissonNames for the list of mean.")}

  all = c( "poisson_pk_ljk", "poisson_pk_lk", "poisson_pk_ljlk", "poisson_p_ljk", "poisson_p_lk", "poisson_p_ljlk")
  propFree  = c( "poisson_pk_ljk", "poisson_pk_ljlk", "poisson_pk_lk")
  propEqual = c( "poisson_p_ljk", "poisson_p_ljlk", "poisson_p_lk")
  meanFree  = c( "poisson_pk_ljk", "poisson_p_ljk")
  meanEqual = c( "poisson_pk_lk", "poisson_p_lk")
  meanProp  = c( "poisson_pk_ljlk", "poisson_p_ljlk")

  res = all;
  if (prop == "free")  { res = intersect(res, propFree);}
  if (prop == "equal") { res = intersect(res, propEqual);}
  if (mean =="Free")  { res = intersect(res, meanFree);}
  if (mean =="equal") { res = intersect(res, meanEqual);}
  if (mean =="proportional") { res = intersect(res, meanProp);}

  res
}

#' check if a vector of Poisson mixture model name is correct.
#' @param names a vector of character
#' @rdname clusterPoissonNames
clusterValidPoissonNames <- function(names)
{
  nb = length(names)
  if ( nb == 0 ) { return(FALSE);}

  all = clusterPoissonNames();
  for (i in 1:nb)
  {  if ( sum(names[i] %in% all) != 1 ) { return(FALSE);}}
  return(TRUE)
}

#' Create a vector of Kernel mixture model names.
#'
#' In a diagonal Kernel mixture model, we can assume that
#' \enumerate{
#'  \item {The proportions can be equal or free.}
#'  \item {The standard deviations can be equal or free for all the clusters.}
#' }
#' This give rise to four models.
#'
#' The model names are summarized in the following array:
#' \tabular{lll}{
#'  Model Name           \tab Proportions \tab s.d. between clusters \cr
#'  kernelGaussian_p_sk  \tab Equal       \tab Free                  \cr
#'  kernelGaussian_p_s   \tab Equal       \tab Equal                 \cr
#'  kernelGaussian_pk_sk \tab Free        \tab Free                  \cr
#'  kernelGaussian_pk_s  \tab Free        \tab Equal                 \cr
#' }
#'
#' @param prop A character string equal to "equal", "free" or "all". Default is "all".
#' @param sdBetweenCluster A character string equal to "equal", "free" or "all". Default is "all".
#'
#' @return A vector of character with the model names.
#' @examples
#' clusterKernelNames()
#' ## same as c("kernelGaussian_p_sk", "kernelGaussian_pk_sk")
#' clusterKernelNames(prop="all", sdBetweenCluster= "free")
#'
#' @rdname clusterKernelNames
#' @export
clusterKernelNames <- function(prop = "all", sdBetweenCluster = "all")
{
  if(sum(prop %in% c("equal","free","all")) != 1)
  { stop("prop is not valid. See ?clusterKernelNames for the list of prop.")}
  if(sum(sdBetweenCluster %in% c("equal","free","all")) != 1)
  { stop("sdBetweenCluster is not valid. See ?clusterKernelNames for the list of sdBetweenCluster.")}

  all = c( "kernelGaussian_pk_sk", "kernelGaussian_pk_s", "kernelGaussian_p_sk", "kernelGaussian_p_s")
  propFree  = c( "kernelGaussian_pk_sk", "kernelGaussian_pk_s")
  propEqual = c( "kernelGaussian_p_sk", "kernelGaussian_p_s")
  sdKFree   = c( "kernelGaussian_pk_sk", "kernelGaussian_p_sk")
  sdKEqual  = c( "kernelGaussian_pk_s", "kernelGaussian_p_s")

  res = all;
  if (prop == "free")  { res = intersect(res, propFree);}
  if (prop == "equal") { res = intersect(res, propEqual);}
  if (sdBetweenCluster =="free")  { res = intersect(res, sdKFree);}
  if (sdBetweenCluster =="equal") { res = intersect(res, sdKEqual);}

  res
}

#' check if a vector of kernel mixture model name is correct.
#' @param names a vector of character
#' @rdname clusterKernelNames
clusterValidKernelNames <- function(names)
{
  nb = length(names)
  if ( nb == 0 ) { return(FALSE);}

  all = clusterKernelNames();
  for (i in 1:nb)
  {  if ( sum(names[i] %in% all) != 1 ) { return(FALSE);}}
  return(TRUE)
}
