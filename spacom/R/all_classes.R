################################################################################
################################################################################
## author Till Junge <till.junge@altermail.ch>                                ##
##                                                                            ##
## Copyright (c) UNIL (Universite de Lausanne)                                ##
## NCCR - LIVES (National Centre of Competence in Research "LIVES -           ##
## Overcoming vulnerability: life course perspectives",                       ##
## <http://www.lives-nccr.ch/>)                                               ##
##                                                                            ##
## spacom is free software: you can redistribute it and/or modify it under    ##
## the terms of the GNU General Public License as published by the Free       ##
## Software Foundation, either version 2 of the License or any later version. ##
##                                                                            ##
## spacom is distributed in the hope that it will be useful, but WITHOUT ANY  ##
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  ##
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      ##
## details, see <http://www.gnu.org/licenses/>.                               ##
################################################################################
################################################################################

## this class defines a kind of hash table, for which multiple keys may refer
## one and the same object without replication of the object
setClass("MultiKeyHash",
         representation(## named list of indices referring to the data list
                        indices="list",

                        ## data list
                        data="list")
         )

## Definition of an object containing everything needed to run a weighted multi-
## level analysis
setClass("MLSpawExactObject",
         representation(## data frame containing the individual level data.
                        ## data in frame must be numeric (logical, integer, or
                        ## plain numeric).
                        ## there may not be any missing values, NaN's, NULL's or
                        ## NA's
                        individual.level.data="data.frame",

                        ## data frame containing the precise contextual data.
                        ## May be NULL, but may not contain any missing values
                        ## like NaN's, NULL's or NA's
                        precise.data="data.frame",

                        ## Name of the dataframe column which contains
                        ## aggregation unit id's
                        context.id="character",

                        ## list of optional design weights at the individual
                        ## formula description of the model
                        formula="formula",

                        ## number of analyses to perform
                        nb.analyses="integer",

                        ## verbose flag
                        verbose="logical"))

setClass("MLSpawExactOutput",
         representation(## the mer model object of the lme4 analysis
                        lme="merMod",

                        ## data frame containing the standardized coefficients
                        beta="numeric"))
setClass("weightsObject",
         representation(## list of weights to apply to contextual variables
                        ## list must have same length as contextual.names
                        ## may contain NULL's for variables which should not be
                        ## weighted at the contextual level
                        distance.matrix="Matrix",

                        ## kernel function used to apply to the distance matrix
                        kernel="function",

                        ## whether the diagonal should be set to zero
                        moran="logical"))


setClass(Class="SpawExactObject",
         ## This object contains more than the minimum necessary to run
         ## SpawExact
         representation(##data frame containing the contextual data on individu-
                        ## al level data in frame must be numeric (logical,
                        ## integer, or plain numeric) there may not be any mis-
                        ## sing values, NaN's, NULL's or NA's
                        contextual.data="data.frame",

                        precise.data="NULL",

                        ## Name of the dataframe column which contains aggrega-
                        ## tion unit id's
                        context.id="character",

                        ## Name of the (perhaps) automatically generated
                        ## numbering of contextual units
                        numeric.context.id="character",

                        ## contexts
                        unique.context.ids="ANY",

                        ## list of contextual variable names to be weighted.
                        ## they can be from either contextual.data or
                        ## precise.data
                        contextual.names="list",

                        ## list of contextual variable names to be weighted
                        ## the part of the contextual.names list which refers
                        ## to contextual.data
                        aggregation.names="list",

                        ## list of contextual variable names to be weighted
                        ## the part of the contextual.names list which refers
                        ## to precise.data
                        precise.names="list",

                        ## number of upper level units
                        nb.area="integer",

                        ## number of weightings to perform in total
                        nb.analyses="integer", # i.e. length of contextual.names

                        ## number of weightings to perform on contextual.data
                        nb.aggregations="integer",

                        ## number of weightings to perform on precise.data
                        nb.precise.weightings="integer",

                        ## list of weights to apply to contextual variables
                        ## list must have same length as contextual.names
                        ## may contain NULL's for variables which should not be
                        ## weighted at the contextual level
                        contextual.weight.matrices="MultiKeyHash",

                        ## list of population weights to be applied
                        population.weight.names="MultiKeyHash",

                        ## verbose flag
                        verbose="logical")
         )







## Definition of an input object for desc when the contextual input is
## at the individual level
setClass(Class="SpawAggregateObject",
         representation(## number of resamples to be evaluated.
                        ## if zero, no boot strapping is performed
                        nb.resamples="integer",

                        ## list of aggregation functions. functions take either
                        ## a) 1 argument in which case the corresponding design
                        ##    weight is NULL
                        ## b) 2 arguments in which case the second argument is
                        ##    taken from the corresponding design weight
                        aggregation.functions="list",

                        ## vector of percentiles and their labels to be
                        ## evaluated
                        percentiles="numeric",

                        ## list of optional design weights at the individual
                        ## level used for aggregation. List must have same
                        ## length as contextual.names. May contain NULL's for
                        ## variables which should not be weighted at the indivi-
                        ## dual level
                        design.weight.names="MultiKeyHash",

                        ## sample seed is one of four things
                        ## a) NULL, in which case whatever the current random
                        ##    seed is is used
                        ## b) an integer, which will be used to set the random
                        ##    seed. this allows reproducible random samples
                        ## c) a saved .Random.seed this allows reproducible
                        ##    random samples as well. The reason why both b) and
                        ##    c) are present is because .Random.seed can be
                        ##    saved a posteriori
                        ## d) a list of samples. The user will most probably ne-
                        ##    ver use this option. It exists in order to avoid
                        ##    to redo sampling
                        sample.seed="ANY",

                        ## list of lists additional argument columns to be
                        ## passed to the aggregation function. e.g. group
                        ## identifiers. If one aggregation function takes
                        ## several additional arguments, pack them in lists,
                        ## e.g. list("group_id", list("group_id", "gender"))
                        ## for one aggregation function taking "group_id" as
                        ## additional argument and one taking both ("group_id"
                        ## and "gender")
                        additional.args="MultiKeyHash"),
         contains=c("SpawExactObject"))






setClass(Class="DescribeExactObjectLocal",
         representation(##data frame containing the contextual data on individu-
                        ## al level data in frame must be numeric (logical,
                        ## integer, or plain numeric) there may not be any mis-
                        ## sing values, NaN's, NULL's or NA's
                        contextual.data="data.frame",

                        ## precise.data="NULL",

                        ## Name of the dataframe column which contains aggrega-
                        ## tion unit id's
                        context.id="character",

                        ## list of contextual variable names to be weighted
                        contextual.names="list",

                        ## number of upper level units
                        nb.area="integer",

                        ## number of analyses to perform
                        nb.analyses="integer",

                        ## list of weights to apply to contextual variables
                        ## list must have same length as contextual.names
                        ## may contain NULL's for variables which should not be
                        ## weighted at the contextual level
                        weight.matrices="list",

                        aggregation.functions="list",

                        design.weight="character",

                        groups.name= "character", # ???
                        groups.number = "integer",
                        groups.id = "ANY"
                        ))



## Definition of an input object for desc when the contextual input is
## at the individual level
setClass(Class="DescribeAggregateObjectLocal",
         representation(## number of resamples to be evaluated
                        nb.resamples="integer",

                        ## list of aggregation functions. functions take either
                        ## a) 1 argument in which case the corresponding design
                        ##    weight is NULL
                        ## b) 2 arguments in which case the second argument is
                        ##    taken from the corresponding design weight
                        percentiles="numeric",
                        sample.seed="ANY"),

         contains=c("DescribeExactObjectLocal"))








## Definition of an output object for SpawAggregate
setClass(Class="SpawAggregateOutput",
         representation(## sample seed to be reused later
                        seed="ANY",

                        ## list of matrices per variable each containing the
                        ## aggregated contexts per sample for a contextual
                        ## variable
                        aggregated.samples="list",

                        ## list of output data.frames
                        frames="list"))









setClass("MoranObject",
         representation(## random effect matrix
                        ranefs="matrix",

                        ## distance matrix,
                        weights.object="weightsObject",

                        ## range of bandwidths to be evaluated
                        bandwidths="numeric",

                        ## number of moran.tests
                        nb.moron="integer",

                        ## number of samples
                        nb.resamples="integer",

                        ## vector of percentiles and their labels to be
                        ## evaluated
                        percentiles="numeric",

                        ## verbose flag
                        verbose="logical"
))

## Definition of an object containing everything needed to run a resampled
## weighted multilevel analysis
setClass("ResampleMLSpawExactObject",
         representation(## number of resamples to be evaluated
                        nb.resamples="integer",

                        ## sample seed is one of four things
                        ## a) NULL, in which case whatever the current random
                        ##    seed is is used
                        ## b) an integer, which will be used to set the random
                        ##    seed. this allows reproducible random samples
                        ## c) a saved .Random.seed this allows reproducible
                        ##    random samples as well. The reason why both b) and
                        ##    c) are present is because .Random.seed can be
                        ##    saved a posteriori
                        ## d) a list of samples. The user will most probably ne-
                        ##    ver use this option. It exists in order to avoid
                        ##    to redo sampling
                        individual.sample.seed="ANY",

                        ## vector of percentiles and their labels to be
                        ## evaluated
                        percentiles="numeric"

                        ),
         contains=c("MLSpawExactObject"))


## Definition of an output object for describeResampledBothContext
setClass(Class="ResampleMLSpawOutput",
         representation(## individual.samples to be reused later
                        individual.sample.seed="integer",

                        ## fixed effect data
                        fixed="data.frame",

                        ## random effect variance data
                        random.var="data.frame",

                        ## model fit data
                        model.fit="data.frame",

                        ## matrix of random effects by area for moron's eye
                        ranefs="matrix",

                        ## number of resamples
                        nb.resamples="integer",

                        ## data frame containing the standardised coefficients
                        ## and descriptives
                        betas="data.frame"))

## Definition of an object containing everything needed to run a resampled
## weighted multilevel analysis
setClass("ResampleMLSpawAggregateObject",
         representation(## The bootstrapped context in form of a
                        ## SpawAggregateOutput object,
                        aggregates="SpawAggregateOutput"),
         contains=c("ResampleMLSpawExactObject"))






 setClass("ResampleMLSpawAggregateresultObject",
         representation(
                        nb.resamples = "integer",
                        precise.frame = "data.frame",
                        individual.data = "data.frame",
                        formula="formula",
                        fixed.effect.formula="formula",
                        context.id = "character",
                        nb.area = "integer",
                        percentiles = "numeric",
                        liste = "list" )
          )










## ## Definition of an output object for describeResampledBothContext
## setClass(Class="ResampledMLSpawExactOutput",
##          representation(## contextual.samples to be reused later
##                         contextual.sample.seed="integer"),
##          contains=c("ResampleMLSpawExactOutput"))



## Definition of an input object for exploratory SRAWE
setClass(Class="ResampleExploreSpawML",
         representation(## data frame containing the individual level data.
                        ## data in frame must be numeric
                        ## (logical, integer, or plain numeric)
                        ## there may not be any missing values, NaN's, NULL's
                        ## or NA's
                        individual.level.data="data.frame",

                        ## list of 1 contextual variable name to be weighted
                        contextual.names="list",

                        ## data frame containing the contextual data on indivi-
                        ## dual level may be NULL, in which case the individual
                        ## data is used data in frame must be numeric (logical,
                        ## integer, or plain numeric)
                        ## there may not be any missing values, NaN's, NULL's or
                        ## NA's
                        contextual.data="ANY",

                        ## number of resamples
                        nb.resamples="integer",

                        ## Name of the dataframe column which contains aggrega-
                        ## tion unit id's
                        context.id="character",

                        ## formula description of the model
                        formula="formula",

                        ## Distance matrix from which weight matrices will be
                        ## computed
                        distance.matrix="Matrix",

                        ## bandwidths for the multilevel analysis
                        multilevel.bandwidths="numeric",

                        ## bandwidths for the moran test
                        moron.bandwidths="ANY",

                        ## optional design weight name
                        design.weight.name="ANY",

                        ## aggregation function
                        aggregation.function="ANY",

                        ## confidence intervals
                        percentiles="numeric",

                        ## number of areas
                        nb.area="integer",

                        ## individual.sample.seed
                        individual.sample.seed="integer",

                        ##dummy slots used for compatibility
                        nb.aggregations="integer",
                        nb.precise.weightings="integer",
                        nb.analyses="integer",
                        aggregation.names="list",
                        precise.names="list",

                        ## contextual sample seed
                        contextual.sample.seed="integer",

                        ## kernel function for proximity weighting
                        kernel="function"))



## Definition of an input object for exploratory SRAWE
setClass(Class="ExploreSpawML",
         representation(## data frame containing the individual level data.
                        ## data in frame must be numeric
                        ## (logical, integer, or plain numeric)
                        ## there may not be any missing values, NaN's, NULL's
                        ## or NA's
                        individual.level.data="data.frame",

                        ## list of 1 contextual variable name to be weighted
                        contextual.names="list",

                        ## data frame containing the precise contextual data.
                        ## May be NULL, but may not contain any missing values
                        ## like NaN's, NULL's or NA's
                        precise.data="ANY",

                        ## data frame containing the contextual data on indivi-
                        ## dual level may be NULL, in which case the individual
                        ## data is used data in frame must be numeric (logical,
                        ## integer, or plain numeric)
                        ## there may not be any missing values, NaN's, NULL's or
                        ## NA's
                        contextual.data="ANY",

                        ## number of resamples
                        nb.resamples="integer",

                        ## Name of the dataframe column which contains aggrega-
                        ## tion unit id's
                        context.id="character",

                        ## formula description of the model
                        formula="formula",

                        ## Distance matrix from which weight matrices will be
                        ## computed
                        distance.matrix="Matrix",

                        ## bandwidths for the multilevel analysis
                        multilevel.bandwidths="numeric",

                        ## optional design weight name
                        design.weight.name="ANY",

                        ## aggregation function
                        aggregation.function="ANY",

                        ## number of areas
                        nb.area="integer",

                        ##dummy slots used for compatibility
                        nb.aggregations="integer",
                        nb.precise.weightings="integer",
                        nb.analyses="integer",
                        aggregation.names="list",
                        precise.names="list",


                        ## kernel function for proximity weighting
                        kernel="function",

                        ## list of lists additional argument columns to be
                        ## passed to the aggregation function. e.g. group
                        ## identifiers. If one aggregation function takes
                        ## several additional arguments, pack them in lists,
                        ## e.g. list("group_id", list("group_id", "gender"))
                        ## for one aggregation function taking "group_id" as
                        ## additional argument and one taking both ("group_id"
                        ## and "gender")
                        additional.args="ANY",


                        ## verbose flag
                        verbose="logical"))



