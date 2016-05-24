#' Full-Scale Simulations of the Fossil Record with Birth, Death and Sampling of Morphotaxa
#'
#' A complete birth-death-sampling branching simulator that captures morphological-taxon identity
#' of lineages, as is typically discussed in models of paleontological data. This function
#' allows for the use of precise point constraints to condition simulation run acceptance and
#' can interpret complex character strings given as rate values for use in modeling
#' complex processes of diversification and sampling.

#' @details
#' \code{simFossilRecord} simulates a birth-death-sampling branching process (ala Foote, 1997,
#'  2000; Stadler, 2009) in which lineages of organisms may branch, go extinct or be sampled
#' thought a continuous time-interval, with the occurrence of these events
#' modeled as Poisson process controlled by some set of instantaneous rates.
#' This model is ultimately based on the birth-death model (Kendall, 1948; Nee, 2006),
#' which is widely implemented in many R packages. Unlike other such typical branching
#' simulators, this function enmeshes the lineage units within explicit models of how
#' lineages are morphologically differentiated (Bapst, 2013). This is key to allow comparison
#' to datasets from the fossil record, as morphotaxa are the basic
#' units of paleontological diversity estimates and phylogenetic analyses. 
#'
#' \emph{Models of Morphological Differentiation and Branching (Cladogenesis and Anagenesis)}
#'
#' These models of morphological differentiation do not involve the direct simulation of
#' morphological traits. Instead, morphotaxon identity is used as a proxy of the
#' distinctiveness of lineages on morphological grounds, as if there was some hypothetical
#' paleontologist attempting to taxonomically sort collections of specimens of these simulated
#' lineages. Two lineages are either identical, and thus share the same morphotaxon identity,
#' or they are distinct, and thus have separate morphotaxon identities. Morphological
#' differentiation is assumed to be an instantaneous process for the purposes of this model,
#' such that no intermediate could be uncovered.
#'
#' Specifically, \code{simFossilRecord} allows for three types of binary branching events
#' (here grouped under the term 'cladogenesis': 'budding cladogenesis', 'bifurcating
#' cladogenesis' and 'cryptic cladogenesis', as well as for a fourth event-type, 'anagenesis'
#' (see Wagner and Erwin, 1995; Foote, 1996, and Bapst, 2013, for further details). Budding,
#' bifurcation and cryptic cladogenetic events all share in common that a single geneological
#' lineage splits into two descendant lineages, but differ in the morphological differentiation
#' of these child lineages relative to their parent. Under budding cladogenesis, only one of the
#' child lineages becomes morphologically distinguishable from the parent, and thus the ancestral
#' morphotaxon persists through the branching event as the un-differtiated child lineage. Under
#' bifurcating cladogenesis, both child lineages become immediately distinct from the ancestor,
#' and thus two new morphotaxa appear while the ancestor terminates in an event known as
#' 'pseudoextinction'. Crytic cladogenesis has no morphological differentiation: both child lineages
#' are presumed to be indistinct from the ancestor and from each other, which means a hypothetical
#" paleontologist would not observe that branching had occurred at all. Anagenesis is morphological
#' differentiation independent of any branching, such that a morphotaxon instanteously transitions
#' to a new morphotaxon identity, resulting in the pseudoextinction of the ancestral morphotaxon
#' and the immediate 'pseudospeciation' of the child morphotaxon. The two morphotaxa do not overlap
#' in time at all, as modeled here (contra to the models described by Ezard et al., 2012). For ease
#' of following these cryptic lineages, cryptic cladogenetic events are treated in terms of data
#' structure similarly to budding cladogenetic events, with one child lineage treated as a
#' persistance of the ancestral lineage, and the other as a new morphologically indistinguishable lineage.
#' This model of cryptic cladogenesis is ultimately based on the hierarchical birth-death model used
#' by many authors for modeling patterns across paraphyletic higher taxa and the lower taxon units within
#' them (e.g. Patzkowsky, 1995; Foote, 2012).
#'
#' The occurrence of the various models is controlled by multiple arguments of \code{simFossilRecord}.
#' The overall instantaneous rate of branching (cladogenesis) is controlled by argument \code{p}, and
#' the proportion of each type of cladogenesis controlled by arguments \code{prop.bifurc} and
#' \code{prop.cryptic}. \code{prop.cryptic} controls the overall probability that any branching event
#' will be cryptic versus involving any morphological differentiation (budding or bifurcating). If
#' \code{prop.cryptic = 1}, all branching events will be cryptic cladogenesis, and if
#' \code{prop.cryptic = 0}, all branching events will involve morphological differentiation and none
#' will be cryptic. \code{prop.bifurc} controls how many branching events that involve morphological
#' differentiation (i.e. the inverse of \code{prop.cryptic}) are bifurcating, as opposed to budding
#' cladogenesis.  If \code{prop.bifurc = 1}, all morphologically-differentiating branching events will
#' be bifurcating cladogenesis, and if \code{prop.bifurc = 0}, all morphologically-differentiating
#' branching events will be budding cladogenesis. Thus, for example, the probability of a given
#' cladogenesis event being budding is given by:
#' 
#'  \code{Prob(budding cladogenesis at a branching event) = (1 - prop.cryptic) * (1 - prop.bifurc)}
#' 
#' By default,
#' \code{prop.cryptic = 0} and \code{prop.bifurc = 0}, so all branching events by default will be
#' instances of budding cladogenesis. Anagenesis is completely independent of these, controlled as its
#' own Poisson process with an instantaneous rated defined by the argument \code{anag.rate}. By default,
#' this rate is set to zero and thus there is no anagenetic events without user intervention.
#'
#' \emph{Stopping Conditions and Acceptance Criteria for Simulations}
#'
#' How forward-time simulations are generated, halted and whether they are accepted
#' or not for output is a critical component of simulation design.
#' Most uses of \code{simFossilRecord} will involve iteratively
#' generating and analyzing multiple simulation runs. Runs are only
#' accepted for output if they meet the conditioning criteria defined
#' in the arguments, either matching point constraints or falling
#' within range constraints. However, this requires separating the processes of
#' halting simulation runs and accepting a run for output, particularly to avoid bias
#' related to statistical sampling issues.
#'
#' Hartmann et al. (2011) recently discovered a potential statistical artifact
#' when branching simulations are conditioned on some number of taxa.
#' Previously within \code{paleotree}, this was accounted for in the deprecated function \code{simFossilTaxa} by
#' a complex arrangement of minimum and maximum constraints, and an (incorrect)
#' presumption that allowing simulations to continue for a short distance after
#' constraints were reached. This strategy is not applied here. Instead,
#' \code{simFossilRecord} applies the General Sampling Algorithm presented
#' by Hartmann et al. (or at least, a close variant). A simulation continues until
#' extinction or some maximum time-constraint is reached, evaluated for intervals
#' that match the set run conditions (e.g. \code{nExtant}, \code{nTotalTime}) and, if some
#' interval or set of intervals matches the run conditions, a date is randomly sampled
#' from within this interval/intervals. The simulation is then cut at this date using
#' the \code{timeSliceFossilRecord} function, and saved as an accepted run.
#' The simulation data is otherwise discarded and then a new simulation initiated
#' (thus, at most, only one simulated dataset is accepted from one simulation run).
#'
#' Thus, accepted simulations runs should reflect unbiased samples of evolutionary
#' histories that precisely match the input constraints, which can be very precise,
#' unlike how stopping and acceptance conditions were handled in the previous (deprecated)
#' \code{simFossilTaxa} function. Of course, selecting very precise constraints that
#' are very unlikely or impossible given other model parameters may take considerable
#' computation time to find acceptable simulation runs, or effectively never find any
#' acceptable simulation runs.
#' 
#' \emph{On Time-Scale Used in Output}
#' 
#' Dates given in the output are on an reversed absolute time-scale; i.e. time
#' decreases going from the past to the future, as is typical in paleontological
#' uses of time (as time before present) and as for most function in package
#' \code{paleotree}. The endpoints of the time-scale are decided by details of the
#' simulation and can be modified by several arguments. By default (with
#' \code{shiftRoot4TimeSlice =} \code{"withExtantOnly"}), any simulation run that is accepted
#' with extant taxa will have zero as the \emph{end-time} (i.e. when those taxa are
#' extant), as zero is the typical time assigned to the modern day in empirical studies.
#' If a simulation ends with all taxa extinct, however, then instead the \emph{start-time}
#' of a run (i.e. when the run initiates with starting taxa) will be maximum value
#' assigned to the conditioning argument \code{totalTime}. If \code{shiftRoot4TimeSlice =} \code{FALSE}
#' then the \emph{start-time} of the run will always be this maximum value for
#' \code{totalTime}, and any extant taxa will stop at some time greater than zero.
#'

#' @param p,q,r,anag.rate These parameters control the instantaneous ('per-capita') rates of branching, extinction,
#' sampling and anagenesis, respectively. These can be given as a number equal to or greater than zero, or as a 
#' character string which will be interpreted as an algebraic equation. These equations can make use of three
#' quantities which will/may change throughout the simulation: the standing richness is \code{N}, the
#' current time passed since the start of the simulation is \code{T}, the present duration of a given still-living
#' lineage since its origination time is code{D}, and the current branching rate is \code{P}
#' (corresponding to the argument name \code{p}).
#' Note that \code{P} cannot be used in equations for the branching rate itself; it is for making other rates
#' relative to the branching rate.
#' By default, the rates \code{r} and \code{anag.rate} are set to zero, so that the default simulator is a birth-death
#' simulator.
#' Rates set to \code{= Inf} are treated as if 0. When a rate is set to 0, this event type will not occur in the simulation.
#' Setting certain processes to zero, like sampling, may increase simulation efficiency, if the goal is a birth-death or
#' pure-birth model.
#' See documentation for argument \code{negRatesAsZero} about the treatment of rates that decrease below zero.
#' Notation of branching, extinction and sampling rates as \code{p, q, r} 
#' follows what is typical for the paleobiology literature (e.g. Foote, 1997), not the  Greek letters \code{lambda, mu, phi}
#' found more typically in the biological literature (e.g. Stadler, 2009; Heath et al., 2014; Gavryushkina et al., 2014).

#' @param totalTime,nTotalTaxa,nExtant,nSamp These arguments represent stopping and
#' acceptance conditions for simulation runs. They are respectively \code{totalTime}, the
#' total length of the simulation in time-units, \code{nTotalTaxa}, the total number of taxa
#' over the past evolutionary history of the clade, \code{nExtant}, the total number of extant taxa at
#' the end of the simulation and \code{nSamp} the total number of sampled taxa (not counting extant
#' taxa sampled at the modern day). These are used to determine when to end simulation runs, and whether to accept
#' or reject them as output. They can be input as a vector of two numbers, representing minimum
#' and maximum values of a range for accepted simulation runs (i.e. the simulation length can be between 0 and
#' 1000 time-steps, by default), or as a single number, representing a point condition (i.e. if
#' \code{nSamp = 100} then the only simulation runs with exactly 100 taxa sampled will be output).
#' Note that it is easy to set combinations of parameters and run conditions that are impossible
#' to produce satisfactory input under, in which case \code{simFossilRecord} would run in a nonstop loop.
#' How cryptic taxa are counted for the sake of these conditions is controlled by argument \code{count.cryptic}.

#' @param negRatesAsZero A logical. Should rates calculated as a negative number cause the simulation to fail
#' with an error message (\code{ = FALSE}) or should these be treated as zero (\code{"= TRUE"}, the default). This
#' is equivalent to saying that the \code{rate.as.used =} \code{max(0, rate.as.given)}.

#' @param prop.cryptic,prop.bifurc These parameters control (respectively) the proportion of branching events that have
#' morphological differentiation, versus those that are cryptic (\code{prop.cryptic}) and the proportion of morphological
#' branching events that are bifurcating, as opposed to budding. Both of these proportions must be a number between 0 and 1.
#' By default, both are set to zero, meaning all branching events are events of budding cladogenesis. See description of
#' the available models of morphological differentiation in the \emph{Description} section.

#' @param tolerance A small number which defines a tiny interval for the sake of placing run-sampling dates before events and
#' for use in determining whether a taxon is extant in simFossilRecordMethods.

#' @param maxStepTime When rates are time-dependent (i.e. when parameters 'D' or 'T' are used in equations input for one of
#' the four rate arguments), then protocol used by \code{simFossilRecord} of drawing waiting times to the next event could
#' produce a serious mismatch of resulting process to the defined model, because the possibility of new events is only
#' considered at the end of these waiting times. Instead, any time a waiting time greater than \code{maxStepTime} is
#' selected, then instead \emph{no} event occurs and a time-step equal to \code{maxStepTime} occurs instead, thus effectively
#' discretizing the progression of time in the simulations run by \code{simFossilRecord}. Decreasing this value will increase
#' accuracy (as the time-scale is effectively more discretized) but increase computation time, as the computer will need
#' to stop and check rates to see if an event happened more often. Users should toggle this value relative to the time-dependent
#' rate equations they input, relative to the rate of change in rates expected in time-dependent rates.

#' @param nruns Number of simulation datasets to accept, save and output. If \code{nruns = 1}, output will be a single
#' object of class \code{fossilRecordSimulation}, and if \code{nruns} is greater than 1, a list will be output composed of
#' \code{nruns} objects of class \code{fossilRecordSimulation}. 

#' @param startTaxa Number of initital taxa to begin a simulation with. All will have the simulation start date
#' listed as their time of origination.

#' @param sortNames If TRUE, output taxonomic lists are sorted by the taxon
#' names (thus sorting cryptic taxa together) rather than by taxon ID number
#' (i.e. the order they were simulated in).

#' @param print.runs If TRUE, prints the proportion of simulations accepted for
#' output to the terminal.

#' @param maxAttempts Number of simulation attempts allowed before the simulation process
#' is halted with an error message. Default is \code{Inf}.

#' @param plot If TRUE, plots the diversity curves of accepted simulations,
#' including both the diversity curve of the true number of taxa and the
#' diversity curve for the 'observed' (sampled) number of taxa.

#' @param count.cryptic If TRUE, cryptic taxa are counted as separate taxa for
#' conditioning limits that count a number of taxon units, such as \code{nTotalTaxa},
#' \code{nExtant} and \code{nSamp}. If FALSE (the default), then each cryptic
#' complex (i.e. each distinguishable morphotaxon) is treated as a single taxon.
#' See examples.

#' @inheritParams simFossilRecordMethods

#' @return
#' \code{simFossilRecord} returns either a single object of class \code{fossilRecordSimulation}
#' or a list of multiple such objects, depending on whether \code{nruns} was 1 or more.
#'
#' An object of class \code{fossilRecordSimulation} consists of a list object composed of multiple
#' elements, each of which is data for 'one taxon'. Each data element for each taxon is itself
#' a list, composed of two elements: the first describes vital information about the taxon unit, 
#' and the second describes the sampling times of each taxon. 
#'
#' The first element of the list (named \code{$taxa.data}) is a distinctive six-element
#' vector composed of numbers (some are nominally integers, but not all, so all are stored
#' as double-precision integers) with the following field names:
#'
#' \describe{

#' \item{\code{taxon.id}}{The ID number of this particular taxon-unit.}

#' \item{\code{ancestor.id}}{The ID number of the ancestral taxon-unit.
#' The initial taxa in a simulation will be listed with \code{NA} as their ancestor.}

#' \item{\code{orig.time}}{True time of origination for a taxon-unit in absolute time.}

#' \item{\code{ext.time}}{True time of extinction for a taxon-unit in absolute time.
#' Extant taxa will be listed with an \code{ext.time} of the run-end time of the
#' simulation run, which for simulations with extant taxa is 0 by default (but this
#' may be modified using argument \code{shiftRoot4TimeSlice}).}

#' \item{\code{still.alive}}{Indicates whether a taxon-unit is 'still alive' or not:
#' '1' indicates the taxon-unit is extant, '0' indicates the taxon-unit is extinct}

#' \item{\code{looks.like}}{The ID number of the first morphotaxon in a dataset that
#' 'looks like' this taxon-unit; i.e. belongs to the same multi-lineage cryptic
#' complex. Taxa that are morphologically distinct from any previous lineage will
#' have their \code{taxon.id} match their \code{looks.like}. Thus, this column
#' is rather uninformative unless cryptic cladogenesis occurred in a simulation.}}
#'
#' The second element for each taxon-unit is a vector of sampling times, creatively
#' named \code{$sampling.times}, with each value representing a data in absolute time
#' when that taxon was sampled in the simulated fossil record. If a taxon was never
#' sampled, this vector is an empty numeric vector of \code{length = 0}.
#'
#' As is typical for paleontological uses of absolute time, absolute time in these
#' simulations is always decreasing toward the modern; i.e. an absolute date of 50
#' means a point in time which is 50 time-units before the present-day, if the
#' present-day is zero (the default, but see argument \code{shiftRoot4TimeSlice}).
#' 
#' Each individual element of a \code{fossilRecordSimulation} list object is named, generally of
#' the form "t1" and "t2", where the number is the \code{taxon.id}. Cryptic taxa are instead
#' named in the form of "t1.2" and "t5.3", where the first number is the taxon which they are a
#' cryptic descendant of (\code{looks.like}) and the second number, after the period, is
#' the order of appearance of lineage units in that cryptic complex. For example, for
#' "t5.3",  the first number is the \code{taxon.id} and the second number communicates
#' that this is the third lineage to appear in this cryptic complex.


#' @seealso
#' \code{\link{simFossilRecordMethods}}
#'
#' This function essentially replaces and adds to all functionality of the
#' deprecated \code{paleotree} functions \code{simFossilTaxa}, \code{simFossilTaxaSRCond},
#' \code{simPaleoTrees}, as well as the combined used of \code{simFossilTaxa}
#' and \code{sampleRanges} for most models of sampling. 

#' @author 
#' David W. Bapst, inspired by code written by Peter Smits.

#' @references
#' Bapst, D. W. 2013. When Can Clades Be Potentially Resolved with
#' Morphology? \emph{PLoS ONE} 8(4):e62312.
#'
#' Ezard, T. H. G., P. N. Pearson, T. Aze, and A. Purvis. 2012. The meaning of birth
#' and death (in macroevolutionary birth-death models). \emph{Biology Letters} 8(1):139-142.
#'
#' Foote, M. 1996 On the Probability of Ancestors in the Fossil
#' Record. \emph{Paleobiology} \bold{22}(2):141--151.
#'
#' Foote, M. 1997. Estimating Taxonomic Durations and Preservation
#' Probability. \emph{Paleobiology} 23(3):278-300.
#'
#' Foote, M. 2000. Origination and extinction components of taxonomic diversity:
#' general problems. Pp. 74-102. In D. H. Erwin, and S. L. Wing, eds. \emph{Deep Time:
#' Paleobiology's Perspective.} The Paleontological Society, Lawrence, Kansas.
#'
#' Foote, M. 2012. Evolutionary dynamics of taxonomic structure. \emph{Biology Letters} 8(1):135-138.
#'
#' Gavryushkina, A., D. Welch, T. Stadler, and A. J. Drummond. 2014. Bayesian Inference
#' of Sampled Ancestor Trees for Epidemiology and Fossil Calibration. \emph{PLoS Comput Biol.}
#' 10(12):e1003919.
#'
#' Hartmann, K., D. Wong, and T. Stadler. 2010 Sampling Trees from Evolutionary
#' Models. \emph{Systematic Biology} \bold{59}(4):465--476.
#'
#' Heath, T. A., J. P. Huelsenbeck, and T. Stadler. 2014. The fossilized birth-death process
#' for coherent calibration of divergence-time estimates. \emph{Proceedings of the National Academy
#' of Sciences} 111(29):E2957-E2966.
#'
#' Kendall, D. G. 1948 On the Generalized "Birth-and-Death" Process. \emph{The
#' Annals of Mathematical Statistics} \bold{19}(1):1--15.
#' 
#' Nee, S. 2006 Birth-Death Models in Macroevolution. \emph{Annual Review of
#' Ecology, Evolution, and Systematics} \bold{37}(1):1--17.
#'
#' Patzkowsky, M. E. 1995. A Hierarchical Branching Model of Evolutionary Radiations.
#' \emph{Paleobiology} 21(4):440-460.
#' 
#' Solow, A. R., and W. Smith. 1997 On Fossil Preservation and the
#' Stratigraphic Ranges of Taxa. \emph{Paleobiology} \bold{23}(3):271--277.
#'
#' Stadler, T. 2009. On incomplete sampling under birth-death models and connections to the
#' sampling-based coalescent. \emph{Journal of Theoretical Biology} 261(1):58-66.
#'
#' Wagner, P. J., and D. H. Erwin. 1995. Phylogenetic patterns as tests of speciation models.
#' Pp. 87-122. In D. H. Erwin, and R. L. Anstey, eds. \emph{New approaches to speciation in the
#' fossil record.} Columbia University Press, New York.

#' @examples
#' 
#' set.seed(2)
#' 
#' # quick birth-death-sampling run with 1 run, 50 taxa
#' 
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=1,
#' 	nTotalTaxa=50, plot=TRUE)
#' 
#' \donttest{ 
#' # examining multiple runs of simulations
#' 
#' #example of repeated pure birth simulations over 50 time-units
#' records <- simFossilRecord(p=0.1, q=0, nruns=10,
#' 	totalTime=50, plot=TRUE)
#' #plot multiple diversity curves on a log scale
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE,plotLogRich=TRUE)
#' #histogram of total number of taxa
#' hist(sapply(records,nrow))
#' 
#' #example of repeated birth-death-sampling simulations over 50 time-units
#' records <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=10,
#' 	totalTime=50, plot=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' 
#' #like above, but conditioned instead on having 10 extant taxa
#' 	# between 1 and 100 time-units
#' set.seed(4)
#' records <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=10,
#' 	totalTime=c(1,300), nExtant=10, plot=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' 
#' ################################################
#' 
#' # How probable were the runs I accepted?
#' 	# The effect of conditions
#' 
#' set.seed(1)
#' 
#' # Let's look at an example of a birth-death process
#' 	# with high extinction relative to branching
#' # use default run conditions (barely any conditioning)
#' # use print.runs to look at acceptance probability
#' records <- simFossilRecord(p=0.1, q=0.8, nruns=10,
#' 	print.runs=TRUE, plot=TRUE)
#' # 10 runs accepted from a total of 10 !
#' 
#' # now let's give much more stringent run conditions
#' 	# require 3 extant taxa at minimum, 5 taxa total minimum
#' records <- simFossilRecord(p=0.1, q=0.8, nruns=10,
#' 	nExtant=c(3,100), nTotalTaxa=c(5,100),
#' 	print.runs=TRUE, plot=TRUE)
#' # thousands of simulations to just obtail 10 accepable runs!
#' 	# most ended in extinction before minimums were hit
#' 
#' # beware analysis of simulated where acceptance conditions 
#' 	# are too stringent: your data will be a 'special case'
#' 	# of the simulation parameters
#' # it will also take you a long time to generate reasonable
#' 	# numbers of replicates for whatever analysis you are doing
#' 
#' # TLDR: You should look at print.runs=TRUE
#' 
#' ######################
#' 
#' # Using the rate equation-input for complex diversification models
#' 
#' # First up... Diversity Dependent Models!
#' # Let's try Diversity-Dependent Branching over 50 time-units
#' 
#' # first, let's write the rate equation
#' # We'll use the diversity dependent rate equation model
#' 	# from Ettienne et al. 2012 as an example here
#' 	# Under this equation, p=q at carrying capacity K
#' # Many others are possible!
#' # Note that we don't need to use max(0,rate) as negative rates
#' 	# are converted to zero by default, as controlled by
#' 	# the argument negRatesAsZero
#' 
#' # From Ettiene et al.
#' #	lambda = lambda0 - (lambda0 - mu)*(n/K)
#' # lambda and mu are branching rate and extinction rate
#' 	# lambda and mu == p and q in paleotree (i.e. Foote convention)
#' # lambda0 is the branching rate at richness=0
#' # K is the carrying capacity
#' # n is the richness
#' 
#' # 'N' is the algebra symbol for standing taxonomic richness 
#' 	# for simFossilRecord's simulation capabilities
#' # also branching rate cannot reference extinction rate
#' # we'll have to set lambda0, mu and K in the rate equation directly
#' 
#' lambda0 <- 0.3	# branching rate at 0 richness in Ltu
#' K <- 40		# carrying capacity 
#' mu <- 0.1		# extinction rate will 0.1 Ltu (= 1/3 of lambda0)
#' 
#' # technically, mu here represents the lambda at richness=K
#' 	# i.e. lambdaK
#' # Ettienne et al. are just implicitly saying that the carrying capacity
#' 	# is the richness at which lambda==mu
#' 
#' # construct the equation programmatically using paste0
#' branchingRateEq<-paste0(lambda0,"-(",lambda0,"-",mu,")*(N/",K,")")
#' # and take a look at it...
#' branchingRateEq
#' # its a thing of beauty, folks
#' 
#' # now let's try it
#' records <- simFossilRecord(p=branchingRateEq, q=mu, nruns=3,
#' 	totalTime=100, plot=TRUE, print.runs=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' # those are some happy little diversity plateaus!
#' 
#' 
#' # now let's do diversity-dependent extinction
#' 
#' # let's slightly modify the model from Ettiene et al.
#' #	mu = mu0 + (mu0 - muK)*(n/K)
#' 
#' mu0<-0.001		# mu at n=0
#' muK<-0.1		# mu at n=K (should be equal to lambda at K)
#' K<-40
#' lambda<-muK		# equal to muK
#' 
#' # construct the equation programmatically using paste0
#' extRateEq<-paste0(mu0,"-(",mu0,"-",muK,")*(N/",K,")")
#' extRateEq
#' 
#' # now let's try it
#' records <- simFossilRecord(p=lambda, q=extRateEq, nruns=3,
#' 	totalTime=100, plot=TRUE, print.runs=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' 
#' # these plateaus looks a little more spiky 
#' 	#( maybe there is more turnover at K? )
#' # also, it took a longer for the rapid rise to occur
#' 
#' # Now let's try an example with time-dependent origination
#' 	# and extinction constrained to equal origination
#' 
#' # Note! Use of time-dependent parameters "D" and "T" may
#' # result in slower than normal simulation run times
#' # as the time-scale has to be discretized; see
#' # info for argument maxTimeStep above
#' 
#' # First, let's define a time-dependent rate equation
#' 	# "T" is the symbol for time passed
#' timeEquation<-"0.4-(0.007*T)"
#' 
#' #in this equation, 0.4 is the rate at time=0
#' 	# and it will decrease by 0.007 with every time-unit
#' 	# at time=50, the final rate will be 0.05
#' # We can easily make it so extinction is always equal to branching rate
#' # "P" is the algebraic equivalent for "branching rate" in simFossilRecord
#'
#' # now let's try it
#' records <- simFossilRecord(p=timeEquation, q="P", nruns=3,
#' 	totalTime=50, plot=TRUE, print.runs=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' # high variability that seems to then smooth out as turnover decreases
#'
#' # And duration what about duration-dependent processes?
#'		# let's do a duration-dep extinction equation:
#' durDepExt<-"0.01+(0.01*D)"
#'
#' # okay, let's take it for a spin
#' records <- simFossilRecord(p=0.1, q=durDepExt, nruns=3,
#' 	totalTime=50, plot=TRUE, print.runs=TRUE)
#' records<-lapply(records,fossilRecord2fossilTaxa)
#' multiDiv(records,plotMultCurves=TRUE)
#' # creates runs full of short lived taxa
#' 
#' ##########################################################
#' 
#' # Speciation Modes
#' # Some examples of varying the 'speciation modes' in simFossilRecord
#' 
#' # The default is pure budding cladogenesis
#' 	# anag.rate = prop.bifurc = prop.cryptic = 0
#' # let's just set those for the moment anyway
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0,
#' 	nruns=1, nTotalTaxa=c(20,30) ,nExtant=0, plot=TRUE)
#' 
#' #convert and plot phylogeny
#' 	# note this will not reflect the 'budding' pattern
#' 	# branching events will just appear like bifurcation
#' 	# its a typical convention for phylogeny plotting
#' converted<-fossilRecord2fossilTaxa(record)
#' tree<-taxa2phylo(converted,plot=TRUE)
#' 
#' #now, an example of pure bifurcation
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=1, prop.cryptic=0,
#' 	nruns=1, nTotalTaxa=c(20,30) ,nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' 
#' # all the short branches are due to ancestors that terminate
#' 	# via pseudoextinction at bifurcation events
#' 
#' # an example with anagenesis = branching
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0.1, prop.bifurc=0, prop.cryptic=0,
#' 	nruns=1, nTotalTaxa=c(20,30) ,nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' # lots of pseudoextinction
#' 
#' # an example with anagenesis, pure bifurcation
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0.1, prop.bifurc=1, prop.cryptic=0,
#' 	nruns=1, nTotalTaxa=c(20,30) ,nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' # lots and lots of pseudoextinction
#' 
#' # an example with half cryptic speciation
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0.5,
#' 	nruns=1, nTotalTaxa=c(20,30), nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' 
#' # notice that the tree has many more than the maximum of 30 tips:
#' 	# that's because the cryptic taxa are not counted as
#' 	# separate taxa by default, as controlled by count.cryptic
#' 
#' # an example with anagenesis, bifurcation, cryptic speciation
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0.1, prop.bifurc=0.5, prop.cryptic=0.5,
#' 	nruns=1, nTotalTaxa=c(20,30), nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#' # note in this case, 50% of branching is cryptic
#' 	# 25% is bifurcation, 25% is budding
#' 
#' # an example with anagenesis, pure cryptic speciation
#' 	# morphotaxon identity will thus be entirely indep of branching!
#' 	# I wonder if this is what is really going on, sometimes...
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0.1, prop.bifurc=0, prop.cryptic=1,
#' 	nruns=1, nTotalTaxa=c(20,30), nExtant=0)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record),plot=TRUE)
#'
#' # merging cryptic taxa when all speciation is cryptic
#' set.seed(1)
#' record <- simFossilRecord(p=0.1,
#' 	q=0.1, r=0.1,
#' 	prop.crypt=1,
#' 	totalTime=50, plot=TRUE)
#' # there looks like there is only a single taxon, but...
#' length(record)	#actual number of cryptic lineages
#' 
#' #############
#' 
#' # playing with count.cryptic with simulations of pure cryptic speciation
#' 
#' #can choose to condition on total morphologically-distinguishable taxa
#'     #or total taxa including cryptic taxa with count.cryptic=FALSE
#' 
#' # an example with pure cryptic speciation with count.cryptic=TRUE
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=1,
#' 	nruns=1, totalTime=50, nTotalTaxa=c(10,100), count.cryptic=TRUE)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record))
#' plot(tree);axisPhylo()
#' # notice how the tip labels indicate all are the same morphotaxon
#' 
#' # we'll replace the # of taxa constraints with a time constraint
#' 	# or else the count.cryptic=FALSE simulation will never end!
#' 
#' # an example with pure cryptic speciation with count.cryptic=FALSE
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=1,
#' 	nruns=1, totalTime=50, count.cryptic=FALSE)
#' tree<-taxa2phylo(fossilRecord2fossilTaxa(record))
#' plot(tree);axisPhylo()
#' 
#' #let's look at numbers of taxa returned when varying count.cryptic
#' 	# with prop.cryptic=0.5
#' 
#' #simple simulation going for 50 total taxa	
#' 
#' #first, count.cryptic=FALSE (default)
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0.5,
#' 	nruns=1, nTotalTaxa=50, count.cryptic=FALSE)
#' taxa<-fossilRecord2fossilTaxa(record)
#' nrow(taxa)                 		#number of lineages (inc. cryptic)
#' length(unique(taxa[,6]))            #number of morph-distinguishable taxa
#' 
#' # and count.cryptic=TRUE
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0.5,
#' 	nruns=1, nTotalTaxa=50, count.cryptic=TRUE)
#' taxa<-fossilRecord2fossilTaxa(record)
#' nrow(taxa)                 		#number of lineages (inc. cryptic)
#' length(unique(taxa[,6]))            #number of morph-distinguishable taxa
#' 
#' # okay...
#' # now let's try with 50 extant taxa
#' 
#' #first, count.cryptic=FALSE (default)
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0.5,
#' 	nruns=1, nExtant=10, totalTime=c(1,100), count.cryptic=FALSE)
#' taxa<-fossilRecord2fossilTaxa(record)
#' sum(taxa[,5])             		  	#number of still-living lineages (inc. cryptic)
#' length(unique(taxa[taxa[,5]==1,6]))	   	#number of still-living morph-dist. taxa
#' 
#' # and count.cryptic=TRUE
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1,
#' 	anag.rate=0, prop.bifurc=0, prop.cryptic=0.5,
#' 	nruns=1, nExtant=10, totalTime=c(1,100), count.cryptic=TRUE)
#' taxa<-fossilRecord2fossilTaxa(record)
#' sum(taxa[,5])             		  	#number of still-living lineages (inc. cryptic)
#' length(unique(taxa[taxa[,5]==1,6]))	   	#number of still-living morph-dist. taxa
#' 
#' #################################################
#'
#' # an example using startTaxa to have more initial taxa
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=1,
#'	nTotalTaxa=100, startTaxa=20, plot=TRUE)
#'
#' ######################################################
#' 
#' # Using run conditions
#' 
#' # Users can generate datasets that meet multiple conditions:
#' 	# such as time, number of total taxa, extant taxa, sampled taxa
#' # These can be set as point conditions or ranges
#' 
#' # let's set time = 10-100 units, total taxa = 30-40, extant = 10
#' 	#and look at acceptance rates with print.run
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=1, 
#' 	totalTime=c(10,100), nTotalTaxa=c(30,40), nExtant=10,
#' 	print.runs=TRUE, plot=TRUE)
#' 
#' # let's make the constraints on totaltaxa a little tighter
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=1, 
#' 	totalTime=c(50,100), nTotalTaxa=30, nExtant=10,
#' 	print.runs=TRUE, plot=TRUE)
#' # still okay acceptance rates
#' 
#' # alright, now let's add a constraint on sampled taxa
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=1, 
#' 	totalTime=c(50,100), nTotalTaxa=30, nExtant=10,
#' 	nSamp=15, print.runs=TRUE, plot=TRUE)
#' # still okay acceptance rates
#'
#' # we can be really odd and condition on having a single taxon
#' set.seed(1)
#' record <- simFossilRecord(p=0.1,
#'	q=0.1, r=0.1, nTotalTaxa=1,
#'	totalTime=c(10,20), plot=TRUE)
#'
#' ########################################################
#' 
#' # Simulations of entirely extinct taxa
#' 
#' #Typically, a user may want to condition on a precise
#' 	# number of sampled taxa in an all-extinct simulation
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=1, 
#' 	nTotalTaxa=c(1,100), nExtant=0, nSamp=20,
#' 	print.runs=TRUE, plot=TRUE)
#'
#' # Note that when simulations don't include
#' # sampling or extant taxa, the plot 
#' # functionality changes
#' record <- simFossilRecord(p=0.1, q=0.1, r=0, nruns=1, 
#' 	nExtant=0, print.runs=TRUE, plot=TRUE)
#' # something similar happens when there is no sampling
#' # and there are extant taxa but they aren't sampled
#' record <- simFossilRecord(p=0.1, q=0.1, r=0, nruns=1, 
#' 	nExtant=10, nTotalTaxa=100, modern.samp.prob=0,
#' 	print.runs=TRUE, plot=TRUE)
#' 
#'	
#' # We can set up a test to make sure that no extant taxa somehow get
#' # returned in many simulations with extinct-only conditioning:
#' res<-simFossilRecord(p=0.1, q=0.1, r=0.1,nTotalTaxa=10,nExtant=0,nruns=1000,plot=TRUE)
#' anyLive<-any(sapply(res,function(z) any(sapply(z,function(x) x[[1]][5]==1))))
#' if(anyLive){
#'	stop("Runs have extant taxa under conditioning for none?")
#'	}
#'
#' }
#'

#' @name simFossilRecord
#' @rdname simFossilRecord
#' @export
simFossilRecord<-function(

	# model parameters
	#
	p, q, r=0, anag.rate=0, prop.bifurc=0, prop.cryptic=0,
	modern.samp.prob=1, startTaxa=1, nruns=1, maxAttempts=Inf,

	# run conditions can be given as vectors of length 1 or length 2 (= min,max)
	#
	totalTime = c(0, 1000), nTotalTaxa = c(1, 1000),
	nExtant = c(0, 1000), nSamp = c(0, 1000),

	#control parameters
	#
	tolerance=10^-4, maxStepTime=0.01, shiftRoot4TimeSlice="withExtantOnly",
	count.cryptic=FALSE, negRatesAsZero=TRUE, print.runs=FALSE, sortNames=FALSE, plot=FALSE){

	#####################################################################################
	
	# NOT USED (but in simFossilTaxa):	min.cond=TRUE

	#################################################################################
	
	#example parameter sets
	#
	# DEFAULTS (without rates set)
		# r=0.1;anag.rate=0;prop.bifurc=0;prop.cryptic=0;startTaxa=1;nruns=1;
		# nTotalTaxa=c(1,1000);totalTime=c(1,1000);nSamp=c(0,1000);nExtant=c(0,1000);
		# tolerance=10^-4; maxStepTime=0.01; shiftRoot4TimeSlice="withExtantOnly";
		# count.cryptic=FALSE; negRatesAsZero=TRUE; print.runs=FALSE; sortNames=FALSE; plot=FALSE
	#
	# BASIC RUN	with diversity-dep extinction
	# p=0.1;q='0.01*N'
		# r=0.1;anag.rate=0;prop.bifurc=0;prop.cryptic=0;startTaxa=1;nruns=1;
		# nTotalTaxa=c(10,200);totalTime=c(1,1000);nSamp=c(0,1000);nExtant=c(0,0);plot=TRUE;
		# tolerance=10^-4; maxStepTime=0.01; shiftRoot4TimeSlice="withExtantOnly";
		# count.cryptic=FALSE; negRatesAsZero=TRUE; print.runs=FALSE; sortNames=FALSE; plot=FALSE
	#
	
	##################################################################################
	#
	#simplified birth-death-sampling simulator for fossil record
	#
	#p is the rate of branching
		#may be either budding (prob = 1-prop.bifurc) or bifurcation (prob = prop.bifurc)
	#anagenesis is a separate rate (anag.rate)
	#q is rate of extinction
	#r is rate of sampling
	#
	#ARGUMENT CHECKING
	#
	# number of starting taxa and runs must be at least 1
	if(nruns<1 | maxAttempts<1){
		stop("nruns and maxAttempts must be at least 1")}
	if(startTaxa<1){
		stop("startTaxa must be at least 1")}
	#nruns, starting taxa must be integer values
	if(!all(sapply(c(nruns,startTaxa),function(x) x==round(x)))){
			stop("nruns and startTaxa must coercible to whole number integers")}
	# check that prop.bifurc, prop.cryptic, modern.samp.prob are greater than 0 and less than 1
	if(any(c(prop.bifurc, prop.cryptic, modern.samp.prob)<0) |
			any(c(prop.bifurc, prop.cryptic, modern.samp.prob)>1)){
		stop("bad parameters input: prop.bifurc, prop.cryptic and modern.samp.prob must be between 0 and 1")}	
	# is prop.bifurc and prop.cryptic consistent?
	if(prop.bifurc>0 & prop.cryptic==1){
		stop("Prop.bifurc greater than 0 when probability of branching being cryptic is 1")}
 	#check that min nSamp isn't higher that 0, if r = 0 or Inf
	if((r==0 | is.infinite(r)) & nSamp[1]>0){
		stop("Minimum number of required sampled taxa is >0 but sampling rate is zero (or infinite)")}
	# check count.cryptic
	if(!count.cryptic){ #if false
		#check that min nTotalTaxa, nExtant, nSamp isn't higher that 1, if prop.cryptic=1
		if((prop.cryptic==1 & anag.rate==0) & ( nTotalTaxa[1]>1 | nExtant[1]>1 | nSamp[1]>1 )){
			stop(paste0("Minimum number of required  nTotalTaxa, nExtant and/or nSamp is >1 but these",
				" constraints cannot be reached as count.cryptic=FALSE, prop.cryptic=1 and anag.rate=0)"))
			}
		}
	#check that count.cryptic,negRatesAsZero,print.runs,sortNames,plot are all logicals
	if(!all(sapply(c(count.cryptic,negRatesAsZero,print.runs,sortNames,plot),is.logical))){
		stop("count.cryptic, negRatesAsZero, print.runs, sortNames, and plot arguments must be logicals")}
	#
	##################################
	# CHECK RUN CONDITIONS
	#
	# nTotalTaxa, nExtant, nSamp must all be integer values
	if(!all(sapply(c(nTotalTaxa,nExtant,nSamp),function(x) x==round(x)))){
			stop("nTotalTaxa, nExtant, nSamp must coercible to whole number integers")}		
	#
	runConditions<-list(totalTime=totalTime,nTotalTaxa=nTotalTaxa,nExtant=nExtant,nSamp=nSamp)
	#check that all are numeric
	if(any(!sapply(runConditions,is.numeric))){
		stop("Run condition arguments must be all of type numeric")
		}
	#are length of 1 or 2
	if(any(sapply(runConditions,length)>2) | any(sapply(runConditions,length)<1)){
		stop("Run condition arguments must be of length 1 or 2")
		}
	# run conditions can be given as vectors of length 1 or 2
		# i.e. a point condition or range
	# turn run conditions of length 1 into vectors of length 2
	runConditions<-lapply(runConditions,function(x)
		if(length(x)==1){c(x,x)}else{x}
		)
	#all values are over or equal to zero
	if(any(!sapply(runConditions,function(x) all(x>=0)))){
		stop("Run Condition values must be equal to or greater than 0")
		}	
	#with minimums less than maximums
	if(any(!sapply(runConditions,function(x) x[1]<=x[2]))){
		stop("Run condition misordered: values given as a range must have the minimum before the maximum")
		}	
	###########################
	#get the basic rate functions
	getBranchRate<-makeParFunct(p,isBranchRate=TRUE)
	getExtRate<-makeParFunct(q,isBranchRate=FALSE)
	getSampRate<-makeParFunct(r,isBranchRate=FALSE)
	getAnagRate<-makeParFunct(anag.rate,isBranchRate=FALSE)
	#
	# check if time-dependent simulation
	isTimeDep<-any(sapply(
		list(getBranchRate,getExtRate,getSampRate,getAnagRate)
		,attr,which="timeDep"))
	#
	##############################################
	#now iterate for nruns
	results<-list()
	ntries<-0
	for(i in 1:nruns){
		accept<-FALSE
		while(!accept){
			ntries<-ntries+1
			#test that haven't exceeded maximum number of attempts
			if(ntries>maxAttempts){
				stop(paste0("Maximum number of attempts (",maxAttempts,
					") has been exceeded with only ",i-1," run(s) successful."))
				}
			#
			#initiate the taxa dataset
			timePassed<-0
			#currentTime is the max time from runConditions
			currentTime<-runConditions$totalTime[2]
			taxa<-initiateTaxa(startTaxa=startTaxa,time=currentTime)
			#
			#get vitals
			startVitals<-getRunVitals(taxa=taxa,count.cryptic=count.cryptic)
			#start vitals table		
			vitalsRecord<-cbind(timePassed=timePassed,t(as.matrix(startVitals)))
			#test to make sure run conditions aren't impossible
			continue<-testContinue(vitals=startVitals,timePassed=timePassed,
				runConditions=runConditions)
			if(!continue){
				stop("Initial starting point already matches given run conditions")
				}
			while(continue){
				#only as long as continue=TRUE
				#
				#timePassed from the initiation of the simulation
				timePassed<-runConditions$totalTime[2]-currentTime
				#
				# get rates, sample new event, have it occur
				#
				# first get durations
				taxaDurations<-getTaxonDurations(taxa,currentTime)
				#
				#get event probability vector
				rateMatrix<-getRateMatrix(taxa=taxa, timePassed=timePassed,
					taxaDurations = taxaDurations,
					getBranchRate=getBranchRate, getExtRate=getExtRate,
					getSampRate=getSampRate, getAnagRate=getAnagRate,
					prop.cryptic=prop.cryptic, prop.bifurc=prop.bifurc,
					negRatesAsZero=negRatesAsZero)
				#
				#get the probabilty for each event in rateMatrix
				eventProbMatrix<-rateMatrix/sum(rateMatrix)
				#
				# from stackoverflow.com/questions/30232740/
					# randomly-sample-entries-of-a-matrix-and-return-the-row-column-indexes-in-r
					# simplified: arrayInd(sample(length(m),1,prob=m),dim(m)) 
				#
				#get event type and which taxon it occurs to
				sampledCell<-sample(length(rateMatrix),1,prob=eventProbMatrix)
				#will return as a 1 row matrix, of row # and col #
				sampledCell<-arrayInd(sampledCell,dim(rateMatrix))
				#
				#draw waiting time to an event (from Peter Smits)
					# exponential with rate = sum of all rates, across all taxa
				changeTime <- rexp(1, rate = sum(rateMatrix))
				#
				# if time dep rates, test if changeTime is greater than 
				if(isTimeDep & changeTime>maxStepTime){
					# redefine changeTime as maxStepTime
					changeTime <- maxStepTime
					#
					# measure time passed
					newTime<- currentTime - changeTime
					newTimePassed<-timePassed+changeTime
					#
					# obviously no event needs to occur...
				}else{
					# what is the event type
					event<-colnames(rateMatrix)[sampledCell[1,2]]
					# who did it happen to (what lineage)
					target<-attr(rateMatrix,"whichExtant")[sampledCell[1,1]]
					#
					# measure time passed
					newTime<- currentTime - changeTime
					newTimePassed<-timePassed+changeTime
					#
					# make the new event so!
					taxa<-eventOccurs(taxa=taxa,target=target,type=event,time=newTime)
					}
				#
				####################################################
				#
				#evaluate run conditions NOW for stopping
				#
				#(1) continue = TRUE until max totalTime, max nTotalTaxa, nSamp or total extinction
					# none of these can REVERSE
				#
				#get vitals
				currentVitals<-getRunVitals(taxa=taxa,count.cryptic=count.cryptic)
				# continue ??
				continue<-testContinue(vitals=currentVitals,timePassed=newTimePassed,
					runConditions=runConditions)
				#
				# Updated vitals table
					#for (2), keep a table that records changes in nTotalTaxa, nExtant, nSamp with timePassed
				#then can quickly evaluate (2)
				currentVitals<-c(timePassed=newTimePassed,t(as.matrix(currentVitals)))
				vitalsRecord<-rbind(vitalsRecord,currentVitals)
				# set new current time
				currentTime<-newTime	
				#
				###############################################################
				# some archived debugging lines for posterity
				#if(newTimePassed>74.5){browser()}
				#if(newTimePassed>120){if(taxa[[4]][[1]][4]<120){browser()}}
				#
				}
			###########################################
			#
			#accepting or rejecting runs
			#
			#discussion with Smits 05/11/15
				#real run condition is max limits / total ext for typical birth-death simulators
				#minimums are just for acceptability of runs when they hit run conditions
			#
			# NEED TO AVOID HARTMANN ET AL. EFFECT ---- simFossilTaxa did it wrong!!
				# sample simulation from intervals where it 'matched' run conditions
			#
			#(1) continue = TRUE until max totalTime, max nTotalTaxa, nSamp or total extinction
				# none of these can REVERSE
			#(2) then go back, find all interval for which run conditions were met
				# if no acceptable intervals, reject run
			#(3) randomly sample within intervals for a single date, apply timeSliceFossilRecord
			#
			###########################################
			#
			# use vitalRecords to identify intervals of acceptable parameter values
			#		
			#is it even worth checking? (were mins reached)
			worthyVitals<-worthCheckingVitalsRecord(vitalsRecord=vitalsRecord,runConditions=runConditions)
			if(worthyVitals){
				#test with testVitalsRecord to get seqVitals
				seqVitals<-testVitalsRecord(vitalsRecord=vitalsRecord,runConditions=runConditions
					,tolerance=tolerance)
				if(all(!is.na(seqVitals))){
					#hey, if its an acceptable simulation!!!!!!
					accept<-TRUE
					}
				}
			}
		#sample the sequences for a date
		passedDate<-sampleSeqVitals(seqVitals=seqVitals)
		#this date is in timePassed units: convert to backwards currentTime
		currentDate<-runConditions$totalTime[2]-passedDate
		#
		class(taxa)<-'fossilRecordSimulation'
		#
		# now time slice
			# if stop and there are extant, evaluate if sampled at modern
			# 0< modern.samp.prob <1 need to randomly sample
		taxa<-timeSliceFossilRecord(fossilRecord=taxa, sliceTime=currentDate,
			shiftRoot4TimeSlice=shiftRoot4TimeSlice, modern.samp.prob=modern.samp.prob)
		#
		##############################################################################
		# FINAL CHECKS
		# test that the produced taxa object actually passed the runConditions
		finalTest<-testFinal(taxa=taxa,timePassed=passedDate,
			runConditions=runConditions,count.cryptic=count.cryptic)
		#are there any non-identical taxa in a simulation with pure cryptic speciation?
		if(anag.rate==0 & prop.cryptic==1 & startTaxa==1){
			taxaIDsTest<-sapply(taxa,function(x) x[[1]][6])
			if(any(!sapply(taxaIDsTest,function(x) all(x==taxaIDsTest)))){
				stop("non-cryptic taxa created in a simulation with pure cryptic speciation?!")
				}
			}
		################################################################################
		#
		#name each normal taxon as t + ID 
			#cryptic taxa are cryptic id + . taxon number within that complex
		names(taxa)<-getTaxaNames(taxa=taxa)
		#sort if sortNames
		if(sortNames){
			taxa<-taxa[order(names(taxa))]
			}
		#
		results[[i]]<-taxa
		if(plot){
			divCurveFossilRecordSim(fossilRecord=taxa)
			if(nruns>1){
				title(paste0("Run Number ",i," of ",nruns))
				}
			}
		}
	if(print.runs){
		message(paste(
			nruns," runs accepted from ",ntries," total runs ("
			,signif(nruns/ntries,2)," Acceptance Probability)",sep=""))
		}
	if(nruns==1){results<-results[[1]]}
	return(results)	
	}	
