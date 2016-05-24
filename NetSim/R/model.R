# TODO: Add comment
# 
# Author: cws
###############################################################################

################################
# Create and update ModelManager

create_model_manager <- function(){
	.Call("create_model_manager", PACKAGE = "NetSim")
}

add_time_model <- function(modelManager, timeModel){
	.Call("add_time_model", modelManager, timeModel, PACKAGE = "NetSim")
}

add_change_model <- function(modelManager, timeModel, changeModel){
	.Call("add_change_model", modelManager, timeModel, changeModel, package = "NetSim")
	
}

add_updater <- function(modelManager, changeModel, updater){
	.Call("add_updater", modelManager, changeModel, updater, package = "NetSim")
}

# RcppExport SEXP add_time_updater(SEXP modelManager, SEXP timeUpdater);
add_time_updater <- function(modelManager, timeUpdater){
	.Call("add_time_updater", modelManager, timeUpdater, PACKAGE = "NetSim")
}

##############################
# Create particular models

# time models

create_poisson_model <- function(param = 1){
	.Call("create_poisson_model", param, PACKAGE = "NetSim")
}

# RcppExport SEXP create_attribute_poisson_model(SEXP attributeIndex);
create_attribute_poisson_model <- function(attributeIndex){
	.Call("create_attribute_poisson_model", attributeIndex, PACKAGE = "NetSim")
}

# change models

# generic function is commented out.. re-activate in later versions?

#create_change_model <- function(type, ...){
#	UseMethod("create_change_model")
#}

# simple model factory
#create_change_model.character <- function(name, ...){
#	
#	# TODO make this working for all change models
#	# and make factories also for time models and all updaters
#	
#	if (name == "jacksonRogers"){
#		type <- structure(name, class="jacksonRogersChangeModel")
#	}
#	else if (name == "netSaom"){
#		type <- structure(name, class="networkMChoice")
#	}
#	else{
#		stop(paste("Unknown change model: ", name, sep=""))
#	}
#	create_change_model(type, ...)
#}

#RcppExport SEXP create_jackson_rogers_change_model(
#		SEXP networkName, SEXP pLinkToParentNode, SEXP pLinkToNeighborNode,
#		SEXP nParentNodes, SEXP nNeighborNodes);
#
#create_change_model.jacksonRogersChangeModel <- function(
#		type,
#		networkId,
#		pLinkToParentNode = 1.0,
#		pLinkToNeigborNode = 1.0,
#		nParentNodes = 1,
#		nNeighborNodes = 1){
#	.Call("create_jackson_rogers_change_model", networkId, pLinkToParentNode, 
#			pLinkToNeigborNode, nParentNodes, nNeighborNodes, PACKAGE = "NetSim")
#	
#}

create_jackson_rogers_change_model <- function(networkIndex, pLinkToParentNode = 1.0, 
		pLinkToNeigborNode = 1.0, nParentNodes = 1, nNeighborNodes = 1){
	.Call("create_jackson_rogers_change_model", networkIndex, pLinkToParentNode, 
			pLinkToNeigborNode, nParentNodes, nNeighborNodes, PACKAGE = "NetSim")
	
}

# RcppExport SEXP create_watts_strogatz_change_model(SEXP networkId);
create_watts_strogatz_change_model <- function(networkIndex){
	.Call("create_watts_strogatz_change_model", networkIndex, PACKAGE = "NetSim")
}

# RcppExport SEXP create_rewire_tie_updater(SEXP networkId_)
create_rewire_tie_updater <- function(networkIndex){
	.Call("create_rewire_tie_updater", networkIndex, PACKAGE = "NetSim");
}

#RcppExport SEXP create_round_based_time_model(
#		SEXP timerIndex, SEXP intervalLength, SEXP startTime);
#
create_round_based_time_model <- function(
		timerIndex, intervalLength = 1.0, startTime = 0.0){
	
	.Call("create_round_based_time_model", timerIndex, 
			intervalLength, startTime, PACKAGE = "NetSim")
}

#RcppExport SEXP create_add_actor_updater();
create_add_actor_updater <- function(){
	.Call("create_add_actor_updater", package = "NetSim")
}

#RcppExport SEXP create_add_ties_from_newborn_actor_updater(
#		SEXP networkIndex);
#};
create_add_ties_from_newborn_actor_updater <- function(
		networkIndex){
	.Call("create_add_ties_from_newborn_actor_updater", networkIndex, 
			PACKAGE = "NetSim")
			
}

# RcppExport SEXP create_set_attribute_of_newborn_actor_updater(
# 		SEXP attributeIndex, SEXP value);
create_set_attribute_of_newborn_actor_updater <- function(
		attributeIndex, value){
	.Call("create_set_attribute_of_newborn_actor_updater", 
			attributeIndex, value, PACKAGE = "NetSim")
}

#RcppExport SEXP create_timer_updater(SEXP timerIndex);
create_timer_updater <- function(timerIndex){
	.Call("create_timer_updater", timerIndex, PACKAGE = "NetSim")
}
		
create_effect_container <- function(){
	.Call("create_effect_container", PACKAGE = "NetSim")
}

add_effect <- function(changeModel, ...){
	UseMethod("add_effect")
}

add_effect.AttributeMultinomialChoiceNetworkChangeModel <- function(changeModel, effect, attributeIndex, ...){
	add_effect_with_index(changeModel, effect, attributeIndex, ...)	
}

add_effect.MultinomialChoiceNetworkChoiceChangeModel <- function(changeModel, ...){
	add_effect_with_parameter(changeModel, ...)
}

#RcppExport SEXP add_effect_with_parameter(SEXP saom, SEXP effect, SEXP parameter);
add_effect_with_parameter <- function(saom, effect, parameter){
	.Call("add_effect_with_parameter", saom, effect, parameter, PACKAGE = "NetSim")
}


#RcppExport SEXP add_effect_with_index(SEXP saom, SEXP effect, SEXP index);
add_effect_with_index <- function(saom, effect, attributeIndex){
	.Call("add_effect_with_index", saom, effect, attributeIndex, PACKAGE = "NetSim")
}

# Generic function to create effects
create_effect <- function(name, ...){
	UseMethod("create_effect")
}

create_effect.character <- function(name, ...){
	typedName = get_effect_type(name);
	create_effect(typedName, ...)
	
}

create_effect.OneModeNetworkEffect <- function(name, networkIndex, ...){
	.Call("create_one_mode_effect", name, networkIndex, PACKAGE = "NetSim")
}


create_effect.AttributeOneModeNetworkEffect <- function(name, attributeIndex, networkIndex, ...){
	.Call("create_attribute_one_mode_effect", name, attributeIndex, networkIndex, PACKAGE = "NetSim")
}

create_effect.SimilarityAttributeOneModeNetworkEffect <- function(name, attributeIndex, networkIndex, meanSimilarityScore, ...){
	.Call("create_similarity_attribute_one_mode_effect", name, attributeIndex, networkIndex, meanSimilarityScore, PACKAGE = "NetSim")
}

create_effect.AttributeEffect <- function(name, attributeIndex, ...){
	.Call("create_attribute_effect", name, attributeIndex, PACKAGE = "NetSim")
}

create_effect.MultiplexNetworkEffect <- function(name, 
		networkIndex1, 
		networkIndex2, ...){
	.Call("create_multiplex_network_effect", name, 
			networkIndex1, networkIndex2, PACKAGE = "NetSim")
}


## error
create_effect.UnknownType <- function(name, ...){
	print(paste("Unknown effect type: ", name, sep=""))
}

get_effect_type <- function(name){
	type <- .Call("get_effect_type", name, PACKAGE = "NetSim")
	structure(name, class=type)
}

add_to_effect_container<- function(effectContainer, effect, parameter){
	.Call("add_to_effect_container", effectContainer, effect, parameter, PACKAGE = "NetSim")
}

# by default implemented as using the tie swap updater
create_multinomial_choice_network_change_model <- function(
		focalActorIndex, networkIndex, effectContainer){
	updater <- create_tie_swap_updater(networkIndex);
	.Call("create_multinomial_choice_network_change_model",
			focalActorIndex, networkIndex, effectContainer, updater,
			PACKAGE = "NetSim")		
}

create_multinomial_choice_behavior_change_model <- function(
		focalActorIndex, attributeIndex, effectContainer){
	.Call("create_multinomial_choice_behavior_change_model", 
			focalActorIndex, attributeIndex, effectContainer,
			PACKAGE = "NetSim")
}

# RcppExport SEXP create_attribute_multinomial_choice_network_change_model(
# 		SEXP networkIndex, SEXP poissonAttributeIndex, SEXP updater);
create_attribute_multinomial_choice_network_change_model <- function(
		networkIndex, poissonAttributeIndex, updater = create_tie_swap_updater(networkIndex)){
	.Call("create_attribute_multinomial_choice_network_change_model",
			networkIndex, poissonAttributeIndex, updater, PACKAGE = "NetSim")
}

# only for OneModeNetworkEffects
create_siena_model_manager <- function(poissonParameter, dependentNetworkIndex, 
		effectNames, effectInitParameters, effectParameters, nActors){
	
	modelManager <- create_model_manager()	
	
	effectContainer <- create_effect_container()
	
	# create effects
	for (i in c(1 : length(effectNames))){
		effectContainer <<- add_to_effect_container(
				effectContainer,
				create_effect(effectNames[i],
						effectInitParameters[i]),
				effectParameters[i]
		)
	}
	
	# create individual models
	for (i in c(0 : (nActors - 1))){
		# Poisson model
		poissonParameter = poissonParameter
		poissonModel <- create_poisson_model(poissonParameter)
		
		#saom
		saomModel <- create_multinomial_choice_network_change_model(
				i,
				dependentNetworkIndex,
				effectContainer
		)
		#tie updater
		tieSwapUpdater <- create_tie_swap_updater(dependentNetworkIndex)
		
		modelManager <<- add_time_model(modelManager, poissonModel)
		modelManager <<- add_change_model(modelManager, poissonModel, saomModel)
		modelManager <<- add_updater(modelManager, saomModel, tieSwapUpdater)
		
	}
	
	return(modelManager)
	
}

# updaters

create_tie_swap_updater <- function(networkIndex){
	.Call("create_tie_swap_updater", networkIndex, PACKAGE = "NetSim")
}

create_actor_attribute_set_updater <- function(attributeIndex, actorIndex){
	.Call("create_actor_attribute_set_updater", attributeIndex, actorIndex, PACKAGE = "NetSim")
}