

.womEffectsTreeCols <- function(name)
{
COLUMNS <- c(effectname=0,nlevels=1,rownumber=2,visible=3,children=4)
as.integer(COLUMNS[name])
}


.womCovTreeCols <- function(name)
{
COLUMNS <- c(effectname=0,parameter=1,rownumber=2,nlevels=3,visible=4,children=5)
as.integer(COLUMNS[name])
}


.womDefinedModelsTreeCols <- function(name)
{
COLUMNS <- c(name=0,hasResults=1,toBeAnalyzed=2,parsOnK=3,parsOnZ=4,parsOnG=5,
				NCovGroups=6,muKMean=7,muKSD=8,muZMean=9,muZSD=10,muGMean=11,muGSD=12,
				invWishartScalar=13,invGammaA0=14,invGammaB0=15,
				iterations=16,epsUpp=17,epsLow=18,leapfrog=19,useMet=20,metropScale=21,metropThin=22,
				accRate=23,burnin=24,DIC=25,pD=26,timeAnalyzed=27,foreground=28,effectiveIterations=29,logLikePostMean=30,accColor=31,DICtext=32)
as.integer(COLUMNS[name])
}



.womActiveModelTab<-function(status=TRUE)
{

	theWidget("modelsPageBox")$setSensitive(status)

}

.womInitializeModelsPageColumns <- function()
{
	.womCreateEffectSelectionColumns()
	.womCreateAllModelColumns()
	.womCreateCovColumns()
	.womCreateModelsDefinedModelsColumns()

}

.womCreateEffectSelectionColumns <- function()
{
	model <- gtkTreeStoreNew("gchararray","gchararray","gchararray","gboolean")
	treeview <- theWidget("modelsEffectsTreeview")

	gtkTreeViewSetModel(treeview, model)

	#Create effects tree view
	if(is.null(gtkTreeViewGetColumn(treeview,0))){
		
		# effect name
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Effect", renderer, 
  								text = .womEffectsTreeCols("effectname"),
								visible = .womEffectsTreeCols("visible"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# effect nlevels
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Levels", renderer, 
  								text = .womEffectsTreeCols("nlevels"),
								visible = .womEffectsTreeCols("visible"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
	}

}

.womCreateModelsDefinedModelsColumns <- function()
{
	
	model <- gtkTreeStoreNew("gchararray",
						 	"gboolean", 
						 	"gboolean", 
						 	"gint", 
						 	"gint",
						 	"gint",
							"gint", 
							"gdouble", 
							"gdouble", 
							"gdouble", 
							"gdouble", 
							"gdouble", 
							"gdouble",
							"gdouble",
							"gdouble",
							"gdouble",
							"gint",
							"gchararray",
							"gchararray",
							"gint",
							"gboolean",
							"gdouble",
							"gint",
							"gdouble",
							"gint",
							"gdouble",
							"gdouble",
							"gchararray",
							"gchararray",
							"gint",
							"gdouble",
							"gchararray",
							"gchararray")
								
	modelsTreeview <- theWidget("modelsDefinedModelsTreeview")
	analysisTreeview <- theWidget("analysisDefinedModelsTreeview")
	diagnosticsTreeview <-theWidget("diagnosticDefinedModelsTreeview")
	resultsTreeview <-theWidget("resultsDefinedModelsTreeview")
	
	resultsSortModel <- gtkTreeModelSortNewWithModel(child.model = model)
	diagnosticsSortModel <- gtkTreeModelSortNewWithModel(child.model = model)
	analysisSortModel <- gtkTreeModelSortNewWithModel(child.model = model)


	gtkTreeViewSetModel(modelsTreeview, model)	
	gtkTreeViewSetModel(analysisTreeview, analysisSortModel)	
	gtkTreeViewSetModel(diagnosticsTreeview, diagnosticsSortModel)	
	gtkTreeViewSetModel(resultsTreeview, resultsSortModel)	
	
	treeview = modelsTreeview
	
	#Create tree view
	if(is.null(gtkTreeViewGetColumn(treeview,0))){
		
		# name
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Name", renderer, 
  								text = .womDefinedModelsTreeCols("name"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# hasResults
		renderer <- gtkCellRendererToggleNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Results?", renderer, 
  								active = .womDefinedModelsTreeCols("hasResults")
								)
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(FALSE)
		
		# params on K
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "K params", 			renderer, text = .womDefinedModelsTreeCols("parsOnK"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# params on Z
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Z params", renderer, 
  								text = .womDefinedModelsTreeCols("parsOnZ"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# params on G
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "G params", renderer, 
  								text = .womDefinedModelsTreeCols("parsOnG"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# num Cov Mats
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Cov. Matrices", renderer, 
  								text = .womDefinedModelsTreeCols("NCovGroups"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
	
		# Prior on K: mean
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "muK Prior Mean", renderer, 
  								text = .womDefinedModelsTreeCols("muKMean"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# Prior on K: SD
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "muK Prior SD", renderer, 
  								text = .womDefinedModelsTreeCols("muKSD"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
	
		# Prior on Z: mean
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "muZ Prior Mean", renderer, 
  								text = .womDefinedModelsTreeCols("muZMean"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# Prior on Z: SD
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "muZ Prior SD", renderer, 
  								text = .womDefinedModelsTreeCols("muZSD"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
	
		# Prior on G: mean
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "muG Prior Mean", renderer, 
  								text = .womDefinedModelsTreeCols("muGMean"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# Prior on G: SD
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "muG Prior SD", renderer, 
  								text = .womDefinedModelsTreeCols("muGSD"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
	
		# Prior on variances: a0
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Var. Prior a0", renderer, 
  								text = .womDefinedModelsTreeCols("invGammaA0"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
	
		# Prior on variances: b0
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Var. Prior b0", renderer, 
  								text = .womDefinedModelsTreeCols("invGammaB0"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
	
		# Prior on covariances: inv Wishart Scalar
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "inv. Wishart Prior Scalar", renderer, 
  								text = .womDefinedModelsTreeCols("invWishartScalar"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
	
	}

}


.womInitialModelSetup<-function()
{
	model <- gtkTreeViewGetModel(theWidget("dataTreeview2"))
	dataset <- wommbatAnalysis$data
	
	# get the selected columns and types
	selectedCols <- .womTreeModelGetNthCol(model,.womColumnsOfInterestTreeCols("columnname"))
	selectedColsType <- .womTreeModelGetNthCol(model,.womColumnsOfInterestTreeCols("type"))	
	respText <- wommbatAnalysis$response
	changeText <- wommbatAnalysis$change 
	setsizeText <- wommbatAnalysis$setsize 
	
	# format them properly
	selectedColsType[selectedColsType=="continuous"] <- "CONTIN"
	selectedColsType[selectedColsType=="categorical"] <- "CATEG"
	SelCols <- wommbatAnalysis$SelCols <- cbind(selectedCols,selectedColsType)

	# Remove all unnecessary columns
	newDat <- wommbatAnalysis$newDat <- .womCreateDat(dataset,respText,changeText,setsizeText,SelCols)
    
	# Get all available effects
	mods <- wommbatAnalysis$mods <- .womListModels(newDat,SelCols)
	
	# Convert all effects into a nice, readable list
	effs <- wommbatAnalysis$effs <- .womNiceListEffects(mods)
	
	.womPopulateEffectsTable(effs)
}

.womPopulateEffectsTable <- function(effs)
{
	treeview <- theWidget("modelsEffectsTreeview")
	model <-  gtkTreeViewGetModel(treeview)
	toplevel <- list()
	allchildren <- list()
	
	maxway <- max(effs[,1])
	
	groupnames <- paste(1:maxway,"way")
	groupnames[1] <- "Main effects" 
	
	for(i in 1:maxway)
	{
		allchildren[[i]] <- list()
		myEffects <- which(effs[,1]==i)
		for(j in 1:length(myEffects))
		{
			allchildren[[i]][[j]] <- list(effs[myEffects[j],4],effs[myEffects[j],3],myEffects[j],TRUE,NULL)
		}
		toplevel[[i]] <- list(groupnames[i],"","",TRUE,allchildren[[i]])
	}
	
	wommbatAnalysis$toplevel = toplevel
	wommbatAnalysis$allchildren = allchildren
	
	sapply(toplevel, function(all)
    {
      group <- all[[.womEffectsTreeCols("children")+1]]

      iter <- model$append(NULL)$iter
      model$set(iter,
	  		  .womEffectsTreeCols("effectname"), all[[1]],
			  .womEffectsTreeCols("nlevels"), all[[2]],
			  .womEffectsTreeCols("rownumber"), all[[3]],
			  .womEffectsTreeCols("visible"), all[[4]])

      # add children
      sapply(group, function(effect)
	{
	  child.iter <- model$append(iter)$iter
	  model$set(child.iter,
			  .womEffectsTreeCols("effectname"), effect[[1]],
			  .womEffectsTreeCols("nlevels"), as.character(as.integer(effect[[2]])),
			  .womEffectsTreeCols("rownumber"), as.character(as.integer(effect[[3]])),
			  .womEffectsTreeCols("visible"), effect[[4]])
	})
    })
	gtkTreeViewExpandAll(treeview)
}

.womCreateAllModelColumns <- function()
{
	.womCreateModelColumns("K")
	.womCreateModelColumns("G")
	.womCreateModelColumns("Z")
}


.womCreateModelColumns <- function(whichPar="K")
{
	model <- gtkTreeStoreNew("gchararray","gchararray","gchararray","gboolean")
	treeview <- theWidget(paste("modelsOn",whichPar,"Treeview",sep=""))

	gtkTreeViewSetModel(treeview, model)

	#Create effects tree view
	if(is.null(gtkTreeViewGetColumn(treeview,0))){
		
		# effect name
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Effect", renderer, 
  								text = .womEffectsTreeCols("effectname"),
								visible = .womEffectsTreeCols("visible"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# effect nlevels
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Levels", renderer, 
  								text = .womEffectsTreeCols("nlevels"),
								visible = .womEffectsTreeCols("visible"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
	}

}

.womCreateCovColumns <- function()
{
	model <- gtkTreeStoreNew("gchararray","gchararray","gchararray","gchararray","gboolean")
	treeview <- theWidget("modelsCovTreeview")

	gtkTreeViewSetModel(treeview, model)

	#Create effects tree view
	if(is.null(gtkTreeViewGetColumn(treeview,0))){
		
		# effect name
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Effect", renderer, 
  								text = .womCovTreeCols("effectname"),
								visible = .womCovTreeCols("visible"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
				
	}

}


.womGetEffectsSelection<-function()
{
	treeview = theWidget("modelsEffectsTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	
	

	return(list(iter=iter))
}

.womGetModelOnParSelection<-function(whichPar="K")
{
	treeview = theWidget(paste("modelsOn",whichPar,"Treeview",sep=""))
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	


	return(list(iter=iter))
}

.womGetCovSelection<-function()
{
	treeview = theWidget("modelsCovTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	
	

	return(list(iter=iter))
}

.womAddToPar = function(whichPar="K")
{
	treeview1 <- theWidget("modelsEffectsTreeview")
	treemodel1 <- gtkTreeViewGetModel(treeview1)
	
	treeview2 <- theWidget(paste("modelsOn",whichPar,"Treeview",sep=""))
	treemodel2 <- gtkTreeViewGetModel(treeview2)
	
	iter = .womGetEffectsSelection()$iter
	
	
	selectednLevels = gtkTreeModelGetValue(treemodel1,iter,.womEffectsTreeCols("nlevels"))$value
	selectedEffect = gtkTreeModelGetValue(treemodel1,iter,.womEffectsTreeCols("effectname"))$value
	selectedRow = gtkTreeModelGetValue(treemodel1,iter,.womEffectsTreeCols("rownumber"))$value
	
	if(is.null(selectednLevels)) return(0)
	if(selectednLevels=="") return(0)
	
	alreadyAdded = .womTreeModelGetNthCol(treemodel2,n=.womEffectsTreeCols("effectname"))

	if(selectedEffect %in% alreadyAdded) return(0)
	
	iter <- treemodel2$append(NULL)$iter
	treemodel2$set(iter, 
			  .womEffectsTreeCols("effectname"), selectedEffect,
			  .womEffectsTreeCols("nlevels"), selectednLevels,
			  .womEffectsTreeCols("rownumber"), selectedRow,
			  .womEffectsTreeCols("visible"), TRUE
			  )

	StateEnv$win$present()
	
}

.clicked_models_addToK <- function(button)
{
	.womAddToPar(whichPar="K")
}

.clicked_models_removeFromK <- function(button)
{
	treeview = theWidget("modelsOnKTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

    if(!is.null(iter)){
		model$remove(iter)
		StateEnv$win$present()
	}	
}

.clicked_models_addToZ <- function(button)
{
	.womAddToPar(whichPar="Z")
}

.clicked_models_removeFromZ <- function(button)
{
	treeview = theWidget("modelsOnZTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

    if(!is.null(iter)){
		model$remove(iter)
		StateEnv$win$present()
	}	
}

.clicked_models_addToG <- function(button)
{
	.womAddToPar(whichPar="G")
}

.clicked_models_removeFromG <- function(button)
{
	treeview = theWidget("modelsOnGTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

    if(!is.null(iter)){
		model$remove(iter)
		StateEnv$win$present()
	}	
}

.toggle_models_useZ <- function(toggle)
{

	checkbox = theWidget("modelsUseZCheckbox")
	toggled = gtkToggleButtonGetActive(checkbox)
	
	if(toggled){
		.womTurnOnZ()
	}else{
		.womTurnOffZ()
	}
}

.womTurnOnZ <-function()
{
	# activate buttons, treeview
	theWidget("modelsAddZButton")$setSensitive(TRUE)
	theWidget("modelsRemoveZButton")$setSensitive(TRUE)
	theWidget("modelsAddZCovButton")$setSensitive(TRUE)
	theWidget("modelsNewZCovButton")$setSensitive(TRUE)
	theWidget("modelsOnZTreeview")$setSensitive(TRUE)
	
	# set forground color on cov treeview
	.womSetZCovActive(TRUE)
	
	# change priors
	theWidget("modelsMuZMeanEntry")$setText(as.character(StateEnv$PrevPriorMuZMean))
	theWidget("modelsMuZSDEntry")$setText(as.character(StateEnv$PrevPriorMuZSD))
	theWidget("modelsMuZMeanEntry")$setSensitive(TRUE)
	theWidget("modelsMuZSDEntry")$setSensitive(TRUE)
}

.womTurnOffZ <- function()
{
	# deactivate buttons, treeview
	theWidget("modelsAddZButton")$setSensitive(FALSE)
	theWidget("modelsRemoveZButton")$setSensitive(FALSE)
	theWidget("modelsAddZCovButton")$setSensitive(FALSE)
	theWidget("modelsNewZCovButton")$setSensitive(FALSE)
	theWidget("modelsOnZTreeview")$setSensitive(FALSE)
	
	# set forground color on cov treeview
	.womSetZCovActive(FALSE)
	
	# change priors
	StateEnv$PrevPriorMuZMean <- theWidget("modelsMuZMeanEntry")$getText()
	StateEnv$PrevPriorMuZSD <- theWidget("modelsMuZSDEntry")$getText()
	theWidget("modelsMuZMeanEntry")$setText("3")
	theWidget("modelsMuZSDEntry")$setText(".05")
	theWidget("modelsMuZMeanEntry")$setSensitive(FALSE)
	theWidget("modelsMuZSDEntry")$setSensitive(FALSE)
	
}

.womSetZCovActive<-function(status, remove=FALSE)
{
	treeview <- theWidget("modelsCovTreeview")
	model <- gtkTreeViewGetModel(treeview)
	
	gtkTreeModelForeach(model,function(model,path,iter,...){
		par = gtkTreeModelGetValue(model,iter,.womCovTreeCols("parameter"))$value
		if(par=="Z")
		{	
			model$set(iter,.womCovTreeCols("visible"),status)
			parentIter = gtkTreeModelIterParent(model, iter)
			if(remove){ 
				model$remove(iter)			
				if(!gtkTreeModelIterHasChild(model, parentIter$iter)) model$remove(parentIter$iter)			
			}
		}
	})

}


.womAddCovElement <- function(whichPar="K", type="new")
{
	treeview1 <- theWidget(paste("modelsOn",whichPar,"Treeview",sep=""))
	treemodel1 <-gtkTreeViewGetModel(treeview1)
	selection1 <- .womGetModelOnParSelection(whichPar)
	iter1 <- selection1$iter
	
	
	treeview2 <- theWidget("modelsCovTreeview")
	treemodel2 <-gtkTreeViewGetModel(treeview2)
	selection2 <- .womGetCovSelection()
	iter2 <- selection2$iter
	
	
	selectednLevels = gtkTreeModelGetValue(treemodel1,iter1,.womEffectsTreeCols("nlevels"))$value
	selectedEffect = gtkTreeModelGetValue(treemodel1,iter1,.womEffectsTreeCols("effectname"))$value
	selectedRow = gtkTreeModelGetValue(treemodel1,iter1,.womEffectsTreeCols("rownumber"))$value
	
	if(is.null(selectednLevels)) return(0)
	
	name <- paste(selectedEffect,"on",whichPar)
	
	alreadyAdded = NULL
	gtkTreeModelForeach(treemodel2,function(model,path,iter,...){
		name = gtkTreeModelGetValue(model,iter,.womCovTreeCols("effectname"))$value
		alreadyAdded <<- c(alreadyAdded,name)
		return(FALSE)
	})
		
	if(name %in% alreadyAdded) return(0)

	
	if(type=="new")
	{
		groupName <- paste(selectednLevels,"levels")		
		iter2 <- treemodel2$append(NULL)$iter
		treemodel2$set(iter2, 
			  .womCovTreeCols("effectname"), groupName,
			  .womCovTreeCols("parameter"), "",
			  .womCovTreeCols("nlevels"), selectednLevels,
			  .womCovTreeCols("rownumber"), "",
			  .womCovTreeCols("visible"), TRUE
			  )
		
		child.iter <- treemodel2$append(iter2)$iter
		treemodel2$set(child.iter, 
			  .womCovTreeCols("effectname"), name,
			  .womCovTreeCols("parameter"), whichPar,
			  .womCovTreeCols("nlevels"), selectednLevels,
			  .womCovTreeCols("rownumber"), selectedRow,
			  .womCovTreeCols("visible"), TRUE
			  )	
	}
	
	if(type=="add")
	{
		currGroupLevels = gtkTreeModelGetValue(treemodel2,iter2,.womCovTreeCols("nlevels"))$value
		if(is.null(currGroupLevels)) return(0)
		if(currGroupLevels != selectednLevels) return(0)
		parentIter = gtkTreeModelIterParent(treemodel2, iter2)
		if(parentIter$retval){
			child.iter <- treemodel2$append(parentIter$iter)$iter
		}else{
			child.iter <- treemodel2$append(iter2)$iter
		}
		treemodel2$set(child.iter, 
			  .womCovTreeCols("effectname"), name,
			  .womCovTreeCols("parameter"), whichPar,
			  .womCovTreeCols("nlevels"), selectednLevels,
			  .womCovTreeCols("rownumber"), selectedRow,
			  .womCovTreeCols("visible"), TRUE
			  )
	}
	gtkTreeViewExpandAll(treeview2)
	StateEnv$win$present()
}

.clicked_models_removeCov <-function(button)
{
	
	treeview = theWidget("modelsCovTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

	if(!is.null(iter)){
		parentIter = gtkTreeModelIterParent(model, iter)
		if(parentIter$retval){
			model$remove(iter)
			if(!gtkTreeModelIterHasChild(model, parentIter$iter)) model$remove(parentIter$iter)
		}else{
			model$remove(iter)
		}
	
		
		StateEnv$win$present()
	}
}



.clicked_models_newCovK <- function(button)
{
	.womAddCovElement(whichPar = "K", type="new")
}

.clicked_models_addCovK <-function(button)
{
	.womAddCovElement(whichPar = "K", type="add")
}


.clicked_models_newCovZ <- function(button)
{
	.womAddCovElement(whichPar = "Z", type="new")
}

.clicked_models_addCovZ <-function(button)
{
	.womAddCovElement(whichPar = "Z", type="add")
}


.clicked_models_newCovG <- function(button)
{
	.womAddCovElement(whichPar = "G", type="new")
}

.clicked_models_addCovG <-function(button)
{
	.womAddCovElement(whichPar = "G", type="add")
}

.womRemoveMissingCovarianceItems <- function()
{
	counter = 0
		
	treeviewCov <- theWidget("modelsCovTreeview")
	treemodelCov <- gtkTreeViewGetModel(treeviewCov)
	treeviewK <- theWidget("modelsOnKTreeview")
	treemodelK <- gtkTreeViewGetModel(treeviewK)
	treeviewZ <- theWidget("modelsOnZTreeview")
	treemodelZ <- gtkTreeViewGetModel(treeviewZ)
	treeviewG <- theWidget("modelsOnGTreeview")
	treemodelG <- gtkTreeViewGetModel(treeviewG)
	
	covNgroups <- gtkTreeModelIterNChildren(treemodelCov,iter=NULL)

	if(covNgroups > 0){
		for(i in 1:covNgroups){
			iter <- gtkTreeModelIterNthChild(treemodelCov, parent = NULL, i-1)$iter
			covNelements <- gtkTreeModelIterNChildren(treemodelCov,iter=iter)
			child.iters = list()
			for(j in 1:covNelements)
			{
				child.iters[[j]] <- gtkTreeModelIterNthChild(treemodelCov, parent = iter, j-1)$iter
			}	
			for(j in 1:covNelements)
			{
				parameter = gtkTreeModelGetValue(treemodelCov,child.iters[[j]],.womCovTreeCols("parameter"))$value
				effectname = gtkTreeModelGetValue(treemodelCov,child.iters[[j]],.womCovTreeCols("effectname"))$value
				effectname = substr(effectname,1,nchar(effectname)-5)
				listedEffects = .womTreeModelGetNthCol(theWidget(paste("modelsOn",parameter,"Treeview",sep=""))$getModel(),n=.womCovTreeCols("effectname"))
				if(!(effectname%in%listedEffects))
				{		
					counter = counter + 1
					treemodelCov$remove(child.iters[[j]])
				}
			}
		}
	}
	
return(counter)

}


.womCheckModelSettings <- function()
{
	# Inverse Gamma:
	IGa0 <- try(as.numeric(theWidget("modelsInvGammaEntry1")$getText()),silent=TRUE)
	IGb0 <- try(as.numeric(theWidget("modelsInvGammaEntry2")$getText()))
	
	# Prior mean/variance on grand means
	MuKMean <- try(as.numeric(theWidget("modelsMuKMeanEntry")$getText()),silent=TRUE)
	MuKSD <- try(as.numeric(theWidget("modelsMuKSDEntry")$getText()),silent=TRUE)
	MuZMean <- try(as.numeric(theWidget("modelsMuZMeanEntry")$getText()),silent=TRUE)
	MuZSD <- try(as.numeric(theWidget("modelsMuZSDEntry")$getText()),silent=TRUE)
	MuGMean <- try(as.numeric(theWidget("modelsMuGMeanEntry")$getText()),silent=TRUE)
	MuGSD <- try(as.numeric(theWidget("modelsMuGSDEntry")$getText()),silent=TRUE)
	
	# Covariance
	invWishartScalar <- try(as.numeric(theWidget("modelsInvWishartEntry")$getText()),silent=TRUE)
	treeviewCov <- theWidget("modelsCovTreeview")
	treemodelCov <- gtkTreeViewGetModel(treeviewCov)
	
	# Model name
	modelName <- theWidget("modelsNameEntry")$getText()
	
	#N effects on each parameter
	treeviewK <- theWidget("modelsOnKTreeview")
	treemodelK <- gtkTreeViewGetModel(treeviewK)
	treeviewZ <- theWidget("modelsOnZTreeview")
	treemodelZ <- gtkTreeViewGetModel(treeviewZ)
	treeviewG <- theWidget("modelsOnGTreeview")
	treemodelG <- gtkTreeViewGetModel(treeviewG)
	
	NumKEffs <- gtkTreeModelIterNChildren(treemodelK,iter=NULL)
	NumZEffs <- gtkTreeModelIterNChildren(treemodelZ,iter=NULL)
	NumGEffs <- gtkTreeModelIterNChildren(treemodelG,iter=NULL)
	
	checkbox = theWidget("modelsUseZCheckbox")
	useZ = gtkToggleButtonGetActive(checkbox)
	if(!useZ) NumZEffs = 0
	
	if(NumKEffs==0 & NumZEffs==0 & NumGEffs==0){
		gWidgets::gmessage("You must specify at least one effect on one parameter.", title="Model definition error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}
	
	if(!is.null(wommbatAnalysis$Models[[modelName]]))
	{
		gWidgets::gmessage(paste("That name is already in use. Please select a unique name for the model."), title="Model definition error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}
	
	if(modelName=="")
	{
		gWidgets::gmessage(paste("Please select a name for the model."), title="Model definition error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}
	
	
	if ( inherits(invWishartScalar, "try-error") | is.na(invWishartScalar))
	{
		gWidgets::gmessage(paste("Could not coerce the inverse Wishart prior scalar parameter to numeric."), title="Prior error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}

	
	if ( inherits(IGa0, "try-error") | inherits(IGb0, "try-error") | is.na(IGa0) | is.na(IGb0))
	{
		gWidgets::gmessage(paste("Could not coerce the inverse gamma prior parameters to numerics."), title="Prior error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}
	
	if( inherits(MuKMean, "try-error") | 
		 inherits(MuGMean, "try-error") | 
		 inherits(MuZMean, "try-error") |
		 inherits(MuKSD, "try-error") |
		 inherits(MuZSD, "try-error") |
		 inherits(MuGSD, "try-error") |
		 is.na(MuKMean) | is.na(MuGMean) | 
		 is.na(MuZMean) | is.na(MuKSD) |
		 is.na(MuZSD) | is.na(MuGSD))
	{
		gWidgets::gmessage(paste("Could not coerce the grand mean prior parameters to numerics."), title="Prior error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}

	if (MuKSD <=0 | MuZSD <=0 | MuGSD <=0) 
	{
		gWidgets::gmessage(paste("Standard deviations of priors on grand means must be greater than 0."), title="Prior error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}
	
	if (IGa0 <=0 | IGb0 <=0) 
	{
		gWidgets::gmessage(paste("Inverse gamma prior parameters must both be greater than 0."), title="Prior error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(FALSE)
	}
	
	# Check Covariance matrices
	removed = .womRemoveMissingCovarianceItems()
	if(removed)
	{
		gWidgets::gmessage(paste("Removed elements of the covariance groups that were missing from the model."), title="Warning",
					icon="warning",toolkit=gWidgets::guiToolkit("RGtk2"))
	}
	
	covNgroups <- gtkTreeModelIterNChildren(treemodelCov,iter=NULL)
	sizes <- obs <- 1:covNgroups * NA
	iters=list()
	removeIter = rep(FALSE,covNgroups)
	
	if(covNgroups > 0){
		for(i in 1:covNgroups){
			iters[[i]] <- gtkTreeModelIterNthChild(treemodelCov, parent = NULL, i-1)$iter
			covNelements <- gtkTreeModelIterNChildren(treemodelCov,iter=iters[[i]])
			if(covNelements<2){
				
				removeIter[i] = TRUE
			}else{
				sizes[i] <- covNelements
				obs[i] <- as.numeric(gtkTreeModelGetValue(treemodelCov,iters[[i]],.womCovTreeCols("nlevels"))$value)
			}
		}
		for(i in 1:covNgroups)
		{
			if(removeIter[i]) treemodelCov$remove(iters[[i]])
		}
	}
	
	
	
	sizes <- sizes[!removeIter]
	obs <- obs[!removeIter]
	covNgroups <- gtkTreeModelIterNChildren(treemodelCov,iter=NULL)
	if(covNgroups > 0)
	{
		minInvWishartScalar <- min(pmax(sizes - obs,0))
		if( invWishartScalar <= minInvWishartScalar ){
			gWidgets::gmessage(paste("For the selected covariance matrices, the inverse Wishart prior scalar parameter must be greater than ",minInvWishartScalar,".",sep=""), title="Prior error",
					icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
			return(FALSE)
		}
	}
	
	
return(TRUE)
}

.womGetCovInformation <- function()
{
	treeviewCov <- theWidget("modelsCovTreeview")
	treemodelCov <- gtkTreeViewGetModel(treeviewCov)

	covNgroups <- gtkTreeModelIterNChildren(treemodelCov,iter=NULL)
	if(covNgroups==0) return(list(sizes=NULL,obs=NULL))
	
	sizes <- obs <- 1:covNgroups * NA
	
	for(i in 1:covNgroups){
		iter <- gtkTreeModelIterNthChild(treemodelCov, parent = NULL, i-1)$iter
		covNelements <- gtkTreeModelIterNChildren(treemodelCov,iter=iter)
		sizes[i] <- as.integer(covNelements)
		obs[i] <- as.integer(as.numeric(gtkTreeModelGetValue(treemodelCov,iter,.womCovTreeCols("nlevels"))$value))	
	}
	
	return(list(sizes=as.integer(sizes),obs=as.integer(obs)))
}

.womInsertNewModelPackage<-function()
{
	# Model name
	modelName <- theWidget("modelsNameEntry")$getText()
	
	wommbatAnalysis$Models[[modelName]] <- list(priors=list(),model=list())
	
	wommbatAnalysis$Models[[modelName]]$modelName = modelName
	
	wommbatAnalysis$Models[[modelName]]$priors$IGa0 <- as.numeric(theWidget("modelsInvGammaEntry1")$getText())
	wommbatAnalysis$Models[[modelName]]$priors$IGb0 <- as.numeric(theWidget("modelsInvGammaEntry2")$getText())
	
	# Prior mean/variance on grand means
	wommbatAnalysis$Models[[modelName]]$priors$muKMean <- as.numeric(theWidget("modelsMuKMeanEntry")$getText())
	wommbatAnalysis$Models[[modelName]]$priors$muKSD <- as.numeric(theWidget("modelsMuKSDEntry")$getText())
	wommbatAnalysis$Models[[modelName]]$priors$muZMean <- as.numeric(theWidget("modelsMuZMeanEntry")$getText())
	wommbatAnalysis$Models[[modelName]]$priors$muZSD <- as.numeric(theWidget("modelsMuZSDEntry")$getText())
	wommbatAnalysis$Models[[modelName]]$priors$muGMean <- as.numeric(theWidget("modelsMuGMeanEntry")$getText())
	wommbatAnalysis$Models[[modelName]]$priors$muGSD <- as.numeric(theWidget("modelsMuGSDEntry")$getText())
	
	#N effects on each parameter
	treeviewK <- theWidget("modelsOnKTreeview")
	treemodelK <- gtkTreeViewGetModel(treeviewK)
	treeviewZ <- theWidget("modelsOnZTreeview")
	treemodelZ <- gtkTreeViewGetModel(treeviewZ)
	treeviewG <- theWidget("modelsOnGTreeview")
	treemodelG <- gtkTreeViewGetModel(treeviewG)
		
	checkbox = theWidget("modelsUseZCheckbox")
	wommbatAnalysis$Models[[modelName]]$model$useZ = gtkToggleButtonGetActive(checkbox)
	
	# reset effects
	wommbatAnalysis$Models[[modelName]]$model$allEffects = wommbatAnalysis$effs
	
	# .womSetSpecifiedEffects modifies eff to add selected models
	gtkTreeModelForeach(treemodelK,.womSetSpecifiedEffects,user.data=list(parameter=1,modelName=modelName))
	if(wommbatAnalysis$Models[[modelName]]$model$useZ)
	{
		gtkTreeModelForeach(treemodelZ,.womSetSpecifiedEffects,user.data=list(parameter=2,modelName=modelName))
	}
	gtkTreeModelForeach(treemodelG,.womSetSpecifiedEffects,user.data=list(parameter=3,modelName=modelName))
	
	effs <- wommbatAnalysis$Models[[modelName]]$model$allEffects
	
	# setup effects representations
	usedRow = apply(effs[,5:7]==1,1,any)
	SelEffs <- wommbatAnalysis$Models[[modelName]]$model$SelEffs <- effs[usedRow,]
	intMods <- wommbatAnalysis$Models[[modelName]]$model$intMods <- effs[usedRow,1:2]
	Lvls <- as.integer(wommbatAnalysis$Models[[modelName]]$model$SelEffs[,3])
	effects  <- as.integer(t(Lvls*effs[usedRow,5:7]))
	dim(effects)=c(3,length(effects)/3)
	wommbatAnalysis$Models[[modelName]]$model$effects <- effects
	
	# setup reduced data sets
	wommbatAnalysis$Models[[modelName]]$model$newDat2Cat  = .womCreateModColsCat(wommbatAnalysis$newDat,wommbatAnalysis$mods,intMods,wommbatAnalysis$SelCols)
    wommbatAnalysis$Models[[modelName]]$model$newDat2Cont = .womCreateModColsCont(wommbatAnalysis$newDat,wommbatAnalysis$mods,intMods,wommbatAnalysis$SelCols)
	wommbatAnalysis$Models[[modelName]]$model$namedDat2   = .womCreateMeaningfulCols(wommbatAnalysis$newDat,wommbatAnalysis$mods,intMods,wommbatAnalysis$SelCols)
	
    incCont <- wommbatAnalysis$Models[[modelName]]$model$incCont <- .womIncludesContinuous(wommbatAnalysis$SelCols,effs,wommbatAnalysis$mods)[usedRow]	
	
	# Covariance
	treeviewCov <- theWidget("modelsCovTreeview")
	treemodelCov <- gtkTreeViewGetModel(treeviewCov)
	wommbatAnalysis$Models[[modelName]]$priors$invWishartScalar <- as.numeric(theWidget("modelsInvWishartEntry")$getText())
	
	covNgroups <- wommbatAnalysis$Models[[modelName]]$model$covNgroups <- as.integer(gtkTreeModelIterNChildren(treemodelCov,iter=NULL))
	if(covNgroups>0)
	{
		.womParseCovMatSettings(incCont,effects,modelName)
	}			
	
	# Design
	if(!is.null(wommbatAnalysis$Ktype))
	{
		wommbatAnalysis$Models[[modelName]]$Ktype = wommbatAnalysis$Ktype
	}else
	{
		if(theWidget("dataDesignMoreyRadio")$getActive()){
			wommbatAnalysis$Models[[modelName]]$Ktype=as.integer(2)
		}
		if(theWidget("dataDesignPashlerRadio")$getActive()){
			wommbatAnalysis$Models[[modelName]]$Ktype=as.integer(1)
		}
		if(theWidget("dataDesignCowanRadio")$getActive()){
			wommbatAnalysis$Models[[modelName]]$Ktype=as.integer(0)
		}
	}
	# Set new models to be analyzed
	wommbatAnalysis$Models[[modelName]]$toBeAnalyzed <- TRUE
	
	return(modelName)
}

.womParseCovMatSettings = function(incCont,effects,modelName)
{
  	treeviewCov <- theWidget("modelsCovTreeview")
	treemodelCov <- gtkTreeViewGetModel(treeviewCov)
	
	n <- wommbatAnalysis$Models[[modelName]]$model$covNgroups
	parnames=c("K","Z","G")

	covInfo <- .womGetCovInformation()
	obs <- wommbatAnalysis$Models[[modelName]]$model$covObs <- covInfo$obs
	sizes <- wommbatAnalysis$Models[[modelName]]$model$covSizes <- covInfo$sizes
	
	pstart=cumsum(c(0,effects))[1:length(effects)]
	inCovMat = effects*0
  
	dim(pstart)=dim(effects)
	parStart = matrix(nrow=n,ncol=max(sizes))
	effSlope = parStart
  	
	incCont <- .womIncludesContinuous(wommbatAnalysis$SelCols,wommbatAnalysis$Models[[modelName]]$model$SelEffs,wommbatAnalysis$mods)
	effs <- wommbatAnalysis$Models[[modelName]]$model$allEffects
	usedRow = apply(effs[,5:7]==1,1,any)
	#usedRowNumbers <- (1:length(incCont))[usedRow]
	usedRowNumbers <- which(usedRow)
	
	for(i in 1:n)
	{
		iter <- gtkTreeModelIterNthChild(treemodelCov, parent = NULL, i-1)$iter
		covNelements <- gtkTreeModelIterNChildren(treemodelCov,iter=iter)
		for(j in 1:covNelements){
			child.iter <- gtkTreeModelIterNthChild(treemodelCov, parent = iter, j-1)$iter
			parameter = gtkTreeModelGetValue(treemodelCov, child.iter, .womCovTreeCols("parameter"))$value
			row = as.numeric(gtkTreeModelGetValue(treemodelCov, child.iter, .womCovTreeCols("rownumber"))$value)
			eff = which(usedRowNumbers==row)
			parNum <- match(parameter,parnames)
			parStart[i,j] = pstart[parNum,eff]
			effSlope[i,j] = incCont[eff]
			inCovMat[parNum,eff] = 1
		}
	}
	
	parStart = parStart + 3
	
	wommbatAnalysis$Models[[modelName]]$model$inCovMat <- inCovMat
	wommbatAnalysis$Models[[modelName]]$model$parStart <- parStart	
	wommbatAnalysis$Models[[modelName]]$model$effSlope <- effSlope

	

	
}


.womSetSpecifiedEffects <- function(model,path,iter,user.data){
	parameter = user.data$parameter
	modelName = user.data$modelName
	row = gtkTreeModelGetValue(model,iter,.womEffectsTreeCols("rownumber"))$value
	wommbatAnalysis$Models[[modelName]]$model$allEffects[row,parameter+4] = 1
	return(FALSE)
}

.clicked_models_delete_selected <- function(button)
{
	treeview = theWidget("modelsDefinedModelsTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

    if(!is.null(iter)){
		modelName = gtkTreeModelGetValue(model,iter,.womDefinedModelsTreeCols("name"))$value
		model$remove(iter)
		if(!is.null(modelName)) wommbatAnalysis$Models[[modelName]]=NULL
		if(length(wommbatAnalysis$Models)==0)
		{
			.womActiveAnalysisTab(status=FALSE)
			.womActiveDiagnosticsTab(status=FALSE)
			.womActiveResultsTab(status=FALSE)
		}else
		{
			if(all(sapply(wommbatAnalysis$Models,function(x){is.null(x$results)})))
			{
				.womActiveDiagnosticsTab(status=FALSE)
				.womActiveResultsTab(status=FALSE)
			}
		}
		StateEnv$win$present()
	}	
}

.clicked_models_load_selected <- function(treeview, path, view_column)
{
	treeview <- theWidget("modelsDefinedModelsTreeview")
	model = gtkTreeViewGetModel(treeview)	
	
	iter = .womGetModelsModelSelection()$iter
	modelName = gtkTreeModelGetValue(model,iter,.womDefinedModelsTreeCols("name"))$value
	
	.womLoadModel(modelName)
}

.womGetModelsModelSelection<-function()
{
	treeview = theWidget("modelsDefinedModelsTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

	#vis <- model$get(iter, .womDefinedModelsTreeCols("visible"))[[1]]
	return(list(iter=iter,visible=NULL))
}


.womLoadModel <- function(modelName)
{
	myModel <- wommbatAnalysis$Models[[modelName]]
	
	
	
	if(myModel$model$useZ){
		.womTurnOnZ()
	}else{
		.womTurnOffZ()
	}
	
	
	theWidget("modelsNameEntry")$setText(modelName)
	
	# Priors on variances/covariances
	theWidget("modelsInvGammaEntry1")$setText(myModel$priors$IGa0)
	theWidget("modelsInvGammaEntry2")$setText(myModel$priors$IGb0)
	theWidget("modelsInvWishartEntry")$setText(myModel$priors$invWishartScalar)

	# Prior mean/variance on grand means
	theWidget("modelsMuKMeanEntry")$setText(myModel$priors$muKMean)
	theWidget("modelsMuKSDEntry")$setText(myModel$priors$muKSD)
	theWidget("modelsMuZMeanEntry")$setText(myModel$priors$muZMean)
	theWidget("modelsMuZSDEntry")$setText(myModel$priors$muZSD)
	theWidget("modelsMuGMeanEntry")$setText(myModel$priors$muGMean)
	theWidget("modelsMuGSDEntry")$setText(myModel$priors$muGSD)

	# Clear models
	theWidget("modelsOnKTreeview")$getModel()$clear()
	theWidget("modelsOnZTreeview")$getModel()$clear()
	theWidget("modelsOnGTreeview")$getModel()$clear()
	theWidget("modelsCovTreeview")$getModel()$clear()
	
	if(myModel$model$useZ){
		.womTurnOnZ()
	}else{
		.womTurnOffZ()
	}

	.womPlaceEffectsInTreeviews(myModel$model$allEffects,myModel)
	if(myModel$model$covNgroups>0)
	{
		.womPlaceCovariancesInTreeview(myModel)
		
	}
	
}

.womPlaceCovariancesInTreeview <- function(myModel)
{
 	treeview <- theWidget("modelsCovTreeview")
	treemodel <- gtkTreeViewGetModel(treeview)

	covNgroups = myModel$model$covNgroups
	parStart = myModel$model$parStart
	covObs = myModel$model$covObs
	covSizes = myModel$model$covSizes

	# Use the parameter start matrix to figure out
	# which parameter the covariance element points to
	x = cumsum(c(3,myModel$model$effects))
	x = x[-length(x)]
	dim(x)=dim(myModel$model$effects)
	x = x*(myModel$model$effects>0)


	for(i in 1:covNgroups){
		groupName <- paste(covObs[i],"levels")		
		parent.iter <- treemodel$append(NULL)$iter
		treemodel$set(parent.iter, 
			.womCovTreeCols("effectname"), groupName,
			.womCovTreeCols("parameter"), "",
			.womCovTreeCols("nlevels"), as.integer(covObs[i]),
			.womCovTreeCols("rownumber"), "",
			.womCovTreeCols("visible"), TRUE
			)
		for(j in 1:covSizes[i]){
			if(!is.na(parStart[i,j])){
				indices = which(x==parStart[i,j], arr.ind = TRUE)
	
				parameter = c("K","Z","G")[indices[1]]
				effectname = myModel$model$SelEffs[indices[2],4]
				row = match(effectname,myModel$model$allEffects[,4])
				effectnameCov = paste(effectname,"on",parameter,sep=" ")	
			
				child.iter <- treemodel$append(parent.iter)$iter
				treemodel$set(child.iter, 
					.womCovTreeCols("effectname"), effectnameCov,
					.womCovTreeCols("parameter"), parameter,
					.womCovTreeCols("nlevels"), as.integer(covObs[i]),
					.womCovTreeCols("rownumber"), as.integer(row),
					.womCovTreeCols("visible"), TRUE
				)	
			}	
		}
	}
	gtkTreeViewExpandAll(treeview)
}

.womPlaceEffectsInTreeviews <- function(effects,myModel)
{
	nEffects = dim(effects)[1]
	treemodelK <- theWidget("modelsOnKTreeview")$getModel()
	treemodelZ <- theWidget("modelsOnZTreeview")$getModel()
	treemodelG <- theWidget("modelsOnGTreeview")$getModel()

	for(i in 1:nEffects)
	{
		if(effects[i,5]) # K
		{
			iter <- treemodelK$append(NULL)$iter
			treemodelK$set(iter, 
			  .womEffectsTreeCols("effectname"), effects[i,4],
			  .womEffectsTreeCols("nlevels"), as.integer(effects[i,3]),
			  .womEffectsTreeCols("rownumber"), i,
			  .womEffectsTreeCols("visible"), TRUE
			  )
		}
		if(effects[i,6] & myModel$model$useZ) # Z
		{
			iter <- treemodelZ$append(NULL)$iter
			treemodelZ$set(iter, 
			  .womEffectsTreeCols("effectname"), effects[i,4],
			  .womEffectsTreeCols("nlevels"), as.integer(effects[i,3]),
			  .womEffectsTreeCols("rownumber"), i,
			  .womEffectsTreeCols("visible"), TRUE
			  )
		}
		if(effects[i,7]) # G
		{
			iter <- treemodelG$append(NULL)$iter
			treemodelG$set(iter, 
			  .womEffectsTreeCols("effectname"), effects[i,4],
			  .womEffectsTreeCols("nlevels"), as.integer(effects[i,3]),
			  .womEffectsTreeCols("rownumber"), i,
			  .womEffectsTreeCols("visible"), TRUE
			  )
		}
	}
}

.clicked_models_define_new<- function(button)
{
	# Check whether Z will be estimated
	checkbox = theWidget("modelsUseZCheckbox")
	toggled = gtkToggleButtonGetActive(checkbox)
	
	# If not, remove all covariance elements including Z parameters
	if(!toggled){
		.womSetZCovActive(FALSE, remove=TRUE)
	}
	
	# Check the model and prior parameters 
	allGood <- .womCheckModelSettings()
	if(allGood){
		modelName = .womInsertNewModelPackage()
		.womAddModelToTreeview(modelName)
	}
}


.womAddModelToTreeview<-function(modelName)
{
						
	modelsTreeview <- theWidget("modelsDefinedModelsTreeview")
	analysisTreeview <- theWidget("analysisDefinedModelsTreeview")
	diagnosticsTreeview <-theWidget("diagnosticDefinedModelsTreeview")
	resultsTreeview <-theWidget("resultsDefinedModelsTreeview")
	
	model = gtkTreeViewGetModel(modelsTreeview)
	#gtkTreeViewSetModel(analysisTreeview,model)
	#gtkTreeViewSetModel(diagnosticsTreeview,model)
	#gtkTreeViewSetModel(resultsTreeview,model)
		
	.womActiveAnalysisTab(status=TRUE)
	.womInitializeAnalysisPageColumns()

	myModel = wommbatAnalysis$Models[[modelName]]
	
	hasResults <- !is.null(myModel$results)
	parsOnK <- as.integer(.womSaneNum(rowSums(myModel$model$effects)[1]))
	parsOnZ <- as.integer(.womSaneNum(rowSums(myModel$model$effects)[2]))
	parsOnG <- as.integer(.womSaneNum(rowSums(myModel$model$effects)[3]))
		
	iter <- model$append(NULL)$iter
	model$set(iter,
	  		  .womDefinedModelsTreeCols("name"), modelName,
			  .womDefinedModelsTreeCols("hasResults"), hasResults,
			  .womDefinedModelsTreeCols("toBeAnalyzed"), myModel$toBeAnalyzed,
			  .womDefinedModelsTreeCols("parsOnK"), parsOnK,
			  .womDefinedModelsTreeCols("parsOnZ"), parsOnZ,
			  .womDefinedModelsTreeCols("parsOnG"), parsOnG,
			  .womDefinedModelsTreeCols("NCovGroups"), as.integer(myModel$model$covNgroups),
			  .womDefinedModelsTreeCols("muKMean"), as.numeric(.womSaneNum(myModel$priors$muKMean,2)),
			  .womDefinedModelsTreeCols("muKSD"), as.numeric(.womSaneNum(myModel$priors$muKSD,2)),
			  .womDefinedModelsTreeCols("muZMean"), as.numeric(.womSaneNum(myModel$priors$muZMean,2)),
			  .womDefinedModelsTreeCols("muZSD"), as.numeric(.womSaneNum(myModel$priors$muZSD,2)),
			  .womDefinedModelsTreeCols("muGMean"), as.numeric(.womSaneNum(myModel$priors$muGMean,2)),
			  .womDefinedModelsTreeCols("muGSD"), as.numeric(.womSaneNum(myModel$priors$muGSD,2)),
			  .womDefinedModelsTreeCols("invWishartScalar"), as.numeric(.womSaneNum(myModel$priors$invWishartScalar,2)),
			  .womDefinedModelsTreeCols("invGammaA0"), as.numeric(.womSaneNum(myModel$priors$IGa0,2)),
			  .womDefinedModelsTreeCols("invGammaB0"), as.numeric(.womSaneNum(myModel$priors$IGb0,2)),
			  .womDefinedModelsTreeCols("iterations"), NULL,
			  .womDefinedModelsTreeCols("epsUpp"), NULL,
			  .womDefinedModelsTreeCols("epsLow"), NULL,
			  .womDefinedModelsTreeCols("leapfrog"), NULL,
			  .womDefinedModelsTreeCols("useMet"), NULL,
			  .womDefinedModelsTreeCols("metropScale"), NULL,
			  .womDefinedModelsTreeCols("metropThin"), NULL,
			  .womDefinedModelsTreeCols("accRate"), NULL,
			  .womDefinedModelsTreeCols("burnin"), NULL,
			  .womDefinedModelsTreeCols("DIC"), NULL,
			  .womDefinedModelsTreeCols("pD"), NULL,
			  .womDefinedModelsTreeCols("timeAnalyzed"), NULL
		)
    
	
}	



.womIncludesContinuous <- function(SelCols,SelEffs,mods){
  ret=1:dim(SelEffs)[1]
  for(i in 1:dim(SelEffs)[1]){
    ret[i]=any(.womCatOrCont(mods[[2]][mods[[1]][[SelEffs[i,1]]][SelEffs[i,2],]],SelCols)=="CONTIN")
  }
  as.integer(ret)
}


.womCreateDat <- function(data,respColumn,changeColumn,setSizeColumn,interestColumns)
{

	names = colnames(data)
	respColumnNum = match(respColumn,names)
	changeColumnNum = match(changeColumn,names)
	setSizeColumnNum = match(setSizeColumn,names)
	interestColumnsNum = match(interestColumns[,1],names)

	newData = data[,c(respColumnNum,changeColumnNum,setSizeColumnNum,interestColumnsNum)]

	maxLen=prod(unlist(lapply(apply(newData[,-(1:2)],2,unique),length)))
	if(maxLen<length(newData[,1])){

		a=as.vector(table(newData))
		dim(a)=c(4,length(a)/4)
		a=t(a[4:1,])	

		
		colnames(a)=c("Hits","Misses","FalseAlarms","CorrectRej")

		names=dimnames(table(newData))[-(1:2)]
		newData=data.frame(a,expand.grid(names))
		newData=newData[rowSums(newData[,1:4])!=0,]
	}else{
		a=as.integer(cbind(newData[,1] & newData[,2],
			!newData[,1] & newData[,2],
			newData[,1] & !newData[,2],
			!newData[,1] & !newData[,2]))	
		dim(a)=c(length(a)/4,4)
		colnames(a)=c("Hits","Misses","FalseAlarms","CorrectRej")
		newData=data.frame(a,newData[,-(1:2)])
		newData=newData[rowSums(newData[,1:4])!=0,]
	}

	for(i in 1:length(interestColumns[,1])){
		if(interestColumns[i,2]=="CATEG"){
			newData[,i+5] = as.factor(as.character(newData[,i+5]))			
		}else if(interestColumns[i,2]=="CONTIN"){
			newData[,i+5] = as.numeric(as.character(newData[,i+5]))
		}
	}

	newData
}


.womListModels <- function(newDat,interestColumns){
	nFactors=dim(newDat)[2]-5
	names=colnames(newDat)[-(1:5)]
	levels=1:nFactors
	for(i in 1:nFactors){
	  if(interestColumns[i,2]=="CONTIN") levels[i]=1
	  if(interestColumns[i,2]=="CATEG") levels[i]=length(unique(newDat[,5+i])) 
	#levels[i] = factMax(newDat[,5+i])
	}
	modsInt=list()

	for(i in 1:nFactors){
		modsInt[[i]]=gtools::combinations(nFactors,i)
	}

	list(modsInt,names,levels)
}

.womNiceListEffects = function(mods){
	ret=NULL
	for(i in 1:length(mods[[1]])){
		for(j in 1:dim(mods[[1]][[i]])[1]){
			cols=mods[[1]][[i]][j,]
			if(length(cols)==1){
				ret=rbind(ret,c(i,j,mods[[3]][mods[[1]][[i]][j]],mods[[2]][mods[[1]][[i]][j]]))
			}else{
				ret=rbind(ret,c(i,j,prod(mods[[3]][mods[[1]][[i]][j,]]),paste(mods[[2]][mods[[1]][[i]][j,]],collapse=" by ")))
			}

		}
	}
	ret=data.frame(ret)
	ret[,1]=as.numeric(as.character(ret[,1]))
	ret[,2]=as.numeric(as.character(ret[,2]))
	ret[,3]=as.numeric(as.character(ret[,3]))
	ret=data.frame(ret,ret[,1]*0,ret[,1]*0,ret[,1]*0)
	colnames(ret)=c("way","effnum","nlevels","name","K","Z","G")
	ret
}


.womCatOrCont0=function(name,SelCols)
SelCols[match(name,SelCols[,1]),2]

.womCatOrCont=Vectorize(.womCatOrCont0,"name")

.womCreateModColsCat <- function(newDat,allMods,intMods,SelCols)
{

	if(is.null(dim(intMods))){
	  nNewCols = length(intMods)
	}else{
	  nNewCols = dim(intMods)[1]
	}
	#nNewCols = dim(intMods)[1]
	newDat2  = newDat[,1:5]

	names=c("Hits","Misses","FalseAlarms","CorrectRej","setSize")

	for(i in 1:nNewCols){
      		myCols=allMods[[1]][[intMods[i,1]]][intMods[i,2],]
      		nLvls=allMods[[3]][myCols]
      		cc=.womCatOrCont(colnames(newDat)[myCols+5],SelCols)
      		if(length(myCols)>1){
			if(sum(cc=="CATEG")==0){
				newDat2=data.frame(newDat2,newDat[,1]*0+1)
	      			names=c(names,paste("column",i,"nocat",sep=""))
			}else if(sum(cc=="CATEG")==1){
				newDat2=data.frame(newDat2,as.integer(as.factor(newDat[,myCols[which(cc=="CATEG")]+5])))
	      			names=c(names,paste(colnames(newDat)[(myCols+5)[cc=="CATEG"]],collapse=".x."))				      				      		
			}else{
				newDat2=data.frame(newDat2,as.integer(as.factor(apply(cbind(newDat[,(myCols+5)[cc=="CATEG"]]),1,paste,collapse="x"))))
	      			names=c(names,paste(colnames(newDat)[(myCols+5)[cc=="CATEG"]],collapse=".x."))				      			
	      		}
		}else{
			if(any(cc=="CATEG")){
				newDat2=data.frame(newDat2,as.integer(as.factor(newDat[,myCols+5])))
      				names=c(names,colnames(newDat)[myCols+5])
      			}else{
				newDat2=data.frame(newDat2,newDat[,1]*0+1)
      				names=c(names,paste("column",i,"nocat",sep=""))      			
      			}
      		}
	}
colnames(newDat2)=names
return(newDat2)
}

.womCreateModColsCont <- function(newDat,allMods,intMods,SelCols){

	if(is.null(dim(intMods))){
	  nNewCols = length(intMods)
	}else{
	  nNewCols = dim(intMods)[1]
	}
	#nNewCols = dim(intMods)[1]
	newDat2  = newDat[,1:5]

	names=c("Hits","Misses","FalseAlarms","CorrectRej","setSize")

	for(i in 1:nNewCols){
      		myCols=allMods[[1]][[intMods[i,1]]][intMods[i,2],]
      		nLvls=allMods[[3]][myCols]
      		cc=.womCatOrCont(colnames(newDat)[myCols+5],SelCols)
      		if(length(myCols)>1){
			if(any(cc=="CONTIN")){
				newDat2=data.frame(newDat2,apply(cbind(newDat[,(myCols+5)[cc=="CONTIN"]]),1,prod))
	      			names=c(names,paste(colnames(newDat)[(myCols+5)[cc=="CONTIN"]],collapse=".x."))
	      		}else{
				newDat2=data.frame(newDat2,newDat[,1]*0+1)
	      			names=c(names,paste("column",i,"nocont",sep=""))      			
	      		}
		}else{
			if(any(cc=="CONTIN")){
				newDat2=data.frame(newDat2,newDat[,myCols+5])
      				names=c(names,colnames(newDat)[myCols+5])
      			}else{
				newDat2=data.frame(newDat2,newDat[,1]*0+1)
      				names=c(names,paste("column",i,"nocont",sep=""))      			
      			}
      		}
	}
colnames(newDat2)=names
return(newDat2)
}




.womCreateMeaningfulCols <- function(newDat,allMods,intMods,SelCols)#,LevelNames)
{

	if(is.null(dim(intMods))){
	  nNewCols = length(intMods)
	}else{
	  nNewCols = dim(intMods)[1]
	}
	#nNewCols = dim(intMods)[1]
	#nNewCols = dim(newDat)[2]-5
	newDat2  = newDat[,1:5]
	namedCols=data.frame(newDat[,-(1:5)])

	#for(i in 1:nNewCols){
	#	if(SelCols[i,2]=="CATEG")
	#	{
	#		#namedCols[,i]=LevelNames[[i]][namedCols[,i]]
	#	}else{
	#		namedCols[,i]=rep("(CATEGORICAL)",length(newDat2[,1]))		
	#	}
	#}


	names=c("Hits","Misses","FalseAlarms","CorrectRej","setSize")

	for(i in 1:nNewCols){
      		myCols=allMods[[1]][[intMods[i,1]]][intMods[i,2],]
      		nLvls=allMods[[3]][myCols]
      		cc=.womCatOrCont(colnames(newDat)[myCols+5],SelCols)
      		if(length(myCols)>1){
				if(any(cc=="CATEG")){
					if(sum(cc=="CATEG")==1)
					{
						newDat2=data.frame(newDat2,as.character(namedCols[,myCols[cc=="CATEG"]]))
					}else
					{
						newDat2=data.frame(newDat2,apply(cbind(apply(namedCols[,myCols[cc=="CATEG"]],c(1,2),as.character)),1,paste,collapse=".x."))
					}
					
				}else{
					newDat2=data.frame(newDat2,rep("(CATEGORICAL)",length(newDat2[,1])))	      			
				}
				names=c(names,paste(colnames(newDat)[myCols+5],collapse=".x."))      			
			}else{
				if(any(cc=="CATEG")){
					newDat2=data.frame(newDat2,as.character(namedCols[,myCols]))
      			}else{
					newDat2=data.frame(newDat2,rep("(CATEGORICAL)",length(newDat2[,1])))
      			}
      			names=c(names,colnames(newDat)[myCols+5]) 
      		}
	}
colnames(newDat2)=names
return(newDat2)
}

