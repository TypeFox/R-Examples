
.womParEstTreeCols <- function(name)
{
COLUMNS <- c(parameter=0,effect=1,level=2,type=3,postMean=4,postSD=5,postMeanErr=6,bgCol=7)
as.integer(COLUMNS[name])
}



.womActiveResultsTab<-function(status=TRUE)
{

	theWidget("resultsPageBox")$setSensitive(status)

}

.resultsSelectionChanged <- function(selection)
{
	selectedRows <- gtkTreeSelectionGetSelectedRows(selection)
	treeview <- theWidget("resultsParEstTreeview")
	
	model <- selectedRows$model
	rows <-selectedRows$retval
	nRows <-length(rows)
	
	if(nRows>0){
		clipboardText =.womFormatResultsSelectionForClipboard(selection,rows,model)
		clipboard <- treeview$getClipboard(GDK_SELECTION_CLIPBOARD)
  
		clipboard$setText(clipboardText) # copy all the text to the clipboard
	}
	
}

.womFormatResultsSelectionForClipboard <- function(selection, rows, model)
{
	nRows = length(rows)
	clipboardText = ""
	rowText=1:7 * NA
	for(i in 1:nRows)
	{
		path = rows[[i]]
		getIter = gtkTreeModelGetIter(model,path)
		if(!is.null(getIter$retval)){
			for(j in 1:length(rowText))
			{
				rowText[j] = gtkTreeModelGetValue(model, getIter$iter, j-1)$value
			}
			clipboardText = paste(clipboardText,paste(rowText,collapse="\t"),"\n",sep="")
			
		}
	}
	return(clipboardText)
}

.womCreateParEstTree <- function(modelName=NULL)
{
	treeview <- theWidget("resultsParEstTreeview")	
	dataframe = wommbatAnalysis$Models[[modelName]]$results$pointEst
	
	colorColumn = .womCharVectorToAlternatingCol(paste(dataframe[,1],dataframe[,2]))
	dataframe = data.frame(dataframe,colorColumn)
	
	model <- rGtkDataFrame(dataframe)
	gtkTreeViewSetModel(treeview, model)
	
	
	#Create tree view
	if(is.null(gtkTreeViewGetColumn(treeview,0))){
		
		
		# name
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Parameter", renderer, 
  								text = .womParEstTreeCols("parameter"),background=.womParEstTreeCols("bgCol"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

		# effect
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Effect", renderer, 
  								text = .womParEstTreeCols("effect"),background=.womParEstTreeCols("bgCol"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

		# level
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Level", renderer, 
  								text = .womParEstTreeCols("level"),background=.womParEstTreeCols("bgCol"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

		# posterior mean
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Post. mean", renderer, 
  								text = .womParEstTreeCols("postMean"),background=.womParEstTreeCols("bgCol"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

		# MCMC error
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Post. Mean MCMC Error", renderer, 
  								text = .womParEstTreeCols("postMeanErr"),background=.womParEstTreeCols("bgCol"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# postSD
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Post. SD", renderer, 
  								text = .womParEstTreeCols("postSD"),background=.womParEstTreeCols("bgCol"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

		#  Type
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 0.5)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Type", renderer, 
  								text = .womParEstTreeCols("type"),background=.womParEstTreeCols("bgCol"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)

		
	}
}

.womCreateResultsDefinedModelsColumns <- function()
{
	#modelsTreeview <- theWidget("modelsDefinedModelsTreeview")
	resultsTreeview <- theWidget("resultsDefinedModelsTreeview")
		
	#model = gtkTreeViewGetModel(modelsTreeview)	
	#gtkTreeViewSetModel(resultsTreeview, model)	
	
	
	treeview = resultsTreeview
	
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
		column$setClickable(TRUE)


		
		# DIC
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "DIC", renderer, 
  								text = .womDefinedModelsTreeCols("DICtext"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("DIC"))
		
		# pD
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "pD", renderer, 
  								text = .womDefinedModelsTreeCols("pD"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("pD"))
		
		# log likelihood posterior mean
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Log likelihood", renderer, 
  								text = .womDefinedModelsTreeCols("logLikePostMean"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("logLikePostMean"))
		
		# params on K
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "K params", renderer, 
  								text = .womDefinedModelsTreeCols("parsOnK"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("parsOnK"))
		
		# params on Z
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Z params", renderer, 
  								text = .womDefinedModelsTreeCols("parsOnZ"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("parsOnZ"))
		
		# params on G
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "G params", renderer, 
  								text = .womDefinedModelsTreeCols("parsOnG"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		gtkTreeViewColumnSetSortColumnId(column,.womDefinedModelsTreeCols("parsOnG"))
		
		# acc rate
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Acc. Rate", renderer, 
  								text = .womDefinedModelsTreeCols("accRate"),
								foreground = .womDefinedModelsTreeCols("accColor"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# iterations
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Iterations", renderer, 
  								text = .womDefinedModelsTreeCols("iterations"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# effective iterations
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Effective Iterations", renderer, 
  								text = .womDefinedModelsTreeCols("effectiveIterations"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
		# burnin
		renderer <- gtkCellRendererTextNew()
		renderer$set(xalign = 1.0)

		col.offset <- treeview$insertColumnWithAttributes(-1, "Burnin", renderer, 
  								text = .womDefinedModelsTreeCols("burnin"))
								
		column <- treeview$getColumn(col.offset - 1)
		column$setClickable(TRUE)
		
	}
}


.womGetResultsModelSelection<-function()
{
	treeview = theWidget("resultsDefinedModelsTreeview")
	model <- gtkTreeViewGetModel(treeview)
	selection <- gtkTreeViewGetSelection(treeview)
	iter <- gtkTreeSelectionGetSelected(selection)$iter	

	return(list(iter=iter))
}


.selected_results_parameter_row<-function(treeview, path, view_column)
{
	selection <- treeview$getSelection()
	selectedRows <- gtkTreeSelectionGetSelectedRows(selection)
	
	model <- selectedRows$model
	rows <-selectedRows$retval
	nRows <-length(rows)
	
	if(nRows!=1) return(0)
	
	row = as.integer(gtkTreePathToString(rows[[1]]))+1
	getIter = gtkTreeModelGetIter(model,rows[[1]])

	if(is.null(getIter$retval)) return(0)
	
	
	modelName = gtkTreeModelGetValue(theWidget("resultsDefinedModelsTreeview")$getModel(),
							.womGetResultsModelSelection()$iter,
							.womDefinedModelsTreeCols("name"))$value
	
	.womCreateMultiComp(modelName, row)
}

.selected_results_model_row <- function(treeview, path, view_column)
{
	treeview <- theWidget("resultsDefinedModelsTreeview")
	model = gtkTreeViewGetModel(treeview)	
	iter = .womGetResultsModelSelection()$iter

	modelName = gtkTreeModelGetValue(model,iter,.womDefinedModelsTreeCols("name"))$value
	
	if(is.null(wommbatAnalysis$Models[[modelName]]$results))
	{
		theWidget("resultsParEstTreeview")$setModel(NULL)
		return()
	}
	.womCreateParEstTree(modelName)
	
}

.womCreateMultiComp<-function(modelName,row)
{
	if(row>3){
		pointEst <- wommbatAnalysis$Models[[modelName]]$results$pointEst
		group <- pointEst[row,2]
		par <- pointEst[row,1]
		myRows <- pointEst[,2]==group & pointEst[,1]==par
		myRows[row]=FALSE # eliminate current row
		myNames = pointEst[myRows,3]
		nRows = sum(myRows)
		
		burnin = wommbatAnalysis$Models[[modelName]]$settings$burninIters
		effective = wommbatAnalysis$Models[[modelName]]$settings$effectiveIters
		goodIterations = (burnin+1):effective
		nGoodIterations = effective - burnin	
		
		chains <- array(wommbatAnalysis$Models[[modelName]]$results$effectChains[goodIterations,myRows],dim=c(nGoodIterations,nRows))
		
		meanDiff = sdDiff = probDiffgt0 = vector(mode = "numeric", length = nRows)
				
		selChain = wommbatAnalysis$Models[[modelName]]$results$effectChains[goodIterations,row]
		for(i in 1:nRows)
		{
			myDiff = selChain - chains[,i]
			meanDiff[i] = mean(myDiff)
			sdDiff[i] = sd(myDiff)
			probDiffgt0[i] = mean(myDiff>0)
		}
		
		# make data.frame
		diffDataFrame = data.frame(myNames,meanDiff,sdDiff,probDiffgt0,probDiffgt0/(1-probDiffgt0))
		colnames(diffDataFrame) = c("Level name","Posterior mean difference","Posterior SD difference", "Posterior prob. difference > 0","Posterior odds difference > 0")
		
		
		window <- gtkWindowNew("toplevel", show=F)
		gtkWindowSetDefaultSize(window, 700, 400)
		window$setTitle("Multiple comparisons")
		gSignalConnect(window, "destroy", .present_main_window_after_destroy)
		
		
		vbox <- gtkVBoxNew(FALSE, 0)
		vbox$setBorderWidth(8)
		label <- gtkLabelNew(paste("Posterior differences from parameter\n",group," on ",par,"\nLevel: ",pointEst[row,3],sep=""))
		gtkLabelSetJustify(label,"center")
		
		vbox$packStart(label, FALSE, FALSE, 0)

		scrolledWin = gtkScrolledWindowNew()
		treeview <- gtkTreeViewNew()
		model <- rGtkDataFrameNew(frame = diffDataFrame)
		gtkTreeViewSetModel(treeview,model)
		
		for(i in 1:dim(diffDataFrame)[2])
		{
			renderer <- gtkCellRendererTextNew()
			renderer$set(xalign = 1.0)
  
			col.offset <- treeview$insertColumnWithAttributes(-1, colnames(diffDataFrame)[i], renderer,
							    text = i-1)
			column <- treeview$getColumn(col.offset - 1)
			column$setClickable(TRUE)
		}
		vbox$packStart(scrolledWin, TRUE, TRUE, 0)
		gtkContainerAdd(scrolledWin,treeview)
		
		# do selection stuff
		selection = treeview$getSelection()
		gtkTreeSelectionSetMode(selection, "GTK_SELECTION_MULTIPLE")
		gtkTreeViewSetRubberBanding(treeview, TRUE)
		gSignalConnect(selection, "changed", .multi_comp_selection_changed)
		
		window$add(vbox)
		window$showAll()
	}else
	{
		return(0)
	}
}


.multi_comp_selection_changed <- function(selection)
{
	selectedRows <- gtkTreeSelectionGetSelectedRows(selection)
	treeview <- gtkTreeSelectionGetTreeView(selection)
	
	model <- selectedRows$model
	rows <-selectedRows$retval
	nRows <-length(rows)
	
	if(nRows>0){
		clipboardText =.womFormatMultiCompSelectionForClipboard(selection)
		clipboard <- treeview$getClipboard(GDK_SELECTION_CLIPBOARD)
  
		clipboard$setText(clipboardText) # copy all the text to the clipboard
	}
		
}


.womFormatMultiCompSelectionForClipboard <- function(selection)
{
	selectedRows <- gtkTreeSelectionGetSelectedRows(selection)
	treeview <- gtkTreeSelectionGetTreeView(selection)
	
	model <- selectedRows$model
	rows <-selectedRows$retval
	nRows <-length(rows)

	clipboardText = ""
	rowText=1:5 * NA
	for(i in 1:nRows)
	{
		path = rows[[i]]
		getIter = gtkTreeModelGetIter(model,path)
		if(!is.null(getIter$retval)){
			for(j in 1:length(rowText))
			{
				rowText[j] = gtkTreeModelGetValue(model, getIter$iter, j-1)$value
			}
			clipboardText = paste(clipboardText,paste(rowText,collapse="\t"),"\n",sep="")
			
		}
	}
	return(clipboardText)

}


