# plot functions to visualize the outputs of AssessMonophyly
# written by Orlando Schwery 2015
PlotMonophyly <-
function(solution, tree, taxlevels=1, plot.type='monophyly', monocoll=FALSE, ladderize=TRUE, PDF=FALSE, PDF_filename='Monophylyplot.pdf', PDF_width='auto', PDF_height='auto', mono.colour='PRGn', tax.colour='rainbow', intrud.colour='rainbow', edge.width=3, cex=0.2, adj.names='auto', adj.tips=0.5, label.offset='auto', type='phylogram', ...) {
# warnings and data preparation
    if (taxlevels == 'ALL') {  # test if taxlevels argumetn has correct format and display error if not
		stop(" 'ALL' is not an option for plotting!")
    }
    if (class(taxlevels) == 'numeric' & taxlevels > length(solution)) {  # check if specified taxlevel is among the ones in the solution object and display error if not
		stop('Requested taxonomic level not available (less levels specified as analysis input)!')
    }
    if (plot.type != 'monophyly' & plot.type != 'monoVStax' & plot.type != 'intruders' & plot.type != 'taxonomy') {  # check if plot type is one of the implemented ones
		stop('Invalid plot.type!')
    }
	if (ladderize == TRUE) {  # ladderizes the tree before starting, if specified
        tree <- ladderize(tree)
    }
    tip.states <- solution[[taxlevels]]$TipStates  # extract tip states data frame of requested taxlevel from solution
    tip.states <- tip.states[match(tree$tip.label, tip.states$Tip), ]  # match/sort order of tip states in data frame with order of tip labels in tree
	row.names(tip.states) <- 1:nrow(tip.states)  # renumber rows of sorted table
# collapsing monophyletic groups to one tip per group
    if (monocoll == TRUE) {  # check if monophyletic groups should be collapsed
		tip.states$Tip <- as.character(tip.states$Tip)  # turn tip labels into characters
		tip.states$Taxon <- as.character(tip.states$Taxon)  # turn taxon names into characers
		tip.states$Status <- as.character(tip.states$Status)  # turn status into characters
		alltaxa <- as.vector(unique(tip.states[, "Taxon"]))  # create vector of all unique taxon names in table
		keeptips <- c()  # create empty vector for tips to keep
		colltipse <- c()  # create empty vector for tips to collapse/discard
		for (icoll in 1:length(alltaxa)) {  # loop through unique taxa
		    if (solution[[taxlevels]]$result[alltaxa[icoll], "Monophyly"] == "Yes") {  # if the current taxon is monophyletic (according to results table in solution)
				matchtips <- which(tip.states[, "Taxon"] == alltaxa[icoll])  # get vector with all tip.states rows that belong (match) to current taxon
				colltips <- c(matchtips[2:length(matchtips)])  # select all but first of these tips...
				colltipse <- c(colltipse, colltips)  # ...and add them to discard vector
				keeptip <- c(matchtips[1])  # select first of these tips...
				keeptips <- c(keeptips, keeptip)  # ...and add it to vector with tips to keep
		    }
		}
		tip.states.temp <- tip.states  # create temporary copy of tip.states table
		for (icollstates in 1:length(colltipse)) {  # loop through vector of tips to discard...
		    tip.states.temp <- tip.states.temp[!(tip.states.temp$Tip == tree$tip.label[colltipse[icollstates]]), ]  # and drop them from temporary table
		}
		row.names(tip.states.temp) <- 1:nrow(tip.states.temp)  # renumber rownames of new table
		for (inameadj in 1:length(tip.states.temp$Tip)) {  #loop through tip names of table...
			if (solution[[taxlevels]]$result[tip.states.temp$Taxon[inameadj], "Monophyly"] == "Yes") {  # ...and check if they belong to a monophyletic group (according to results table in solution)...
				tip.states.temp$Tip[inameadj] <- tip.states.temp$Taxon[inameadj]  # ...and then replace tip name with taxon name
			}
		}
		newlabels <- tree$tip.label  # create temporary vector of tip labels
		for (itips in 1:length(keeptips)) {  # loop through vector of tips to keep
		    newlabels[keeptips[itips]] <- as.character(tip.states[keeptips[itips], "Taxon"])  # replace tip names of tips to keep with the name of their taxon
		}
		tree$tip.label <- newlabels  # replace tip labels of tree with modified ones
		tree <- drop.tip(tree, colltipse)  # drop discard tips from tree
		tip.states <- tip.states.temp  # replace tip states table with modified one
    }
##########################
# Reconstruct Monophyly
    if (plot.type == 'monophyly' | plot.type == 'monoVStax' | plot.type == 'intruders') {  # for the plot variants that require monophyly status reconstructed
        mono.tree <- tree  # assigns tree specifically to reconstruction to avoid conflicts
# assign numbers to tip status
        tipdata <- as.character(tip.states[, "Status"])  # extracting monophyly status of tips from table
        tipdata[tipdata == "Monophyletic"] <- 2  # number-coding monophyly status: monophyletic
        tipdata[tipdata == "Non-Monophyletic"] <- 3  # number-coding monophyly status: non-monophyletic
        tipdata[tipdata == "Intruder"] <- 4  # number-coding monophyly status: intruder
        tipdata[tipdata == "Outlier"] <- 4  # number-coding monophyly status: outlier
		tipdata[tipdata == "unknown"] <- 1  # number-coding monophyly status: unknown
		tipdata <- as.numeric(tipdata)  # format number coded tip states as numeric
		names(tipdata) <- tip.states[, "Tip"]  # assign tip names to number coded tip states
# run reco model
        monophyly.reco <- fastAnc(mono.tree, tipdata, vars=FALSE, CI=FALSE)  # run reconstruction of numeric monophyly states
        edgestates <- c()  # empty vector for reconstructed states at edges
        for (i in 1:length(mono.tree$edge[, 2])) {  # loop through edges of tree
            if (mono.tree$edge[i, 2] > length(mono.tree$tip.label)) {  # if internal edge
                edgestates.i <- as.vector(monophyly.reco[mono.tree$edge[i, 2] - (length(mono.tree$tip.label))])  # use reconstructed edge state
            } else {  # if terminal branch
                edgestates.i <- as.vector(tipdata[mono.tree$edge[i, 2]])  # use tip state
            }
            edgestates <- c(edgestates, edgestates.i)  # add to vector
        }
        mono.tree$edge <- cbind(mono.tree$edge, round(as.numeric(edgestates), digits=0))  # round values to full digits (i.e. initial states) and add to edges of tree
    }
##########################
# reconstruct intruders in monophyly
    if (plot.type == 'intruders') {  # for the intruder plot
        int.tree <- tree  # assigns tree specifically to reconstruction to avoid conflicts
# assign numbers to tip status
        tipdataI <- as.character(tip.states[, "Status"])  # extracting monophyly status of tips from table
        tipdataI[tipdataI == "Monophyletic"] <- 2  # number-coding monophyly status: monophyletic
        tipdataI[tipdataI == "Non-Monophyletic"] <- 3  # number-coding monophyly status: non-monophyletic
        tipdataI[tipdataI == "Intruder"] <- 4  # number-coding monophyly status: intruder
        tipdataI[tipdataI == "Outlier"] <- 4  # number-coding monophyly status: outlier
		tipdataI[tipdataI == "unknown"] <- 1  # number-coding monophyly status: unknown
        tipdataI <- as.numeric(tipdataI)  # format number coded tip states as numeric
# assign numbers to tip genus
        tipdataII <- as.character(tip.states[, "Taxon"])  # vector with genus names
        taxai <- c()  # create empty vector for intruder/outlier state taxa
        for (i in 1:length(tipdataI)){  # loop through numeric tip states
            if (tipdataI[i] == 4){  # pick taxa of tips classified as intruders or outliers...
                taxai <- c(taxai, tipdataII[i])  # ...and add them to vector
            }
        }
        if (length(taxai) == 0) {  # stop and display warning if no intruders/outliers are present
            stop('No intruders present to be plotted! Get on with it!')
        }
		taxaI <- as.vector(unique(taxai))  # create vector of intruder taxa (without doubles)
        taxaII <- c()  # create empty vector to be filled with taxon number code
        for (i in 1:length(taxaI)) {  # loop through unique taxon names of intruders/outliers
            taxaII[i] <- i - 1  # create vector of numbers assigned to intruding taxon (starting at 0)
        }
        taxaIII <- matrix(c(taxaI, taxaII), nrow=length(taxaI)) #create matrix consisting of intruder taxon names and their number
        tipdataIII <- c()  # create empty vector for tip states to be used in reconstruction
        for (i in 1:length(tipdataI)) {  # loop through numeric tip states...
            if (tipdataI[i] == 4) {  # ...and for all which are intruders/outliers...
                tipdataIII <- c(tipdataIII,(as.numeric(tipdataI[i]) + as.numeric(taxaIII[(which(taxaIII[, 1] == tipdataII[i])), 2])))  # ...add previously assigned taxon number (leading intruders to be coded from 4 onwards)
            } else {  # if not intruder/outlier...
                tipdataIII <- c(tipdataIII, as.numeric(tipdataI[i]))  # ...keep number 1, 2 or 3 respectively
            }
        }
		tipdataIII <- as.numeric(tipdataIII)  # make final intruder/outlier taxon tip state vector numeric
		names(tipdataIII) <- tip.states[, "Tip"]  # assign tip names to states vector
# run reco model
        int.reco <- fastAnc(int.tree, tipdataIII, vars=FALSE, CI=FALSE)  # run reconstruction of numeric intruder/outlier states
        edgestatesI <- c()  # empty vector for reconstructed states at edges
        for (i in 1:length(int.tree$edge[, 2])) {  # loop through tree edges...
            if (int.tree$edge[i, 2] > length(int.tree$tip.label)) {  # ...if internal edge...
                edgestatesI.i <- as.vector(int.reco[int.tree$edge[i, 2] - (length(int.tree$tip.label))])  # ...use reconstructed edge state...
            } else {  # ...if terminal branch...
                edgestatesI.i <- as.vector(tipdataIII[int.tree$edge[i, 2]])   # ...use tip state
            }
            edgestatesI <- c(edgestatesI, edgestatesI.i)  # ... and add to vector
        }
        edgestates.monoround <- c(round(as.numeric(edgestates), digits=0)) # round values to full digits (i.e. initial states)
        edgestatesII <- edgestates.monoround  # assign rounded edge states new name
        for (i in 1:length(edgestates)) {  # loop through monophyly edgestates...
        	if (edgestates.monoround[i] == 4) {  # ...if they are intruders/outliers...
        		edgestatesII[i] <- edgestatesI[i]  # ...replace their number with the one from the intruder reconstruction
        	}
        }
		int.tree$edge <- cbind(int.tree$edge, as.numeric(edgestatesII))  # add to edges of tree
    }
##########################
# reconstruct taxonomy
    if (plot.type == 'monoVStax' | plot.type == 'taxonomy') {  # for the plot variants that require taxonomy status reconstructed
        tax.tree <- tree  # assigns tree specifically to reconstruction to avoid conflicts
        tipdataT <- as.character(tip.states[, "Taxon"])  # extract taxon of tip from table
        taxaT <- as.vector(unique(tip.states[, "Taxon"])) # vector with taxon names (no doubles)
        for (i in 1:length(taxaT)) {  # loop through unique taxa
            tipdataT[tipdataT == taxaT[i]] <- i  # translate taxon associated to tip into taxon specific number
        }
        tipdataT <- as.numeric(tipdataT)  # format tip taxon coding to be numeric
# run reco model
        taxa.reco <- fastAnc(tax.tree, tipdataT, vars=FALSE, CI=FALSE)  # run reconstruction of numeric intruder/outlier states
        edgestatesT <- c()  # empty vector for reconstructed states at edges
        for (i in 1:length(tax.tree$edge[, 2])) {  # loop through tree edges...
            if (tax.tree$edge[i, 2] > length(tax.tree$tip.label)) {  # ...if internal edge...
                edgestatesT.i <- as.vector(taxa.reco[tax.tree$edge[i, 2] - (length(tax.tree$tip.label))])  # ...use reconstructed edge state...
            } else {  # ...if terminal branch...
                edgestatesT.i <- as.vector(tipdataT[tax.tree$edge[i, 2]])    # ...use tip state
            }
            edgestatesT <- c(edgestatesT, edgestatesT.i)  # add picked state to vector
        }
        tax.tree$edge <- cbind(tax.tree$edge, as.numeric(edgestatesT))  # add edge state vector to edges of tree
    }
##########################
# plotting itself
    if (PDF == TRUE) {  # if output should be printed to PDF
		if (PDF_width == 'auto') {  # if PDF width is atuomatically selected
			if (type == "fan" | type == "radial") {  # for circular tree shapes...
				pdf_width <- 2.5 * sqrt(length(tree$tip.label)) / pi  # ...create PDF with width adjusted to tree size (square shaped for round trees)
			} else {  # for linear tree shapes...
				pdf_width <- (9 - (3 - (3 * (2 ^ -(length(tree$tip.label) / 100)))))  # ...create PDF with width adjusted to tree size (rectangular for straight trees)
			}
		} else {  # if PDF width is customized...
		    pdf_width <- PDF_width  # ...paste specified value
		}
		if (PDF_height == 'auto') {  # if PDF heigth is atuomatically selected
		    if (type == "fan" | type == "radial") {  # for circular tree shapes...
				pdf_height <- 2.5 * sqrt(length(tree$tip.label)) / pi  # ...create PDF with lenght adjusted to tree size (square shaped for round trees)
		    } else {  # for linear tree shapes...
				pdf_height <- (length(tree$tip.label) / 10)  # ...create PDF with lenght adjusted to tree size (rectangular for straight trees)
			}
		} else {  # if PDF width is customized...
		    pdf_height <- PDF_height  # ...paste specified value
		}
		if (pdf_height > 200) {  # if manual or automatic height exceeds 200in...
			pdf_height <- 200  # ... set it back to 200in and display warning
			print('Warning: pdf_height too large, capped to 200in. If output not satisfying, consider different plotting type or tree slicing.')
		}
		if (pdf_width > 200) {  # if manual or automatic width exceeds 200in...
		    pdf_width <- 200  # ... set it back to 200in and display warning
		    print('Warning: pdf_width too large, capped to 200in. If output not satisfying, consider different plotting type or tree slicing.')
		}
    }
    if (adj.names == 'auto') {  # if horizontal adjustment of tip.labels is automatic...
		if (plot.type == 'monophyly' | plot.type == 'intruders' | plot.type == 'taxonomy') {  # ...for single tree plot types...
			adj.name <- 0  # set it to left adjusted
		} else if (plot.type == 'monoVStax') {  # ...for mirrored tree plot...
			adj.name <- 0.5  # set it to centered
		}
    } else {  # if horizontal adjustment of tip.labels is customised...
		adj.name <- adj.names  # ...paste specified value
    }
    if (label.offset == 'auto') {  # if label offset is automatic...
		if (plot.type == 'monophyly' | plot.type == 'intruders' | plot.type == 'taxonomy') {  # for single tree plot types...
			labelOffset <- 50000 / (length(tree$tip.label) ^ 2)  # determine it based on the number of tips
		} else if (plot.type == 'monoVStax') {  # for mirrored tree plot...
			 labelOffset <- 200000 / (length(tree$tip.label) ^ 2)  # determine it based on the number of tips
		}
	} else {  # if label offset is customized...
		labelOffset <- label.offset  # ...paste specified value
    }
# monophyly plot
    if (plot.type == 'monophyly') {  # for monophyly plot
        if (mono.colour == 'PRGn') {  # use colours from colourblind friendly palettes of RColorBrewer
			co <- c('gray', '#5aae61', '#c2a5cf', '#762a83') #green and purple, which is the standard colour scheme
		} else if (mono.colour == 'RdBu') {
			co <- c('gray', '#4393c3', '#f4a582', '#b2182b') #red and blue
		} else if (mono.colour == 'PuOr') {
		    co <- c('gray', '#542788', '#fdb863', '#b35806') #orange and purple
		} else if (mono.colour == 'PiYG') {
		    co <- c('gray', '#4d9221','#f1b6da', '#8e0152') #green and pink
		} else if (mono.colour == 'BrBG') {
		    co <- c('gray', '#35978f','#dfc27d', '#543005') #petrol and brown
		} else {  # if none of the predefined palettes chosen...
            co <- mono.colour  # ...use custom colours
        }
        if (PDF == TRUE) {  # if output should be printed to PDF
	    pdf(PDF_filename, width=pdf_width, height=pdf_height)  # create PDF frame using width and height (see above)
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], show.tip.label=TRUE, cex=cex, adj=adj.name, label.offset=labelOffset, edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction and all other parameters as specified (see above)
            tiplabels(pch=22, bg=co[tipdata], cex=1, adj=adj.tips)  # add square tip labels with tip state colour
            dev.off()  # close PDF device
        } else {  # if not plotted to PDF but R plot window
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], show.tip.label = TRUE, cex=cex, adj = adj.name, label.offset=labelOffset, edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction and all other parameters as specified (see above)
            tiplabels(pch=22, bg=co[tipdata], cex=1, adj=adj.tips)  # add square tip labels with tip state colour
        }
    }
# taxonomy plot
    if (plot.type == 'taxonomy') {  # for taxonomy plot
        if (tax.colour == 'rainbow') {  # if default rainbow colours are chosen...
            coTax <- rainbow(length(taxaT))  # ...define rainbow colours based on number of taxa
        } else {  #if customized colours are specified...
            coTax <- tax.colour  # ... assign custom taxonomy colours
        }
        names(coTax) <- 1:length(taxaT)  # assign taxon numbers as names
        if (PDF == TRUE) {  # if output should be printed to PDF
			pdf(PDF_filename, width=pdf_width, height=pdf_height)  # create PDF frame using width and height (see above)
            plot(tax.tree, edge.col=coTax[as.numeric(tax.tree$edge[, 3])], show.tip.label=TRUE, cex=cex, adj=adj.name, label.offset=labelOffset, edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction and all other parameters as specified (see above)
            tiplabels(pch=22, bg=coTax[tipdataT], cex=1, adj=adj.tips)  # add square tip labels with tip state colour
            dev.off()  # close PDF device
        } else {  # if not plotted to PDF but R plot window
            plot(tax.tree, edge.col=coTax[as.numeric(tax.tree$edge[, 3])], show.tip.label=TRUE, cex=cex, adj=adj.name, label.offset=labelOffset, edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction and all other parameters as specified (see above)
            tiplabels(pch=22, bg=coTax[tipdataT], cex=1, adj=adj.tips)  # add square tip labels with tip state colour
        }
    }
# intruders plot
    if (plot.type == 'intruders') {  # for intruder plot
	    if (intrud.colour == 'rainbow') {  # if default rainbow colours are chosen...
	        coInt <- c('gray82', 'gray50', 'black', rainbow(length(taxaI)))  # ...use default colours (gray for monophyletic, black for invaded and rainbow colours for invading taxa)
	    } else {  #if customized colours are specified...
	        coInt <- intrud.colour  # ...assign custom intruder colours
	    }
        names(coInt) <- 1:(length(taxaI) + 2)  #assign taxon numbers as names (adjusted for the non-intruder ones)
        if (PDF == TRUE) {  # if output should be printed to PDF
			pdf(PDF_filename, width=pdf_width, height=pdf_height)  # create PDF frame using width and height (see above)
			plot(int.tree, edge.col=coInt[as.numeric(int.tree$edge[, 3])], show.tip.label=TRUE, cex=cex, adj=adj.name, label.offset=labelOffset, edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction and all other parameters as specified (see above)
            tiplabels(pch=22, bg=coInt[tipdataIII], cex=1, adj=adj.tips)  # add square tip labels with tip state colour
            dev.off()  # close PDF device
        } else {  # if not plotted to PDF but R plot window
            plot(int.tree, edge.col=coInt[as.numeric(int.tree$edge[, 3])], show.tip.label=TRUE, cex=cex, adj=adj.name, label.offset=labelOffset, edge.width=edge.width, type=type, ...)  # plot tree with edge colours according to reconstruction and all other parameters as specified (see above)
            tiplabels(pch=22, bg=coInt[tipdataIII], cex=1, adj=adj.tips)  # add square tip labels with tip state colour
        }
    }
# monophyly vs taxonomy mirror tree plot
    if (plot.type == 'monoVStax') {  # for monophyly vs. taxonomy mirror plot
        if (type == 'fan' | type == 'unrooted' | type == "radial") {  # if tree type is unsuitable to be plotted as mirror plot, stop and display error
			stop("Phylogeny types 'fan', 'radial' and 'unrooted' aren't suited for plot.type 'monoVStax', use types 'phylogram' or 'cladogram' instead!")
        }
        if (mono.colour == 'PRGn') {  # use colours from colourblind friendly palettes of RColorBrewer
			co <- c('gray', '#5aae61', '#c2a5cf', '#762a83') #green and purple, which is the standard colour scheme
		} else if (mono.colour == 'RdBu') {
		    co <- c('gray', '#4393c3', '#f4a582', '#b2182b') #red and blue
		} else if (mono.colour == 'PuOr') {
		    co <- c('gray', '#542788','#fdb863', '#b35806') #orange and purple
		} else if (mono.colour == 'PiYG') {
		    co <- c('gray', '#4d9221', '#f1b6da', '#8e0152') #green and pink
		} else if (mono.colour == 'BrBG') {
		    co <- c('gray', '#35978f', '#dfc27d', '#543005') #petrol and brown
		} else {  # if none of the predefined palettes chosen...
            co <- mono.colour  # ...assign custom monophyly colours
        }
        if (tax.colour == 'rainbow') {  # if default rainbow colours are chosen...
            coTax <- rainbow(length(taxaT))  # ...use default rainbow colours based on number of taxa
        } else {  # if none of the predefined palettes chosen...
            coTax <- tax.colour  # ...assign custom taxonomy colours
        }
        names(coTax) <- 1:length(taxaT)  # assign taxon numbers as names
        if (PDF == TRUE) {  # if output should be printed to PDF
            pdf(PDF_filename, width=pdf_width, height=pdf_height)  # create PDF frame using width and height (see above)
            par(oma=c(1, 1, 1, 1), mar=c(0, 0, 0, 0))  # set up plotting margins
            layout(matrix(c(1, 2), 1, 2, byrow=TRUE), widths=c(1.5, 1))  # set up for plotting two trees next to each other
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], show.tip.label=TRUE, cex=cex, adj=adj.name, label.offset=labelOffset, edge.width=edge.width, no.margin=TRUE, type=type, ...)  # plot monophyly tree with edge colours according to reconstruction and all other parameters as specified (see above)
            tiplabels(pch=22, bg=co[tipdata], cex=1, adj=adj.tips)  # add square tip labels with tip state colour
            plot(tax.tree, edge.col=coTax[as.numeric(tax.tree$edge[, 3])], show.tip.label=FALSE, cex=cex, label.offset=0, edge.width=edge.width, direction="leftwards", no.margin=TRUE, type=type, ...)  # plot mirrored taxonomy tree with edge colours according to reconstruction and all other parameters as specified (see above)
            tiplabels(pch=22, bg=coTax[tipdataT], cex=1, adj=adj.tips)  # add square tip labels with tip state colour
            dev.off()  # close PDF device
        } else {  # if not plotted to PDF but R plot window
            par(oma=c(1, 1, 1, 1), mar=c(0, 0, 0, 0))  # set up plotting margins
            layout(matrix(c(1, 2), 1, 2, byrow=TRUE), widths=c(1.5, 1))  # set up for plotting two trees next to each other
            plot(mono.tree, edge.col=co[as.numeric(mono.tree$edge[, 3])], show.tip.label=TRUE, cex=cex, adj=adj.name, label.offset=labelOffset, edge.width=edge.width, no.margin=TRUE, type=type, ...)  # plot monophyly tree with edge colours according to reconstruction and all other parameters as specified (see above)
            tiplabels(pch=22, bg=co[tipdata], cex=1, adj=adj.tips)  # add square tip labels with tip state colour
            plot(tax.tree, edge.col=coTax[as.numeric(tax.tree$edge[, 3])], show.tip.label=FALSE, cex=cex, label.offset=0, edge.width=edge.width, direction="leftwards", no.margin=TRUE, type=type, ...)  # plot mirrored taxonomy tree with edge colours according to reconstruction and all other parameters as specified (see above)
            tiplabels(pch=22, bg=coTax[tipdataT], cex=1, adj=adj.tips)  # add square tip labels with tip state colour
        }
    }
}
