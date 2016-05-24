# The NanoStringNorm package is copyright (c) 2013 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

Plot.NanoStringNorm.gvis <- function(x, plot.type = c("gene.norm", "sample"), save.plot = FALSE, path.to.mongoose = "web", output.directory = "NanoStringNorm_gvis_plots") {

	if (!suppressPackageStartupMessages(requireNamespace("googleVis"))) {
		stop ("Plot.NanoStringNorm.gvis:  googleVis is not available.");
		}

	for (plot.item in plot.type) {

		if (plot.item == "gene.norm") {
			NSN.output.name <- "gene.summary.stats.norm";
			idvar <- "Gene";
			}
		else if (plot.item == "gene.raw") {
			NSN.output.name <- "gene.summary.stats.raw";
			idvar <- "Gene";
			}
		else if (plot.item == "sample") {
			NSN.output.name <- "sample.summary.stats.norm";
			idvar <- "Sample";
			}
		else {
			stop(paste("Plot.NanoStringNorm.gvis:  Unrecognized plot.type", plot.item))
			};

		data.to.plot <- x[[NSN.output.name]];
		
		# add the annotation and dummy time variable for classification in plotting
		if (grepl("gene", plot.item) ) { 
			data.to.plot <- data.frame(Gene = rownames(data.to.plot), time = 1, Code.Class = x$raw.data$Code.Class, data.to.plot);
			}
		else if (grepl("sample", plot.item)) {
			data.to.plot <- data.to.plot[rownames(x$traits),];
			data.to.plot <- data.frame(Sample = rownames(data.to.plot), time = 1, data.to.plot, x$traits);
			#colnames(data.to.plot)[1] <-  'Sample';
			# add a prefix to the sample names because they sometimes cause errors
			data.to.plot$Sample <- paste(1:nrow(data.to.plot), data.to.plot$Sample, sep = "-" );
			}

		data.to.plot$time <- 1;
		
		# take the -log10P for scaling.  This should only be done on gene.norm
		pval.columns <- colnames(data.to.plot)[grepl("P_",colnames(data.to.plot))];
		for (pval in pval.columns) {
			data.to.plot[,pval] <- round(-log10(data.to.plot[,pval]),2);
			}

		# the intial plotting parameters
		initial.plotting.parameters <- '{"iconKeySettings":[],"stateVersion":3,"time":"notime","xAxisOption":"_NOTHING","playDuration":15,"iconType":"BUBBLE","sizeOption":"_NOTHING","xZoomedDataMin":null,"xZoomedIn":false,"duration":{"multiplier":1,"timeUnit":"none"},"yZoomedDataMin":null,"xLambda":1,"colorOption":"_NOTHING","nonSelectedAlpha":0.4,"dimensions":{"iconDimensions":[]},"yZoomedIn":false,"yAxisOption":"_NOTHING","yLambda":1,"yZoomedDataMax":null,"showTrails":false,"xZoomedDataMax":null};';
		# "iconKeySettings":[{"key":{"dim0":"Cyp1b1"},"trailStart":"1901"}]

		# call googlevis and make motionChart
		plot.motion = googleVis::gvisMotionChart(
			data = data.to.plot, 
			idvar = idvar, 
			timevar = "time", 
			options = list(
				gvis.editor="Editor",
				height=700,
				width=900,
				showChartButtons = TRUE,
				showHeader = TRUE,
				showSelectComponent = TRUE,
				showSidePanel = TRUE,
				showMetrixPicker = TRUE,
				showYMetricPicker = TRUE,
				showXScalePicker = TRUE,
				showYScalePicker = TRUE,
				showAdvancedPanel = TRUE,
				state = initial.plotting.parameters
				)
			);
		
		# create a data table containing the plotted data
		data.to.plot$time <- NULL;

		plot.table <- googleVis::gvisTable(
			data = data.to.plot,
			options=list(width=600, height=700)
			);

		# merge the chart and table
		plot.merge <- googleVis::gvisMerge(plot.motion, plot.table, horizontal = TRUE);

		if (save.plot == TRUE) {
			# don't plot just save

			# create a directory to dump the html files
			if (!file.exists(output.directory))  {
				dir.create(output.directory);
				}	
			
			# copy the mongoose embedded web server.  
			# note: mongoose was written under MIT licence http://code.google.com/p/mongoose/
			# if no path to mongoose then download the mongoose executables and put them in the report directory
			
			if ((path.to.mongoose) == "web") {
				download.file(url = "http://mongoose.googlecode.com/files/mongoose-3.0.exe", destfile = paste(output.directory,"/mongoose.exe",sep = ""));
				download.file(url = "http://mongoose.googlecode.com/files/mongoose-3.0.tgz", destfile = paste(output.directory, "/mongoose.tgz", sep = ""));
				cat("Plot.NanoStringNorm.gvis: Note that only the source code was downloaded for non windows systems.  You will have to untar and compile the code.  See the docs.\n");
				}
			else if (path.to.mongoose != "none") {
				if (file.exists(paste(path.to.mongoose, "/mongoose", sep = ""))) { 
					file.copy(from = paste(path.to.mongoose, "/mongoose"), to = output.directory);
					}
				else {
					cat("Plot.NanoStringNorm: No linux mongoose binary found.");
					}
				if (file.exists(paste(path.to.mongoose, "/mongoose.exe", sep = ""))) { 
					file.copy(from = paste(path.to.mongoose, "/mongoose.exe"), to = output.directory);
					}
				else {
					cat("Plot.NanoStringNorm: No windows mongoose binary found.\n");
					}
				}
			else {
				cat("Plot.NanoStringNorm: No mongoose binary found.  You will need to either download this or \n\t use an alternate method to display the interactive googleVis plots in the future.\n");
				}

			# Create Google Gadget
			#cat(createGoogleGadget(plot.merge), file = paste("NanoStringNorm_gvis_", plot.type ,"_summary.html", sep = ""))
			
			cat(plot.merge$html$chart, file=paste(output.directory, "/NanoStringNorm_gvis_", plot.item ,"_summary.html", sep = ""));

			# print message about links
			cat("Plot.NanoStringNorm.gvis: First run the mongoose binary found in the NanoStringNorm_gvis_plots and \n\t then navigate to http://127.0.01:8080 in your browser to view the plots\n");
			}
		else {

			# resets the browser setting
#			if (getOption("browser") == "/usr/bin/open" | getOption("browser") == "") options(browser = "xdg-open");
			if (getOption("browser", "/usr/bin/open") == "/usr/bin/open") options(browser = "xdg-open");

			# plot using internal R webserver
			plot(plot.merge);
			}


		}
	}
