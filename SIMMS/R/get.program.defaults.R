get.program.defaults <- function(networks.database = "default") {

	# template file to search for
	entry.secret = "pathway_based_sub_networks.txt";	

	# make a list of potential locations for the datasets file
	program.data.dirs <- paste(.libPaths(), "/SIMMS/programdata/networkdb/", networks.database, "/", sep = "");

	# then search all locations
	file.checks <- file.exists( paste(program.data.dirs, "/", entry.secret, sep = "") );

	# check to see if the file was actually found
	if (any(file.checks)) {
		program.data.dir <- program.data.dirs[ order(file.checks, decreasing = TRUE)[1] ];
		}
	else {
		stop("Unable to find pathway_based_sub_networks.txt file");
		}
	
	return (
		list(
			"program.data.dir" = program.data.dir,
			"subnets.file" = paste(program.data.dir, "/pathway_based_sub_networks.txt", sep=""),
			"subnets.file.flattened" = paste(program.data.dir, "/pathway_based_networks_flattened.txt", sep=""),
			"subnets.file.all" = paste(program.data.dir, "/pathway_based_sub_networks_all.txt", sep=""), # only used by BL.pipeline.SIMMS
			"test.data.dir" = paste(program.data.dir, "../../testdata/", sep="")
			)
		);
	}
