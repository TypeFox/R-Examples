# The bedr package is copyright (c) 2014 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

.onAttach <- function(libname, pkgname) {

	if (file.exists(system.file("config/config.yml", package = "bedr"))) {
		config.bedr <- yaml::yaml.load_file(input = system.file("config/config.yml", package = "bedr"));
		}
	else {
		config.bedr <- list();
		}
	if (file.exists(paste0(Sys.getenv("HOME"), "/bedr/config.yml"))) {
		config.bedr.user <- yaml::yaml.load_file(input = paste0(Sys.getenv("HOME"), "/bedr/config.yml"));
		config.bedr <- modifyList2(config.bedr, config.bedr.user);
		}

	packageStartupMessage(
		paste0("\n\n######################\n"),
		paste0("#### bedr v", utils::packageVersion("bedr"), " ####\n"),
		paste0("######################\n\n"),
		"checking binary availability...\n",
		paste(utils::capture.output(invisible(check.binary("bedtools"))), collapse = "\n"), "\n",
		paste(utils::capture.output(invisible(check.binary("bedops"))), collapse = "\n"), "\n",
		paste(utils::capture.output(invisible(check.binary("tabix"))), collapse = "\n"), "\n",
		"tests and examples will be skipped on R CMD check if binaries are missing\n"
		);
	}

