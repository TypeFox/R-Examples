# == title
# Wrapper of the Perl module ``Getopt::Long`` in R
#
# == param
# -...     specification of options. The value should be a vector having even number of elements.
# -help     whether to add help option
# -version  whether to add version option
# -envir    user's enrivonment where `GetoptLong` will look for default values and export variables
# -argv_str command-line arguments, only for testing purpose
#
# == details
# Following shows a simple example. Put following code at the beginning of your script (e.g. ``foo.R``):
#
#     library(GetoptLong)
#     cutoff = 0.05
#     GetoptLong(
#         "number=i", "Number of items, integer, mandatory option",
#         "cutoff=f", "cutoff to filter results, optional, default (0.05)",
#         "verbose",  "print messages"
#     )
#
# Then you can call the script from command line either by:
#
#     ~\> Rscript foo.R --number 4 --cutoff 0.01 --verbose
#     ~\> Rscript foo.R -n 4 -c 0.01 -v
#     ~\> Rscript foo.R -n 4 --verbose
#
# In above example, ``number`` is a mandatory option and should only be integer mode. ``cutoff``
# is optional and already has a default value. ``verbose`` is a logical option. If parsing is
# successful, two variables with name ``number`` and ``verbose`` will be imported into the working
# environment with specified values, and value for ``cutoff`` will be updated if it is specified in
# command-line argument.
#
# For advanced use of this function, please go to the vignette.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
GetoptLong = function(..., help = TRUE, version = TRUE, envir = parent.frame(), argv_str = NULL) {
	
	spec = list(...)

	# a vector or a two-column matrix
	if(length(spec) == 1) {
		spec = spec[[1]]
	} else {
		spec = unlist(spec)
	}

	if(is.null(get_scriptname())) {
		.IS_UNDER_COMMAND_LINE = FALSE
	} else {
		.IS_UNDER_COMMAND_LINE = TRUE
	}

	if(is.null(argv_str)) {
		argv_str = GetoptLong.options("__argv_str__")
	}
	on.exit(GetoptLong.options("__argv_str__" = NULL))

	# to test whether the script is run under command-line or in R interactive environment
	if(.IS_UNDER_COMMAND_LINE || is.null(argv_str)) {
		OUT = stderr()
	} else {
		OUT = stdout()  # message from STDOUT is important under testing mode
	}
	
	# get the path of binary perl
	# it will look in PATH and also user's command-line argument
	perl_bin = find_perl_bin(con = OUT, from_command_line = .IS_UNDER_COMMAND_LINE)
	
	# check whether `perl` is real Perl,
	# in fact, this step is not necessary
	if(!check_perl(perl_bin = perl_bin)) {
		qqcat("Error when testing Perl: @{perl_bin}.\n", file = OUT)
		if(.IS_UNDER_COMMAND_LINE) {
			q(save = "no", status = 127)
		} else if(!is.null(argv_str)) {  # under test
			return(invisible(NULL))
		} else {
			stop("You have an error.\n")
		}
	}
	
	# check whether Getopt::Long is in @INC
	# normally, it is shippped with standard Perl distributions
	if(!check_perl("Getopt::Long", perl_bin = perl_bin)) {
		cat("Cannot find Getopt::Long module in your Perl library.\n", file = OUT)
		if(.IS_UNDER_COMMAND_LINE) {
			q(save = "no", status = 127)
		} else if(!is.null(argv_str)) {  # under test
			return(invisible(NULL))
		} else {
			stop("You have an error.\n")
		}
	}
	
	# check whether JSON is in @INC
	if(!check_perl("JSON", inc = qq("@{system.file('extdata', package='GetoptLong')}/perl_lib"), perl_bin = perl_bin)) {
		cat("Cannot find JSON module in your Perl library.\n", file = OUT)
		if(.IS_UNDER_COMMAND_LINE) {
			q(save = "no", status = 127)
		} else if(!is.null(argv_str)) {  # under test
			return(invisible(NULL))
		} else {
			stop("You have an error.\n")
		}
	}

	# check first argument
	# it should be a matrix with two columns or a vector with even number of elements
	if(is.matrix(spec)) {
		if(ncol(spec) != 2) {
			stop("If your specification is a matrix, it should be a two-column matrix.\n")
		}
	} else {
		if(is.vector(spec)) {
			if(length(spec) %% 2) {
				stop("Since your specification is a vector, it should have even number of elements.\n")
			} else {
				spec = matrix(spec, ncol = 2, byrow = TRUE)
			}
		} else {
			stop("Wrong specification class.\n")
		}
	}
	
	# we use another way to implement 'optional' options
	if(any(detect_optional(spec[, 1]))) {
		stop("type :[isfo] is not allowed, use =[isfo] instead.\n")
	}
	
	# get arguments string
	if(is.null(argv_str)) {
		ARGV = commandArgs(TRUE)
		ARGV_string = combine_and_escape_ARGV(ARGV)
	} else {
		ARGV_string = argv_str
	}
	
	# first name in each options
	long_name = extract_first_name(spec[, 1])
	if(help && long_name %in% "help") {
		stop("`help` is reserved as a default option, please do not use it.\n")
	}
	
	if(version && long_name %in% "version") {
		stop("`version` is reserved as default option, please do not use it.\n")
	}
	
	# test whether first name in option name is a valid R variable name
	test_long_name = grepl("^[a-zA-Z_\\.][a-zA-Z0-9_\\.]*$", long_name) 
		
	if(!all(test_long_name)) {
		cat("First name in option names can only be valid R variable names which only\nuse numbers, letters,'.' and '_' (It should match\n/^[a-zA-Z_\\.][a-zA-Z0-9_\\.]+$/).", file = OUT)
		qqcat(" Following option name@{ifelse(sum(test_long_name)==1, ' is', 's are')}\nnot valid:\n\n", file = OUT)
		for(k in seq_along(test_long_name)) {
			if(!test_long_name[k]) qqcat("  @{spec[k, 1]}\n", file = OUT)
		}
		cat("\n", file = OUT)
			
		if(.IS_UNDER_COMMAND_LINE) {
			q(save = "no", status = 127)
		} else if(!is.null(argv_str)) {  # under test
			return(invisible(NULL))
		} else {
			stop("You have an error.\n")
		}
	}
	
	json_file = tempfile(fileext = ".json")
	spec2 = spec
	if(help) {
		spec2 = rbind(spec2, c("help", ""))
	}
	if(version) {
		spec2 = rbind(spec2, c("version", ""))
	}
	perl_script = generate_perl_script(spec2, json_file)
	
	cmd = qq("\"@{perl_bin}\" \"@{perl_script}\" @{ARGV_string}")
	res = system(cmd, intern = TRUE)
	res = as.vector(res)
	
	# if you specified wrong arguments
	if(length(res)) {
		qqcat("@{res}\n", file = OUT)
		
		if(.IS_UNDER_COMMAND_LINE) {
			print_help_msg(spec, file = OUT, help = TRUE, version = TRUE)
		}

		ow = options("warn")[[1]]
		options(warn = -1)
		file.remove(json_file)
		file.remove(perl_script)
		options(warn = ow)
		
		if(.IS_UNDER_COMMAND_LINE) {
			q(save = "no", status = 127)
		} else if(!is.null(argv_str)) {  # under test
			return(invisible(NULL))
		} else {
			stop("You have an error.\n")
		}
	}
	
	# if arguments are correct, values for options will be stored in .json file
	opt = fromJSON(file = json_file)
	file.remove(json_file)
	file.remove(perl_script)
	
	# if detect user has specified --help or --version
	# basically, !is.null(opt$help) measn opt$help == 1
	if(!is.null(opt$help) && opt$help) {
		print_help_msg(spec, file = OUT, help = help, version = version)
		
		if(.IS_UNDER_COMMAND_LINE) {
			q(save = "no", status = 127)
		} else if(!is.null(argv_str)) {  # under test
			return(invisible(NULL))
		} else {
			stop("You have an error.\n")
		}
	}
	
	if(!is.null(opt$version) && opt$version) {
		print_version_msg(envir, file = OUT)
		
		if(.IS_UNDER_COMMAND_LINE) {
			q(save = "no", status = 127)
		} else if(!is.null(argv_str)) {  # under test
			return(invisible(NULL))
		} else {
			stop("You have an error.\n")
		}
	}
	
	# check mandatory options
	# note here `spec` does not contain `help`` or `version`
	is_mandatory = detect_mandatory(spec[, 1])
	for(i in seq_len(nrow(spec))) {
		# if variable not defined, or defined as a function
		if(is.null(opt[[ long_name[i]] ]) && is_mandatory[i] &&
		   (!exists(long_name[i], envir = envir) || class(get(long_name[i], envir = envir)) == "function")) {
			qqcat("@{long_name[i]} is mandatory, please specify it.\n", file = OUT)
			if(.IS_UNDER_COMMAND_LINE) {
				print_help_msg(spec, file = OUT, help = help, version = version)
			}

			if(.IS_UNDER_COMMAND_LINE) {
				q(save = "no", status = 127)
			} else if(!is.null(argv_str)) {  # under test
				return(invisible(NULL))
			} else {
				stop("You have an error.\n")
			}
		}
	}
	
	# check default values
	for(i in seq_len(nrow(spec))) {
		if(is_mandatory[i] && exists(long_name[i], envir = envir)) {
			tmp = get(long_name[i], envir = envir)
			
			# if mode is function, skip it
			if(mode(tmp) == "function") {
				next
			}
			
			if(is.null(tmp)) {
				next
			}
			
			# if option is specified as a list (ss=%)
			if(grepl("%$", spec[i, 1])) {
				# if default value is not a list
				if(!is.list(tmp)) {
					qqcat("@{long_name[i]} is mandatory, and also detect in evoking environment you have already \ndefined `@{long_name[i]}`. Since it is defined as a named option, please\nmake sure default value of `@{long_name[i]}` is a list.\n", file = OUT)
					if(.IS_UNDER_COMMAND_LINE) {
						print_help_msg(spec, file = OUT, help = help, version = version)
					}

					if(.IS_UNDER_COMMAND_LINE) {
						q(save = "no", status = 127)
					} else if(!is.null(argv_str)) {  # under test
						return(invisible(NULL))
					} else {
						stop("You have an error.\n")
					}
				} else if(is.null(names(tmp))) {
					qqcat("@{long_name[i]} is mandatory, and also detect in evoking environment you have already \ndefined `@{long_name[i]}`. Since it is defined as a named option, please\nmake sure default value of `@{long_name[i]}` is a list with names.\n", file = OUT)
					if(.IS_UNDER_COMMAND_LINE) {
						print_help_msg(spec, file = OUT, help = help, version = version)
					}

					if(.IS_UNDER_COMMAND_LINE) {
						q(save = "no", status = 127)
					} else if(!is.null(argv_str)) {  # under test
						return(invisible(NULL))
					} else {
						stop("You have an error.\n")
					}
				} else if(!(all(sapply(tmp, is_simple_vector)))) {
					qqcat("@{long_name[i]} is mandatory, and also detect in evoking environment you have already \ndefined `@{long_name[i]}`. Since it is defined as a named option, please\nmake sure default value of `@{long_name[i]}` is a list containing simple vectors.\n", file = OUT)
					if(.IS_UNDER_COMMAND_LINE) {
						print_help_msg(spec, file = OUT, help = help, version = version)
					}

					if(.IS_UNDER_COMMAND_LINE) {
						q(save = "no", status = 127)
					} else if(!is.null(argv_str)) {  # under test
						return(invisible(NULL))
					} else {
						stop("You have an error.\n")
					}
				}
			} else if(! is_simple_vector(tmp)) {
				qqcat("@{long_name[i]} is mandatory, and also detect in evoking environment you have already \ndefined `@{long_name[i]}`. Please make sure default value of `@{long_name[i]}` should only be a simple vector or NULL.\n", file = OUT)
				if(.IS_UNDER_COMMAND_LINE) {
					print_help_msg(spec, file = OUT, help = help, version = version)
				}

				if(.IS_UNDER_COMMAND_LINE) {
					q(save = "no", status = 127)
				} else if(!is.null(argv_str)) {  # under test
					return(invisible(NULL))
				} else {
					stop("You have an error.\n")
				}
			}
		}
	}
	
	# export to envir
	export_parent_env(opt, envir = envir)
	
	return(invisible(opt))
}

# == title
# Wrapper of the Perl module ``Getopt::Long`` in R
#
# == param
# -... pass to `GetoptLong`
#
# == details
# This function is the same as `GetoptLong`. It is just to make it consistent as the ``GetOptions()`` 
# subroutine in ``Getopt::Long`` module in Perl.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
GetOptions = function(...) GetoptLong(...)

is_simple_vector = function(obj) {
	mode(obj) %in% c("numeric", "character", "NULL", "logical", "complex")
}

# this function will be improved later
# it will check whether there are blacks and quotes in users' values
combine_and_escape_ARGV = function(ARGV) {
	paste(ARGV, collapse = " ") 
}

generate_perl_script = function(spec, json_file) {
	perl_script = tempfile(fileext = ".pl")
	#perl_script = "tmp.pl"
	
	long_name = extract_first_name(spec[, 1])

	var_type = detect_var_type(spec[, 1])  # which is scalar, array and hash
	opt_type = detect_opt_type(spec[, 1])  # which is integer, numeric, character, ...
	
	# construct perl code
	perl_code = NULL
	perl_code = c(perl_code, qq("#!/usr/bin/perl"))
	perl_code = c(perl_code, qq(""))
	perl_lib = qq("@{system.file('extdata', package='GetoptLong')}/perl_lib")
	perl_code = c(perl_code, qq("BEGIN { push (@INC, '@{perl_lib}'); }"))
	perl_code = c(perl_code, qq(""))
	perl_code = c(perl_code, qq("use strict;"))
	
	config = NULL
	if(!is.null(options("GetoptLong.Config")[[1]]) && is.null(GetoptLong.options("config"))) {
		config = options("GetoptLong.Config")[[1]]
	} else {
		config = GetoptLong.options("config")
	}
	if(is.null(config)) {
		perl_code = c(perl_code, qq("use Getopt::Long;"))
	} else {
		perl_code = c(perl_code, qq("use Getopt::Long qw(:config @{paste(config, collapse = ' ')});"))
	}
	
	perl_code = c(perl_code, qq("use JSON;"))
	perl_code = c(perl_code, qq("use Data::Dumper;"))
	perl_code = c(perl_code, qq(""))
	
	# declare variables according to variable types
	for (i in seq_len(nrow(spec))) {
	
		perl_code = c(perl_code, qq("my @{perl_sigil(var_type[i])}opt_@{i};    # var_type = @{var_type[i]}, opt_type = @{opt_type[i]}"))
		
	}
	
	perl_code = c(perl_code, qq("*STDERR = *STDOUT;")) # all write to STDOUT
	perl_code = c(perl_code, qq("my $is_successful = GetOptions("))
	
	for (i in seq_len(nrow(spec))) {
		perl_code = c(perl_code, qq("    '@{spec[i, 1]}' => \\@{perl_sigil(var_type[i])}opt_@{i},"))
	}
	perl_code = c(perl_code, qq(");"))
	
	perl_code = c(perl_code, qq("if(!$is_successful) {"))
	perl_code = c(perl_code, qq("    exit;"))
	perl_code = c(perl_code, qq("}"))
	
	# if var_type == integer or numberic, value should be forced ensured
	for (i in seq_len(nrow(spec))) {
		if(opt_type[i] %in% c("integer", "numeric")) {
			
			if(var_type[i] == "scalar") {
			
				perl_code = c(perl_code, qq("if(defined(@{perl_sigil(var_type[i])}opt_@{i})) {"))
				perl_code = c(perl_code, qq("    $opt_@{i} += 0;"))
				perl_code = c(perl_code, "}")
				
			} else if(var_type[i] == "array") {
				
				# if array is defined
				perl_code = c(perl_code, qq("if(@opt_@{i}) {"))
				perl_code = c(perl_code, qq("    foreach (@opt_@{i}) { $_ += 0; }"))
				perl_code = c(perl_code, "}")
				
			} else if(var_type[i] == "hash") {
				
				# if hash is defined
				perl_code = c(perl_code, qq("if(%opt_@{i}) {"))
				perl_code = c(perl_code, qq("    foreach (keys %opt_@{i}) { $opt_@{i}{$_} += 0; }"))
				perl_code = c(perl_code, "}")
				
			}
			
		}
	}
	perl_code = c(perl_code, qq(""))
	
	perl_code = c(perl_code, qq("my $all_opt = {"))
	
	for (i in seq_len(nrow(spec))) {
	
		if(opt_type[i] == "logical") {
		
			perl_code = c(perl_code, qq("    '@{long_name[i]}' => $opt_@{i} ? JSON::true : JSON::false,"))

		} else if(var_type[i] == "scalar") {
			
			perl_code = c(perl_code, qq("    '@{long_name[i]}' => $opt_@{i},"))

		} else {
			
			# in scalar content, empty list will be 0
			perl_code = c(perl_code, qq("    '@{long_name[i]}' => scalar(@{perl_sigil(var_type[i])}opt_@{i}) ? \\@{perl_sigil(var_type[i])}opt_@{i} : undef,"))

		}
	}
	
	perl_code = c(perl_code, qq("};"))
	perl_code = c(perl_code, qq(""))

	perl_code = c(perl_code, qq("open JSON, '>@{json_file}' or die 'Cannot create temp file: @{json_file}\\n';"))
	perl_code = c(perl_code, qq("print JSON to_json($all_opt, {pretty => 1});"))
	perl_code = c(perl_code, qq("close JSON;"))
	#perl_code = c(perl_code, qq("print Dumper $all_opt;"))
	
	writeLines(perl_code, perl_script)
	return(perl_script)
}

perl_sigil = function(type) {
	if(type == "scalar") {
		return("$")
	} else if(type == "array") {
		return("@")
	} else if(type == "hash") {
		return("%")
	} else {
		return("$")
	}
}

print_help_msg = function(spec, file = stderr(), help = TRUE, version = TRUE) {
	
	# add help and version options in `spec`
	if(help) {
		spec = rbind(spec, c("help", "Print help message and exit"))
	}
	if(version) {
		spec = rbind(spec, c("version", "Print version information and exit"))
	}
	
	startingMsg = NULL
	if(!is.null(options("GetoptLong.startingMsg")[[1]]) && is.null(GetoptLong.options("startingMsg"))) {
		startingMsg = options("GetoptLong.startingMsg")[[1]]
	} else {
		startingMsg = GetoptLong.options("startingMsg")
	}
	
	if(!is.null(startingMsg)) {
		cat(startingMsg, file = file)
	}
	
    script_name = get_scriptname()
    if(is.null(script_name)) {
    	script_name = "foo.R"
    } else {
    	script_name = basename(script_name)
    }
    qqcat("Usage: Rscript @{script_name} [options]\n\n", file = file)
    	
	for(i in seq_len(nrow(spec))) {
		print_single_option(spec[i, 1], spec[i, 2], file = file)
	}
	
	endingMsg = NULL
	if(!is.null(options("GetoptLong.endingMsg")[[1]]) && is.null(GetoptLong.options("endingMsg"))) {
		endingMsg = options("GetoptLong.endingMsg")[[1]]
	} else {
		endingMsg = GetoptLong.options("endingMsg")
	}
	
	if(!is.null(endingMsg)) {
		cat(endingMsg, file = file)
	}
}

print_single_option = function(opt, desc, file = stderr()) {
	var_type = detect_var_type(opt)
	opt_type = detect_opt_type(opt)
	
	opt = gsub("[$@%]$", "", opt)
	opt = gsub("\\{.*\\}$", "", opt)
	opt = gsub("[=:][siof]$", "", opt)
	opt = gsub("[!+]$", "", opt)
	
	choices = strsplit(opt, "\\|")[[1]]
	
	cat("  ", file = file)
	for(i in seq_along(choices)) {
		qqcat("@{ifelse(nchar(choices[i]) == 1, '-', '--')}@{choices[i]}@{ifelse(i == length(choices), '', ', ')}", file = file)
	}
	cat(" ", file = file)
	if(var_type == "scalar" && opt_type == "extended_integer") {
		cat("extended_integer", file = file)
	} else if(var_type == "scalar" && opt_type == "logical") {
		cat("", file = file)
	} else if(var_type == "scalar") {
		cat(opt_type, file = file)
	} else if(var_type == "array") {
		qqcat("[ @{opt_type}, ... ]", file = file)
	} else if(var_type == "hash") {
		qqcat("{ name=@{opt_type}, ... }", file = file)
	}
	
	cat("\n", file = file)
	
	cat_format_line(desc, prefix = "    ", file = file)

}

print_version_msg = function(envir, file = stderr()) {
	if(exists("VERSION", envir = envir)) {
		cat(get("VERSION", envir = envir), file = file)
	} else {
		cat("No version information is found in source code.\n", file = file)
	}
	cat("\n", file = file)
}

cat_format_line = function(text, prefix = "", max.width = 70, file = stderr()) {
	words = strsplit(text, "\\s+")[[1]]
	
	i_width = nchar(prefix)
	cat(prefix, file = file)
	for(i in seq_along(words)) {
		if(i_width + 1 + nchar(words[i]) > max.width) {
			cat("\n", file = file)
			cat(prefix, file = file)
			cat(words[i], file = file)
			i_width = nchar(prefix) + nchar(words[i])
		} else {
			cat(ifelse(i == 1, "", " "), file = file)
			qqcat("@{words[i]}", file = file)
			i_width = i_width + nchar(prefix)
		}
	}
	cat("\n\n", file = file)
}

detect_var_type = function(opt) {
	sapply(opt, function(x) {
		if (grepl("\\$$", x)) {
			return("scalar")
		} else if (grepl("@$", x)) {
			return("array")
		} else if (grepl("%$", x)) {
			return("hash")
		} else if (grepl("\\{\\d?,?\\d?\\}$", x)) {
			return("array")
		} else {
			return("scalar")
		}
	}, USE.NAMES = FALSE)
}

detect_opt_type = function(opt) {
	
	opt = gsub("[$@%]$", "", opt)
	opt = gsub("\\{.*\\}$", "", opt)
	
	sapply(opt, function(x) {
		if (grepl("[=:]s$", x)) {
			return("character")
		} else if (grepl("[=:]i$", x)) {
			return("integer")
		} else if (grepl("[=:]o$", x)) {
			return("extended_integer")
		} else if (grepl("[=:]f$", x)) {
			return("numeric")
		} else if (grepl("!$", x)) {
			return("logical")
		} else if (grepl("\\+$", x)) {
			return("integer")
		} else {
			return("logical")
		}
	}, USE.NAMES = FALSE)
}

detect_mandatory = function(opt) {
	opt = gsub("[$@%]$", "", opt)
	opt = gsub("\\{.*\\}$", "", opt)
	
	grepl("=[siof]$", opt)
}

detect_optional = function(opt) {
	opt = gsub("[$@%]$", "", opt)
	opt = gsub("\\{.*\\}$", "", opt)
	
	grepl(":[siof]$", opt)
}

extract_first_name = function(opt) {
	opt = gsub("[$@%]$", "", opt)
	opt = gsub("\\{.*\\}$", "", opt)
	opt = gsub("[=:][siof]$", "", opt)
	opt = gsub("[!+]$", "", opt)

	first_name = sapply(strsplit(opt, "\\|"), function(x) x[1])
	return(first_name)
}

export_parent_env = function(opt, envir = parent.frame()) {
	
	opt_name = names(opt)
	
	# specified from command line
	# if option is not specified and has default values, it is NULL by Getopt::Long
	specified_opt_name = opt_name[ !sapply(opt, is.null) ]
	# export to global environment
	for(o in specified_opt_name) {
		if(o == "help" || o == "version") {
			next
		}

		# if the options is a named list, those default elements that are not specified on command-line will be retained
		if(is.list(opt[[o]]) && exists(o, envir = envir)) {
			ol = get(o, envir = envir)  # default value
			ol[ names(opt[[o]]) ] = opt[[o]]
			assign(o, ol, envir = envir)
		} else {    # for simple vector, just overwrite it
			assign(o, opt[[o]], envir = envir)
		}
	}

	# defined with default values while not specified in command line
	#specified_parent_opt_name = intersect(opt_name[ sapply(opt, is.null) ], parent_opt_name)
	# already have, do nothing

	return(invisible(NULL))
	
}

# == title
# Full path of current script
#
# == value
# If the R script is not run under command-line, it return ``NULL``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
get_scriptname = function() {
	args = commandArgs()
	
	if(length(args) == 1) {
		return(NULL)
	}
	i_arg = which(args == "--args")
	if(length(i_arg) == 0) {
		return(NULL)
	}
	i_arg = i_arg[1]
	args = args[seq_len(i_arg)]
    f = grep("^--file=", args, value = TRUE)
    if(length(f)) {
    	f = gsub("^--file=(.*)$", "\\1", f[1])
    	return(f)	
    } else {
    	return(GetoptLong.options("__script_name__"))
    }  
}

# find path of binary perl
find_perl_bin = function(con = stderr(), from_command_line = TRUE) {

	# first look at user's options
	args = commandArgs()
	i = which(args == "--")
	if(length(i) && length(args) > i) {
		perl_bin = args[i + 1]
		
		if(!file.exists(perl_bin)) {
			qqcat("Cannot find @{perl_bin}\n", file = con)
			if(from_command_line) {
				q(save = "no", status = 127)
			} else {
				return(invisible(NULL))
			}
		}
		
		if(!file.info(perl_bin)$isdir) {
			qqcat("@{perl_bin} should be a file, not a directory.\n", file = con)
			if(from_command_line) {
				q(save = "no", status = 127)
			} else {
				return(invisible(NULL))
			}
		}
		
	} else {  # look at PATH
		perl_bin = Sys.which("perl")
		if(perl_bin == "") {
			cat("cannot find Perl in PATH.\n", file = con)
			if(from_command_line) {
				q(save = "no", status = 127)
			} else {
				return(invisible(NULL))
			}
		}
	}
	
	return(perl_bin)
}

# check whether perl can be called
# check whether perl has certain module
check_perl = function(module = NULL, inc = NULL, perl_bin = "perl") {
	
	if(is.null(module)) {
		cmd = qq("\"@{perl_bin}\" -v")
	} else if(!is.null(module) && is.null(inc)) {
		cmd = qq("\"@{perl_bin}\" -M@{module} -e \"use @{module}\"")
	} else if(!is.null(module) && !is.null(inc)) {
		cmd = qq("\"@{perl_bin}\" \"-I@{inc}\" -M@{module} -e \"use @{module}\"")
	}
	
	OS = Sys.info()["sysname"]
	if(OS == "Windows") {
		res = system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE, show.output.on.console = FALSE)
	} else {
		res = system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
	}

	return(ifelse(res, FALSE, TRUE))
}


is.dir = function(dir) {
	sapply(dir, function(x) file.exists(x) && file.info(x)[1, "isdir"])
}

#  title
# Read R source code with arguments
#
# == param
# -file file name
# -... pass to `base::source`
# -argv a string which contains command line arguments
#
# == details
# This function overwrites the default base::`base::source`.
#
# Command-line arguments can be used when sourcing an R file, just like running R script in the command-line.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
source = function(file, ..., argv = NULL) {
	GetoptLong.options("__argv_str__" = argv)
	GetoptLong.options("__script_name__" = file)
	base::source(file, ...)
	GetoptLong.options("__script_name__" = NULL)
	GetoptLong.options("__argv_str__" = NULL)
}

# stop = function(msg) {
# 	base::stop(strwrap(msg), call. = FALSE)
# }
