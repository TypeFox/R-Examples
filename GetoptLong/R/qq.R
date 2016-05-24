# == title
# Simple variable interpolation in texts
#
# == param
# -text         text string in which variables are marked with certain rules
# -envir        environment where to look for variables. By default it is the environment
#               where `qq` is envoked. It can also be a list in which element names are
#               the variable names to be interpolated.
# -code.pattern pattern of marks for the variables. By default it is ``@\\\\{CODE\\\\}`` which means
#               you can write your variable as ``@{variable}``.
# -collapse     If variables return vector of length larger than one, whether collapse into one string
#               or return a vector
#
# == details
# I like variable interpolation in Perl. But in R, if you want to concatenate plain text and variables,
# you need to use functions such as `base::paste`. However, if there are so many variables, quotes, braces 
# in the string you want to construct, it would be painful.
#   
# This function allows you to construct strings as in Perl style. Variables are
# marked in the text with certain rule. `qq` will look up these variables in user's
# environment and replace the variable marks with their real values.
#
# For more explaination of this function, please refer to vignette.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# a = 1
# b = "text"
# qq("a = @{a}, b = '@{b}'")
#
# a = 1:2
# qq("a = @{a}, b = '@{b}'")
# qq("a = @{a}, b = '@{b}'", collapse = FALSE)
#
# a = 1
# qq("a = `a`, b = '`b`'", code.pattern = "`CODE`")
#
qq = function(text, envir = parent.frame(), code.pattern = NULL, collapse = TRUE) {
	
	if(is.null(code.pattern)) {
		if(!is.null(options("code.pattern")[[1]])) {
			code.pattern = options("code.pattern")[[1]]
		} else {
			code.pattern = qq.options("code.pattern")
		}
	}
	
    if(length(text) != 1) {
        stop("Now only support text with length of 1.\n")
    }
	
    if(!is.null(envir)) {
		if(is.environment(envir)) {
			e = envir
		} else {
			e = as.environment(envir)
		}
    } else {
        e = .GlobalEnv
    }
	
#    for (i in 1:length(text)) {
 
        # check wether there are code replacements
        fc = find_code(code.pattern, text)
        template = fc[[1]]
        code = fc[[2]]
 
        if(length(template)) {   # if there is code replacement
           
            # replace the code with its value
            return_value = lapply(code, function(c) {
				x = eval(parse(text = c), envir = e)
                if(is.null(x)) {
                    return("")
				} else if(length(x) == 0) {
					return("")
				} else {
					if(is.factor(x)) {
						x = as.vector(x)
					}
					return(x)
				}
			})  # anony function is the first level parent
			
			is_return_value_vector = sapply(return_value, function(r) is.atomic(r))
			if(! all(is_return_value_vector)) {
				stop("All your codes should only return simple atomic vectors.\n")
			}
			
            # length of the return value
            # need to test it since not all code returns a scalar
            return_value_length = sapply(return_value, function(x) length(x))
 
            if(max(return_value_length) > 1) {
            # if some piece of codes return a vector
            # recycle to `max(return_value_length)`
                current = rep(text, max(return_value_length))
                
                for(ai in 1:max(return_value_length)) {
                    for(iv in 1:length(template)) {
                        current[ai] = gsub(template[iv],
                                           return_value[[iv]][(ai-1) %% length(return_value[[iv]]) + 1],
                                           current[ai], fixed = TRUE)
                    }
                }
                
				if(collapse) {
					text = paste(current, collapse = "")
				} else {
					text = current
				}
                
            } else if(max(return_value_length == 1)) {   # all variable returns a scalar
            
                current = text
                
                for(iv in 1:length(template)) {
                    current = gsub(template[iv], return_value[[iv]], current, fixed = TRUE)
                }
                
                text = current
                
            }
        }
#    }
 	
	return(text)
}
 
find_code = function(m, text) {
 
    if(length(text) != 1) {
        stop("text must be length of 1.")
    }
 
    m2 = gsub("CODE", ".+?", m)
 
    reg = gregexpr(m2, text)[[1]]
    v1 = character(0)
    if(reg[1] > -1) {
        v1 = sapply(1:length(reg), function(i) substr(text, as.numeric(reg)[i], as.numeric(reg)[i]+ attr(reg, "match.length")[i] - 1))
    }
    edge = strsplit(m, "CODE")[[1]]
    v2 = gsub(paste("^", edge[1], "|", edge[2], "$", sep=""), "", v1)
    
    return(list(template=v1, code=v2))
}


# == title
# Print a string which has been intepolated with variables
#
# == param
# -text         text string in which variables are marked with certain rules
# -envir          environment where to look for those variables
# -code.pattern pattern of marks for the variables
# -file        pass to `base::cat`
# -sep         pass to `base::cat`
# -fill        pass to `base::cat`
# -labels      pass to `base::cat`
# -append      pass to `base::cat`
# -cat_prefix  prefix string. It is prior than ``qq.options(cat_prefix)``.
#
# == details
# This function is a shortcut of
#
#     cat(qq(text, envir, code.pattern), ...)
#
# Additionally, you can add global prefix:
#
#     qq.options("cat_prefix" = "[INFO] ")
#     qq.options("cat_prefix" = function(x) format(Sys.time(), "[\%Y-\%m-\%d \%H:\%M:\%S] "))
#     qq.options("cat_prefix" = NULL)
#
# You can also add local prefix by specifying ``cat_prefix`` in `qqcat`.
#
#     qqcat(text, cat_prefix = "[INFO] ")
#
# Please refer to `qq` to find more details.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# a = 1
# b = "text"
# qqcat("a = @{a}, b = '@{b}'\n")
# qqcat("a = `a`, b = '`b`'\n", code.pattern = "`CODE`")
#
# qq.options("cat_prefix" = function(x) format(Sys.time(), "[\%Y-\%m-\%d \%H:\%M:\%S] "))
# qqcat("a = @{a}, b = '@{b}'\n")
# Sys.sleep(2)
# qqcat("a = @{a}, b = '@{b}'\n")
# qq.options(RESET = TRUE)
qqcat = function(text, envir = parent.frame(), code.pattern = NULL, file = "",
    sep = " ", fill = FALSE, labels = NULL, append = FALSE, cat_prefix = NULL) {
	cat(qq(text, envir, code.pattern), file = file, sep = sep, fill = fill, labels = labels, append = append, cat_prefix = cat_prefix)
}

cat = function(... , file = "", sep = " ", fill = FALSE, labels = NULL, append = FALSE, cat_prefix = NULL) {
    
	if(!is.null(options("cat_verbose")[[1]])) {
		if(!options("cat_verbose")[[1]]) {
			return(invisible(NULL))
		}
	}
	
	if(!qq.options("cat_verbose")) {
		return(invisible(NULL))
	}
	
	if(is.null(cat_prefix)) {
		if(!is.null(options("cat_prefix")[[1]])) {
			cat_prefix = options("cat_prefix")[[1]]
		} else {
			cat_prefix = qq.options("cat_prefix")
		}		
	}
    
	if(is.function(cat_prefix)) {
		cat_prefix = cat_prefix()
	}
	
    base::cat(cat_prefix, file = file, sep = sep, fill = fill, labels = labels, append = append)
    base::cat(... , file = file, sep = sep, fill = fill, labels = labels, append = append)
}
