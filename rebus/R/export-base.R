# From rebus.base alternation.R

#' Alternation
#' 
#' See \code{\link[rebus.base]{or}}.
#' @name or
#' @aliases or1 %|%
#' @export or
#' @export or1
#' @export %|%
NULL 


# From rebus.base backreferences.R

#' Backreferences
#' 
#' See \code{\link[rebus.base]{Backreferences}}.
#' @name Backreferences
#' @aliases REF1 REF2 REF3 REF4 REF5 REF6 REF7 REF8 REF9 ICU_REF1 ICU_REF2 ICU_REF3 ICU_REF4 ICU_REF5 ICU_REF6 ICU_REF7 ICU_REF8 ICU_REF9
#' @export REF1
#' @export REF2
#' @export REF3
#' @export REF4
#' @export REF5
#' @export REF6
#' @export REF7
#' @export REF8
#' @export REF9
#' @export ICU_REF1
#' @export ICU_REF2
#' @export ICU_REF3
#' @export ICU_REF4
#' @export ICU_REF5
#' @export ICU_REF6
#' @export ICU_REF7
#' @export ICU_REF8
#' @export ICU_REF9
NULL 

#' Make the regular expression recursive.
#' 
#' See \code{\link[rebus.base]{recursive}}.
#' @name recursive
#' @export recursive
NULL


# From rebus.base capture.R

#' Capture a token, or not
#' 
#' See \code{\link[rebus.base]{capture}}.
#' @name capture
#' @aliases capture group engroup
#' @export capture
#' @export group
#' @export engroup
NULL 


# From rebus.base class-groups.R

#' Character classes
#' 
#' See \code{\link[rebus.base]{ClassGroups}}.
#' @name ClassGroups
#' @aliases alnum alpha blank cntrl digit graph lower printable punct space upper hex_digit any_char dgt wrd spc not_dgt not_wrd not_spc ascii_digit ascii_lower ascii_upper ascii_alpha ascii_alnum char_range
#' @export alnum
#' @export alpha
#' @export blank
#' @export cntrl
#' @export digit
#' @export graph
#' @export lower
#' @export printable
#' @export punct
#' @export space
#' @export upper
#' @export hex_digit
#' @export any_char
#' @export dgt
#' @export wrd
#' @export spc
#' @export not_dgt
#' @export not_wrd
#' @export not_spc
#' @export ascii_digit
#' @export ascii_lower
#' @export ascii_upper
#' @export ascii_alpha
#' @export ascii_alnum
#' @export char_range
NULL 


# From rebus.base compound-constants.R, constants.R

#' Special characters
#' 
#' See \code{\link[rebus.base]{SpecialCharacters}}.
#' @name SpecialCharacters
#' @aliases BACKSLASH CARET DOLLAR DOT PIPE QUESTION STAR PLUS OPEN_PAREN CLOSE_PAREN OPEN_BRACKET OPEN_BRACE
#' @export BACKSLASH
#' @export CARET
#' @export DOLLAR
#' @export DOT
#' @export PIPE
#' @export QUESTION
#' @export STAR
#' @export PLUS
#' @export OPEN_PAREN
#' @export CLOSE_PAREN
#' @export OPEN_BRACKET
#' @export OPEN_BRACE
NULL 

#' The start or end of a string
#' 
#' See \code{\link[rebus.base]{Anchors}}.
#' @name Anchors
#' @aliases START END
#' @export START
#' @export END
NULL 

#' Class Constants
#' 
#' See \code{\link[rebus.base]{CharacterClasses}}.
#' @name CharacterClasses
#' @aliases ALNUM ALPHA BLANK CNTRL DIGIT GRAPH LOWER PRINT PUNCT SPACE UPPER HEX_DIGIT ANY_CHAR DGT WRD SPC NOT_DGT NOT_WRD NOT_SPC ASCII_DIGIT ASCII_LOWER ASCII_UPPER ASCII_ALPHA ASCII_ALNUM UNMATCHABLE
#' @export ALNUM
#' @export ALPHA
#' @export BLANK
#' @export CNTRL
#' @export DIGIT
#' @export GRAPH
#' @export LOWER
#' @export PRINT
#' @export PUNCT
#' @export SPACE
#' @export UPPER
#' @export HEX_DIGIT
#' @export ANY_CHAR
#' @export DGT
#' @export WRD
#' @export SPC
#' @export NOT_DGT
#' @export NOT_WRD
#' @export NOT_SPC
#' @export ASCII_DIGIT
#' @export ASCII_LOWER
#' @export ASCII_UPPER
#' @export ASCII_ALPHA
#' @export ASCII_ALNUM
#' @export UNMATCHABLE
NULL 

#' Word boundaries
#' 
#' See \code{\link[rebus.base]{WordBoundaries}}.
#' @name WordBoundaries
#' @aliases BOUNDARY NOT_BOUNDARY
#' @export BOUNDARY
#' @export NOT_BOUNDARY
NULL 


# From rebus.base concatenation.R

#' Combine strings together
#' 
#' See \code{\link[rebus.base]{Concatenation}}.
#' @name Concatenation
#' @aliases %R% %c%
#' @export %R%
#' @export %c%
NULL 

#' Escape special characters
#' 
#' See \code{\link[rebus.base]{escape_special}}.
#' @name escape_special
#' @export escape_special
NULL 


# From rebus.base grouping-and-repetition.R

#' A range or char_class of characters
#' 
#' See \code{\link[rebus.base]{char_class}}.
#' @name char_class
#' @aliases negated_char_class
#' @export char_class
#' @export negated_char_class
NULL 

#' Repeat values
#' 
#' See \code{\link[rebus.base]{repeated}}.
#' @name repeated
#' @aliases optional lazy zero_or_more one_or_more
#' @export repeated
#' @export optional
#' @export lazy
#' @export zero_or_more
#' @export one_or_more
NULL 


# From rebus.base lookaround.R

#' Lookaround
#' 
#' See \code{\link[rebus.base]{lookahead}}.
#' @name lookahead
#' @aliases negative_lookahead lookbehind negative_lookbehind
#' @export lookahead 
#' @export negative_lookahead
#' @export lookbehind
#' @export negative_lookbehind
NULL 


# From rebus.base misc.R

#' Make a regex exact
#' 
#' See \code{\link[rebus.base]{exactly}}.
#' @name exactly
#' @export exactly
NULL 

#' Treat part of a regular expression literally
#' 
#' See \code{\link[rebus.base]{literal}}.
#' @name literal
#' @export literal
NULL 


# From rebus.base mode-modifiers.R

#' Apply mode modifiers
#' 
#' See \code{\link[rebus.base]{modify_mode}}.
#' @name modify_mode
#' @aliases case_insensitive free_spacing single_line multi_line duplicate_group_names no_backslash_escaping
#' @export modify_mode
#' @export case_insensitive
#' @export free_spacing
#' @export single_line
#' @export multi_line 
#' @export duplicate_group_names
#' @export no_backslash_escaping
NULL 


# From rebus.base regex-methods.R

#' Convert or test for regex objects
#' 
#' See \code{\link[rebus.base]{as.regex}}.
#' @name as.regex
#' @aliases is.regex
#' @export as.regex
#' @export is.regex
NULL 

#' Create a regex
#' 
#' See \code{\link[rebus.base]{regex}}.
#' @name regex
#' @export regex
NULL 

#' Print or format regex objects
#' 
#' See \code{\link[rebus.base]{format.regex}}.
#' @name format.regex
#' @aliases print.regex
NULL 
