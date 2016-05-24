#' @rdname CharacterClasses
#' @include regex-methods.R
#' @include concatenation.R
#' @include class-groups.R
#' @export
ASCII_ALPHA <- ASCII_LOWER %R% ASCII_UPPER

#' @rdname CharacterClasses
#' @export
ASCII_ALNUM <- ASCII_ALPHA %R% ASCII_DIGIT

#' @rdname CharacterClasses
#' @export
UNMATCHABLE <- END %R% "a"
