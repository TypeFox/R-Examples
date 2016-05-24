#' @title Estimate body surface area of an adult
#' @description \code{bsa_adult} Estimate body surface area of an
#'   adult using \code{sqrt(wt*ht)/6} TODO: reference for this.
#' @template heightm
#' @template weightkg
#' @return numeric vector
#' @examples
#' bsa_adult(2, 80)
#' bsa_adult(1.5, 80)
#' @export
bsa_adult <- function(heightm, weightkg)
  sqrt(heightm * weightkg) / 6

#' @title ideal weight for adults
#' @description \code{ideal_weight_adult} gives the ideal weight using default
#'   adult algorithm, Devine.
#'
#'   TODO: param to switch to different method without knowing the function
#'   name.
#' @template heightm
#' @template male
#' @template dots
#' @rdname ideal_weight
#' @examples
#' ideal_weight_adult(1.7, TRUE)
#' ideal_weight_adult(1.7, FALSE)
#' ideal_weight_adult(6 * 12 * 2.54, TRUE)
#' @export
ideal_weight_adult <- function(heightm, male, ...)
  ideal_weight_Devine(heightm, male, ...)

#' @title ideal weight for adults
#' @description \code{ideal_weight_} gives the ideal weight using default
#'   paediatric algorithm. TODO: account for Down's/other well calibrated
#'   developmental differences. Age specifications are mutually exclusive, and
#'   an error will be generated if more than on age i specified per patient.
#' @template heightm
#' @param age.years numeric vector, age(s) in years
#' @param age.months numeric vector, age(s) in months
#' @param age.days  numeric vector, age(s) in days
#' @template warn
#' @rdname ideal_weight
#' @export
ideal_weight_child <- function(heightm,
                               age.years = NULL,
                               age.months = NULL,
                               age.days = NULL,
                               warn = FALSE)
  ideal_weight_Straub(heightm, age.years, age.months, age.days, warn)

#' @title ideal weight for child per Straub
#' @description http://www.ncbi.nlm.nih.gov/pubmed/6823980 2.396e0.01863(ht),
#'   where height is in cm. There is an argument for using another package to
#'   capture durations, of which age is a special case. However, I am resisting
#'   bringing in external dependencies, and for almost all use-cases I can
#'   imagine, the age will be captured as a single number of one type, not a mix
#'   of types. Note that gender does not appear to be important in this
#'   relationship.
#'
#'   See package AGD for CDC growth chart data.
#' @inheritParams ideal_weight_child
#' @source \url{http://www.ncbi.nlm.nih.gov/pubmed/6823980}
#' @examples
#' # will warn if given age is not in validate range from publication:
#' \dontrun{
#'   ideal_weight_child(0.5, age.years = 0, warn = TRUE)
#'   ideal_weight_child(0.8, age.months = 11, warn = TRUE)
#'   ideal_weight_child(0.5, age.days = 25, warn = TRUE)
#' }
#'   ideal_weight_child(0.5, age.days = 25, warn = FALSE)
#'   ideal_weight_child(1, age.years = 2)
#'   ideal_weight_child(0.75, age.months = 15)
#' @export
ideal_weight_Straub <- function(heightm,
                                age.years = NULL,
                                age.months = NULL,
                                age.days = NULL,
                                warn = FALSE) {
  stopifnot(sum(is.null(age.years),
                is.null(age.months),
                is.null(age.days)) == 2)
  if (warn) {
    if (any(!is.null(age.years) & (age.years < 1 | age.years > 17)))
      warning("age < 1 year or age > 17 year not validated from Straub formula")
    if (any(!is.null(age.months) & (age.months < 24 | age.months > 12 * 17 + 11)))
      warning("age < 1 year or age > 17 year not validated from Straub formula")
    if (any(!is.null(age.months) & (age.days < 364 | age.days > 364 * 12 * 18)))
      warning("age < 1 year or age > 17 year not validated from Straub formula")
  }
  #
  # 2.396e0.01863(ht), where height is in cm
  2.396 ^ (1.863 * heightm)
}

#' @title ideal weight by Devine method
#' @description Devine method is the default and most widely used. Normally
#'   stated in inches. Male: 50kg + 2.3kg * inches over 5ft. Female: 45.5kg +
#'   2.3kg * inches over 5ft. (from 1974 genatamicin paper - see Lemmens for
#'   ref.)
#' @rdname ideal_weight
#' @export
ideal_weight_Devine <- function(heightm, male, ...)
  ideal_weight_linear(heightm, male, 60, 50, 45.5, 2.3, 2.3, ...)

#' @title ideal weight by Robinson method
#' @rdname ideal_weight
#' @export
#' @description Robinson's method for ideal weight: different linear
#'   relationship. (Robinson JD, Lupkiewicz SM, Palenik L et al. Determination
#'   of ideal body weight for drug dosage calculations. Am J Hosp Pharm 1983;
#'   40: 1016-9.)
ideal_weight_Robinson <- function(heightm, male, ...)
  ideal_weight_linear(heightm, male, 60, 52, 49, 1.9, 1.7, ...)

#' @title ideal weight by Miller
#' @export
#' @rdname ideal_weight
#' @description Miller's method for ideal weight: different linear relationship.
#'   (Miller DR, Carlson JD, Loyd BJ et al. Determining ideal body weight.
#'   (Letter). Am J Hosp Pharm 1983; 40: 1622.)
ideal_weight_Miller <- function(heightm, male, ...)
  ideal_weight_linear(heightm, male, 60, 56.2, 53.1, 1.41, 1.36, ...)

#' @title ideal weight by Broca
#' @description Calculate ideal weight based on Broca (1871) Height in cm -100
#'   for women, -105 for men Broca PP. Memoires d'anthropologie. Paris 1871 /
#'   1877.
#' @rdname ideal_weight
#' @export
ideal_weight_Broca <- function(heightm, male, ...)
  ideal_weight_linear(heightm, male, 0, -100, -105, 2.54, 2.54, ...)

#' @title ideal weight by Lemmens
#' @description Lemmens merhod assumes BMI 22 as ideal (Obesity Surgery 2005)
#' TODO: verbose height bounds check
#' @rdname ideal_weight
#' @export
ideal_weight_Lemmens <- function(heightm)
  22 * heightm ^ 2

#' @title ideal weight by gender, offset and gradient
#' @description generic internal function to handle linear ideal weight
#'   calculations. Unofrtunately mixes inches and meters at present.
#' @param heightmininch, height above which to start scaling; consider as a
#'   minimum (for at least some algorithms)
#' @param male_min_kg, y intercept i.e. ideal weight at minimum height for males
#' @param female_min_kg, y intercept i.e. ideal weight at minimum height for
#'   females
#' @param male_kg_per_inch, slope for males
#' @param female_kg_per_inch, slope for females
#' @template warn
#' @rdname ideal_weight
#' @keywords internal
ideal_weight_linear <- function(heightm, male,
                                heightmininch,
                                male_min_kg, female_min_kg,
                                male_kg_per_inch, female_kg_per_inch,
                                warn = FALSE) {
  #NA height or NA maleness are both allowed, and should give NA.

  if (length(heightm) != length(male))
    stop("height (%d) and sex (%d) vectors are different lengths",
         length(heightm), length(male))

  heightinch <- heightm * 100 / 2.54

  f2mintercept <- male_min_kg - female_min_kg
  f2mgradient <- male_kg_per_inch - female_kg_per_inch

  # TODO: vectorize errors and result!
  if (any(heightinch < 0.75 * heightmininch, na.rm = TRUE))
    if (warn) warning(sprintf(
      "calculating ideal weight based on some very low height of %.2fm inches",
      heightm[which(heightinch < 0.75 * heightmininch)]))

  if (any(heightinch < heightmininch, na.rm = TRUE))
    if (warn) warning(sprintf(
      "calculating ideal Weight based on some low height of %.3fm inches",
      heightm[which(heightinch < heightmininch)]))
  if (any(heightinch > 9 * 12, na.rm = TRUE))
    if (warn) warning(
      "calculating ideal_weight_ based on some big heights of %.3fm inches",
      heightm[which(heightinch > 9 * 12)])
  if (any(heightinch > 8 * 12, na.rm = TRUE))
    if (warn) warning(
      "calculating ideal_weight_ based on some big heights of %.3fm inches",
      heightm[which(heightinch > 8 * 12)])

  female_min_kg + f2mintercept * male +
    (heightinch - heightmininch) * (female_kg_per_inch + f2mgradient * male)
}

#' @title Estimate Blood Volume
#' @rdname bloodvol
#' @description estimate blood volume according to the classic 1960s paper by
#'   Nadler. Surgery. 1962 Feb;51(2):224-32. Prediction of blood volume in
#'   normal human adults. Nadler SB, Hidalgo JH, Bloch T.
#' @inheritParams ideal_weight_adult
#' @template weightkg
#' @template warn
#' @examples
#' blood_vol_Nadler(1.8, 80, male = TRUE)
#' blood_vol_Nadler(1.8, 160, male = TRUE)
#' blood_vol_Nadler(1.8, 80, male = FALSE)
#' @export
blood_vol_Nadler <- function(heightm, weightkg, male,
                             warn = FALSE) {

  if (!is.numeric(heightm)) stop("need numeric height input")
  if (any(is.na(heightm)) & warn) warning("need non-NA height input")
  if (!is.numeric(weightkg)) stop("need numeric weight input")
  if (any(is.na(weightkg)) & warn) warning("need non-NA weight input")

  if (length(heightm) != length(weightkg) |
        length(male) != length(heightm)) {
    stop("length(heightm)=%d", length(heightm))
    stop("length(weightkg)=%d", length(weightkg))
    stop("length(male)=%d", length(male))
    stop("NadlerBloodVolume requires that the height weight and \
         male vectors are all the same length.")
  }
  if (warn) {
    if (any(heightm <  0.1, na.rm = TRUE))
      warning("some heights are less than a 10cm!")
    if (any(heightm >  3,   na.rm = TRUE))
      warning("some heights are greater than 3m")
    if (any(weightkg < 0.1, na.rm = TRUE))
      warning("some weights are less than 100g")
    if (any(weightkg > 400, na.rm = TRUE))
      warning("some weights are greater than 400kg")
  }

  nadler <- (0.3669 - (0.3669 - 0.3561) * !male) * heightm ^ 3 +
    (0.03219 - (0.03219 - 0.03308) * !male) * weightkg +
    (0.6041 - (0.6041 - 0.1833) * !male)
}

#' @title Blood volume by Lemmens et al, 2006
#' @rdname bloodvol
#' @description This effectively reverses engineers an ideal weight from BMI of
#'   22, then use the sqaure root of its ratio to actual body weight to adjust
#'   the 70ml/kg of an ideal weight person. Age-dependent regression equations
#'   for indexed blood volume (InBV) at ideal body weight. (No adjustment made
#'   in obesity by Lemmens.) InBV = 90-0.4 X age (males) InBV = 85-0.4 X age
#'   (females). Sounds like he is saying either they are slim and old or younger
#'   and obese. he doesn't attempt to integrate the formulae.
#'
#'   TODO: include age as cut-off butween the use of differing formulae.
#' @return numeric vector
#' @examples
#' blood_vol_lemmens_sedentary(1.8, 80)
#' blood_vol_lemmens_sedentary(1.8, 160)
#' @export
blood_vol_lemmens_sedentary <- function(heightm, weightkg)
  weightkg * blood_vol_lemmens_indexed(heightm, weightkg)

#' @rdname bloodvol
#' @examples
#' blood_vol_lemmens_indexed(1.8, 80)
#' blood_vol_lemmens_indexed(1.8, 160)
#' @export
blood_vol_lemmens_indexed <- function(heightm, weightkg) {
  stopifnot(length(heightm) == length(weightkg))
  70 / sqrt(weightkg / (22 * heightm ^ 2))
}

#' @title Blood volume for non-obese adult, per Lemmens
#' @rdname bloodvol
#' @description applies to slim adults, but note that the age-related decline is
#'   not seen if high degree of physical activity is maintained. TODO: check BMI
#'   not elevated
#' @details Davy KP, Seals DR. Total blood volume in healthy young and older
#'   men. J Appl Physiol 1994; 76: 2059-62.
#'
#'   Parker-Jones P, Davy KP, DeSouza CA et al. Absence of agerelated decline in
#'   total blood volume in physically active females. Am J Physiol 1997; 272:
#'   H2534-40.
#' @param age years
#' @param male logical
#' @examples
#'   blood_vol_lemmens_non_obese(80, age = 25, male = TRUE)
#'   blood_vol_lemmens_non_obese(80, age = 75, male = TRUE)
#' @export
blood_vol_lemmens_non_obese <- function(weightkg, age, male)
  ifelse(male,
         weightkg * (90 - (0.4 * age)),
         weightkg * (85 - (0.4 * age))
  )

#' @title adjusted body weight
#' @description returns ideal weight + 40% of difference between ideal and
#'   actual weights. Ideal weight is calculated using default algorithm. TODO:
#'   is downward adjustment valid?
#' @inheritParams ideal_weight_adult
#' @param weightkg weight in kg, may be a vector
#' @examples
#' adj_weight_adult(1.6, 120, male = FALSE)
#' @export
adj_weight_adult <- function(heightm, weightkg, male) {
  stopifnot(length(heightm) == length(weightkg))
  stopifnot(length(male) == length(weightkg))
  #TODO: is downward adjustment valid?
  0.6 * ideal_weight_adult(heightm, male) + 0.4 * weightkg
}

#' @title Body Mass Index (BMI) for adults
#' @rdname bmi
#' @template heightm
#' @template weightkg
#' @examples
#' bmi_adult(1.6, 120)
#' bmi_adult(2, 75)
#' @export
bmi_adult <- function(heightm, weightkg)
  weightkg / (heightm ^ 2)

#' @rdname bmi
#' @param heightin height in inches
#' @param weightlb weight in pounds
#' @examples
#' bmi_adult_ins_lbs(72, 200)
#' @export
bmi_adult_ins_lbs <- function(heightin, weightlb)
  bmi_adult(heightin * 0.0254, weightlb * 2.20462)
