#' @include imports.R
NULL

# European Scripts
#' @rdname Unicode
#' @export
ARMENIAN <- as.regex("\u0530-\u058F")

#' @rdname Unicode
#' @export
ARMENIAN_LIGATURES <- as.regex("\uFB00-\uFB4F")

#' @rdname Unicode
#' @export
CAUCASIAN_ALBANIAN <- as.regex("\U10530-\U1056F")

#' @rdname Unicode
#' @export
CYPRIOT_SYLLABARY <- as.regex("\U10800-\U1083F")

#' @rdname Unicode
#' @export
CYRILLIC <- as.regex("\u0400-\u04FF")

#' @rdname Unicode
#' @export
CYRILLIC_SUPPLEMENT <- as.regex("\u0500-\u052F")

#' @rdname Unicode
#' @export
CYRILLIC_EXTENDED_A <- as.regex("\u2DE0-\u2DFF")

#' @rdname Unicode
#' @export
CYRILLIC_EXTENDED_B <- as.regex("\uA640-\uA69F")

#' @rdname Unicode
#' @export
ELBASAN <- as.regex("\U10500-\U1052F")

#' @rdname Unicode
#' @export
GEORGIAN <- as.regex("\u10A0-\u10FF")

#' @rdname Unicode
#' @export
GEORGIAN_SUPPLEMENT <- as.regex("\u2D00-\u2D2F")

#' @rdname Unicode
#' @export
GLAGOLITIC <- as.regex("\u2C00-\u2C5F")

#' @rdname Unicode
#' @export
GOTHIC <- as.regex("\U10330-\U1034F")

#' @rdname Unicode
#' @export
GREEK_AND_COPTIC <- as.regex("\u0370-\u03FF")

#' @rdname Unicode
#' @export
GREEK_EXTENDED <- as.regex("\u1F00-\u1FFF")

#' @rdname Unicode
#' @export
LATIN <- as.regex("\u0020-\u007F") #This includes punctuation and space, not control chars.

#' @rdname Unicode
#' @export
LATIN_1_SUPPLEMENT <- as.regex("\u0080-\u00FF")

#' @rdname Unicode
#' @export
LATIN_EXTENDED_A <- as.regex("\u0100-\u017F")

#' @rdname Unicode
#' @export
LATIN_EXTENDED_B <- as.regex("\u0180-\u024F")

#' @rdname Unicode
#' @export
LATIN_EXTENDED_C <- as.regex("\u2C60-\u2C7F")

#' @rdname Unicode
#' @export
LATIN_EXTENDED_D <- as.regex("\uA720-\uA7FF")

#' @rdname Unicode
#' @export
LATIN_EXTENDED_E <- as.regex("\uAB30-\uAB6F")

#' @rdname Unicode
#' @export
LATIN_EXTENDED_ADDITIONAL <- as.regex("\u1E00-\u1EFF")

#' @rdname Unicode
#' @export
LATIN_LIGATURES <- as.regex("\uFB00-\uFB4F")

#' @rdname Unicode
#' @export
LINEAR_A <- as.regex("\U10600-\U1077F")

#' @rdname Unicode
#' @export
LINEAR_B_SYLLABARY <- as.regex("\U10000-\U1007F")

#' @rdname Unicode
#' @export
LINEAR_B_IDEOGRAMS <- as.regex("\U10080-\U100FF")

#' @rdname Unicode
#' @export
OGHAM <- as.regex("\u1680-\u169F")

#' @rdname Unicode
#' @export
OLD_ITALIC <- as.regex("\U10300-\U1032F")

#' @rdname Unicode
#' @export
OLD_PERMIC <- as.regex("\U10350-\U1037F")

#' @rdname Unicode
#' @export
PHAISTOS_DISC <- as.regex("\U101D0-\U101FF")

#' @rdname Unicode
#' @export
RUNIC <- as.regex("\u16A0-\u16FF")

#' @rdname Unicode
#' @export
SHAVIAN <- as.regex("\U10450-\U1047F")


# Phonetic & Shorthand Symbols
#' @rdname Unicode
#' @export
DUPLOYAN <- as.regex("\U1BC00-\U1BC9F")

#' @rdname Unicode
#' @export
SHORTHAND_FORMAT_CONTROLS <- as.regex("\U1BCA0-\U1BCAF")

#' @rdname Unicode
#' @export
IPA_EXTENSIONS <- as.regex("\u0250-\u02AF")

#' @rdname Unicode
#' @export
PHONETIC_EXTENSIONS <- as.regex("\u1D00-\u1D7F")

#' @rdname Unicode
#' @export
PHONETIC_EXTENSIONS_SUPPLEMENT <- as.regex("\u1D80-\u1DBF")

#' @rdname Unicode
#' @export
MODIFIER_TONE_LETTERS <- as.regex("\uA700-\uA71F")

#' @rdname Unicode
#' @export
SPACING_MODIFIER_LETTERS <- as.regex("\u02B0-\u02FF")

#' @rdname Unicode
#' @export
SUPERSCRIPTS_AND_SUBSCRIPTS <- as.regex("\u2070-\u209F")


# Combining Diacritics
#' @rdname Unicode
#' @export
COMBINING_DIACRITIC_MARKS <- as.regex("\u0300-\u036F")

#' @rdname Unicode
#' @export
COMBINING_DIACRITIC_SUPPLEMENT <- as.regex("\u1DC0-\u1DFF")

#' @rdname Unicode
#' @export
COMBINING_DIACRITIC_EXTENDED <- as.regex("\u1AB0-\u1AFF")

#' @rdname Unicode
#' @export
COMBINING_HALF_MARKS <- as.regex("\uFE20-\uFE2F")


#African Scripts
#' @rdname Unicode
#' @export
BAMUN <- as.regex("\uA6A0-\uA6FF")

#' @rdname Unicode
#' @export
BAMUN_SUPPLEMENT <- as.regex("\U16800-\U16A3F")

#' @rdname Unicode
#' @export
BASSA_VAH <- as.regex("\U16AD0-\U16AFF")

#' @rdname Unicode
#' @export
COPTIC <- as.regex("\u2C80-\u2CFF")

#' @rdname Unicode
#' @export
COPTIC_EPACT_NUMBERS <- as.regex("\U102E0-\U102FF")

#' @rdname Unicode
#' @export
EGYPTIAN_HIEROGLYPHS <- as.regex("\U13000-\U1342F")

#' @rdname Unicode
#' @export
ETHIOPIC <- as.regex("\u1200-\u137F")

#' @rdname Unicode
#' @export
ETHIOPIC_SUPPLEMENT <- as.regex("\u1380-\u139F")

#' @rdname Unicode
#' @export
ETHIOPIC_EXTENDED <- as.regex("\u2D80-\u2DDF")

#' @rdname Unicode
#' @export
ETHIOPIC_EXTENDED_A <- as.regex("\uAB00-\uAB2F")

#' @rdname Unicode
#' @export
MENDE_KIKAKUI <- as.regex("\U1E800-\U1E8DF")

#' @rdname Unicode
#' @export
MEROITIC_CURSIVE <- as.regex("\U109A0-\U109FF")

#' @rdname Unicode
#' @export
MEROITIC_HIEROGLYPHS <- as.regex("\U10980-\U1099F")

#' @rdname Unicode
#' @export
NKO <- as.regex("\u07C0-\u07FF")

#' @rdname Unicode
#' @export
OSMANYA <- as.regex("\U10480-\U104AF")

#' @rdname Unicode
#' @export
TIFINAGH <- as.regex("\u2D30-\u2D7F")

#' @rdname Unicode
#' @export
VAI <- as.regex("\uA500-\uA63F")


# Middle Eastern Scripts
#' @rdname Unicode
#' @export
ARABIC <- as.regex("\u0600-\u06FF")

#' @rdname Unicode
#' @export
ARABIC_SUPPLEMENT <- as.regex("\u0750-\u077F")

#' @rdname Unicode
#' @export
ARABIC_EXTENDED_A <- as.regex("\u08A0-\u08FF")

#' @rdname Unicode
#' @export
ARABIC_PRESENTATION_FORMS_A <- as.regex("\uFB50-\uFDFF")

#' @rdname Unicode
#' @export
ARABIC_PRESENTATION_FORMS_B <- as.regex("\uFE70-\uFEFF")

#' @rdname Unicode
#' @export
IMPERIAL_ARAMAIC <- as.regex("\U10840-\U1085F")

#' @rdname Unicode
#' @export
AVESTAN <- as.regex("\U10B00-\U10B3F")

#' @rdname Unicode
#' @export
CARIAN <- as.regex("\U102A0-\U102DF")

#' @rdname Unicode
#' @export
CUNEIFORM <- as.regex("\U12000-\U123FF")

#' @rdname Unicode
#' @export
CUNEIFORM_NUMBERS_AND_PUNCTUATION <- as.regex("\U12400-\U1247F")

#' @rdname Unicode
#' @export
OLD_PERSIAN <- as.regex("\U103A0-\U103DF")

#' @rdname Unicode
#' @export
UGARITIC <- as.regex("\U10380-\U1039F")

#' @rdname Unicode
#' @export
HEBREW <- as.regex("\u0590-\u05FF")

#' @rdname Unicode
#' @export
LYCIAN <- as.regex("\U10280-\U1029F")

#' @rdname Unicode
#' @export
LYDIAN <- as.regex("\U10920-\U1093F")

#' @rdname Unicode
#' @export
MANDAIC <- as.regex("\u0840-\u085F")

#' @rdname Unicode
#' @export
NABATAEAN <- as.regex("\U10880-\U108AF")

#' @rdname Unicode
#' @export
OLD_NORTH_ARABIAN <- as.regex("\U10A80-\U10A9F")

#' @rdname Unicode
#' @export
OLD_SOUTH_ARABIAN <- as.regex("\U10A60-\U10A7F")

#' @rdname Unicode
#' @export
PAHLAVI_INSCRIPTIONAL <- as.regex("\U10B60-\U10B7F")

#' @rdname Unicode
#' @export
PAHLAVI_PSALTER <- as.regex("\U10B80-\U10BAF")

#' @rdname Unicode
#' @export
PALMYRENE <- as.regex("\U10860-\U1087F")

#' @rdname Unicode
#' @export
PHOENICIAN <- as.regex("\U10900-\U1091F")

#' @rdname Unicode
#' @export
SAMARITAN <- as.regex("\u0800-\u083F")

#' @rdname Unicode
#' @export
SYRIAC <- as.regex("\u0700-\u074F")


# Central Asian Scripts
#' @rdname Unicode
#' @export
MANICHAEAN <- as.regex("\U10AC0-\U10AFF")

#' @rdname Unicode
#' @export
MONGOLIAN <- as.regex("\u1800-\u18AF")

#' @rdname Unicode
#' @export
OLD_TURKIC <- as.regex("\U10C00-\U10C4F")

#' @rdname Unicode
#' @export
PHAGS_PA <- as.regex("\uA840-\uA87F")

#' @rdname Unicode
#' @export
TIBETAN <- as.regex("\u0F00-\u0FFF")


# South Asian Scripts
#' @rdname Unicode
#' @export
BENGALI_AND_ASSAMESE <- as.regex("\u0980-\u09FF")

#' @rdname Unicode
#' @export
BRAHMI <- as.regex("\U11000-\U1107F")

#' @rdname Unicode
#' @export
CHAKMA <- as.regex("\U11100-\U1114F")

#' @rdname Unicode
#' @export
DEVANAGARI <- as.regex("\u0900-\u097F")

#' @rdname Unicode
#' @export
DEVANAGARI_EXTENDED <- as.regex("\uA8E0-\uA8FF")

#' @rdname Unicode
#' @export
GRANTHA <- as.regex("\U11300-\U1137F")

#' @rdname Unicode
#' @export
GUJARATI <- as.regex("\u0A80-\u0AFF")

#' @rdname Unicode
#' @export
GURMUKHI <- as.regex("\u0A00-\u0A7F")

#' @rdname Unicode
#' @export
KAITHI <- as.regex("\U11080-\U110CF")

#' @rdname Unicode
#' @export
KANNADA <- as.regex("\u0C80-\u0CFF")

#' @rdname Unicode
#' @export
KHAROSHTHI <- as.regex("\U10A00-\U10A5F")

#' @rdname Unicode
#' @export
KHOJKI <- as.regex("\U11200-\U1124F")

#' @rdname Unicode
#' @export
KHUDAWADI <- as.regex("\U112B0-\U112FF")

#' @rdname Unicode
#' @export
LEPCHA <- as.regex("\u1C00-\u1C4F")

#' @rdname Unicode
#' @export
LIMBU <- as.regex("\u1900-\u194F")

#' @rdname Unicode
#' @export
MAHAJANI <- as.regex("\U11150-\U1117F")

#' @rdname Unicode
#' @export
MALAYALAM <- as.regex("\u0D00-\u0D7F")

#' @rdname Unicode
#' @export
MEETEI_MAYEK <- as.regex("\uABC0-\uABFF")

#' @rdname Unicode
#' @export
MEETEI_MAYEK_EXTENSIONS <- as.regex("\uAAE0-\uAAFF")

#' @rdname Unicode
#' @export
MODI <- as.regex("\U11600-\U1165F")

#' @rdname Unicode
#' @export
MRO <- as.regex("\U16A40-\U16A6F")

#' @rdname Unicode
#' @export
OL_CHIKI <- as.regex("\u1C50-\u1C7F")

#' @rdname Unicode
#' @export
ORIYA <- as.regex("\u0B00-\u0B7F")

#' @rdname Unicode
#' @export
SAURASHTRA <- as.regex("\uA880-\uA8DF")

#' @rdname Unicode
#' @export
SHARADA <- as.regex("\U11180-\U111DF")

#' @rdname Unicode
#' @export
SIDDHAM <- as.regex("\U11580-\U115FF")

#' @rdname Unicode
#' @export
SINHALA <- as.regex("\u0D80-\u0DFF")

#' @rdname Unicode
#' @export
SINHALA_ARCHAIC_NUMBERS <- as.regex("\U111E0-\U111FF")

#' @rdname Unicode
#' @export
SORA_SOMPENG <- as.regex("\U110D0-\U110FF")

#' @rdname Unicode
#' @export
SYLOTI_NAGRI <- as.regex("\uA800-\uA82F")

#' @rdname Unicode
#' @export
TAKRI <- as.regex("\U11680-\U116CF")

#' @rdname Unicode
#' @export
TAMIL <- as.regex("\u0B80-\u0BFF")

#' @rdname Unicode
#' @export
TELUGU <- as.regex("\u0C00-\u0C7F")

#' @rdname Unicode
#' @export
THAANA <- as.regex("\u0780-\u07BF")

#' @rdname Unicode
#' @export
TIRHUTA <- as.regex("\U11480-\U114DF")

#' @rdname Unicode
#' @export
VEDIC_EXTENSIONS <- as.regex("\u1CD0-\u1CFF")

#' @rdname Unicode
#' @export
WARANG_CITI <- as.regex("\U118A0-\U118FF")


# Southeast Asian Scripts
#' @rdname Unicode
#' @export
CHAM <- as.regex("\uAA00-\uAA5F")

#' @rdname Unicode
#' @export
KAYAH_LI <- as.regex("\uA900-\uA92F")

#' @rdname Unicode
#' @export
KHMER <- as.regex("\u1780-\u17FF")

#' @rdname Unicode
#' @export
KHMER_SYMBOLS <- as.regex("\u19E0-\u19FF")

#' @rdname Unicode
#' @export
LAO <- as.regex("\u0E80-\u0EFF")

#' @rdname Unicode
#' @export
MYANMAR <- as.regex("\u1000-\u109F")

#' @rdname Unicode
#' @export
MYANMAR_EXTENDED_A <- as.regex("\uAA60-\uAA7F")

#' @rdname Unicode
#' @export
MYANMAR_EXTENDED_B <- as.regex("\uA9E0-\uA9FF")

#' @rdname Unicode
#' @export
NEW_TAI_LUE <- as.regex("\u1980-\u19DF")

#' @rdname Unicode
#' @export
PAHAWH_HMONG <- as.regex("\U16B00-\U16B8F")

#' @rdname Unicode
#' @export
PAU_CIN_HAU <- as.regex("\U11AC0-\U11AFF")

#' @rdname Unicode
#' @export
TAI_LE <- as.regex("\u1950-\u197F")

#' @rdname Unicode
#' @export
TAI_THAM <- as.regex("\u1A20-\u1AAF")

#' @rdname Unicode
#' @export
TAI_VIET <- as.regex("\uAA80-\uAADF")

#' @rdname Unicode
#' @export
THAI <- as.regex("\u0E00-\u0E7F")


# Indonesia & Oceania Scripts
#' @rdname Unicode
#' @export
BALINESE <- as.regex("\u1B00-\u1B7F")

#' @rdname Unicode
#' @export
BATAK <- as.regex("\u1BC0-\u1BFF")

#' @rdname Unicode
#' @export
BUGINESE <- as.regex("\u1A00-\u1A1F")

#' @rdname Unicode
#' @export
BUHID <- as.regex("\u1740-\u175F")

#' @rdname Unicode
#' @export
HANUNOO <- as.regex("\u1720-\u173F")

#' @rdname Unicode
#' @export
JAVANESE <- as.regex("\uA980-\uA9DF")

#' @rdname Unicode
#' @export
REJANG <- as.regex("\uA930-\uA95F")

#' @rdname Unicode
#' @export
SUNDANESE <- as.regex("\u1B80-\u1BBF")

#' @rdname Unicode
#' @export
SUNDANESE_SUPPLEMENT <- as.regex("\u1CC0-\u1CCF")

#' @rdname Unicode
#' @export
TAGALOG <- as.regex("\u1700-\u171F")

#' @rdname Unicode
#' @export
TAGBANWA <- as.regex("\u1760-\u177F")


# East Asian Scripts
#' @rdname Unicode
#' @export
BOPOMOFO <- as.regex("\u3100-\u312F")

#' @rdname Unicode
#' @export
BOPOMOFO_EXTENDED <- as.regex("\u31A0-\u31BF")

#' @rdname Unicode
#' @export
CJK_UNIFIED_IDEOGRAPHS <- as.regex("\u4E00-\u9FCC")

#' @rdname Unicode
#' @export
CJK_UNIFIED_IDEOGRAPHS_EXTENSION_A <- as.regex("\u3400-\u4DB5")

#' @rdname Unicode
#' @export
CJK_UNIFIED_IDEOGRAPHS_EXTENSION_B <- as.regex("\U20000-\U2A6D6")

#' @rdname Unicode
#' @export
CJK_UNIFIED_IDEOGRAPHS_EXTENSION_C <- as.regex("\U2A700-\U2B734")

#' @rdname Unicode
#' @export
CJK_UNIFIED_IDEOGRAPHS_EXTENSION_D <- as.regex("\U2B740-\U2B81D")

#' @rdname Unicode
#' @export
CJK_COMPATIBILITY_IDEOGRAPHS <- as.regex("\uF900-\uFAFF")

#' @rdname Unicode
#' @export
CJK_COMPATIBILITY_IDEOGRAPHS_SUPPLEMENT <- as.regex("\U2F800-\U2FA1F")

#' @rdname Unicode
#' @export
KANGXI_RADICALS <- as.regex("\u2F00-\u2FDF")

#' @rdname Unicode
#' @export
KANGXI_RADICALS_SUPPLEMENT <- as.regex("\u2E80-\u2EFF")

#' @rdname Unicode
#' @export
CJK_STROKES <- as.regex("\u31C0-\u31EF")

#' @rdname Unicode
#' @export
CJK_IDEOGRAPHIC_DESCRIPTION_CHARACTERS <- as.regex("\u2FF0-\u2FFF")

#' @rdname Unicode
#' @export
HANGUL_JAMO <- as.regex("\u1100-\u11FF")

#' @rdname Unicode
#' @export
HANGUL_JAMO_EXTENDED_A <- as.regex("\uA960-\uA97F")

#' @rdname Unicode
#' @export
HANGUL_JAMO_EXTENDED_B <- as.regex("\uD7B0-\uD7FF")

#' @rdname Unicode
#' @export
HANGUL_COMPATIBILITY_JAMO <- as.regex("\u3130-\u318F")

#' @rdname Unicode
#' @export
HANGUL_SYLLABLES <- as.regex("\uAC00-\uD7AF")

#' @rdname Unicode
#' @export
HIRAGANA <- as.regex("\u3040-\u309F")

#' @rdname Unicode
#' @export
KATAKANA <- as.regex("\u30A0-\u30FF")

#' @rdname Unicode
#' @export
KATAKANA_PHONETIC_EXTENSIONS <- as.regex("\u31F0-\u31FF")

#' @rdname Unicode
#' @export
KANA_SUPPLEMENT <- as.regex("\U1B000-\U1B0FF")

#' @rdname Unicode
#' @export
KANBUN <- as.regex("\u3190-\u319F")

#' @rdname Unicode
#' @export
LISU <- as.regex("\uA4D0-\uA4FF")

#' @rdname Unicode
#' @export
MIAO <- as.regex("\U16F00-\U16F9F")

#' @rdname Unicode
#' @export
YI_SYLLABLES <- as.regex("\uA000-\uA48F")

#' @rdname Unicode
#' @export
YI_RADICALS <- as.regex("\uA490-\uA4CF")


# American Scripts
#' @rdname Unicode
#' @export
CHEROKEE <- as.regex("\u13A0-\u13FF")

#' @rdname Unicode
#' @export
DESERET <- as.regex("\U10400-\U1044F")

#' @rdname Unicode
#' @export
UNIFIED_CANADIAN_ABORIGINAL_SYLLABICS <- as.regex("\u1400-\u167F")

#' @rdname Unicode
#' @export
UNIFIED_CANADIAN_ABORIGINAL_SYLLABICS_EXTENDED <- as.regex("\u18B0-\u18FF")


# Other
#' @rdname Unicode
#' @export
ALPHABETIC_PRESENTATION_FORMS <- as.regex("\uFB00-\uFB4F")

#' @rdname Unicode
#' @export
HALFWIDTH_AND_FULLWIDTH_FORMS <- as.regex("\uFF00-\uFFEF")


# Punctuation
#' @rdname Unicode
#' @export
GENERAL_PUNCTUATION <- as.regex("\u2000-\u206F")

#' @rdname Unicode
#' @export
LATIN_1_PUNCTUATION <- as.regex("\u00A1-00BF")

#' @rdname Unicode
#' @export
SMALL_FORM_VARIANTS <- as.regex("\uFE50-\uFE6F")

#' @rdname Unicode
#' @export
SUPPLEMENTAL_PUNCTUATION <- as.regex("\u2E00-\u2E7F")

#' @rdname Unicode
#' @export
CJK_SYMBOLS_AND_PUNCTUATION <- as.regex("\u3000-\u303F")

#' @rdname Unicode
#' @export
CJK_COMPATIBILITY_FORMS <- as.regex("\uFE30-\uFE4F")

#' @rdname Unicode
#' @export
FULLWIDTH_ASCII_PUNCTUATION <- as.regex("\uFF01-\uFF60")

#' @rdname Unicode
#' @export
VERTICAL_FORMS <- as.regex("\uFE10-\uFE1F")


# Alphanumeric Symbols
#' @rdname Unicode
#' @export
LETTERLIKE_SYMBOLS <- as.regex("\u2100-\u214F")

#' @rdname Unicode
#' @export
ANCIENT_SYMBOLS <- as.regex("\U10190-\U101CF")

#' @rdname Unicode
#' @export
MATHEMATICAL_ALPHANUMERIC_SYMBOLS <- as.regex("\U1D400-\U1D7FF")

#' @rdname Unicode
#' @export
ARABIC_MATHEMATICAL_ALPHANUMERIC_SYMBOLS <- as.regex("\U1EE00-\U1EEFF")

#' @rdname Unicode
#' @export
ENCLOSED_ALPHANUMERICS <- as.regex("\u2460-\u24FF")

#' @rdname Unicode
#' @export
ENCLOSED_ALPHANUMERIC_SUPPLEMENT <- as.regex("\U1F100-\U1F1FF")

#' @rdname Unicode
#' @export
ENCLOSED_CJK_LETTERS_AND_MONTHS <- as.regex("\u3200-\u32FF")

#' @rdname Unicode
#' @export
ENCLOSED_IDEOGRAPHIC_SUPPLEMENT <- as.regex("\U1F200-\U1F2FF")

#' @rdname Unicode
#' @export
CJK_COMPATIBILITY <- as.regex("\u3300-\u33FF")


# Technical Symbols
#' @rdname Unicode
#' @export
MISCELLANEOUS_TECHNICAL <- as.regex("\u2300-\u23FF")

#' @rdname Unicode
#' @export
CONTROL_PICTURES <- as.regex("\u2400-\u243F")

#' @rdname Unicode
#' @export
OPTICAL_CHARACTER_RECOGNITION <- as.regex("\u2440-\u245F")


# Combining Diacritics
#' @rdname Unicode
#' @export
COMBINING_DIACRITIC_MARKS_FOR_SYMBOLS <- as.regex("\u20D0-\u20FF")


# Numbers and Digits
#' @rdname Unicode
#' @export
AEGEAN_NUMBERS <- as.regex("\U10100-\U1013F")

#' @rdname Unicode
#' @export
ANCIENT_GREEK_NUMBERS <- as.regex("\U10140-\U1018F")

#' @rdname Unicode
#' @export
FULLWIDTH_ASCII_DIGITS <- as.regex("\uFF10-\uFF19")

#' @rdname Unicode
#' @export
COMMON_INDIC_NUMBER_FORMS <- as.regex("\uA830-\uA83F")

#' @rdname Unicode
#' @export
COPTIC_EPACT_NUMBERS <- as.regex("\U102E0-\U102FF")

#' @rdname Unicode
#' @export
COUNTING_ROD_NUMERALS <- as.regex("\U1D360-\U1D37F")

#' @rdname Unicode
#' @export
NUMBER_FORMS <- as.regex("\u2150-\u218F")

#' @rdname Unicode
#' @export
RUMI_NUMERAL_SYMBOLS <- as.regex("\U10E60-\U10E7F")

#' @rdname Unicode
#' @export
SINHALA_ARCHAIC_NUMBERS <- as.regex("\U111E0-\U111FF")


# Mathematical Symbols
#' @rdname Unicode
#' @export
MATH_ARROWS <- as.regex("\u2190-\u21FF")

#' @rdname Unicode
#' @export
SUPPLEMENTAL_ARROWS_A <- as.regex("\u27F0-\u27FF")

#' @rdname Unicode
#' @export
SUPPLEMENTAL_ARROWS_A <- as.regex("\u2900-\u297F")

#' @rdname Unicode
#' @export
SUPPLEMENTAL_ARROWS_A <- as.regex("\U1F800-\U1F8FF")

#' @rdname Unicode
#' @export
ADDITIONAL_ARROWS <- as.regex("\u2B00-\u2BFF")


# Mathematical Operators
#' @rdname Unicode
#' @export
SUPPLEMENTAL_MATHEMATICAL_OPERATORS <- as.regex("\u2A00-\u2AFF")

#' @rdname Unicode
#' @export
MISCELLANEOUS_MATHEMATICAL_SYMBOLS_A <- as.regex("\u27C0-\u27EF")

#' @rdname Unicode
#' @export
MISCELLANEOUS_MATHEMATICAL_SYMBOLS_B <- as.regex("\u2980-\u29FF")

#' @rdname Unicode
#' @export
FLOORS_AND_CEILINGS <- as.regex("\u2308-\u230B")

#' @rdname Unicode
#' @export
INVISIBLE_OPERATORS <- as.regex("\u2061-\u2064")

#' @rdname Unicode
#' @export
GEOMETRIC_SHAPES <- as.regex("\u25A0-\u25FF")

#' @rdname Unicode
#' @export
BOX_DRAWING <- as.regex("\u2500-\u257F")

#' @rdname Unicode
#' @export
BLOCK_ELEMENTS <- as.regex("\u2580-\u259F")

#' @rdname Unicode
#' @export
GEOMETRIC_SHAPES_EXTENDED <- as.regex("\U1F780-\U1F7FF")


# Other Symbols
#' @rdname Unicode
#' @export
ALCHEMICAL_SYMBOLS <- as.regex("\U1F700-\U1F77F")

#' @rdname Unicode
#' @export
BRAILLE_PATTERNS <- as.regex("\u2800-\u28FF")

#' @rdname Unicode
#' @export
CURRENCY_SYMBOLS <- as.regex("\u20A0-\u20CF")


# Dingbats
#' @rdname Unicode
#' @export
DINGBATS <- as.regex("\u2700-\u27BF")

#' @rdname Unicode
#' @export
ORNAMENTAL_DINGBATS <- as.regex("\U1F650-\u1F67")


# Emoticons
#' @rdname Unicode
#' @export
EMOTICONS <- as.regex("\U1F600-\U1F64F")


# Game Symbols
#' @rdname Unicode
#' @export
CHESS_CHECKERS_DRAUGHTS <- as.regex("\u2654-\u265F\u26C0-\u26C3")

#' @rdname Unicode
#' @export
DOMINO_TILES <- as.regex("\U1F030-\U1F09F")

#' @rdname Unicode
#' @export
JAPANESE_CHESS <- as.regex("\u2616-\u2617")

#' @rdname Unicode
#' @export
MAHJONG_TILES <- as.regex("\U1F000-\U1F02F")

#' @rdname Unicode
#' @export
PLAYING_CARDS <- as.regex("\U1F0A0-\U1F0FF")

#' @rdname Unicode
#' @export
CARD_SUITS <- as.regex("\u2660-\u2667")


# Miscellanous Symbols
#' @rdname Unicode
#' @export
MISCELLANEOUS_SYMBOLS_AND_PICTOGRAPHS <- as.regex("\U1F300-\U1F5FF")


# Musical Symbols
#' @rdname Unicode
#' @export
MUSICAL_SYMBOLS <- as.regex("\U1D100-\U1D1FF")

#' @rdname Unicode
#' @export
ANCIENT_GREEK_MUSICAL_NOTATION <- as.regex("\U1D200-\U1D24F")

#' @rdname Unicode
#' @export
BYZANTINE_MUSICAL_SYMBOLS <- as.regex("\U1D000-\U1D0FF")


# Transport and Map Symbols
#' @rdname Unicode
#' @export
TRANSPORT_AND_MAP_SYMBOLS <- as.regex("\U1F680-\U1F6FF")


# Yijing Symbols
#' @rdname Unicode
#' @export
YIJING_MONO_DI_AND_TRIGRAMS <- as.regex("\u268A-\u268F\u2630-\u2637")

#' @rdname Unicode
#' @export
YIJING_HEXAGRAM_SYMBOLS <- as.regex("\u4DC0-\u4DFF")

#' @rdname Unicode
#' @export
TAI_XUAN_JING_SYMBOLS <- as.regex("\U1D300-\U1D35F")


# Specials
#' @rdname Unicode
#' @export
SPECIALS <- as.regex("\uFFF0-\uFFFD")  #roxygen chokes on the non-characters \uFFFE and \uFFFF


# Tags
#' @rdname Unicode
#' @export
TAGS <- as.regex("\UE0000-\UE007F")


# Variation Selectors
#' @rdname Unicode
#' @export
VARIATION_SELECTORS <- as.regex("\uFE00-\uFE0F")

#' @rdname Unicode
#' @export
VARIATION_SELECTORS_SUPPLEMENT <- as.regex("\UE0100-\UE01EF")


# Private Use
#' @rdname Unicode
#' @export
PRIVATE_USE_AREA <- as.regex("\uE000-\uF8FF")

#' @rdname Unicode
#' @export
SUPPLEMENTARY_PRIVATE_USE_AREA_A <- as.regex("\UF0000-\UFFFFF")

#' @rdname Unicode
#' @export
SUPPLEMENTARY_PRIVATE_USE_AREA_B <- as.regex("\U100000-\U10FFFD")
