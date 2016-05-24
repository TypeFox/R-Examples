library("tools")
Rd2txt_options(
  width = 80L,
  minIndent = 2L,
  extraIndent = 2L,
  sectionIndent = 2L,
  sectionExtra = 2L,
  underline_titles = FALSE
)

Rd2txt("NEWS.Rd", out="NEWS")
