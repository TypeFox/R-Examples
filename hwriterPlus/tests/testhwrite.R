require(hwriterPlus)
### Text in a heading tag
hwrite("Some text", heading = 1)

### Text in a div tag
hwrite("Some text", div = TRUE)

### A link to an external website
hwrite("Auckland University", link = "http://www.auckland.ac.nz")

### A heading with a name so that linking to it is possible
hwrite(hwrite("A Level 2 Heading", name = 'a named heading'),
       heading = 2, center = FALSE, br = TRUE)
