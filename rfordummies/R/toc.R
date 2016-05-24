tableOfContents <- "Introduction 1

Part I: R You Ready? 7

Chapter 1: Introducing R: The Big Picture 9

Chapter 2: Exploring R 15

Chapter 3: The Fundamentals of R 31

Part II: Getting Down to Work in R 43

Chapter 4: Getting Started with Arithmetic 45

Chapter 5: Getting Started with Reading and Writing 71

Chapter 6: Going on a Date with R 93

Chapter 7: Working in More Dimensions 103

Part III: Coding in R 137

Chapter 8: Putting the Fun in Functions 139

Chapter 9: Controlling the Logical Flow 159

Chapter 10: Debugging Your Code 179

Chapter 11: Getting Help 193

Part IV: Making the Data Talk 203

Chapter 12: Getting Data into and out of R 205

Chapter 13: Manipulating and Processing Data 219

Chapter 14: Summarizing Data 253

Chapter 15: Testing Differences and Relations 275

Part V: Working with Graphics 299

Chapter 16: Using Base Graphics 301

Chapter 17: Creating Faceted Graphics with Lattice 317

Chapter 18: Looking At ggplot2 Graphics 333

Part VI: The Part of Tens 347

Chapter 19: Ten Things You Can Do in R That You Would've Done in Microsoft Excel 349

Chapter 20: Ten Tips on Working with Packages 359

Appendix: Installing R and RStudio 365

Index 371"



#' Print table of contents.
#' 
#' @export
#' @examples
#' toc()
toc <- function(){
  cat(tableOfContents)
  invisible(tableOfContents)
}