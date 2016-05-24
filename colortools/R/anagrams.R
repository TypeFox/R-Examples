#'@title anagram colors
#'
#'@description
#'This function returns the anagrams of a given color.  Anagrams are made by
#'rearranging the hexadecimal letter pairs to form another color.
#'
#'@details
#'An anagram is a word that is made by rearranging the letters of another word.
#'In this case, we have color names instead of words, and hexadecimal
#'letter-pairs instead of single letters. For example, the color "#CD1091" is
#'composed by three pairs "CD", "10", and "91". By rearranging these three
#'pairs, we could obtain an anagram like "#CD9110" or "#10CD91"
#'
#'@param color an R color name or a color in hexadecimal notation
#'@return A character vector with the anagram colors of a given color in
#'hexadecimal notation
#'@author Gaston Sanchez
#'@seealso \code{\link{pizza}}
#'@export
#'@examples
#'
#' # anagrams of 'tomato'
#' ana_tom <- anagrams("tomato")
#'
#' # plot tomato anagrams in a color 'pizza' wheel
#' pizza(ana_tom, bg = "gray20")
#'
anagrams <-
function(color)
{
	# convert to HSV
	col_hsv = rgb2hsv(col2rgb((color)))[,1]
	# get hue, saturation, and value
	hue = col_hsv[1]
	sat = col_hsv[2]
	val = col_hsv[3]
	# transparency
	alpha = 1
	if (substr(color, 1, 1) == "#" && nchar(color) == 9)
		alpha = substr(color, 8, 9)
	# make sure it is in HEX code
	hex = hsv(hue, sat, val)
	# break down hex in pairs
	par1 = substr(hex, 2, 3)
	par2 = substr(hex, 4, 5)
	par3 = substr(hex, 6, 7)
	# put pars in a vector
	hp = c(par1, par2, par3)
	# create vector to store results
	anacols = rep("", 27)
	# generate anagrams
	aux = 1
	for (i in 1:3)
	{
		for (j in 1:3)
		{
			for (k in 1:3)
			{
				anacols[aux] = paste("#", hp[i], hp[j], hp[k], sep="")
				aux = aux + 1
			}
		}
	}
	if (alpha != 1)
		anacols = paste(anacols, alpha, sep="")
	# result
	anacols
}
