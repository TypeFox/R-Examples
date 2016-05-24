### confusionMatrix.R  (2014-11-13)
###
###    Compute Confusion Matrix 
###
### Copyright 2014  Korbinian Strimmer
###
###
### This file is part of the `crossval' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


confusionMatrix = function(actual, predicted, negative="control") 
{
  idx.null = (actual==negative)
  idx.alternative = !idx.null
  
  # True Positives ("hit")
  # example: predicted is cancer, and actual is cancer
  TP = sum( predicted[idx.alternative] != negative ) 

  # True Negatives ("correct rejection")
  # example:predict is control, and actual is control
  TN = sum( predicted[idx.null] == negative )

  # False Positives ("false alarm", Type I error)
  # example: predicted is cancer, and actual is control  
  FP = sum( predicted[idx.null] != negative ) 

  # False Negatives ("miss", Type II error)
  # example: predicted is control, and actual is cancer  
  FN = sum( predicted[idx.alternative] == negative )

  cm = c(FP, TP, TN, FN)  
  names(cm) = c("FP", "TP", "TN", "FN")
  attr(cm, "negative") = negative

  return ( cm )
}

