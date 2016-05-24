### diagnosticErrors.R  (2014-11-13)
###
###    Diagnostic Errors: Accuracy, Sensitivity, Specificity, PPV, NPV, LOR 
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


diagnosticErrors = function(cm)
{
  FP = cm["FP"]
  TP = cm["TP"]
  TN = cm["TN"]
  FN = cm["FN"]

  significant = FP+TP       
  notSignificant = TN+FN  
  null = FP+TN              # not interesting 
  alternative = TP+FN       # interesting

  acc = (TP+TN)/(null+alternative)  # accuracy, recognition rate
# err = 1-acc                       # error rate, misclassification rate, 0-1 risk

  sens = TP/alternative     # power, recall, TPR
# alphaII = 1-sens          # FNR, miss rate

  spec = TN/null            # TNR
# alphaI = 1-spec           # FPR, false alarm rate

  ppv = TP/significant      # positive predictive value, TDR, precision
# fdr = 1-ppv

  npv = TN/notSignificant   # negative predictive value, TNDR
# fndr = 1-npv

  lor = log(TP)+log(TN)-log(FN)-log(FP) # log odds ratio
                            # = log ( sens/(1-sens)*spec/(1-spec) )
                            # = log ( ppv/(1-ppv)*npv/(1-npv) )


  err = c(acc, sens, spec, ppv, npv, lor)
  names(err) = c("acc", "sens", "spec", "ppv", "npv", "lor")
  attr(err, "negative") = attr(cm, "negative")

  return ( err )
}






