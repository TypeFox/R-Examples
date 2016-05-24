# --------------------------------------------------------------------------------
##' Aligned Rank Transform for Nonparametric Factorial Analysis 
##'
##' The function computes a separate aligned response variable for each effect of an user-specified model,
##' transform it into a ranking, and applies a separate ANOVA to every resulting ranked aligned response to
##' check the significance of the corresponding effect.
##' @title Aligned Rank Transform procedure
##' @param formula A formula indicating the model to be fitted.
##' @param data A data frame containing the input data. The name of the columns should match the names used in
##' the user-specified \code{formula} of the model that will be fitted.
##' @param perform.aov Optional: whether separate ANOVAs should be run on the Ranked aligned responses or not.
##' In case it should not, only the ranked aligned responses will be returned. Defaults to \code{TRUE}.
##'	@param SS.type A string indicating the type of sums of squares to be used in the ANOVA on the aligned responses.
##'	Must be one of "I", "II", "III". If \code{perform.aov} was set to \code{FALSE}, the value of \code{SS.type} will be ignored.
##'	Please note SS types coincide when the design is balanced (equal number of observations per cell) but differ otherwise. 
##' Refer to Shaw and Mitchell-Olds (1993) or Fox (1997) for further reading and recomentations on how to conduct ANOVA analyses with unbalanced designs.
##'	@param ... Other arguments passed to \link{lm} when computing effect estimates via ordinary least squares for the alignment.
##' @return A tagged list with the following elements:
##'	\itemize{
##'		\item \code{$aligned}: a data frame with the input data and additional columns to the right, containing the aligned 
##'		and the ranked aligned responses for each model effect.
##'		\item \code{$significance}: (only when \code{perform.aov = TRUE}) the ANOVA table that collects every unique meaningful row of
##'		each of the separate ANOVA tables obtained from the ranked aligned responses.
##' }
##' @author Pablo J. Villacorta Iglesias
##' @references Higgins, J. J., Blair, R. C. and Tashtoush, S. (1990). The aligned rank transform procedure. Proceedings of the Conference on Applied Statistics in Agriculture. 
##'	Manhattan, Kansas: Kansas State University, pp. 185-195.
##' @references Higgins, J. J. and Tashtoush, S. (1994). An aligned rank transform test for interaction. Nonlinear World 1 (2), pp. 201-211. 
##' @references Mansouri, H. (1999). Aligned rank transform tests in linear models. Journal of Statistical Planning and Inference 79, pp. 141 - 155.
##' @references Wobbrock, J.O., Findlater, L., Gergle, D. and Higgins, J.J. (2011). The Aligned Rank Transform for nonparametric factorial analyses using only ANOVA procedures. 
##'	Proceedings of the ACM Conference on Human Factors in Computing Systems (CHI '11). New York: ACM Press, pp. 143-146.
##' @references Higgins, J.J. (2003). Introduction to Modern Nonparametric Statistics. Cengage Learning.
##' @references Shaw, R.G. and Mitchell-Olds, T. (1993). Anova for Unbalanced Data: An Overview. Ecology 74, 6, pp. 1638 - 1645.
##' @references Fox, J. (1997). Applied Regression Analysis, Linear Models, and Related Methods. SAGE Publications.
##' @references ARTool R package, for full models only. \url{http://cran.r-project.org/package=ARTool}
##' @seealso \link{lm}
##' @examples
##' # Input data contained in the Higgins1990-Table1.csv file distributed with ARTool
##'	# The data were used in the 1990 paper cited in the References section 
##'	data(higgins1990, package = "ART"); 
##' # Two-factor full factorial model that will be fitted to the data
##' art.results = aligned.rank.transform(Response ~ Row * Column, data = data.higgins1990);
##' print(art.results$aligned, digits = 4);
##' print(art.results$significance);
aligned.rank.transform<-function(formula, data, perform.aov = TRUE, SS.type = c("III", "II", "I"), ...){

  SS.type = match.arg(SS.type);
  
  ## ------- FORMULA PROCESSING ----------
  
  termsObject = terms(formula, keep.order = FALSE, simplify = FALSE, allowDotAsName = FALSE); 
  responseIndexInFormula = attr(termsObject, "response"); # index of the response within the vector "variables"
  factorMatrix = attr(termsObject, "factors");
  variables = as.character(attr(termsObject, "variables"))[-1]; # first element is always the string "list"
  intercept = attr(termsObject, "intercept");
  
  ## -------------------------------------
  
  frameNames = names(data);
  ncols = length(frameNames); # including the response variable
  responseColumn = NULL;
  
  frameNames = gsub(":", "_", frameNames); # strip character ":" from the column names to avoid confusion with interaction names
  
  if(length(factorMatrix) == 0){
    stop("ERROR: the formula does not contain any variable, but only an intercept constant");
  }  

  for(i in 1:length(variables)){
    if(!(variables[[i]] %in% frameNames)){
      stop("ERROR: variable ",variables[[i]]," not found in the frame data");
    }
    if(i == responseIndexInFormula){ # index of the response variable in the formula
      responseColumn = match(variables[[i]], frameNames);
    }
  }
  
  ## Delete rows with NA in the response
  isdata = !is.na(data[[responseColumn]]);
  data = data[isdata,];
  nrows = nrow(data);  
  
  # --------------------------------------------------------------------------
  # First step: compute estimated effects using ordinary least squares with function lm (linear model)
  # and strip those effects from the response variable 
  # --------------------------------------------------------------------------
  matrixEffectNames = colnames(factorMatrix);

  old.contrasts = getOption("contrasts"); # get contrasts settings ...
  options(contrasts=c("contr.sum", "contr.sum"));  # ... modify them temporarily...

  # --------------- LEAST SQUARES FITTING ---------------
  mymodel = lm(formula = formula, x = TRUE, data = data, ...); # get the model matrix x
  coefficients = mymodel$coefficients;
  coefficients[is.na(coefficients)] = 0;
  # -----------------------------------------------------
  
  model.matrix = mymodel$x;  
  options(contrasts = old.contrasts); # ... and restore original contrast settings
    
  chopped.model.effects = strsplit(x = colnames(model.matrix), ":");
  chopped.effectNames = strsplit(x = matrixEffectNames, ":");
  
  aligned_col_names = sapply(X = matrixEffectNames, FUN = function(x) gsub(":","_", paste("Aligned", x, sep="_")));
  ranked_col_names = sapply(X = matrixEffectNames, FUN = function(x) gsub(":","_", paste("Ranks", x, sep="_")));
  
  aligned.matrix = data.frame(matrix(NA, nrow = nrows, ncol = ncol(factorMatrix)));
  ranked.aligned.matrix = data.frame(matrix(NA, nrow = nrows, ncol = ncol(factorMatrix)));
  colnames(aligned.matrix) = aligned_col_names;
  colnames(ranked.aligned.matrix) = ranked_col_names; 

  for(i in 1:length(chopped.effectNames)){
    # Alignment for effect matrixEffectNames[i]
    # appearances[j] = TRUE iff the i-th effect or interaction appears in the j-th column of the model matrix
    appearances = rep(FALSE, length(chopped.model.effects)); # column 1 of the model is the intercept
    subtraction.matrix = model.matrix;

    for(j in 1:length(chopped.model.effects)){ # check all column names of the model matrix
      if(length(chopped.effectNames[[i]]) == length(chopped.model.effects[[j]])){
        # Both have the same number of effects. Check they are the same
        bool.mat = sapply(X = chopped.effectNames[[i]], FUN = grepl, x = chopped.model.effects[[j]]);

        if( ( is.matrix(bool.mat) && sum(colSums(bool.mat) == 0) == 0 ) || # all effects (e.g. "a", "b") appear at least once in the effect combination (e.g. "a1:b2")
            ( !is.matrix(bool.mat) && sum(bool.mat) == length(chopped.effectNames[[i]]) ) ){ 

					# Delete the effect name that matches from the effect combinations vector and make sure the residual are just numbers
					# that encode the levels (e.g. delete "Row" and "Column" from the combination "Row2:Column1") and check the residuals are "2" and "1".
					replacement.mat = sapply(X = chopped.effectNames[[i]], FUN = sub, x = chopped.model.effects[[j]], replacement = "");
			
					replacement.mat[replacement.mat == ""] = 1; # correction (any number) for residual "" (returned by function sub when there is a perfect match)
					if(is.matrix(replacement.mat)){					
				
            valid = suppressWarnings( !is.na(matrix(as.integer(replacement.mat), nrow = nrow(replacement.mat), ncol=ncol(replacement.mat))) );
            correct = sum(colSums(valid) == 0) == 0;
					}
					else{
            all.numbers = suppressWarnings(as.integer(replacement.mat[bool.mat]));
            correct = sum(is.na(all.numbers)) == 0 # all the conversions were successful -> all residuals are names
					}
					if(correct){
						appearances[j] = TRUE; # found: this column of the model matrix must NOT be stripped out for this alignment
          }
        }
      } 
    }
    subtraction.matrix[,appearances] = 0; 

    # Align the response variable with respect to effect matrixEffectNames[i]    

    aligned.matrix[,i] = data[,responseColumn] - subtraction.matrix %*% coefficients;
    ranked.aligned.matrix[,i] = rank(aligned.matrix[,i], ties.method = "average");
  }
  
  data = data.frame(data, aligned.matrix, ranked.aligned.matrix);
  
  # -------------------------------------------------------------------------------
  # Final step: apply conventional ANOVA to the ranked aligned responses separately
  # -------------------------------------------------------------------------------
  if(perform.aov){

    final.rows = ncol(factorMatrix);
    final.table = data.frame(matrix(nrow=final.rows, ncol = 4)); # Columns of any usual ANOVA table: Sum Sq, Df, F value,  Pr(>F)
    old.contrasts = getOption("contrasts");
    options(contrasts=c("contr.sum", "contr.sum")); 
    
    for(k in 1:length(ranked_col_names)){
      expr = update.formula(formula, as.formula(paste(ranked_col_names[[k]], "~ . ")));
      aux.model = lm(expr, data);    
      if(SS.type == "I" || sum(mymodel$residuals)==0){  # for a perfect fit, we turn to stats::anova() regardless the SS.type
            my.table = anova(aux.model)[,-3];     } # delete the third column (mean sums of squares)
      else{ my.table = Anova(aux.model, type=SS.type, singular.ok = TRUE);    } # Anova() from car package        
      final.table[k,] = my.table[matrixEffectNames[[k]],];
      if(k == 1){ 
        names(final.table) = names(my.table); 
      }
    }
    row.names(final.table) = matrixEffectNames;
    options(contrasts = old.contrasts);
    return(list(aligned = data, significance = final.table));
  }
  else{ # No ANOVA required after the ART 
    return(list(aligned = data));
  }
}