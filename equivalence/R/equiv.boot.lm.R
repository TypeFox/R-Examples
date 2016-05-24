# $Id: equiv.boot.lm.R,v 1.4 2005/10/10 10:14:43 andrewr Exp $

"equiv.boot.lm" <-
    function(compare, index, rel.int.int=0.25, rel.int.slope=0.25,
                           b0.ii.absolute = FALSE) {
        observation <- compare[,1]
        prediction <- compare[,2]
	m.pred <- mean(prediction, na.rm=TRUE)
        val.model <- lm(observation[index] ~ I(prediction[index]-m.pred))
        intercept.hat <- coef(val.model)[1]
        intercept <- m.pred
        slope.hat <- coef(val.model)[2]
        inter.up <- ifelse(b0.ii.absolute,
                           intercept + rel.int.int,
                           intercept * (1 + rel.int.int))
        inter.down <- ifelse(b0.ii.absolute,
                             intercept - rel.int.int,
                             intercept * (1 - rel.int.int))
        validate <- c(intercept.hat <= inter.down, 
                      (intercept.hat <= inter.up) & (intercept.hat > 
                         inter.down), intercept.hat > inter.up,
                      slope.hat <= 1 - rel.int.slope, 
                      (slope.hat <= 1 + rel.int.slope) &
                      (slope.hat > 1 - rel.int.slope), 
                      slope.hat > 1 + rel.int.slope, intercept.hat, slope.hat)
      }
