"gamlist" <-
function(...)
{
	gl <- list(...)
	oldClass(gl) <- c("gamlist", "glmlist")
	gl
}
