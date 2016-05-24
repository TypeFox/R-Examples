"initModellist" <-
function (m) 
{
        m@parnames <- sort(intersect(
	          setdiff(slotNames(theta()), c("prel","drel")), 
	          slotNames(m)))
        m <- initOneModel(m)

}