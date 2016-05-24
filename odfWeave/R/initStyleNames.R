initStyleNames <- function(node, env)
{
   sfun <- function(s) xmlGetAttr(s, 'style:name', 'ERROR')
   stylenames <- unlist(treeapply(node, 'style:style', sfun, onlyFirst=FALSE, rooted=FALSE))

   for (sname in stylenames)
   {
      assign(sname, sname, pos=env)
   }
}
