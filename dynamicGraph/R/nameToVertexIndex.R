"nameToVertexIndex" <-
function (vertexnames, vertices) 
{
    return(match(vertexnames, Names(vertices)))
}
