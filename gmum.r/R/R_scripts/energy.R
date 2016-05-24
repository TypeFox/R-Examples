standard_entropy <- function(cluster_points)
{
    dimension <- dim(cluster_points)[2]
    cluster_cov_mat <- cov_mat(cluster_points)
    det_cluster_cov_mat <- det(cluster_cov_mat)
    if(det_cluster_cov_mat == 0)
    {
        det_cluster_cov_mat <- 1.0e-32
    }
    return( dimension/2 * log(2 * pi * exp(1)) + log(det_cluster_cov_mat) / 2 )
}

sphere_entropy <- function(cluster_points)
{
    dimension <- dim(cluster_points)[2]
    cluster_cov_mat_trace <- cov_mat_trace(cluster_points)
    if(cluster_cov_mat_trace == 0)
    {
        cluster_cov_mat_trace <- 1.0e-32
    }
    return ( dimension/2 * log(2 * pi * exp(1) / dimension) + dimension / 2 * log(cluster_cov_mat_trace) )
}

diagonal_entropy <- function(cluster_points)
{
    dimension <- dim(cluster_points)[2]
    cluster_cov_mat <- cov_mat(cluster_points)
    det_cluster_cov_mat <- prod(diag(cluster_cov_mat))
    if(det_cluster_cov_mat == 0)
    {
        det_cluster_cov_mat <- 1.0e-32
    }    
    return ( dimension/2 * log(2 * pi * exp(1)) + log(det_cluster_cov_mat) / 2 )
}

cluster_energy <- function(cluster_entropy, cluster_npoints, npoints)
{
    p <- cluster_npoints / npoints
    return( p * (cluster_entropy - log(p)) )
}

cov_mat <- function(cluster_points)
{
    npoints <- dim(cluster_points)[1]
    dimension <- dim(cluster_points)[2]
    mean <- as.vector(colMeans(cluster_points))
    result <- matrix(nrow = dimension, ncol = dimension, data = 0)
    for(i in 1:npoints)
    {
        p <- as.matrix(cluster_points[i, ] - mean)
        result <- result + (p %*% t(p)) / npoints        
    }
    return(result)    
}

cov_mat_trace <- function(cluster_points)
{
    npoints <- dim(cluster_points)[1]
    mean <- as.vector(colMeans(cluster_points))
    result <- 0.0
    for(i in 1:npoints)
    {
        p <- cluster_points[i, ] - mean
        result <- result + (p %*% p)       
    }
    result <- result / npoints
    return(result)   
}

cec_energy <- function(dataset, clustering, entropy_func)
{
    dimension <- ncol(dataset)
    npoints <- dim(dataset)[1]
    energy <- 0
    for (i in unique(clustering))
    {
        cluster_points <- dataset[clustering == i,] 
        cluster_npoints <- dim(cluster_points)[1]
        curr_cluster_entropy <- entropy_func(cluster_points)
        curr_cluster_energy <- cluster_energy(curr_cluster_entropy, cluster_npoints, npoints)
        energy <- energy + curr_cluster_energy
    }
    return(as.numeric(energy))
}
