apply.pca <-
function (matrix.signal) 
{
    if (dim(matrix.signal)[2] <= 2) {
        m <- apply(matrix.signal, MARGIN = 1, FUN = mean)
        return(m)
    }
    pca <- prcomp(matrix.signal, scale = TRUE)$x[, 1]
    return(pca/sd(pca))
}
