CNV.fitModel <-
function (ncomp, nind, hyp = "H0", data, logit.offset, design.matrix.mean, 
    design.matrix.variance, design.matrix.disease, pi.model = 0, 
    mix.model = 10, control = list(tol = 1e-05, max.iter = 3000, 
        min.freq = 4)) 
{
    pi.mod <- as.integer(pi.model)
    mix.mod <- as.integer(mix.model)
    initial <- sort(unique(data$strats.var))
    data$strats.var <- as.integer(rank(initial)[match(x = data$strats.var, 
        table = initial)])
    initial <- sort(unique(data$strats.mean))
    data$strats.mean <- as.integer(rank(initial)[match(x = data$strats.mean, 
        table = initial)])
    initial <- sort(unique(data$association.strata))
    data$association.strata <- as.integer(rank(initial)[match(x = data$association.strata, 
        table = initial)])
    if (mix.model < 10) {
        stop("Specification of mix.model incorrect\n")
    }
    res <- .Call("C_fitmodel", as.integer(ncomp), as.integer(nind), 
        hyp, data, logit.offset, as.matrix(design.matrix.mean[, 
            -1]), as.matrix(design.matrix.variance[, -1]), as.matrix(design.matrix.disease[, 
            -1]), control, mix.mod, pi.mod)
    new.data <- res[[1]]
    data$posterior <- new.data[, 1]
    data$mean <- new.data[, 2]
    data$var <- new.data[, 3]
    data$pr <- new.data[, 4]
    data$alpha <- new.data[, 5]
    data$pdc <- new.data[, 6]
    data$nu <- new.data[, 7]
    return(list(data = data, status = res[[2]]))
}
