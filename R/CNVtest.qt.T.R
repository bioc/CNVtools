CNVtest.qt.T <-
function (signal, batch, sample = NULL, qt = NULL, ncomp, n.H0 = 5, 
    n.H1 = 0, model.mean = "~ strata(cn)", model.var = "~ strata(cn)", 
    model.qt = "~ cn", beta.estimated = NULL, start.mean = NULL, 
    start.var = NULL, control = list(tol = 1e-05, max.iter = 3000, 
        min.freq = 4)) 
{
    model.nu <- model.var
    start.nu <- NULL
    nind <- length(signal)
    batch <- factor(batch)
    if (length(batch) != nind) {
        stop("Specification of batches does not have the same length as the input signal\n")
    }
    if (length(qt) != nind) {
        stop("Specification of qt does not have the same length as the input signal\n")
    }
    if (is.null(sample)) {
        sample <- paste(batch, c(1:nind), sep = "_")
    }
    else {
        if (length(sample) != nind) {
            stop("The sample names you have specified do not match the number of data points\n")
        }
    }
    if (is.null(qt)) {
        if (n.H1 > 0) {
            stop("Must specify quantitative trait under H1")
        }
        qt <- rep(0, nind)
    }
    sample.split <- split(x = sample, f = batch)
    signal.split <- split(x = signal, f = batch)
    qt.split <- split(x = qt, f = batch)
    batch.split <- split(x = batch, f = batch)
    ncohort <- length(levels(batch))
    data.for.Ccode <- ExpandData(batch = batch.split, trait = qt.split, 
        names = sample.split, signal = signal.split, ncomp = ncomp)
    if (!missing(start.mean) && !is.null(start.mean)) {
        if ((length(start.mean) != ncomp) || (class(start.mean) != 
            "numeric")) {
            stop("You provided invalid starting values for the mean")
        }
    }
    else {
        start.mean <- NULL
    }
    if (!is.null(start.var) && !missing(start.var)) {
        if ((length(start.var) != ncomp) || (class(start.var) != 
            "numeric")) {
            stop("You provided invalid starting values for the variances")
        }
    }
    else {
        start.var <- NULL
    }
    if (!is.null(start.nu) && !missing(start.nu)) {
        if ((length(start.nu) != ncomp) || (class(start.nu) != 
            "numeric")) {
            stop("You provided invalid starting values for the nu")
        }
    }
    else {
        start.nu <- NULL
    }
    design.matrix.mean <- as.matrix(0)
    design.matrix.variance <- as.matrix(0)
    design.matrix.qt <- model.matrix(data = data.for.Ccode, object = as.formula(model.qt))
    model.spec <- get.model.spec("T", model.mean, model.var, 
        model.nu, design.matrix.mean, design.matrix.variance, 
        ncomp, ncohort)[[1]]
    status.H0 <- ""
    res <- list()
    offset <- rep(0, dim(data.for.Ccode)[1])
    if (n.H0 > 0) {
        best.model.H0 <- 0
        best.lnL.H0 <- -Inf
        best.status.H0 <- ""
        for (i in 1:n.H0) {
            cat("Iteration", i, "under H0\n")
            data.for.Ccode <- EM.starting.point(data.for.Ccode)
            if (!is.null(start.mean)) 
                data.for.Ccode$mean = start.mean[data.for.Ccode$cn]
            if (!is.null(start.var)) 
                data.for.Ccode$var = start.var[data.for.Ccode$cn]
            if (!is.null(start.nu)) 
                data.for.Ccode$nu = start.nu[data.for.Ccode$cn]
            final.frame <- CNV.fitModel(ncomp, nind, hyp = "H0", 
                data.for.Ccode, logit.offset = offset, design.matrix.mean = matrix(0), 
                design.matrix.variance = matrix(0), design.matrix.qt, 
                pi.model = 2, mix.model = model.spec, control = control)
            if (final.frame[["status"]] == "F") 
                lnL <- -Inf
            else lnL <- getparams(final.frame[["data"]])$lnL
            if (lnL > best.lnL.H0) {
                best.status.H0 <- final.frame[["status"]]
                best.lnL.H0 <- lnL
                best.model.H0 <- final.frame[["data"]]
            }
        }
        res$model.H0 = getparams(best.model.H0)
        res$posterior.H0 = compact.data.frame(best.model.H0)
        res$posterior.H0 <- res$posterior.H0[match(sample, res$posterior.H0$subject), 
            ]
        res[["status.H0"]] = best.status.H0
        if (best.status.H0 == "C" && test.posterior(res[["posterior.H0"]], 
            ncomp) == TRUE) 
            res[["status.H0"]] = "P"
    }
    if (n.H1 > 0) {
        best.model.H1 <- 0
        best.lnL.H1 <- -1e+09
        best.status.H1 <- ""
        status.H1 <- ""
        if (!is.null(beta.estimated)) 
            offset <- beta.estimated * data.for.Ccode$cn
        for (i in 1:n.H1) {
            cat("Iteration", i, "under H1\n")
            if ((i == 1) & (n.H0 > 0)) {
                data.for.Ccode <- best.model.H0
            }
            else {
                data.for.Ccode <- EM.starting.point(data.for.Ccode)
                if (!is.null(start.mean)) 
                  data.for.Ccode$mean = start.mean[data.for.Ccode$cn]
                if (!is.null(start.var)) 
                  data.for.Ccode$var = start.var[data.for.Ccode$cn]
            }
            final.frame <- CNV.fitModel(ncomp, nind, hyp = "H1", 
                data = data.for.Ccode, logit.offset = offset, 
                design.matrix.mean, design.matrix.variance, design.matrix.qt, 
                pi.model = 2, mix.model = model.spec, control = control)
            if (final.frame[["status"]] == "F") 
                lnL <- -1e+09
            else lnL <- getparams(final.frame[["data"]])$lnL
            if (lnL > best.lnL.H1) {
                best.status.H1 <- final.frame[["status"]]
                best.lnL.H1 <- lnL
                best.model.H1 <- final.frame[["data"]]
            }
        }
        res[["model.H1"]] = getparams(best.model.H1)
        res[["posterior.H1"]] = compact.data.frame(best.model.H1)
        res[["status.H1"]] <- best.status.H1
    }
    return(res)
}
