CNVtest.binary.T <-
function (signal, batch, sample = NULL, disease.status = NULL, 
    ncomp, n.H0 = 5, n.H1 = 0, output = "compact", model.mean = "~ strata(batch, cn)", 
    model.var = "~ strata(batch, cn)", model.disease = "~ cn", 
    beta.estimated = NULL, start.mean = NULL, start.var = NULL, 
    control = list(tol = 1e-05, max.iter = 3000, min.freq = 4)) 
{
    model.nu <- model.var
    start.nu <- NULL
    nind <- length(signal)
    batch <- factor(batch)
    if (length(batch) != nind) {
        stop("Specification of batches does not have the same length as the input signal\n")
    }
    cat("Attempting to cluster the data with", nind, "individuals and", 
        ncomp, "components\n")
    if ((sd(signal) > 3) || (sd(signal) < 0.33)) {
        cat("Warning: The signal you provided to the algorithm has a standard deviation significantly different from 1. \n\tTo maximize the stability of the underlying clustering algorithm we recommend you normalize this signal to make sure that the standard deviation is indeed 1.\n")
    }
    if (is.null(sample)) {
        sample <- paste(batch, c(1:nind), sep = "_")
    }
    else {
        if (length(sample) != nind) {
            stop("The sample names you have specified do not match the number of data points\n")
        }
    }
    if (is.null(disease.status)) {
        if (n.H1 > 0) {
            stop("Must specify disease.status under H1")
        }
        disease.status <- rep(0, nind)
    }
    sample.split <- split(x = sample, f = batch)
    signal.split <- split(x = signal, f = batch)
    disease.split <- split(x = disease.status, f = batch)
    batch.split <- split(x = batch, f = batch)
    trait.sample.list <- list(sample[which(disease.status == 
        0)], sample[which(disease.status == 1)])
    ncohort <- length(levels(batch))
    data.for.Ccode <- ExpandData(batch = batch.split, trait = disease.split, 
        names = sample.split, signal = signal.split, ncomp = ncomp)
    for (start.values in c("mean", "var")) {
        arg.start <- get(paste("start", start.values, sep = "."))
        if (!is.null(arg.start)) {
            if (!(class(arg.start) %in% c("numeric", "matrix"))) 
                stop("If you provided starting values for ", 
                  start.values, ", these have to be either numeric of matrix")
            if (class(arg.start) == "numeric") 
                my.l <- length(arg.start)
            if (class(arg.start) == "matrix") 
                my.l <- dim(arg.start)[2]
            if (my.l != ncomp) 
                stop("The size of starting.values for ", start.values, 
                  " (", my.l, ") does not fit the number of components (", 
                  ncomp, ")")
            if (class(arg.start) == "numeric") 
                assign(x = paste("start", start.values, sep = "."), 
                  value = matrix(data = rep(arg.start, max(n.H0, 
                    n.H1)), nrow = max(n.H0, n.H1), byrow = TRUE))
        }
    }
    model.spec <- get.model.spec("T", model.mean, model.var, 
        model.nu, ncomp = ncomp, nbatch = ncohort)[[1]]
    status.H0 <- ""
    res <- list()
    offset <- rep(0, dim(data.for.Ccode)[1])
    design.matrix.disease <- model.matrix(data = data.for.Ccode, 
        object = as.formula(model.disease))
    if (n.H0 > 0) {
        best.model.H0 <- 0
        best.lnL.H0 <- -Inf
        best.status.H0 <- ""
        for (i in 1:n.H0) {
            cat("Iteration", i, "under H0\n")
            data.for.Ccode <- EM.starting.point(data.for.Ccode)
            for (start.values in c("mean", "var", "nu")) {
                arg.start <- get(paste("start", start.values, 
                  sep = "."))
                if (!is.null(arg.start)) {
                  line.arg.start <- arg.start[i, ]
                  if (sum(is.na(line.arg.start) == 0)) 
                    data.for.Ccode[, start.values] <- line.arg.start[data.for.Ccode$cn]
                }
            }
            final.frame <- CNV.fitModel(ncomp, nind, hyp = "H0", 
                data.for.Ccode, logit.offset = offset, design.matrix.mean = matrix(0), 
                design.matrix.variance = matrix(0), design.matrix.disease, 
                mix.model = model.spec, control = control)
            if (final.frame$status == "F") 
                lnL <- -Inf
            else lnL <- getparams(final.frame$data)$lnL
            if (lnL > best.lnL.H0) {
                best.status.H0 <- final.frame$status
                best.lnL.H0 <- lnL
                best.model.H0 <- final.frame$data
            }
        }
        res$model.H0 <- getparams(best.model.H0)
        if (output == "compact") {
            res$posterior.H0 <- compact.data.frame(best.model.H0)
            res$posterior.H0 <- res$posterior.H0[match(sample, 
                res$posterior.H0$subject), ]
            res$status.H0 = best.status.H0
            if (best.status.H0 == "C" && test.posterior(frame = res$posterior.H0, 
                ncomp = ncomp) == TRUE) 
                res$status.H0 = "P"
        }
        else res$posterior.H0 <- best.model.H0
    }
    if (n.H1 > 0) {
        best.model.H1 <- 0
        best.lnL.H1 <- -Inf
        best.status.H1 <- ""
        status.H1 <- ""
        if (!is.null(beta.estimated)) 
            offset <- beta.estimated * data.for.Ccode$cn
        for (i in 1:n.H1) {
            cat("Iteration", i, "under H1\n")
            if ((i == 1) & (n.H0 > 0)) 
                data.for.Ccode <- best.model.H0
            else {
                data.for.Ccode <- EM.starting.point(data.for.Ccode)
                for (start.values in c("mean", "var", "nu")) {
                  arg.start <- get(paste("start", start.values, 
                    sep = "."))
                  if (!is.null(arg.start)) {
                    line.arg.start <- arg.start[i, ]
                    if (sum(is.na(line.arg.start) == 0)) 
                      data.for.Ccode[, start.values] <- line.arg.start[data.for.Ccode$cn]
                  }
                }
            }
            final.frame <- CNV.fitModel(ncomp, nind, hyp = "H1", 
                data = data.for.Ccode, logit.offset = offset, 
                design.matrix.mean = matrix(0), design.matrix.variance = matrix(0), 
                design.matrix.disease, mix.model = model.spec, 
                control = control)
            if (final.frame$status == "F") 
                lnL <- -Inf
            else lnL <- getparams(final.frame$data)$lnL
            if (lnL > best.lnL.H1) {
                best.status.H1 <- final.frame$status
                best.lnL.H1 <- lnL
                best.model.H1 <- final.frame$data
            }
        }
        res$model.H1 <- getparams(best.model.H1)
        res$posterior.H1 <- compact.data.frame(best.model.H1)
        res$posterior.H1 <- res$posterior.H1[match(sample, res$posterior.H1$subject), 
            ]
        res$status.H1 <- best.status.H1
        if (best.status.H1 == "C" && test.posterior(frame = res$posterior.H1, 
            ncomp = ncomp, samples.by.disease = trait.sample.list) == 
            TRUE) 
            res$status.H1 <- "P"
    }
    return(res)
}
