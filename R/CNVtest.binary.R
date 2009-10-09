CNVtest.binary <-
function (signal, batch, sample = NULL, disease.status = NULL, 
    ncomp, n.H0 = 5, n.H1 = 0, output = "compact", model.mean = "~ strata(batch, cn)", 
    model.var = "~ strata(batch, cn)", model.disease = "~ cn", 
    association.test.strata = NULL, beta.estimated = NULL, start.mean = NULL, 
    start.var = NULL, control = list(tol = 1e-05, max.iter = 3000, 
        min.freq = 4)) 
{
    nind <- length(signal)
    batch <- factor(batch)
    if (length(batch) != nind) {
        stop("Specification of batches does not have the same length as the input signal\n")
    }
    if (!is.null(association.test.strata)) {
        if (length(association.test.strata) != nind) 
            stop("Specification of strata for association test does not have the same length as the input signal\n")
        association.test.strata <- as.numeric(factor(association.test.strata))
    }
    message(paste("Attempting to cluster the data with", nind, 
        "individuals and", ncomp, "components"))
    if ((sd(signal) > 3) || (sd(signal) < 0.33)) {
        warning("The signal you provided to the algorithm has a standard deviation significantly different from 1. \n\tTo maximize the stability of the underlying clustering algorithm we recommend you normalize this signal to make sure that the standard deviation is indeed 1.")
    }
    if (is.null(sample)) {
        sample <- paste(batch, c(1:nind), sep = "_")
    }
    else {
        if (length(sample) != nind) {
            stop("The sample names you have specified do not match the number of data points\n")
        }
        if (sum(make.unique(as.character(sample)) != as.character(sample)) > 
            0) {
            warning("Sample names not unique. CNVtools will modify them for uniqueness.")
            sample <- make.unique(as.character(sample))
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
    if (!is.null(association.test.strata)) 
        association.test.strata.split <- split(x = association.test.strata, 
            f = batch)
    else association.test.strata.split <- NULL
    trait.sample.list <- list(sample[which(disease.status == 
        0)], sample[which(disease.status == 1)])
    ncohort <- length(levels(batch))
    data.for.Ccode <- ExpandData(batch = batch.split, trait = disease.split, 
        names = sample.split, signal = signal.split, ncomp = ncomp, 
        association.strata = association.test.strata.split)
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
    design.matrix.disease <- model.matrix(data = data.for.Ccode, 
        object = as.formula(model.disease))
    design.matrix.mean <- NULL
    design.matrix.var <- NULL
    special <- c("strata")
    for (design in c("var", "mean")) {
        my.formula <- as.formula(get(paste("model", design, sep = ".")))
        Terms <- terms(my.formula, special, data = data.for.Ccode)
        strats <- attr(Terms, "specials")$strata
        if (!is.null(strats)) {
            m <- list()
            m[[1]] <- as.name("model.frame")
            m[[2]] <- Terms
            names(m)[2] <- "formula"
            m[["data"]] <- data.for.Ccode
            m <- as.call(c(as.list(m), list(na.action = as.symbol("na.omit"))))
            m <- eval(m)
            temps <- untangle.specials(Terms, "strata", 1)
            data.for.Ccode[, paste("strats", design, sep = ".")] <- as.integer(strata(m[, 
                temps$vars], shortlabel = TRUE))
            nstrata <- length(unique(data.for.Ccode[, paste("strats", 
                design, sep = ".")]))
            if (nstrata == 1) 
                assign(x = paste("design.matrix", design, sep = "."), 
                  value = model.matrix(data = data.for.Ccode, 
                    object = as.formula(" ~ 1")))
            else assign(x = paste("design.matrix", design, sep = "."), 
                value = matrix(0))
        }
        else {
            assign(x = paste("design.matrix", design, sep = "."), 
                value = model.matrix(data = data.for.Ccode, object = my.formula))
        }
    }
    model.spec <- 10
    status.H0 <- ""
    res <- list()
    offset <- rep(0, dim(data.for.Ccode)[1])
    if (n.H0 > 0) {
        best.model.H0 <- 0
        best.lnL.H0 <- -Inf
        best.status.H0 <- ""
        for (i in 1:n.H0) {
            message(paste("Iteration", i, "under H0"))
            data.for.Ccode <- EM.starting.point(data.for.Ccode)
            for (start.values in c("mean", "var")) {
                arg.start <- get(paste("start", start.values, 
                  sep = "."))
                if (!is.null(arg.start)) {
                  line.arg.start <- arg.start[i, ]
                  if (sum(is.na(line.arg.start) == 0)) 
                    data.for.Ccode[, start.values] <- line.arg.start[data.for.Ccode$cn]
                }
            }
            final.frame <- CNV.fitModel(ncomp, nind, hyp = "H0", 
                data.for.Ccode, logit.offset = offset, design.matrix.mean, 
                design.matrix.var, design.matrix.disease, mix.model = model.spec, 
                control = control)
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
        res$model.H0$nu <- 0
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
            message(paste("Iteration", i, "under H1"))
            if ((i == 1) & (n.H0 > 0)) 
                data.for.Ccode <- best.model.H0
            else {
                data.for.Ccode <- EM.starting.point(data.for.Ccode)
                for (start.values in c("mean", "var")) {
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
                design.matrix.mean, design.matrix.var, design.matrix.disease, 
                mix.model = model.spec, control = control)
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
        if (output == "compact") {
            res$posterior.H1 <- compact.data.frame(best.model.H1)
            res$posterior.H1 <- res$posterior.H1[match(sample, 
                res$posterior.H1$subject), ]
        }
        else {
            res$posterior.H1 <- best.model.H1
        }
        res$status.H1 <- best.status.H1
        if (best.status.H1 == "C" && test.posterior(frame = res$posterior.H1, 
            ncomp = ncomp, samples.by.disease = trait.sample.list) == 
            TRUE) 
            res$status.H1 <- "P"
        res$model.H1$nu <- 0
    }
    return(res)
}
