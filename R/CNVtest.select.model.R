CNVtest.select.model <-
function (signal, batch, sample = NULL, n.H0 = 3, method = "BIC", 
    v.ncomp = 1:6, v.model.component = rep("gaussian", 6), v.model.mean = rep("~ strata(cn)", 
        6), v.model.var = rep("~1", 6), control = list(tol = 1e-05, 
        max.iter = 500, min.freq = 4)) 
{
    start.values <- list()
    v.model.nu <- v.model.var
    y <- c(length(v.ncomp), length(v.model.component), length(v.model.mean), 
        length(v.model.var), length(v.model.nu))
    if (sum(ifelse(y == y[1], 1, 0)) != 5) {
        stop("ncomp, model.mean, model.var, model.nu, must have same length")
    }
    bics <- vector(mode = "numeric", length = 0)
    aics <- vector(mode = "numeric", length = 0)
    nind <- length(signal)
    batch <- factor(batch)
    if (length(batch) != nind) {
        stop("Specification of batches does not have the same length as the input signal\n")
    }
    if (is.null(sample)) {
        sample <- paste(batch, c(1:nind), sep = "_")
    }
    else {
        if (length(sample) != nind) {
            stop("The sample names you have specified do not match the number of data points\n")
        }
    }
    sample.split <- split(x = sample, f = batch)
    signal.split <- split(x = signal, f = batch)
    disease.split <- split(x = rep(0, nind), f = batch)
    batch.split <- split(x = batch, f = batch)
    nbatches <- length(levels(batch))
    ncohort <- length(levels(batch))
    nmodel <- length(v.ncomp)
    res <- list()
    res[["model"]] <- vector("list", nmodel)
    res[["BIC"]] <- vector("numeric", nmodel)
    res[["AIC"]] <- vector("numeric", nmodel)
    res[["status"]] <- vector("character", nmodel)
    res[["np"]] <- vector("numeric", nmodel)
    res[["posteriors"]] <- vector("list", nmodel)
    res[["models"]] <- vector("list", nmodel)
    res[["selected"]] <- 0
    for (m in c(1:nmodel)) {
        message(paste("Fitting model ", m))
        status <- ""
        ncomp <- v.ncomp[m]
        model.component <- v.model.component[m]
        model.mean <- v.model.mean[m]
        model.var <- v.model.var[m]
        model.nu <- v.model.nu[m]
        data.for.Ccode <- ExpandData(batch = batch.split, trait = disease.split, 
            names = sample.split, signal = signal.split, ncomp = ncomp, 
            association.strata = NULL)
        design.matrix.disease <- model.matrix(data = data.for.Ccode, 
            object = as.formula("~ 1"))
        if (model.component == "gaussian") {
            special <- c("strata")
            nparameters <- ncomp - 1
            my.model <- 10
            for (design in c("var", "mean")) {
                my.formula <- as.formula(get(paste("model", design, 
                  sep = ".")))
                Terms <- terms(my.formula, special, data = data.for.Ccode)
                strats <- attr(Terms, "specials")$strata
                if (!is.null(strats)) {
                  mm <- list()
                  mm[[1]] <- as.name("model.frame")
                  mm[[2]] <- Terms
                  names(mm)[2] <- "formula"
                  mm[["data"]] <- data.for.Ccode
                  mm <- as.call(c(as.list(mm), list(na.action = as.symbol("na.omit"))))
                  mm <- eval(mm)
                  temps <- untangle.specials(Terms, "strata", 
                    1)
                  data.for.Ccode[, paste("strats", design, sep = ".")] <- as.integer(strata(mm[, 
                    temps$vars], shortlabel = TRUE))
                  nstrata <- length(unique(data.for.Ccode[, paste("strats", 
                    design, sep = ".")]))
                  nparameters <- nparameters + nstrata
                  if (nstrata == 1) 
                    assign(x = paste("design.matrix", design, 
                      sep = "."), value = model.matrix(data = data.for.Ccode, 
                      object = as.formula(" ~ 1")))
                  else assign(x = paste("design.matrix", design, 
                    sep = "."), value = matrix(0))
                }
                else {
                  assign(x = paste("design.matrix", design, sep = "."), 
                    value = model.matrix(data = data.for.Ccode, 
                      object = my.formula))
                  nparameters <- nparameters + ncol(get(paste("design.matrix", 
                    design, sep = ".")))
                }
            }
        }
        else {
            design.matrix.mean <- as.matrix(0)
            design.matrix.var <- as.matrix(0)
            model.spec <- get.model.spec(model.component, model.mean, 
                model.var, model.nu, design.matrix.mean, design.matrix.var, 
                ncomp, ncohort)
            my.model <- model.spec[[1]]
            nparameters <- model.spec[[2]]
        }
        best.model.H0 <- 0
        best.lnL.H0 <- -Inf
        best.status <- "F"
        offset <- rep(0, dim(data.for.Ccode)[1])
        suc <- 0
        start.mean <- NULL
        start.var <- NULL
        start.nu <- NULL
        mean.label <- paste("M", ncomp, sep = "")
        if (mean.label %in% names(start.values)) 
            start.mean <- start.values[[mean.label]]
        for (i in 1:n.H0) {
            message(paste("\t Iteration", i))
            data.for.Ccode <- EM.starting.point(data.for.Ccode)
            for (start.v in c("mean", "var", "nu")) {
                arg.start <- get(paste("start", start.v, sep = "."))
                if (!is.null(arg.start)) {
                  line.arg.start <- arg.start[i, ]
                  if (sum(is.na(line.arg.start) == 0)) 
                    data.for.Ccode[, start.v] <- line.arg.start[data.for.Ccode$cn]
                }
            }
            final.frame <- CNV.fitModel(ncomp, nind, hyp = "H0", 
                data = data.for.Ccode, logit.offset = offset, 
                design.matrix.mean, design.matrix.var, design.matrix.disease, 
                mix.model = my.model, control = control)
            if (final.frame[["status"]] == "F") 
                next
            suc <- suc + 1
            lnL <- getparams(final.frame[["data"]])$lnL
            if (lnL > best.lnL.H0) {
                best.lnL.H0 <- lnL
                best.model.H0 <- final.frame[["data"]]
                best.status <- final.frame[["status"]]
            }
        }
        if (suc > 0) {
            bic <- -2 * best.lnL.H0 + nparameters * log(nind)
            aic <- -2 * best.lnL.H0 + 2 * nparameters
            bics <- append(bics, bic)
            aics <- append(aics, aic)
            res[["np"]][[m]] <- nparameters
            res[["model"]][[m]] <- list(ncomp = ncomp, model.mean = model.mean, 
                model.var = model.var)
            res[["BIC"]][[m]] <- bic
            res[["AIC"]][[m]] <- aic
            res[["lnL"]][[m]] <- best.lnL.H0
            res[["posteriors"]][[m]] <- compact.data.frame(best.model.H0)
            res[["models"]][[m]] <- getparams(best.model.H0)
            if (best.status == "C" && ncomp != 1 && test.posterior(frame = compact.data.frame(best.model.H0), 
                ncomp = ncomp) == TRUE) 
                res[["status"]][[m]] = "P"
            else res[["status"]][[m]] = best.status
            if (method == "BIC") {
                res[["selected"]] <- order(bics)[1]
            }
            else {
                res[["selected"]] <- order(aics)[1]
            }
        }
    }
    return(res)
}
