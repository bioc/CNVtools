get.model.spec <-
function (model.component, model.mean, model.var, model.nu, design.matrix.mean = NULL, 
    design.matrix.variance = NULL, ncomp, nbatch) 
{
    nparam.alpha <- (ncomp - 1)
    model.spec <- 0
    if (model.component != "gaussian" && !(model.component %in% 
        c("t", "T"))) {
        cat("Did not specify valid model.component - going with gaussians\n")
        model.component <- "gaussian"
    }
    model.mean <- gsub(model.mean, pattern = " ", replacement = "")
    model.var <- gsub(model.var, pattern = " ", replacement = "")
    full.model <- c("~strata(batch,cn)", "~strata(cn,batch)")
    if (model.component == "gaussian") {
        model.spec <- 10
        if (!is.null(design.matrix.mean)) 
            nparam.mean <- ncol(design.matrix.mean)
        if (!is.null(design.matrix.variance)) 
            nparam.var <- ncol(design.matrix.variance)
        if (model.mean %in% full.model) 
            nparam.mean <- nbatch * ncomp
        if (model.mean == "~strata(cn)") 
            nparam.mean <- ncomp
        if (model.var == "~ 1") 
            nparameter.var <- 1
        if (model.var == "~strata(batch)") 
            nparameter.var <- nbatch
        if (model.var == "~strata(cn)") 
            nparameter.var <- ncomp
        if (model.var %in% full.model) 
            nparameter.var <- nbatch * ncomp
        nparameter <- nparam.alpha + nparam.mean + nparam.var
    }
    if (model.component %in% c("T", "t")) {
        nparameter <- nparam.alpha
        model.nu <- gsub(model.nu, pattern = " ", replacement = "")
        if (model.mean == "~strata(cn)") {
            model.spec <- 2300
            nparameter <- nparameter + ncomp
        }
        else if (model.mean %in% full.model) {
            model.spec <- 2400
            nparameter <- nparameter + nbatch * ncomp
        }
        else {
            cat("Mean formula not recognized, going with ~ strata(batch, cn)\n")
            model.spec <- 2400
            nparameter <- nparameter + nbatch * ncomp
        }
        if (model.var == "~ 1") {
            model.spec <- model.spec + 10
            nparameter <- nparameter + 1
        }
        else if (model.var == "~strata(batch)") {
            model.spec <- model.spec + 20
            nparameter <- nparameter + nbatch
        }
        else if (model.var == "~strata(cn)") {
            model.spec <- model.spec + 30
            nparameter <- nparameter + ncomp
        }
        else if (model.var %in% full.model) {
            model.spec <- model.spec + 40
            nparameter <- nparameter + ncomp * nbatch
        }
        else {
            cat("Variance formula not recognized, going with ~ strata(batch, cn)\n")
            model.spec <- model.spec + 40
            nparameter <- nparameter + ncomp * nbatch
        }
        if (model.nu == "~ 1") {
            model.spec <- model.spec + 1
            nparameter <- nparameter + 1
        }
        else if (model.nu == "~strata(batch)") {
            model.spec <- model.spec + 2
            nparameter <- nparameter + nbatch
        }
        else if (model.nu == "~strata(cn)") {
            model.spec <- model.spec + 3
            nparameter <- nparameter + ncomp
        }
        else if (model.nu %in% full.model) {
            model.spec <- model.spec + 4
            nparameter <- nparameter + nbatch * ncomp
        }
        else {
            cat("constrained variance fitting for t not implemented, going with ~strata(batch,cn)\n")
            model.spec <- model.spec + 4
            nparameter <- nparameter + nbatch * ncomp
        }
    }
    return(list(model.spec, nparameter))
}
