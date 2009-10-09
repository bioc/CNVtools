ExpandData <-
function (batch, trait, names, signal, ncomp, association.strata = NULL) 
{
    ncohort <- length(batch)
    cn <- vector(mode = "numeric", length = 0)
    trait.f <- vector(mode = "numeric", length = 0)
    batch.f <- vector(mode = "character", length = 0)
    subject <- vector(mode = "character", length = 0)
    signal.f <- vector(mode = "numeric", length = 0)
    a.strata.f <- vector(mode = "numeric", length = 0)
    for (nc in 1:ncomp) {
        for (ns in 1:ncohort) {
            nind <- length(signal[[ns]])
            cn <- append(cn, rep(nc, nind))
            trait.f <- append(trait.f, trait[[ns]])
            batch.f <- append(batch.f, as.character(batch[[ns]]))
            subject <- append(subject, as.character(names[[ns]]))
            signal.f <- append(signal.f, signal[[ns]])
            if (!is.null(association.strata)) 
                a.strata.f <- append(a.strata.f, association.strata[[ns]])
            else a.strata.f <- append(a.strata.f, rep(1, nind))
        }
    }
    DataFrame <- data.frame(batch = batch.f, cn = cn, trait = trait.f, 
        subject = subject, signal = signal.f, alpha = 0, mean = 0, 
        var = 0, nu = 0, pr = 0, posterior = 0, pdc = 0, strats.var = as.integer(0), 
        strats.mean = as.integer(0), association.strata = a.strata.f)
    return(DataFrame)
}
