## Copyright (C) 2012 Marius Hofert and Valerie Chavez-Demoulin
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


## *G*eneralized *A*dditive *M*odeling for *E*xtremes
##
## Remarks:
## 1) Notation:
##    paper | Valerie's thesis | Valerie's former code | Coles' book + ismev
##      xi        -kappa                kappa                    xi
##    beta         sigma                sigma                 sigma
##      nu            nu                  eta                    --
## 2) GPD(xi in IR, beta > 0) distribution function (same as in ismev + EKM; up to notation):
##    G_{xi,beta}(x) = 1-(1+xi*x/beta)^(-1/xi) if xi!=0
##                   = 1-exp(-x/beta)          if xi =0
##    x>=0 when xi>=0 and x in [0,-beta/xi] when xi<0
##    Note: in EKM, x ~> (x-nu)/beta (with the same meaning of beta)
## 3) GPD(xi, beta) density for xi>0:
##    g_{xi,beta}(x) = (1+xi*x/beta)^(-(1+1/xi))/beta if x>0
## 4) We need xi > -1 in order for beta (= exp(nu) / (1+xi)) > 0. In comparison
##    to Chavez-Demoulin, Embrechts, Hofert (2016), the below code was updated
##    in 2021 ('QRM' version 0.4-33) to also include a reparameterization in xi
##    (= exp(eta) - 1).


### Auxiliary functions ########################################################

##' Reparameterized log-likelihood l^r(eta, nu; y_1,..,y_n) = l(xi, beta;
##' y_1,..,y_n) = l(exp(eta)-1, exp(nu)/(1+xi); y_1,..,y_n) = l(exp(eta)-1, exp(nu-eta);
##' y_1,..,y_n), where l(xi, beta; y_1,..,y_n) denotes the log-likelihood of a
##' GPD(xi, beta) distribution
##'
##' @title Reparameterized log-Likelihood l^r
##' @param y vector of excesses (over a high threshold u)
##' @param eta vector of GPD(xi, beta) parameters eta (= log(1+xi))
##' @param nu vector of GPD(xi, beta) parameters nu (= log(beta*(1+xi)) = log(beta)+eta);
##'        orthogonal to eta in the Fisher information metric: mixed derivative
##'        (w.r.t. eta and nu) of the reparameterized log-likelihood is 0
##' @return reparametrized log-likelihood l^r
##' @author Marius Hofert
rlogL <- function(y, eta, nu, verbose = TRUE)
{
    stopifnot((n <- length(y)) > 0, length(eta)==n, length(nu)==n)
    ## Before reparametrization of xi
    ## ii <- xi <= -1
    ## if(all(ii)) stop("Can't adjust xi <= -1 since there are no xi > -1") # beta = exp(nu)/(1+xi)
    ## if(any(ii)) { # not a valid reparameterization either (for all xi <= -1) but we make it one
    ##     perc <- 100*sum(ii)/length(ii) # percentage of xi <= -1
    ##     if(verbose) warning(round(perc,2), "% of all xi are <= -1 and are adjusted to be > -1")
    ##     xi[ii] <- mean(xi[!ii])
    ## }
    xi <- expm1(eta)
    res <- rep(-Inf, n)
    ii <- (xi < 0 & (0 <= y & y < -exp(nu)/(xi*(1+xi)))) | (xi > 0 & y >= 0)
    if(any(ii)) {
        xi. <- xi[ii]
        nu. <- nu[ii]
        res[ii] <- log1p(xi.)-nu.-(1+1/xi.)*log1p(xi.*(1+xi.)*exp(-nu.)*y[ii])
    }
    ii <- xi == 0
    if(any(ii)) res[ii] <- -(nu[ii]+exp(-nu[ii])*y[ii])
    sum(res)
}

##' @title Function to Adjust Derivatives
##' @param x vector of values (derivatives)
##' @param order order of the derivatives to be adjusted (1, 2)
##' @param verbose logical indicating whether warnings about adjustments of
##'        the derivatives are printed
##' @return adjusted derivatives
##' @author Marius Hofert
##' Note: This is an auxiliary function of DrlogL()
adjustD <- function(x, order, verbose = TRUE)
{
    stopifnot(order==1 || order==2)

    ## order == 1
    ii <- is.finite(x) # indices of value which are fine (not NA, NaN, Inf, or -Inf)
    res <- x
    if(any(!ii)) {
        if(sum(ii)==0) stop("Can't adjust derivatives, there are no finite values")
        if(verbose) warning(sum(!ii)," (",
                            round(100*sum(!ii)/length(res), 2),
                            "%) non-finite derivatives adjusted")
        res[!ii] <- mean(res[ii])
    }

    ## order == 2 => additionally check that 2nd order derivatives are negative
    ii <- res < 0
    if(order == 2 && any(!ii)) {
        if(sum(ii)==0) stop("Can't adjust second-order derivatives, there are no negative values")
        if(verbose) warning(sum(!ii)," (",
                            round(100*sum(!ii)/length(res), 2),
                            "%) non-negative second-order derivatives adjusted")
        res[!ii] <- mean(res[ii])
    }

    ## return
    res
}

##' @title Compute Derivatives of the Reparameterized log-Likelihood
##' @param y vector of excesses (over a high threshold u)
##' @param eta vector of GPD(xi, beta) parameters eta (= log(1+xi))
##' @param nu vector of GPD(xi, beta) parameters nu (= log(beta*(1+xi)) = log(beta)+eta)
##' @param adjust logical indicating whether non-real values of the derivatives are adjusted
##' @param verbose logical indicating whether modified arguments are printed
##' @return (n x 4) matrix containing the partial derivatives of the
##'         reparameterized log-likelihood l^r where
##'         column 1: derivative of the reparameterized log-likelihood w.r.t. eta
##'         column 2: derivative of the reparameterized log-likelihood w.r.t. nu
##'         column 3: 2nd derivative of the reparameterized log-likelihood w.r.t. eta
##'         column 4: 2nd derivative of the reparameterized log-likelihood w.r.t. nu
##' @author Marius Hofert
##' Note: Column 3 and 4 have different sign than in the old code
DrlogL <- function(y, eta, nu, adjust = TRUE, verbose = TRUE)
{
    stopifnot((n <- length(y)) > 0, length(eta)==n, length(nu)==n)

    ## Main
    res <- matrix(0, nrow = n, ncol = 4,
                  dimnames = list(NULL, c("rl.eta", "rl.nu", "rl.etaeta", "rl.nunu")))
    xi <- expm1(eta)
    ii <- (xi < 0 & (0 <= y & y < -exp(nu)/(xi*(1+xi)))) | (xi > 0 & y >= 0)
    if(any(ii)) {
        ## auxiliary terms
        eta. <- eta[ii]
        xi. <- xi[ii]
        nu. <- nu[ii]
        y.  <-  y[ii]
        a <- 1+xi.*(1+xi.)*exp(-nu.)*y.
        la <- log1p(xi.*(1+xi.)*exp(-nu.)*y.)
        a. <- exp(-nu.)*y.*(1+2*xi.)
        ## main
        res[ii,"rl.eta"] <- exp(eta.)*(1/(1+xi.) + la/xi.^2 - (1+1/xi.)*a./a) # formerly: 1/(1+xi.) + la/xi.^2 - (1+1/xi.)*a./a
        res[ii,"rl.etaeta"] <- exp(eta.)*(-1/(1+xi.)^2 - 2*la/xi.^3 + 2*a./(a*xi.^2) -
                                          (1+1/xi.)*exp(-nu.)*y. * (2*a - (1+2*xi.)^2 * exp(-nu.)*y.) / a^2) +
            exp(eta.)*(1/(1+xi.) + la/xi.^2 - (1+1/xi.)*a./a) # formerly: -1/(1+xi.)^2 - 2*la/xi.^3 + 2*a./(a*xi.^2) - (1+1/xi.)*exp(-nu.)*y. * (2*a - (1+2*xi.)^2 * exp(-nu.)*y.) / a^2
        res[ii,"rl.nu"] <- (-1+(1+xi.)*exp(-nu.)*y.) / a
        res[ii,"rl.nunu"] <- -(1+xi.)^2*y.*exp(-nu.)/a^2
    }
    ii <- xi == 0
    if(any(ii)) {
        ## auxiliary terms
        ynu <- y[ii]*exp(-nu[ii])
        res[ii,"rl.eta"] <- 0
        res[ii,"rl.etaeta"] <- 0
        res[ii,"rl.nu"] <- -1+ynu
        res[ii,"rl.nunu"] <- -ynu
    }

    ## replace non-finite derivatives by mean of all finite ones
    if(adjust) {
        res[,"rl.eta"]    <- adjustD(res[,"rl.eta"],    order = 1, verbose = verbose)
        res[,"rl.etaeta"] <- adjustD(res[,"rl.etaeta"], order = 2, verbose = verbose)
        res[,"rl.nu"]     <- adjustD(res[,"rl.nu"],     order = 1, verbose = verbose)
        res[,"rl.nunu"]   <- adjustD(res[,"rl.nunu"],   order = 2, verbose = verbose)
    }

    ## return
    res
}

##' @title Compute Update (one Iteration) in gamGPDfit()
##' @param y data.frame containing the excesses over the threshold in a column
##'        labeled yname
##' @param eta.nu 2-column matrix of GPD parameters (eta, nu) to be updated
##' @param etaFrhs right-hand side of the formula for eta in the gam() call
##'        for fitting eta
##' @param nuFrhs right-hand side of the formula for nu in the gam() call
##'        for fitting nu
##' @param yname string containing the name of the column of y which contains
##'        the excesses
##' @param adjust logical indicating whether non-real values of the derivatives are adjusted
##' @param verbose logical indicating whether warnings about adjustments of
##'        the derivatives, wrong arguments in DrlogL() and failed gam() calls are printed
##' @param ... additional arguments passed to gam()
##' @return a list of length four containing
##'         element 1 (eta): object of class gamObject for eta as returned by mgcv::gam()
##'         element 2 (nu): object of class gamObject for nu as returned by mgcv::gam()
##'         element 3 (eta.weights): weights associated with eta
##'         element 4 (nu.weights): weights associated with nu
##'         or list() (in case gam() or the Newton step failed)
##' @author Marius Hofert
##' Note: That's a helper function of gamGPDfit()
gamGPDfitUp <- function(y, eta.nu, etaFrhs, nuFrhs, yname, adjust=TRUE, verbose=TRUE, ...)
{
    stopifnot(is.data.frame(y), (dim. <- dim(y))[2] >= 1,
              length(which(colnames(y)==yname))==1,
              (n <- dim.[1]) > 0, dim(eta.nu)==c(n, 2))

    ## pick out eta and nu (for readability)
    eta <- eta.nu[,1]
    nu  <- eta.nu[,2]

    ## compute one Newton step in xi
    DrLL <- tryCatch(DrlogL(y[,yname], eta = eta, nu = nu, adjust = adjust, verbose = verbose), # (n1,4) matrix
                     error = function(e) e)
    if(is(DrLL, "simpleError")) return(list())
    rl.eta <- DrLL[,"rl.eta"] # score in eta
    rl.etaeta. <- DrLL[,"rl.etaeta"] # -weight
    Newton.eta <- eta - rl.eta / rl.etaeta. # Newton step

    ## concatenate Newton.eta and rl.etaeta. to y, build formula, and estimate eta
    if("Newton.eta" %in% colnames(y))
        stop("y is not allowed to have a column named 'Newton.eta'")
    y. <- cbind(y, Newton.eta = Newton.eta, rl.etaeta. = rl.etaeta.)
    eta.formula <- update(etaFrhs, Newton.eta~.) # build formula Newton.eta ~ etaFrhs
    ## note: the following tryCatch() was used to avoid the error
    ##       "no valid set of coefficients has been found: please supply starting values"
    ##       when using gamGPDboot() is called with a small sample size
    eta.obj <- tryCatch(gam(eta.formula, data = y., weights = -rl.etaeta., ...),
                       error = function(e) e) # updated eta object of type gamObject
    ## Update on 2020-07-07: eta.formula's 'by' argument needs to be a factor! Otherwise:
    ## Error in smoothCon(split$smooth.spec[[i]], data, knots, absorb.cons, scale.penalty = scale.penalty, : Can't find by variable
    ## => see also https://stackoverflow.com/questions/45832928/gam-model-error
    if(is(eta.obj, "simpleError")) return(list())
    ## warning("gam() produced the error:", conditionMessage(eta.obj), " when fitting eta; propagate list() as result")

    ## build fitted (eta) object and check
    eta.fit <- fitted(eta.obj)

    ## compute one Newton step in nu (for given new eta)
    DrLL <- tryCatch(DrlogL(y[,yname], eta = eta.fit, nu = nu, adjust = adjust, verbose = verbose), # (n1,4) matrix
                     error = function(e) e)
    if(is(DrLL, "simpleError")) return(list())
    rl.nu <- DrLL[,"rl.nu"] # score in nu
    rl.nunu. <- DrLL[,"rl.nunu"] # -weight
    Newton.nu <- nu - rl.nu / rl.nunu. # Newton step

    ## concatenate Newton.nu and rl.nunu. to y, build formula, and estimate nu
    if("Newton.nu" %in% colnames(y))
        stop("y is not allowed to have a column named 'Newton.nu'")
    y. <- cbind(y, Newton.nu=Newton.nu, rl.nunu.=rl.nunu.)
    nu.formula <- update(nuFrhs, Newton.nu~.) # build formula Newton.nu ~ nuFrhs
    nu.obj <- tryCatch(gam(nu.formula, data = y., weights = -rl.nunu., ...),
                       error = function(e) e) # updated nu object of type gamObject
    if(is(nu.obj, "simpleError")) return(list())
    ## warning("gam() produced the error:", conditionMessage(nu.obj), " when fitting nu; propagate list() as result")

    ## return list of two gamObject objects (for eta, nu)
    list(eta = eta.obj, nu = nu.obj, eta.weights = -rl.etaeta., nu.weights = -rl.nunu.)
    ## naming (eta, nu) not ideal but guarantees colnames(param.old) = colnames(param.new) in gamGPDfit()
}


### gamGPDfit() ################################################################

##' @title Semi-parametric Estimation of GPD Parameters via Penalized
##'        Maximum Likelihood Estimation Based on Reweighted Least Squares (Backfitting)
##' @param x data.frame containing the losses (all other columns are treated
##'        as covariates)
##' @param threshold POT threshold above which losses are considered
##' @param nextremes number of excesses
##' @param datvar name of the data column which contains the data to be modeled,
##'        for example, the losses
##' @param etaFrhs right-hand side of the formula for eta in the gam() call
##'        for fitting eta
##' @param nuFrhs right-hand side of the formula for nu in the gam() call
##'        for fitting nu
##' @param init bivariate vector containing initial values for (xi, beta)
##' @param niter maximal number of iterations in the backfitting algorithm
##' @param include.updates logical indicating whether updates for eta and nu are
##'        returned as well
##' @param eps.eta epsilon for stopping criterion for eta
##' @param eps.nu epsilon for stopping criterion for nu
##' @param progress logical indicating whether progress information is displayed
##' @param adjust logical indicating whether non-real values of the derivatives are adjusted
##' @param verbose logical passed to gamGPDfitUp() (thus to DrlogL() and
##'        adjustD())
##' @param ... additional arguments passed to gam() (called by gamGPDfitUp())
##' @return a list (see below) or list() (in case gam() or the Newton step in gamGPDfitUp() failed)
##' @author Marius Hofert
gamGPDfit <- function(x, threshold, nextremes = NULL, datvar, etaFrhs, nuFrhs,
                      init = fit.GPD(x[,datvar], threshold=threshold, type="pwm", verbose=FALSE)$par.ests,
                      niter = 32, include.updates = FALSE, eps.eta = 1e-5, eps.nu = 1e-5,
                      progress = TRUE, adjust = TRUE, verbose = FALSE, ...)
{
    ## checks
    stopifnot(is.data.frame(x), length(init) == 2, niter >= 1, eps.eta > 0, eps.nu > 0)
    has.threshold <- !missing(threshold)
    has.nextr <- !is.null(nextremes)
    if(has.threshold && has.nextr)
        warning("Only one of 'threshold' and 'nextremes' is allowed -- will take 'threshold'") # both threshold and nextremes given
    if(!has.threshold && !has.nextr)
        stop("Provide either 'threshold' or 'nextremes'") # none of threshold or nextremes given
    dim. <- dim(x) # dimension of x
    stopifnot((n <- dim.[1])>=1, dim.[2]>=2) # there should at least be one observation (actually, quite a bit more) and one column of covariates
    if(has.nextr) { # nextremes given but no threshold
        stopifnot(0 < nextremes, nextremes <= n)
        threshold <- quantile(x[,datvar], probs=1-nextremes/n, names=FALSE)
    } # => now we can work with threshold

    ## determine excesses
    y. <- x[x[,datvar]>threshold,] # pick out excesses; note: y. still contains covariates
    y.[,datvar] <- y.[,datvar] - threshold # replace excesses by excess sizes
    n.ex <- nrow(y.) # number of excesses

    ## initial values for xi, eta, beta, and nu (nu = reparameterization of beta)
    xi.init <- init[1]
    beta.init <- init[2]
    nu.init <- log((1+xi.init) * beta.init)
    eta.init <- log(1+xi.init)

    ## iteration ###############################################################

    iter <- 1 # iteration number
    updates <- list() # (empty) list of update (gamGPDfitUp) objects
    while(TRUE){

        ## update/fit parameters (eta, nu)
        param.old <- if(iter==1)
            matrix(rep(c(eta.init, nu.init), each=n.ex), ncol=2,
                       dimnames=list(rownames(y.), c("eta", "nu"))) # (n.ex,2)-matrix with cols "eta" and "nu"
        else param.new
        updates[[iter]] <- gamGPDfitUp(y., eta.nu = param.old,
                                       etaFrhs = etaFrhs, nuFrhs = nuFrhs,
                                       yname = datvar, adjust = adjust, verbose = verbose, ...) # returns a list of two gam() objects containing the fitted (eta, nu) (or list() in case a gam() call fails)

        ## check whether gam() calls failed
        if(!length(updates[[iter]])) {
            warning("gam() call(s) failed, will return list()")
            return(list()) # one or both gam() calls in gamGPDfitUp() failed
        }
        ## check convergence status of gam() in gamGPDfitUp()
        conv <- sapply(updates[[iter]][c("eta", "nu")], function(x) x$converged)
        if(any(!conv))
            warning("gam() in gamGPDfitUp() did not converge for ", paste(c("eta","nu")[!conv], collapse=", "))
        ## check parameter dimensions (param.old and param.new have same rownames/colnames)
        param.new <- sapply(updates[[iter]][c("eta", "nu")], fitted)
        if(any(dim(param.new)!=dim(param.old)))
            stop("gamGPDfitUp() returned an updated gamObject object of wrong dimension")

        ## tracing
        MRD.. <- colMeans(abs((param.old-param.new)/param.old)) # mean relative distance for (eta, nu)
        MRD. <- c(iter=iter, MRD..[1], MRD..[2])
        MRD <- if(iter==1) MRD. else rbind(MRD, MRD., deparse.level=0)
        if(progress){
            cat("Mean relative differences in iteration ", iter,
                " (eta, nu): (", sprintf("%g", MRD..[1]), ", ",
                sprintf("%g", MRD..[2]), ")\n", sep="")
        }

        ## check for "convergence"
        conv <- c(MRD..[1] <= eps.eta, MRD..[2] <= eps.nu)
        if(all(conv)) break
        if(iter >= niter){
            if(progress) warning("Reached 'niter' without the required precision for ",
                                 paste(c("eta","nu")[!conv], collapse=", "))
            break
        }
        iter <- iter + 1 # update iteration

    }

    ## finish and return #######################################################

    ## fitted parameters
    eta <- param.new[,"eta"] # fitted eta's
    xi <- expm1(eta)
    nu <- param.new[,"nu"] # fitted nu's
    beta <- exp(nu)/(1+xi) # corresponding (fitted) beta (= exp(nu-eta))

    ## check beta for being a valid reparameterization (otherwise pGPD() fails)
    ## note: nu can be too small (=> beta = 0) => set corresponding beta to NA
    ii <- is.finite(beta) & beta > 0
    beta[!ii] <- NA

    ## compute the residuals
    resi <- mapply(function(x, xi, beta) if(!is.na(beta)) -log1p(-pGPD(x, xi=xi, beta=beta)) else NA,
                   x=y.[,datvar], xi=xi, beta=beta)

    ## log-likelihood
    logL <- rlogL(y=y.[,datvar], eta=eta, nu=nu) # reparameterized log-likelihood (reparameterization doesn't matter, it's the same *value* as the original log-likelihood)

    ## fitted gamObjects for eta and nu and standard errors
    ## (contained in updates but for convenience we treat this additionally)
    etaObj <- updates[[iter]]$eta
    nuObj  <- updates[[iter]]$nu
    se.eta <- updates[[iter]]$eta.weights
    se.nu  <- updates[[iter]]$nu.weights

    ## determine *all* combinations of the provided covariates
    ## eta
    eta.covars <- names(etaObj$var.summary)
    x.eta <- x[,eta.covars, drop = FALSE] # all cols with covariates used for fitting eta
    all.covars.eta <- lapply(seq_len(ncol(x.eta)), function(j) sort(unique(x.eta[,j]))) # determine levels
    ## => we only can determine those which appear at least in one column
    names(all.covars.eta) <- colnames(x.eta) # put in names
    ## nu
    nu.covars <- names(nuObj$var.summary)
    x.nu <- x[,nu.covars, drop=FALSE] # all cols with covariates used for fitting nu
    all.covars.nu <- lapply(seq_len(ncol(x.nu)), function(j) sort(unique(x.nu[,j]))) # determine levels
    names(all.covars.nu) <- colnames(x.nu) # put in names

    ## build list
    res <- list(xi=xi, # estimated xi
                beta=beta, # estimated beta
                eta=eta, # estimated eta
                nu=nu, # estimated nu
                se.eta=se.eta, # standard error for eta
                se.nu=se.nu, # standard error for nu
                eta.covar=all.covars.eta, # (unique) covariates for eta
                nu.covar=all.covars.nu, # (unique) covariates for nu
                covar=y.[,union(eta.covars, nu.covars), drop=FALSE], # *available* (not necessarily all) covariate combinations used for fitting beta (= eta *and* nu)
                y=y.[,datvar], # excesses
                res=resi, # residuals
                MRD=MRD, # mean relative distances between old/new (eta, nu) for all iterations
                logL=logL, # log-likelihood at the estimated parameters
                etaObj=etaObj, # gamObject for estimated eta (return object of mgcv::gam())
                nuObj=nuObj) # gamObject for estimated nu (return object of mgcv::gam())
    if(include.updates) res <- c(res, etaUpdates=lapply(updates, `[[`, "eta"), # updates for eta for each iteration (list of gamObject objects); contains etaObj as last element
                                 nuUpdates=lapply(updates, `[[`, "nu")) # updates for nu for each iteration (list of gamObject objects); contains nuObj as last element

    ## return
    res
}


### gamGPDboot() ###############################################################

##' @title Post-blackend Bootstrap of Chavez-Demoulin and Davison (2005) for gamGPDfit()
##' @param x see gamGPDfit()
##' @param B number of bootstrap replications
##' @param threshold see gamGPDfit()
##' @param nextremes see gamGPDfit()
##' @param datvar see gamGPDfit()
##' @param etaFrhs see gamGPDfit()
##' @param nuFrhs see gamGPDfit()
##' @param init see gamGPDfit()
##' @param niter see gamGPDfit()
##' @param include.updates see gamGPDfit()
##' @param eps.eta see gamGPDfit()
##' @param eps.nu see gamGPDfit()
##' @param boot.progress logical indicating whether progress information is displayed
##' @param progress see gamGPDfit() (only used if progress==TRUE)
##' @param adjust logical indicating whether non-real values of the derivatives are adjusted
##' @param verbose see gamGPDfit() (only used if progress==TRUE)
##' @param debug logical indicating whether initial fit is saved
##' @param ... see gamGPDfit()
##' @return a list of length B+1, the first component being the fitted object
##'         as returned by gamGPDfit(); the other components contain similar
##'         objects based on the B bootstrap replications, but can be list() in case
##'         gamGPDfit() reports it (see there)
##' @author Marius Hofert
gamGPDboot <- function(x, B, threshold, nextremes=NULL, datvar, etaFrhs, nuFrhs,
                       init=fit.GPD(x[,datvar], threshold=threshold, type="pwm", verbose=FALSE)$par.ests,
                       niter=32, include.updates=FALSE, eps.eta=1e-5, eps.nu=1e-5,
                       boot.progress=TRUE, progress=FALSE, adjust=TRUE, verbose=FALSE,
                       debug=FALSE, ...)
{
    ## progress
    if(boot.progress) {
         if(progress) cat("\nStarting initial fit:\n") else {
             pb <- txtProgressBar(max=B+1, style=if(isatty(stdout())) 3 else 1) # setup progress bar
             on.exit(close(pb)) # close progress bar
        }
    }

    ## (major) fit using gamGPDfit()
    fit <- gamGPDfit(x=x, threshold=threshold, nextremes=nextremes, datvar=datvar,
                     etaFrhs=etaFrhs, nuFrhs=nuFrhs, init=init, niter=niter,
                     include.updates=include.updates, eps.eta=eps.eta, eps.nu=eps.nu,
                     progress=if(!boot.progress) FALSE else progress, adjust=adjust,
                     verbose=if(!boot.progress) FALSE else verbose, ...)
    ## => can be list() (if gam() in gamGPDfitUp() failed)
    if(!length(fit))
        stop("gamGPDfit() reported that gam() in gamGPDfitUp() failed in the major fit (already)")

    ## progress
    if(boot.progress) if(progress) cat("\n") else setTxtProgressBar(pb, 1)

    ## pick out fitted values
    eta <- fit$eta # fitted eta
    nu <- fit$nu # fitted nu
    xi <- expm1(eta) # fitted xi
    beta <- exp(nu)/(1+xi) # fitted beta

    ## for debugging
    if(debug) save(fit, file="gamGPDboot_debug.rda")

    ## post-blackened bootstrap; see Chavez-Demoulin and Davison (2005)
    rnum <- 1 # iteration number
    bfit <- lapply(1:B, function(b){
        ## resample residuals within each group of (same) covariates
        ## => in particular, the *number* of excesses remains the same
        ##    for each covariate combination
        rr <- ave(fit$res, fit$covar,
                  FUN=function(r) if(length(r)==1) r else
                        sample(r, size=length(r), replace=TRUE))

        ## compute corresponding excesses
        y. <- mapply(qGPD, p=-expm1(-rr), xi=xi, beta=beta) # reconstruct excesses from residuals
        x. <- data.frame(fit$covar, y=y.) # add excesses/covariates

        ## progress
        if(boot.progress && progress) cat("Starting fit in bootstrap run ",
                                          b, " of ", B, ":\n", sep="")

        ## call gamGPDfit()
        ## note: threshold=0 (since the data is x.!)
        ##       => we discard those excesses which are equal to 0
        bfitobj <- gamGPDfit(x=x., threshold=0, nextremes=nextremes, datvar="y",
                             etaFrhs=etaFrhs, nuFrhs=nuFrhs,
                             init=fit.GPD(x.[,"y"], threshold=0, type="pwm", verbose=FALSE)$par.ests,
                             niter=niter, include.updates=include.updates,
                             eps.eta=eps.eta, eps.nu=eps.nu,
                             progress=if(!boot.progress) FALSE else progress,
                             adjust=adjust,
                             verbose=if(!boot.progress) FALSE else verbose, ...)
        ## => can be list() (if gam() in gamGPDfitUp() failed)

        ## progress
        if(boot.progress) if(progress) cat("\n") else setTxtProgressBar(pb, b+1)

        ## return fit
        bfitobj
    })

    ## return list of length B+1 containing all the gamGPDfit() objects
    c(fit=list(fit), bfit=bfit)
}


### Functions to extract fitted results and for predicting #####################

##' @title Extract Fits in a minimalized form
##' @param x object as returned by gam()
##' @return a list with components
##'         covar: containing the ('minimalized') covariate combinations
##'         fit:   corresponding fitted values
##' @author Marius Hofert
get.gam.fit <- function(x)
{
    ## check x
    stopifnot(inherits(x, "gam"))

    ## determine the minimal grid which contains all combinations of covars used for fitting the model
    covar.nms <- if(length(x$var.summary)>0) names(x$var.summary) else "1" # names of covariates used for fitting (right-hand side of the formula in gam()); NULL if no covariates are used
    covars <- if(length(x$var.summary)>0) x$model[,covar.nms, drop=FALSE] else NULL # build 'minimal' grid of covariates used for fitting

    ## determine fitted values for each covariate combination in (the 'minimalized') covars
    ## => Only combinations of the covariates used for fitting are given and thus
    ##    there are no fitted values returned for each excess.
    y <- cbind(covars, value=x$fitted.values) # value = lambda or rho
    frml <- as.formula(paste("value ~", paste(rev(covar.nms), collapse=" + ")))
    covar.hat <- aggregate(frml, data=y, function(z) z[1]) # pick out only first value (they are equal anyways)
    covar.hat <- covar.hat[,rev(names(covar.hat)), drop=FALSE] # revert again to keep column order [now value is in the *first* column]

    ## return
    list(covar = if(ncol(covar.hat) > 1) covar.hat[,-1] else NULL, # covariate combinations used for fitting
         fit   = covar.hat[,1]) # fitted value
}

##' @title Compute Predicted lambda or rho and Pointwise Asymptotic Two-Sided 1-alpha Confidence Intervals
##' @param x object as returned by gam()
##' @param newdata 'newdata' object as required by predict(); named data.frame
##'        of type expand.grid(covar1=, covar2=) with at least the covariates
##'        used for fitting with gam(); if more are provided, predict() returns
##'        values which are equal uniformly over all of these additional
##'        covariates. Each covariate which appears when fitting with gam() can
##'        have more values than were actually used in gam() (for example,
##'        half-years). In this case predict() 'interpolates' correctly with the
##'        fitted model.
##' @param alpha significance level
##' @return a list with components
##'         covar:   containing the covariate combinations as provided by newdata
##'         predict: the predicted lambda or rho
##'         CI.low:  lower CI [based on *predicted* values]
##'         CI.up:   upper CI [based on *predicted* values]
##' @author Marius Hofert
gam.predict <- function(x, newdata=NULL, alpha=0.05, value=c("lambda", "rho"))
{
    ## check x
    stopifnot(inherits(x, "gam"))
    value <- match.arg(value)

    ## default for newdata
    if(is.null(newdata)) { # choose useful default (exactly the covariates used for fitting but *all* combis of such)
	if(length(x$var.summary) > 0) {
            covar.nms <- names(x$var.summary) # names of covariates used for fitting (right-hand side of the formula in gam())
            ## determine the minimal grid which contains all combinations of covars used for fitting the model
            covars <- x$model[,covar.nms, drop=FALSE] # build 'minimal' grid of covariates used for fitting
            covars. <- lapply(1:ncol(covars), function(j) sort(unique(covars[,j]))) # determine levels
            names(covars.) <- colnames(covars) # put in names
            newdata <- expand.grid(rev(covars.))[,rev(seq_len(length(covars.))), drop=FALSE] # expand grid with reverted covariates (and reverting back) to guarantee the same sorting order of rows as covars
	## check newdata
            if(!all(covar.nms %in% colnames(newdata)))
                stop("'newdata' requires at least the covariates used for fitting ", value, " via gam()")
        } else {
            newdata <- data.frame(dummy.covar=1)
        }
    }

    ## predict (note: can contain half-years etc.)
    pred <- predict(x, newdata=newdata, se.fit=TRUE) # predict object
    pred. <- as.numeric(pred$fit) # predicted values (lambda/rho)
    se. <- as.numeric(pred$se.fit) # standard error
    qa2 <- qnorm(1-alpha/2) # 1-alpha/2 quantile of N(0,1)

    ## return (exp() is due to family=poisson used for lambda;
    ##         rho comes from a logistic regression => not needed there)
    list(covar   = if(length(x$var.summary) > 0) newdata else NULL, # covariate combinations as specified by newdata
         predict = if(value=="lambda") exp(pred.) else pred., # predicted lambda/rho
         CI.low  = if(value=="lambda") exp(pred.-qa2*se.) else pred.-qa2*se., # lower CI [based on *predicted* values]
         CI.up   = if(value=="lambda") exp(pred.+qa2*se.) else pred.+qa2*se.) # upper CI [based on *predicted* values]
}

##' @title Compute Fitted GPD Parameters xi and beta and Bootstrapped Pointwise Two-Sided
##'        1-alpha Confidence Intervals
##' @param x object as returned by gamGPDboot()
##' @param alpha significance level
##' @return a list with components
##'         xi:   a list with components
##'               covar:  a data.frame (possibly data.frame()) containing the 'minimal'
##'                       covariate combinations for the covariates used for fitting xi
##'               fit:    corresponding fitted xi's
##'               CI.low: corresponding lower CIs
##'               CI.up:  corresponding upper CIs
##'               boot:   a matrix containing the corresponding bootstrapped xi's
##'         beta: same as xi, just for beta.
##' @author Marius Hofert
##' Note: Standard errors as for lambda would be available via predict(..., se.fit=TRUE),
##'       but only for xi and nu, not for beta. That's why we need bootstrapped values
##'       here (and thus only have CIs for the fitted values)
get.GPD.fit <- function(x, alpha=0.05)
{
    if(length(x[[1]])==0)
        stop("length(x[[1]]) must be > 0; otherwise, already the fitting step in gamGPDboot() failed")
    ## => x[[1]] serves as a role model (for the minimized grid of covariate combinations)
    x. <- x[sapply(x, length) > 0] # pick out non-empty (= non-failed gam()) sublists
    l. <- length(x.)
    if(l. == 0) stop("no non-failed fitted object available")
    l <- length(x)
    if(l. != l)
        warning("found ", l-l.," list() in x (failed fits); those (resampled) cases will be discarded")
    ## => guarantees that xi.mat. etc. are indeed matrices

    ## Note: 1) The matrices xi.mat. and beta.mat. below typically have a huge number of rows nr.
    ##          For each combination of covariates, they contain the same fitted value for
    ##          each excess. This is 'overhead' we don't want. We therefore now pick
    ##          out the fitted values for each unique combination of covariates which was
    ##          originally used for fitting.
    ##       2) We first consider beta, since that consists of the largest number of covars
    ##          (by having those of xi *and* nu)

    ## determine the minimal grid which only contains all *available* combinations of
    ## covars used for fitting xi and nu (and thus *beta*) [in particular not for each excess]

    ## x.[[1]] is the role model; this is the 'long' version (covariate combination for each excesses)
    covars. <- x.[[1]]$covar
    if(length(covars.) == 0) { # neither xi nor nu (thus beta) depend on covariates
        covar.index <- 1
        covars <- data.frame()
    } else { # at least one of xi or nu has covariates => minimize to unique covariate combinations
        covar.nms <- names(covars.) # names of covariates used for fitting xi and nu (thus beta)
        y <- cbind(x.[[1]]$covar, index=seq_len(nrow(covars.))) # append index ('long' version)
        frml <- as.formula(paste("index ~", paste(rev(covar.nms), collapse=" + "))) # formula [with reverted covariates to guarantee the same sorting order of rows (cols are then reverted but we don't care here)]
        covar.index <- aggregate(frml, data=y, FUN=function(z) z[1])[,"index"] # pick out only first value (they are equal anyways)
        covars <- y[covar.index, covar.nms, drop=FALSE] # minimal version
        rownames(covars) <- NULL
    }

    ## now consider the eta covariates (and thus those of xi)
    xi.covars. <- x.[[1]]$eta.covar
    if(length(xi.covars.) == 0) {
        xi.covar.index <- 1
        xi.covars <- data.frame()
    } else { # xi has at least one covariate => minimize to unique covariate combinations
        xi.covar.nms <- names(xi.covars.) # names of covariates used for fitting xi/eta
        xi.frml <- as.formula(paste("index ~", paste(rev(xi.covar.nms), collapse=" + ")))
        xi.covar.index <- aggregate(xi.frml, data=y, FUN=function(z) z[1])[,"index"] # pick out only first value (they are equal anyways)
        xi.covars <- y[xi.covar.index, xi.covar.nms, drop=FALSE] # minimal version
        rownames(xi.covars) <- NULL
    }

    ## determine fitted values for each covariate combination in (the 'minimalized') covars
    ## pick out all (max. B+1) vectors of estimated xi/beta; (nr, <= B+1) matrix
    xi.mat. <- sapply(x., `[[`, "xi")
    xi.mat <- xi.mat.[xi.covar.index, , drop=FALSE]
    rownames(xi.mat) <- NULL; colnames(xi.mat) <- NULL
    beta.mat. <- sapply(x., `[[`, "beta")
    beta.mat <- beta.mat.[covar.index, , drop=FALSE]
    rownames(beta.mat) <- NULL; colnames(beta.mat) <- NULL

    ## compute CIs (derived from fitted values + bootstrapped ones)
    ## 'na.rm=TRUE' is for those covariate combinations where there is no data
    ## (=> the whole row is NA)
    xi.CI   <- t(apply(xi.mat, 1, quantile, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE, names=FALSE)) # CIs (low, up) for xi; (nrow(xi.mat), 2) matrix
    beta.CI <- t(apply(beta.mat, 1, quantile, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE, names=FALSE)) # CIs (low, up) for beta; (nrow(beta.mat), 2) matrix

    ## result (note: beta depends on xi and nu => on all covariates)
    list(xi   = list(covar  = xi.covars, # covariate combinations used for fitting (or data.frame())
                     fit    = xi.mat[,1], # fitted xi
                     CI.low = xi.CI[,1], # lower CI
                     CI.up  = xi.CI[,2], # upper CI
                     boot   = if(ncol(xi.mat) > 1) xi.mat[,-1] else NULL), # bootstrapped xi's
         beta = list(covar  = covars, # covariate combinations used for fitting (or data.frame())
                     fit    = beta.mat[,1], # fitted beta
                     CI.low = beta.CI[,1], # lower CI
                     CI.up  = beta.CI[,2], # upper CI
                     boot   = if(ncol(beta.mat) > 1) beta.mat[,-1] else NULL)) # bootstrapped beta's
}

##' @title Compute Predicted GPD Parameters xi and beta
##' @param x object as returned by gamGPDboot()
##' @param xi.newdata 'newdata' object for xi as required by predict(); named data.frame
##'        of type expand.grid(covar1=, covar2=) with at least the covariates
##'        used for fitting xi with gam(); if more are provided, predict() returns
##'        values which are equal uniformly over all of these additional
##'        covariates. Each covariate which appears when fitting with gam() can
##'        have more values than were actually used in gam() (for example,
##'        half-years). In this case predict() 'interpolates' correctly with the
##'        fitted model.
##' @param beta.newdata similar to xi.newdata. Since beta is estimated based on xi
##'        (and nu), beta.newdata must contain at least the covariates of xi.newdata.
##' @return a list with components
##'         xi:   a list with components
##'               covar:   a data.frame containing the covariate combinations as provided by xi.newdata
##'               predict: the predicted xi's
##'         beta: same as xi, just for beta; predicted values are based on beta.newdata
##' @author Marius Hofert
GPD.predict <- function(x, xi.newdata = NULL, beta.newdata = NULL)
{
    ## default for xi.newdata and beta.newdata
    if(is.null(xi.newdata)) { # choose useful default (covariates used for fitting but *all* combis of such)
        xi.newdata <- if(length(x[[1]]$eta.covar)>0)
            expand.grid(rev(x[[1]]$eta.covar))[,rev(seq_len(length(x[[1]]$eta.covar))), drop=FALSE]
        else data.frame(xi.dummy.covar=1)
    }
    if(is.null(beta.newdata)) {
        fulllist <- c(x[[1]]$eta.covar, x[[1]]$nu.covar) # some may be duplicate
        sublist <- fulllist[unique(names(fulllist))] # pick out maximal unique subset
        beta.newdata <- if(length(x[[1]]$eta.covar)>0 || length(x[[1]]$nu.covar)>0)
            expand.grid(rev(sublist))[,rev(seq_len(length(sublist))), drop=FALSE] # build grid
        else data.frame(beta.dummy.covar=1)
    }

    ## check xi.newdata and beta.newdata (true for the default)
    if(!all(names(x[[1]]$eta.covar) %in% colnames(xi.newdata)))
        stop("'xi.newdata' requires at least the covariates used for fitting xi via gamGPDfit()")
    if(!all(names(x[[1]]$covar) %in% colnames(beta.newdata)))
       stop("'beta.newdata' requires at least the covariates used for fitting beta via gamGPDfit()")
    ## => implicitly also checks that beta.newdata contains the covariates of both eta and nu

    ## predict
    ## note: - *.newdata can be of different length than fitted values
    ##         (can contain half-years etc.)
    ##       - eta.newdata and beta.newdata can be of different length
    ##         (depending on the covariates used for fitting, for example)
    eta.pred      <- as.numeric(predict(x[[1]]$etaObj, newdata=xi.newdata)) # predict eta on xi.newdata
    eta.pred.beta <- as.numeric(predict(x[[1]]$etaObj, newdata=beta.newdata)) # predict eta on beta.newdata
    nu.pred.beta <- as.numeric(predict(x[[1]]$nuObj, newdata=beta.newdata)) # predict nu on beta.newdata
    ## note: - there is no betaObj so we can't predict beta directly (only through nu)
    ##       - we predict eta and nu both on beta.newdata since that's what we need
    ##         to predict beta (see below)

    ## result
    list(xi   = list(covar   = if(length(x[[1]]$eta.covar)>0) xi.newdata else NULL, # covariate combinations as specified by newdata
                     predict = expm1(eta.pred)), # predicted xi
         beta = list(covar   = if(length(x[[1]]$eta.covar)>0 || length(x[[1]]$nu.covar)>0)
                                   beta.newdata else NULL, # covariate combinations as specified by newdata
                     predict = exp(nu.pred.beta - eta.pred.beta))) # predicted beta
}


### Computing risk measures ####################################################

##' @title Compute Value-at-Risk or Expected Shortfall for a GPD
##' @param x matrix with three columns containing rho (estimate for \bar{F}(u)
##'        depending on covariates; obtained via logistic regression),
##'        corresponding xi and beta
##' @param alpha confidence level
##' @param u threshold
##' @param method either "VaR" for Value-at-Risk or "ES" for expected shortfall
##' @return Value-at-Risk or expected shortfall
risk.measure <- function(x, alpha, u, method = c("VaR", "ES"))
{
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    stopifnot(ncol(x)==3)
    rho <- x[,1] # rho (estimate for \bar{F}(u) depending on covariates)
    xi <- x[,2] # corresponding xi
    beta <- x[,3] # corresponding beta
    method <- match.arg(method)
    switch(method,
           "VaR"={
               ## The theory:
               ## The number of exceedances until t is Poi(Lambda(t)) distributed,
               ## hence N_t ~ Poi(Lambda(t)); in our case N_{x, t} ~
               ## Poi(Lambda(x, t)). The number of exceedances in [t, t+1]
               ## (1y here) is thus Poi(Lambda(x, [t,t+1])) distributed, with
               ## Lambda(x, [t,t+1]) = \int_t^{t+1} lambda(x, s) ds. Since our
               ## lambda(x, .) is constant in [t,t+1] (only depending on
               ## covariates x other than time), we obtain Lambda(x, [t,t+1]) =
               ## lambda(x, t) * (t+1 - t) = lambda(x, t). So the number of
               ## exceedances in [t,t+1] is Poi(lambda(x, t)) distributed and
               ## thus the expected number of exceedances in [t,t+1] is
               ## lambda(x, t); we do this additional approximation of
               ## considering the *expected* number to bring lambda into play
               ## and thus obtain a parametric estimators (depending on covariates).
               ## The denominator in the formula for VaR has to be
               ## \bar{F}(u) ~= N_u/n (s. QRM book, p. 283) ~= E[N_u]/n (our
               ## additional approximation) ~= lambda(x, t)/n_{x, t}, where
               ## n_{x, t} is the number of overall losses for covariate x in
               ## [t, t+1].
               ##
               ## The problems:
               ## 1) lambda(x, t)/n_{x, t} does not have to be in [0,1]
               ##    (but it is for large enough thresholds u; if this is not
               ##     the case (lambda [= expected number of exceedances] > n_{x, t}
               ##     [= total number of events]), this clearly indicates that u
               ##     needs to be larger).
               ## 2) If t is in the future, how would we predict n_{x, t}?
               ##
               ## The/A solution:
               ## We directly approximate \bar{F}(u) by a *logistic* regression
               ## => predict() available
               u + (beta/xi) * (( (1-alpha) / rho )^(-xi) - 1)
           },
           "ES"={
               ES <- (risk.measure(x, alpha=alpha, u=u, method="VaR")+beta-xi*u)/(1-xi)
               ## adjust to be Inf if xi >= 1 (i.e., ES < 0)
               ## that's a convention, see p. 79 Coles (2001)
               if(any(xi >= 1)) ES[xi >= 1] <- Inf
               ES
           },
           stop("wrong method"))
}
