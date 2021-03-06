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
## Remark:
## 1) Notation:
##    paper | Valerie's thesis | Valerie's former code | Coles' book + ismev
##      xi        -kappa                kappa                    xi
##    beta         sigma                sigma                 sigma
##      nu            nu                  eta                    --
## 2) GPD(xi in IR, beta > 0) distribution function (same as in ismev; up to notation):
##    G_{xi,beta}(x) = 1-(1+xi*x/beta)^(-1/xi) if xi!=0
##                   = 1-exp(-x/beta)          if xi =0
##    x>=0 when xi>=0 and x in [0,-beta/xi] when xi<0
## 3) GPD(xi, beta) density for xi>0:
##    g_{xi,beta}(x) = (1+xi*x/beta)^(-(1+1/xi))/beta if x>0


### Auxiliary functions ########################################################

##' Reparameterized log-likelihood l^r(xi, nu; y_1,..,y_n) = l(xi, exp(nu)/(1+xi);
##' y_1,..,y_n), where l(xi, beta; y_1,..,y_n) denotes the log-likelihood of a
##' GPD(xi, beta) distribution
##'
##' @title Reparameterized log-Likelihood l^r
##' @param y vector of excesses (over a high threshold u)
##' @param xi GPD(xi,beta) vector of parameters xi
##' @param nu GPD(xi, beta) vector of parameters nu (orthogonal in the Fisher
##'        information metric to xi)
##' @return reparametrized log-likelihood l^r
##' @author Marius Hofert
rlogL <- function(y, xi, nu)
{
    ## checks
    stopifnot((n <- length(y)) > 0, length(xi)==n, length(nu)==n)
    if(any(ii <- xi <= -1)) { # not a valid reparameterization (but we make it one)
        perc <- 100*sum(ii)/length(ii) # percentage of xi <= -1
        stop(round(perc,2), "% of all xi are <= -1; not a valid reparameterization since corresponding beta are not > 0") # beta = exp(nu)/(1+xi)
    }

    ## set up result
    res <- rep(-Inf, n)

    ## main
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

##' @title Compute Derivatives of the Reparameterized log-Likelihood
##' @param y vector of excesses (over a high threshold u)
##' @param xi vector of GPD(xi,beta) parameters xi
##' @param nu vector of (orthogonal in the Fisher information metric)
##'        GPD(xi, beta) parameters nu (> 0)
##' @param adjust logical indicating whether non-real values of the derivatives are adjusted
##' @param verbose logical indicating whether modified arguments are printed
##' @return (n x 4) matrix containing the partial derivatives of the
##'         reparameterized log-likelihood l^r where
##'         column 1: derivative of the reparameterized log-likelihood w.r.t. xi
##'         column 2: derivative of the reparameterized log-likelihood w.r.t. nu
##'         column 3: 2nd derivative of the reparameterized log-likelihood w.r.t. xi
##'         column 4: 2nd derivative of the reparameterized log-likelihood w.r.t. nu
##' @author Marius Hofert
##' Note: Column 3 and 4 have different sign than in the old code
DrlogL <- function(y, xi, nu, adjust=TRUE, verbose=TRUE)
{
    ## checks
    stopifnot((n <- length(y)) > 0, length(xi)==n, length(nu)==n)
    if(any(ii <- xi <= -1)) { # not a valid reparameterization (but we make it one)
        perc <- 100*sum(ii)/length(ii) # percentage of xi <= -1
        if(verbose) warning(round(perc,2), "% of all xi are <= -1; truncated at -1 (=> beta = Inf)") # beta = exp(nu)/(1+xi)
        xi <- pmax(xi, -1)
    }

    ## set up result
    res <- matrix(0, nrow=n, ncol=4,
                  dimnames=list(NULL, c("rl.xi", "rl.nu", "rl.xixi", "rl.nunu")))

    ## main
    ii <- (xi < 0 & (0 <= y & y < -exp(nu)/(xi*(1+xi)))) | (xi > 0 & y >= 0)
    if(any(ii)) {
        ## auxiliary terms
        xi. <- xi[ii]
        nu. <- nu[ii]
        y.  <-  y[ii]
        a <- 1+xi.*(1+xi.)*exp(-nu.)*y.
        la <- log1p(xi.*(1+xi.)*exp(-nu.)*y.)
        a. <- exp(-nu.)*y.*(1+2*xi.)
        ## main
        res[ii,"rl.xi"] <- 1/(1+xi.) + la/xi.^2 - (1+1/xi.)*a./a
        res[ii,"rl.xixi"] <- -1/(1+xi.)^2 - 2*la/xi.^3 + 2*a./(a*xi.^2) -
            (1+1/xi.)*exp(-nu.)*y. * (2*a - (1+2*xi.)^2 * exp(-nu.)*y.) / a^2
        res[ii,"rl.nu"] <- (-1+(1+xi.)*exp(-nu.)*y.) / a
        res[ii,"rl.nunu"] <- -(1+xi.)^2*y.*exp(-nu.)/a^2
    }
    ii <- xi == 0
    if(any(ii)) {
        ## auxiliary terms
        ynu <- y[ii]*exp(-nu[ii])
        res[ii,"rl.xi"] <- 0
        res[ii,"rl.xixi"] <- 0
        res[ii,"rl.nu"] <- -1+ynu
        res[ii,"rl.nunu"] <- -ynu
    }

    ## replace non-finite derivatives by mean of all finite ones
    if(adjust) {
        res[,"rl.xi"]   <- adjustD(res[,"rl.xi"],   order=1, verbose=verbose)
        res[,"rl.xixi"] <- adjustD(res[,"rl.xixi"], order=2, verbose=verbose)
        res[,"rl.nu"]   <- adjustD(res[,"rl.nu"],   order=1, verbose=verbose)
        res[,"rl.nunu"] <- adjustD(res[,"rl.nunu"], order=2, verbose=verbose)
    }

    ## return
    res
}

##' @title Function to Adjust Derivatives
##' @param x vector of values (derivatives)
##' @param order order of the derivatives to be adjusted
##' @param verbose logical indicating whether warnings about adjustments of
##'        the derivatives are printed
##' @return adjusted derivatives
##' @author Marius Hofert
##' Note: This is an auxiliary function of DrlogL()
adjustD <- function(x, order, verbose=TRUE)
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

##' @title Compute Update (one Iteration) in gamGPDfit()
##' @param y data.frame containing the excesses over the threshold in a column
##'        labeled yname
##' @param xi.nu 2-column matrix of GPD parameters (xi,nu) to be updated where
##'        nu is orthogonal to xi in the Fisher information metric
##' @param xiFrhs right-hand side of the formula for xi in the gam() call
##'        for fitting xi
##' @param nuFrhs right-hand side of the formula for nu in the gam() call
##'        for fitting nu
##' @param yname string containing the name of the column of y which contains
##'        the excesses
##' @param adjust logical indicating whether non-real values of the derivatives are adjusted
##' @param verbose logical indicating whether warnings about adjustments of
##'        the derivatives and wrong arguments in DrlogL() are printed
##' @param ... additional arguments passed to gam()
##' @return a list of length four containing
##'         element 1 (xi): object of class gamObject for xi as returned by mgcv::gam()
##'         element 2 (nu): object of class gamObject for nu as returned by mgcv::gam()
##'         element 3 (xi.weights): weights associated with xi
##'         element 4 (nu.weights): weights associated with nu
##' @author Marius Hofert
##' Note: That's a helper function of gamGPDfit()
gamGPDfitUp <- function(y, xi.nu, xiFrhs, nuFrhs, yname, adjust=TRUE, verbose=TRUE, ...)
{
    stopifnot(is.data.frame(y), (dim. <- dim(y))[2] >= 1,
              length(which(colnames(y)==yname))==1,
              (n <- dim.[1]) > 0, dim(xi.nu)==c(n, 2))

    ## pick out xi and nu (for readability)
    xi <- xi.nu[,1]
    nu <- xi.nu[,2]

    ## compute one Newton step in xi
    DrLL <- DrlogL(y[,yname], xi=xi, nu=nu, adjust=adjust, verbose=verbose) # (n1,4) matrix
    rl.xi <- DrLL[,"rl.xi"] # score in xi
    rl.xixi. <- DrLL[,"rl.xixi"] # -weight
    Newton.xi <- xi - rl.xi / rl.xixi. # Newton step

    ## concatenate Newton.xi and rl.xixi. to y, build formula, and estimate xi
    if("Newton.xi" %in% colnames(y))
        stop("y is not allowed to have a column named 'Newton.xi'")
    y. <- cbind(y, Newton.xi=Newton.xi, rl.xixi.=rl.xixi.)
    xi.formula <- update(xiFrhs, Newton.xi~.) # build formula Newton.xi ~ xiFrhs
    xi.obj <- gam(xi.formula, data=y., weights=-rl.xixi., ...) # updated xi object of type gamObject

    ## build fitted (xi) object and check
    xi.fit <- fitted(xi.obj)
    if((n. <- length(xi.fit)) != n) stop("After introducing adjustD(), this error should not appear anymore")

    ## compute one Newton step in nu (for given new xi)
    DrLL <- DrlogL(y[,yname], xi=xi.fit, nu=nu, adjust=adjust, verbose=verbose) # (n1,4) matrix
    rl.nu <- DrLL[,"rl.nu"] # score in nu
    rl.nunu. <- DrLL[,"rl.nunu"] # -weight
    Newton.nu <- nu - rl.nu / rl.nunu. # Newton step

    ## concatenate Newton.nu and rl.nunu. to y, build formula, and estimate nu
    if("Newton.nu" %in% colnames(y))
        stop("y is not allowed to have a column named 'Newton.nu'")
    y. <- cbind(y, Newton.nu=Newton.nu, rl.nunu.=rl.nunu.)
    nu.formula <- update(nuFrhs, Newton.nu~.) # build formula Newton.nu ~ nuFrhs
    nu.obj <- gam(nu.formula, data=y., weights=-rl.nunu., ...) # updated nu object of type gamObject

    ## return list of two gamObject objects (for xi, nu)
    list(xi=xi.obj, nu=nu.obj, xi.weights=-rl.xixi., nu.weights=-rl.nunu.) # note: the naming (xi, nu) is not ideal but guarantees that colnames(param.old) = colnames(param.new) in gamGPDfit
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
##' @param xiFrhs right-hand side of the formula for xi in the gam() call
##'        for fitting xi
##' @param nuFrhs right-hand side of the formula for nu in the gam() call
##'        for fitting nu
##' @param init bivariate vector containing initial values for (xi, beta)
##' @param niter maximal number of iterations in the backfitting algorithm
##' @param include.updates logical indicating whether updates for xi and nu are
##'        returned as well
##' @param epsxi epsilon for stop criterion for xi
##' @param epsnu epsilon for stop criterion for nu
##' @param progress logical indicating whether progress information is displayed
##' @param adjust logical indicating whether non-real values of the derivatives are adjusted
##' @param verbose logical passed to gamGPDfitUp() (thus to DrlogL() and
##'        adjustD())
##' @param ... additional arguments passed to gam() (called by gamGPDfitUp())
##' @return a list; see below
##' @author Marius Hofert
gamGPDfit <- function(x, threshold, nextremes = NULL, datvar, xiFrhs, nuFrhs,
                      init = fit.GPD(x[,datvar], threshold=threshold, type="pwm", verbose=FALSE)$par.ests,
                      niter = 32, include.updates = FALSE, epsxi = 1e-5, epsnu = 1e-5,
                      progress = TRUE, adjust = TRUE, verbose = FALSE, ...)
{
    ## checks
    stopifnot(is.data.frame(x), length(init)==2, niter>=1, epsxi>0, epsnu>0)
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

    ## initial values for xi, beta, and nu (nu = reparameterization of beta)
    xi.init <- init[1]
    beta.init <- init[2]
    nu.init <- log((1+xi.init) * beta.init)

    ## iteration ###############################################################

    iter <- 1 # iteration number
    updates <- list() # (empty) list of update (gamGPDfitUp) objects
    while(TRUE){

        ## update/fit parameters (xi, nu)
        param.old <- if(iter==1){
            matrix(rep(c(xi.init, nu.init), each=n.ex), ncol=2,
                       dimnames=list(rownames(y.), c("xi", "nu"))) # (n.ex,2)-matrix with cols "xi" and "nu"
        } else{
            param.new
        }
        updates[[iter]] <- gamGPDfitUp(y., xi.nu=param.old,
                                       xiFrhs=xiFrhs, nuFrhs=nuFrhs,
                                       yname=datvar, adjust=adjust, verbose=verbose, ...) # returns a list of two gam() objects containg the fitted (xi, nu)
        param.new <- sapply(updates[[iter]][c("xi", "nu")], fitted)
        ## note: param.old and param.new have the same rownames/colnames

        ## check
        if(any(dim(param.new)!=dim(param.old))) stop("gamGPDfitUp() returned an updated gamObject object of wrong dimension")
        if(progress){
            conv <- sapply(updates[[iter]][c("xi", "nu")], function(x) x$converged) # convergence status for (xi, nu)
            if(any(!conv)) warning("gam() in gamGPDfitUp() did not converge for ", paste(c("xi","nu")[!conv], collapse=", "))
        }

        ## tracing
        MRD.. <- colMeans(abs((param.old-param.new)/param.old)) # mean relative distance for (xi, nu)
        MRD. <- c(iter=iter, MRD..[1], MRD..[2])
        MRD <- if(iter==1) MRD. else rbind(MRD, MRD., deparse.level=0)
        if(progress){
            cat("Mean relative differences in iteration ", iter,
                " (xi, nu): (", sprintf("%g", MRD..[1]), ", ",
                sprintf("%g", MRD..[2]), ")\n", sep="")
        }

        ## check for "convergence"
        conv <- c(MRD..[1] <= epsxi, MRD..[2] <= epsnu)
        if(all(conv)) break
        if(iter >= niter){
            if(progress) warning("Reached 'niter' without the required precision for ",
                                 paste(c("xi","nu")[!conv], collapse=", "))
            break
        }
        iter <- iter + 1 # update iteration

    }

    ## finish and return #######################################################

    ## fitted parameters
    xi <- param.new[,"xi"] # fitted xi's
    nu <- param.new[,"nu"] # fitted nu's
    beta <- exp(nu)/(1+xi) # corresponding (fitted) beta

    ## check beta for being a valid reparameterization (otherwise pGPD() fails)
    ## note: nu can be too small (=> beta = 0) or xi <= -1
    ##       => set corresponding beta to NA
    ii <- is.finite(beta) & beta > 0
    beta[!ii] <- NA

    ## compute the residuals
    resi <- mapply(function(x, xi, beta) if(!is.na(beta)) -log1p(-pGPD(x, xi=xi, beta=beta)) else NA,
                   x=y.[,datvar], xi=xi, beta=beta)

    ## log-likelihood
    logL <- rlogL(y=y.[,datvar], xi=xi, nu=nu) # reparameterized log-likelihood (reparameterization doesn't matter, it's the same *value* as the original log-likelihood)

    ## fitted gamObjects for xi and nu and standard errors
    ## (contained in updates but for convenience we treat this additionally)
    xiObj <- updates[[iter]]$xi
    nuObj <- updates[[iter]]$nu
    se.xi <- updates[[iter]]$xi.weights
    se.nu <- updates[[iter]]$nu.weights

    ## determine *all* combinations of the provided covariates
    ## xi
    xi.covars <- names(xiObj$var.summary)
    x.xi <- x[,xi.covars, drop=FALSE] # all cols with covariates used for fitting xi
    all.covars.xi <- lapply(seq_len(ncol(x.xi)), function(j) sort(unique(x.xi[,j]))) # determine levels
    ## => we only can determine those which appear at least in one column
    names(all.covars.xi) <- colnames(x.xi) # put in names
    ## nu
    nu.covars <- names(nuObj$var.summary)
    x.nu <- x[,nu.covars, drop=FALSE] # all cols with covariates used for fitting nu
    all.covars.nu <- lapply(seq_len(ncol(x.nu)), function(j) sort(unique(x.nu[,j]))) # determine levels
    names(all.covars.nu) <- colnames(x.nu) # put in names

    ## build list
    res <- list(xi=xi, # estimated xi
                beta=beta, # estimated beta
                nu=nu, # estimated nu
                se.xi=se.xi, # standard error for xi
                se.nu=se.nu, # standard error for nu
                xi.covar=all.covars.xi, # (unique) covariates for xi
                nu.covar=all.covars.nu, # (unique) covariates for nu
                covar=y.[,union(xi.covars, nu.covars), drop=FALSE], # *available* (not necessarily all) covariate combinations used for fitting beta (= xi *and* nu)
                y=y.[,datvar], # excesses
                res=resi, # residuals
                MRD=MRD, # mean relative distances between old/new (xi, nu) for all iterations
                logL=logL, # log-likelihood at the estimated parameters
                xiObj=xiObj, # gamObject for estimated xi (return object of mgcv::gam())
                nuObj=nuObj) # gamObject for estimated nu (return object of mgcv::gam())
    if(include.updates) res <- c(res, xiUpdates=lapply(updates, `[[`, "xi"), # updates for xi for each iteration (list of gamObject objects); contains xiObj as last element
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
##' @param xiFrhs see gamGPDfit()
##' @param nuFrhs see gamGPDfit()
##' @param init see gamGPDfit()
##' @param niter see gamGPDfit()
##' @param include.updates see gamGPDfit()
##' @param epsxi see gamGPDfit()
##' @param epsnu see gamGPDfit()
##' @param boot.progress logical indicating whether progress information is displayed
##' @param progress see gamGPDfit() (only used if progress==TRUE)
##' @param adjust logical indicating whether non-real values of the derivatives are adjusted
##' @param verbose see gamGPDfit() (only used if progress==TRUE)
##' @param debug logical indicating whether initial fit is saved
##' @param ... see gamGPDfit()
##' @return a list of length B+1, the first component being the fitted object
##'         as returned by gamGPDfit(); the other components contain similar
##'         objects based on the B bootstrap replications.
##' @author Marius Hofert
gamGPDboot <- function(x, B, threshold, nextremes=NULL, datvar, xiFrhs, nuFrhs,
                       init=fit.GPD(x[,datvar], threshold=threshold, type="pwm", verbose=FALSE)$par.ests,
                       niter=32, include.updates=FALSE, epsxi=1e-5, epsnu=1e-5,
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
                     xiFrhs=xiFrhs, nuFrhs=nuFrhs, init=init, niter=niter,
                     include.updates=include.updates, epsxi=epsxi, epsnu=epsnu,
                     progress=if(!boot.progress) FALSE else progress, adjust=adjust,
                     verbose=if(!boot.progress) FALSE else verbose, ...)

    ## progress
    if(boot.progress) if(progress) cat("\n") else setTxtProgressBar(pb, 1)

    ## pick out fitted values
    xi <- fit$xi # fitted xi
    nu <- fit$nu # fitted nu
    beta <- exp(nu)/(1+xi) # fitted beta

    ## for debugging
    if(debug) save(fit, file="gamGPDboot_debug.rda")

    ## post-blackened bootstrap; see Chavez-Demoulin and Davison (2005)
    rnum <- 1 # iteration number
    bfit <- lapply(1:B, function(b){
        ## resample residuals within each group of (same) covariates,
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
                             xiFrhs=xiFrhs, nuFrhs=nuFrhs,
                             init=fit.GPD(x.[,"y"], threshold=0, type="pwm", verbose=FALSE)$par.ests,
                             niter=niter, include.updates=include.updates,
                             epsxi=epsxi, epsnu=epsnu,
                             progress=if(!boot.progress) FALSE else progress,
                             adjust=adjust,
                             verbose=if(!boot.progress) FALSE else verbose, ...)

        ## progress
        if(boot.progress) if(progress) cat("\n") else setTxtProgressBar(pb, b+1)

        ## return fit
        bfitobj
    })

    ## return list of length B+1 containing all the gamGPDfit() objects
    c(fit=list(fit), bfit=bfit)
}


### Functions to extract fitted results and for predicting #####################

##' @title Compute Fitted lambda
##' @param x object as returned by gam()
##' @return a list with components
##'         covar: containing the ('minimalized') covariate combinations
##'         fit:   corresponding fitted values of lambda
##' @author Marius Hofert
get.lambda.fit <- function(x)
{
    ## check x
    stopifnot(inherits(x, "gam"))

    ## determine the minimal grid which contains all combinations of covars used for fitting the model
    covar.nms <- if(length(x$var.summary)>0) names(x$var.summary) else "1" # names of covariates used for fitting (right-hand side of the formula in gam()); NULL if no covariates are used
    covars <- if(length(x$var.summary)>0) x$model[,covar.nms, drop=FALSE] else NULL # build 'minimal' grid of covariates used for fitting

    ## determine fitted values for each covariate combination in (the 'minimalized') covars
    ## => Only combinations of the covariates used for fitting are given and thus
    ##    there are no fitted values returned for each excess.
    y <- cbind(covars, lambda=x$fitted.values)
    frml <- as.formula(paste("lambda ~", paste(rev(covar.nms), collapse=" + ")))
    covar.lam.hat <- aggregate(frml, data=y, function(z) z[1]) # pick out only first value (they are equal anyways)
    covar.lam.hat <- covar.lam.hat[,rev(names(covar.lam.hat)), drop=FALSE] # revert again to keep column order [now lambda is in the *first* column]

    ## return
    list(covar = if(ncol(covar.lam.hat) > 1) covar.lam.hat[,-1] else NULL, # covariate combinations used for fitting
         fit   = covar.lam.hat[,1]) # fitted lambda
}

##' @title Compute Predicted lambda and Pointwise Asymptotic Two-Sided 1-alpha Confidence Intervals
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
##'         predict: the predicted lambda
##'         CI.low:  lower CI [based on *predicted* values]
##'         CI.up:   upper CI [based on *predicted* values]
##' @author Marius Hofert
lambda.predict <- function(x, newdata=NULL, alpha=0.05)
{
    ## check x
    stopifnot(inherits(x, "gam"))

    ## default for newdata
    if(is.null(newdata)) { # choose useful default (exactly the covariates used for fitting but *all* combis of such)
	if(length(x$var.summary) > 0) {
            covar.nms <- names(x$var.summary) # names of covariates used for fitting (right-hand side of the formula in gam())
            ## determine the minimal grid which contains all combinations of covars used for fitting the model
            covars <- x$model[,covar.nms, drop=FALSE] # build 'minimal' grid of covariates used for fitting
            lam.covars <- lapply(1:ncol(covars), function(j) sort(unique(covars[,j]))) # determine levels
            names(lam.covars) <- colnames(covars) # put in names
            newdata <- expand.grid(rev(lam.covars))[,rev(seq_len(length(lam.covars))), drop=FALSE] # expand grid with reverted covariates (and reverting back) to guarantee the same sorting order of rows as covars
	## check newdata
            if(!all(covar.nms %in% colnames(newdata)))
                stop("'newdata' requires at least the covariates used for fitting lambda via gam()")
        } else {
            newdata <- data.frame(lambda.dummy.covar=1)
        }
    }

    ## predict (note: can contain half-years etc.)
    pred <- predict(x, newdata=newdata, se.fit=TRUE) # predict object
    lam.pred <- as.numeric(pred$fit) # predicted values for lambda
    lam.se <- as.numeric(pred$se.fit) # standard error
    qa2 <- qnorm(1-alpha/2) # 1-alpha/2 quantile of N(0,1)

    ## return
    list(covar   = if(length(x$var.summary) > 0) newdata else NULL, # covariate combinations as specified by newdata
         predict = exp(lam.pred), # predicted lambda
         CI.low  = exp(lam.pred-qa2*lam.se), # lower CI [based on *predicted* values]
         CI.up   = exp(lam.pred+qa2*lam.se)) # upper CI [based on *predicted* values]
}

##' @title Compute Fitted GPD Parameters xi and beta and Bootstrapped Pointwise Two-Sided
##'        1-alpha Confidence Intervals
##' @param x object as returned by gamGPDboot()
##' @param alpha significance level
##' @return a list with components
##'         xi:   a list with components
##'               covar:  a data.frame containing the 'minimal' covariate combinations
##'                       for the covariates used for fitting xi
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
    ## basic check (B = number of bootstrap replicates; nr = number of rows of the data set
    ## provided to gamGPDboot())
    xi.mat.   <- sapply(x, `[[`, "xi") # pick out all B+1 vectors of xi; (nr, B+1) matrix
    beta.mat. <- sapply(x, `[[`, "beta") # pick out all B+1 vectors of beta; (nr, B+1) matrix
    stopifnot(dim(xi.mat.)==dim(beta.mat.))

    ## Note: Now these matrices typically have a huge number of rows nr. For each
    ##       combination of covariates, they contain the same fitted value for
    ##       each loss. This is 'overhead' we don't want. We therefore now pick
    ##       out the fitted values for each unique combination of covariates which was
    ##       originally used for fitting.

    ## determine the minimal grid which contains all *available* combinations of
    ## covars used for fitting xi and nu (and beta) but not more
    ## (in particular not for each excess)
    covars. <- x[[1]]$covar # 'long' version (covariate combination for each excesses, too)
    covar.nms <- if(length(covars.)>0) names(covars.) else "1" # names of covariates used for fitting xi and nu (or "1" if not depending on covariates)
    y <- cbind(covars., index=seq_len(nrow(covars.))) # dummy data set ('long' version)
    frml <- as.formula(paste("index ~", paste(rev(covar.nms), collapse=" + "))) # formula [with reverted covariates to guarantee the same sorting order of rows (cols are then reverted but we don't care here)]
    covar.index <- aggregate(frml, data=y, FUN=function(z) z[1])[,"index"] # = 1 if y is list()
    covars <- if(length(covar.index)==1) NULL else y[covar.index, covar.nms, drop=FALSE] # 'minimal' version (or NULL for no covariates)

    ## determine further minimalized version for xi (possibly less combinations than beta)
    xi.covars. <- x[[1]]$xi.covar
    xi.covar.nms <- if(length(xi.covars.)>0) names(xi.covars.) else "1" # names of covariates used for fitting xi (or "1" if not depending on covariates)
    xi.frml <- as.formula(paste("index ~", paste(rev(xi.covar.nms), collapse=" + "))) # formula [see above]
    xi.covar.index <- aggregate(xi.frml, data=y, FUN=function(z) z[1])[,"index"] # pick out only first value (they are equal anyways)
    xi.covars <- if(length(xi.covar.index)==1) NULL else y[xi.covar.index, xi.covar.nms, drop=FALSE] # 'minimal' version (or NULL for no covariates)

    ## determine fitted values for each covariate combination in (the 'minimalized') covars
    xi.mat <- xi.mat.[xi.covar.index,, drop=FALSE] # minimal number of rows
    rownames(xi.mat) <- NULL
    beta.mat <- beta.mat.[covar.index,, drop=FALSE] # minimal number of rows
    rownames(beta.mat) <- NULL
    ## compute CIs (derived from fitted values + bootstrapped ones)
    ## 'na.rm=TRUE' is for those covariate combinations where there is no data (=> the whole row is NA)
    xi.CI   <- t(apply(xi.mat, 1, quantile, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE, names=FALSE)) # CIs (low, up) for xi; (nrow(xi.mat), 2) matrix
    beta.CI <- t(apply(beta.mat, 1, quantile, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE, names=FALSE)) # CIs (low, up) for beta; (nrow(beta.mat), 2) matrix

    ## result
    list(xi   = list(covar  = xi.covars, # covariate combinations used for fitting (or NULL)
                     fit    = xi.mat[,1], # fitted xi
                     CI.low = xi.CI[,1], # lower CI
                     CI.up  = xi.CI[,2], # upper CI
                     boot   = xi.mat[,-1]), # bootstrapped xi's
         beta = list(covar  = covars, # covariate combinations used for fitting (or NULL)
                     fit    = beta.mat[,1], # fitted beta
                     CI.low = beta.CI[,1], # lower CI
                     CI.up  = beta.CI[,2], # upper CI
                     boot   = beta.mat[,-1])) # bootstrapped beta's
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
GPD.predict <- function(x, xi.newdata=NULL, beta.newdata=NULL)
{
    ## default for xi.newdata and beta.newdata
    if(is.null(xi.newdata)) { # choose useful default (covariates used for fitting but *all* combis of such)
        xi.newdata <- if(length(x[[1]]$xi.covar)>0) {
            expand.grid(rev(x[[1]]$xi.covar))[,rev(seq_len(length(x[[1]]$xi.covar))), drop=FALSE]
        } else {
            data.frame(xi.dummy.covar=1)
        }
    }
    if(is.null(beta.newdata)) {
        fulllist <- c(x[[1]]$xi.covar, x[[1]]$nu.covar) # some may be duplicate
        sublist <- fulllist[unique(names(fulllist))] # pick out maximal unique subset
        beta.newdata <- if(length(x[[1]]$xi.covar)>0 || length(x[[1]]$nu.covar)>0) {
            expand.grid(rev(sublist))[,rev(seq_len(length(sublist))), drop=FALSE] # build grid
        } else {
            data.frame(beta.dummy.covar=1)
        }
    }

    ## check xi.newdata and beta.newdata (true for the default)
    if(!all(names(x[[1]]$xi.covar) %in% colnames(xi.newdata)))
        stop("'xi.newdata' requires at least the covariates used for fitting xi via gamGPDfit()")
    if(!all(names(x[[1]]$covar) %in% colnames(beta.newdata)))
       stop("'beta.newdata' requires at least the covariates used for fitting beta via gamGPDfit()")
    ## => implicitly also checks that beta.newdata contains the covariates of both xi and nu

    ## predict
    ## note: - *.newdata can be of different length than fitted values
    ##         (can contain half-years etc.)
    ##       - xi.newdata and beta.newdata can be of different length
    ##         (depending on the covariates used for fitting, for example)
    xi.pred      <- as.numeric(predict(x[[1]]$xiObj, newdata=xi.newdata)) # predict xi (for its own)
    xi.pred.beta <- as.numeric(predict(x[[1]]$xiObj, newdata=beta.newdata)) # predict xi on beta.newdata
    nu.pred.beta <- as.numeric(predict(x[[1]]$nuObj, newdata=beta.newdata)) # predict nu on beta.newdata
    ## note: - there is no betaObj so we can't predict beta directly (only through nu)
    ##       - we predict xi and nu both on beta.newdata since that's what we need
    ##         to predict beta (see below)

    ## result
    list(xi   = list(covar   = if(length(x[[1]]$xi.covar)>0) xi.newdata else NULL, # covariate combinations as specified by newdata
                     predict = xi.pred), # predicted xi
         beta = list(covar   = if(length(x[[1]]$xi.covar)>0 || length(x[[1]]$nu.covar)>0)
                                   beta.newdata else NULL, # covariate combinations as specified by newdata
                     predict = exp(nu.pred.beta)/(1+xi.pred.beta))) # predicted beta
}


### Computing risk measures ####################################################

##' @title Compute Value-at-Risk or Expected Shortfall
##' @param x matrix with three columns containing lambda, xi, and beta
##' @param alpha confidence level
##' @param u threshold
##' @param method either "VaR" for Value-at-Risk or "ES" for expected shortfall
##' @return Value-at-Risk or expected shortfall
risk.measure <- function(x, alpha, u, method=c("VaR", "ES")){
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    stopifnot(ncol(x)==3)
    lambda <- x[,1]
    xi <- x[,2]
    beta <- x[,3]
    method <- match.arg(method)
    switch(method,
           "VaR"={
               u+(beta/xi)*(((1-alpha)/lambda)^(-xi)-1)
           },
           "ES"={
               ES <- (risk.measure(cbind(lambda, xi, beta),
                                   alpha=alpha, u=u, method="VaR")+beta-xi*u)/(1-xi)
               ## adjust to be Inf if xi > 1 (i.e., ES < 0)
               ## that's a convention, see p. 79 Coles (2001)
               if(any(xi > 1)) ES[xi > 1] <- Inf
               ES
           },
           stop("wrong method"))
}
