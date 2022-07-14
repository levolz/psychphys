#' Fit the ternary indecision model to detection or discrimination data from
#' a dual-presentation task
#'
#' @description Further details:  Garc?a-P?rez, M.A. & Alcal?-Quintana, R. (2017).
#' The indecision model of psychophysical performance in
#' dual-presentation tasks: Parameter estimation and
#' comparative analysis of response formats.
#' Frontiers in Psychology, 8:1142. <https://doi.org/10.3389/fpsyg.2017.01142>
#'
#' @param data 7-by-N array with data from the dual-presentaion task. The
#' first row gives the N levels of the test stimulus; rows 2-4 respectively
#' give the count of "first", "undecided", and "second" responses at each
#' level when the test stimulus was presented in the first interval/position;
#' rows 5-7 analogously give the count of "first", "undecided", and "second"
#' responses at each level when the test stimulus was presented in the second
#' interval/position. If \code{format <- "2AFC"} or
#' \code{format <- "equality"}, array data still has the same number of rows
#' but its content is subject to the constraints listed under argument format
#' below.
#'
#' @param standard Level of the standard stimulus. For detection data (i.e.,
#' with a nominally null standard), set \code{standard <- -Inf}; for
#' discrimination data (i.e., with a non-null standard), set standard to the
#' appropriate level, in the same units in which levels of the test stimulus
#' are given in the first row of array data.
#'
#' @param format Out of \code{c("ternary", "2AFC", or "equality")} to indicate
#' that data were collected with a ternary response format (first, second, or
#' undecided) a 2AFC response format (forcing first or second responses by
#' guessing when uncertain), or the same-different response format (both stimuli
#' are subjective equal or they are subjectively different).
#' When \code{format <- "2AFC"}, rows 3 and 6 of array data must be filled with
#' zeros (as no "undecided" responses are given under this format);
#' When \code{format <- "equality"}, rows 2 and 5 of array data must be filled with zeros
#' (so "different" responses are stored as "second" responses).
#'
#' @param alpha_bounds 2-element vector with the lower and upper bounds on alpha_t. Reals.
#'
#' @param beta_bounds 2-element vector with the lower and upper bounds on beta_t. Non-negative reals.
#'
#' @param delta1_bounds 2-element vector with the lower and upper bounds on Delta_1. Reals.
#'
#' @param width_bounds 2-element vector with the lower and upper bounds on the width of the central
#' region in decision space (i.e., Delta_2 - Delta_1). Reals equal to or greater than 0.
#'
#' @param alpha_start Starting value(s) for alpha_t. Scalar or vector of reals.
#'
#' @param beta_start Starting value(s) for beta_t. Scalar or vector of non-negative reals.
#'
#' @param delta1_start Starting value(s) for delta_1. Scalar or vector of reals.
#'
#' @param width_start Starting value(s) for width. Scalar or vector of non-negative reals.
#'
#' @param eps_start Starting value(s) for epsilon parameters. Scalar or vector of reals in \eqn{[0, 1]}.
#'
#' @param kappa_start Starting value(s) for kappa parameters. Scalar or vector of reals in \eqn{[0, 1]}.
#'
#' @param model Choice of model to be fitted. An integer scalar in the range between 0 and 7
#' implying different assumptions about error parameters (see the table below)
#' that will be fitted to both presentation orders (i.e., test first and test
#' second), a 2-element vector with components in the same range (describing the
#' error model to be fitted to each presentation order), or the string "best"
#' (case insensitive) to find the best-fitting error model for each presentation
#' order by the criterion of choice (see the next input argument).
#' When \code{format <- "2AFC"}, models 0, 3, 5, and 7 are not permitted because
#' e_U must be in the fitted model;
#' when \code{format <- "equality"}, models 0, 4, 6, and 7 are not permitted
#' because e_F must be in the fitted model.
#' Setting \code{model <- "best"} takes these constraints into account when
#' searching for the best-fitting model.
#'
#' \tabular{rccc}{
#' \strong{model} \tab \strong{e_F} \tab \strong{e_U} \tab \strong{e_S} \cr
#' 0 \tab x \tab x \tab x \cr
#' 1 \tab o \tab o \tab o \cr
#' 2 \tab o \tab o \tab x \cr
#' 3 \tab o \tab x \tab o \cr
#' 4 \tab x \tab o \tab o \cr
#' 5 \tab o \tab x \tab x \cr
#' 6 \tab x \tab o \tab x \cr
#' 7 \tab x \tab x \tab o \cr
#' }
#'          x : parameter is not included in the model
#'          o : parameter is included in the model
#'
#' @param criterion Choice of criterion when \code{model <- "best"}.
#' Options are \code{"LogL"} (to use the value of the negative log-likelihood)
#' or \code{"BIC"} (to use the value of the Bayesian Information criterion).
#' Not used when argument model is not the string \code{"best"}.
#'
#' @param type Choice of \code{c("same")} (when the same psychophysical function holds
#' for test and standard stimuli) or \code{c("diff")} (when different psychophysical
#' functions are assumed to hold for test and standard stimuli).
#' Not used for detection data (i.e., when \code{standard <- -Inf}).
#'
#' @param plot Option to plot data and fitted functions. Logical scalar.
#'
#' @param disp Option to issue warnings or display progress information. Logical scalar.
#'
#' @return List of includes parameter estimates and goodness-of-fit measures and p-values.
#'
#' @importFrom stats optim qnorm pchisq
#'
#' @export
fit_TIM_2P <- function(data, standard, format,
                       alpha_bounds, beta_bounds, delta1_bounds, width_bounds,
                       alpha_start, beta_start, delta1_start, width_start, eps_start, kappa_start,
                       model, criterion, type, plot_graphs = FALSE, disp = FALSE){

    if (is.character(model)) {
        model <- tolower(model)
    }

    ### Sort boundary vectors
    alpha_bounds <- sort(alpha_bounds)
    beta_bounds <- sort(beta_bounds)
    delta1_bounds <- sort(delta1_bounds)
    width_bounds <- sort(width_bounds)

    ### Check conditions ###
    if (missing(data)){
        stop("No data supplied")
    } else if (length(dim(data))!=2 || nrow(data)!=7 || !is.numeric(data)) {
        stop("Invalid data (must be a 2D numeric matrix with 7 rows)")}
    if (missing(standard)) {
        stop("No standard argument supplied")
    } else if (!(is.numeric(standard) || length(standard)!=1 || (is.infinite(standard) && standard>0))) {
        stop("Invalid standard: must be a real scalar (including -Inf)")
    }
    if (missing(format)){
        stop("No format supplied")
    } else if (is.character(format)) {
        format  <- tolower(format)
        if (!any(identical(format,"ternary"),
                 identical(format,"2afc"),
                 identical(format,"equality"))) {
            stop("Wrong format (must be 'ternary', '2AFC', or 'equality', case insensitive)")
        }
    }
    if (!is.numeric(alpha_bounds) || length(alpha_bounds)!=2 || identical(alpha_bounds[1], alpha_bounds[2])) {
        stop("Invalid alpha_bounds: must be vector of two different real numbers")}
    if (any(!is.numeric(beta_bounds), beta_bounds<0) || length(beta_bounds)!=2 || identical(beta_bounds[1], beta_bounds[2])) {
        stop("Invalid beta_bounds:  must be vector of two different non-negative real numbers")}
    if (!is.numeric(delta1_bounds) || length(delta1_bounds)!=2 || identical(delta1_bounds[1], delta1_bounds[2])) {
        stop("Invalid delta1_bounds: must be vector of two different real numbers)")}

    if (any(!is.numeric(width_bounds), width_bounds<0) || length(width_bounds)!=2) {
        stop("Invalid width_bounds: must be vector of two non-negative reals")
    } else if (!identical(format,"2afc") && identical(width_bounds[1], width_bounds[2])){
        stop("Invalid contents of width_bounds: components must have different values")
    } else if (identical(format,"2afc") && identical(width_bounds[1], width_bounds[2]) && width_bounds[1]!=0){
        stop("Invalid contents of width_bounds: for 2AFC data, components must have different values or be zero)")
    }

    if (!is.numeric(alpha_start) || any(alpha_start>alpha_bounds[2]) || any(alpha_start<alpha_bounds[1])) {
        stop("Invalid alpha_start: must be real number within supplied alpha_bounds")}
    if (!is.numeric(beta_start) || any(beta_start<0) || any(beta_start>beta_bounds[2]) || any(beta_start<beta_bounds[1])) {
        stop("Invalid beta_start: must be non-negative real number within supplied beta_bounds)")}
    if (!is.numeric(delta1_start) || any(delta1_start>delta1_bounds[2]) || any(delta1_start<delta1_bounds[1])) {
        stop("Invalid delta1_start: must be real number within supplied delta1_bounds")}

    if (any(!is.numeric(width_start), width_start<0) || any(width_start>width_bounds[2]) || any(width_start<width_bounds[1])) {
        stop("Invalid width_start: must be non-negative reals within supplied width_bounds")}

    if (any(!is.numeric(eps_start), eps_start<0, eps_start>1)) {
        stop("Invalid eps_start (must be reals in [0,1])")}
    if (any(!is.wholeNumber(kappa_start), kappa_start<0, kappa_start>1)) {
        stop("Invalid kappa_start (must be reals in [0,1])")}

    if (any(data[2:7,]<0)) {
        stop("Negative counts in rows 2-7 of data ")}
    if (identical(format,"2afc") && any(data[c(3,6),]>0)) {
        stop("Rows 3 and 6 of array data must be filled with zeros when format='2AFC'")}
    if (identical(format,"2afc") && is.numeric(model) &&
        (any(model %in% c(0, 3, 5, 7)))){
        stop("model cannot be set to 0, 3, 5, or 7 when format='2AFC'")}
    if (identical(format,"equality") && any(data[c(2,5),]>0)){
        stop("Rows 2 and 5 of array data must be filled with zeros when format='equality'")}
    if (identical(format,"equality") && is.numeric(model) &&
        (any(model %in% c(0, 4, 6, 7)))) {
        stop("model cannot be set to 0, 4, 6, or 7 when format='equality'")}
    if (is.character(model)) {
        if (identical(model,"best")){
            if (missing(criterion) || (!identical(toupper(criterion),"BIC") && !identical(toupper(criterion),"LOGL"))) {
                stop("criterion must be either of c('LogL','BIC')")
            }
        } else {
            stop("Invalid value for model (the only valid string is 'best')")
        }
    } else if (!is.integer(model) || length(model)>2 || any(model<0) || any(model>7)) {
        stop('Invalid value for model: components must be integers between 0 and 7 (or "best")')
    } else {
        criterion = "Not applicable"
    }
    if (!identical(type,"same") && !identical(type,"diff")) {
        stop("Wrong string for type (must be 'same' or 'diff', case insensitive)")}
    if (!is.logical(plot) || length (plot)!=1) {
        plot = FALSE
        warning("Invalid value for plot (must be a logical scalar); 'plot' was set to FALSE")}
    if (!is.logical(disp) || length (disp)!=1) {
        disp = FALSE
        warning("Invalid value for disp (must be a logical scalar); 'disp' was set to FALSE")}

    ### Main ###
    output = fit_model(data, standard, format,
                       alpha_bounds, beta_bounds, delta1_bounds, width_bounds,
                       alpha_start, beta_start, delta1_start, width_start, eps_start, kappa_start,
                       model, criterion, type, plot_graphs, disp)

    if (plot_graphs){
        #plot_results()
    }
    return(output)
}

fit_model = function(data, standard, format,
                     alpha_bounds, beta_bounds, delta1_bounds, width_bounds,
                     alpha_start, beta_start, delta1_start, width_start, eps_start, kappa_start,
                     model, criterion, type, plot, disp){
    useBIC <- identical(toupper(criterion),"BIC")
    same_mu <- isTRUE(identical (type,"same") || is.infinite(standard))
    ter <- identical(format,"ternary")
    bin <- identical(format,"2afc")
    equ <- identical(format,"equality")
    det <- is.infinite(standard)
    noU <- bin && identical(width_bounds[1], width_bounds[2])
    logs <- c(same_mu, ter, bin, equ, det, noU)
    Lev <- data[1,]
    F1 <- data[2,]; U1 <- data[3,]; S1 <- data[4,]
    F2 <- data[5,]; U2 <- data[6,]; S2 <- data[7,]
    kita <- 1
    if (same_mu) kita <- kita+1
    if (!ter) kita <- kita+2
    upperwidth <- width_bounds[2]
    if (noU) {kita <- kita+3; upperwidth <- 0.0001; width_start <- upperwidth}
    NforBIC <- sum(colSums(data[2:4,]>0)) + sum(colSums(data[5:7,]>0))
    Bounds <- rbind(c(rep(alpha_bounds[1], times=2), rep(beta_bounds[1], times=2), delta1_bounds[1], width_bounds[1], rep(0, times=18)),
                    c(rep(alpha_bounds[2], times=2), rep(beta_bounds[2], times=2), delta1_bounds[2], upperwidth, rep(1, times=18)))

    # choice of algorithm
    if (all(model==0)) {
        eps_start <- 0
        kappa_start <- 0
    }
    if (identical(model,"best")) {
        if (ter) {models1 <- seq(0,7)}
        if (bin) {models1 <- c(1,2,4,6)}
        if (equ) {models1 <- c(1,2,3,5)}
        models2 <- models1
    } else if(length(model)==1) {
        models1 <- model; models2 <- model
    } else {
        models1 <- model[1]; models2 <- model[2]
    }
    num_models <- length(models1)*length(models2)
    num_iter <- length(alpha_start)*length(beta_start)*length(delta1_start)*
        length(width_start)*length(eps_start)*length(kappa_start)
    chosen_model <- NaN
    objective_function_min <- Inf
    CritMin <- .Machine$double.xmax
    nm <- 0
    for(modelLoop1 in models1) {
        for(modelLoop2 in models2) {
            models <- c(modelLoop1, modelLoop2)
            nm <- nm+1
            bounds <- Bounds
            stringm <- paste("Tried ",as.character(nm)," of ",as.character(num_models)," models (",
                             as.character(modelLoop1,width=1),",",as.character(modelLoop2,1), "); ", sep="")
            # proceed through initial values
            niter <- 0
            out.par <- matrix(NA,num_iter,50)
            out.value <- numeric(NA, num_iter)
            out.counts <-matrix(NA,num_iter,2)
            out.convg <- numeric(NA, num_iter)
            out.msg <- character(num_iter)
            OnBound <- array(NaN, c(num_iter,50,2))
            for (AlphaInit in alpha_start) {
                for (BetaInit in beta_start) {
                    for (Delta1Init in delta1_start) {
                        for (WidthInit in width_start) {
                            for (EpsInit in eps_start) {
                                for (KappaInit in kappa_start) {
                                    niter <- niter+1
                                    string <- paste("Tried ",as.character(niter),"/",as.character(num_iter)," sets of initial values; ")
                                    Initial <- c(AlphaInit, AlphaInit, BetaInit, BetaInit, Delta1Init, WidthInit)
                                    pstn <- c() ; kita2 <- 0
                                    for (i in c(1,2)) {
                                        if (models[i]==1) {
                                            Initial <- c(Initial, EpsInit, EpsInit, EpsInit, KappaInit, KappaInit, KappaInit)
                                            j <- length(Initial); pstn <- c(pstn, j-5, j-4, j-3)
                                            if (bin) {Initial[j-4] <- 0.9995; bounds[1,j-4] <- 0.999
                                            Initial[j-2] <- 0.0005; bounds[2,j-2] <- 0.001
                                            Initial[j] <- 0.9995; bounds[1,j] <- 0.999}
                                            if (equ) {Initial[j-5] <- 0.9995; bounds[1,j-5] <- 0.999
                                            Initial[j-1] <- 0.0005; bounds[2,j-1] <- 0.001
                                            Initial[j] <- 0.0005; bounds[2,j] <- 0.001}
                                            if (bin || equ) {kita2 <- kita2+2}
                                        } else if (models[i]==2) {
                                            Initial <- c(Initial, EpsInit, EpsInit,          KappaInit, KappaInit           )
                                            j <- length(Initial); pstn <- c(pstn, j-3, j-2)
                                            if (bin) {Initial[j-2] <- 0.9995; bounds[1,j-2] <- 0.999
                                            Initial[j-1] <- 0.0005; bounds[2,j-1] <- 0.001}
                                            if (equ) {Initial[j-3] <- 0.9995; bounds[1,j-3] <- 0.999
                                            Initial[j] <- 0.0005; bounds[2,j] <- 0.001}
                                            if (bin || equ) {kita2 <- kita2+1}
                                        } else if (models[i]==3) {
                                            Initial <- c(Initial, EpsInit,          EpsInit, KappaInit,            KappaInit)
                                            j <- length(Initial); pstn <- c(pstn,j-3,j-2)
                                            if (equ) {Initial[j-3] <- 0.9995; bounds[1,j-3] <- 0.999
                                            Initial[j] <- 0.0005; bounds[2,j] <- 0.001}
                                            if (equ) {kita2 <- kita2+1}
                                        } else if (models[i]==4) {
                                            Initial <- c(Initial,          EpsInit, EpsInit,            KappaInit, KappaInit)
                                            j <- length(Initial); pstn <- c(pstn, j-3, j-2)
                                            if (bin) {Initial[j-3] <- 0.9995; bounds[1,j-3] <- 0.999
                                            Initial[j] <- 0.9995; bounds[1,j] <- 0.999}
                                            if (bin) {kita2 <- kita2+1}
                                        } else if (models[i]==5) {
                                            Initial <- c(Initial, EpsInit,                   KappaInit                      )
                                            j <- length (Initial); pstn <- c(pstn, j-1)
                                            if (equ) {Initial[j-1] <- 0.9995; bounds[1,j-1] <- 0.999}
                                        } else if (models[i]==6) {
                                            Initial <- c(Initial,          EpsInit,                     KappaInit           )
                                            j <- length(Initial); pstn <- c(pstn, j-1)
                                            if (bin) {Initial[j-1] <- 0.9995; bounds[1,j-1] <- 0.999}
                                        } else if (models[i]==7) {
                                            Initial <- c(Initial,                   EpsInit,                       KappaInit)
                                            j <- length(Initial); pstn <- c(pstn, j-1)
                                        }
                                    }
                                    npars <- length(Initial)
                                    Initial <- log((Initial-bounds[1,1:npars])/(bounds[2,1:npars]-Initial))
                                    Initial[which(Initial==Inf)] <- 20
                                    Initial[which(Initial==-Inf)] <- -20
                                    # estimate parameters
                                    optimset <- list(trace=0, maxit=500, lmm=70, REPORT=1)
                                    out <- optim(Initial, LogLikelihood, gr=NULL, method="L-BFGS-B", control=optimset,
                                                 F1, U1, S1, F2, U2, S2, standard, Lev, models, logs, bounds[,1:npars],
                                                 lower=-Inf, upper=Inf)
                                    out.par[niter,1:npars] <- out$par
                                    out.value[niter] <- out$value
                                    out.counts[niter,] <-out$counts
                                    out.convg[niter] <- out$convergence
                                    out.msg[niter] <- out$message
                                    OnBound[niter,1:6,1] <- out$par[1:6] <= -20
                                    OnBound[niter,1:6,2] <- out$par[1:6] >= 20
                                    if (num_models==1 && disp) {
                                        if (out$value < objective_function_min) {
                                            objective_function_min <- out$value
                                            cat(paste(string,"-2LogL = ",as.character(2*out$value),"; best thus far ..."), sep=" ","\n")
                                        } else {
                                            cat(paste(string,"-2LogL = ",as.character(2*out$value)), sep=" ","\n")}
                                    }
                                }
                            }
                        }
                    }
                }
            }
            objective_function_min <- min(out.value)
            minIndex <- which.min(out.value)
            temp.pars <- out.par[minIndex,1:npars]
            temp.counts <- out.counts[minIndex,]
            temp.convg <- out.convg[minIndex]
            temp.msg <- out.msg[minIndex]
            temp.OnBound <- rbind(OnBound[minIndex,,1], OnBound[minIndex,,2])
            LogL <- 2*objective_function_min
            BIC <- 2*objective_function_min + (npars-kita-kita2)*log(NforBIC)
            if (!det) BIC2 <- 2*objective_function_min + (npars-kita-kita2-1)*log(NforBIC)
            ifelse(useBIC, Crit <- BIC, Crit <- LogL)
            doit <- isTRUE(useBIC || (nm==1 || !any(temp.pars[pstn]==0)))
            if (Crit<CritMin && doit) {
                CritMin <- Crit
                finalBIC <- BIC
                if (!det) finalBIC <- c(BIC, BIC2)
                finalLogL <- LogL
                chosen_model <- c(modelLoop1, modelLoop2)
                final.npars <- npars
                final.pars <- temp.pars
                final.counts <- temp.counts
                final.convg <- temp.convg
                final.msg <- temp.msg
                final.OnBound <- temp.OnBound
                if (num_models>1 && disp) {
                    if(useBIC){
                        cat(paste(stringm,"BIC =",as.character(Crit),"; best thus far ..."), sep=" ","\n")
                    } else {
                        cat(paste(stringm,"-2LogL =",as.character(Crit),"; best thus far ..."), sep=" ","\n")}
                }
            } else if (num_models>1 && disp) {
                if(useBIC){
                    cat(paste(stringm,"BIC =",as.character(Crit)), sep=" ","\n")
                } else {
                    cat(paste(stringm,"-2LogL =",as.character(Crit)), sep=" ","\n")}
            }
        }
    }
    npars <- final.npars; pars <- final.pars; BIC <- finalBIC; LogL <- finalLogL
    # extract estimated parameters, fix misleading kappas, and check boundary conditions
    all.pars <- GetParVec(pars, chosen_model, logs, Bounds)
    alpha_st <- all.pars[1]; alpha_t <- all.pars[2]
    beta_st <- all.pars[3]; beta_t <- all.pars[4]
    Delta_1 <- all.pars[5]; Delta_2 <- all.pars[6]
    lapses <- all.pars[-1:-6]
    warn <-c()
    if (final.OnBound[1,2]) {warn <- c(warn, " lower_Alpha ")}
    if (final.OnBound[2,2]) {warn <- c(warn, " upper_Alpha ")}
    if (final.OnBound[1,4]) {warn <- c(warn, " lower_Beta ")}
    if (final.OnBound[2,4]) {warn <- c(warn, " upper_Beta ")}
    if (final.OnBound[1,5]) {warn <- c(warn, " lower_Delta1 ")}
    if (final.OnBound[2,5]) {warn <- c(warn, " upper_Delta1 ")}
    if (final.OnBound[1,6] && width_bounds[1]>0) {warn <- c(warn, " lower_Width ")}
    if (final.OnBound[2,6]) {warn <- c(warn, " upper_Width ")}
    if (is.null(warn)){ warn <- " none "}
    # compute PSE and threshold, if applicable
    if (det) {
        PSE <- "not applicable"
        DL <- "not applicable"
        Thr_84 <- muinv(alpha_t,beta_t,qnorm(0.84,0,1)*sqrt(2))
    } else {
        PSE <- muinv(alpha_t,beta_t,mu(standard,alpha_st,beta_st))
        DL <- muinv(alpha_t,beta_t,mu(standard,alpha_st,beta_st)+qnorm(0.75,0,1)*sqrt(2)) - PSE
        Thr_84 <- "not applicable"
    }
    ### Goodness of fit test ###
    # observed frequencies
    O <- c(F1, U1, S1, F2, U2, S2)
    # number of observations per order of presentation and stimulus level
    n_1 <- F1 + U1 + S1
    n_2 <- F2 + U2 + S2
    # theoretical probabilities
    p.output <- Psy(standard, Lev, alpha_st, beta_st, alpha_t, beta_t, Delta_1, Delta_2, lapses)
    P_F1 <- p.output[1,]; P_U1 <- p.output[2,]; P_S1 <- p.output[3,]
    P_F2 <- p.output[4,]; P_U2 <- p.output[5,]; P_S2 <- p.output[6,]
    # expectations
    E <- c(P_F1*n_1, P_U1*n_1, P_S1*n_1, P_F2*n_2, P_U2*n_2, P_S2*n_2)
    # degrees of freedom
    ifelse (ter, factor <- 3, factor <- 2)
    dof <- (factor-1)*sum(n_1>0) + (factor-1)*sum(n_2>0) - (npars-kita-kita2)
    if (!det) dof <- c(dof, dof+1)
    # index of fit and significance
    EE <- E[E>0]; OO <- O[E>0]
    G2 <- 2*sum(OO*log(OO/EE),na.rm=T)
    X2 <- sum((OO-EE)^2/EE,na.rm=T)
    dof2 <- dof; dof2[dof<=0] <- -1; oldw <- getOption("warn"); options(warn = -1)
    pG2 <- pchisq(G2, df=dof2, lower.tail=F)
    pX2 <- pchisq(X2, df=dof2, lower.tail=F)
    options(warn = oldw)
    params <- c(alpha_st, alpha_t, beta_st, beta_t, Delta_1, Delta_2, lapses)
    if (plot) {plot_results(params, standard, Lev, F1, U1, S1, F2, U2, S2, PSE, Thr_84, logs)}
    # create report
    # rearrange parameters to return
    if (chosen_model[1]==0) {
        params[7:15] <- NA
    } else if (chosen_model[1]==2) {
        params[c(9, 14:15)] <- NA
    } else if (chosen_model[1]==3) {
        params[c(8, 12:13)] <- NA
    } else if (chosen_model[1]==4) {
        params[c(7, 10:11)] <- NA
    } else if (chosen_model[1]==5) {
        params[c(8:9, 12:15)] <- NA
    } else if (chosen_model[1]==6) {
        params[c(7, 9:11, 14:15)] <- NA
    } else if (chosen_model[1]==7) {
        params[c(7:8, 10:13)] <- NA
    }
    if (chosen_model[2]==0) {
        params[16:24] <- NA
    } else if (chosen_model[2]==2) {
        params[c(18, 23:24)] <- NA
    } else if (chosen_model[2]==3) {
        params[c(17, 21:22)] <- NA
    } else if (chosen_model[2]==4) {
        params[c(16, 19:20)] <- NA
    } else if (chosen_model[2]==5) {
        params[c(17:18, 21:24)] <- NA
    } else if (chosen_model[2]==6) {
        params[c(16, 18:20, 23:24)] <- NA
    } else if (chosen_model[2]==7) {
        params[c(16:17, 19:22)] <- NA
    }
    if (noU) params[c(6, 8, 12, 13, 17, 21, 22)] <- NA
    ifelse (det, type2 <- "same (detection)", type2 <- type)
    if (ter) {
        format2 <- paste(format, " (respond first, second, or uncertain)")
    } else if (bin) {
        format2 <- paste(format, " (respond first or second, guessing when uncertain)")
    } else if (equ) {
        format2 <- paste(format, " (respond same or different)")
    }
    numfree <- npars-kita-kita2
    if (!det) numfree <- c(numfree, numfree-1)
    report <- list(Problem = "Fit of the indecision model to dual-presentation data; L-BFGS-B",
                   optim_output = c(final.convg, final.counts[1]),
                   data = data,
                   format = format2,
                   standard = standard,
                   Usermodel = model,
                   fitted_model = chosen_model,
                   criterion = criterion,
                   type = type2,
                   num_free_parameters = numfree,
                   num_cells = factor*sum(n_1>0) + factor*sum(n_2>0),
                   num_exp_below5 = sum(E<5),
                   num_exp_below5_obs_above0 = sum((E<5)*(O>0)),
                   num_exp_below1 = sum(E<1),
                   num_exp_below1_obs_above0 = sum((E<1)*(O>0)),
                   df = dof,
                   ChiSquareStatistic = X2,
                   ChiSquareP_value = pX2,
                   likelihood_ratio_statistic = G2,
                   likelihood_ratio_p_value = pG2,
                   BIC = BIC,
                   Neg2LogL = LogL,
                   alpha_bounds = alpha_bounds,
                   beta_bounds = beta_bounds,
                   delta1_bounds = delta1_bounds,
                   width_bounds = width_bounds,
                   boundaries_reached = warn,
                   alpha_st = params[1],
                   alpha_t = params[2],
                   beta_st = params[3],
                   beta_t = params[4],
                   mu_s_at_standard_level = mu(standard,params[1],params[3]),
                   mu_t_at_standard_level = mu(standard,params[2],params[4]),
                   delta_1 = params[5],
                   delta_2 = params[6],
                   epsilon_F_1 = params[7],
                   epsilon_U_1 = params[8],
                   epsilon_S_1 = params[9],
                   kappa_FintoU_1 = params[10],
                   kappa_FintoS_1 = params[11],
                   kappa_UintoF_1 = params[12],
                   kappa_UintoS_1 = params[13],
                   kappa_SintoF_1 = params[14],
                   kappa_SintoU_1 = params[15],
                   epsilon_F_2 = params[16],
                   epsilon_U_2 = params[17],
                   epsilon_S_2 = params[18],
                   kappa_FintoU_2 = params[19],
                   kappa_FintoS_2 = params[20],
                   kappa_UintoF_2 = params[21],
                   kappa_UintoS_2 = params[22],
                   kappa_SintoF_2 = params[23],
                   kappa_SintoU_2 = params[24],
                   same = logs[1],
                   ter = logs[2],
                   bin = logs[3],
                   equ = logs[4],
                   det = logs[5],
                   lapses = params[7:24],
                   PSE = PSE,
                   DL = DL,
                   Threshold_84 = Thr_84
    )

    class(report) <- "fit"
    return(report)
}
