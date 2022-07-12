#' Fit the ternary indecision model to detection or discrimination data from
#' a dual-presentation task
#'
#' @description Further details:  Garc?a-P?rez, M.A. & Alcal?-Quintana, R. (2017).
#' The indecision model of psychophysical performance in
#' dual-presentation tasks: Parameter estimation and
#' comparative analysis of response formats.
#' Frontiers in Psychology, 8:1142. <https://doi.org/10.3389/fpsyg.2017.01142>
#'
#' @param Data 7-by-N array with data from the dual-presentaion task. The
#' first row gives the N levels of the test stimulus; rows 2-4 respectively
#' give the count of "first", "undecided", and "second" responses at each
#' level when the test stimulus was presented in the first interval/position;
#' rows 5-7 analogously give the count of "first", "undecided", and "second"
#' responses at each level when the test stimulus was presented in the second
#' interval/position. If \code{Format <- "2AFC"} or
#' \code{Format <- "equality"}, array Data still has the same number of rows
#' but its content is subject to the constraints listed under argument Format
#' below.
#'
#' @param Standard Level of the standard stimulus. For detection data (i.e.,
#' with a nominally null standard), set \code{Standard <- -Inf}; for
#' discrimination data (i.e., with a non-null standard), set Standard to the
#' appropriate level, in the same units in which levels of the test stimulus
#' are given in the first row of array Data.
#'
#' @param Format Out of \code{c("ternary", "2AFC", or "equality")} to indicate
#' that data were collected with a ternary response format (first, second, or
#' undecided) a 2AFC response format (forcing first or second responses by
#' guessing when uncertain), or the same-different response format (both stimuli
#' are subjective equal or they are subjectively different).
#' When \code{Format <- "2AFC"}, rows 3 and 6 of array Data must be filled with
#' zeros (as no "undecided" responses are given under this format);
#' When \code{Format <- "equality"}, rows 2 and 5 of array Data must be filled with zeros
#' (so "different" responses are stored as "second" responses).
#'
#' @param AlphaBounds 2-element vector with the lower and upper bounds on Alpha_t. Reals.
#'
#' @param BetaBounds 2-element vector with the lower and upper bounds on Beta_t. Non-negative reals.
#'
#' @param Delta1Bounds 2-element vector with the lower and upper bounds on Delta_1. Reals.
#'
#' @param WidthBounds 2-element vector with the lower and upper bounds on the width of the central
#' region in decision space (i.e., Delta_2 - Delta_1). Reals equal to or greater than 0.
#'
#' @param AlphaStart Starting value(s) for Alpha_t. Scalar or vector of reals.
#'
#' @param BetaStart Starting value(s) for Beta_t. Scalar or vector of non-negative reals.
#'
#' @param Delta1Start Starting value(s) for Delta_1. Scalar or vector of reals.
#'
#' @param WidthStart Starting value(s) for width. Scalar or vector of non-negative reals.
#'
#' @param EpsStart Starting value(s) for epsilon parameters. Scalar or vector of reals in \eqn{[0, 1]}.
#'
#' @param KappaStart Starting value(s) for kappa parameters. Scalar or vector of reals in \eqn{[0, 1]}.
#'
#' @param Model Choice of model to be fitted. An integer scalar in the range between 0 and 7
#' implying different assumptions about error parameters (see the table below)
#' that will be fitted to both presentation orders (i.e., test first and test
#' second), a 2-element vector with components in the same range (describing the
#' error model to be fitted to each presentation order), or the string "best"
#' (case insensitive) to find the best-fitting error model for each presentation
#' order by the criterion of choice (see the next input argument).
#' When \code{Format <- "2AFC"}, models 0, 3, 5, and 7 are not permitted because
#' e_U must be in the fitted model;
#' when \code{Format <- "equality"}, models 0, 4, 6, and 7 are not permitted
#' because e_F must be in the fitted model.
#' Setting \code{Model <- "best"} takes these constraints into account when
#' searching for the best-fitting model.
#'
#' \tabular{rccc}{
#' \strong{Model} \tab \strong{e_F} \tab \strong{e_U} \tab \strong{e_S} \cr
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
#' @param Criterion Choice of criterion when \code{Model <- "best"}.
#' Options are \code{"LogL"} (to use the value of the negative log-likelihood)
#' or \code{"BIC"} (to use the value of the Bayesian Information Criterion).
#' Not used when argument Model is not the string \code{"best"}.
#'
#' @param Type Choice of \code{c("same")} (when the same psychophysical function holds
#' for test and standard stimuli) or \code{c("diff")} (when different psychophysical
#' functions are assumed to hold for test and standard stimuli).
#' Not used for detection data (i.e., when \code{Standard <- -Inf}).
#'
#' @param Plot Option to plot data and fitted functions. Logical scalar.
#'
#' @param Disp Option to issue warnings or display progress information. Logical scalar.
#'
#' @return List of includes parameter estimates and goodness-of-fit measures and p-values.
#'
#' @importFrom stats optim qnorm pchisq
#'
#' @export
fit_TIM_2P <- function(Data, Standard, Format,
                       AlphaBounds, BetaBounds, Delta1Bounds, WidthBounds,
                       AlphaStart, BetaStart, Delta1Start, WidthStart, EpsStart, KappaStart,
                       Model = "best", Criterion = "LogL", Type, Plot = FALSE, Disp = FALSE){

    ### Convert strings to lower case ###
    if (is.character(Format))  Format  <- tolower(Format)
    if (is.character(Model)) Model <- tolower(Model)

    ### Sort boundary vectors
    AlphaBounds <- sort(AlphaBounds)
    BetaBounds <- sort(BetaBounds)
    Delta1Bounds <- sort(Delta1Bounds)
    WidthBounds <- sort(WidthBounds)

    ### Check conditions ###
    if (length(dim(Data))!=2 || nrow(Data)!=7 || !is.double(Data)) {
        stop("Invalid Data (must be a 2D numeric matrix with 7 rows)")}
    if (!is.double(Standard) || length(Standard)!=1 || (is.infinite(Standard) && Standard>0)) {
        stop("Invalid Standard (must be a real scalar, including -Inf)")}
    if (!is.character(Format) || !any(identical(Format,"ternary"),
                                      identical(Format,"2afc"),
                                      identical(Format,"equality"))) {
        stop("Wrong Format (must be 'ternary', '2AFC', or 'equality', case insensitive)")}
    if (!is.double(AlphaBounds)) {
        stop("Invalid AlphaBounds (must be reals)")}
    if (length(AlphaBounds)!=2 ) {
        stop("Invalid size of AlphaBounds (must have two components only)")}
    if (identical(AlphaBounds[1], AlphaBounds[2])) {
        stop("Invalid contents of AlphaBounds (components must have different values)")}
    if (any(!is.double(BetaBounds), BetaBounds<0)) {
        stop("Invalid BetaBounds (must be non-negative reals)")}
    if (length(BetaBounds)!=2 ) {
        stop("Invalid size of BetaBounds (must have two components only)")}
    if (identical(BetaBounds[1], BetaBounds[2])) {
        stop("Invalid contents of BetaBounds (components must have different values)")}
    if (!is.double(Delta1Bounds)) {
        stop("Invalid Delta1Bounds (must be reals)")}
    if (length(Delta1Bounds)!=2 ) {
        stop("Invalid size of Delta1Bounds (must have two components only)")}
    if (identical(Delta1Bounds[1], Delta1Bounds[2])) {
        stop("Invalid contents of Delta1Bounds (components must have different values)")}
    if (any(!is.double(WidthBounds), WidthBounds<0)) {
        stop("Invalid WidthBounds (must be non-negative reals)")}
    if (length(WidthBounds)!=2) {
        stop("Invalid size of WidthBounds (must have two components only)")}
    if (!identical(Format,"2afc") && identical(WidthBounds[1], WidthBounds[2])) {
        stop("Invalid contents of WidthBounds (components must have different values)")}
    if (identical(Format,"2afc") && identical(WidthBounds[1], WidthBounds[2]) && WidthBounds[1]!=0) {
        stop("Invalid contents of WidthBounds (for 2AFC data, components must have different values or zeros)")}
    if (!is.double(AlphaStart)) {
        stop("Invalid AlphaStart (must be reals)")}
    if (any(AlphaStart>AlphaBounds[2]) || any(AlphaStart<AlphaBounds[1])) {
        stop("Invalid AlphaStart (one or more values are not within AlphaBounds)")}
    if (any(!is.double(BetaStart), BetaStart<0)) {
        stop("Invalid BetaStart (must be non-negative reals)")}
    if (any(BetaStart>BetaBounds[2]) || any(BetaStart<BetaBounds[1])) {
        stop("Invalid BetaStart (one or more values are not within BetaBounds)")}
    if (!is.double(Delta1Start)) {
        stop("Invalid Delta1Start (must be reals)")}
    if (any(Delta1Start>Delta1Bounds[2]) || any(Delta1Start<Delta1Bounds[1])) {
        stop("Invalid Delta1Start (one or more values are not within Delta1Bounds)")}
    if (any(!is.double(WidthStart), WidthStart<0)) {
        stop("Invalid WidthStart (must be non-negative reals)")}
    if (any(WidthStart>WidthBounds[2]) || any(WidthStart<WidthBounds[1])) {
        stop("Invalid WidthStart (one or more values are not within WidthBounds)")}
    if (any(!is.double(EpsStart), EpsStart<0, EpsStart>1)) {
        stop("Invalid EpsStart (must be reals in [0,1])")}
    if (any(!is.double(KappaStart), KappaStart<0, KappaStart>1)) {
        stop("Invalid KappaStart (must be reals in [0,1])")}
    if (any(Data[2:7,]<0)) {
        stop("Negative counts in rows 2-7 of Data ")}
    if (identical(Format,"2afc") && any(Data[c(3,6),]>0)) {
        stop("Rows 3 and 6 of array Data must be filled with zeros when Format='2AFC'")}
    if (identical(Format,"2afc") && is.double(Model) &&
        (any(Model==0)||any(Model==3) || any(Model==5) || any(Model==7))){
        stop("Model cannot be set to 0, 3, 5, or 7 when Format='2AFC'")}
    if (identical(Format,"equality") && any(Data[c(2,5),]>0)){
        stop("Rows 2 and 5 of array Data must be filled with zeros when Format='equality'")}
    if (identical(Format,"equality") && is.double(Model) &&
        (any(Model==0)||any(Model==4) || any(Model==6) || any(Model==7))) {
        stop("Model cannot be set to 0, 4, 6, or 7 when Format='equality'")}
    if (is.character(Model) && !identical(Model,"best")) {
        stop("Invalid value for Model (the only valid string is 'best')")}
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
    if (is.double(Model) &&
        (length(Model)>2 || !is.wholenumber(Model) || any(Model<0) || any(Model>7))) {
        stop("Invalid value for Model (components must be integers between 0 and 7)")}
    if (!identical(toupper(Criterion),"BIC") && !identical(toupper(Criterion),"LOGL")) {
        stop("Wrong string for Type (must be 'LogL' or 'BIC', case insensitive)")}
    if (!identical(Type,"same") && !identical(Type,"diff")) {
        stop("Wrong string for Type (must be 'same' or 'diff', case insensitive)")}
    if (!is.logical(Plot) || length (Plot)!=1) {
        stop("Invalid value for Plot (must be a logical scalar)")}
    if (!is.logical(Disp) || length (Disp)!=1) {
        stop("Invalid value for Disp (must be a logical scalar)")}

    ### Main ###
    useBIC <- identical(toupper(Criterion),"BIC")
    same_mu <- isTRUE(identical (Type,"same") || is.infinite(Standard))
    ter <- identical(Format,"ternary")
    bin <- identical(Format,"2afc")
    equ <- identical(Format,"equality")
    det <- is.infinite(Standard)
    noU <- bin && identical(WidthBounds[1], WidthBounds[2])
    logs <- c(same_mu, ter, bin, equ, det, noU)
    Lev <- Data[1,]
    F1 <- Data[2,]; U1 <- Data[3,]; S1 <- Data[4,]
    F2 <- Data[5,]; U2 <- Data[6,]; S2 <- Data[7,]
    kita <- 1
    if (same_mu) kita <- kita+1
    if (!ter) kita <- kita+2
    upperwidth <- WidthBounds[2]
    if (noU) {kita <- kita+3; upperwidth <- 0.0001; WidthStart <- upperwidth}
    NforBIC <- sum(colSums(Data[2:4,]>0)) + sum(colSums(Data[5:7,]>0))
    Bounds <- rbind(c(rep(AlphaBounds[1], times=2), rep(BetaBounds[1], times=2), Delta1Bounds[1], WidthBounds[1], rep(0, times=18)),
                    c(rep(AlphaBounds[2], times=2), rep(BetaBounds[2], times=2), Delta1Bounds[2], upperwidth, rep(1, times=18)))

    # choice of algorithm
    if (all(Model==0)) {
        EpsStart <- 0
        KappaStart <- 0
    }
    if (identical(Model,"best")) {
        if (ter) {Models1 <- seq(0,7)}
        if (bin) {Models1 <- c(1,2,4,6)}
        if (equ) {Models1 <- c(1,2,3,5)}
        Models2 <- Models1
    } else {
        if (length(Model)==1) {
            Models1 <- Model; Models2 <- Model
        } else {
            Models1 <- Model[1]; Models2 <- Model[2]}
    }
    NumModels <- length(Models1)*length(Models2)
    NumIter <- length(AlphaStart)*length(BetaStart)*length(Delta1Start)*
        length(WidthStart)*length(EpsStart)*length(KappaStart)
    ChosenModel <- NaN
    ObjectiveFunctionMin <- Inf
    CritMin <- .Machine$double.xmax
    nm <- 0
    for(ModelLoop1 in Models1) {
        for(ModelLoop2 in Models2) {
            Models <- c(ModelLoop1, ModelLoop2)
            nm <- nm+1
            bounds <- Bounds
            stringm <- paste("Tried ",as.character(nm)," of ",as.character(NumModels)," models (",
                             as.character(ModelLoop1,width=1),",",as.character(ModelLoop2,1), "); ", sep="")
            # proceed through initial values
            niter <- 0
            out.par <- matrix(NA,NumIter,50)
            out.value <- vector("double",NumIter)+NA
            out.counts <-matrix(NA,NumIter,2)
            out.convg <- vector("double",NumIter)+NA
            out.msg <- vector("character",NumIter)
            OnBound <- array(NaN, c(NumIter,50,2))
            for (AlphaInit in AlphaStart) {
                for (BetaInit in BetaStart) {
                    for (Delta1Init in Delta1Start) {
                        for (WidthInit in WidthStart) {
                            for (EpsInit in EpsStart) {
                                for (KappaInit in KappaStart) {
                                    niter <- niter+1
                                    string <- paste("Tried ",as.character(niter),"/",as.character(NumIter)," sets of initial values; ")
                                    Initial <- c(AlphaInit, AlphaInit, BetaInit, BetaInit, Delta1Init, WidthInit)
                                    pstn <- c() ; kita2 <- 0
                                    for (i in c(1,2)) {
                                        if (Models[i]==1) {
                                            Initial <- c(Initial, EpsInit, EpsInit, EpsInit, KappaInit, KappaInit, KappaInit)
                                            j <- length(Initial); pstn <- c(pstn, j-5, j-4, j-3)
                                            if (bin) {Initial[j-4] <- 0.9995; bounds[1,j-4] <- 0.999
                                            Initial[j-2] <- 0.0005; bounds[2,j-2] <- 0.001
                                            Initial[j] <- 0.9995; bounds[1,j] <- 0.999}
                                            if (equ) {Initial[j-5] <- 0.9995; bounds[1,j-5] <- 0.999
                                            Initial[j-1] <- 0.0005; bounds[2,j-1] <- 0.001
                                            Initial[j] <- 0.0005; bounds[2,j] <- 0.001}
                                            if (bin || equ) {kita2 <- kita2+2}
                                        } else if (Models[i]==2) {
                                            Initial <- c(Initial, EpsInit, EpsInit,          KappaInit, KappaInit           )
                                            j <- length(Initial); pstn <- c(pstn, j-3, j-2)
                                            if (bin) {Initial[j-2] <- 0.9995; bounds[1,j-2] <- 0.999
                                            Initial[j-1] <- 0.0005; bounds[2,j-1] <- 0.001}
                                            if (equ) {Initial[j-3] <- 0.9995; bounds[1,j-3] <- 0.999
                                            Initial[j] <- 0.0005; bounds[2,j] <- 0.001}
                                            if (bin || equ) {kita2 <- kita2+1}
                                        } else if (Models[i]==3) {
                                            Initial <- c(Initial, EpsInit,          EpsInit, KappaInit,            KappaInit)
                                            j <- length(Initial); pstn <- c(pstn,j-3,j-2)
                                            if (equ) {Initial[j-3] <- 0.9995; bounds[1,j-3] <- 0.999
                                            Initial[j] <- 0.0005; bounds[2,j] <- 0.001}
                                            if (equ) {kita2 <- kita2+1}
                                        } else if (Models[i]==4) {
                                            Initial <- c(Initial,          EpsInit, EpsInit,            KappaInit, KappaInit)
                                            j <- length(Initial); pstn <- c(pstn, j-3, j-2)
                                            if (bin) {Initial[j-3] <- 0.9995; bounds[1,j-3] <- 0.999
                                            Initial[j] <- 0.9995; bounds[1,j] <- 0.999}
                                            if (bin) {kita2 <- kita2+1}
                                        } else if (Models[i]==5) {
                                            Initial <- c(Initial, EpsInit,                   KappaInit                      )
                                            j <- length (Initial); pstn <- c(pstn, j-1)
                                            if (equ) {Initial[j-1] <- 0.9995; bounds[1,j-1] <- 0.999}
                                        } else if (Models[i]==6) {
                                            Initial <- c(Initial,          EpsInit,                     KappaInit           )
                                            j <- length(Initial); pstn <- c(pstn, j-1)
                                            if (bin) {Initial[j-1] <- 0.9995; bounds[1,j-1] <- 0.999}
                                        } else if (Models[i]==7) {
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
                                                 F1, U1, S1, F2, U2, S2, Standard, Lev, Models, logs, bounds[,1:npars],
                                                 lower=-Inf, upper=Inf)
                                    out.par[niter,1:npars] <- out$par
                                    out.value[niter] <- out$value
                                    out.counts[niter,] <-out$counts
                                    out.convg[niter] <- out$convergence
                                    out.msg[niter] <- out$message
                                    OnBound[niter,1:6,1] <- out$par[1:6] <= -20
                                    OnBound[niter,1:6,2] <- out$par[1:6] >= 20
                                    if (NumModels==1 && Disp) {
                                        if (out$value < ObjectiveFunctionMin) {
                                            ObjectiveFunctionMin <- out$value
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
            ObjectiveFunctionMin <- min(out.value)
            minIndex <- which.min(out.value)
            temp.pars <- out.par[minIndex,1:npars]
            temp.counts <- out.counts[minIndex,]
            temp.convg <- out.convg[minIndex]
            temp.msg <- out.msg[minIndex]
            temp.OnBound <- rbind(OnBound[minIndex,,1], OnBound[minIndex,,2])
            LogL <- 2*ObjectiveFunctionMin
            BIC <- 2*ObjectiveFunctionMin + (npars-kita-kita2)*log(NforBIC)
            if (!det) BIC2 <- 2*ObjectiveFunctionMin + (npars-kita-kita2-1)*log(NforBIC)
            ifelse(useBIC, Crit <- BIC, Crit <- LogL)
            doit <- isTRUE(useBIC || (nm==1 || !any(temp.pars[pstn]==0)))
            if (Crit<CritMin && doit) {
                CritMin <- Crit
                FinalBIC <- BIC
                if (!det) FinalBIC <- c(BIC, BIC2)
                FinalLogL <- LogL
                ChosenModel <- c(ModelLoop1, ModelLoop2)
                Final.npars <- npars
                Final.pars <- temp.pars
                Final.counts <- temp.counts
                Final.convg <- temp.convg
                Final.msg <- temp.msg
                Final.OnBound <- temp.OnBound
                if (NumModels>1 && Disp) {
                    if(useBIC){
                        cat(paste(stringm,"BIC =",as.character(Crit),"; best thus far ..."), sep=" ","\n")
                    } else {
                        cat(paste(stringm,"-2LogL =",as.character(Crit),"; best thus far ..."), sep=" ","\n")}
                }
            } else if (NumModels>1 && Disp) {
                if(useBIC){
                    cat(paste(stringm,"BIC =",as.character(Crit)), sep=" ","\n")
                } else {
                    cat(paste(stringm,"-2LogL =",as.character(Crit)), sep=" ","\n")}
            }
        }
    }
    npars <- Final.npars; pars <- Final.pars; BIC <- FinalBIC; LogL <- FinalLogL
    # extract estimated parameters, fix misleading kappas, and check boundary conditions
    all.pars <- GetParVec(pars, ChosenModel, logs, Bounds)
    Alpha_st <- all.pars[1]; Alpha_t <- all.pars[2]
    Beta_st <- all.pars[3]; Beta_t <- all.pars[4]
    Delta_1 <- all.pars[5]; Delta_2 <- all.pars[6]
    lapses <- all.pars[-1:-6]
    warn <-c()
    if (Final.OnBound[1,2]) {warn <- c(warn, " lower_Alpha ")}
    if (Final.OnBound[2,2]) {warn <- c(warn, " upper_Alpha ")}
    if (Final.OnBound[1,4]) {warn <- c(warn, " lower_Beta ")}
    if (Final.OnBound[2,4]) {warn <- c(warn, " upper_Beta ")}
    if (Final.OnBound[1,5]) {warn <- c(warn, " lower_Delta1 ")}
    if (Final.OnBound[2,5]) {warn <- c(warn, " upper_Delta1 ")}
    if (Final.OnBound[1,6] && WidthBounds[1]>0) {warn <- c(warn, " lower_Width ")}
    if (Final.OnBound[2,6]) {warn <- c(warn, " upper_Width ")}
    if (is.null(warn)){ warn <- " none "}
    # compute PSE and threshold, if applicable
    if (det) {
        PSE <- "not applicable"
        DL <- "not applicable"
        Thr_84 <- muinv(Alpha_t,Beta_t,qnorm(0.84,0,1)*sqrt(2))
    } else {
        PSE <- muinv(Alpha_t,Beta_t,mu(Standard,Alpha_st,Beta_st))
        DL <- muinv(Alpha_t,Beta_t,mu(Standard,Alpha_st,Beta_st)+qnorm(0.75,0,1)*sqrt(2)) - PSE
        Thr_84 <- "not applicable"
    }
    ### Goodness of fit test ###
    # observed frequencies
    O <- c(F1, U1, S1, F2, U2, S2)
    # number of observations per order of presentation and stimulus level
    n_1 <- F1 + U1 + S1
    n_2 <- F2 + U2 + S2
    # theoretical probabilities
    p.output <- Psy(Standard, Lev, Alpha_st, Beta_st, Alpha_t, Beta_t, Delta_1, Delta_2, lapses)
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
    params <- c(Alpha_st, Alpha_t, Beta_st, Beta_t, Delta_1, Delta_2, lapses)
    if (Plot) {plot_results(params, Standard, Lev, F1, U1, S1, F2, U2, S2, PSE, Thr_84, logs)}
    # create report
    # rearrange parameters to return
    if (ChosenModel[1]==0) {params[7:15] <- NA
    } else if (ChosenModel[1]==2) {params[9] <- NA; params[14:15] <- NA
    } else if (ChosenModel[1]==3) {params[8] <- NA; params[12:13] <- NA
    } else if (ChosenModel[1]==4) {params[7] <- NA; params[10:11] <- NA
    } else if (ChosenModel[1]==5) {params[8:9] <- NA; params[12:15] <- NA
    } else if (ChosenModel[1]==6) {params[7] <- NA; params[9] <- NA
    params[10:11] <- NA; params[14:15] <- NA
    } else if (ChosenModel[1]==7) {params[7:8] <- NA; params[10:13] <- NA
    }
    if (ChosenModel[2]==0) {params[16:24] <- NA
    } else if (ChosenModel[2]==2) {params[18] <- NA; params[23:24] <- NA
    } else if (ChosenModel[2]==3) {params[17] <- NA; params[21:22] <- NA
    } else if (ChosenModel[2]==4) {params[16] <- NA; params[19:20] <- NA
    } else if (ChosenModel[2]==5) {params[17:18] <- NA; params[21:24] <- NA
    } else if (ChosenModel[2]==6) {params[16] <- NA; params[18] <- NA
    params[19:20] <- NA; params[23:24] <- NA
    } else if (ChosenModel[2]==7) {params[16:17] <- NA; params[19:22] <- NA
    }
    if (noU) params[c(6, 8, 12, 13, 17, 21, 22)] <- NA
    ifelse (det, Type2 <- "same (detection)", Type2 <- Type)
    ifelse (NumModels==1, Criterion2 <- "not applicable", Criterion2 <- Criterion)
    if (ter) {
        Format2 <- paste(Format, " (respond first, second, or uncertain)")
    } else if (bin) {
        Format2 <- paste(Format, " (respond first or second, guessing when uncertain)")
    } else if (equ) {
        Format2 <- paste(Format, " (respond same or different)")
    }
    numfree <- npars-kita-kita2
    if (!det) numfree <- c(numfree, numfree-1)
    report <- list(Problem = "Fit of the indecision model to dual-presentation data; L-BFGS-B",
                   optim_output = c(Final.convg, Final.counts[1]),
                   Data = Data,
                   Format = Format2,
                   Standard = Standard,
                   UserModel = Model,
                   FittedModel = ChosenModel,
                   Criterion = Criterion2,
                   Type = Type2,
                   NumFreeParameters = numfree,
                   NumCells = factor*sum(n_1>0) + factor*sum(n_2>0),
                   NumExpBelow5 = sum(E<5),
                   NumExpBelow5_ObsAbove0 = sum((E<5)*(O>0)),
                   NumExpBelow1 = sum(E<1),
                   NumExpBelow1_ObsAbove0 = sum((E<1)*(O>0)),
                   DegreesOfFreedom = dof,
                   ChiSquareStatistic = X2,
                   ChiSquareP_value = pX2,
                   LikelihoodRatioStatistic = G2,
                   LikelihoodRatioP_value = pG2,
                   BIC = BIC,
                   Neg2LogL = LogL,
                   ALphaBounds = AlphaBounds,
                   BetaBounds = BetaBounds,
                   Delta1Bounds = Delta1Bounds,
                   WidthBounds = WidthBounds,
                   BoundariesReached = warn,
                   Alpha_t = params[2],
                   Beta_t = params[4],
                   Mu_s_at_StandardLevel = mu(Standard,params[1],params[3]),
                   Mu_t_at_StandardLevel = mu(Standard,params[2],params[4]),
                   Delta_1 = params[5],
                   Delta_2 = params[6],
                   Epsilon_F_1 = params[7],
                   Epsilon_U_1 = params[8],
                   Epsilon_S_1 = params[9],
                   Kappa_FintoU_1 = params[10],
                   Kappa_FintoS_1 = params[11],
                   Kappa_UintoF_1 = params[12],
                   Kappa_UintoS_1 = params[13],
                   Kappa_SintoF_1 = params[14],
                   Kappa_SintoU_1 = params[15],
                   Epsilon_F_2 = params[16],
                   Epsilon_U_2 = params[17],
                   Epsilon_S_2 = params[18],
                   Kappa_FintoU_2 = params[19],
                   Kappa_FintoS_2 = params[20],
                   Kappa_UintoF_2 = params[21],
                   Kappa_UintoS_2 = params[22],
                   Kappa_SintoF_2 = params[23],
                   Kappa_SintoU_2 = params[24],
                   PSE = PSE,
                   DL = DL,
                   Threshold_84 = Thr_84)
    class(report) <- "fit"
    return(report)
}
