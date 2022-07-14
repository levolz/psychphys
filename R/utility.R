#' whole number check
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

#' Get parameters out of vector of estimates
GetParVec <- function (p, model, logs, bounds) {
    npars <- length(p)
    p1 <- p>700
    p <- colSums(rbind((1-p1)*(bounds[1,1:npars]+bounds[2,1:npars]*exp(p))/(1+exp(p)),
                       p1*bounds[2,1:npars]), na.rm=T)
    if (npars>6) {
        p[which(p[7:npars] < 0.001)+6] <- 0
        p[which(p[7:npars] > 0.999)+6] <- 1
    }
    alpha_st <- p[1]; alpha_t <- p[2]; Beta_st <- p[3]; beta_t <- p[4]; delta_1 <- p[5]; delta_2 <- p[5]+p[6]
    if (logs[1]) {alpha_st <- alpha_t; Beta_st <- beta_t}
    if (logs[6]) delta_2 <- delta_1
    lapses <- c(); add <- 0
    for (i in c(1,2)) {
        if (model[i]==0) {
            lapses <- c(lapses, rep(0,times=9)); add <- 0
        } else if (model[i]==1) {
            lapses <- c(lapses, p[add+(7:9)], p[add+10], 1-p[add+10], p[add+11], 1-p[add+11], p[add+12], 1-p[add+12]); add <- 6
        } else if (model[i]==2) {
            lapses <- c(lapses, p[add+(7:8)], 0, p[add+9], 1-p[add+9], p[add+10], 1-p[add+10], 0, 0); add <- 4
        } else if (model[i]==3) {
            lapses <- c(lapses, p[add+7], 0, p[add+8], p[add+9], 1-p[add+9], 0, 0, p[add+10], 1-p[add+10]); add <- 4
        } else if (model[i]==4) {
            lapses <- c(lapses, 0, p[add+7], p[add+8], 0, 0, p[add+9], 1-p[add+9], p[add+10], 1-p[add+10]); add <- 4
        } else if (model[i]==5) {
            lapses <- c(lapses, p[add+7], 0, 0, p[add+8], 1-p[add+8], 0, 0, 0, 0); add <- 2
        } else if (model[i]==6) {
            lapses <- c(lapses, 0, p[add+7], 0, 0, 0, p[add+8], 1-p[add+8], 0, 0); add <- 2
        } else if (model[i]==7) {
            lapses <- c(lapses, 0, 0, p[add+7], 0, 0, 0, 0, p[add+8], 1-p[add+8]); add <- 2
        }
    }
    if (logs[3]) {
        lapses[2] <- 1; lapses[11] <- 1
        lapses[4] <- 0; lapses[5] <- 1; lapses[13] <- 0; lapses[14] <- 1
        lapses[8] <- 1; lapses[9] <- 0; lapses[17] <- 1; lapses[18] <- 0
    }
    if (logs[4]) {
        lapses[1] <- 1; lapses[10] <- 1
        lapses[6] <- 0; lapses[7] <- 1; lapses[15] <- 0; lapses[16] <- 1
        lapses[8] <- 0; lapses[9] <- 1; lapses[17] <- 0; lapses[18] <- 1
    }
    all.pars <- c(alpha_st, alpha_t, Beta_st, beta_t, delta_1, delta_2, lapses)
    return(all.pars)
}

#### statistical Utility Functions ####

#' Psychophysical function
mu <- function(x, a, b) {
    x1 <- (x-a)/b
    if (x1<700) {
        return(log(1+2*exp(x1)))
    }
    return(log(2)+x1)
    #y <- colSums(rbind(x2*log(1+2*exp(x1)), (1-x2)*(log(2)+x1)),na.rm=T)
    #return(y)
}


#' Inverse psychophysical function
muinv <- function(a, b, y) {
    if (y<700){
        return(a+b*log((exp(y)-1)/2))
    }
    return((y-log(2))*b+a)
    #x <- colSums(rbind(y1*(a+b*log((exp(y)-1)/2)), (1-y1)*((y-log(2))*b+a)), na.rm=T)
    #return(x)
}


#' Psychometric functions
#' @importFrom stats pnorm
Psy <- function(st, t, a_st, b_st, a_t, b_t, d_1, d_2, l) {
    mu_t <- mu(t,a_t,b_t)
    mu_st <- mu(st,a_st,b_st)
    p <- matrix(NA, nrow=6, ncol=length(t))
    # test first
    l1 <- l[1]; lU <- l[2]; l2 <- l[3]
    k1U <- l[4]; k12 <- l[5]; kU1 <- l[6]; kU2 <- l[7]; k21 <- l[8]; k2U <- l[9]
    mudiff <- mu_st - mu_t
    p_l <- pnorm(d_1,mudiff,sqrt(2),lower.tail=T)
    p_r <- pnorm(d_2,mudiff,sqrt(2),lower.tail=F)
    p_c <- 1-p_l-p_r
    p_l <- pmin(pmax(p_l, 0), 1)
    p_r <- pmin(pmax(p_r, 0), 1)
    p_c <- pmin(pmax(p_c, 0), 1)
    p[1,] <- (1-l1)*p_l + lU*kU1*p_c + l2*k21*p_r
    p[2,] <- l1*k1U*p_l + (1-lU)*p_c + l2*k2U*p_r
    p[3,] <- l1*k12*p_l + lU*kU2*p_c + (1-l2)*p_r
    # test second
    l1 <- l[10]; lU <- l[11]; l2 <- l[12]
    k1U <- l[13]; k12 <- l[14]; kU1 <- l[15]; kU2 <- l[16]; k21 <- l[17]; k2U <- l[18]
    mudiff <- mu_t - mu_st
    p_l <- pnorm(d_1,mudiff,sqrt(2),lower.tail=T)
    p_r <- pnorm(d_2,mudiff,sqrt(2),lower.tail=F)
    p_c <- 1-p_l-p_r
    p_l <- pmin(pmax(p_l, 0), 1)
    p_r <- pmin(pmax(p_r, 0), 1)
    p_c <- pmin(pmax(p_c, 0), 1)
    p[4,] <- (1-l1)*p_l + lU*kU1*p_c + l2*k21*p_r
    p[5,] <- l1*k1U*p_l + (1-lU)*p_c + l2*k2U*p_r
    p[6,] <- l1*k12*p_l + lU*kU2*p_c + (1-l2)*p_r
    return(p)
}

#' Likelihood function
LogLikelihood <- function(x, F1, U1, S1, F2, U2, S2, standard, Lev, models, logs, bounds) {
    all.pars <- GetParVec(x, models, logs, bounds)
    lapses <- all.pars[-1:-6]
    if (any(lapses<0) || any (lapses>1)) {
        L <- 1.0e10
    } else if (any(!is.finite(all.pars)) || any (is.complex(all.pars))) {
        L <- 1.0e10
    } else {
        p.output <- Psy(standard, Lev, all.pars[1], all.pars[3], all.pars[2], all.pars[4],
                        all.pars[5], all.pars[6], lapses)
        P_F1 <- p.output[1,]; P_U1 <- p.output[2,]; P_S1 <- p.output[3,]
        P_F2 <- p.output[4,]; P_U2 <- p.output[5,]; P_S2 <- p.output[6,]
        L <- -sum(c(F1*log(P_F1), U1*log(P_U1), S1*log(P_S1),
                    F2*log(P_F2), U2*log(P_U2), S2*log(P_S2)), na.rm=T)
        if (!is.finite(L) || is.complex (L)) {L <- 1.0e10}
    }
    return(L)
}
