#' plot psychophysical function
plot_psychphys <- function(model){
    alpha_st <- model$alpha_st
    alpha_t <- model$alpha_t
    beta_st <- model$beta_st
    beta_t <- model$beta_t
    same <- model$same
    det <- model$det
    xinf <- min(model$data[1,])
    xsup <- max(model$data[1,])
    thr_84 <- model$Threshold_84
    standard <- model$standard
    pse <- model$PSE

    step <- (xsup-xinf)/800
    xval <- seq(xinf,xsup,step)
    mu_st <- mu(standard, alpha_st, beta_st)
    yval <- mu(xval, alpha_t, beta_t)
    xy <- c(xinf, xsup, min(yval), max(yval))

    plot_data = data.frame(x = xval, y = yval)

    psychphys <- ggplot(data = plot_data, aes(x = x, y = y)) +
        geom_line() +
        theme_minimal() +
        ggtitle("Psychophysical function") +
        xlab("Stimulus level") +
        ylab("Subjective level")

    if(det){
        mu_thr <- mu(thr_84, alpha_t, beta_t)
        psychphys <- psychphys +
            geom_segment(aes(x = xinf, y = mu_thr,
                             xend = thr_84, yend = mu_thr,
                             linetype = 3)) +
            geom_segment(aes(x = thr_84, y = mu_thr,
                             xend = thr_84, yend = y[1],
                             linetype = 3)) +
            geom_point(aes(x = thr_84, y = mu_thr), color = "blue") +
            scale_linetype_identity() +
            annotate("text", label = paste("\u03B1 = ", round(alpha_t, 3),
                                           "\n\u03B2 = ", round(beta_t, 3),
                                           "\n\u03B8 = ", round(thr_84, 3), sep = ""),
                     y = (xy[4] - 0.05 * (xy[4]-xy[3])),
                     x = (xy[2] - 0.975 * (xy[2]-xy[1])),
                     hjust = 0)
    } else {
        psychphys <- psychphys + #add 'Standard' text
            geom_segment(aes(x = xinf, y = mu_st,
                             xend = standard, yend = mu_st,
                             linetype = 3)) +
            geom_segment(aes(x = standard, y = mu_st,
                             xend = standard, yend = y[1],
                             linetype = 3)) +
            annotate("text", label = "Standard",
                     x = standard, y = mu_st + 0.025*(xy[4]-xy[3])) +
            geom_point(aes(x = standard, y = mu_st), color = "blue") +
            scale_linetype_identity()

        if(is.numeric(pse) & !same){
            psychphys <- psychphys +
                geom_segment(aes(x = pse, y = mu_st,
                                 xend = pse, yend = y[1],
                                 linetype = 3)) +
                annotate("text", label = "PSE",
                         x = pse - 0.025*(xy[2]-xy[1]), y = mu_st + 0.025*(xy[4]-xy[3])) +
                geom_point(aes(x = pse, y = mu_st), color = "red") +
                annotate("text", label = paste("\u03B2 = ", round(beta_t, 3),
                                               "\nPOE - PSE = ", round(standard - pse, 3), sep=""),
                         y = (xy[4] - 0.05 * (xy[4]-xy[3])),
                         x = (xy[2] - 0.975 * (xy[2]-xy[1])),
                         hjust = 0)
        } else {
            psychphys <- psychphys +
            annotate("text", label = paste("\n\u03B2 = ", round(beta_t, 3)),
                     y = (xy[4] - 0.05 * (xy[4]-xy[3])),
                     x = (xy[2] - 0.975 * (xy[2]-xy[1])),
                     hjust = 0)
        }
    }

    return(psychphys)
}

#' plot decision boundary function
#' @import ggplot
plot_boundary <- function(model){

    d1 = model$delta_1; d2 = model$delta_2
    xsup <- max(abs(c(8, d1, d2))); # take absmax of decision bounds for xlims
    xsup <- ceiling(xsup); # round up this value to whole number
    xinf <- -xsup # set negative version for lower bound
    step <- (xsup-xinf)/800 # abs sum of both over 800 to define intervals
    x <- seq(xinf,xsup,step) # set 800 steps between limits
    y <- dnorm(x,0,sqrt(2)) # normal distribution (u=0,sd=2^0.5)
    data <- data.frame(x=x, y=y)

    p <- ggplot(data, aes(x,y)) +
        geom_line() +
        geom_vline(xintercept = d1) +
        geom_vline(xintercept = d2) +
        geom_segment(aes(x=0 , y=0, xend=0, yend=max(y)), linetype=2) +
        xlab('Decision Variable') +
        theme_minimal() +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank())

    return(p)
}

#' plot psychometric
plot_psychmet <- function(model){


}


#' plotting wrapper function
#' @import ggplot
#' @importFrom cowplot plot_grid
plot_results <- function(model,to_plot = 1:3){
    if(1 %in% to_plot){
        p1 <- plot_psychphys(model)
    }
    if(2 %in% to_plot){
        p2 <- plot_boundary(model)
    }
    if(3 %in% to_plot){
        p3 <- plot_psychmet(model)
    }

    cowplot::plot_grid(plotlist=mget(paste0("p", to_plot)))
}
