# plot psychophysical function
plot_psychphys <- function(){

    #alpha_st, beta_st, lapses,

}

# plot decision boundary function
plot_boundary <- function(model){

    d1 = model$Delta_1; d2 = model$Delta_2
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

# plot psychometric
plot_psychmet <- function(params, standard, lev, f1, u1, s1,
                          f2, u2, s2, pse, thr_84, logs){


}


# wrapper function
plot_results <- function(plot_id = 1:4){

    if(1 %in% plot_id){
        p1 <- plot_psychphys()
    }
    if(2 %in% plot_id){
        p2 <- plot_boundary()
    }
    if(3 %in% plot_id){
        p3 <- plot_psychmet()
    }

}
