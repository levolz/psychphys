# plot psychophysical function
plot_psychphys <- function(pars, Standard, Lev, F1, U1, S1,
                           F2, U2, S2, PSE, Thr_84, logs){


}

# plot decision boundary function
plot_boundary <- function(pars, Standard, Lev, F1, U1, S1,
                           F2, U2, S2, PSE, Thr_84, logs){


}

# plot psychometric
plot_psychmet <- function(pars, Standard, Lev, F1, U1, S1,
                          F2, U2, S2, PSE, Thr_84, logs){


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
