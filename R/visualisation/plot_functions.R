# plot psychophysical function
plot_psychphys <- function(pars, standard, lev, f1, u1, s1,
                           f2, u2, s2, pse, thr_84, logs){


}

# plot decision boundary function
plot_boundary <- function(pars, standard, lev, f1, u1, s1,
                           f2, u2, s2, pse, thr_84, logs){


}

# plot psychometric
plot_psychmet <- function(pars, standard, lev, f1, u1, s1,
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
