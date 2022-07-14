#' Print the results
#'
#' @importFrom utils read.table
#'
#' @export
print.fit <- function(output){
    cat ("----------------------------------------------------------------", "\n")
    cat (output$problem,"\n")
    cat ("----------------------------------------------------------------", "\n")
    rownames (output$data) <- c("Lev", " F1", " U1", " S1"," F2", " U2", " S2")
    colnames (output$data) <- rep("",times=dim(output$data)[2])
    print.table (output$data); cat("   ")
    cat (paste (names(output[-c(1,3)]),":", " ",output[-c(1,3)], "\n", sep=""))
    invisible (output)
}

#' Plot the results
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics lines par points text title
#' @importFrom stats dnorm
#'
#' @export
plot_results <- function(pars, standard, Lev, F1, U1, S1, F2, U2, S2, PSE, Thr_84, logs) {
    alpha_st <- pars[1]; alpha_t <- pars[2]; beta_st <- pars[3]; beta_t <- pars[4]
    delta_1 <- pars[5]; delta_2 <- pars[6]; lapses <- pars[7:24]
    same <- logs[1]; ter <- logs[2]; bin <- logs[3]; equ <- logs[4]; det <- logs[5]
    xinf <- min(Lev);  xsup <- max(Lev)
    step <- (xsup-xinf)/800
    xval <- seq(xinf,xsup,step)
    psycho <- Psy(standard, xval, alpha_st, beta_st, alpha_t, beta_t, delta_1, delta_2, lapses)
    P_F1 <- psycho[1,]; P_U1 <- psycho[2,]; P_S1 <- psycho[3,]
    P_F2 <- psycho[4,]; P_U2 <- psycho[5,]; P_S2 <- psycho[6,]
    #windows(5,5, rescale='R') # Only works under windows OS
    dev.new(width=5, height=4)
    par(mai=c(0.8,0.5,0.8,0.1), cex=1.05, mex=0.8, xaxs='i', yaxs='i', tck=-0.03, xpd=T,
        family='sans', pty='s',las=1, pch=19)
    lw <- 3 # line width#
    n_1 <- F1 + U1 + S1
    n_2 <- F2 + U2 + S2
    plot(c(xinf,xsup), c(0,1), type="n",
         xlab="Stimulus level", ylab="Probability of response", main = "")
    if (ter) {
        title(main =
                  expression(bold(paste("Test ", 1^st,": ", phantom(paste("F U S   Test ", 2^nd,": ", "F U S")))))
              , cex=0.9)
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": ")), "F", phantom(paste("U S   Test ", 2^nd,": F U S")))))
              , cex=0.9, col.main="blue")
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": F")), "U", phantom(paste("S   Test ", 2^nd,": F U S")))))
              , cex=0.9, col.main="black")
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": F U")), "S", phantom(paste("   Test ", 2^nd,": F U S")))))
              , cex=0.9, col.main="red")
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": F U S")), "   Test ", 2^nd,": ", phantom("F U S"))))
              , cex=0.9, col.main="black")
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": F U S   Test ",2^nd,": ")),"F", phantom("U S"))))
              , cex=0.9, col.main="pink")
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": F U S   Test ",2^nd,": F")),"U", phantom("S"))))
              , cex=0.9, col.main="gray")
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": F U S   Test ",2^nd,": F U")),"S")))
              , cex=0.9, col.main="cyan")
    } else if (bin) {
        title(main =
                  expression(bold(paste("Test ", 1^st,": ", phantom(paste("F   Test ", 2^nd,": S")))))
              , cex=0.9)
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": ")), "F", phantom(paste("   Test ", 2^nd,": S")))))
              , cex=0.9, col.main="blue")
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": F")), "   Test ", 2^nd,": ", phantom("S"))))
              , cex=0.9, col.main="black")
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": F   Test ", 2^nd,": ")), "S")))
              , cex=0.9, col.main="cyan")
    } else if (equ) {
        title(main =
                  expression(bold(paste("Test ", 1^st,": same   Test ", 2^nd,": ", phantom("same"))))
              , cex=0.9)
        title(main =
                  expression(bold(paste(phantom(paste("Test ", 1^st,": same   Test ", 2^nd,": ")), "same")))
              , cex=0.9, col.main="gray")
    }
    if (!equ) {
        lines(xval, P_F1, type='l', col='blue', lwd=lw)
        lines(xval, P_S2, type='l', col='cyan', lwd=lw)
        points(Lev, F1/n_1, type='p', col='blue', lwd=lw)
        points(Lev, S2/n_2, type='p', col='cyan', lwd=lw)
        if (ter) {
            lines(xval, P_S1, col='red', lwd=lw)
            lines(xval, P_F2, col='pink', lwd=lw)
            points(Lev, S1/n_1, type='p', col='red', lwd=lw)
            points(Lev, F2/n_2, type='p', col='pink', lwd=lw)
        }
    }
    if (!bin) {
        lines(xval, P_U1, type='l', col='black', lwd=lw)
        lines(xval, P_U2, type='l', col='gray', lwd=lw)
        points(Lev, U1/n_1, type='p', col='black', lwd=lw)
        points(Lev, U2/n_2, type='p', col='gray', lwd=lw)
    }
    if (!det) {
        lines(c(standard, standard),c(0, 1),lty=3, col='black')
        if (is.numeric(PSE)) {lines(c(PSE, PSE),c(0, 1),col='black')
        } else {
            lines(c(Thr_84, Thr_84),c(0, 1),col='black')}
    }

    #### second figure
    #Part 1
    #windows(18,6, rescale='R') # Only works under windows OS
    dev.new(width=18, height=6)
    par(mfrow=c(1,3), mai=c(0.8,0.5,0.8,0.1), cex=1.05, mex=0.8, xaxs='i', yaxs='i', tck=-0.03, xpd=T,
        family='sans', pty='s',las=1, pch=16)
    mu_st <- mu(standard,alpha_st,beta_st)
    yval <- mu(xval,alpha_t,beta_t)
    xy <- c(xinf, xsup, min(yval), max(yval))
    plot(xy[1:2], xy[3:4], type="n", adj=0.5, axes=T, cex=0.9,
         xlab='Stimulus level', ylab='Subjective level',
         main='Psychophysical function')
    lines (xval, yval,type='l',col='black', lwd=lw/2)
    if (!det) {
        lines(standard, mu_st, type='p', col='blue', lwd=lw)
        text(standard,mu_st,labels='   standard', pos=4, offset=0.25,cex=0.75)
        lines(c(xinf, standard, standard),c(mu_st, mu_st, xy[3]),lty=3,col='black')
        if (is.numeric(PSE) && !same) {
            lines(c(xinf, PSE, PSE),c(mu_st, mu_st, xy[3]),lty=3,col='black')
            text(PSE,xy[3]+0.08*(xy[4]-xy[3]),labels=paste(' PSE = ',round(PSE, digits=3)),
                 pos=4, offset=0.25, cex=0.75)
        }
    } else {
        mu_thr <- mu(Thr_84,alpha_t,beta_t)
        lines(c(xinf, Thr_84, Thr_84),c(mu_thr, mu_thr, xy[3]),lty=3,col='black')
        text(Thr_84,xy[3]+0.08*(xy[4]-xy[3]),
             labels=substitute(paste(theta,"=",label),list(label=Thr_84)),
             pos=4, offset=0.5, cex=0.75)
    }
    xpos <- xinf+0.05*(xsup-xinf)
    ypos <-  xy[3]+0.9*(xy[4]-xy[3])
    text(xpos,ypos,labels=substitute(paste(alpha[t],"=",label),list(label=alpha_t)), pos=4, cex=0.75)
    xpos <- xinf+0.05*(xsup-xinf)
    ypos <- xy[3]+0.8*(xy[4]-xy[3])
    text(xpos,ypos,labels=substitute(paste(beta[t],"=",label),list(label=beta_t)),pos=4, cex=0.75)
    #Part 2
    xsup2 <- max(abs(c(8, delta_1, delta_2))); xsup2 <- ceiling(xsup2); xinf2 <- -xsup2
    step2 <- (xsup2-xinf2)/800
    xval2 <- seq(xinf2,xsup2,step2)
    yval <- dnorm(xval2,0,sqrt(2))
    plot(c(xinf2,xsup2),c(0, 1.4*dnorm(0,0,sqrt(2))), type="n", adj=0.5, axes=F,cex=0.9,
         xlab='Decision variable', ylab='',
         main='Decision boundaries')
    lines (xval2,yval,type='l',col='black',lwd=lw)
    lines (c(0,0),c(0,dnorm(0,0,sqrt(2))),lty=3,col='black')
    yval2 <- 1.2*dnorm(0,0,sqrt(2))
    lines(c(delta_1,delta_1),c(0,yval2),lty=2, col='black', lwd=lw/2)
    lines(c(delta_2,delta_2),c(0,yval2),lty=2, col='black', lwd=lw/2)
    lines(c(xinf2,xsup2),c(0,0),lty=1, col='black', lwd=lw/2)
    text(delta_1, yval2, pos=2,cex=0.75,
         labels=substitute(paste(delta[1],"=",label),list(label=delta_1)))
    text(delta_2, yval2, pos=4, cex=0.75,
         labels=substitute(paste(delta[2],"=",label),list(label=delta_2)))
    text(0,0,pos=1,labels='0')
    #Part 3
    psycho_0 <- Psy(standard, xval, alpha_st, beta_st, alpha_t, beta_t, delta_1, delta_2, lapses*0)
    plot(xy[1:2], c(0,1), type="n", cex=0.9,
         xlab='Stimulus level', ylab='Probability of response',
         main=expression(bold(paste('Ternary functions with all ', epsilon, ' = 0'))))
    lines(xval,psycho_0[1,],col='blue', type='l', lwd=lw)
    lines(xval,psycho_0[2,],col='black', type='l',lwd=lw)
    lines(xval,psycho_0[3,],col='red', type='l', lwd=lw)
    lines(xval,psycho_0[4,],col='pink', type='l', lwd=lw)#[1 0.4 0.6]
    lines(xval,psycho_0[5,],col='gray', type='l', lwd=lw)#[0.5 0.5 0.5]
    lines(xval,psycho_0[6,],col='cyan', type='l', lwd=lw)
    if (!det) {
        lines(c(standard, standard),c(0,1),col='black', lty=3)
        if (is.numeric(PSE)) lines(c(PSE, PSE),c(0, 1),col='black')
    } else {
        lines(c(Thr_84, Thr_84),c(0,1),col='black')
    }
}
