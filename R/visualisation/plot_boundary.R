Delta_1 = -1.32; Delta_2 = 2.91
xsup2 <- max(abs(c(8, Delta_1, Delta_2))); # take absmax of decision bounds for xlims
xsup2 <- ceiling(xsup2); # round up this value to whole number
xinf2 <- -xsup2 # set negative version for lower bound
step2 <- (xsup2-xinf2)/800 # abs sum of both over 800 to define intervals
xval2 <- seq(xinf2,xsup2,step2) # set 800 steps between limits
yval <- dnorm(xval2,0,sqrt(2)) # normal distribution (u=0,sd=2^0.5)

x = c(xinf2,xsup2)
y = c(0, 1.4*dnorm(0,0,sqrt(2)))
step2 = step2,

data <- data.frame(
    xval2 = xval2,
    yval = yval
)


ggplot(data, aes(xval2,yval)) +
    geom_line() +
    geom_vline(xintercept = Delta_1) +
    geom_vline(xintercept = Delta_2) +
    geom_segment(aes(x=0 , y=0, xend=0, yend=max(yval)), linetype=2) +
    xlab('Decision Variable') +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank())


    geom_text(mapping = aes(x = vals,
                            y = 0,
                            label = Ref,
                            hjust = -1,
                            vjust = -1),
              data = cuts)



plot(c(xinf2,xsup2),c(0, 1.4*dnorm(0,0,sqrt(2))), type="n", adj=0.5, axes=F,cex=0.9,
     xlab='Decision variable', ylab='',
     main='Decision boundaries')
lines (xval2,yval,type='l',col='black',lwd=lw)
lines (c(0,0),c(0,dnorm(0,0,sqrt(2))),lty=3,col='black')
yval2 <- 1.2*dnorm(0,0,sqrt(2))
lines(c(Delta_1,Delta_1),c(0,yval2),lty=2, col='black', lwd=lw/2)
lines(c(Delta_2,Delta_2),c(0,yval2),lty=2, col='black', lwd=lw/2)
lines(c(xinf2,xsup2),c(0,0),lty=1, col='black', lwd=lw/2)
text(Delta_1, yval2, pos=2,cex=0.75,
     labels=substitute(paste(delta[1],"=",label),list(label=Delta_1)))
text(Delta_2, yval2, pos=4, cex=0.75,
     labels=substitute(paste(delta[2],"=",label),list(label=Delta_2)))
text(0,0,pos=1,labels='0')


p <- ggplot()
