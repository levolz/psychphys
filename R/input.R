#' Data reading from raw data
#'
#' @param code
#'
#' @returns Data matrix with fitting structure to be used in \code{fit_TIM_2P()}.
#'
#' @export
data.prep <- function (code) {
    f.name <- paste('data/', code, '_AMOCS.txt', sep="")
    raw <- read.table(f.name,
                      col.names=c("staircase", "trial", "standard", "test", "order", "response"))
    clean <- matrix(as.double(table(raw[,c("test", "order","response")])), c(6, 10), byrow=T)
    std <- unique(raw$standard)
    rearrange <- c(3,1,5,4,2,6) # To get the right order: F1, U1, S1, F2, U2, S2
    clean <- clean[rearrange,]
    rownames (clean) <- c(" F1", " U1", " S1"," F2", " U2", " S2")
    Lev <- sort(unique(raw$test))
    freq <-rbind(Lev, clean)
    return (freq)
}
