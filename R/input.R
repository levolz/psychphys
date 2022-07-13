#' Reading In Raw Data
#'
#' @description  Convenience function for reading in raw data (only certain formats)
#'
#' @param file name of file to be cleaned
#'
#' @param format raw data format - must be one of c("AMOCS")
#'
#' @param folder optional folder in which data file is located
#'
#' @returns Data matrix with fitting structure to be used in \code{fit_TIM_2P()}.
#'
#' @export
data_prep <- function (file, format, folder) {
    if (!format %in% c("AMOCS")){
        stop("No valid raw data format supplied")
    }
    if (missing(file)) {
        stop("Missing file name")
    } else if (!missing(folder)){
        f.name = paste(folder, "/", file, sep = "")
    } else {
        f.name <- paste(file, sep = "")
    }
    if (format == "AMOCS"){
        raw <- read.table(f.name,
                          col.names=c("staircase", "trial", "standard", "test", "order", "response"))
        tabled = table(raw[,c("test", "order","response")])
        dim_tabled = dim(tabled)
        clean <- matrix(as.double(tabled), c(dim_tabled[2] * dim_tabled[3], dim_tabled[1]), byrow=T)
        # std <- unique(raw$standard)
        rearrange <- c(3,1,5,4,2,6) # To get the right order: F1, U1, S1, F2, U2, S2
        clean <- clean[rearrange,]
        rownames (clean) <- c(" F1", " U1", " S1"," F2", " U2", " S2")
        Lev <- sort(unique(raw$test))
        freq <-rbind(Lev, clean)
    }
    return (freq)
}


## #' Simulate Data from Underlying Model
## #'
## #' @export
## simulate_data = function(){
##
## }
