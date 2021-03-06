% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{fit_TIM_2P}
\alias{fit_TIM_2P}
\title{Fit the ternary indecision model to detection or discrimination data from
a dual-presentation task}
\usage{
fit_TIM_2P(
  Data,
  Standard,
  Format,
  AlphaBounds,
  BetaBounds,
  Delta1Bounds,
  WidthBounds,
  AlphaStart,
  BetaStart,
  Delta1Start,
  WidthStart,
  EpsStart,
  KappaStart,
  Model = "best",
  Criterion = "LogL",
  Type,
  Plot = FALSE,
  Disp = FALSE
)
}
\arguments{
\item{Data}{7-by-N array with data from the dual-presentaion task. The
first row gives the N levels of the test stimulus; rows 2-4 respectively
give the count of "first", "undecided", and "second" responses at each
level when the test stimulus was presented in the first interval/position;
rows 5-7 analogously give the count of "first", "undecided", and "second"
responses at each level when the test stimulus was presented in the second
interval/position. If \code{Format <- "2AFC"} or
\code{Format <- "equality"}, array Data still has the same number of rows
but its content is subject to the constraints listed under argument Format
below.}

\item{Standard}{Level of the standard stimulus. For detection data (i.e.,
with a nominally null standard), set \code{Standard <- -Inf}; for
discrimination data (i.e., with a non-null standard), set Standard to the
appropriate level, in the same units in which levels of the test stimulus
are given in the first row of array Data.}

\item{Format}{Out of \code{c("ternary", "2AFC", or "equality")} to indicate
that data were collected with a ternary response format (first, second, or
undecided) a 2AFC response format (forcing first or second responses by
guessing when uncertain), or the same-different response format (both stimuli
are subjective equal or they are subjectively different).
When \code{Format <- "2AFC"}, rows 3 and 6 of array Data must be filled with
zeros (as no "undecided" responses are given under this format);
When \code{Format <- "equality"}, rows 2 and 5 of array Data must be filled with zeros
(so "different" responses are stored as "second" responses).}

\item{AlphaBounds}{2-element vector with the lower and upper bounds on Alpha_t. Reals.}

\item{BetaBounds}{2-element vector with the lower and upper bounds on Beta_t. Non-negative reals.}

\item{Delta1Bounds}{2-element vector with the lower and upper bounds on Delta_1. Reals.}

\item{WidthBounds}{2-element vector with the lower and upper bounds on the width of the central
region in decision space (i.e., Delta_2 - Delta_1). Reals equal to or greater than 0.}

\item{AlphaStart}{Starting value(s) for Alpha_t. Scalar or vector of reals.}

\item{BetaStart}{Starting value(s) for Beta_t. Scalar or vector of non-negative reals.}

\item{Delta1Start}{Starting value(s) for Delta_1. Scalar or vector of reals.}

\item{WidthStart}{Starting value(s) for width. Scalar or vector of non-negative reals.}

\item{EpsStart}{Starting value(s) for epsilon parameters. Scalar or vector of reals in \eqn{[0, 1]}.}

\item{KappaStart}{Starting value(s) for kappa parameters. Scalar or vector of reals in \eqn{[0, 1]}.}

\item{Model}{Choice of model to be fitted. An integer scalar in the range between 0 and 7
implying different assumptions about error parameters (see the table below)
that will be fitted to both presentation orders (i.e., test first and test
second), a 2-element vector with components in the same range (describing the
error model to be fitted to each presentation order), or the string "best"
(case insensitive) to find the best-fitting error model for each presentation
order by the criterion of choice (see the next input argument).
When \code{Format <- "2AFC"}, models 0, 3, 5, and 7 are not permitted because
e_U must be in the fitted model;
when \code{Format <- "equality"}, models 0, 4, 6, and 7 are not permitted
because e_F must be in the fitted model.
Setting \code{Model <- "best"} takes these constraints into account when
searching for the best-fitting model.

\tabular{rccc}{
\strong{Model} \tab \strong{e_F} \tab \strong{e_U} \tab \strong{e_S} \cr
0 \tab x \tab x \tab x \cr
1 \tab o \tab o \tab o \cr
2 \tab o \tab o \tab x \cr
3 \tab o \tab x \tab o \cr
4 \tab x \tab o \tab o \cr
5 \tab o \tab x \tab x \cr
6 \tab x \tab o \tab x \cr
7 \tab x \tab x \tab o \cr
}
x : parameter is not included in the model
o : parameter is included in the model}

\item{Criterion}{Choice of criterion when \code{Model <- "best"}.
Options are \code{"LogL"} (to use the value of the negative log-likelihood)
or \code{"BIC"} (to use the value of the Bayesian Information Criterion).
Not used when argument Model is not the string \code{"best"}.}

\item{Type}{Choice of \code{c("same")} (when the same psychophysical function holds
for test and standard stimuli) or \code{c("diff")} (when different psychophysical
functions are assumed to hold for test and standard stimuli).
Not used for detection data (i.e., when \code{Standard <- -Inf}).}

\item{Plot}{Option to plot data and fitted functions. Logical scalar.}

\item{Disp}{Option to issue warnings or display progress information. Logical scalar.}
}
\value{
List of includes parameter estimates and goodness-of-fit measures and p-values.
}
\description{
Further details:  Garc?a-P?rez, M.A. & Alcal?-Quintana, R. (2017).
The indecision model of psychophysical performance in
dual-presentation tasks: Parameter estimation and
comparative analysis of response formats.
Frontiers in Psychology, 8:1142. \url{https://doi.org/10.3389/fpsyg.2017.01142}
}
