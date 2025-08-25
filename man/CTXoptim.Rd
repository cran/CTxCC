\name{CTXoptim}
\alias{CTXoptim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimum Design of a CTX Control Chart
}
\description{
Computes the optimal weight vector and control limit for a CTX control chart that minimizes the overall out-of-control Average Run Length (ARL) for a given set of shift vectors.
}
\usage{
CTXoptim(delta, R, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
 list of numeric vectors, each representing a mean shift scenario to be optimized for (out-of-control cases).
}
  \item{R}{
A positive-definite correlation matrix for the monitored variables.
}
  \item{alpha}{
Target in-control false alarm rate (e.g., 1/500) corresponds to an in-control ARL of 500).
}
}
\details{
The optimization is performed over proxy weights (log-transformed)
to ensure positivity. These weights are normalized to sum to one,
and a quasi-Newton method (`BFGS`) is used to minimize the sum
of out-of-control ARLs computed for each specified shift scenario.
The control limit is computed to satisfy the desired in-control ARL
(via \code{alpha} = 1/ARL0).
}
\value{
  \item{weights}{Numeric vector of optimal weights that sum to one.}
  \item{control_limit}{Numeric value of the CTX control chart limit that satisfies the specified alpha}.

}
\references{
Paper (pending complete bibliographic reference when available).

}
\author{
Dr. Burcu Aytaçoğlu (\email{burcuaytacoglu@gmail.com}),
Dr. Diana Barraza-Barraza (\email{diana.barraza@ujed.mx}),
Dr. Víctor G. Tercero-Gómez (\email{victor.tercero@tec.mx}),
Dr. A. Eduardo Cordero-Franco (\email{lalo.cordero@gmail.com})
}

\examples{
\dontrun{
# Define shift scenarios
delta1 <- c(1, 0, 0)
delta2 <- c(0, 3, 0)
delta3 <- c(0, 0, 2)
delta <- list(delta1, delta2, delta3)

 # Correlation matrix and false alarm rate
 R <- diag(3)
 alpha <- 1/500

 # Compute optimal weights and control limit
 CTXoptim(delta = delta, R = R, alpha = alpha)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
