\name{j2}
\alias{j2}

\title{
  MCMC estimates for the j2 model
}
\description{
  Estimates j2 model parameters as described in Zijlstra (2017) <doi:10.1080/0022250X.2017.1387858>.
}
\usage{
j2(net, sender = NULL, receiver = NULL , density = NULL, reciprocity = NULL, 
burnin = NULL, sample = NULL, adapt= NULL, center = NULL, seed = NULL) 
}

\arguments{
  \item{net}{Directed dichotomous n*n network (digraph).}
  \item{sender}{Optional sender covariates of lenght n.}
  \item{receiver}{Optinal receiver covariates of length n.}
  \item{density}{Optional density covariates of dimensions n*n.}
  \item{reciprocity}{Optional symmetric reciprocity covariates of dimensions n*n.}
  \item{burnin}{Optional specification of number of burn-in iterations (default is 10000).}
  \item{sample}{Optional specification of number of MCMC samples (default is 40000).}
  \item{adapt}{Optional number of adaptive sequenses (default is 100).}
  \item{center}{Optional boolean argument for centering predictors (default is TRUE).}
  \item{seed}{Optional specification of random seed (delfault is 1).}
}

\value{
Returns a matrix with MCMC means, standard deviations, quantiles and effective sample sizes for j2 parameters.
}
\references{
Zijlstra, B.J.H. (2017). Regression of directed graphs on independent effects for density and reciprocity. \emph{The Journal of Mathematical Sociology} 41 (4).
}
\author{
Bonne J.H. Zijlstra \email{b.j.h.zijlstra@uva.nl}
}

\examples{

# create a very small network with covariates for illustrative purposes
S <- c(1,0,1,0,1,1,0,1,0,1)
REC <- c(0,0,1,1,0,0,1,1,0,0)
D1 <- matrix(c(0,1,0,1,0,1,0,1,0,0,
               0,0,1,1,0,1,0,1,0,1,
               1,1,0,0,1,0,0,0,0,0,
               1,1,1,0,1,0,0,0,0,1,
               1,0,1,0,0,1,1,0,1,1,
               0,0,0,0,0,0,1,1,1,1,
               0,0,0,0,0,1,0,1,0,1,
               1,0,0,0,0,1,1,0,1,1,
               0,1,0,1,0,1,0,1,0,0,
               0,0,1,1,1,0,0,0,0,0), ncol=10)
D2 <- abs(matrix(rep(S,10), byrow = FALSE, ncol= 10) -
            matrix(rep(REC,10), byrow = TRUE, ncol= 10))
R <- D1*t(D1)
Y <- matrix(c(0,0,1,1,1,1,0,0,1,1,
              0,0,0,1,1,1,0,0,1,0,
              1,1,0,1,1,1,0,0,1,1,
              0,1,1,0,1,1,0,1,1,0,
              1,1,1,1,0,1,1,0,1,1,
              0,1,1,1,1,0,1,1,1,0,
              1,0,1,0,1,1,0,1,0,1,
              0,1,1,1,0,1,1,0,1,1,
              1,0,1,0,1,0,1,1,0,1,
              1,1,1,0,0,1,1,1,1,0), ncol=10) 

# estimate j2 model
j2(Y,sender= ~ S, receiver =  ~ REC, density = ~ D1 + D2, reciprocity= ~ R,
   burnin = 100, sample = 400, adapt = 10)
# notice: burn-in, sample size and number of adaptive sequenses are 
# much smaller than recommended to keep computation time low.
# recommended code: 
\dontrun{
j2(Y,sender= ~ S, receiver =  ~ REC, density = ~ D1 + D2, reciprocity= ~ R)
}
}
\keyword{j2 model}

