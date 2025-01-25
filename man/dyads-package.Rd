\name{dyads-package}
\alias{dyads-package}
\alias{dyads}
\docType{package}
\title{
dyads
}
\description{
Package for Dyadic Network Analysis.
}
\details{


\packageDESCRIPTION{dyads}
\packageIndices{dyads}
Includes functions for estimation of the (multilevel) p2 model (van Duijn, Snijders and Zijlstra (2004) <doi:10.1046/j.0039-0402.2003.00258.x>), more specifically the adaptive random walk algorithm (Zijlstra, van Duijn and Snijders (2009) <doi:10.1348/000711007X255336>), for the estimation of the j2 model (Zijlstra (2017) <doi:10.1080/0022250X.2017.1387858>), and for their bidirectional counterpart, b2.
}
\author{
Bonne J.H. Zijlstra
Maintainer: Bonne J.H. Zijlstra <B.J.H.Zijlstra@uva.nl>
}
\references{
Zijlstra, B.J.H., Duijn, M.A.J. van, and  Snijders, T.A.B. (2009). MCMC estimation for the $p_2$ network regression model with crossed random effects. \emph{British Journal of Mathematical and Statistical Psychology}, 62, 143-166.
Zijlstra, B.J.H. (2017). Regression of directed graphs on independent effects for density and reciprocity. \emph{Journal of Mathematical Sociology}, 41(4), 185-192.
}


\concept{j2 model}
\concept{p2 model}
\concept{b2 model}


\examples{

# create a very small network with covariates for illustrative purposes
S <- c(1,0,1,0,1,1,0,1,0,1)
REC <- (S*-1)+1
D1 <- matrix(c(0,1,0,1,0,1,0,1,0,1,
              0,0,0,1,0,1,0,1,0,1,
              1,1,0,0,1,0,0,0,0,0,
              1,1,1,0,1,0,0,0,0,1,
              1,0,1,0,0,1,1,0,1,0,
              0,0,0,0,0,0,1,1,1,1,
              0,0,0,0,0,1,0,1,0,1,
              1,0,0,0,0,1,1,0,1,1,
              0,1,0,1,0,1,0,1,0,0,
              1,0,1,1,1,0,0,0,0,0), ncol=10)
D2 <- abs(matrix(rep(S,10), byrow = FALSE, ncol= 10) -
            matrix(rep(REC,10), byrow = TRUE, ncol= 10))
R <- D1*t(D1)
Y <- matrix(c(0,1,1,1,1,1,0,0,1,1,
              0,0,0,1,1,1,0,0,1,0,
              1,1,0,1,1,1,0,0,1,1,
              1,1,1,0,1,1,0,1,1,0,
              1,1,1,1,0,1,1,0,1,1,
              0,1,1,1,1,0,1,1,1,0,
              1,0,1,0,1,1,0,1,0,1,
              0,1,1,1,0,1,1,0,1,1,
              1,0,1,0,1,0,1,1,0,1,
              1,1,1,0,0,1,1,1,1,0), ncol=10) 

# estimate p2 model
p2(Y,sender= ~ S, receiver =  ~ REC, density = ~ D1 + D2, reciprocity= ~ R,
   burnin = 100, sample = 400, adapt = 10)
# Notice: burn-in, sample size and number of adaptive sequenses are 
# much smaller than recommended to keep computation time low.
# recommended code: 
\dontrun{
p2(Y,sender= ~ S, receiver =  ~ REC, density = ~ D1 + D2, reciprocity= ~ R)
}
}
