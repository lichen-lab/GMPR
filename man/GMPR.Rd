\name{GMPR}
\alias{GMPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geometric Mean of Pairwise Ratios (GMPR) of microbime data
}
\description{
Wrapper function of GMPR with appropriate datatype check and exception handling.
    A simple but effective normalization method - for zeroinflated
    sequencing data such as microbiome-Seq data.
}
\usage{
GMPR(OTUmatrix, min_ct = 2, intersect_no = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OTUmatrix}{
An OTU table matrix, where OTUs arranged in columns and samples in rows.
}
  \item{min_ct}{
The minimal number of OTUs counts. Default 2.
}
  \item{intersect_no}{
The minimal number of shared OTUs between samples. Default 4.
}
}
\details{
Normalization is the first and a critical step in microbiome sequencing data analysis to account for variable library sizes. Microbiome data contains a vast number of zeros due to the physical absence or under-sampling of the microbes. Here we propose GMPR - a simple but effective normalization method - for zero-inflated sequencing data such as microbiome data. Simulation studies and real datasets analyses demonstrate that the proposed method is more robust than competing methods, leading to more powerful detection of differentially abundant taxa and higher reproducibility of the relative abundances of taxa.

GMPR normalizes data following two steps:

First, to calculate rij, the median count ratio of nonzero counts between samples:
  rij=median(cki/ckj) (k in 1:OTU_number and cki, ckj is the non-zero count of the kth OTU)

Secondly, to caculate the size factor si for a given sample i:
  si=geometric_mean(rij)
}
\value{
A vector of geometric mean of pairwise ratios (GMPR) for each sample.
}
\references{
Li Chen, James Reeve, Lujun Zhang, Shenbing Huang, and Jun
Chen. GMPR: A robust normalization method for zero-inflated count data with application to microbiome sequencing data
}
\author{
Chen Jun; Zhang Lujun.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(GMPR)
data(otutable)
GMPR(t(otutable))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
