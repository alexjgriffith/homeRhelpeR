% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifUtils.r
\name{reverseIUPAC}
\alias{reverseIUPAC}
\title{reverse IUPAC}
\usage{
reverseIUPAC(stringIn)
}
\arguments{
\item{stringIn}{a string of iupac characters}
}
\value{
returns a string equal in length to stringIn
}
\description{
Takes one set of iupac characters and returns the reverse compliment
}
\examples{
motif<-"CANN"
reverseIUPAC("CANN")
# > NNTG
}

