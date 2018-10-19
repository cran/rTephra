\name{VILL}
\alias{VILL}
\docType{data}
\title{Villarrica DEM}
\description{
Digital Elevation Model (DEM) of Villarrica Volcano, Chile.
}
\format{
  List with following elements:
  \describe{
    \item{ x   }{ Easting from vent (m)}
    \item{ y  }{ Northing from vent (m)}
    \item{ z  }{ Elevation above sea level (m)}
    \item{ lat }{ Latitude (degrees)}
    \item{ lon }{ Longitude (degrees)}
  }
}

\examples{
data(VILL)

contour(VILL)
}
\keyword{datasets}
