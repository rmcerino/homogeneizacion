\name{valor_actualizado}
\alias{valor_actualizado}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Actualización de valores muestrales a un determinado valor de tipo de cambio.
}
\description{
Paso 2 de 3:

Con las elasticidades y la base de datos obtenidas en el paso 1, se procede a la actualización de valores muestrales a un determinado valor de tipo de cambio.
}
\usage{
valor_actualizado(tc_act, datos, elasticidad)
}
\arguments{
  \item{tc_act}{Valor de tipo de cambio al cual se desea actualizar}
  \item{datos}{Base de datos muestral obtenida de la función "elasticidades", es importante que tenga definidos los cuantiles de cada observación}
  \item{elasticidad}{Data frame obtenido de la función "elasticidades", la cual posee las elasticidades por cada cuantil}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
La función arroja la base de datos con dos nuevas variables:

elasticidad = Elasticidad de tipo de cambio en cada observación

valor_actualizado = Es el valor en pesos por metro cuadrado llevado al tipo de cambio de referencia

\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Este es el segundo de tres pasos para actualizar y homogeneizar la base muestral.

Primero: elasticidades(datos, cuantiles, form, dist_lw)
Segundo: valor_actualizado(tc_act, datos, elasticidad)
Tercero: func_homog(form, datos, dist_lw, p_valor, parcelas)

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
