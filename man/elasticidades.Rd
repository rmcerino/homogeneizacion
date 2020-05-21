\name{elasticidades}
\alias{elasticidades}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cálculo de elasticidades del tipo de cambio
}
\description{
Paso 1 de 3:
Estimación de las elasticidades del tipo de cambio, de acuerdo a una determinada cantidad de cuantiles. Que serán luego utilizadas para la actualización de valores muestrales.
}
\usage{
elasticidades(datos, cuantiles, form, dist_lw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datos}{Base de datos muestrales (tiene que ser tipo sf)}
  \item{cuantiles}{Número de cuantiles que se desean estimar (mayor o igual a 1)}
  \item{form}{Fórmula para la estimación de las elasticidades.

  Si cuantiles = 1, la fórmula debe ser:

  log(valor_m2) ~ log(tc) + otras variables

  Si cuantiles > 1, la fórmula debe ser:

  log(valor_m2) ~ log(tc) + log(tc):quant + otras variables}

  \item{dist_lw}{Distancia para generar la matriz W de vecinos}


}
\details{
Es importante que la variable tc sea el tipo de cambio en el momento que se relevó el dato.
Y es importante llamar a esa variable "tc" en la base muestral, sino no funciona.
}

\value{

Devuelve un data frame llamado elasticidad, donde se presenta la elasticidad estimada para cada cuantil.

Devuelve una base (llamada datos), donde a la base datos original se le ha incorporado las variable quant (cuantil al que pertenece la observación)
}

\references{
%% ~put references to the literature/web site here ~
}

\author{
%%  ~~who you are~~
}

\note{

Este es el primer paso de tres para actualizar y homogeneizar la base muestral:

Primero: elasticidades(datos, cuantiles, form, dist_lw)
Segundo: valor_actualizado(tc_act, datos, elasticidad)
Tercero: func_homog(form, datos, dist_lw, p_valor, parcelas)

}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
