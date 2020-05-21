\name{func_homog}
\alias{func_homog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Función de homogeneización}
\description{
Paso 3 de 3:

Cálculo de los coeficientes de homogeneización (a nivel muestral y parcelario) y de los valores homogeneizados por metro cuadrado de la base muestral.
}
\usage{
func_homog(form, datos, dist_lw, p_valor, parcelas)
}
\arguments{
  \item{form}{Fórmula con la que se calcularán los coeficientes. Es importante que siga la siguiente estructura:

  log(valor_actualizado) ~ log(p_sup) + log(largo_frente) + forma + ubicacion_cuadra + p_tipodevalor + p_sj + otras variables
    }
  \item{datos}{Base de datos muestrales con valores actualizados del paso 2. Es importante que contenga además las variables superficie, frente, forma, ubicación en la cuadra esquina, interno, salida a dos calles, tipo de valor, situación jurídica y el resto de las variables definidas en la fórmula}
  \item{dist_lw}{Distancia para generar la matriz W de vecinos}
  \item{p_valor}{P_valor referente al nivel de significación que se utilizará para admitir las variables en las regresiones}
  \item{parcelas}{Base parcelaria que será utilizada para la predicción. Es importante que contenga las variables superficie, frente, forma, ubicación en la cuadra esquina, interno, salida a dos calles}
    }
\details{
En la fórmula es necesario que se respeten los nombres de las variables y su contenido:

p_sup = Es la superficie en metros cuadrados (variable numérica)
largo_frente = Es el largo (o ancho) del frente en metros (variable numérica)
forma = Debe ser de clase "factor" donde 0 indique forma regular y 1 irregular
ubicacion_cuadra = Debe ser de clase "factor" donde 0 indique ubicación en la cuadra medial, 1 esquina, 2 interno, 3 salida a dos o más calles
p_tipodevalor = Debe ser de clase "factor" donde 0 indique venta y 1 oferta
p_sj = Debe ser de clase "factor" donde 0 indique con escritura y 1 sin escritura

NO CAMBIAR LOS NOMBRES, DE LO CONTRARIO LA FUNCIÓN ARROJARÁ ERROR

}
\value{

En el directorio especificado se creará una carpeta con la siguiente información:

Valores de los parámetros superficie, frente, forma, ubicación en la cuadra esquina, interno, salida a dos calles, tipo de valor y situación jurídica si existiesen.

Valor de la mediana de la superficie y el frente del parcelario.

Una base de datos denominada "datos_coef" en .Rda que contendrá los coeficientes de homogeneización y los valores homogeneizados para cada observación muestral.

Una base denominada "parcelas_coef" en .Rda, donde a la base parcelaria se le agregan los coeficientes de homogeneización de cada observación.

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Este es el tercer y último paso para actualizar y homogeneizar la base muestral.

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
